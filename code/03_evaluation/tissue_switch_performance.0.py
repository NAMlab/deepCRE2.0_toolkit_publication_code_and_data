"""
Purpose:
    Evaluate how accurately both SSR and MSR leaf/root models classify genes whose expression
    class flips between tissues against each model's own tissue labels (accuracy/MCC at
    threshold 0.5 and threshold-free auROC across balanced draws).

Manuscript (Peleke et al., 2026):
    Generates `results/evaluation/tissue_cross/tissue_switch_performance.tsv` plotted in
    Supplementary Figure 2m,o.

Usage & I/O (run from repo root):
    python3 code/03_evaluation/tissue_switch_performance.0.py <ssr_pred_dir> <msr_pred_dir> <out.tsv>

How to change or extend:
    Imports loaders from `tissue_cross_eval.1.py` (`T`) and `tissue_cross_msr.0.py` (`M`);
    edit metric thresholds or model sets inside `main()`.
"""
import importlib.util
import os
import re
import sys
import zlib

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location("tcross", os.path.join(HERE, "tissue_cross_eval.1.py"))
T = importlib.util.module_from_spec(spec)
spec.loader.exec_module(T)
specm = importlib.util.spec_from_file_location("tmsr", os.path.join(HERE, "tissue_cross_msr.0.py"))
M = importlib.util.module_from_spec(specm)
specm.loader.exec_module(M)
ORDER = T.ORDER
MSR_RE = re.compile(r"^MSR_dPm_all_(leaf|root)_seed2026_([A-Z][a-z]{3})_dPm_(leaf|root)_"
                    r"msr_train_models_\d{6}_\d{6}\.h5$")


def main(ssr_dir, msr_dir, out):
    rows = []
    for code in ORDER:
        shared = T.shared_genes(ssr_dir, code)
        lab = {t: T.labels(code, t) for t in ("leaf", "root")}
        j = pd.concat([lab["leaf"].rename("leaf"), lab["root"].rename("root")], axis=1, join="inner")
        j = j[j.leaf.isin([0, 1]) & j.root.isin([0, 1])]
        j = j[j.index.isin(shared)]
        flip = j[j.leaf != j.root]
        for model_tissue in ("leaf", "root"):
            y = flip[model_tissue].to_numpy().astype(int)
            preds = {}
            d = T.predictions(ssr_dir, code, model_tissue, shared).set_index("genes")
            col = [c for c in d.columns if c.startswith(code) and c.endswith(".h5")]
            assert len(col) == 1, (code, model_tissue)
            preds["SSR"] = d[col[0]]
            dm = M.msr_predictions(msr_dir, code, model_tissue, shared).set_index("genes")
            colm = [c for c in dm.columns if MSR_RE.match(c) and MSR_RE.match(c).group(2) == code]
            assert len(colm) == 1, (code, model_tissue, colm)
            preds["MSR (all 13 species)"] = dm[colm[0]]
            for kind, p in preds.items():
                pv = p.reindex(flip.index).to_numpy()
                keep = ~np.isnan(pv)
                r = dict(species=code, model_type=kind, model_tissue=model_tissue,
                         n_flipped=int(keep.sum()),
                         n_high_in_this_tissue=int((y[keep] == 1).sum()),
                         n_low_in_this_tissue=int((y[keep] == 0).sum()),
                         mean_P=round(float(np.nanmean(pv)), 4))
                r.update(T.balanced(y[keep], pv[keep],
                                    zlib.crc32(f"flip|{kind}|{code}|{model_tissue}".encode())))
                rows.append(r)
        print(f"  {code}: {len(flip)} flipped genes scored for both model types and tissues")
    d = pd.DataFrame(rows)
    d["_s"] = d.species.map(ORDER.index)
    d = d.sort_values(["model_type", "model_tissue", "_s"]).drop(columns="_s")
    d.to_csv(out, sep="\t", index=False)
    print("\nPerformance on the flipped genes, scored against the model's own tissue "
          "(chance = 0.5 accuracy, 0.5 auROC, 0 MCC):")
    for kind in ("SSR", "MSR (all 13 species)"):
        for t in ("leaf", "root"):
            s = d[(d.model_type == kind) & (d.model_tissue == t)]
            print(f"  {kind:22s} {t:5s} median accuracy {s.accuracy.median():.3f}  "
                  f"auROC {s.auROC.median():.3f}  MCC {s.MCC.median():+.3f}  "
                  f"(auROC > 0.5 in {int((s.auROC > 0.5).sum())}/{len(s)} species)")
    print("\nwritten:", out)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2], sys.argv[3])

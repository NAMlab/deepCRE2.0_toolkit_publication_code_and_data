#!/usr/bin/env python3
"""
Purpose:
    Core publication plotting module providing data readers, UPGMA/Mantel distance trees,
    cladograms, SSR/MSR dot plots, heatmap grids, and standalone split figure builders across all
    13 species (including *H. vulgare*). Note: `build_figure1_setup()` expects `EXTERNAL/figure_inputs/fig1_ab_crop.png`.

Manuscript (Peleke et al., 2026):
    Shared computational and graphical library imported by `make_fig2_comprehensive.{0,2}.py`,
    `make_suppfig1_root.2.py`, `make_suppfig2_tissue.7.py`, and `build_source_data_fig2.5.py`.

Usage & I/O (run from repo root; defaults to `results/evaluation/`):
    FIG_OUT_DIR=/path/to/out python3 code/04_figures/make_pub_figures_split.7.py
    - Reads `results/evaluation/final/` and `results/evaluation/msr/`
    - Writes split figures and companion `_values.tsv` / `_legend.txt` into `$FIG_OUT_DIR`

How to change or extend:
    Modify `ORDER`, `LABELS`, `TOPO_ALL`, `TOPO_DIC`, `TOPO_GRA`, `CLADES`, color map `CMAP`, or
    font sizes cleanly in one place; all downstream Figure 2 and Supp. Fig. 1–2 scripts inherit them.
"""
import csv
import os
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from PIL import Image, ImageDraw, ImageFont

WORK_DIR = Path("EXTERNAL/figure_output/figures_working_20261015")  # EXTERNAL (figure_output): see EXTERNAL_INPUTS.md
BASE = Path("results/evaluation")
FINAL = BASE / "final"
MSR_DIR = BASE / "msr"
OUT = Path(os.environ.get("FIG_OUT_DIR", str(WORK_DIR / "figures")))
OUT.mkdir(parents=True, exist_ok=True)

MM = 1.0 / 25.4
ORDER = ["Atha", "Csat", "Bole", "Cucu", "Gmax", "Dcar", "Slyc", "Cqui", "Bvul", "Osat", "Hvul", "Sbic", "Zmay"]
LABELS = dict(zip(
    ORDER,
    ["A. thaliana", "Ca. sativa", "B. oleracea", "Cu. sativus", "G. max", "D. carota",
     "S. lycopersicum", "C. quinoa", "B. vulgaris", "O. sativa", "H. vulgare", "S. bicolor", "Z. mays"]
))
COLLAB = {c: (c if c != "Dcar" else "Dcau") for c in ORDER}
TOPO_ALL = ((((("Atha", "Csat"), "Bole"), ("Cucu", "Gmax")), ("Dcar", "Slyc")), ("Cqui", "Bvul")), (("Osat", "Hvul"), ("Sbic", "Zmay"))
TOPO_DIC = TOPO_ALL[0]
TOPO_GRA = TOPO_ALL[1]
CLADES = [(0, 5, "p"), (5, 7, "D"), (7, 9, "s"), (9, 13, "^")]

CMAP = LinearSegmentedColormap.from_list("cyan_white_orange", ["#00CDCD", "#FFFFFF", "#EE9A00"], N=256)
LIGHT, MID, DARK, GRID = "#D4D4D4", "#7A7A7A", "#1A1A1A", "#9A9A9A"
XMIN = 0.5
FS, FS_S, FS_L = 5.8, 5.2, 8.0
CELL = 33.0 / 13  # mm per heatmap cell (13 cells = 33.0 mm, same heatmap width as v5)

# Horizontal layout (mm) across 183 mm Nature double-column width
X_TREE = (0.0, 24.5)
X_AUROC = (26.0, 18.0)    # 26.0 .. 44.0 (+ value label up to ~49.5)
X_ACC = (51.0, 18.0)      # 51.0 .. 69.0 (+ value label up to ~74.5)
X_HEAT0 = 76.0            # 76.0 .. 109.0 (33.0 mm) and 111.5 .. 144.5 (33.0 mm)
X_CB = 146.3              # colorbar at 146.3 .. 148.3 (+ label up to ~155.5)
X_PT = 156.5              # UPGMA performance similarity tree at 156.5 .. 183.0
FIG_W = 183.0


def setup_rc():
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Liberation Sans", "DejaVu Sans", "Arial"],
        "font.size": FS,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.linewidth": 0.5,
        "lines.linewidth": 0.6,
    })


# ---------------------------------------------------------------------------
# Data loading & pure-NumPy UPGMA / Mantel helpers
# ---------------------------------------------------------------------------
def read_tsv_dicts(path: Path):
    with open(path, "r", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def read_matrix_tsv(path: Path):
    with open(path, "r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        cols = header[1:]
        data = {}
        for row in reader:
            if not row:
                continue
            rname = row[0]
            data[rname] = {c: (float(v) if v not in ("", "nan", "NA") else np.nan)
                           for c, v in zip(cols, row[1:])}
    return data, cols


def reindex_matrix(data_dict, rows, cols):
    arr = np.full((len(rows), len(cols)), np.nan, dtype=float)
    for i, r in enumerate(rows):
        row_d = data_dict.get(r, {})
        for j, c in enumerate(cols):
            arr[i, j] = row_d.get(c, np.nan)
    return arr


def write_tsv(path: Path, rows, fieldnames=None):
    if not rows:
        return
    if fieldnames is None:
        fieldnames = list(rows[0].keys())
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        for r in rows:
            out_r = {}
            for k in fieldnames:
                v = r.get(k, "")
                if isinstance(v, (float, np.floating)):
                    out_r[k] = "" if np.isnan(v) else f"{float(v):.4f}"
                else:
                    out_r[k] = v
            w.writerow(out_r)


def tree_distance(topo):
    """Topological distance (edges) between all tips of a nested-tuple tree."""
    paths = {}

    def walk(node, path):
        if isinstance(node, str):
            paths[node] = path
        else:
            for ch in node:
                walk(ch, path + (id(node),))

    walk(topo, ())

    def d(a, b):
        pa, pb = paths[a], paths[b]
        n = 0
        while n < min(len(pa), len(pb)) and pa[n] == pb[n]:
            n += 1
        return (len(pa) - n + 1) + (len(pb) - n + 1) if a != b else 0

    return np.array([[d(a, b) for b in ORDER] for a in ORDER], dtype=float)


def mantel(d1, d2, n_perm=9999, seed=2026):
    iu = np.triu_indices_from(d1, 1)
    r0 = float(np.corrcoef(d1[iu], d2[iu])[0, 1])
    rng = np.random.default_rng(seed)
    hits = 0
    for _ in range(n_perm):
        p = rng.permutation(len(d1))
        hits += int(np.corrcoef(d1[p][:, p][iu], d2[iu])[0, 1] >= r0)
    return r0, (hits + 1) / (n_perm + 1)


class ClusterNode:
    def __init__(self, nid, left=None, right=None, dist=0.0):
        self.id = nid
        self.left = left
        self.right = right
        self.dist = dist

    def is_leaf(self):
        return self.left is None and self.right is None

    def get_left(self):
        return self.left

    def get_right(self):
        return self.right

    def pre_order(self):
        if self.is_leaf():
            return [self.id]
        return self.left.pre_order() + self.right.pre_order()


def upgma_to_tree(D_square: np.ndarray) -> ClusterNode:
    """Pure-NumPy average linkage (UPGMA) returning a ClusterNode root."""
    n = D_square.shape[0]
    D = D_square.astype(float).copy()
    np.fill_diagonal(D, np.inf)
    nodes = {i: ClusterNode(i) for i in range(n)}
    members = {i: [i] for i in range(n)}
    next_id = n
    for _ in range(n - 1):
        keys = sorted(nodes.keys())
        best_d = np.inf
        best_pair = None
        for ia in range(len(keys)):
            a = keys[ia]
            for ib in range(ia + 1, len(keys)):
                b = keys[ib]
                if D[a, b] < best_d:
                    best_d = D[a, b]
                    best_pair = (a, b)
        a, b = best_pair
        na, nb = nodes.pop(a), nodes.pop(b)
        ma, mb = members.pop(a), members.pop(b)
        new_m = ma + mb
        nodes[next_id] = ClusterNode(next_id, left=na, right=nb, dist=float(best_d))
        members[next_id] = new_m
        new_D = np.full((next_id + 1, next_id + 1), np.inf)
        new_D[:next_id, :next_id] = D
        for k, mk in members.items():
            if k == next_id:
                continue
            d_ak = float(np.mean([D_square[i, j] for i in new_m for j in mk]))
            new_D[next_id, k] = new_D[k, next_id] = d_ak
        D = new_D
        next_id += 1
    return nodes[next_id - 1]


def compute_msr_summary(full_rows, pf_rows, tissue):
    """
    Compute per-species MSR and SSR reference metrics for both auROC and accuracy.
    v7: H. vulgare included.
    """
    # all rows (v7: H. vulgare included)
    msr_rows = [
        r for r in full_rows
        if r["kind"] == "MSR"
    ]
    oth_rows = [
        r for r in full_rows
        if r["kind"] == "SSR_deployed" and r["status"] == "other_species"
    ]
    dep_map = {}
    for r in pf_rows:
        if r["tissue"] == tissue and str(r["deployed"]).lower() == "true" and r["code"] in ORDER:
            dep_map[r["code"]] = r

    summary = {}
    for code in ORDER:
        dep = dep_map[code]
        oth_c = [r for r in oth_rows if r["target"] == code]
        r_out = {
            "target": code,
            "SSR_deployed_own_auROC": float(dep["auROC_v8"]),
            "SSR_deployed_own_accuracy": float(dep["accuracy_v8"]),
            "SSR_deployed_own_MCC": float(dep["MCC_v8"]),
            "SSR_other_species_median_auROC": float(np.median([float(x["auROC"]) for x in oth_c])),
            "SSR_other_species_median_accuracy": float(np.median([float(x["accuracy"]) for x in oth_c])),
            "SSR_other_species_median_MCC": float(np.median([float(x["MCC"]) for x in oth_c])),
        }
        for key in ["all", "dicots", "monocots"]:
            m = [x for x in msr_rows if x["target"] == code and x["set"] == key]
            h = [x for x in m if x["status"] == "held_out"]
            o = [x for x in m if x["status"] == "outside"]
            if h:
                assert len(h) == 1, (code, key, len(h))
                r_out[f"{key}_kind"] = "held_out"
                r_out[f"{key}_auROC"] = float(h[0]["auROC"])
                r_out[f"{key}_accuracy"] = float(h[0]["accuracy"])
                r_out[f"{key}_MCC"] = float(h[0]["MCC"])
                r_out[f"{key}_n_models"] = 1
            elif o:
                # 4 grass models for dicots outside_median, 9 eudicot models for monocots outside_median
                expected_n = 4 if key == "monocots" else 9
                assert len(o) == expected_n, f"Expected {expected_n} models for {code}/{key}, got {len(o)}"
                r_out[f"{key}_kind"] = "outside_median"
                r_out[f"{key}_auROC"] = float(np.median([float(x["auROC"]) for x in o]))
                r_out[f"{key}_accuracy"] = float(np.median([float(x["accuracy"]) for x in o]))
                r_out[f"{key}_MCC"] = float(np.median([float(x["MCC"]) for x in o]))
                r_out[f"{key}_n_models"] = len(o)
        summary[code] = r_out
    return summary, oth_rows


# ---------------------------------------------------------------------------
# Plotting primitives
# ---------------------------------------------------------------------------
def draw_tree(ax, topo, rows):
    y = {c: i for i, c in enumerate(rows)}
    depth = lambda n: 0 if isinstance(n, str) else 1 + max(depth(c) for c in n)
    D = depth(topo)

    def draw(n):
        if isinstance(n, str):
            return 1.0, y[n]
        kids = [draw(c) for c in n]
        x = 1.0 - depth(n) * 0.9 / max(D, 1)
        ys = [k[1] for k in kids]
        ax.plot([x, x], [min(ys), max(ys)], color=DARK, lw=0.7)
        for kx, ky in kids:
            ax.plot([x, kx], [ky, ky], color=DARK, lw=0.7, solid_capstyle="butt")
        yy = ys[1] if len(kids) == 3 else (min(ys) + max(ys)) / 2
        return x, yy

    x, yy = draw(topo)
    ax.plot([0.0, x], [yy, yy], color=DARK, lw=0.7)
    for c in rows:
        ax.text(1.05, y[c], LABELS[c], va="center", style="italic", fontsize=FS)
    ax.set_ylim(len(rows) - 0.5, -0.5)
    ax.set_xlim(0, 2.25)
    ax.axis("off")


def metric_axis(ax, xlabel="auROC"):
    ticks = [0.5, 0.75, 1.0]
    ax.set_xlim(XMIN, 1.0)
    ax.set_xticks(ticks)
    ax.set_xticklabels(["0.5", "0.75", "1"], fontsize=FS_S)
    ax.set_xticks([0.625, 0.875], minor=True)
    ax.set_yticks([])
    for sp in ["top", "right", "left"]:
        ax.spines[sp].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.5)
    ax.tick_params(axis="x", which="both", length=1.5, width=0.5, pad=1)
    ax.set_xlabel(xlabel, fontsize=FS, labelpad=1)


def value_label(ax, i, x):
    ax.text(1.015, i, f"{x:.2f}", fontsize=4.9, va="center", ha="left",
            transform=ax.get_yaxis_transform(), clip_on=False)


def ssr_dots(ax, pf_rows, tissue, rows, metric_col, xlabel, values_out, panel_letter):
    """
    Dot plot for single-species models (auROC or Accuracy):
    - Grey dots: individual converged chromosome folds
    - Black dot: deployed web-tool model (selected by peak validation auROC) + printed value
    - Open circle + horizontal bar below: mean +/- s.d. of the folds
    """
    rng = np.random.default_rng(1 if "auROC" in metric_col else 11)
    for i, code in enumerate(rows):
        s = [r for r in pf_rows if r["tissue"] == tissue and r["code"] == code]
        vv = np.array([float(r[metric_col]) for r in s if r[metric_col] not in ("", "nan", "NA")], dtype=float)
        dep = [float(r[metric_col]) for r in s if str(r["deployed"]).lower() == "true"]
        assert len(dep) == 1, (tissue, code)
        if metric_col == "auROC_v8":
            assert np.isclose(dep[0], vv.max()), (tissue, code)
        m = float(vv.mean())
        sd = float(vv.std(ddof=1)) if len(vv) > 1 else 0.0
        ax.scatter(vv, i - 0.12 + rng.uniform(-0.1, 0.1, len(vv)),
                   s=2.2, color=MID, zorder=3, linewidths=0)
        ax.scatter(dep, [i - 0.12], s=9, color=DARK, zorder=4, linewidths=0)
        ax.plot([m - sd, m + sd], [i + 0.28] * 2, color=DARK, lw=0.6, solid_capstyle="butt", zorder=3)
        ax.scatter([m], [i + 0.28], s=5, facecolor="white", edgecolor=DARK, linewidths=0.5, zorder=4)
        value_label(ax, i, dep[0])
        assert vv.min() >= XMIN, f"{tissue} {code} {metric_col} min={vv.min()} < {XMIN}"
        values_out.append({
            "panel": panel_letter,
            "tissue": tissue,
            "species": code,
            "metric": xlabel,
            "web_tool_val": dep[0],
            "fold_mean": m,
            "fold_sd": sd,
            "fold_min": float(vv.min()),
            "fold_max": float(vv.max()),
            "n_folds": len(vv),
        })
    ax.set_ylim(len(rows) - 0.5, -0.5)
    metric_axis(ax, xlabel=xlabel)


def msr_dots(ax, summary, oth_rows, rows, key, cross, metric_name, xlabel, values_out, panel_letter):
    """
    Dot plot for multi-species models (auROC or Accuracy):
    - Light grey dots: single-species models of the other 12 species
    - Open diamond: web-tool single-species model of this species
    - Grey square (clade-only blocks): outside-clade multi-species models (median)
    - Black dot: multi-species model on its held-out species + printed value
    """
    rng = np.random.default_rng(2 if metric_name == "auROC" else 12)
    for i, code in enumerate(rows):
        r = summary[code]
        assert r[f"{key}_kind"] == "held_out", (code, key)
        o = np.array([float(x[metric_name]) for x in oth_rows if x["target"] == code], dtype=float)
        ax.scatter(o, i + rng.uniform(-0.22, 0.22, len(o)),
                   s=1.8, color=LIGHT, edgecolor=MID, linewidths=0.2, zorder=2)
        own_val = r[f"SSR_deployed_own_{metric_name}"]
        ax.scatter([own_val], [i], s=8, marker="D", facecolor="white",
                   edgecolor=DARK, linewidths=0.6, zorder=3)
        if cross is not None:
            assert r[f"{cross}_kind"] == "outside_median", (code, cross)
            ax.scatter([r[f"{cross}_{metric_name}"]], [i], s=9, marker="s",
                       facecolor=MID, edgecolor=MID, linewidths=0, zorder=3)
        ho_val = r[f"{key}_{metric_name}"]
        ax.scatter([ho_val], [i], s=11, color=DARK, zorder=4, linewidths=0)
        value_label(ax, i, ho_val)
        values_out.append({
            "panel": panel_letter,
            "species": code,
            "msr_set": key,
            "metric": xlabel,
            "msr_held_out_val": ho_val,
            "other_clade_median_val": r[f"{cross}_{metric_name}"] if cross else np.nan,
            "web_tool_own_val": own_val,
            "single_species_other_median_val": float(np.median(o)),
            "msr_held_out_auROC": r[f"{key}_auROC"],
            "msr_held_out_accuracy": r[f"{key}_accuracy"],
            "msr_held_out_MCC": r[f"{key}_MCC"],
        })
    ax.set_ylim(len(rows) - 0.5, -0.5)
    metric_axis(ax, xlabel=xlabel)


def draw_heat(ax, m_dict, rows, cols, title):
    a = reindex_matrix(m_dict, rows, cols)
    assert np.isfinite(a).all(), f"empty heatmap cells in {title}"
    ax.imshow(a, cmap=CMAP, vmin=0, vmax=1, aspect="auto", interpolation="nearest")
    nr, nc = len(rows), len(cols)
    for k in range(nr + 1):
        ax.plot([-0.5, nc - 0.5], [k - 0.5, k - 0.5], color=GRID, lw=0.3)
    for k in range(nc + 1):
        ax.plot([k - 0.5, k - 0.5], [-0.5, nr - 0.5], color=GRID, lw=0.3)
    if all(r in cols for r in rows):
        c = [cols.index(r) for r in rows]
        ax.plot([c[0] - 0.5, c[-1] + 0.5], [-0.5, nr - 0.5], color="#505050", lw=0.7)
    for a0, a1, _ in CLADES:
        cl = ORDER[a0:a1]
        ri = [i for i, r in enumerate(rows) if r in cl]
        ci = [j for j, c in enumerate(cols) if c in cl]
        if ri and ci:
            ax.add_patch(Rectangle((ci[0] - 0.5, ri[0] - 0.5), len(ci), len(ri),
                                   fill=False, ec="white", lw=0.7, ls=(0, (2, 1.2))))
    ax.set_xticks(range(nc))
    ax.set_xticklabels([COLLAB[c] for c in cols], rotation=90, fontsize=FS_S)
    ax.set_yticks([])
    ax.tick_params(length=0, pad=1)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_title(title, fontsize=FS, loc="left", pad=2)
    ax.set_xlim(-0.5, nc - 0.5)
    ax.set_ylim(nr - 0.5, -0.5)


def draw_brackets(ax, cols):
    names = dict(zip([0, 5, 7, 9], ["rosids", "asterids", "Caryophyllales", "grasses"]))
    present = []
    for a0, a1, _ in CLADES:
        ci = [j for j, c in enumerate(cols) if c in ORDER[a0:a1]]
        if ci:
            present.append((a0, ci[0], ci[-1]))
    has_ast = any(a0 == 5 for a0, _, _ in present)
    for a0, j0, j1 in present:
        lab = names[a0]
        ax.plot([j0 - 0.42, j0 - 0.42, j1 + 0.42, j1 + 0.42], [0.0, 0.3, 0.3, 0.0], color=MID, lw=0.5)
        low = lab == "Caryophyllales" and has_ast
        ax.text((j0 + j1) / 2, 0.82 if low else 0.36, lab, fontsize=4.8, ha="center", va="top")
    eud = [j for j, c in enumerate(cols) if c in ORDER[:9]]
    if len([p for p in present if p[0] < 9]) > 1:
        ax.plot([eud[0] - 0.42, eud[0] - 0.42, eud[-1] + 0.42, eud[-1] + 0.42],
                [1.3, 1.55, 1.55, 1.3], color=MID, lw=0.5)
        ax.text((eud[0] + eud[-1]) / 2, 1.63, "eudicots", fontsize=4.8, ha="center", va="top")
    ax.set_xlim(-0.5, len(cols) - 0.5)
    ax.set_ylim(2.05, -0.1)
    ax.axis("off")


def draw_perf_tree(ax, m0_dict, m1_dict):
    a0 = reindex_matrix(m0_dict, ORDER, ORDER)
    a1 = reindex_matrix(m1_dict, ORDER, ORDER)
    Xp = np.hstack([a0, a1, a0.T, a1.T])
    n = len(ORDER)
    D_sq = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(n):
            D_sq[i, j] = float(np.linalg.norm(Xp[i] - Xp[j]))
    root = upgma_to_tree(D_sq)
    pos = lambda nd: (np.mean(nd.pre_order()), min(nd.pre_order()))
    order = []

    def lay(nd):
        if nd.is_leaf():
            order.append(nd.id)
            return
        for ch in sorted([nd.get_left(), nd.get_right()], key=pos):
            lay(ch)

    lay(root)
    # v7: with 13 species the top leaf of the UPGMA performance tree is no longer A. thaliana in
    # every tissue (it is data, not a layout error), so this is reported instead of asserted.
    # The full leaf order is returned below and written into the figure legend.
    if ORDER[order[0]] != "Atha":
        print(f"  note: UPGMA performance tree starts with {ORDER[order[0]]}, not Atha")
    y = {leaf: i for i, leaf in enumerate(order)}
    H = root.dist

    def draw(nd):
        if nd.is_leaf():
            return 0.0, y[nd.id]
        kids = [draw(c) for c in (nd.get_left(), nd.get_right())]
        ys = [k[1] for k in kids]
        ax.plot([nd.dist, nd.dist], [min(ys), max(ys)], color=DARK, lw=0.7)
        for kx, ky in kids:
            ax.plot([nd.dist, kx], [ky, ky], color=DARK, lw=0.7, solid_capstyle="butt")
        return nd.dist, float(np.mean(ys))

    x, yy = draw(root)
    ax.plot([x, x * 1.06], [yy, yy], color=DARK, lw=0.7)
    for i, leaf in enumerate(order):
        code = ORDER[leaf]
        ax.text(-0.05 * H, i, LABELS[code], style="italic", fontsize=FS, va="center", ha="left")
    ax.set_xlim(1.06 * H, -1.25 * H)
    ax.set_ylim(len(ORDER) - 0.5, -0.5)
    ax.axis("off")
    ax.set_title("Performance similarity", fontsize=FS, loc="left", pad=2)
    r, p = mantel(D_sq, tree_distance(TOPO_ALL))
    return dict(order=[ORDER[i] for i in order], mantel_r=r, mantel_p=p)


def draw_legend(fig, items, x_mm, y_mm, FW, FH):
    for j, (hd, lab) in enumerate(items):
        y = (y_mm - j * 3.0) * MM / FH
        x = x_mm * MM / FW
        if hd is None:
            continue
        if "marker" in hd:
            fig.add_artist(Line2D([x + 1.5 * MM / FW], [y], ls="none", transform=fig.transFigure, **hd))
        elif "meansd" in hd:
            fig.add_artist(Line2D([x, x + 3.0 * MM / FW], [y, y], color=DARK, lw=0.6, transform=fig.transFigure))
            fig.add_artist(Line2D([x + 1.5 * MM / FW], [y], marker="o", ms=2.2, mfc="white",
                                  mec=DARK, mew=0.5, ls="none", transform=fig.transFigure))
        elif "line" in hd:
            fig.add_artist(Line2D([x + 1.5 * MM / FW] * 2, [y - 1.0 * MM / FH, y + 1.0 * MM / FH],
                                  color=DARK, lw=0.8, transform=fig.transFigure))
        else:
            fig.add_artist(Rectangle((x, y - 0.9 * MM / FH), 3.0 * MM / FW, 1.8 * MM / FH,
                                     transform=fig.transFigure, **hd))
        fig.text(x + 4.2 * MM / FW, y, lab, fontsize=FS, va="center")


def build_figure(blocks, out_prefix, legend_items, legend_block, legend_x_mm=None):
    """
    Build a multi-block figure where each block has:
      tree | auROC dot plot | Accuracy dot plot | [optional heatmaps] | [optional extra panel]
    """
    HEAD, FOOT, GAP = 11.0, 16.0, 3.0
    LEG = 0.0 if legend_block is not None else 4.0 + 3.0 * ((len(legend_items) + 1) // 2)
    FW_mm = FIG_W
    FH_mm = sum(HEAD + len(b["rows"]) * CELL + FOOT for b in blocks) + GAP * (len(blocks) - 1) + LEG
    FW, FH = FW_mm * MM, FH_mm * MM
    fig = plt.figure(figsize=(FW, FH))
    ax_at = lambda x0, w, y0, h: fig.add_axes([x0 / FW_mm, y0 / FH_mm, w / FW_mm, h / FH_mm])
    top, extra, chk = FH_mm, {}, []

    for k, b in enumerate(blocks):
        has_heats = b.get("show_heats", True)
        heats = b.get("heats")
        if has_heats and heats is None:
            heats = [
                (b["m0"], ORDER, "True lowly expressed genes", None),
                (b["m1"], ORDER, "True highly expressed genes", None),
            ]
        bh = len(b["rows"]) * CELL
        top -= HEAD
        y0 = top - bh

        ax_tree = ax_at(X_TREE[0], X_TREE[1], y0, bh)
        ax_auroc = ax_at(X_AUROC[0], X_AUROC[1], y0, bh)
        ax_acc = ax_at(X_ACC[0], X_ACC[1], y0, bh)

        draw_tree(ax_tree, b["topo"], b["rows"])
        b["panel_auroc"](ax_auroc)
        b["panel_acc"](ax_acc)

        A = [ax_tree, ax_auroc, ax_acc]
        xs = [X_TREE[0], X_AUROC[0], X_ACC[0]]

        if has_heats and heats:
            x = X_HEAT0
            for j, (m_dict, cols, title, group) in enumerate(heats):
                hw = len(cols) * CELL
                ax_h = ax_at(x, hw, y0, bh)
                draw_heat(ax_h, m_dict, b["rows"], cols, title)
                draw_brackets(ax_at(x, hw, y0 - FOOT + 0.5, 8.0), cols)
                if group:
                    fig.text(x / FW_mm, (top + 7.6) / FH_mm, group, fontsize=FS, va="bottom")
                A.append(ax_h)
                # Record x for panel letter according to b["heat_letter_mode"]
                if b.get("heat_letter_mode", "first_only") == "each":
                    xs.append(x)
                elif b.get("heat_letter_mode", "first_only") == "pair_first" and j % 2 == 0:
                    xs.append(x)
                elif b.get("heat_letter_mode", "first_only") == "first_only" and j == 0:
                    xs.append(x)
                gap = 2.5 if j % 2 == 0 else 6.0
                x += max(hw, 10.0 if b.get("heats") else 30.0) + gap

        for xx, lab in zip(xs, b["letters"]):
            if lab:
                fig.text(xx / FW_mm, (top + 3.2) / FH_mm, lab, fontsize=FS_L, fontweight="bold", va="bottom")

        fig.text(X_TREE[0] / FW_mm, (top + 7.6) / FH_mm, b["name"],
                 fontsize=FS + 0.5, fontweight="bold", va="bottom")

        if legend_block is not None and k == legend_block:
            lx = legend_x_mm if legend_x_mm is not None else X_CB
            draw_legend(fig, legend_items, lx, y0 + bh - 1.5, FW, FH)

        if "extra" in b:
            axp = ax_at(X_PT, FW_mm - X_PT, y0, bh)
            extra = b["extra"](axp)
            fig.text(X_PT / FW_mm, (top + 3.2) / FH_mm, b["extra_letter"],
                     fontsize=FS_L, fontweight="bold", va="bottom")

        if k == 0 and has_heats:
            cb_x = b.get("cb_x_mm", max(X_CB, x - gap + 2.0))
            axcb = ax_at(cb_x, 2.0, y0 + bh * 0.2, bh * 0.6)
            cb = fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(0, 1), cmap=CMAP),
                              cax=axcb, ticks=[0, 0.5, 1])
            cb.ax.tick_params(labelsize=FS_S, length=1.5, width=0.5, pad=1)
            cb.outline.set_linewidth(0.4)
            cb.set_label("Mean predicted probability\nof high expression", fontsize=FS_S, labelpad=2)

        chk.append((A, len(b["rows"])))
        top = y0 - FOOT - GAP

    if legend_block is None:
        half = (len(legend_items) + 1) // 2
        draw_legend(fig, legend_items[:half], 0.0, LEG - 2.5, FW, FH)
        draw_legend(fig, legend_items[half:], 92.0, LEG - 2.5, FW, FH)

    for axs, n in chk:
        rr = np.arange(n)
        ys = [a.transData.transform(np.c_[np.zeros_like(rr), rr])[:, 1] for a in axs]
        assert max(float(np.abs(y - ys[0]).max()) for y in ys[1:]) < 0.5, "rows not aligned"

    fig.savefig(str(out_prefix) + ".pdf")
    fig.savefig(str(out_prefix) + ".png", dpi=600)
    plt.close(fig)
    print(f"Saved {out_prefix} ({FW_mm:.0f} x {FH_mm:.0f} mm)")
    return extra


# ---------------------------------------------------------------------------
# Figure 1 (Panels a & b: Platform setup & deepCRE architecture)
# ---------------------------------------------------------------------------
def build_figure1_setup():
    """
    Create standalone Figure 1 (panels a & b) from the round-2 high-resolution crop:
      - Figure1_setup_original.pdf / .png (untouched panels a & b)
      - Figure1.pdf / .png (with log2Max(TPM+1) -> log10(mean TPM+1) and dDR -> gDR)
    """
    crop_path = Path("EXTERNAL/figure_inputs/fig1_ab_crop.png")  # EXTERNAL (figure_inputs): see EXTERNAL_INPUTS.md
    if not crop_path.exists():
        raise FileNotFoundError(f"Expected {crop_path} (see EXTERNAL_INPUTS.md)")

    im_orig = Image.open(crop_path).convert("RGB")
    w_px, h_px = im_orig.size
    # Save original crop as 183 mm wide PDF and PNG
    fw_in = FIG_W * MM
    fh_in = fw_in * (h_px / w_px)

    for tag, img in [("Figure1_setup_original", im_orig)]:
        png_out = OUT / f"{tag}.png"
        pdf_out = OUT / f"{tag}.pdf"
        img.save(png_out, dpi=(300, 300))
        fig = plt.figure(figsize=(fw_in, fh_in), dpi=300)
        ax = fig.add_axes([0, 0, 1, 1])
        ax.imshow(np.array(img))
        ax.axis("off")
        fig.savefig(pdf_out, dpi=300)
        plt.close(fig)
        print(f"Saved {OUT / tag} ({FIG_W:.0f} x {fh_in / MM:.0f} mm)")

    # Create updated Figure 1 with dDR -> gDR and log2Max(TPM+1) -> log10(mean TPM+1)
    arr = np.array(im_orig).copy()
    # 1. Replace 'd' in dDR (x=1468..1488, y=804..840) with 'g' from gUR (x=1084..1104, y=804..840)
    arr[804:840, 1468:1488] = arr[804:840, 1084:1104]
    im_upd = Image.fromarray(arr)
    draw = ImageDraw.Draw(im_upd)
    # 2. Replace 'log2Max(TPM+1)' (x=1355..1630, y=415..464) with 'log10(mean TPM+1)'
    draw.rectangle([1355, 415, 1630, 464], fill=(255, 255, 255))
    font_path = "/usr/share/fonts/opentype/urw-base35/URWGothic-Book.otf"
    font = ImageFont.truetype(font_path, 33)
    # Render onto a temporary image so we can match the 92% horizontal condensation of the original label
    tmp = Image.new("RGB", (340, 52), (255, 255, 255))
    tdraw = ImageDraw.Draw(tmp)
    tdraw.text((0, 6), "log10(mean TPM+1)", fill=(54, 53, 53), font=font)
    tmp_scaled = tmp.resize((int(340 * 0.90), 52), resample=Image.BICUBIC)
    im_upd.paste(tmp_scaled, (1358, 415))

    png_out = OUT / "Figure1.png"
    pdf_out = OUT / "Figure1.pdf"
    im_upd.save(png_out, dpi=(300, 300))
    fig = plt.figure(figsize=(fw_in, fh_in), dpi=300)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.imshow(np.array(im_upd))
    ax.axis("off")
    fig.savefig(pdf_out, dpi=300)
    plt.close(fig)
    print(f"Saved {OUT / 'Figure1'} ({FIG_W:.0f} x {fh_in / MM:.0f} mm)")

    with open(OUT / "Figure1_legend.txt", "w", encoding="utf-8") as fh:
        fh.write(
            "Figure 1 | Overview of the deepCRE web toolkit and convolutional neural network workflow.\n"
            "a, User workflow on the deepCRE web platform: users select a reference organism, upload a gene list, "
            "and choose a pre-trained single-species or multi-species deepCRE model (Inputs); the platform extracts "
            "flanking regulatory sequences and runs inference, saliency, and in silico mutagenesis (Processing / Modelling); "
            "and returns interactive gene-level expression class predictions, nucleotide-resolution saliency maps, and "
            "mutation impact profiles (Outputs).\n"
            "b, Schematic of deepCRE sequence extraction and model training: genes are ranked by log10(mean TPM + 1) "
            "to define lowly (bottom 25th percentile) and highly (top 25th percentile) expressed classes; 1.0 kb upstream (gUR) "
            "+ 0.5 kb 5' transcribed region (gTUR) around the transcription start site (TSS) and 0.5 kb 3' transcribed region (gTDR) "
            "+ 1.0 kb downstream (gDR) around the transcription termination site (TTS) are one-hot encoded with a 20 bp N-spacer "
            "[1, 3020, 4] and passed through three convolutional/pooling/dropout blocks to predict the probability of high expression P in [0, 1].\n"
        )


def main():
    setup_rc()
    # Note: Figure 1 (updated with live website screenshots and cleaned Panel b)
    # is built by build_updated_fig1.py so we do not overwrite it here.

    pf_all = read_tsv_dicts(FINAL / "S2v8_per_fold.tsv")
    pf = [
        r for r in pf_all
        if r["status"] in ("converged", "converged_rescue")
        and r["code"] in ORDER
        and r["auROC_v8"] not in ("", "nan", "NA")
    ]

    L_SSR = [
        (dict(marker="o", ms=3.0, mfc=DARK, mec=DARK), "web-tool model (best fold; value)"),
        (dict(marker="o", ms=1.7, mfc=MID, mec=MID), "chromosome folds"),
        (dict(meansd=True), "mean \u00b1 s.d. of folds"),
    ]
    L_MSR = [
        (dict(marker="o", ms=3.2, mfc=DARK, mec=DARK), "multi-species model, held-out (value)"),
        (dict(marker="D", ms=2.6, mfc="white", mec=DARK, mew=0.6), "web-tool model of this species"),
        (dict(marker="o", ms=1.8, mfc=LIGHT, mec=MID, mew=0.2), "single-species, other species"),
    ]
    L_CROSS = [
        (dict(marker="s", ms=2.8, mfc=MID, mec=MID), "other clade's models (median)"),
    ]

    for tissue in ["leaf", "root"]:
        full_rows = read_tsv_dicts(MSR_DIR / f"msr_eval_{tissue}_full.tsv")
        summary, oth_rows = compute_msr_summary(full_rows, pf_all, tissue)

        e0, _ = read_matrix_tsv(FINAL / f"crossspecies_mean_{tissue}_class0.tsv")
        e1, _ = read_matrix_tsv(FINAL / f"crossspecies_mean_{tissue}_class1.tsv")
        M = {}
        for s in ["all", "dicots", "monocots"]:
            for k in (0, 1):
                M[(s, k)], _ = read_matrix_tsv(MSR_DIR / f"suppfig_msr_mean_{tissue}_{s}_class{k}.tsv")

        # -------------------------------------------------------------------
        # 1) Figure 2 (panels a-g: SSR a-e + MSR f, g)
        # -------------------------------------------------------------------
        vals_fig2 = []
        ssr_block = dict(
            name="Single-species models",
            topo=TOPO_ALL,
            rows=ORDER,
            letters=["a", "b", "c", "d"],
            heat_letter_mode="first_only",
            show_heats=True,
            m0=e0,
            m1=e1,
            panel_auroc=lambda ax: ssr_dots(ax, pf, tissue, ORDER, "auROC_v8", "auROC", vals_fig2, "b"),
            panel_acc=lambda ax: ssr_dots(ax, pf, tissue, ORDER, "accuracy_v8", "Accuracy", vals_fig2, "c"),
            extra=lambda ax: draw_perf_tree(ax, e0, e1),
            extra_letter="e",
        )
        msr_block_ag = dict(
            name="Multi-species models (row: held-out species)",
            topo=TOPO_ALL,
            rows=ORDER,
            letters=["", "f", "g"],
            show_heats=False,
            panel_auroc=lambda ax: msr_dots(ax, summary, oth_rows, ORDER, "all", None, "auROC", "auROC", vals_fig2, "f"),
            panel_acc=lambda ax: msr_dots(ax, summary, oth_rows, ORDER, "all", None, "accuracy", "Accuracy", vals_fig2, "g"),
        )
        t_info = build_figure(
            [ssr_block, msr_block_ag],
            OUT / f"Figure2_{tissue}",
            L_SSR + [(None, None)] + L_MSR,
            legend_block=1,
            legend_x_mm=X_HEAT0,
        )
        write_tsv(OUT / f"Figure2_{tissue}_values.tsv", vals_fig2)
        with open(OUT / f"Figure2_{tissue}_legend.txt", "w", encoding="utf-8") as fh:
            fh.write(
                f"Figure 2 ({tissue}): Single-species and leave-one-species-out multi-species model performance across the 13 plant species (9 eudicots, 4 grasses).\n"
                f"a, Cladogram of the 13 species.\n"
                f"b, c, Single-species chromosome-fold validation auROC (b) and accuracy (c) (grey dots = individual chromosome folds; black dot = deployed web-tool model with value printed at right; open circle and horizontal bar = mean +/- s.d. across folds).\n"
                f"d, Cross-species prediction heatmaps showing the mean predicted probability of high expression for true lowly expressed genes (left) and true highly expressed genes (right) across all 13 x 13 model-target species pairs.\n"
                f"e, UPGMA tree of Euclidean distances between the 52-value cross-prediction profiles of the 13 species (rows and columns of d). Mantel test against the topological distances of the cladogram in a: r = {t_info['mantel_r']:.3f}, p = {t_info['mantel_p']:.4f} (Pearson, 9,999 permutations, 78 species pairs). Leaf order: {', '.join(t_info['order'])}.\n"
                f"f, g, Leave-one-species-out multi-species model auROC (f) and accuracy (g) on each held-out species (black dot with value at right), compared against the held-out species' own single-species web-tool model (open diamond) and single-species models transferred from the other 12 species (light grey dots).\n"
            )

        # -------------------------------------------------------------------
        # 1b) v7: Figure 2 a-g + eudicot-only (h-l) and grass-only (m-q) MSR blocks (as v6)
        # -------------------------------------------------------------------
        vals_fig2c = []
        Dq, Gq = ORDER[:9], ORDER[9:]
        blocks_c = [
            dict(ssr_block,
                 panel_auroc=lambda ax: ssr_dots(ax, pf, tissue, ORDER, "auROC_v8", "auROC", vals_fig2c, "b"),
                 panel_acc=lambda ax: ssr_dots(ax, pf, tissue, ORDER, "accuracy_v8", "Accuracy", vals_fig2c, "c")),
            dict(msr_block_ag, name="Multi-species models, all species (row: held-out species)",
                 panel_auroc=lambda ax: msr_dots(ax, summary, oth_rows, ORDER, "all", None, "auROC", "auROC", vals_fig2c, "f"),
                 panel_acc=lambda ax: msr_dots(ax, summary, oth_rows, ORDER, "all", None, "accuracy", "Accuracy", vals_fig2c, "g")),
            dict(name="Eudicot-only models (row: held-out eudicot)", topo=TOPO_DIC, rows=Dq,
                 letters=["h", "i", "j", "k", "l"], heat_letter_mode="pair_first", show_heats=True,
                 heats=[(M[("dicots", 0)], Dq, "True low", "on eudicots (own clade)"),
                        (M[("dicots", 1)], Dq, "True high", None),
                        (M[("dicots", 0)], Gq, "True low", "on grasses"),
                        (M[("dicots", 1)], Gq, "True high", None)],
                 panel_auroc=lambda ax: msr_dots(ax, summary, oth_rows, Dq, "dicots", None, "auROC", "auROC", vals_fig2c, "i"),
                 panel_acc=lambda ax: msr_dots(ax, summary, oth_rows, Dq, "dicots", None, "accuracy", "Accuracy", vals_fig2c, "j")),
            dict(name="Grass-only models (row: held-out grass)", topo=TOPO_GRA, rows=Gq,
                 letters=["m", "n", "o", "p", "q"], heat_letter_mode="pair_first", show_heats=True,
                 heats=[(M[("monocots", 0)], Gq, "True low", "on grasses (own clade)"),
                        (M[("monocots", 1)], Gq, "True high", None),
                        (M[("monocots", 0)], Dq, "True low", "on eudicots"),
                        (M[("monocots", 1)], Dq, "True high", None)],
                 panel_auroc=lambda ax: msr_dots(ax, summary, oth_rows, Gq, "monocots", None, "auROC", "auROC", vals_fig2c, "n"),
                 panel_acc=lambda ax: msr_dots(ax, summary, oth_rows, Gq, "monocots", None, "accuracy", "Accuracy", vals_fig2c, "o")),
        ]
        # Omit "other clade's models (median)" square in the combined Figure 2 legend.
        t_c = build_figure(blocks_c, OUT / f"Figure2_{tissue}_with_clade_MSR",
                           L_SSR + [(None, None)] + L_MSR, legend_block=1, legend_x_mm=X_HEAT0)
        write_tsv(OUT / f"Figure2_{tissue}_with_clade_MSR_values.tsv", vals_fig2c)
        with open(OUT / f"Figure2_{tissue}_with_clade_MSR_legend.txt", "w", encoding="utf-8") as fh:
            fh.write(
                f"Figure 2 | Flowering plant deepCRE model performance for single-species, multi-species and clade-specific training ({tissue}).\n"
                f"a, A broad phylogenetic spectrum of flowering plants and crop species was selected, encompassing rosids (Arabidopsis thaliana (thale cress; Atha), Camelina sativa (false flax; Csat), Brassica oleracea (cabbage; Bole), Cucumis sativus (cucumber; Cucu), Glycine max (soy; Gmax)), asterids (Daucus carota (carrot; Dcar), Solanum lycopersicum (tomato; Slyc)), Caryophyllales (Chenopodium quinoa (quinoa; Cqui), Beta vulgaris (sugar beet; Bvul)), and grasses/cereals (Oryza sativa (rice; Osat), Hordeum vulgare (barley; Hvul), Sorghum bicolor (sorghum; Sbic), Zea mays (maize; Zmay)).\n"
                f"b, c, Area under the receiver operating characteristic curve (auROC, b) and accuracy (c) of deepCRE models for {tissue} tissue show excellent predictive performance. Following the deepCRE single-species reference training strategy, cross-validations were performed on individual chromosomes. Individual grey dots show the performance per species on hold-out chromosome validation sets; the best performing model per species (selected by maximum validation auROC) is provided on the web toolkit and shown as a black dot with its value printed at the right. The mean and standard deviation across all folds are displayed as an open circle and a horizontal line, respectively (n = number of chromosome folds per species, see Supplementary Table 1).\n"
                f"d, Cross-species comparison of single-species deepCRE models, showing the mean predicted probability (PP) of high expression for genes labelled as lowly expressed (left) and for genes labelled as highly expressed (right). A mean PP near 0 indicates a gene is predicted to be lowly expressed, a PP near 1 highly expressed. The black diagonal line marks models predicting their own native species (including training data). Dashed rectangles emphasise cross-prediction within phylogenetically closely related groups. In the left heatmap, blue cells indicate that more than 50% of lowly expressed genes were predicted correctly, while orange cells show the opposite.\n"
                f"e, UPGMA tree of Euclidean distances between the 52-value cross-prediction profiles of the 13 species (rows and columns of d). Mantel test against the topological distances of the cladogram in a: r = {t_c['mantel_r']:.3f}, p = {t_c['mantel_p']:.4f} (Pearson, 9,999 permutations, 78 species pairs) shows that relatedness correlates with cross-species prediction performance. Leaf order: {', '.join(t_c['order'])}.\n"
                f"f, g, Multi-species reference training performance on the hold-out species using auROC (f) and accuracy (g), for models trained on all 13 species (black dot with value at right), in comparison to the single-species deepCRE model of that species provided with the toolkit (open diamond) and to the single-species models of the other 12 species (light grey dots).\n"
                f"h-l, Eudicot-only multi-species models (rows: held-out eudicot): cladogram of the 9 eudicot species (h); hold-out auROC (i) and accuracy (j) as in f and g; cross-species comparison as in d for eudicot-only models evaluated within the eudicots (k) and across clades on the grasses (l).\n"
                f"m-q, Grass-only multi-species models (rows: held-out grass), shown as in h-l: cladogram of the 4 grass species (m); hold-out auROC (n) and accuracy (o); cross-species comparison within the grasses (p) and across clades on the eudicots (q).\n"
            )

        # -------------------------------------------------------------------
        # 2) Figure 2 with MSR heatmaps (panels a-h)
        #    For leaf: this IS the main Figure 2 (Figure2_leaf / Figure2_with_MSR_heatmaps_leaf)
        #    For root: this is the Root Supplementary Figure (SuppFig_root_SSR_and_allMSR)
        # -------------------------------------------------------------------
        vals_fig2_ah = []
        ssr_block_ah = dict(
            ssr_block,
            panel_auroc=lambda ax: ssr_dots(ax, pf, tissue, ORDER, "auROC_v8", "auROC", vals_fig2_ah, "b"),
            panel_acc=lambda ax: ssr_dots(ax, pf, tissue, ORDER, "accuracy_v8", "Accuracy", vals_fig2_ah, "c"),
        )
        msr_block_ah = dict(
            name="Multi-species models (row: held-out species)",
            topo=TOPO_ALL,
            rows=ORDER,
            letters=["", "f", "g", "h"],
            heat_letter_mode="first_only",
            show_heats=True,
            m0=M[("all", 0)],
            m1=M[("all", 1)],
            panel_auroc=lambda ax: msr_dots(ax, summary, oth_rows, ORDER, "all", None, "auROC", "auROC", vals_fig2_ah, "f"),
            panel_acc=lambda ax: msr_dots(ax, summary, oth_rows, ORDER, "all", None, "accuracy", "Accuracy", vals_fig2_ah, "g"),
        )
        t_info_ah = build_figure(
            [ssr_block_ah, msr_block_ah],
            OUT / f"Figure2_with_MSR_heatmaps_{tissue}",
            L_SSR + [(None, None)] + L_MSR,
            legend_block=1,
            legend_x_mm=X_CB,
        )
        write_tsv(OUT / f"Figure2_with_MSR_heatmaps_{tissue}_values.tsv", vals_fig2_ah)
        with open(OUT / f"Figure2_with_MSR_heatmaps_{tissue}_legend.txt", "w", encoding="utf-8") as fh:
            fh.write(
                f"Figure ({tissue}, panels a-h): Single-species and leave-one-species-out multi-species model performance and cross-species prediction heatmaps across the 13 plant species (9 eudicots, 4 grasses).\n"
                f"a, Cladogram of the 13 species.\n"
                f"b, c, Single-species chromosome-fold validation auROC (b) and accuracy (c) (grey dots = individual chromosome folds; black dot = deployed web-tool model with value printed at right; open circle and horizontal bar = mean +/- s.d. across folds).\n"
                f"d, Single-species cross-species prediction heatmaps showing the mean predicted probability of high expression for true lowly expressed genes (left) and true highly expressed genes (right) across all 13 x 13 model-target species pairs.\n"
                f"e, UPGMA tree of Euclidean distances between the 52-value cross-prediction profiles of the 13 species (rows and columns of d). Mantel test against the topological distances of the cladogram in a: r = {t_info_ah['mantel_r']:.3f}, p = {t_info_ah['mantel_p']:.4f} (Pearson, 9,999 permutations, 78 species pairs). Leaf order: {', '.join(t_info_ah['order'])}.\n"
                f"f, g, Leave-one-species-out multi-species model auROC (f) and accuracy (g) on each held-out species (black dot with value at right), compared against the held-out species' own single-species web-tool model (open diamond) and single-species models transferred from the other 12 species (light grey dots).\n"
                f"h, Multi-species prediction heatmaps showing the mean predicted probability of high expression for true lowly expressed genes (left) and true highly expressed genes (right) across all 13 leave-one-species-out models (rows) and 13 target species (columns; grey diagonal = held-out species).\n"
            )
        if tissue == "root":
            # Also save explicitly as SuppFig_root_SSR_and_allMSR for direct supplementary use
            build_figure(
                [ssr_block_ah, msr_block_ah],
                OUT / "SuppFig_root_SSR_and_allMSR",
                L_SSR + [(None, None)] + L_MSR,
                legend_block=1,
                legend_x_mm=X_CB,
            )
            write_tsv(OUT / "SuppFig_root_SSR_and_allMSR_values.tsv", vals_fig2_ah)
            with open(OUT / "SuppFig_root_SSR_and_allMSR_legend.txt", "w", encoding="utf-8") as fh:
                fh.write(open(OUT / "Figure2_with_MSR_heatmaps_root_legend.txt", "r", encoding="utf-8").read())

        # -------------------------------------------------------------------
        # 3) Clade-only Supplementary MSR Figure (Eudicot-only & Grass-only only,
        #    with All-species MSR removed since it is already in Figure 2 / Root SuppFig):
        #    SuppFig_MSR_clades_only_{leaf,root} (panels a-j)
        # -------------------------------------------------------------------
        vals_clades = []
        D, G = ORDER[:9], ORDER[9:]
        s_dic_only = dict(
            name="Eudicot-only models (row: held-out eudicot)",
            topo=TOPO_DIC,
            rows=D,
            letters=["a", "b", "c", "d", "e"],
            heat_letter_mode="pair_first",
            show_heats=True,
            heats=[
                (M[("dicots", 0)], D, "True low", "on eudicots (own clade)"),
                (M[("dicots", 1)], D, "True high", None),
                (M[("dicots", 0)], G, "True low", "on grasses"),
                (M[("dicots", 1)], G, "True high", None),
            ],
            panel_auroc=lambda ax: msr_dots(ax, summary, oth_rows, D, "dicots", "monocots", "auROC", "auROC", vals_clades, "b"),
            panel_acc=lambda ax: msr_dots(ax, summary, oth_rows, D, "dicots", "monocots", "accuracy", "Accuracy", vals_clades, "c"),
        )
        s_gra_only = dict(
            name="Grass-only models (row: held-out grass)",
            topo=TOPO_GRA,
            rows=G,
            letters=["f", "g", "h", "i", "j"],
            heat_letter_mode="pair_first",
            show_heats=True,
            heats=[
                (M[("monocots", 0)], G, "True low", "on grasses (own clade)"),
                (M[("monocots", 1)], G, "True high", None),
                (M[("monocots", 0)], D, "True low", "on eudicots"),
                (M[("monocots", 1)], D, "True high", None),
            ],
            panel_auroc=lambda ax: msr_dots(ax, summary, oth_rows, G, "monocots", "dicots", "auROC", "auROC", vals_clades, "g"),
            panel_acc=lambda ax: msr_dots(ax, summary, oth_rows, G, "monocots", "dicots", "accuracy", "Accuracy", vals_clades, "h"),
        )
        build_figure(
            [s_dic_only, s_gra_only],
            OUT / f"SuppFig_MSR_clades_only_{tissue}",
            L_MSR + L_CROSS,
            legend_block=None,
        )
        write_tsv(OUT / f"SuppFig_MSR_clades_only_{tissue}_values.tsv", vals_clades)
        with open(OUT / f"SuppFig_MSR_clades_only_{tissue}_legend.txt", "w", encoding="utf-8") as fh:
            fh.write(
                f"Supplementary Figure ({tissue}, clade-specific multi-species models; All-species MSR block removed to avoid duplicating Figure 2): "
                f"Leave-one-species-out multi-species generalization within eudicots (9 species; a-e) and within grasses (4 species: O. sativa, H. vulgare, S. bicolor, Z. mays; f-j). "
                f"a, Cladogram of the 9 eudicot species. b, c, Eudicot-only leave-one-species-out auROC (b) and accuracy (c) on each held-out eudicot (black dot with value), compared with the median of the 4 grass-only models transferred across clades (grey square), the species' own single-species web-tool model (open diamond), and single-species models of the other 12 species (light grey dots). "
                f"d, e, Mean predicted probability of high expression for eudicot-only models evaluated within eudicots (d) and across clades on grasses (e). "
                f"f-j, Corresponding cladogram (f), held-out auROC (g), held-out accuracy (h), within-grass heatmaps (i), and cross-clade on-eudicot heatmaps (j) for grass-only leave-one-species-out models.\n"
            )

        # -------------------------------------------------------------------
        # 4) Full 3-block Supplementary MSR Figure (All-species + Eudicot-only + Grass-only)
        #    (useful for root if all root MSR blocks are shown together)
        # -------------------------------------------------------------------
        vals_supp = []
        s_all = dict(
            name="All-species multi-species models (row: held-out species)",
            topo=TOPO_ALL,
            rows=ORDER,
            letters=["a", "b", "c", "d"],
            heat_letter_mode="first_only",
            show_heats=True,
            m0=M[("all", 0)],
            m1=M[("all", 1)],
            panel_auroc=lambda ax: msr_dots(ax, summary, oth_rows, ORDER, "all", None, "auROC", "auROC", vals_supp, "b"),
            panel_acc=lambda ax: msr_dots(ax, summary, oth_rows, ORDER, "all", None, "accuracy", "Accuracy", vals_supp, "c"),
        )
        s_dic = dict(
            name="Eudicot-only models (row: held-out eudicot)",
            topo=TOPO_DIC,
            rows=D,
            letters=["e", "f", "g", "h", "i"],
            heat_letter_mode="pair_first",
            show_heats=True,
            heats=[
                (M[("dicots", 0)], D, "True low", "on eudicots (own clade)"),
                (M[("dicots", 1)], D, "True high", None),
                (M[("dicots", 0)], G, "True low", "on grasses"),
                (M[("dicots", 1)], G, "True high", None),
            ],
            panel_auroc=lambda ax: msr_dots(ax, summary, oth_rows, D, "dicots", "monocots", "auROC", "auROC", vals_supp, "f"),
            panel_acc=lambda ax: msr_dots(ax, summary, oth_rows, D, "dicots", "monocots", "accuracy", "Accuracy", vals_supp, "g"),
        )
        s_gra = dict(
            name="Grass-only models (row: held-out grass)",
            topo=TOPO_GRA,
            rows=G,
            letters=["j", "k", "l", "m", "n"],
            heat_letter_mode="pair_first",
            show_heats=True,
            heats=[
                (M[("monocots", 0)], G, "True low", "on grasses (own clade)"),
                (M[("monocots", 1)], G, "True high", None),
                (M[("monocots", 0)], D, "True low", "on eudicots"),
                (M[("monocots", 1)], D, "True high", None),
            ],
            panel_auroc=lambda ax: msr_dots(ax, summary, oth_rows, G, "monocots", "dicots", "auROC", "auROC", vals_supp, "k"),
            panel_acc=lambda ax: msr_dots(ax, summary, oth_rows, G, "monocots", "dicots", "accuracy", "Accuracy", vals_supp, "l"),
        )
        build_figure(
            [s_all, s_dic, s_gra],
            OUT / f"SuppFig_MSR_{tissue}",
            L_MSR + L_CROSS,
            legend_block=None,
        )
        write_tsv(OUT / f"SuppFig_MSR_{tissue}_values.tsv", vals_supp)
        with open(OUT / f"SuppFig_MSR_{tissue}_legend.txt", "w", encoding="utf-8") as fh:
            fh.write(
                f"Supplementary Figure ({tissue}): Multi-species (MSR) leave-one-species-out generalization across all 13 species (a-d), within eudicots (9 species; e-i), and within grasses (4 species: O. sativa, H. vulgare, S. bicolor, Z. mays; j-n).\n"
            )


        # -------------------------------------------------------------------
        # 5) Combined Leaf + Root Clade-Only Supplementary Figure (panels a-t)
        # -------------------------------------------------------------------
    # After both tissues are processed, build a combined Leaf + Root clade-only MSR figure
    full_leaf = read_tsv_dicts(MSR_DIR / "msr_eval_leaf_full.tsv")
    sum_leaf, oth_leaf = compute_msr_summary(full_leaf, pf_all, "leaf")
    full_root = read_tsv_dicts(MSR_DIR / "msr_eval_root_full.tsv")
    sum_root, oth_root = compute_msr_summary(full_root, pf_all, "root")
    M_both = {}
    for tis in ["leaf", "root"]:
        for s in ["dicots", "monocots"]:
            for k in (0, 1):
                M_both[(tis, s, k)], _ = read_matrix_tsv(MSR_DIR / f"suppfig_msr_mean_{tis}_{s}_class{k}.tsv")
    D, G = ORDER[:9], ORDER[9:]
    vals_comb = []
    b_leaf_dic = dict(
        name="Leaf — Eudicot-only models (row: held-out eudicot)",
        topo=TOPO_DIC,
        rows=D,
        letters=["a", "b", "c", "d", "e"],
        heat_letter_mode="pair_first",
        show_heats=True,
        heats=[
            (M_both[("leaf", "dicots", 0)], D, "True low", "on eudicots (own clade)"),
            (M_both[("leaf", "dicots", 1)], D, "True high", None),
            (M_both[("leaf", "dicots", 0)], G, "True low", "on grasses"),
            (M_both[("leaf", "dicots", 1)], G, "True high", None),
        ],
        panel_auroc=lambda ax: msr_dots(ax, sum_leaf, oth_leaf, D, "dicots", "monocots", "auROC", "auROC", vals_comb, "b"),
        panel_acc=lambda ax: msr_dots(ax, sum_leaf, oth_leaf, D, "dicots", "monocots", "accuracy", "Accuracy", vals_comb, "c"),
    )
    b_leaf_gra = dict(
        name="Leaf — Grass-only models (row: held-out grass)",
        topo=TOPO_GRA,
        rows=G,
        letters=["f", "g", "h", "i", "j"],
        heat_letter_mode="pair_first",
        show_heats=True,
        heats=[
            (M_both[("leaf", "monocots", 0)], G, "True low", "on grasses (own clade)"),
            (M_both[("leaf", "monocots", 1)], G, "True high", None),
            (M_both[("leaf", "monocots", 0)], D, "True low", "on eudicots"),
            (M_both[("leaf", "monocots", 1)], D, "True high", None),
        ],
        panel_auroc=lambda ax: msr_dots(ax, sum_leaf, oth_leaf, G, "monocots", "dicots", "auROC", "auROC", vals_comb, "g"),
        panel_acc=lambda ax: msr_dots(ax, sum_leaf, oth_leaf, G, "monocots", "dicots", "accuracy", "Accuracy", vals_comb, "h"),
    )
    b_root_dic = dict(
        name="Root — Eudicot-only models (row: held-out eudicot)",
        topo=TOPO_DIC,
        rows=D,
        letters=["k", "l", "m", "n", "o"],
        heat_letter_mode="pair_first",
        show_heats=True,
        heats=[
            (M_both[("root", "dicots", 0)], D, "True low", "on eudicots (own clade)"),
            (M_both[("root", "dicots", 1)], D, "True high", None),
            (M_both[("root", "dicots", 0)], G, "True low", "on grasses"),
            (M_both[("root", "dicots", 1)], G, "True high", None),
        ],
        panel_auroc=lambda ax: msr_dots(ax, sum_root, oth_root, D, "dicots", "monocots", "auROC", "auROC", vals_comb, "l"),
        panel_acc=lambda ax: msr_dots(ax, sum_root, oth_root, D, "dicots", "monocots", "accuracy", "Accuracy", vals_comb, "m"),
    )
    b_root_gra = dict(
        name="Root — Grass-only models (row: held-out grass)",
        topo=TOPO_GRA,
        rows=G,
        letters=["p", "q", "r", "s", "t"],
        heat_letter_mode="pair_first",
        show_heats=True,
        heats=[
            (M_both[("root", "monocots", 0)], G, "True low", "on grasses (own clade)"),
            (M_both[("root", "monocots", 1)], G, "True high", None),
            (M_both[("root", "monocots", 0)], D, "True low", "on eudicots"),
            (M_both[("root", "monocots", 1)], D, "True high", None),
        ],
        panel_auroc=lambda ax: msr_dots(ax, sum_root, oth_root, G, "monocots", "dicots", "auROC", "auROC", vals_comb, "q"),
        panel_acc=lambda ax: msr_dots(ax, sum_root, oth_root, G, "monocots", "dicots", "accuracy", "Accuracy", vals_comb, "r"),
    )
    build_figure(
        [b_leaf_dic, b_leaf_gra, b_root_dic, b_root_gra],
        OUT / "SuppFig_MSR_clades_only_leaf_and_root",
        L_MSR + L_CROSS,
        legend_block=None,
    )
    write_tsv(OUT / "SuppFig_MSR_clades_only_leaf_and_root_values.tsv", vals_comb)


if __name__ == "__main__":
    main()

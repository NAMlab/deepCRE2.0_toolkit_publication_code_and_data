"""
Purpose:
    Extend `train_models_seed_restart.1.py` with `--only A,B` filtering so interrupted multi-fold
    or multi-species configs can resume only specific validation entities (`<fold>` or `<species_name>`)
    using identical deterministic seeds (`--seed 2026`), plateau restart logic, and restart logging.

Manuscript (Peleke et al., 2026):
    Used to finish interrupted all-species MSR tasks (`job_finish_msr_all_dPmean.0.sh`), run all
    shuffled-sequence controls (`job_ssc_control_dPmean.0.sh`, `--max_attempts 1`), and train the
    *A. thaliana* TAIR10 exception (`job_train_atha_tair10.0.sh`).

Usage & I/O (run from `code/06_benchmarking/workflows/deepcre/src/deepCRE/` with TF 2.10 env):
    python train_models_seed_restart.2.py --input <config.json> \
        [--seed 2026] [--max_attempts 5] [--plateau_loss 0.68] [--first_attempt 0] [--only A,B]

How to change or extend:
    Pass `--only <item1,item2>` to target specific failed chromosomes or MSR hold-out species,
    or `--max_attempts 1` for single-shot controls (`model_case="ssc"`).
"""
import argparse
import os
import random
import sys

import numpy as np
import tensorflow as tf

import train_models as tm

ARGS = None
_original_deep_cre = tm.deep_cre
_original_train_deep_cre = tm.train_deep_cre


def seed_all(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    tf.random.set_seed(seed)


def train_deep_cre_seeded(*args, **kwargs):
    val = kwargs.get("val_chromosome")
    if ARGS.only and val not in ARGS.only:
        print(f"[seed_restart] skipping {val} (not in --only)")
        return (float("nan"),) * 4
    seed_all(ARGS.seed)  # fixes the random class balancing of this fold
    return _original_train_deep_cre(*args, **kwargs)


def deep_cre_with_restart(x_train, y_train, x_val, y_val, output_name, model_case, chrom, time_stamp, test=False):
    log = f"restart_log_{os.path.basename(str(output_name))}.tsv"
    new_log = not os.path.exists(log)
    output = None
    for attempt in range(ARGS.first_attempt, ARGS.first_attempt + ARGS.max_attempts):
        seed = ARGS.seed + attempt
        seed_all(seed)
        tf.keras.backend.clear_session()
        output = _original_deep_cre(x_train, y_train, x_val, y_val, output_name, model_case, chrom, time_stamp, test)
        loss = float(output[0])
        accepted = loss <= ARGS.plateau_loss
        with open(log, "a") as fh:
            if new_log:
                fh.write("output_name\tfold\tattempt\tseed\tval_loss\tval_accuracy\tval_auROC\tval_auPR\taccepted\n")
                new_log = False
            fh.write(f"{output_name}\t{chrom}\t{attempt + 1}\t{seed}\t" + "\t".join(f"{v:.6f}" for v in output) + f"\t{accepted}\n")
        print(f"[seed_restart] fold {chrom} attempt {attempt + 1} seed {seed}: val_loss {loss:.4f} -> {'accepted' if accepted else 'plateau, restarting'}")
        if accepted:
            return output
    print(f"[seed_restart] WARNING fold {chrom}: still on the plateau after {ARGS.max_attempts} attempts; last checkpoint kept")
    return output


def main():
    global ARGS
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", "-i", required=True)
    p.add_argument("--seed", type=int, default=2026)
    p.add_argument("--max_attempts", type=int, default=5)
    p.add_argument("--plateau_loss", type=float, default=0.68)
    p.add_argument("--first_attempt", type=int, default=0)
    p.add_argument("--only", type=lambda x: set(x.split(",")), default=None)
    ARGS = p.parse_args()
    tm.deep_cre = deep_cre_with_restart
    tm.train_deep_cre = train_deep_cre_seeded
    sys.argv = [sys.argv[0], "--input", ARGS.input]  # train_models parses --input itself
    tm.main()


if __name__ == "__main__":
    main()

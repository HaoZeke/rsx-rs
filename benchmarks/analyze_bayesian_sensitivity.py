#!/usr/bin/env python3
"""
Prior sensitivity analysis for the Bayesian triage layer.

Sweeps π (baseline sex-linkage probability) and p_sex on the literature
candidate tables and reports how the number of posterior > 0.9 markers
and the final biological calls change.
"""

import argparse
from pathlib import Path
import pandas as pd
import sys

# Reuse the exact posterior function from the existing analysis
sys.path.insert(0, str(Path(__file__).parent))
from analyze_bayesian_evidence import posterior_sex_linked

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--triage-csv", type=Path, default=Path("results/literature_candidate_triage.csv"))
    parser.add_argument("--out", type=Path, default=Path("results/prior_sensitivity.csv"))
    args = parser.parse_args()

    df = pd.read_csv(args.triage_csv)

    # We need raw counts. The triage CSV has Group1_Present, Group2_Present, Group1_Total, Group2_Total
    # Filter to the rows that were actually evaluated (they have the counts)
    required_cols = ["dataset", "min_depth", "group1_present", "group2_present", "group1_total", "group2_total"]
    if not all(c in df.columns for c in required_cols):
        print("Triage CSV does not have the required count columns.")
        print("Available:", list(df.columns)[:10])
        return

    # Parameters to sweep (common values used in literature)
    pis = [0.001, 0.005, 0.01, 0.02, 0.05, 0.1]
    p_sexes = [0.8, 0.85, 0.9, 0.95]

    records = []

    for (dataset, depth), group in df.groupby(["dataset", "min_depth"]):
        for pi in pis:
            for p_sex in p_sexes:
                n_post = 0
                for _, row in group.iterrows():
                    g1 = int(row["group1_present"])
                    g2 = int(row["group2_present"])
                    t1 = int(row["group1_total"])
                    t2 = int(row["group2_total"])
                    post = posterior_sex_linked(g1, g2, t1, t2, pi, p_sex)
                    if post > 0.9:
                        n_post += 1

                records.append({
                    "dataset": dataset,
                    "min_depth": depth,
                    "pi": pi,
                    "p_sex": p_sex,
                    "n_posterior_gt_0_9": n_post,
                })

    out_df = pd.DataFrame(records)
    out_df.to_csv(args.out, index=False)
    print(f"Wrote {len(out_df)} rows to {args.out}")
    print(out_df.head(20))

if __name__ == "__main__":
    main()

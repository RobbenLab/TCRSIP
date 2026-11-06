#!/usr/bin/env python3

# Computes Accuracy/ROC_AUC/PR_AUC/PPV/NPV for the DL pipeline's FullRun comparison.
# predictions.csv already carries both the true label ("response") and a raw binder
# probability ("Prob") per sample, plus a Set column (train/val) - it just hasn't been
# reduced to a summary yet, unlike the DL pipeline's older 5-fold CV runs (which already
# write their own results_summary.csv). Only the held-out "val" rows are scored.

import os

import pandas as pd

from ScorePredictions import compute_binary_metrics

PREDICTIONS_CSV = "./data/DL/FullRun/predictions.csv"
OUT_FILE = "./data/DL/FullRun/fullrun_summary.csv"


def main():
    preds = pd.read_csv(PREDICTIONS_CSV)
    val = preds[preds["Set"] == "val"]

    metrics = compute_binary_metrics(val["response"], val["Prob"])
    metrics["Model"] = "ResNet"
    metrics["Comparison"] = "core_immrep_split"

    summary = pd.DataFrame([metrics])[["Model", "Comparison", "n", "n_positive", "Accuracy", "ROC_AUC", "PR_AUC", "PPV", "NPV"]]
    summary.to_csv(OUT_FILE, index=False)
    print(summary.to_string(index=False))
    print(f"Saved to {OUT_FILE}")


if __name__ == "__main__":
    main()

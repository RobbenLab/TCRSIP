#!/usr/bin/env python3

# Computes Accuracy/ROC_AUC/PR_AUC/PPV/NPV for ERGO-II's two model variants (Auto,
# LSTM) against the full immrep test set. ERGO's own output files
# (ImmpredAutoresults.csv/ImmpredLMSTresults.csv) only carry a raw binder probability
# ("Score") per sample - no ground-truth label at all - but they're in the exact same
# row order as the immrep rows of Annotation.csv (confirmed: the Peptide column matches
# 1:1 in order), so the true binder label for each row is attached positionally rather
# than via any shared id column.

import os

import pandas as pd

from ScorePredictions import compute_binary_metrics

ANNOTATION_CSV = "./data/ML/Annotation.csv"
ERGO_DIR = "./data/ERGO"
ERGO_RESULTS = {
    "ERGO-Auto": "ImmpredAutoresults.csv",
    "ERGO-LSTM": "ImmpredLMSTresults.csv",
}
OUT_FILE = os.path.join(ERGO_DIR, "ergo_summary.csv")


def main():
    annotation = pd.read_csv(ANNOTATION_CSV)
    immrep = annotation[annotation["source"] == "immrep"].reset_index(drop=True)

    rows = []
    for model_name, filename in ERGO_RESULTS.items():
        preds = pd.read_csv(os.path.join(ERGO_DIR, filename))
        if len(preds) != len(immrep) or not (preds["Peptide"].values == immrep["Peptide"].values).all():
            raise ValueError(
                f"{filename} no longer lines up row-for-row with Annotation.csv's immrep rows - "
                "can't attach the true binder label this way anymore."
            )
        metrics = compute_binary_metrics(immrep["response"], preds["Score"])
        metrics["Model"] = model_name
        # Scored against the complete immrep test set (not just the structured subset
        # ./dat covers) - unlike our own pipelines, ERGO's own training data is external
        # and not one of our defined train/test splits, so "Comparison" here only
        # describes what it was evaluated against.
        metrics["Comparison"] = "immrep_full"
        rows.append(metrics)

    summary = pd.DataFrame(rows)[["Model", "Comparison", "n", "n_positive", "Accuracy", "ROC_AUC", "PR_AUC", "PPV", "NPV"]]
    summary.to_csv(OUT_FILE, index=False)
    print(summary.to_string(index=False))
    print(f"Saved to {OUT_FILE}")


if __name__ == "__main__":
    main()

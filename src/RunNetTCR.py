#!/usr/bin/env python3

# Runs NetTCR-2.0 (chain=ab) on each comparison directory produced by PrepNetTCR.py,
# and scores the resulting predictions against the true binder labels.

import os
import subprocess
import sys
from datetime import datetime

import pandas as pd

from ScorePredictions import compute_binary_metrics

COMPARISONS_DIR = "./data/NetTCR/Comparisons"
NETTCR_DIR = "./data/NetTCR/NetTCR-2.0"
EPOCHS = 100
CHAIN = "ab"


def run_nettcr(comparison_dir):
    train_file = os.path.abspath(os.path.join(comparison_dir, "train.csv"))
    test_file = os.path.abspath(os.path.join(comparison_dir, "test.csv"))
    out_file = os.path.abspath(os.path.join(comparison_dir, "predictions.csv"))

    print(f"Running NetTCR on {comparison_dir} ...")
    subprocess.run(
        [
            sys.executable, "nettcr.py",
            "--trainfile", train_file,
            "--testfile", test_file,
            "--chain", CHAIN,
            "--epochs", str(EPOCHS),
            "--outfile", out_file,
        ],
        cwd=NETTCR_DIR,
        check=True,
    )
    return out_file


def score(comparison_dir, out_file):
    preds = pd.read_csv(out_file)
    train = pd.read_csv(os.path.join(comparison_dir, "train.csv"))

    metrics = compute_binary_metrics(preds["binder"], preds["prediction"])
    n_test, n_test_positive = metrics.pop("n"), metrics.pop("n_positive")
    return {
        "n_train": len(train),
        "n_train_positive": int(train["binder"].sum()),
        "n_test": n_test,
        "n_test_positive": n_test_positive,
        **metrics,
    }


def write_report(summary, report_file):
    # A plain-text report meant to be opened and read directly, as opposed to the CSV
    # (meant for further processing) - one block per comparison, ranked by ROC_AUC so
    # the best-generalizing comparison is immediately visible at the top.
    ranked = summary.sort_values("ROC_AUC", ascending=False)
    lines = [
        "NetTCR-2.0 comparison summary",
        f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        "=" * 60,
        "",
    ]
    for _, row in ranked.iterrows():
        lines.append(f"{row['Comparison']}")
        lines.append(f"  Train: {row['n_train']} samples ({row['n_train_positive']} positive)")
        lines.append(f"  Test:  {row['n_test']} samples ({row['n_test_positive']} positive)")
        lines.append(
            f"  Accuracy={row['Accuracy']:.3f}  ROC_AUC={row['ROC_AUC']:.3f}  PR_AUC={row['PR_AUC']:.3f}"
            f"  PPV={row['PPV']:.3f}  NPV={row['NPV']:.3f}"
        )
        lines.append("")

    with open(report_file, "w") as f:
        f.write("\n".join(lines))


def main(skip_existing=True):
    comparisons = sorted(
        d for d in os.listdir(COMPARISONS_DIR)
        if os.path.isdir(os.path.join(COMPARISONS_DIR, d))
    )
    results = []
    for name in comparisons:
        comparison_dir = os.path.join(COMPARISONS_DIR, name)
        out_file = os.path.join(comparison_dir, "predictions.csv")
        if skip_existing and os.path.exists(out_file):
            print(f"Skipping {name} - {out_file} already exists")
        else:
            out_file = run_nettcr(comparison_dir)
        metrics = score(comparison_dir, out_file)
        metrics["Comparison"] = name
        results.append(metrics)
        print(f"  {name}: {metrics}")

    columns = ["Comparison", "n_train", "n_train_positive", "n_test", "n_test_positive", "Accuracy", "ROC_AUC", "PR_AUC", "PPV", "NPV"]
    summary = pd.DataFrame(results)[columns]
    summary_file = os.path.join(COMPARISONS_DIR, "nettcr_summary.csv")
    summary.to_csv(summary_file, index=False)

    report_file = os.path.join(COMPARISONS_DIR, "nettcr_summary.txt")
    write_report(summary, report_file)

    print()
    print(summary.to_string(index=False))
    print(f"\nSaved summary to {summary_file} and {report_file}")


if __name__ == "__main__":
    main()

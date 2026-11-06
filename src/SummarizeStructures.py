#!/usr/bin/env python3

# End-of-run summary over a set of ./dat/<name>/ structures: how many, their MHC class
# mix, how many look misfolded, and - reading whatever Inference.<Model> keys InferML.py/
# InferDL.py have written into each structure's json - each model's binder calls plus an
# ensemble across every model present (mean probability). Also writes results.csv, one row
# per structure: name, peptide, CDR3a/b, MHC class, misfold score, and every model's/the
# ensemble's Prob/Pred. Pure read of already-written per-structure json, no recomputation.
#
#   python src/SummarizeStructures.py --dat-dir ./dat --names-file ./dat/new_structures.txt

import argparse
import json
import os

import pandas as pd

from TCRclassifier import getMisfolded

DAT_DIR = "./dat"
MISFOLD_THRESHOLD = 0.8
THRESHOLD = 0.5


def cdr3_seq(js, chain, region_key):
    """Slice the chain's full sequence by its CDR3 region - an inclusive [start, end]
    residue index pair into Sequence.<chain>.Seq, as SetupData.py stores it."""
    try:
        seq = js["Sequence"][chain]["Seq"]
        start, end = js["Region"][chain][region_key]
        return seq[start:end + 1]
    except (KeyError, TypeError, IndexError):
        return None


def collect_rows(dat_dir, subdirs, misfolded, threshold):
    rows = []
    n_missing_json = 0
    for name in subdirs:
        json_path = os.path.join(dat_dir, name, f"{name}.json")
        if not os.path.isfile(json_path):
            n_missing_json += 1
            continue
        with open(json_path) as f:
            js = json.load(f)

        row = {
            "name": name,
            "Peptide": js.get("Sequence", {}).get("Peptide", {}).get("Seq"),
            "CDR3a": cdr3_seq(js, "TCRa", "CDR3a"),
            "CDR3b": cdr3_seq(js, "TCRb", "CDR3b"),
            "MHC_class": js.get("Features", {}).get("Sequence", {}).get("MHC_class"),
            "Misfold_Score": js.get("Scoring", {}).get("Final"),
            "Misfolded": int(name in misfolded),
        }

        probs = []
        for model_name, result in js.get("Inference", {}).items():
            prob = result.get("Prob")
            row[f"{model_name}_Prob"] = prob
            row[f"{model_name}_Pred"] = result.get("Pred")
            if prob is not None:
                probs.append(prob)
        if probs:
            ens_prob = sum(probs) / len(probs)
            row["Ensemble_Prob"] = ens_prob
            row["Ensemble_Pred"] = int(ens_prob >= threshold)

        rows.append(row)
    return rows, n_missing_json


def print_summary(dat_dir, rows, total, n_missing_json, restricted):
    mhc_class_counts = {}
    inference_by_model = {}
    ensemble_binders = 0
    n_with_inference = 0
    n_misfolded = 0

    for row in rows:
        mhc_class_counts[row["MHC_class"]] = mhc_class_counts.get(row["MHC_class"], 0) + 1
        n_misfolded += row["Misfolded"]
        row_had_inference = False
        for key, value in row.items():
            if not key.endswith("_Prob") or key == "Ensemble_Prob":
                continue
            model_name = key[:-len("_Prob")]
            if value is None:
                continue
            bucket = inference_by_model.setdefault(model_name, {"n": 0, "n_binder": 0})
            bucket["n"] += 1
            bucket["n_binder"] += int(row.get(f"{model_name}_Pred") == 1)
            row_had_inference = True
        if row_had_inference:
            n_with_inference += 1
            ensemble_binders += int(row.get("Ensemble_Pred") == 1)

    print("=" * 60)
    scope = f"{total} structure(s) in {dat_dir}" + (" (restricted set)" if restricted else "")
    print(f"Run summary: {scope}")
    if n_missing_json:
        print(f"  {n_missing_json} name(s) had no {dat_dir}/<name>/<name>.json")

    if total:
        print("\nMHC class mix:")
        for cls, count in sorted(mhc_class_counts.items(), key=lambda kv: (kv[0] is None, kv[0])):
            label = f"Class {cls}" if cls is not None else "Unknown"
            print(f"  {label}: {count} ({100 * count / total:.1f}%)")

        print(f"\nPotential misfolds (Scoring.Final < {MISFOLD_THRESHOLD}): "
             f"{n_misfolded} ({100 * n_misfolded / total:.1f}%)")

    if inference_by_model:
        print("\nBinding model predictions:")
        for model_name, bucket in sorted(inference_by_model.items()):
            n, n_binder = bucket["n"], bucket["n_binder"]
            print(f"  {model_name}: {n_binder}/{n} predicted binder ({100 * n_binder / n:.1f}%)")
        print(f"  Ensemble (mean of available model probabilities): "
             f"{ensemble_binders}/{n_with_inference} predicted binder "
             f"({100 * ensemble_binders / n_with_inference:.1f}%)")
    else:
        print("\nNo model predictions recorded yet (run with --ml/--dl first).")
    print("=" * 60)


def summarize(dat_dir, names=None, misfold_threshold=MISFOLD_THRESHOLD, threshold=THRESHOLD, output=None):
    subdirs = sorted(d for d in os.listdir(dat_dir) if os.path.isdir(os.path.join(dat_dir, d)))
    if names is not None:
        subdirs = [d for d in subdirs if d in names]

    misfolded = set(getMisfolded(dat_dir, threshold=misfold_threshold))
    rows, n_missing_json = collect_rows(dat_dir, subdirs, misfolded, threshold)
    print_summary(dat_dir, rows, len(subdirs), n_missing_json, names is not None)

    output = output or os.path.join(dat_dir, "results.csv")
    pd.DataFrame(rows).to_csv(output, index=False)
    print(f"Wrote {len(rows)} structures' results to {output}")


def main():
    ap = argparse.ArgumentParser(description="Summarize a set of processed structures and their model predictions")
    ap.add_argument("--dat-dir", default=DAT_DIR, help="./dat/<name>/<name>.json folders")
    ap.add_argument("--names-file", default=None,
                    help="Optional text file of structure names (one per line) to restrict the "
                         "summary to, instead of every structure in --dat-dir")
    ap.add_argument("--misfold-threshold", type=float, default=MISFOLD_THRESHOLD)
    ap.add_argument("--threshold", type=float, default=THRESHOLD, help="probability cutoff for the ensemble call")
    ap.add_argument("-o", "--output", default=None, help="results CSV path (default: <dat-dir>/results.csv)")
    args = ap.parse_args()

    names = None
    if args.names_file:
        with open(args.names_file) as f:
            names = {line.strip() for line in f if line.strip()}

    summarize(args.dat_dir, names, args.misfold_threshold, args.threshold, args.output)


if __name__ == "__main__":
    main()

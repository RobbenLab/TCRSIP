#!/usr/bin/env python3

# Reuses the features ExtractFeaturestoJSON.py already computed and stored in each
# <Dir>/<name>/<name>.json's "Features" key, flattening them into the same CSV format
# ExtractFeaturestoJSON.py writes - without recomputing anything (no biopython/pyrosetta
# calls). For when --Dir already has every structure's features extracted and only the
# combined CSV needs rebuilding, e.g. Run.sh's feature-extraction bypass when no raw
# structures directory is given.
#
#   python src/CompileFeaturesCSV.py --Dir ./dat --Out ./dat/new_features.csv

import argparse
import json
import os

import flatdict
import pandas as pd


def uniformDict(d):
    return dict(flatdict.FlatDict(d))


def compile_features(dat_dir, names=None):
    subdirs = sorted(d for d in os.listdir(dat_dir) if os.path.isdir(os.path.join(dat_dir, d)))
    if names is not None:
        subdirs = [d for d in subdirs if d in names]

    rows, missing = [], []
    for name in subdirs:
        json_path = os.path.join(dat_dir, name, f"{name}.json")
        if not os.path.isfile(json_path):
            missing.append(name)
            continue
        with open(json_path) as f:
            js = json.load(f)
        feat = js.get("Features")
        if not feat:
            missing.append(name)
            continue
        row = {"name": name}
        row.update(uniformDict(feat.get("Sequence", {})))
        row.update(uniformDict(feat.get("Structure", {})))
        row.update(uniformDict(feat.get("Physics", {})))
        rows.append(row)

    if missing:
        print(f"Note: {len(missing)} structure(s) in {dat_dir} have no extracted Features yet, "
             f"skipped: {missing[:10]}{'...' if len(missing) > 10 else ''}")
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(
        description="Compile already-extracted per-structure Features into one CSV, without recomputing them"
    )
    ap.add_argument("-d", "--Dir", required=True, help="dat/ directory of <name>/<name>.json folders")
    ap.add_argument("-o", "--Out", required=True, help="output features CSV path")
    ap.add_argument("--names-file", default=None,
                    help="Optional text file of structure names (one per line) to restrict to, "
                         "instead of every subdirectory in --Dir")
    args = ap.parse_args()

    names = None
    if args.names_file:
        with open(args.names_file) as f:
            names = {line.strip() for line in f if line.strip()}

    df = compile_features(args.Dir, names)
    df.to_csv(args.Out, index=False)
    print(f"Compiled {len(df)} structures' already-extracted features -> {args.Out}")


if __name__ == "__main__":
    main()

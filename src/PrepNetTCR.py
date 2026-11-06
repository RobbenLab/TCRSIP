#!/usr/bin/env python3

# Builds CDR3a/CDR3b/peptide/binder train-test CSVs for NetTCR-2.0 from this project's
# Annotation.csv plus the per-structure JSONs in ./dat/, for a fixed set of model
# comparisons:
#   core_80_20         : core (ground+exp+fake) split 80/20, stratified by binder
#   core_vs_immrep     : train = all of core, test = all of immrep
#   core_vs_scan       : train = all of core, test = all of scan
#   core_immrep_80_20  : core+immrep combined, split 80/20, stratified by binder
#   full_run           : core+immrep, train/val split fixed to match the DL pipeline's
#                        own FullRun comparison (data/DL/FullRun/split.csv) exactly
#   paper_train_vs_immrep_full : NetTCR-2.0's own published training set (not derived
#                        from our Annotation.csv/dat at all), tested on the COMPLETE
#                        immrep test set (not just the structured subset - see
#                        build_full_immrep_table)
# "core" = the samples with folded/annotated structures (ground, exp, fake sources);
# immrep has sequence-only entries (no folded structure) for the vast majority of rows,
# so only the subset of immrep with a structure in ./dat/ can be used for the other
# immrep comparisons above.

import os
import json
import pandas as pd
from sklearn.model_selection import train_test_split

ANNOTATION_CSV = "./data/ML/Annotation.csv"
STRUCTURE_DIR = "./dat"
STANDARD_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")
MAX_PEPTIDE_LEN = 9  # NetTCR-2.0 hard-codes a 9-residue peptide buffer (class I only);
# longer peptides (e.g. class II, 12-25 aa) crash its encoder outright rather than
# truncating, so they're excluded entirely here rather than silently cut down to 9.
MAX_CDR3_LEN = 30  # NetTCR-2.0 hard-codes a 30-residue CDR3 buffer too. A real CDR3 is
# essentially always well under that (our data's 75th percentile is ~14); anything
# longer here is an IMGT-numbering/Turn-annotation failure upstream (an entire chain
# mislabeled as "cdr3"), not a real CDR3 - drop it rather than feed NetTCR garbage.
OUTPUT_DIR = "./data/NetTCR/Comparisons"
RANDOM_STATE = 42
TEST_SIZE = 0.2
FULLRUN_SPLIT_CSV = "./data/DL/FullRun/split.csv"  # fixed train/val assignment to match the DL pipeline's FullRun comparison
NATIVE_TRAIN_CSV = "./data/NetTCR/NetTCR-2.0/data/train_ab_95_alphabeta.csv"  # NetTCR-2.0's own published training set
# ERGO-II's raw-prediction output on the immrep test set (see ScoreERGO.py) happens to
# also carry the original TRA/TRB/Peptide for every immrep row, in the exact same order
# as Annotation.csv's immrep rows (verified: Peptide matches 1:1 positionally) - reused
# here purely as a source of CDR3 sequences, so the COMPLETE immrep test set can be used
# (most immrep rows have no folded structure in ./dat, so extract_cdr3 alone only
# recovers a small subset of it).
ERGO_IMMREP_CSV = "./data/ERGO/ImmpredAutoresults.csv"


def extract_cdr3(name):
    # Pulls CDR3a/CDR3b from <STRUCTURE_DIR>/<name>/<name>.json: the Consensus "Turn"
    # array labels each residue of the structure-derived sequence (PDBseq) with its
    # IMGT region, so the CDR3 is simply every residue labeled "cdr3". Returns
    # (None, None) if the structure/JSON is missing or IMGT numbering failed for
    # either chain.
    json_path = os.path.join(STRUCTURE_DIR, name, f"{name}.json")
    if not os.path.exists(json_path):
        return None, None
    try:
        with open(json_path) as f:
            data = json.load(f)
        cdr3s = {}
        for chain in ("TCRa", "TCRb"):
            seq = data["Sequence"][chain]
            consensus = seq.get("Consensus")
            if not consensus:
                return None, None
            cdr3 = "".join(aa for aa, turn in zip(seq["PDBseq"], consensus["Turn"]) if turn == "cdr3")
            if not cdr3:
                return None, None
            cdr3s[chain] = cdr3
        return cdr3s["TCRa"], cdr3s["TCRb"]
    except Exception:
        return None, None


def is_standard_aa_sequence(seq):
    return set(seq) <= STANDARD_AMINO_ACIDS


def build_nettcr_table(annotation):
    # Adds CDR3a/CDR3b to annotation (dropping rows where extraction failed) and
    # renames to NetTCR's expected column names.
    annotation = annotation.copy()
    cdr3 = annotation["name"].apply(extract_cdr3)
    annotation["CDR3a"] = [c[0] for c in cdr3]
    annotation["CDR3b"] = [c[1] for c in cdr3]
    annotation = annotation.dropna(subset=["CDR3a", "CDR3b"])

    # NetTCR's BLOSUM50 encoder has no entry for non-standard residues (e.g. "X" for an
    # unresolved residue from the original structure) and aborts the whole run on one -
    # drop any row with one in CDR3a/CDR3b/peptide rather than let a handful of rows
    # crash the entire comparison.
    valid = (
        annotation["CDR3a"].apply(is_standard_aa_sequence)
        & annotation["CDR3b"].apply(is_standard_aa_sequence)
        & annotation["Peptide"].apply(is_standard_aa_sequence)
    )
    dropped = (~valid).sum()
    if dropped:
        print(f"Dropping {dropped} rows with a non-standard amino acid in CDR3a/CDR3b/peptide")
    annotation = annotation[valid]

    too_long_pep = annotation["Peptide"].str.len() > MAX_PEPTIDE_LEN
    if too_long_pep.sum():
        print(f"Dropping {too_long_pep.sum()} rows with a peptide longer than {MAX_PEPTIDE_LEN} residues (likely MHC class II)")
    annotation = annotation[~too_long_pep]

    too_long_cdr3 = (annotation["CDR3a"].str.len() > MAX_CDR3_LEN) | (annotation["CDR3b"].str.len() > MAX_CDR3_LEN)
    if too_long_cdr3.sum():
        print(f"Dropping {too_long_cdr3.sum()} rows with a CDR3a/CDR3b longer than {MAX_CDR3_LEN} residues (IMGT annotation failure)")
    annotation = annotation[~too_long_cdr3]

    table = annotation.rename(columns={"Peptide": "peptide", "response": "binder"})
    return table[["name", "CDR3a", "CDR3b", "peptide", "binder", "source"]]


def build_full_immrep_table(annotation):
    # The vast majority of immrep rows have no folded structure in ./dat, so
    # extract_cdr3 (which reads ./dat/<name>/<name>.json) only recovers CDR3 sequences
    # for a small subset of them. ERGO_IMMREP_CSV happens to carry the original
    # TRA/TRB/Peptide for every immrep row, in the exact same order as Annotation.csv's
    # immrep rows - reused here purely as a source of CDR3 sequences, so the complete
    # immrep test set can be used rather than just the structured subset.
    immrep = annotation[annotation["source"] == "immrep"].reset_index(drop=True)
    ergo = pd.read_csv(ERGO_IMMREP_CSV)
    if len(ergo) != len(immrep) or not (ergo["Peptide"].values == immrep["Peptide"].values).all():
        raise ValueError(
            f"{ERGO_IMMREP_CSV} no longer lines up row-for-row with Annotation.csv's immrep "
            "rows - can't recover CDR3 sequences for the full immrep test set this way anymore."
        )
    table = immrep[["name", "response"]].rename(columns={"response": "binder"})
    table["CDR3a"] = ergo["TRA"].values
    table["CDR3b"] = ergo["TRB"].values
    table["peptide"] = ergo["Peptide"].values

    valid = (
        table["CDR3a"].apply(is_standard_aa_sequence)
        & table["CDR3b"].apply(is_standard_aa_sequence)
        & table["peptide"].apply(is_standard_aa_sequence)
    )
    dropped = (~valid).sum()
    if dropped:
        print(f"Dropping {dropped} rows with a non-standard amino acid in CDR3a/CDR3b/peptide (full immrep set)")
    table = table[valid]

    too_long_pep = table["peptide"].str.len() > MAX_PEPTIDE_LEN
    if too_long_pep.sum():
        print(f"Dropping {too_long_pep.sum()} rows with a peptide longer than {MAX_PEPTIDE_LEN} residues (full immrep set)")
    table = table[~too_long_pep]

    too_long_cdr3 = (table["CDR3a"].str.len() > MAX_CDR3_LEN) | (table["CDR3b"].str.len() > MAX_CDR3_LEN)
    if too_long_cdr3.sum():
        print(f"Dropping {too_long_cdr3.sum()} rows with a CDR3a/CDR3b longer than {MAX_CDR3_LEN} residues (full immrep set)")
    table = table[~too_long_cdr3]

    return table[["name", "CDR3a", "CDR3b", "peptide", "binder"]]


def build_native_train_table():
    # NetTCR-2.0's own published training set - not derived from our Annotation.csv/dat
    # at all, already in the CDR3a/CDR3b/peptide/binder format write_split expects, just
    # with an extra "partition" column (the paper's own CV fold assignment, irrelevant
    # here since the whole file is used as a single training set).
    train = pd.read_csv(NATIVE_TRAIN_CSV)
    return train[["CDR3a", "CDR3b", "peptide", "binder"]]


def write_split(train_df, test_df, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    cols = ["CDR3a", "CDR3b", "peptide", "binder"]
    train_df[cols].to_csv(os.path.join(out_dir, "train.csv"), index=False)
    test_df[cols].to_csv(os.path.join(out_dir, "test.csv"), index=False)
    print(
        f"{out_dir}: train={len(train_df)} (binder={train_df['binder'].value_counts().to_dict()}), "
        f"test={len(test_df)} (binder={test_df['binder'].value_counts().to_dict()})"
    )


def main():
    annotation = pd.read_csv(ANNOTATION_CSV)
    table = build_nettcr_table(annotation)
    print(f"Usable rows with CDR3a/CDR3b: {len(table)} of {len(annotation)}")
    print(table.groupby("source")["binder"].value_counts().unstack(fill_value=0))
    print()

    core = table[table["source"].isin(["ground", "exp", "fake"])]
    immrep = table[table["source"] == "immrep"]
    scan = table[table["source"] == "scan"]
    core_immrep = pd.concat([core, immrep], ignore_index=True)

    train, test = train_test_split(
        core, test_size=TEST_SIZE, stratify=core["binder"], random_state=RANDOM_STATE
    )
    write_split(train, test, os.path.join(OUTPUT_DIR, "core_80_20"))

    write_split(core, immrep, os.path.join(OUTPUT_DIR, "core_vs_immrep"))

    write_split(core, scan, os.path.join(OUTPUT_DIR, "core_vs_scan"))

    train, test = train_test_split(
        core_immrep, test_size=TEST_SIZE, stratify=core_immrep["binder"], random_state=RANDOM_STATE
    )
    write_split(train, test, os.path.join(OUTPUT_DIR, "core_immrep_80_20"))

    # Use the DL pipeline's own fixed train/val assignment rather than re-splitting, so
    # this comparison is scored against the identical holdout set.
    fullrun_split = pd.read_csv(FULLRUN_SPLIT_CSV)
    fullrun_table = core_immrep.merge(fullrun_split[["name", "Set"]], on="name", how="inner")
    train = fullrun_table[fullrun_table["Set"] == "train"]
    test = fullrun_table[fullrun_table["Set"] == "val"]
    write_split(train, test, os.path.join(OUTPUT_DIR, "full_run"))

    # NetTCR-2.0's own published training set, tested on the complete immrep test set
    # (not just the subset with a folded structure in ./dat).
    native_train = build_native_train_table()
    full_immrep = build_full_immrep_table(annotation)
    write_split(native_train, full_immrep, os.path.join(OUTPUT_DIR, "paper_train_vs_immrep_full"))

    # Scanning-dataset generalization: everything we have (core+immrep combined) vs.
    # NetTCR-2.0's own published training set, both tested on the full scanning dataset
    # (all 26 A6/B7 WT + point-mutant samples - all of them have a folded structure in
    # ./dat, so no immrep-style recovery trick is needed here).
    write_split(core_immrep, scan, os.path.join(OUTPUT_DIR, "all_vs_scan"))
    write_split(native_train, scan, os.path.join(OUTPUT_DIR, "paper_vs_scan"))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3

# Builds an ERGO-ready input dataframe for the scanning dataset (the 26 A6/B7 WT +
# point-mutant samples), matching the column layout of data/ERGO/ImmpredAutoresults.csv
# (an existing ERGO *output* file) minus its Score column, which is ERGO's own output,
# not something we provide as input.
#
# ERGO only strictly needs TRA/TRB (CDR3a/CDR3b) and Peptide - unlike NetTCR-2.0, it has
# no fixed-length buffer that requires excluding long peptides/CDR3s, so none of
# PrepNetTCR.py's MAX_PEPTIDE_LEN/MAX_CDR3_LEN filtering applies here. TRAV/TRBV and MHC
# are included too since we have them; TRAJ/TRBJ and T-Cell-Type aren't tracked anywhere
# in this project's data, so they're left blank (same convention the existing
# tcr_test_ergo.csv/tcr_train_ergo.csv files already use for the columns they don't have
# either).

import os

import pandas as pd

from PrepNetTCR import ANNOTATION_CSV, STANDARD_AMINO_ACIDS, extract_cdr3, is_standard_aa_sequence

OUT_FILE = "./data/ERGO/scan_ergo.csv"
ERGO_COLUMNS = ["TRA", "TRB", "TRAV", "TRAJ", "TRBV", "TRBJ", "T-Cell-Type", "Peptide", "MHC"]


def main():
    annotation = pd.read_csv(ANNOTATION_CSV)
    scan = annotation[annotation["source"] == "scan"].copy()

    cdr3 = scan["name"].apply(extract_cdr3)
    scan["TRA"] = [c[0] for c in cdr3]
    scan["TRB"] = [c[1] for c in cdr3]
    missing = scan["TRA"].isna() | scan["TRB"].isna()
    if missing.any():
        print(f"Dropping {missing.sum()} scan samples with no extractable CDR3a/CDR3b: {scan.loc[missing, 'name'].tolist()}")
    scan = scan[~missing]

    valid = (
        scan["TRA"].apply(is_standard_aa_sequence)
        & scan["TRB"].apply(is_standard_aa_sequence)
        & scan["Peptide"].apply(is_standard_aa_sequence)
    )
    if (~valid).any():
        print(f"Dropping {(~valid).sum()} scan samples with a non-standard amino acid in CDR3a/CDR3b/Peptide")
    scan = scan[valid]

    scan["TRAJ"] = ""
    scan["TRBJ"] = ""
    scan["T-Cell-Type"] = ""
    # ImmpredAutoresults.csv's MHC column has no "HLA-" prefix (e.g. "A*02:01"), unlike
    # Annotation.csv's Allele column ("HLA-A*02:01") - stripped to match ERGO's expected
    # format exactly, not just carry the metadata over as-is.
    scan["MHC"] = scan["Allele"].str.removeprefix("HLA-")

    table = scan[["name"] + ERGO_COLUMNS]
    table[ERGO_COLUMNS].to_csv(OUT_FILE, index=False)
    print(f"Wrote {len(table)} rows to {OUT_FILE}")
    print(table.to_string(index=False))


if __name__ == "__main__":
    main()

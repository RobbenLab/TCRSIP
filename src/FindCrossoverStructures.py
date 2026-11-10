#!/usr/bin/env python3

# Scans every structure in ./dat for TCR "crossover" topology - IMGT 110-127 of TCRa
# (the end of the alpha variable domain, going into the elbow region toward the
# constant domain) sitting physically closer to IMGT 1-15 of TCRb (the very start of
# the beta chain) than to IMGT 1-15 of TCRa itself (the same region on its own chain).
# Distance is measured centroid-to-centroid (mean CA position of each residue group),
# a simple, robust proxy for "is this segment nearer the other chain's start than its
# own" rather than being driven by a single close/far atom.
#
# Each sample's own per-sample JSON already carries a literal IMGT number per residue
# (Sequence[chain]['Consensus']['IMGT'], parallel to PDBseq) - no need to reconstruct
# IMGT numbering from the Region dict's CDR/framework boundaries like elsewhere in this
# project. Residue coordinates are recovered the same way as everywhere else here:
# position-in-PDBseq-order maps to position-in-file-order within that chain (NOT resseq
# value, which doesn't correspond 1:1 to sequence index in general).
#
# Second pass: for every structure flagged as crossover, scores how similar its TCRa/
# TCRb sequences are to the XL-MS structures' own (9GV6 chain D/E, 9GV7 chain D/E) via
# global pairwise alignment (BLOSUM62) - the idea being to find a good comparative
# structure for what a crossover conformation might look like on a similar TCR.
#
# Run as a standalone script - scanning ~2600 structures is not something to redo every
# notebook execution. Writes two CSVs to ./data/Results/.

import json
import os

import numpy as np
import pandas as pd
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.PDB import PDBParser

DAT_DIR = "./dat"
XLMS_DIR = "/mnt/d/Lab/TCRBinding/XLMS/Ecolifolded"
XLMS_STRUCTURES = {
    "9GV6": f"{XLMS_DIR}/structures_wlink/alphfold/9gv6.pdb",
    "9GV7": f"{XLMS_DIR}/structures_wlink/alphfold/9gv7.pdb",
}
TCRA_CHAIN, TCRB_CHAIN = "D", "E"  # this project's structure convention (A=MHCa, C=Peptide, D=TCRa, E=TCRb)

CROSSOVER_RANGE = (110, 127)  # TCRa IMGT range being tested
NTERM_RANGE = (1, 15)         # TCRa/TCRb IMGT range compared against

OUT_DIR = "./data/Results"
SCAN_FILE = os.path.join(OUT_DIR, "crossover_scan_all.csv")
SIMILARITY_FILE = os.path.join(OUT_DIR, "crossover_vs_xlms_similarity.csv")

parser = PDBParser(QUIET=True)


def chain_ca_coords_in_order(pdb_path, chain_id):
    """CA coordinates for one chain, in file order (same order as that chain's own
    PDBseq/IMGT arrays in the sample JSON)."""
    structure = parser.get_structure("x", pdb_path)
    model = next(structure.get_models())
    coords = []
    for residue in model[chain_id]:
        if residue.id[0] != " ":
            continue
        if "CA" not in residue:
            coords.append(None)
            continue
        coords.append(residue["CA"].get_coord())
    return coords


def imgt_indices_in_range(imgt_array, lo, hi):
    return [i for i, v in enumerate(imgt_array) if v is not None and lo <= v <= hi]


def centroid(coords, indices):
    pts = [coords[i] for i in indices if i < len(coords) and coords[i] is not None]
    if not pts:
        return None
    return np.mean(np.vstack(pts), axis=0)


def scan_structure(name):
    """Returns a result dict, or a dict with 'skip_reason' set if this structure can't
    be evaluated (missing data, empty IMGT ranges, etc.)."""
    json_path = os.path.join(DAT_DIR, name, f"{name}.json")
    pdb_path = os.path.join(DAT_DIR, name, "Structure", "TCR.pdb")
    if not os.path.exists(json_path) or not os.path.exists(pdb_path):
        return {"name": name, "skip_reason": "missing json or Structure/TCR.pdb"}

    try:
        with open(json_path) as f:
            data = json.load(f)
        tcra = data["Sequence"]["TCRa"]
        tcrb = data["Sequence"]["TCRb"]
        imgt_a, imgt_b = tcra["Consensus"]["IMGT"], tcrb["Consensus"]["IMGT"]
        if imgt_a is None or imgt_b is None:
            return {"name": name, "skip_reason": "no IMGT numbering for TCRa/TCRb"}

        a_cross_idx = imgt_indices_in_range(imgt_a, *CROSSOVER_RANGE)
        a_nterm_idx = imgt_indices_in_range(imgt_a, *NTERM_RANGE)
        b_nterm_idx = imgt_indices_in_range(imgt_b, *NTERM_RANGE)
        if not a_cross_idx or not a_nterm_idx or not b_nterm_idx:
            return {"name": name, "skip_reason": "one or more IMGT ranges empty for this structure"}

        coords_a = chain_ca_coords_in_order(pdb_path, TCRA_CHAIN)
        coords_b = chain_ca_coords_in_order(pdb_path, TCRB_CHAIN)
        if len(coords_a) != len(imgt_a) or len(coords_b) != len(imgt_b):
            return {"name": name, "skip_reason": "chain residue count doesn't match JSON IMGT array length"}

        cross_centroid = centroid(coords_a, a_cross_idx)
        a_nterm_centroid = centroid(coords_a, a_nterm_idx)
        b_nterm_centroid = centroid(coords_b, b_nterm_idx)
        if cross_centroid is None or a_nterm_centroid is None or b_nterm_centroid is None:
            return {"name": name, "skip_reason": "no resolved CA atoms in one of the residue groups"}

        dist_same = float(np.linalg.norm(cross_centroid - a_nterm_centroid))
        dist_cross = float(np.linalg.norm(cross_centroid - b_nterm_centroid))

        return {
            "name": name,
            "skip_reason": None,
            "dist_to_own_Nterm_A": dist_same,
            "dist_to_other_Nterm_B": dist_cross,
            "crossover": dist_cross < dist_same,
            "TCRa_len": len(tcra["PDBseq"]),
            "TCRb_len": len(tcrb["PDBseq"]),
            "TCRa_seq": tcra["PDBseq"],
            "TCRb_seq": tcrb["PDBseq"],
        }
    except Exception as exc:
        return {"name": name, "skip_reason": f"error: {exc}"}


def load_xlms_chain_seq(pdb_path, chain_id):
    structure = parser.get_structure("x", pdb_path)
    model = next(structure.get_models())
    from Bio.PDB.Polypeptide import three_to_index, index_to_one
    seq = []
    for residue in model[chain_id]:
        if residue.id[0] != " ":
            continue
        try:
            seq.append(index_to_one(three_to_index(residue.get_resname())))
        except KeyError:
            seq.append("X")
    return "".join(seq)


def alignment_identity(seq1, seq2, aligner):
    if not seq1 or not seq2:
        return 0.0
    alignment = aligner.align(seq1, seq2)[0]
    a1, a2 = alignment[0], alignment[1]
    matches = sum(1 for x, y in zip(a1, a2) if x == y and x != "-")
    return matches / min(len(seq1), len(seq2))


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    annotation = pd.read_csv("./data/ML/Annotation.csv")[["name", "source"]].set_index("name")["source"]

    names = sorted(os.listdir(DAT_DIR))
    print(f"Scanning {len(names)} structures under {DAT_DIR} ...")
    results = [scan_structure(name) for name in names]
    scan_df = pd.DataFrame(results)
    scan_df["source"] = scan_df["name"].map(annotation)
    scan_df.to_csv(SCAN_FILE, index=False)

    n_skipped = scan_df["skip_reason"].notna().sum()
    n_crossover = (scan_df["crossover"] == True).sum() if "crossover" in scan_df else 0
    print(f"{n_skipped} skipped, {n_crossover} crossover of {len(scan_df) - n_skipped} evaluated")
    print(f"Full scan saved to {SCAN_FILE}")

    crossover_df = scan_df[scan_df["crossover"] == True].copy()
    if crossover_df.empty:
        print("No crossover structures found - skipping XL-MS similarity ranking.")
        return

    xlms_seqs = {
        label: {"TCRa": load_xlms_chain_seq(pdb, TCRA_CHAIN), "TCRb": load_xlms_chain_seq(pdb, TCRB_CHAIN)}
        for label, pdb in XLMS_STRUCTURES.items()
    }

    aligner = PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    aligner.mode = "global"

    for label, seqs in xlms_seqs.items():
        crossover_df[f"TCRa_identity_{label}"] = crossover_df["TCRa_seq"].apply(
            lambda s: alignment_identity(s, seqs["TCRa"], aligner)
        )
        crossover_df[f"TCRb_identity_{label}"] = crossover_df["TCRb_seq"].apply(
            lambda s: alignment_identity(s, seqs["TCRb"], aligner)
        )
        crossover_df[f"mean_identity_{label}"] = crossover_df[[f"TCRa_identity_{label}", f"TCRb_identity_{label}"]].mean(axis=1)
        crossover_df[f"TCRa_len_diff_{label}"] = (crossover_df["TCRa_len"] - len(seqs["TCRa"])).abs()
        crossover_df[f"TCRb_len_diff_{label}"] = (crossover_df["TCRb_len"] - len(seqs["TCRb"])).abs()

    crossover_df = crossover_df.sort_values("mean_identity_9GV6", ascending=False)
    crossover_df.drop(columns=["TCRa_seq", "TCRb_seq"]).to_csv(SIMILARITY_FILE, index=False)
    print(f"Crossover-vs-XLMS similarity saved to {SIMILARITY_FILE}")

    for label in XLMS_STRUCTURES:
        top = crossover_df.sort_values(f"mean_identity_{label}", ascending=False).head(5)
        print(f"\nTop 5 crossover structures most similar to {label}:")
        print(top[["name", "source", f"TCRa_identity_{label}", f"TCRb_identity_{label}",
                    f"mean_identity_{label}", "TCRa_len", "TCRb_len"]].to_string(index=False))


if __name__ == "__main__":
    main()

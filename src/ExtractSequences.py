#!/usr/bin/env python3

# Takes a directory of .pdb files, extracts the residue sequence of each chain (A-E),
# and exports a CSV with one row per structure and one column per chain.

import os
import argparse
import warnings

import pandas as pd
from tqdm import tqdm
from Bio.PDB import PDBParser
from Bio import BiopythonWarning

warnings.simplefilter('ignore', BiopythonWarning)

seq_dict = {'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
            'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
            'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
            'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'}

CHAIN_MAP = {"A": "MHCa", "B": "MHCb", "C": "Peptide", "D": "TCRa", "E": "TCRb"}


def getChainSequence(chain):
    #Non-standard/unknown residues are kept and marked as X rather than skipped
    return "".join(seq_dict.get(res.get_resname(), "X") for res in chain)


def extractSequences(pdb_file):
    parser = PDBParser(QUIET=True)
    name = os.path.basename(pdb_file).split(".")[0][:4]
    structure = parser.get_structure(name, pdb_file)
    model = structure[0]

    row = {"ID": name}
    for chain_id, col_name in CHAIN_MAP.items():
        row[col_name] = getChainSequence(model[chain_id]) if chain_id in model else ""
    return row


def extractDirectory(input_dir):
    files = [f for f in os.listdir(input_dir) if f.lower().endswith(".pdb")]
    rows = [extractSequences(os.path.join(input_dir, file)) for file in tqdm(files)]
    return pd.DataFrame(rows, columns=["ID"] + list(CHAIN_MAP.values()))


if __name__ == "__main__":
    argparser = argparse.ArgumentParser(description="Extract per-chain residue sequences from PDB files into a CSV")
    argparser.add_argument("input", help="Directory containing .pdb files")
    argparser.add_argument("-o", "--output", default="./sequences.csv", help="Path to output CSV file")
    args = argparser.parse_args()

    df = extractDirectory(args.input)
    df.to_csv(args.output, index=False)

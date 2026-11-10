#!/usr/bin/env python3

# Refines the 9GV6/9GV7 AlphaFold structures (/mnt/d/Lab/TCRBinding/XLMS/Ecolifolded/)
# against their XL-MS crosslink restraints using PyRosetta: AtomPairConstraints on each
# crosslinked pair's CA-CA distance, using a FlatHarmonicFunc centered at 0 - zero
# penalty up to MAX_DISTANCE Angstroms, harmonic beyond - since a crosslinker is a
# MAXIMUM reach, not a fixed target distance, so only genuinely-too-long pairs should be
# pulled in. Minimized in torsion (dihedral) space rather than Cartesian space: bond
# lengths/angles aren't degrees of freedom there at all, so they stay exactly fixed
# throughout - "minimize crosslink distances while keeping all other bond lengths equal"
# is what torsion-space minimization does by construction, not as an approximation.
# Inter-chain crosslinks (e.g. TCR alpha to TCR beta) also need the pose's rigid-body
# jumps between chains marked flexible, since intra-chain torsion moves alone can't
# bring two separate chains closer together.
#
# MAX_DISTANCE: sulfo-SDA's own spacer-arm length (reactive-group to reactive-group) is
# ~3.9 A, but that excludes the lysine side chain's own reach (~6.5 A) and the
# photoreactive group's reach on the other end - a CA-CA restraint needs to account for
# both, so this uses 25 A, a commonly-cited effective CA-CA cutoff for SDA-type
# crosslinkers in the XL-MS structural literature. Adjust below if a different value is
# wanted - this is the single most consequential modeling choice in this script.
#
# Run as a standalone script (not from the notebook): PyRosetta refinement takes real
# time (not something to redo every notebook execution), and the notebook just loads the
# refined PDB this script writes out to visualize it.

import csv
import os

import pandas as pd
import pyrosetta
from pyrosetta.rosetta.core.id import AtomID
from pyrosetta.rosetta.core.kinematics import MoveMap
from pyrosetta.rosetta.core.scoring import ScoreType, get_score_function
from pyrosetta.rosetta.core.scoring.constraints import AtomPairConstraint
from pyrosetta.rosetta.core.scoring.func import FlatHarmonicFunc
from pyrosetta.rosetta.protocols.minimization_packing import MinMover

XLMS_DIR = "/mnt/d/Lab/TCRBinding/XLMS/Ecolifolded"
CHAIN_SUFFIX_MAP = {1: "A", 2: "B", 3: "C", 4: "D", 5: "E"}

STRUCTURES = [
    {
        "name": "9GV6",
        "pdb": f"{XLMS_DIR}/structures_wlink/alphfold/9gv6.pdb",
        "xvis_csv": f"{XLMS_DIR}/results/TCRDxvis.csv",
        "out_pdb": f"{XLMS_DIR}/structures_wlink/alphfold/9gv6_xlrefined.pdb",
        "report_csv": f"{XLMS_DIR}/structures_wlink/alphfold/9gv6_xlrefined_report.csv",
    },
    {
        "name": "9GV7",
        "pdb": f"{XLMS_DIR}/structures_wlink/alphfold/9gv7.pdb",
        "xvis_csv": f"{XLMS_DIR}/results/TCRExvis.csv",
        "out_pdb": f"{XLMS_DIR}/structures_wlink/alphfold/9gv7_xlrefined.pdb",
        "report_csv": f"{XLMS_DIR}/structures_wlink/alphfold/9gv7_xlrefined_report.csv",
    },
]

MAX_DISTANCE = 25.0   # Angstroms - see note above
RESTRAINT_SD = 5.0    # FlatHarmonicFunc width beyond MAX_DISTANCE
CONSTRAINT_WEIGHT = 5.0  # atom_pair_constraint score term weight
N_CYCLES = 3           # minimization rounds


def parse_crosslinks(xvis_csv):
    """
    Same parsing as the notebook's parse_crosslinks (duplicated here so this script has
    no dependency on the notebook): the structure prefix in Protein1/Protein2 (9GV6 vs
    9GV7) is an annotation artifact and is ignored, only the chain-number suffix is
    used. AbsPos is 0-indexed - +1 to get the PDB's own resseq.
    """
    df = pd.read_csv(xvis_csv)

    def parse_chain_num(value):
        return int(value.split("_", 1)[1].split("|")[0])

    links = set()
    for _, row in df.iterrows():
        pos1, pos2 = int(row["AbsPos1"]) + 1, int(row["AbsPos2"]) + 1
        chain1 = CHAIN_SUFFIX_MAP[parse_chain_num(row["Protein1"])]
        chain2 = CHAIN_SUFFIX_MAP[parse_chain_num(row["Protein2"])]
        links.add(tuple(sorted([(chain1, pos1), (chain2, pos2)])))
    return sorted(links)


def refine_with_crosslinks(pdb_path, crosslinks, max_distance=MAX_DISTANCE, restraint_sd=RESTRAINT_SD,
                            constraint_weight=CONSTRAINT_WEIGHT, n_cycles=N_CYCLES):
    """
    Refines pdb_path against `crosslinks` (a list of ((chain, resseq), (chain, resseq))
    pairs). Returns (pose, report) - report is a list of dicts with before/after CA-CA
    distances per crosslink.
    """
    pose = pyrosetta.pose_from_pdb(pdb_path)
    pdb_info = pose.pdb_info()

    def ca_resnum(chain, resseq):
        resnum = pdb_info.pdb2pose(chain, resseq)
        if resnum == 0:
            raise ValueError(f"Residue {chain}{resseq} not found in {pdb_path}")
        return resnum

    report = []
    constrained_pairs = []
    for (c1, r1), (c2, r2) in crosslinks:
        resnum1, resnum2 = ca_resnum(c1, r1), ca_resnum(c2, r2)
        atom_id1 = AtomID(pose.residue(resnum1).atom_index("CA"), resnum1)
        atom_id2 = AtomID(pose.residue(resnum2).atom_index("CA"), resnum2)
        before = pose.residue(resnum1).xyz("CA").distance(pose.residue(resnum2).xyz("CA"))

        func = FlatHarmonicFunc(0.0, restraint_sd, max_distance)
        pose.add_constraint(AtomPairConstraint(atom_id1, atom_id2, func))

        report.append({"chain1": c1, "resseq1": r1, "chain2": c2, "resseq2": r2, "before_A": before})
        constrained_pairs.append((resnum1, resnum2))

    scorefxn = get_score_function()
    scorefxn.set_weight(ScoreType.atom_pair_constraint, constraint_weight)

    movemap = MoveMap()
    movemap.set_bb(True)
    movemap.set_chi(True)
    for jump in range(1, pose.num_jump() + 1):
        movemap.set_jump(jump, True)

    min_mover = MinMover()
    min_mover.movemap(movemap)
    min_mover.score_function(scorefxn)
    min_mover.min_type("lbfgs_armijo_nonmonotone")
    min_mover.tolerance(0.0001)

    for _ in range(n_cycles):
        min_mover.apply(pose)

    for row, (resnum1, resnum2) in zip(report, constrained_pairs):
        row["after_A"] = pose.residue(resnum1).xyz("CA").distance(pose.residue(resnum2).xyz("CA"))

    return pose, report


def main():
    pyrosetta.init("-mute all -ex1 -ex2aro")

    for spec in STRUCTURES:
        print(f"=== {spec['name']} ===")
        crosslinks = parse_crosslinks(spec["xvis_csv"])
        print(f"{len(crosslinks)} crosslinks parsed from {spec['xvis_csv']}")

        pose, report = refine_with_crosslinks(spec["pdb"], crosslinks)

        print(f"{'chain1':7}{'res1':6}{'chain2':7}{'res2':6}{'before (A)':12}{'after (A)':12}")
        for row in report:
            print(f"{row['chain1']:7}{row['resseq1']:<6}{row['chain2']:7}{row['resseq2']:<6}"
                  f"{row['before_A']:<12.1f}{row['after_A']:<12.1f}")

        violated_before = sum(1 for r in report if r["before_A"] > MAX_DISTANCE)
        violated_after = sum(1 for r in report if r["after_A"] > MAX_DISTANCE)
        print(f"Violations (> {MAX_DISTANCE} A): {violated_before} before -> {violated_after} after")

        pose.dump_pdb(spec["out_pdb"])
        print(f"Refined structure written to {spec['out_pdb']}")

        with open(spec["report_csv"], "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["chain1", "resseq1", "chain2", "resseq2", "before_A", "after_A"])
            writer.writeheader()
            writer.writerows(report)
        print(f"Report saved to {spec['report_csv']}")
        print()


if __name__ == "__main__":
    main()

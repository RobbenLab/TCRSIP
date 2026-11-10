#!/usr/bin/env python3

# Builds the 2D residue-pair input maps for the TCR ResNet from processed structure
# directories (the same layout ExtractFeaturestoJSON.py uses):
#
#   <in_dir>/<name>/<name>.json         sequences + regions
#   <in_dir>/<name>/Structure/Complex.pdb
#
# The channels (CHANNELS) and trimming rules (TRIM) are defined below (the load-time
# canvas, CANVAS, in tcr_dataset.py) and documented in Input_def.txt; each channel name maps to a feature function
# (@provides). Each complex is trimmed (MHC peptide-binding domains, peptide, TCR V
# domains), then written to <out_dir>/<name>.npz with:
#   x             float16 (C, L, L) map of the trimmed complex, unpadded, physical units
#   channels      channel names (C,)
#   groups        channel groups (C,), e.g. for holding out "physics"
#   seq           one-letter trimmed sequence in canonical chain order
#   chain         (L,) chain index 0-4 = MHCa, MHCb, Peptide, TCRa, TCRb
#   anchor_index  index of the last peptide residue (centered on the canvas at load time)
#   trim          (5, 2) kept [start, end) per chain in original chain offsets
#   flags         trimming warnings (e.g. no_J_a, short_tcr_b, cdr_from_anchors_a)
#   name, mhc_class, peptide, cdr3a, cdr3b   identifiers derived from the JSON/structure
# plus Input_structures.csv - one row per .npz in out_dir (File relative to the CSV, Name,
# MHC class, peptide, CDR3s, lengths, flags), rebuilt from the files after every run; it
# holds no labels, so it works the same for training and inference sets - and
# build_log.txt (skipped/failed structures).
#
#   conda activate tcr_dl
#   python src/build_inputs.py -i data/structures -o data/input --workers 32
#   python src/build_inputs.py -o data/input --index-only        # just rewrite the CSV

import argparse
import json
import os
import re
import tempfile
import traceback
import warnings
from multiprocessing import Pool

import numpy as np
import pandas as pd
from Bio.Align import substitution_matrices
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
from scipy.spatial import cKDTree
from tqdm import tqdm

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

CHAINS = ["MHCa", "MHCb", "Peptide", "TCRa", "TCRb"]
CHAIN_IDS = {"MHCa": "A", "MHCb": "B", "Peptide": "C", "TCRa": "D", "TCRb": "E"}
AA_ORDER = "ARNDCQEGHILKMFPSTWYV"
BACKBONE = {"N", "CA", "C", "O", "OXT"}
# Region channel -> JSON Region names it covers
REGION_CHANNELS = {"MHCha": ["MHCha"], "MHChb": ["MHChb"], "CDR1": ["CDR1a", "CDR1b"],
                   "CDR2": ["CDR2a", "CDR2b"], "CDR3": ["CDR3a", "CDR3b"]}
# Residue classes from getDistanceFeatures in ExtractFeaturestoJSON.py; one group = both
# residues in it, two groups = one residue from each
RESIDUE_TYPE_PAIRS = {
    "polar": ["STNQ"],
    "nonpolar": ["AVILMFYW"],
    "charged": ["RHK", "DE"],
    "aromatic": ["FYW"],
}
IMGT_CDRS = {1: (27, 38), 2: (56, 65), 3: (105, 117)}
CONTACT_CUTOFF = 4.5
SEQSEP_MAX = 64
FP16_MAX = 65504.0


##############################
###### Channel definitions and settings
##############################
# Input_def.txt documents these channels and rules; keep it in sync when changing them.

ROSETTA_TERMS = ["fa_atr", "fa_rep", "fa_sol", "lk_ball_wtd", "fa_elec",
                 "hbond_sr_bb", "hbond_lr_bb", "hbond_bb_sc", "hbond_sc", "dslf_fa13"]

# (name, group, symmetry) in channel order
CHANNELS = (
    [(f"seq_blosum_diff_{a}", "sequence", "anti") for a in AA_ORDER]
    + [("seq_blosum_score", "sequence", "sym")]
    + [(f"pos_chain_{c}", "position", "sym") for c in CHAINS]
    + [("pos_seqsep", "position", "sym")]
    + [(f"pos_region_{r}", "position", "sym") for r in REGION_CHANNELS]
    + [("str_ca_dist", "structure", "sym"),
       ("str_cb_dist", "structure", "sym"),
       ("str_sc_sc_contacts", "structure", "sym"),
       ("str_sc_bb_contacts", "structure", "asym"),
       ("str_bb_bb_contacts", "structure", "sym")]
    + [(f"str_{t}_dist", "structure", "sym") for t in RESIDUE_TYPE_PAIRS]
    + [(f"phys_{t}", "physics", "sym") for t in ROSETTA_TERMS]
)

# Trimming rules (see Input_def.txt for the reasoning)
TRIM = {
    "mhc1_keep": 180,            # class I MHCa: HLA residues 1-180 (alpha1 + alpha2)
    "mhc1_helix_a": (49, 83),    # class I helix [start, end) offsets from HLA residue 1
    "mhc1_helix_b": (138, 175),
    "mhc2_max_len": 100,         # flag class II chains longer than this
    "tcr_pre_cys": 22,           # TCR start = Cys23 - 22
    "tcr_post_j": 10,            # TCR end = F118 + 10 (IMGT 128)
    "tcr_fallback_len": 120,     # TCR length when the J motif is missing
    "tcr_min_len": 60,           # flag shorter (fragment) TCR chains
}

# Load-time canvas (last peptide residue at the anchor) lives with the loader so the
# loader stays standalone; imported here only to check that built maps fit
from tcr_dataset import CANVAS  # noqa: E402


# Registry: channel name -> function computing a dict {name: (L, L) array} for one complex.
# A function can provide several channels (e.g. all 20 BLOSUM columns) and is run at most
# once per complex.
FEATURES = {}


def provides(*names):
    def register(fn):
        for n in names:
            FEATURES[n] = fn
        return fn
    return register


##############################
###### Structure loading
##############################

class SkipComplex(Exception):
    """Structure intentionally not built (e.g. no peptide / non-classical MHC)."""


def tcr_anchors(seq, st):
    """Trim range and conserved anchors for a TCR chain (rules in Input_def.txt)."""
    flags = []
    m = re.search(r"C.{9,17}W[YF]", seq[:70])
    cys23 = m.start() if m else None
    trp41 = m.end() - 2 if m else None
    start = max(0, cys23 - st["tcr_pre_cys"]) if m else 0
    if not m:
        flags.append("no_cys23")
    j = [x.start() for x in re.finditer(r"[FW]G.G", seq[:170])]
    f118 = j[-1] if j else None
    if j:
        end = min(len(seq), f118 + st["tcr_post_j"])
    else:
        end = min(len(seq), start + st["tcr_fallback_len"])
        flags.append("no_J")
    return dict(start=start, end=end, cys23=cys23, trp41=trp41, f118=f118, flags=flags)


def mhc1_anchors(seq, st):
    m = re.search(r"H[ST][LM]RY", seq[:40])
    start = max(0, m.start() - 2) if m else 0
    return dict(start=start, end=min(len(seq), start + st["mhc1_keep"]),
                flags=[] if m else ["no_mhc1_motif"])


class Complex:
    """A TCR-pMHC complex with residues in canonical chain order (MHCa, MHCb, Peptide,
    TCRa, TCRb), validated against the JSON sequences and trimmed with the rules in
    Input_def.txt."""

    def __init__(self, struct_dir, settings=TRIM, pdb_name="Complex.pdb"):
        self.name = os.path.basename(os.path.normpath(struct_dir))
        self.dir = struct_dir
        with open(os.path.join(struct_dir, f"{self.name}.json")) as f:
            self.js = json.load(f)
        self.mhc_class = self.js.get("Features", {}).get("Sequence", {}).get("MHC_class")
        if self.mhc_class not in (1, 2) or not self.js["Sequence"]["Peptide"]["Seq"]:
            raise SkipComplex(f"no peptide or MHC class (class={self.mhc_class})")
        # Some sources name it complex.pdb, so match the file name case-insensitively
        sdir = os.path.join(struct_dir, "Structure")
        match = [f for f in os.listdir(sdir) if f.lower() == pdb_name.lower()] if os.path.isdir(sdir) else []
        if not match:
            raise FileNotFoundError(f"No {pdb_name} in {sdir}")
        self.pdb_file = os.path.join(sdir, match[0])
        structure = PDBParser(QUIET=True).get_structure(self.name, self.pdb_file)
        model = next(iter(structure))
        st = settings
        self.flags = []

        residues, chain_idx, seq, self.trim, regions = [], [], [], {}, {}
        for c, chain_name in enumerate(CHAINS):
            js_seq = self.js["Sequence"][chain_name]["Seq"]
            cid = CHAIN_IDS[chain_name]
            chain_res = [r for r in model[cid] if r.id[0] == " "] if cid in model else []
            pdb_seq = "".join(seq1(r.get_resname()) for r in chain_res)
            if pdb_seq != js_seq:
                raise ValueError(f"{self.name} {chain_name} (chain {cid}): PDB sequence does not match JSON\n"
                                 f"  PDB:  {pdb_seq}\n  JSON: {js_seq}")

            # Trim range [a, b) in chain coordinates and regions in chain coordinates
            a, b = 0, len(js_seq)
            if chain_name == "MHCa" and self.mhc_class == 1:
                an = mhc1_anchors(js_seq, st)
                a, b = an["start"], an["end"]
                self.flags += an["flags"]
                regions["MHCha"] = (cid, a + st["mhc1_helix_a"][0], a + st["mhc1_helix_a"][1])
                regions["MHChb"] = (cid, a + st["mhc1_helix_b"][0], a + st["mhc1_helix_b"][1])
            elif chain_name in ("MHCa", "MHCb") and self.mhc_class == 2:
                if len(js_seq) > st["mhc2_max_len"]:
                    self.flags.append("mhc2_long")
                helix = "MHCha" if chain_name == "MHCa" else "MHChb"
                r = self.js["Region"].get(chain_name, {}).get(helix)
                if r and r[1] > r[0]:
                    regions[helix] = (cid, r[0], r[1])
            elif chain_name in ("TCRa", "TCRb"):
                an = tcr_anchors(js_seq, st)
                a, b = an["start"], an["end"]
                sfx = chain_name[-1]
                self.flags += [f"{f}_{sfx}" for f in an["flags"]]
                if b - a < st["tcr_min_len"]:
                    self.flags.append(f"short_tcr_{sfx}")
                cdrs, from_anchors = self.cdr_ranges(chain_name, js_seq, an)
                if from_anchors:
                    self.flags.append(f"cdr_from_anchors_{sfx}")
                for k, (lo, hi) in cdrs.items():
                    regions[f"CDR{k}{sfx}"] = (cid, lo, hi)

            self.trim[chain_name] = (a, b)
            kept = chain_res[a:b]
            residues += kept
            chain_idx += [c] * len(kept)
            seq.append(js_seq[a:b])

        self.residues = residues
        self.chain = np.array(chain_idx, dtype=np.int8)
        self.seq = "".join(seq)
        self.L = len(residues)
        self.chain_start = {n: int(np.argmax(self.chain == c)) if (self.chain == c).any() else None
                            for c, n in enumerate(CHAINS)}
        # Index of the last peptide residue in the trimmed map (the centering anchor)
        self.anchor_index = int(np.where(self.chain == CHAINS.index("Peptide"))[0][-1])
        # Regions -> boolean masks in trimmed map coordinates
        self.regions = {}
        for name, (cid, lo, hi) in regions.items():
            chain_name = [k for k, v in CHAIN_IDS.items() if v == cid][0]
            a, b = self.trim[chain_name]
            lo, hi = max(lo, a) - a, min(hi, b) - a
            mask = np.zeros(self.L, dtype=bool)
            if hi > lo:
                s = self.chain_start[chain_name]
                mask[s + lo:s + hi] = True
            self.regions[name] = mask

    def cdr_ranges(self, chain_name, seq, an):
        """CDR [start, end) chain offsets from the JSON IMGT numbering, else from anchors.
        Returns ({1: (lo, hi), 2: ..., 3: ...}, from_anchors)."""
        cons = self.js["Sequence"][chain_name].get("Consensus") or {}
        imgt = cons.get("IMGT") if isinstance(cons, dict) else None
        if imgt and len(imgt) == len(seq):
            num = np.array([x if isinstance(x, (int, float)) else -1 for x in imgt])
            out = {}
            for k, (lo, hi) in IMGT_CDRS.items():
                idx = np.where((num >= lo) & (num <= hi))[0]
                if len(idx):
                    out[k] = (int(idx[0]), int(idx[-1]) + 1)
            # The JSON numbering runs CDR3 into the J region for ~25% of chains; CDR3 ends
            # at IMGT 117, right before the J motif's F118
            if 3 in out and an["f118"] is not None and out[3][1] > an["f118"]:
                out[3] = (out[3][0], an["f118"])
            if len(out) == 3:
                return out, False
        out = {}
        if an["cys23"] is not None:
            out[1] = (an["cys23"] + 4, an["trp41"] - 3)
            out[2] = (an["trp41"] + 15, an["trp41"] + 25)
        if an["f118"] is not None:
            cys104 = seq.rfind("C", 0, an["f118"])
            if cys104 >= 0:
                out[3] = (cys104 + 1, an["f118"])
        return out, True

    def region_mask(self, region_names):
        mask = np.zeros(self.L, dtype=bool)
        for r in region_names:
            if r in self.regions:
                mask |= self.regions[r]
        return mask

    def write_trimmed_pdb(self, path):
        from Bio.PDB import PDBIO, Select
        # PDBIO copies the model, so match residues by (chain, residue id), not identity
        keep = {(r.get_parent().id, r.id) for r in self.residues}

        class Keep(Select):
            def accept_residue(self, residue):
                return (residue.get_parent().id, residue.id) in keep

        io = PDBIO()
        io.set_structure(self.residues[0].get_parent().get_parent())
        io.save(path, Keep())

    def atom_coord(self, res, atom):
        return res[atom].coord if atom in res else np.full(3, np.nan, dtype=np.float32)


def pairwise_dist(coords):
    return np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=-1)


##############################
###### Sequence features
##############################

BLOSUM62 = substitution_matrices.load("BLOSUM62")


def blosum_row(aa):
    aa = aa if aa in AA_ORDER else "X"
    return np.array([BLOSUM62[aa][k] for k in AA_ORDER], dtype=np.float32)


@provides(*[f"seq_blosum_diff_{a}" for a in AA_ORDER])
def seq_blosum_diff(cx):
    sq = np.stack([blosum_row(a) for a in cx.seq]) ** 2  # (L, 20)
    diff = sq[:, None, :] - sq[None, :, :]               # (L, L, 20)
    return {f"seq_blosum_diff_{a}": diff[:, :, k] for k, a in enumerate(AA_ORDER)}


@provides("seq_blosum_score")
def seq_blosum_score(cx):
    seq = [a if a in AA_ORDER else "X" for a in cx.seq]
    score = np.array([[BLOSUM62[a][b] for b in seq] for a in seq], dtype=np.float32)
    return {"seq_blosum_score": score}


##############################
###### Position features
##############################

@provides(*[f"pos_chain_{c}" for c in CHAINS])
def pos_chain(cx):
    out = {}
    for c, name in enumerate(CHAINS):
        on = cx.chain == c
        out[f"pos_chain_{name}"] = (on[:, None] | on[None, :]).astype(np.float32)
    return out


@provides("pos_seqsep")
def pos_seqsep(cx):
    idx = np.arange(cx.L)
    sep = np.minimum(np.abs(idx[:, None] - idx[None, :]), SEQSEP_MAX).astype(np.float32)
    sep[cx.chain[:, None] != cx.chain[None, :]] = SEQSEP_MAX + 1
    return {"pos_seqsep": sep}


@provides(*[f"pos_region_{r}" for r in REGION_CHANNELS])
def pos_region(cx):
    out = {}
    for name, regions in REGION_CHANNELS.items():
        m = cx.region_mask(regions).astype(np.float32)
        out[f"pos_region_{name}"] = m[:, None] + m[None, :]
    return out


##############################
###### Structure features
##############################

def virtual_cb(n, ca, c):
    # Ideal CB position from backbone atoms (trRosetta), used for glycine
    b = ca - n
    c_ = c - ca
    a = np.cross(b, c_)
    return -0.58273431 * a + 0.56802827 * b - 0.54067466 * c_ + ca


@provides("str_ca_dist")
def str_ca_dist(cx):
    ca = np.stack([cx.atom_coord(r, "CA") for r in cx.residues])
    return {"str_ca_dist": pairwise_dist(ca)}


@provides("str_cb_dist")
def str_cb_dist(cx):
    cb = []
    for r in cx.residues:
        if "CB" in r:
            cb.append(r["CB"].coord)
        else:
            cb.append(virtual_cb(*(cx.atom_coord(r, a) for a in ("N", "CA", "C"))))
    return {"str_cb_dist": pairwise_dist(np.stack(cb))}


@provides(*[f"str_{t}_dist" for t in RESIDUE_TYPE_PAIRS])
def str_type_dist(cx):
    # Side-chain heavy-atom centroids; glycine (no side chain) is in none of the classes
    cen = []
    for r in cx.residues:
        sc = [a.coord for a in r if a.element != "H" and a.get_id() not in BACKBONE]
        cen.append(np.mean(sc, axis=0) if sc else cx.atom_coord(r, "CA"))
    dist = pairwise_dist(np.stack(cen))
    seq = np.array(list(cx.seq))
    out = {}
    for name, groups in RESIDUE_TYPE_PAIRS.items():
        g1 = np.isin(seq, list(groups[0]))
        g2 = np.isin(seq, list(groups[-1]))
        pair = (g1[:, None] & g2[None, :]) | (g2[:, None] & g1[None, :])
        np.fill_diagonal(pair, False)
        out[f"str_{name}_dist"] = np.where(pair, dist, 0).astype(np.float32)
    return out


@provides("str_sc_sc_contacts", "str_sc_bb_contacts", "str_bb_bb_contacts")
def str_contacts(cx):
    coords, res_idx, is_sc = [], [], []
    for i, r in enumerate(cx.residues):
        for atom in r:
            if atom.element == "H":
                continue
            coords.append(atom.coord)
            res_idx.append(i)
            is_sc.append(atom.get_id() not in BACKBONE)
    coords = np.array(coords)
    res_idx = np.array(res_idx)
    is_sc = np.array(is_sc)

    pairs = cKDTree(coords).query_pairs(CONTACT_CUTOFF, output_type="ndarray")
    a, b = pairs[:, 0], pairs[:, 1]
    keep = res_idx[a] != res_idx[b]
    a, b = a[keep], b[keep]
    ra, rb = res_idx[a], res_idx[b]
    sa, sb = is_sc[a], is_sc[b]

    L = cx.L
    sc_sc = np.zeros((L, L), dtype=np.float32)
    sc_bb = np.zeros((L, L), dtype=np.float32)
    bb_bb = np.zeros((L, L), dtype=np.float32)
    m = sa & sb
    np.add.at(sc_sc, (ra[m], rb[m]), 1)
    np.add.at(sc_sc, (rb[m], ra[m]), 1)
    m = ~sa & ~sb
    np.add.at(bb_bb, (ra[m], rb[m]), 1)
    np.add.at(bb_bb, (rb[m], ra[m]), 1)
    # sc_bb[i, j] counts side-chain atoms of i near backbone atoms of j
    m = sa & ~sb
    np.add.at(sc_bb, (ra[m], rb[m]), 1)
    m = ~sa & sb
    np.add.at(sc_bb, (rb[m], ra[m]), 1)
    return {"str_sc_sc_contacts": sc_sc, "str_sc_bb_contacts": sc_bb, "str_bb_bb_contacts": bb_bb}


##############################
###### Physics features (Rosetta)
##############################

_rosetta = {}


def init_rosetta():
    if not _rosetta:
        import pyrosetta
        pyrosetta.init("-mute all", silent=True)
        _rosetta["pyrosetta"] = pyrosetta
        scorefxn = pyrosetta.get_fa_scorefxn()
        # Backbone-backbone H-bonds are otherwise only added to the pose total, not to the
        # residue-pair energy graph, so hbond_sr_bb/hbond_lr_bb maps would be all zero
        opts = scorefxn.energy_method_options()
        opts.hbond_options().decompose_bb_hb_into_pair_energies(True)
        scorefxn.set_energy_method_options(opts)
        _rosetta["scorefxn"] = scorefxn
    return _rosetta["pyrosetta"], _rosetta["scorefxn"]


@provides(*[f"phys_{t}" for t in ROSETTA_TERMS])
def phys_rosetta(cx):
    pyrosetta, scorefxn = init_rosetta()
    from pyrosetta.rosetta.core.scoring import EMapVector, ScoreType
    from pyrosetta.rosetta.core.scoring.methods import LongRangeEnergyType

    # Score the trimmed complex so energies only reflect the residues in the map
    with tempfile.NamedTemporaryFile(suffix=".pdb") as tmp:
        cx.write_trimmed_pdb(tmp.name)
        pose = pyrosetta.pose_from_pdb(tmp.name)
    scorefxn(pose)
    weights = scorefxn.weights()
    terms = [(t, getattr(ScoreType, t)) for t in ROSETTA_TERMS]

    # Map pose residue numbers to canonical indices via PDB chain/number/insertion code
    canon = {(r.get_parent().id, r.id[1], r.id[2].strip()): i for i, r in enumerate(cx.residues)}
    info = pose.pdb_info()
    pose2canon = {}
    for p in range(1, pose.total_residue() + 1):
        key = (info.chain(p), info.number(p), info.icode(p).strip())
        if key in canon:
            pose2canon[p] = canon[key]
    if len(pose2canon) != cx.L:
        raise ValueError(f"{cx.name}: only {len(pose2canon)}/{cx.L} residues matched between Rosetta pose and PDB")

    maps = {t: np.zeros((cx.L, cx.L), dtype=np.float32) for t in ROSETTA_TERMS}
    energies = pose.energies()

    def add(p1, p2, emap):
        i, j = pose2canon.get(p1), pose2canon.get(p2)
        if i is None or j is None:
            return
        for t, st in terms:
            e = emap[st] * weights[st]
            if e:
                maps[t][i, j] += e
                maps[t][j, i] += e

    # Short-range two-body terms live on the energy graph edges
    n = pose.total_residue()
    graph = energies.energy_graph()
    for p1 in range(1, n + 1):
        for p2 in range(p1 + 1, n + 1):
            edge = graph.find_energy_edge(p1, p2)
            if edge:
                add(p1, p2, edge.fill_energy_map())

    # Disulfides are a long-range two-body term stored in their own container
    lrec = energies.long_range_container(LongRangeEnergyType.fa_disulfide_energy)
    if lrec is not None and not lrec.empty():
        for p1 in range(1, n + 1):
            if not lrec.any_upper_neighbors_for_residue(p1):
                continue
            it = lrec.const_upper_neighbor_iterator_begin(p1)
            end = lrec.const_upper_neighbor_iterator_end(p1)
            while it != end:
                emap = EMapVector()
                it.retrieve_energy(emap)
                add(p1, it.upper_neighbor_id(), emap)
                it.pre_increment()

    return {f"phys_{t}": m for t, m in maps.items()}


##############################
###### Build
##############################

def build_one(struct_dir, channels, settings, out_dir, pdb_name):
    cx = Complex(struct_dir, settings, pdb_name)
    names = [c[0] for c in channels]
    computed = {}
    for fn in dict.fromkeys(FEATURES[n] for n in names):  # each function once, in channel order
        computed.update(fn(cx))

    x = np.stack([computed[n] for n in names]).astype(np.float32)
    n_nan = int(np.isnan(x).sum())
    if n_nan:
        bad = [n for n, ch in zip(names, x) if np.isnan(ch).any()]
        warnings.warn(f"{cx.name}: {n_nan} NaN values (missing atoms?) in {bad}, set to 0")
        x = np.nan_to_num(x, nan=0.0)
    n_big = int((np.abs(x) > FP16_MAX).sum())
    if n_big:
        warnings.warn(f"{cx.name}: {n_big} values exceed the float16 range and were clipped")
        x = np.clip(x, -FP16_MAX, FP16_MAX)

    def region_seq(name):
        m = cx.regions.get(name)
        return "".join(np.array(list(cx.seq))[m]) if m is not None else ""

    out_file = os.path.join(out_dir, f"{cx.name}.npz")
    np.savez_compressed(
        out_file,
        x=x.astype(np.float16),
        channels=np.array(names),
        groups=np.array([c[1] for c in channels]),
        name=np.array(cx.name),
        seq=np.array(cx.seq),
        chain=cx.chain,
        anchor_index=cx.anchor_index,
        mhc_class=cx.mhc_class,
        peptide=np.array(region_seq("Peptide") or "".join(np.array(list(cx.seq))[cx.chain == 2])),
        cdr3a=np.array(region_seq("CDR3a")),
        cdr3b=np.array(region_seq("CDR3b")),
        trim=np.array([cx.trim[c] for c in CHAINS]),  # kept [start, end) per chain, original offsets
        flags=np.array(cx.flags, dtype=str),
    )
    return index_row(out_file)


def index_row(npz_file):
    """One Input_structures.csv row describing a built .npz (reads only the small arrays,
    not the map, except for files built before peptide/CDR3 were stored)."""
    with np.load(npz_file) as d:
        seq = str(d["seq"])
        chain = d["chain"]
        anchor = int(d["anchor_index"])
        pep_mask = chain == CHAINS.index("Peptide")
        if "cdr3a" in d:
            name, mhc_class = str(d["name"]), int(d["mhc_class"])
            peptide, cdr3a, cdr3b = str(d["peptide"]), str(d["cdr3a"]), str(d["cdr3b"])
        else:
            name = os.path.basename(npz_file)[:-4]
            mhc_class = 2 if (chain == CHAINS.index("MHCb")).any() else 1
            peptide = "".join(np.array(list(seq))[pep_mask])
            region = np.diag(d["x"][list(d["channels"]).index("pos_region_CDR3")]) == 2
            cdr3a = "".join(np.array(list(seq))[region & (chain == CHAINS.index("TCRa"))])
            cdr3b = "".join(np.array(list(seq))[region & (chain == CHAINS.index("TCRb"))])
        n_channels = len(d["channels"])
        flags = ";".join(d["flags"].tolist())
    row = {"File": os.path.basename(npz_file), "Name": name, "MHC_class": mhc_class,
           "Peptide": peptide, "CDR3a": cdr3a, "CDR3b": cdr3b, "Length": len(seq),
           "Anchor_index": anchor, "TCR_side": len(seq) - anchor - 1, "Channels": n_channels}
    for c, n in enumerate(CHAINS):
        row[f"Len_{n}"] = int((chain == c).sum())
    row["Flags"] = flags
    return row


def write_index(out_dir, workers=1):
    """Rebuild <out_dir>/Input_structures.csv from every .npz in out_dir, so it always
    matches the files on disk. File paths are relative to the CSV."""
    files = sorted(os.path.join(out_dir, f) for f in os.listdir(out_dir) if f.endswith(".npz"))
    if workers > 1:
        with Pool(workers) as pool:
            rows = pool.map(index_row, files, chunksize=16)
    else:
        rows = [index_row(f) for f in files]
    df = pd.DataFrame(rows)
    out = os.path.join(out_dir, "Input_structures.csv")
    df.to_csv(out, index=False)
    return df, out


def _build_worker(job):
    struct_dir, channels, settings, out_dir, pdb_name = job
    name = os.path.basename(struct_dir)
    try:
        return build_one(struct_dir, channels, settings, out_dir, pdb_name), None
    except SkipComplex as e:
        return None, (name, "skipped", str(e))
    except Exception:
        return None, (name, "failed", traceback.format_exc())


def main():
    p = argparse.ArgumentParser(description="Build 2D residue-pair input maps for the TCR ResNet")
    p.add_argument("-i", "--in-dir", default=os.path.join(PROJECT_DIR, "data", "structures"),
                   help="directory of structure folders (<name>/<name>.json, <name>/Structure/)")
    p.add_argument("-o", "--out-dir", default=os.path.join(PROJECT_DIR, "data", "input"))
    p.add_argument("--pdb", default="Complex.pdb", help="structure file inside each Structure/ folder")
    p.add_argument("--workers", type=int, default=1)
    p.add_argument("--overwrite", action="store_true", help="rebuild structures that already have an .npz")
    p.add_argument("--names-file", default=None,
                   help="Optional text file of structure names (one per line) to restrict "
                        "building to, instead of every subdirectory in --in-dir")
    p.add_argument("--index-only", action="store_true",
                   help="skip building, just rewrite Input_structures.csv from the .npz files")
    args = p.parse_args()
    if args.index_only:
        df, out = write_index(args.out_dir, args.workers)
        print(f"Indexed {len(df)} structures -> {out}")
        return

    channels, settings = CHANNELS, TRIM
    missing = [c[0] for c in channels if c[0] not in FEATURES]
    if missing:
        raise SystemExit(f"No feature function for channels: {missing}")
    groups = dict.fromkeys(c[1] for c in channels)
    print(f"{len(channels)} channels: "
          + ", ".join(f"{g} {sum(c[1] == g for c in channels)}" for g in groups))

    os.makedirs(args.out_dir, exist_ok=True)
    dirs = sorted(d for d in os.listdir(args.in_dir) if os.path.isdir(os.path.join(args.in_dir, d)))
    if args.names_file is not None:
        with open(args.names_file) as f:
            wanted = [line.strip() for line in f if line.strip()]
        missing = sorted(set(wanted) - set(dirs))
        if missing:
            print(f"Warning: {len(missing)} names from {args.names_file} have no subdirectory "
                  f"in {args.in_dir}: {missing[:10]}")
        dirs = [d for d in dirs if d in wanted]
        print(f"Restricted to {len(dirs)} of {len(wanted)} names from {args.names_file}")
    if not args.overwrite:
        dirs = [d for d in dirs if not os.path.exists(os.path.join(args.out_dir, f"{d}.npz"))]
    jobs = [(os.path.join(args.in_dir, d), channels, settings, args.out_dir, args.pdb) for d in dirs]
    print(f"Building {len(jobs)} structures from {args.in_dir}")

    rows, failed, skipped = [], [], []
    if args.workers > 1:
        with Pool(args.workers) as pool:
            results = list(tqdm(pool.imap_unordered(_build_worker, jobs), total=len(jobs)))
    else:
        results = [_build_worker(j) for j in tqdm(jobs)]
    for row, err in results:
        if err is None:
            rows.append(row)
        elif err[1] == "skipped":
            skipped.append(err)
        else:
            failed.append(err)

    # Write problems to a log rather than the terminal; there can be hundreds
    log = os.path.join(args.out_dir, "build_log.txt")
    with open(log, "w") as f:
        for name, kind, msg in sorted(skipped + failed):
            f.write(f"=== {name}: {kind}\n{msg.rstrip()}\n")
    if failed:
        reasons = pd.Series([e[2].strip().splitlines()[-1][:120] for e in failed]).value_counts()
        print("Failure reasons:\n" + reasons.to_string())

    # Check everything fits the canvas
    if rows:
        df_new = pd.DataFrame(rows)
        size, anchor = CANVAS["canvas_size"], CANVAS["anchor"]
        over = df_new[(df_new["Anchor_index"] > anchor) | (df_new["TCR_side"] > size - anchor - 1)]
        print(f"pMHC side max {df_new['Anchor_index'].max() + 1} (canvas allows {anchor + 1}), "
              f"TCR side max {df_new['TCR_side'].max()} (allows {size - anchor - 1})")
        if len(over):
            print(f"WARNING: {len(over)} structures do not fit the canvas: {over['Name'].tolist()[:10]}")

    df, index = write_index(args.out_dir, args.workers)
    print(f"Done: {len(rows)} built, {len(skipped)} skipped, {len(failed)} failed. "
          f"{len(df)} structures in {index}, problems: {log}")

if __name__ == "__main__":
    main()

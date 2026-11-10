#!/usr/bin/env python3

# This file will be used to setup the data from the D://Lab/TCRBinding/ directories.
# We want to go through all the specified subdirectories, run the setup function on each structure and save the results in the ./dat directory and create a tsv which contains originating information


import argparse
import difflib
import json
import os
import shutil
import traceback
import pandas as pd
import numpy as np
from tqdm import tqdm

from tempfile import gettempdir
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import biotite
import biotite.database.entrez as entrez
import biotite.database.rcsb as rcsb
import biotite.sequence as bseq
import biotite.sequence.graphics as graphics
import biotite.sequence.io.genbank as gb
import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx
from immunum import Annotator


import MDAnalysis as mda
from MDAnalysis.analysis.dssp import DSSP
import MDAnalysis.transformations as tr

from scipy.spatial.transform import Rotation
from scipy.spatial import distance
from scipy.optimize import linear_sum_assignment

import plotly.graph_objects as go
import plotly.express as px
import random
import matplotlib.pyplot as plt
from plotly.tools import mpl_to_plotly
from Bio.PDB import PDBParser, Select, Structure, Model
from Bio.PDB import PDBIO
import warnings
from Bio import BiopythonWarning


from pdbfixer import PDBFixer
from openmm.app import PDBFile

# Suppress all Biopython-specific warnings globally
warnings.simplefilter("ignore", BiopythonWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)


class ChainSelect(Select):
    def __init__(self, accepted_chains):
        self.accepted_chains = accepted_chains

    # def accept_model(self, model):
    #         return 1  # Required: accept or reject models

    def accept_chain(self, chain):
            # chain.get_id() returns the chain identifier (e.g., 'A', 'B')
            if chain.get_id() in self.accepted_chains:
                return 1
            return 0
    
    # def accept_residue(self, residue):
    #     return 1  # Required: accept or reject residues

    # def accept_atom(self, atom):
    #     return 1  # Required: accept or reject atoms

# class SelectChains(Select):
#     def __init__(self, chain_ids):
#         self.chain_ids = chain_ids

#     def accept_chain(self, chain):
#         return chain.get_id() in self.chain_ids
    

    

######################################
########   Functions #################
#######################################


def setupData(input, output, seqfile, prefix=None):
    # This function takes in the input directory which contains a bunch of .pdb files and then runs the setup function on each one, saving the output to the output directory, appending the prefix if not None
    # Returns the list of structure names successfully processed (i.e. the names of the
    # ./<output>/<name>/ folders this call wrote), so a caller can tell which names in a
    # shared/accumulating output directory are new from this call rather than reprocessing
    # (or re-extracting features/inputs for) every structure already there.
    files = os.listdir(input)
    seqdf = pd.read_csv(seqfile, sep=",")
    processed = []
    for file in tqdm(files):
        try:
            # We need to setup a directory for each file if it doesn't exist yet and then create a Structure directory in each one and copy the pdb file there
            name = os.path.splitext(file)[0]
            # print(name)
            # remove everything after the first "." in the name
            name = name.split(".")[0]
            if prefix is not None:
                name = f"{prefix}_{name}"
            if not os.path.exists(os.path.join(output, name)):
                os.makedirs(os.path.join(output, name))
            structureDir = os.path.join(output, name, "Structure")
            os.makedirs(structureDir, exist_ok=True)
            #First copy the original pdb file to the structure directory
            shutil.copy2(os.path.join(input, file), os.path.join(structureDir, f"{name}.pdb"))

            # Now we need to run the json setup function so that we construct and save a JSON in each folder
            # First, let's pull from the sequence file by name and get the row for that structure
            row = seqdf[seqdf["ID"] == name].fillna("")
            # Add hydrogens and normalize the structures
            add_hydrogens(
                os.path.join(input, file), os.path.join(structureDir, f"hydro.pdb")
            )
            # Now save use the renameChains function to rename the chains to match the sequence and save to the structure directory
            renameChains(
                os.path.join(structureDir, f"hydro.pdb"),  # os.path.join(input, file),
                os.path.join(structureDir, f"Complex.pdb"),
                row,
            )
            #Also split up the different components and save to the same folder
            splitComponents(os.path.join(structureDir, f"Complex.pdb"), structureDir)
            #We should have 4 pdb files, the named pdb (protonated), complex pdb (protonated and renamed), TCR.pdb, pMHC.pdb
            # Now we need to turn every column but the first into a dictionary so that we can pass it to the setupJSON function
            sequence = {k: {"Seq": v} for k, v in row.iloc[0, 1:].to_dict().items()}
            jsdat, u = setupJSON(sequence, os.path.join(structureDir, f"Complex.pdb"),os.path.join(structureDir, f"{name}.pdb"))
            # Now write json to the output directory
            with open(os.path.join(output, name, f"{name}.json"), "w") as f:
                json.dump(jsdat, f, cls=NumpyEncoder)
            # Plot the transformed component COMs to output/plots/
            plotDir = os.path.join(output, name, "Plots")
            os.makedirs(plotDir, exist_ok=True)
            plotComponents(jsdat, u, plotDir)
            processed.append(name)
        except Exception as e:
            print(f"Error processing {file}. Skipping. Error: {e}")
            #traceback
            # traceback.print_exc()
    return processed


def splitComponents(infile, outdir):
    class ChainSelect:
        def __init__(self, chain_ids):
            self.chain_ids = chain_ids
        def accept_chain(self, chain):
            return chain.id in self.chain_ids
    parser = PDBParser()
    structure = parser.get_structure("structure", infile)
    # We need to split the structure into a TCR (chain D,E) and pMHC (chain A,B,C) and save each one to the outdir
    tcr_chains = ["D", "E"]
    pmhc_chains = ["A", "B", "C"]

    for comp_name, chains in [("TCR", tcr_chains), ("pMHC", pmhc_chains)]:
        temp = subset_chains(structure.copy(), chains)
        io = PDBIO()
        io.set_structure(temp)
        io.save(os.path.join(outdir, f"{comp_name}.pdb"))


def renameChains(infile, outfile, row):
    parser = PDBParser()
    row = row.iloc[0]
    structure = parser.get_structure(row.iloc[0], infile)

    # columns 1:6 are the reference sequences for the 5 possible chain roles
    # (MHCa, MHCb, Peptide, TCRa, TCRb -> A, B, C, D, E); "" means that role
    # isn't present in this structure.
    role_seqs = {
        name: seq for name, seq in zip(["A", "B", "C", "D", "E"], row.iloc[1:6]) if seq != ""
    }
    namekey = matchChainsBySequence(structure, role_seqs)

    normalize_structure(structure, namekey)  # Normalize using the computed names
    io = PDBIO()
    io.set_structure(structure)
    io.save(outfile)


def chainSequence(chain):
    # One-letter sequence of a chain's standard (non-hetero) residues, for matching
    # against a reference sequence. "X" stands in for any residue seq_dict doesn't
    # recognize (e.g. a modified residue).
    return "".join(
        seq_dict.get(residue.get_resname(), "X")
        for residue in chain
        if residue.id[0] == " "
    )


def matchChainsBySequence(structure, role_seqs):
    # Matches each chain in the structure to whichever reference sequence in role_seqs
    # (role name -> sequence) it most closely resembles, instead of relying on the raw
    # PDB chain ids/order - those aren't consistent across structures and don't reliably
    # line up with the intended MHCa/MHCb/Peptide/TCRa/TCRb roles. Similarity is scored
    # with difflib's SequenceMatcher ratio, and scipy's Hungarian algorithm
    # (linear_sum_assignment) finds the globally best one-to-one pairing across all
    # chains at once, rather than greedily matching chain by chain, which could lock in
    # a worse pairing for a chain considered later. Returns {original chain id: role name}.
    chains = [c for c in structure[0] if chainSequence(c)]
    chain_seqs = [chainSequence(c) for c in chains]
    role_names = list(role_seqs.keys())
    role_sequences = [role_seqs[name] for name in role_names]

    cost = np.array(
        [
            [1.0 - difflib.SequenceMatcher(None, cseq, rseq).ratio() for rseq in role_sequences]
            for cseq in chain_seqs
        ]
    )
    chain_idx, role_idx = linear_sum_assignment(cost)
    namekey = {chains[i].id: role_names[j] for i, j in zip(chain_idx, role_idx)}

    if len(namekey) != len(role_seqs):
        print("Output structure does not match rename chains")

    return namekey


def normalize_structure(structure, namekey):
    # Renames each chain in the structure's first model per namekey (original chain id -> new chain id),
    # then renumbers each chain's residues sequentially from 1 to len(chain). Both passes go through a
    # temporary id first, since Biopython's Model/Chain key their children by id and won't tolerate two
    # children transiently sharing one (e.g. renaming A->D while D->A, or shifting residue numbers into a
    # range another residue in the chain still occupies).
    model = structure[0]

    rename_pairs = [(chain, namekey[chain.id]) for chain in model]
    if len(rename_pairs) != len(namekey):
        print("Output structure does not match rename chains")
    for i, (chain, _) in enumerate(rename_pairs):
        chain.id = f"_tmp{i}"
    for chain, new_id in rename_pairs:
        chain.id = new_id

    for chain in model:
        residues = list(chain)
        for idx, residue in enumerate(residues):
            residue.id = (residue.id[0], -(idx + 1), residue.id[2])
        for idx, residue in enumerate(residues, start=1):
            residue.id = (residue.id[0], idx, " ")

    return structure


def add_hydrogens(file, outfile):
    # 1. Load the PDB file
    fixer = PDBFixer(filename=file)

    # 2. Find and add missing structural heavy atoms first (required step)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()

    # 3. Add hydrogen atoms at a specified pH (default is 7.0)
    fixer.addMissingHydrogens(pH=7.0)

    # 4. Save the finalized structure
    with open(outfile, "w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)


class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, np.generic):
            return obj.item()
        return super(NumpyEncoder, self).default(obj)

def subset_chains(structure, chains):
    chain_ids = [chain.id for chain in structure.get_chains()]
    model = structure[0]
    struct = Structure.Structure("1")
    mod = Model.Model("1")
    struct.add(mod)
    for c in chains:
        if c in chain_ids:
            mod.add(model[c])
    return struct



########################################################
#############Secondary Structure#########################
###########################################################


def drawHelix(fig, begin, end, bbox, name):
    # Approx. 1 turn per 3.6 residues to resemble natural helix
    n_turns = np.ceil((end - begin + 1) / 3.6)
    x_val = np.linspace(0, n_turns * 2 * np.pi, 1000)
    # Curve ranges from 0.3 to 0.7
    y_val = (-0.4 * np.sin(x_val) + 1) / 2

    # Transform values for correct location in feature map
    x_val *= bbox["width"] / (n_turns * 2 * np.pi)
    x_val += bbox["x0"]
    y_val *= bbox["height"]
    y_val += bbox["y0"]
    fig.add_shape(
        type="rect",
        x0=bbox["x0"],
        y0=bbox["y0"],
        x1=bbox["x0"] + bbox["width"],
        y1=bbox["y0"] + bbox["height"],  # Define corners
        line=dict(color="white"),
        fillcolor="white",
        layer="below",
    )
    fig.add_trace(
        go.Scatter(
            x=x_val,
            y=y_val,
            mode="lines",
            line_width=2.5,
            line_color="green",
            name=name,
        )
    )


def drawSheet(fig, begin, end, bbox, name):
    fig.add_trace(
        go.Scatter(
            x=[
                bbox["x0"],
                bbox["x0"],
                bbox["x0"] + (bbox["width"] * 0.6),
                bbox["x0"] + (bbox["width"] * 0.6),
            ],
            y=[8, 12, 12, 8],
            fill="toself",
            fillcolor="green",
            mode="lines",
            name=name,
            line=dict(color="green"),
        )
    )
    # fig.add_shape(
    #     type="rect",
    #     x0=bbox['x0'], y0=8, x1=bbox['x0'] + (bbox['width']*0.5), y1=12, # Define corners
    #     line=dict(color="green"),
    #     fillcolor="green",
    #     layer = "above",
    #     name = name)
    fig.add_trace(
        go.Scatter(
            x=[
                bbox["x0"] + (bbox["width"] * 0.5),
                bbox["x0"] + (bbox["width"] * 0.5),
                bbox["x0"] + bbox["width"],
                bbox["x0"] + (bbox["width"] * 0.5),
            ],
            y=[7, 13, 10, 7],
            fill="toself",
            mode="lines",
            fillcolor="green",
            line_color="green",
            name=name,
        )
    )


def viz_sec(annotation, sse, region):
    # Build a plotly figure where the
    fig = go.Figure()
    h = 20  # (math.ceil((len(sse)/8)/10) * 10)
    l = len(sse) + 10
    fig.update_layout(
        width=800,
        height=200,
        margin=dict(l=20, r=20, t=20, b=20),
        yaxis=dict(rangemode="tozero"),
        yaxis_range=[0, h],
        xaxis_range=[0, l],
        template="plotly_white",
        showlegend=False,
    )
    fig.update_xaxes(showgrid=False)
    fig.update_yaxes(showgrid=False, showticklabels=False)
    # fig.add_hline(y=2, line_width=2,line_color="black")
    fig.add_shape(
        type="rect",
        x0=0,
        y0=(h / 2) - 0.2,
        x1=l,
        y1=(h / 2) + 0.2,  # Define corners
        line=dict(color="grey"),
        fillcolor="grey",
        layer="below",
    )
    ns, nh = (0, 0)
    for a in annotation:
        bbox = {
            "x0": a["location"][0],
            "y0": (h / 2) - 7,
            "width": a["location"][1] - a["location"][0],
            "height": 14,
        }
        if a["Type"] == "helix":
            # nh = nh + 1
            # name = f'helix {nh}'
            drawHelix(fig, a["location"][0], a["location"][1], bbox, a["Name"])
        if a["Type"] == "sheet":
            # ns = ns + 1
            # name = f'sheet {ns}'
            drawSheet(fig, a["location"][0], a["location"][1], bbox, a["Name"])

    # add region annotation (line and label)
    palette = px.colors.qualitative.Plotly
    for r in region:
        if r[0] == 0 and r[1] == 0:
            continue
        col = random.choice(palette)
        if r in ["TCRca", "TCRcb", "TCRva", "TCRvb", "MHCha", "MHChb"]:
            yline = 5
            ylabel = 2
        else:
            yline = 15
            ylabel = 18
        fig.add_shape(
            type="line",
            x0=region[r][0],
            y0=yline,
            x1=region[r][1],
            y1=yline,
            line=dict(color=col, width=4),
        )
        dist = region[r][1] - region[r][0]
        fig.add_annotation(
            x=region[r][0] + dist / 2,
            y=ylabel,
            text=r,
            showarrow=False,
            font=dict(color=col, size=10),
        )
    fig.update_layout()
    return fig


def get_annotation(sse, first_id, chain="all"):
    def _add_sec_str(annotation, first, last, str_type):
        # if str_type == "a":
        if str_type == "H":
            str_type = "helix"
        # elif str_type == "b":
        elif str_type == "E":
            str_type = "sheet"
        else:
            # coil
            return
        feature = bseq.Feature(
            "SecStr", [bseq.Location(first, last)], {"sec_str_type": str_type}
        )
        annotation.add_feature(feature)

    # Find the intervals for each secondary structure element
    # and add to annotation
    annotation = bseq.Annotation()
    curr_sse = None
    curr_start = None
    for i in range(len(sse)):  # walk through the secondary structure and add a

        if curr_start is None:
            curr = sse[i]
            err = 0
            curr_start = i

        elif i > 1:
            last = [sse[i - 1]]
            if curr not in last:
                err = err + 1
        if err > 1:
            _add_sec_str(annotation, curr_start + first_id, i - 1 + first_id, curr)
            curr_start = i
            curr = sse[i]
            err = 0

    # Add last secondary structure element to annotation
    _add_sec_str(annotation, curr_start + first_id, i - 1 + first_id, curr_sse)

    return annotation


###################################################
def annotation_to_dict(annotation):
    annotation_dict = []
    for feature in annotation:
        feature_data = {
            "Type": feature.qual["sec_str_type"],
            "location": (
                feature.get_location_range()[0].item(),
                feature.get_location_range()[1].item(),
            ),
        }
        annotation_dict.append(feature_data)
    # We want to also sort the this annotation dict object by the lowest feature so that we can name them in order
    annotation_dict = sorted(annotation_dict, key=lambda d: d["location"][0])
    # loop through and add a name key based on sheet/helix
    h, s = 0, 0
    for i in range(len(annotation_dict)):
        if annotation_dict[i]["Type"] == "helix":
            annotation_dict[i]["Name"] = f"Helix {h+1}"
            h = h + 1
        if annotation_dict[i]["Type"] == "sheet":
            annotation_dict[i]["Name"] = f"Sheet {s+1}"
            s = s + 1
    return annotation_dict


##################################################
# Create a function that will turn the dict back into an annotation object for plotting
def dict_to_annotation(annotation_dict):
    annotation = b.Annotation()
    for element in annotation_dict:
        feature = bseq.Feature(
            "SecStr",
            [bseq.Location(element["location"][0], element["location"][1])],
            {"sec_str_type": element["Type"]},
        )
        annotation.add_feature(feature)
    return annotation


def getSecondaryStructures(u):
    cum_score = []
    chains = np.unique(u.atoms.chainIDs)
    if "A" in chains and "B" in chains:
        cls = 2
    else:
        cls = 1
    sec_dict = {}
    chain_dict = {"A": "MHCa", "B": "MHCb", "C": "Peptide", "D": "TCRa", "E": "TCRb"}
    for c in u.segments:
        name = chain_dict[np.unique(c.atoms.chainIDs)[0]]

        dssp_analysis = DSSP(c.atoms).run()
        sse = dssp_analysis.results.dssp[0]
        annotation = get_annotation(sse, 1)
        if name in ["Peptide","MHCb"]:
            score = {
                "n_helix": 0,
                "n_sheet": 0,
                "exp_helix": 0,
                "exp_sheet": 0,
                "score": 0,
            }
        else:
            score = scoreSecondary(annotation, sse, name, cls)

        cum_score.append(score["score"])
        sec_dict[name] = {
            "Name": name,
            "Annotation": annotation_to_dict(annotation),
            "Structure": sse,
            "Score": score,
        }
        # sec_dict[name] = {'Annotation': annotation, 'Structure': sse, 'Score': score}
        # visualize_secondary_structure(annotation,sse,title=name) #Don't visualize here
    sec_dict["Cumulative score"] = np.sum(cum_score)/3
    return sec_dict


# Score the secondary structure
def scoreSecondary(annotation, sse, name, cls):
    # if cls == 1:
    #     exp = {"MHCa": (3, 8), "TCRa": (3, 16), "TCRb": (4, 20)}
    # else:
    #     exp = {"MHCa": (1, 4), "MHCb": (2, 4), "TCRa": (3, 16), "TCRb": (4, 20)}

    score = 0
    length = len(sse)
    r_helix = (sse == "H").sum()
    r_sheet = (sse == "E").sum()
    r_inter = (sse == "-").sum()
    p_helix = r_helix / length * 100 if length > 0 else 0
    p_sheet = r_sheet / length * 100 if length > 0 else 0
    p_inter = r_inter / length * 100 if length > 0 else 0
    # print(f"name: {name}, r_helix: {r_helix}, r_sheet: {r_sheet}, r_inter: {r_inter}, length: {length}, %_helix: {r_helix/length*100:.2f}, %_sheet: {r_sheet/length*100:.2f}, %_inter: {r_inter/length*100:.2f}")
    # n_helix = 0
    # n_sheet = 0
    # # exp_h = exp[name][0]
    # # exp_s = exp[name][1]
    # for a in annotation.get_features():
    #     feat = a.qual["sec_str_type"]
    #     if feat == "helix":
    #         n_helix = n_helix + 1
    #     if feat == "sheet":
    #         n_sheet = n_sheet + 1
    # print(f"name: {name}, length: {length}, n_helix: {n_helix}, n_sheet: {n_sheet}, %_helix: {n_helix/length*100:.2f}, %_sheet: {n_sheet/length*100:.2f}")
    if name == "MHCa":
        scoring = [(27<= p_helix <= 45), (25 <= p_sheet <= 35), (p_inter <= 70)]
    elif name == "TCRa":
        scoring = [(p_helix <= 30), (5 <= p_sheet <= 55), (p_inter <= 70)]
    elif name == "TCRb":
        scoring = [(3 <= p_helix <= 30), (30 <= p_sheet <= 55), (p_inter <= 70)]
    else: 
        return {}
    score = {
        # "n_helix": n_helix,
        # "n_sheet": n_sheet,
        # "exp_helix": exp_h,
        # "exp_sheet": exp_s,
        "r_helix": r_helix,
        "r_sheet": r_sheet,
        "r_inter": r_inter,
        "length": length,
        "score": sum(scoring)/3,
    }
    return score


def extract_helix_stretches(features):
    starts = []
    stops = []
    current_start = None
    current_stop = None
    for item in features:
        feature_type = item["Type"]
        start, stop = item["location"]
        if feature_type == "helix":
            # Start a new helix stretch if not already in one
            if current_start is None:
                current_start = start
            current_stop = stop
        elif feature_type == "sheet":
            # End current helix stretch when interrupted by a sheet
            if current_start is not None:
                starts.append(current_start)
                stops.append(current_stop)
                current_start = None
                current_stop = None
    # Catch a helix stretch that continues to the end
    if current_start is not None:
        starts.append(current_start)
        stops.append(current_stop)
    return starts, stops


##################################################################
#####################Regions######################################
#################################################################

# Predict regions (based on the expected conserved secondary structure or conserved resid)
def predictRegions(u, secondary, sequences):
    regions = {}
    chains = np.unique(u.atoms.chainIDs)
    chain_dict = {"A": "MHCa", "B": "MHCb", "C": "Peptide", "D": "TCRa", "E": "TCRb"}
    # First lets do the TCR since we can always expect there to be a TCRab

    # First the MHC
    # try:
    mhca = secondary["MHCa"]["Structure"]
    mhca_helix_indices = [i for i, x in enumerate(mhca) if x == "H"]
    if "B" in chains:
        mhcb = secondary["MHCb"]["Structure"]
        #We want to find the min and max extents of helix ([mhca] == "H") within the defined ranges so first need to get the indices of all elements that are helices
        mhcb_helix_indices = [i for i, x in enumerate(mhcb) if x == "H"]
        #Now get a subset of indices between a min and max value
        mhcb_helix = [i for i in mhcb_helix_indices if i >= 40 and i <= 85]
        mhca_helix = [i for i in mhca_helix_indices if i >= 40 and i <= 85]
        # helix = [x for x in mhca if x["Type"] == "helix"]
        # MHCha = (helix[0]["location"][0], helix[-1]["location"][1])
        # helix = [x for x in mhca if x["Type"] == "helix"]
        # MHChb = (helix[0]["location"][0], helix[-1]["location"][1])
        regions["MHCa"] = {"MHCha": (min(mhca_helix), max(mhca_helix))}
        regions["MHCb"] = {"MHChb": (min(mhcb_helix), max(mhcb_helix))}
    else:
        #
        # starts, stops = extract_helix_stretches(mhca)
        # MHCha = (starts[0], stops[0])
        # MHChb = (starts[-1], stops[-1])
        # helix = [x for x in mhca if x['Type'] == 'helix']
        # MHCha = (helix[0]['location'][0],helix[-3]['location'][1])
        # MHChb = (helix[-2]['location'][0],helix[-1]['location'][1])
        mhca_helix = [i for i in mhca_helix_indices if i >= 40 and i <= 85]
        mhcb_helix = [i for i in mhca_helix_indices if i >= 138 and i <= 175]
        regions["MHCa"] = {"MHCha": (min(mhca_helix), max(mhca_helix)), "MHChb": (min(mhcb_helix), max(mhcb_helix))}
    # except:
    #     regions["MHCa"] = {"MHCha": (0, 0), "MHChb": (0, 0)}    

    # First do the peptide since it is just chain c
    peptide = u.select_atoms("segid C")
    peptide_length = len(peptide.residues)
    regions["Peptide"] = {"Peptide": (1, peptide_length)}
    # Do CDR1-3 based on IMGT annotation
    try:
        TCRa_anno = sequences["TCRa"]["Consensus"]
        tcrva_idx = [i for i, x in enumerate(TCRa_anno["IMGT"]) if x != -1]
        tcrca_idx = [i for i, x in enumerate(TCRa_anno["IMGT"]) if x == -1]
        TCRva = (min(tcrva_idx), max(tcrva_idx))
        if len(tcrca_idx) > 0:
            TCRca = (min(tcrca_idx), max(tcrca_idx))
        else:
            TCRca = (0, 0)

        cdr1a_idx = [i for i, x in enumerate(TCRa_anno["Turn"]) if x == "cdr1"]
        cdr2a_idx = [i for i, x in enumerate(TCRa_anno["Turn"]) if x == "cdr2"]
        cdr3a_idx = [i for i, x in enumerate(TCRa_anno["Turn"]) if x == "cdr3"]
        CDR1a = (min(cdr1a_idx), max(cdr1a_idx))
        CDR2a = (min(cdr2a_idx), max(cdr2a_idx))
        CDR3a = (min(cdr3a_idx), max(cdr3a_idx))
        regions["TCRa"] = {
            "TCRva": TCRva,
            "TCRca": TCRca,
            "CDR1a": CDR1a,
            "CDR2a": CDR2a,
            "CDR3a": CDR3a,
        }
    except:
        regions["TCRa"] = {
            "TCRva": (0, 0),
            "TCRca": (0, 0),
            "CDR1a": (0, 0),
            "CDR2a": (0, 0),
            "CDR3a": (0, 0),
        }

    try:
        TCRb_anno = sequences["TCRb"]["Consensus"]
        tcrvb_idx = [i for i, x in enumerate(TCRb_anno["IMGT"]) if x != -1]
        tcrcb_idx = [i for i, x in enumerate(TCRb_anno["IMGT"]) if x == -1]
        TCRvb = (min(tcrvb_idx), max(tcrvb_idx))
        if len(tcrcb_idx) > 0:
            TCRcb = (min(tcrcb_idx), max(tcrcb_idx))
        else:
            TCRcb = (0, 0)

        cdr1b_idx = [i for i, x in enumerate(TCRb_anno["Turn"]) if x == "cdr1"]
        cdr2b_idx = [i for i, x in enumerate(TCRb_anno["Turn"]) if x == "cdr2"]
        cdr3b_idx = [i for i, x in enumerate(TCRb_anno["Turn"]) if x == "cdr3"]
        CDR1b = (min(cdr1b_idx), max(cdr1b_idx))
        CDR2b = (min(cdr2b_idx), max(cdr2b_idx))
        CDR3b = (min(cdr3b_idx), max(cdr3b_idx))
        regions["TCRb"] = {
            "TCRvb": TCRvb,
            "TCRcb": TCRcb,
            "CDR1b": CDR1b,
            "CDR2b": CDR2b,
            "CDR3b": CDR3b,
        }
    except:
        regions["TCRb"] = {
            "TCRvb": (0, 0),
            "TCRcb": (0, 0),
            "CDR1b": (0, 0),
            "CDR2b": (0, 0),
            "CDR3b": (0, 0),
        }

    return regions


#########################################################
############Components and Bounding Boxes#########################
###########################################################


class BoundingEllipse3D:
    """
    Represents a 3D bounding ellipsoid for a given set of points.
    Calculates the center, axis lengths, rotation matrix, and volume.
    """

    def __init__(self, points):
        # Convert points to a numpy array for vector operations
        P = np.asarray(points)
        if P.shape[1] != 3:
            raise ValueError("Points must be 3-dimensional coordinates.")
        if P.shape[0] < 4:
            raise ValueError("At least 4 points are required to define a 3D shape.")

        # 1. Calculate centroid
        self.center = np.mean(P, axis=0)

        # 2. Calculate the covariance matrix
        centered_points = P - self.center
        covariance_matrix = np.cov(centered_points, rowvar=False)

        # 3. Perform SVD/Eigen decomposition
        # The eigenvalues represent squared semi-axis lengths
        # The eigenvectors represent the principal axes of rotation
        eigenvalues, eigenvectors = np.linalg.eigh(covariance_matrix)

        # Multiply by 3 for scaling (99% coverage roughly fits points inside)
        # For a minimum bounding ellipsoid, use the Khachiyan Algorithm.
        # This approach approximates the bounding shape using PCA.
        self.radii = np.sqrt(eigenvalues * 3.0)
        self.rotation = eigenvectors

        # 4. Calculate the total bounding volume
        # V = (4/3) * pi * r1 * r2 * r3
        self.volume = (4.0 / 3.0) * np.pi * np.prod(self.radii)

    def get_axis_lengths(self):
        """Returns the semi-axis lengths (r1, r2, r3)."""
        return tuple(self.radii)


# Function that will extract the structural components
def getComponents(u, regions):
    chains = np.unique(u.atoms.chainIDs)
    ################################
    #######COM######################
    ################################
    COM = np.mean(u.atoms.positions, axis=0)
    TCR_COM = np.mean(u.select_atoms("segid D E").positions, axis=0)
    TCRa = u.select_atoms("segid D")
    TCRb = u.select_atoms("segid E")
    TCRa_COM = np.mean(TCRa.positions, axis=0)
    TCRb_COM = np.mean(TCRb.positions, axis=0)
    MHC_COM = np.mean(u.select_atoms("segid A B").positions, axis=0)
    MHCa = u.select_atoms("segid A")
    MHCa_COM = np.mean(MHCa.positions, axis=0)
    Pep_COM = np.mean(u.select_atoms("segid C").positions, axis=0)
    COMs = {
        "All": {"COM": COM},
        "TCR": {"COM": TCR_COM},
        "MHC": {"COM": MHC_COM},
        "MHCa": {"COM": MHCa_COM},
        "Peptide": {"COM": Pep_COM},
        "TCRa": {"COM": TCRa_COM},
        "TCRb": {"COM": TCRb_COM},
    }

    if "B" in chains:
        MHCb = u.select_atoms("segid B")
        MHCb_COM = np.mean(MHCb.positions, axis=0)
        COMs["MHCb"] = {"COM": MHCb_COM}
    else:
        MHCb_COM = MHCa_COM
        COMs["MHCb"] = {"COM": np.array([np.nan, np.nan, np.nan])}
    # Instead of how we just did it, we want to take the "region" input which is a dict of the region and subregion and calculate the COM for alpha carbons between those ranges
    chain_dict = {"A": "MHCa", "B": "MHCb", "C": "Peptide", "D": "TCRa", "E": "TCRb"}
    for c in u.atoms.segments:
        name = chain_dict[np.unique(c.atoms.chainIDs)[0]]
        # print(np.unique(c.atoms.chainIDs), len(c.residues))
        for r in regions[name]:
            if r not in ["TCR", "TCRa", "TCRb", "MHC", "MHCa", "MHCb", "Peptide"]:
                res_range = regions[name][r]
                res_com = np.mean(
                    c.atoms.select_atoms(
                        f"resid {res_range[0]}:{res_range[1]}"
                    ).positions,
                    axis=0,
                )
                if r == "MHChb":
                    # print(name, res_range, res_com)
                    pass
                COMs[name][r] = res_com

    ################################
    ##########Segments##############
    ################################
    # Now make a line segment for each region connecting the 1st and last residue CA
    Segments = {"MHCa": {}, "MHCb": {}, "Peptide": {}, "TCRa": {}, "TCRb": {}}
    for c in u.atoms.segments:
        name = chain_dict[np.unique(c.atoms.chainIDs)[0]]
        for r in regions[name]:
            if r not in [
                "TCR",
                "TCRa",
                "TCRb",
                "MHC",
                "MHCa",
                "MHCb",
                "TCRca",
                "TCRcb",
                "TCRva",
                "TCRvb",
            ]:
                res_range = regions[name][r]
                res_seg = c.atoms.select_atoms(
                    f"resid {res_range[0]} {res_range[1]} and name CA"
                ).positions
                Segments[name][r] = res_seg

    ################################
    ############Bounding***********
    #   ################################
    # Now we want to make a bounding box and ellipsoid for each region based on the CA atoms in that region
    bound = calcBounding(u)
    TCRbound = calcBounding(u.select_atoms("segid D E"))
    MHCbound = calcBounding(u.select_atoms("segid A B"))
    TCRabound = calcBounding(u.select_atoms("segid D"))
    TCRbbound = calcBounding(u.select_atoms("segid E"))

    Bounding = {
        "All": bound,
        "TCR": TCRbound,
        "MHC": MHCbound,
        "TCRa": {"All": TCRabound},
        "TCRb": {"All": TCRbbound},
        "MHCa": {},
        "MHCb": {},
        "Peptide": {},
    }
    for c in u.atoms.segments:
        name = chain_dict[np.unique(c.atoms.chainIDs)[0]]
        for r in regions[name]:
            # if r not in ['TCR','TCRa','TCRb','MHC','MHCa','MHCb']:
            res_range = regions[name][r]
            if (res_range[0] - res_range[1]) == 0:
                continue
            bound = calcBounding(
                c.atoms.select_atoms(f"resid {res_range[0]}:{res_range[1]}")
            )
            Bounding[name][r] = bound

    ##################################
    ########Compactness###############
    ##################################
    #use a function to calculate how many residues of a chain are within 10,20,30,40 of the COM
    Compactness = {}
    for c in u.atoms.segments:
        name = chain_dict[np.unique(c.atoms.chainIDs)[0]]
        com = COMs[name]["COM"]
        Compactness[name] = {}
        res_seg = c.atoms.select_atoms(f"segid {c.segid} and name CA").positions
        #Turn the list of CA atom positions into a numpy array for distance from com
        dists = np.linalg.norm(res_seg - com, axis=1)
        for t in [10, 20, 30, 40]:
            within = np.sum(dists <= t)/len(dists)
            Compactness[name][t] = within

    ###################################
    ########Features and score#########
    ###################################
    scores = scoreComponents(u, COMs, Segments, Bounding, Compactness)
    # try:
    #     scores = scoreComponents(u, COMs, Segments, Bounding)
    # except Exception as e:
    #     print(f"Error occurred while scoring components: {e}")
    #     scores = {}
    structure = {
        "COM": COMs,
        "Segments": Segments,
        "Bounding": Bounding,
        "Compactness": Compactness,
        "Scores": scores
    }
    return structure


def calcBounding(u):
    # This function will calculate the bounding box and ellipsoid for the entire structure
    try:
        box = u.atoms.select_atoms("name CA").bbox()
    except:
        box = np.array([[np.nan, np.nan, np.nan], [np.nan, np.nan, np.nan]])
    try:
        sphere = u.atoms.select_atoms("name CA").bsphere()
    except:
        sphere = (np.array([np.nan, np.nan, np.nan]), np.nan)
    # bounding = {'Box': {'Min': box[0],'Max': box[1],'Volume': np.prod(box[1]-box[0])},'Sphere': {'Center': sphere[0],'Radius': sphere[1],'Volume': (4/3)*np.pi*(sphere[1]**3)}}
    bounding = {"Box": box, "Sphere": sphere}
    return bounding


def _direction_alignment_rotation(y=None, x=None):
    # Builds a rotation that aligns y (peptide->TCR COM direction) to the +y axis
    # and x (MHCa->MHCb COM direction) to the +x axis. y is aligned first and
    # exactly, if given. x is then aligned using only the rotational freedom left
    # after y - i.e. a twist about the y-axis - so it never disturbs y's
    # alignment. If y is not given, x alone is aligned exactly instead.
    rot = Rotation.identity()
    primary_target = None

    for vec, target in [(y, np.array([0.0, 1.0, 0.0])), (x, np.array([1.0, 0.0, 0.0]))]:
        if vec is None:
            continue
        vec = np.asarray(vec, dtype=float)
        vec_norm = np.linalg.norm(vec)
        if not np.isfinite(vec_norm) or vec_norm < 1e-8:
            # Unusable direction (NaN/Inf, e.g. from a chain missing in this
            # structure, or a zero-length vector) - treat it as not given.
            continue
        current = rot.apply(vec / vec_norm)

        if primary_target is None:
            step, _ = Rotation.align_vectors([target], [current])
            rot = step * rot
            primary_target = target
        else:
            axis = primary_target
            # Only the component perpendicular to the already-fixed axis can be
            # moved by a rotation about that axis.
            current_perp = current - np.dot(current, axis) * axis
            target_perp = target - np.dot(target, axis) * axis
            norm_c = np.linalg.norm(current_perp)
            norm_t = np.linalg.norm(target_perp)
            if not np.isfinite(norm_c) or norm_c < 1e-8 or norm_t < 1e-8:
                # vec is (anti)parallel to the already-fixed axis - no twist can align it
                continue
            current_perp /= norm_c
            target_perp /= norm_t
            angle = np.arctan2(
                np.dot(np.cross(current_perp, target_perp), axis),
                np.dot(current_perp, target_perp),
            )
            twist = Rotation.from_rotvec(angle * axis)
            rot = twist * rot

    return rot


def transformCOMs(COM, t, y=None, x=None):
    # Translates every point by t, then applies the y/x axis-alignment rotation
    # from _direction_alignment_rotation (see there for details).
    rot = _direction_alignment_rotation(y=y, x=x)

    def transform(value):
        if isinstance(value, dict):
            return {key: transform(v) for key, v in value.items()}
        return rot.apply(np.array(value) + t)

    return transform(COM)


def transformStructure(u, t, y=None, x=None):
    # Same translation/rotation as transformCOMs, applied to every atom position.
    rot = _direction_alignment_rotation(y=y, x=x)
    u.atoms.positions = rot.apply(u.atoms.positions + t)
    return u


def scoreComponents(u, comps, seg, box, compact):
    ###########################################
    ########## Calc basic orientations#########
    ###########################################

    # 1. normalize the structure so that the peptide is at the origin, the TCR is in the positive y direction, and the MHCa->MHCb direction is in the positive x direction
    pep = comps["Peptide"]["COM"]
    tcr = comps["TCR"]["COM"]
    mhcha = comps["MHCa"]["MHCha"]
    try:
        mhchb = comps["MHCa"]["MHChb"]
    except:
        mhchb = comps["MHCb"]["MHChb"]
    xdir = mhchb - mhcha
    transCOM = transformCOMs(comps, pep, y = tcr - pep, x = xdir) # t, y, x

    # 2. determine if the peptide is between the MHC and TCR (The peptide should be between the MHC and TCR in the y direction so the MHC COM should be negative y and the TCR COM should be positive y)
    score = {}
    mhc_com = transCOM["MHC"]["COM"]
    tcr_com = transCOM["TCR"]["COM"]
    pep_com = transCOM["Peptide"]["COM"]
    if mhc_com[1] < pep_com[1] and tcr_com[1] > pep_com[1]:
        score["Peptide position"] = 0
    else:
        score["Peptide position"] = 1
    # 3. determine if the peptide is between the MHC arms (The peptide should be between the MHC arms in the x direction so the MHCa COM should be negative x and the MHCb COM should be positive x or vice versa)
    mhcha_com = transCOM["MHCa"]["MHCha"]
    mhchb_com = (
        transCOM["MHCa"]["MHChb"]
        if "MHChb" in transCOM["MHCa"]
        else transCOM["MHCb"]["MHChb"]
    )
    if mhcha_com[0] < 0 and mhchb_com[0] > 0:
        score["Peptide horizontal"] = 0
    elif mhcha_com[0] > 0 and mhchb_com[0] < 0:
        score["Peptide horizontal"] = 0
    else:
        score["Peptide horizontal"] = 1
    # 4. determine if the MHC is oriented facing the peptide (MHCcom to the average of the MHCa and MHCb COMs should be pointing in the positive y direction/towards peptide)
    MHCorient = (
        np.mean(np.array([mhcha_com, mhchb_com]), axis=0) - transCOM["MHC"]["COM"]
    )
    MHCorient = MHCorient / np.linalg.norm(MHCorient)
    if MHCorient[1] > 0:
        score["MHC orientation"] = 0
    else:
        score["MHC orientation"] = 1
    #Scoring the compactness of the complex using the compact 
    compact_sum = np.sum([
        (compact["MHCa"][20] > 0.7), (compact["MHCa"][30] > 0.7), (compact["MHCa"][40] > 0.7),
        (compact["TCRa"][20] > 0.25), (compact["TCRa"][30] > 0.7), (compact["TCRa"][40] > 0.9),
        (compact["TCRb"][20] > 0.25), (compact["TCRb"][30] > 0.7), (compact["TCRb"][40] > 0.9)
    ])
    score["Compactness"] = compact_sum/9

    return score


#########################################################
###############Sequence################################
########################################################
# Get a dict of sequences from universe, in the future this will be from an input where we have the sequences in the JSON and need to add other information (like the PDB file sequence/IMGT numbers)
seq_dict = {
    "ALA": "A",
    "CYS": "C",
    "ASP": "D",
    "GLU": "E",
    "PHE": "F",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LYS": "K",
    "LEU": "L",
    "MET": "M",
    "ASN": "N",
    "PRO": "P",
    "GLN": "Q",
    "ARG": "R",
    "SER": "S",
    "THR": "T",
    "VAL": "V",
    "TRP": "W",
    "TYR": "Y",
}


def exactMatch(chaine1, chaine2):
    return [c1 != c2 for c1, c2 in zip(chaine1, chaine2)]


def getSequences(sequences, u):
    chain_dict = {"A": "MHCa", "B": "MHCb", "C": "Peptide", "D": "TCRa", "E": "TCRb"}
    for c in u.segments:
        name = chain_dict[np.unique(c.atoms.chainIDs)[0]]
        # print(np.unique(c.atoms.chainIDs)[0])
        seq = "".join(seq_dict.get(res, "X") for res in c.residues.resnames)
        full = sequences[name]["Seq"]
        exact = exactMatch(full, seq)
        sscore = 1-(sum(exact) / len(full))
        if name in ["TCRa", "TCRb"]:
            anno = annotateChain(seq)
            score = scoreConsensus(anno)
        else:
            anno = None
            score = 0
        sequences[name]["PDBseq"] = seq
        sequences[name]["Exact"] = exact
        sequences[name]["seqscore"] = sscore
        sequences[name]["Consensus"] = anno
        sequences[name]["ConsensusScore"] = score

    return sequences


# Numbering and lettering based on:
# https://www.imgt.org/IMGTindex/numbering.php
# https://www.sciencedirect.com/science/article/pii/S0145305X02000393?via=ihub
def annotateChain(sequence):
    numbering = {}
    numbers = [-1] * len(
        sequence
    )  # First we will make a dummy list and assign the IMGT number based on sequence
    turn = ["none"] * len(sequence)
    # print(sequence)
    annotator = Annotator(chains=["A", "B", "D", "G"], scheme="imgt")
    anno = annotator.number(sequence)
    seg = annotator.segment(sequence)
    # print(anno)
    if anno.numbering == None:
        return None

    for pos in anno.numbering:
        pos = "".join(filter(str.isdigit, pos))
        for s in range(len(sequence)):
            if anno.numbering[pos] == sequence[s] and numbers[s] == -1:
                numbers[s] = int(pos)
                break
            elif numbers[s] == -1:
                numbers[s] = 0
        # numbering[pos] = anno.numbering[pos]
    lab = "none"
    for s in range(len(numbers)):
        if numbers[s] > 0:
            lab = returnLabel(int(numbers[s]))
        if numbers[s] < 0:
            lab = "none"
        turn[s] = lab
    numbering["IMGT"] = numbers
    numbering["Confidence"] = anno.confidence
    numbering["Segments"] = seg.as_dict()
    numbering["Turn"] = turn
    return numbering


def returnLabel(pos):
    if pos <= 15:
        return "A"
    if pos <= 26:
        return "B"
    if pos <= 38:
        return "cdr1"
    if pos <= 46:
        return "C"
    if pos <= 56:
        return "C+"
    if pos <= 66:
        return "cdr2"
    if pos <= 74:
        return "C++"
    if pos <= 84:
        return "D"
    if pos <= 96:
        return "E"
    if pos <= 104:
        return "F"
    if pos <= 117:
        return "cdr3"
    return "G"


def scoreConsensus(anno):
    if anno == None:
        return 0
    score = 0
    turns = ["A", "B", "C", "D", "E", "G", "cdr1", "cdr2", "cdr3"]
    for t in turns:
        if t not in anno["Turn"]:
            score = score + 1
    return 1 - (score / len(turns))


################################################################################
#########pLDDT extraction
##############################################################################


def getpLDDT(u, region,raw):
    pscore = {}
    f = mda.Universe(raw, dt=1)
    chains = np.unique(u.atoms.chainIDs)
    chain_dict = {"A": "MHCa", "B": "MHCb", "C": "Peptide", "D": "TCRa", "E": "TCRb"}
    for c in chains:
        name = chain_dict[c]
        chain_atoms = f.select_atoms(f"chainID {c}")
        plddt = np.array(chain_atoms.tempfactors)
        avg = np.mean(plddt)
        regions = region[name]
        pscore[name] = {"Average": avg}
        for r in regions:
            pscore[name][r] = np.mean(plddt[regions[r][0] : regions[r][1]])
    return pscore


# score the extracted sequences
def extract_nested_values(data_dict, target_key):
    results = []
    # Iterating through keys directly
    for key in data_dict:
        nested_dict = data_dict[key]
        # Ensure the value is a dict before trying to access a key
        if isinstance(nested_dict, dict) and target_key in nested_dict:
            results.append(nested_dict[target_key])
    return np.array(results)


####################################################
############    Scores #############################
##################################################


def buildScores(seq, sec, reg, coms, plddt):
    #With this function we will add a new key to the dict that compiles the scored factors as to the structural validity
    # The number for each will be a 1 if perfectly fit and 0 if not perfect
    score = {}
    #Calculate sequence score each chain has a seqscore which is a number 0-1 for the amount of correct residues
    score["Sequence"] = np.mean(extract_nested_values(seq, "seqscore"))
    #Calculate the consensus score which is the sum of individual consensus scores for each chain
    score["Consensus"] = np.mean([seq["TCRa"]["ConsensusScore"],seq["TCRb"]["ConsensusScore"]])
    #Calc the secondary scores
    score["Secondary"] = sec["Cumulative score"]
    # Calc region scores
    # we want to go through region keys and flatten all the key names into one single list
    detected = []
    for i in reg:
        for k,v in reg[i].items():
            if v[-1] != 0 and v[-1] != np.nan:
                detected.append(k.upper())
    required = [
        "MHCha",
        "MHChb",
        "Peptide",
        "TCRva",
        "TCRvb",
        "cdr1a",
        "cdr2a",
        "cdr3a",
        "cdr1b",
        "cdr2b",
        "cdr3b",
    ]
    found = [r for r in required if r.upper() in detected]
    score["Region"] = len(found)/len(required)
    #Calc component scores (the number of structural components that have a bad score)
    pos_score = (3 - (coms["Scores"]["Peptide position"] + coms["Scores"]["Peptide horizontal"] + coms["Scores"]["MHC orientation"]))/3
    score["Components"] = np.mean(np.array([pos_score,coms["Scores"]["Compactness"]]))
    # Score the pLDDT (Average of variable region pLDDT values)
    values = []
    variable = [
        "MHCha",
        "MHChb",
        "Peptide",
        "cdr1a",
        "cdr2a",
        "cdr3a",
        "cdr1b",
        "cdr2b",
        "cdr3b",
    ]
    for i in plddt:
        for j in plddt[i]:
            if j in variable:
                values.append(plddt[i][j])
    score["pLDDT"] = min(np.mean(values),1)
    #Get a final average score (mean of all individual scores)
    score["Final"] = np.mean(list(score.values()))
    return score




#################################################
##############Function to setup the dict#########
#################################################


def setupJSON(sequence, pdb_file,raw):
    name = os.path.splitext(os.path.basename(pdb_file))[0]
    u = mda.Universe(pdb_file, dt=1)
    seq = getSequences(sequence, u)
    sec = getSecondaryStructures(u)
    reg = predictRegions(u, sec, seq)  # Gets variable/interacting regions based on
    coms = getComponents(u, reg)
    plddt = getpLDDT(u, reg,raw)
    scores = buildScores(seq, sec, reg, coms, plddt)
    d = {
        "Name": name,
        "Sequence": seq,
        "Region": reg,
        "Secondary": sec,
        "Components": coms,
        "pLDDT": plddt,
        "Scoring": scores,
    }
    return d, u


#################################################
##############Plotting###########################
#################################################


def flattenPoints(d, prefix=""):
    # Recursively walks a nested dict of points (arbitrary depth) and returns a flat {label: point} dict,
    # joining nested keys with "_" to build each label
    points = {}
    for key, value in d.items():
        label = key if not prefix else f"{prefix}_{key}"
        if isinstance(value, dict):
            points.update(flattenPoints(value, label))
        else:
            points[label] = value
    return points


def getChainLines(u, t, y, x):
    # Transforms every CA atom position the same way as transformCOMs (via transformStructure) and
    # groups the resulting backbone trace by chain
    chain_dict = {"A": "MHCa", "B": "MHCb", "C": "Peptide", "D": "TCRa", "E": "TCRb"}
    utrans = transformStructure(u.copy(), t, y, x)
    lines = {}
    for chain_id, name in chain_dict.items():
        ca = utrans.select_atoms(f"segid {chain_id} and name CA")
        if len(ca) > 0:
            lines[name] = ca.positions
    return lines


def plotComponents(jsdat, u, plotDir):
    # Takes every COM point (not just the top-level ['COM'] entries) from jsdat['Components']['COM'],
    # and every chain's CA backbone, and transforms them the same way scoreComponents does (peptide at
    # origin, TCR pointing +y, MHCa->MHCb direction pointing +x). Saves a standalone 2D (x,y) line plot of
    # the backbone (colored per chain) as structure_flat.jpg, and a COM scatter plot to plotDir with the
    # same backbone lines overlaid at 30% opacity.
    comps = jsdat["Components"]["COM"]
    seg = jsdat["Components"]["Segments"]

    pep = comps["Peptide"]["COM"]
    tcr = comps["TCR"]["COM"]
    mhcha = comps["MHCa"]["MHCha"]
    try:
        mhchb = comps["MHCa"]["MHChb"]
    except:
        mhchb = comps["MHCb"]["MHChb"]
    xdir = mhchb - mhcha
    transCOM = transformCOMs(comps, -pep, y = tcr - pep, x = xdir)
    points = flattenPoints(transCOM)
    lines = getChainLines(u, -pep, y = tcr - pep, x = xdir)

    # Standalone line plot of the transformed CA backbone, colored per chain
    fig, ax = plt.subplots()
    for name, pos in lines.items():
        ax.plot(pos[:, 0], pos[:, 1], label=name, linewidth=4)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.legend()
    fig.savefig(os.path.join(plotDir, "structure_flat.jpg"))
    plt.close(fig)

    # COM scatter plot with the same backbone lines overlaid at 30% opacity
    fig, ax = plt.subplots()
    for name, pos in lines.items():
        ax.plot(pos[:, 0], pos[:, 1], alpha=0.3, label=name, linewidth=4)
    for label, point in points.items():
        if np.any(np.isnan(point)):
            continue
        ax.scatter(point[0], point[1])
        ax.text(point[0], point[1], label)

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.legend()

    fig.savefig(os.path.join(plotDir, f"COM_flat.png"))
    plt.close(fig)


######################################
#########   Main #####################
######################################

# Converts a directory of raw TCR-pMHC .pdb structures into the per-structure ./dat/<name>/
# subdirectories (Structure/Complex.pdb + <name>.json) the rest of the pipeline
# (ExtractFeaturestoJSON.py, the DL/ML dataset builders, PrepNetTCR.py, ...) reads from.
# seqfile is a CSV with one row per structure (column "ID", matching the pdb filename
# stem) and one reference sequence column per chain role, in order: MHCa, MHCb, Peptide,
# TCRa, TCRb (blank/omit a role that isn't present in that structure) - used to identify
# each raw chain by sequence similarity rather than relying on its original chain id.
# ExtractSequences.py can generate a seqfile of this shape from a directory of structures
# that already have their chains labeled A-E.
if __name__ == "__main__":
    argparser = argparse.ArgumentParser(
        description="Set up ./dat/<name>/ structure+JSON folders from raw TCR-pMHC PDB files"
    )
    argparser.add_argument("input", help="Directory containing raw .pdb structures")
    argparser.add_argument("-s", "--seqfile", required=True,
                           help="CSV with reference chain sequences per structure "
                                "(columns: ID, MHCa, MHCb, Peptide, TCRa, TCRb)")
    argparser.add_argument("-o", "--output", default="./dat",
                           help="Output directory for the per-structure folders (default: ./dat)")
    argparser.add_argument("-p", "--prefix", default=None,
                           help="Optional prefix prepended to every structure's folder/name "
                                "(e.g. so a second batch doesn't collide with an existing one)")
    args = argparser.parse_args()

    if not os.path.isdir(args.input):
        raise FileNotFoundError(f"Input directory {args.input} does not exist.")
    if not os.path.isfile(args.seqfile):
        raise FileNotFoundError(f"Sequence file {args.seqfile} does not exist.")

    os.makedirs(args.output, exist_ok=True)
    processed = setupData(args.input, args.output, args.seqfile, args.prefix)

    # Manifest of just this call's structures, so a caller (e.g. Run.sh's feature
    # extraction/inference steps) can act on only what's new even when --output is a
    # shared, accumulating directory like the default ./dat.
    manifest = os.path.join(args.output, "new_structures.txt")
    with open(manifest, "w") as f:
        f.write("\n".join(processed) + ("\n" if processed else ""))
    print(f"Processed {len(processed)} of {len(os.listdir(args.input))} input files "
         f"-> {args.output} (new_structures.txt lists the {len(processed)} successes)")

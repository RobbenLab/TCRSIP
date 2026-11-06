#!/usr/bin/env python3

###################################
########Imports################
#################################
import os
import argparse
from pathlib import Path
from struct import Struct
import sys
import numpy as np
import pandas as pd
import json
import shutil
import time
import scipy
import json

from math import nan, sqrt
from math import pow

from Bio.PDB.kdtrees import KDTree
from Bio.PDB import PDBIO

import collections
import math
from collections.abc import MutableMapping

import Bio
import Bio.PDB
import Bio.SeqRecord
from Bio.PDB import (
    PDBParser,
    MMCIFIO,
    Superimposer,
    SASA,
    Structure,
    Model,
    Entity,
    NeighborSearch,
)
from Bio.PDB.PDBIO import PDBIO
from Bio import SeqIO

# from pymol import cmd

from skspatial.objects import Line, Points, Plane
from skspatial.plotting import plot_3d

import MDAnalysis as mda
from MDAnalysis.analysis import align

# from MDAnalysis.tests.datafiles import CRD, PSF, DCD, DCD2
# import nglview as nv
from tqdm import tqdm
import flatdict

from pdbfixer import PDBFixer
from openmm.app import PDBFile

from pyrosetta import *
# from pyrosetta.toolbox.atom_pair_energy import _atom_pair_energy_table as pair_energy
from pyrosetta.toolbox.atom_pair_energy import print_atom_pair_energy_table as pair_energy
from pyrosetta.toolbox.atom_pair_energy import _reisude_pair_energies as all_energies
init(extra_options="-mute all", silent=True)
# init()
scorefxn = get_fa_scorefxn()


###################################
######## Description################
###################################

##TCR Analysis script
##This python script can be run on a standard directory to generate all of the analysis steps for TCR folding

##Usage:
# python ExtractFeaturestoJSON.py -i <input_directory> -o <output_csv_file>
#python ./src/ExtractFeaturestoJSON.py -i ./dat -o ./data/Features/features.csv 

##The script requires a file structure of:

# /dat
#   /pdb1
#       /Plots
#       /Structure
#           /pdb1.pdb
#       /Contacts
#       /pdb1.json
#   /pdb2
#   .....

## The JSON is organized into 6 different keys
# Sequence: Contains the raw sequence as input to the folding, the output sequence from the folding, and the predicted IMGT numbering
# Region: Contains the residue positions for the predicted regions of the TCR
# Secondary: Contains the secondary structure predictions for the TCR
# Components: Contains the positions of the Center of Masses, Segment lines, and the bounding box/spheres of the TCR
# pLDDT: Contains the predicted lDDT scores for the TCR (If from AlphaFold, otherwise meaningless)
# Scoring: Cumulative scores from the other sections

# The following script will create a new key in the JSON called "features" which will contain the following subkeys:

# Sequence
# MHC_class
# Num_residues: Number of residues per region and chains of the TCR
# Num_hydrophobic: Number of hydrophobic residues per region and chains of the TCR
# Num_negative: Number of negatively charged residues per region and chains of the TCR
# Num_positive: Number of positively charged residues per region and chains of the TCR
# Num_polar: Number of polar residues per region and chains of the TCR
# Num_aromatic: Number of aromatic residues per region and chains of the TCR
# Num_basic: Number of basic residues per region and chains of the TCR
# Num_acidic: Number of acidic residues per region and chains of the TCR
# Num_aliphatic: Number of aliphatic residues per region and chains of the TCR
# Num_special: Number of special residues per region and chains of the TCR
# Num_hydrophilic: Number of hydrophilic residues per region and chains of the TCR
# Mean_Hydropathy_Index: The mean of the hydropathy index of the region or chain
# Total_mass: Total molecular mass of the region or chain

# Structure (Will also save a contact map under /pdb1/Contact/)
# r: Distance from pMHC COM to TCR COM
# tilt: Angle of TCR over pMHC
# twist: angle of TCR on top of pMHC
# Buried/Solvent accessible Surface area: of each region and chain
# Min distance: between regions
# Max distance: between regions
# Min dist residue charge/polarity/hydrophobicity etc.... (can only use bit encoded)
# Min max distance to COM: minimum/maximum distance between region residue and COM
# region_angle: The angle formed by the segments made by two regions
# within 5-10-15-20

# Physics
# Num_clashes: pairwise between regions and chains
# Num_hbonds: pairwise between regions and chains (main chain-side chain, main-main, side-side?)
# Num_disulfide_bonds: pairwise between regions and chains
# Num_hydrophobic
# Num_ionic: positive and negative charged interaction (salt bridge)
# Num aro-aro: Number of interactions between aromatic-aromatic groups
# Num aro-sul: Number of interactions between aromatic and sulfur groups
# Num catpi: Number of cation-pi interactions
# Rosetta_Energy


###############################################
##########  Constants   #######################
###############################################

residue_prop = pd.DataFrame(
    {
        "Letter": [
            "A",
            "R",
            "N",
            "D",
            "C",
            "Q",
            "E",
            "G",
            "H",
            "I",
            "L",
            "K",
            "M",
            "F",
            "P",
            "S",
            "T",
            "W",
            "Y",
            "V",
        ],
        "Abbr": [
            "ALA",
            "ARG",
            "ASN",
            "ASP",
            "CYS",
            "GLN",
            "GLU",
            "GLY",
            "HIS",
            "ILE",
            "LEU",
            "LYS",
            "MET",
            "PHE",
            "PRO",
            "SER",
            "THR",
            "TRP",
            "TYR",
            "VAL",
        ],
        "Class": [
            "Aliphatic",
            "Positive",
            "Amide",
            "Negative",
            "Thiol",
            "Amide",
            "Negative",
            "Aliphatic",
            "Positive",
            "Aliphatic",
            "Aliphatic",
            "Positive",
            "Thioether",
            "Aromatic",
            "Cyclic",
            "Hydroxylic",
            "Hydroxylic",
            "Aromatic",
            "Aromatic",
            "Aliphatic",
        ],  # Class
        "Charge": [
            0,
            1,
            0,
            -1,
            0,
            0,
            -1,
            0,
            1,
            0,
            0,
            1,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
        ],  # Charge
        "Polarity": [
            "Nonpolar",
            "Polar",
            "Polar",
            "Base",
            "Acid",
            "Polar",
            "Base",
            "Nonpolar",
            "Acid",
            "Nonpolar",
            "Nonpolar",
            "Acid",
            "Nonpolar",
            "Nonpolar",
            "Nonpolar",
            "Polar",
            "Polar",
            "Nonpolar",
            "Acid",
            "Nonpolar",
        ],  # chem polarity
        "Hydropathy": [
            1.8,
            -4.5,
            -3.5,
            -3.5,
            2.5,
            -3.5,
            -3.5,
            -0.4,
            -3.2,
            4.5,
            3.8,
            -3.9,
            1.9,
            2.8,
            -1.6,
            -0.8,
            -0.7,
            -0.9,
            -1.3,
            4.2,
        ],
        "Molarmass": [
            89.094,
            174.203,
            132.119,
            133.104,
            121.154,
            146.146,
            147.131,
            75.067,
            155.156,
            131.175,
            131.175,
            146.189,
            149.208,
            165.192,
            115.132,
            105.093,
            119.119,
            204.228,
            181.191,
            117.148,
        ],
    }
)

amino = {
    "hydrophobic": ["ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "TYR"],
    "disulphide": ["CYS"],
    "ionic": ["ARG", "LYS", "HIS", "ASP", "GLU"],
    "aroaro": ["PHE", "TRP", "TYR"],
    "arosul": ["PHE", "TRP", "TYR", "CYS", "MET"],
    "cationpi": ["LYS", "ARG", "PHE", "TRP", "TYR"],
    "all": [
        "ALA",
        "GLY",
        "ILE",
        "LEU",
        "PRO",
        "VAL",
        "PHE",
        "TRP",
        "TYR",
        "ASP",
        "GLU",
        "ARG",
        "HIS",
        "LYS",
        "SER",
        "THR",
        "CYS",
        "MET",
        "ASN",
        "GLN",
    ],
}

# vdW radii taken from:
# https://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page)
#
# Radii for CL, K, NA, etc are _not_ ionic radii.
#
# References:
# A. Bondi (1964). "van der Waals Volumes and Radii".
# M. Mantina, A.C. et al., J. Phys. Chem. 2009, 113, 5806.
ATOMIC_RADII: MutableMapping[str, float] = collections.defaultdict(lambda: 2.0)
ATOMIC_RADII.update(
    {
        "H": 1.200,
        "HE": 1.400,
        "C": 1.700,
        "N": 1.550,
        "O": 1.520,
        "F": 1.470,
        "NA": 2.270,
        "MG": 1.730,
        "P": 1.800,
        "S": 1.800,
        "CL": 1.750,
        "K": 2.750,
        "CA": 2.310,
        "NI": 1.630,
        "CU": 1.400,
        "ZN": 1.390,
        "SE": 1.900,
        "BR": 1.850,
        "CD": 1.580,
        "I": 1.980,
        "HG": 1.550,
    }
)
atom_radii = {
    "H": 1.20,  # Who cares about hydrogen?? I do
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "S": 1.80,
    "F": 1.47,
    "P": 1.80,
    "CL": 1.75,
    "MG": 1.73,
}

##########################
####Features        #####
###########################


def CalcSequenceFeatures(js):
    # This function takes in a JSON file and outputs a dict which contains multiple features based on sequence
    outdict = {}
    # save class
    if js["Sequence"]["MHCb"]["Seq"] == "":
        outdict["MHC_class"] = 1
    else:
        outdict["MHC_class"] = 2
    Total_dict = {
        "Total": js["Sequence"]["MHCa"]["Seq"]
        + js["Sequence"]["MHCb"]["Seq"]
        + js["Sequence"]["Peptide"]["Seq"]
        + js["Sequence"]["TCRa"]["Seq"]
        + js["Sequence"]["TCRb"]["Seq"],
        "Total_MHCa": js["Sequence"]["MHCa"]["Seq"],
        "Total_MHCb": js["Sequence"]["MHCb"]["Seq"],
        "Total_Peptide": js["Sequence"]["Peptide"]["Seq"],
        "Total_TCRa": js["Sequence"]["TCRa"]["Seq"],
        "Total_TCRb": js["Sequence"]["TCRb"]["Seq"],
    }
    region_seqs = {}
    for r in js["Region"]:
        region = js["Region"][r]
        for k, v in region.items():
            if k not in ["TCRva", "TCRca", "TCRvb", "TCRcb"]:
                region_seqs[k] = js["Sequence"][r]["Seq"][v[0] : v[1]]

    all_seqs = Total_dict | region_seqs
    # Now run on each sequence of dicts
    for i, s in all_seqs.items():
        outdict[i] = getSeqData(s)
    # When combining into a single CSV, make sure to calc the perc of each region for each measurment using Num_residues
    return outdict


def CalcStructureFeatures(js, struct, contact):
    # Now calculate all structure features and return as a dict
    all_regions = getAllRegions(js)
    structdict = {}
    chainnames = ["MHCa", "MHCb", "Peptide", "TCRa", "TCRb"]
    chainlength = [len(js["Sequence"][k]["Seq"]) for k in chainnames]
    if np.array(chainlength).sum() != contact.shape[0]:
        raise ValueError(
            "Sum of chain lengths does not match contact map shape, feature extraction failed. Sum of chain lengths: {}, contact map shape: {}".format(
                np.array(chainlength).sum(), contact.shape
            )
        )

    # r: Distance from pMHC COM to TCR COM
    # tilt: Angle of TCR over pMHC
    # twist: angle of TCR on top of pMHC
    r, theta, phi = calcCOM(js)
    structdict["r"] = r
    structdict["tilt"] = theta
    structdict["twist"] = phi

    # Buried/Solvent accessible Surface area: of each region and chain
    structdict["SASA"] = getSASA(struct, js, all_regions)
    # Min-max residue distance: between regions
    structdict["distance"] = getDistanceFeatures(js, struct, contact, all_regions)
    # Min dist residue charge/polarity/hydrophobicity etc.... (can only use bit encoded)
    # Dist between classes? per region?
    structdict["polar_distance"] = getDistanceFeatures(
        js, struct, contact, all_regions, allow=[["S", "T", "N", "Q"]]
    )
    structdict["nonpolar_distance"] = getDistanceFeatures(
        js,
        struct,
        contact,
        all_regions,
        allow=[["A", "V", "I", "L", "M", "F", "Y", "W"]],
    )
    structdict["p_np_distance"] = getDistanceFeatures(
        js,
        struct,
        contact,
        all_regions,
        allow=[["S", "T", "N", "Q"], ["A", "V", "I", "L", "M", "F", "Y", "W"]],
    )
    structdict["charged_distance"] = getDistanceFeatures(
        js, struct, contact, all_regions, allow=[["R", "H", "K"], ["D", "E"]]
    )
    structdict["aromatic_distance"] = getDistanceFeatures(
        js, struct, contact, all_regions, allow=[["F", "Y", "W"]]
    )
    structdict["aro_charge_distance"] = getDistanceFeatures(
        js,
        struct,
        contact,
        all_regions,
        allow=[["F", "Y", "W"], ["R", "H", "K", "D", "E"]],
    )
    # region_angle: The angle formed by the segments made by two regions
    structdict["region_angle"] = getAngleFeatures(js, struct, all_regions)
    # bounding volume: The volume of each chain/region bounding sphere
    # structure volume: Volume taken up by the structure (convex hull)
    # fill: ratio bounding/structure volume
    structdict["volume"] = getVolumeFeatures(js, struct, all_regions)
    # within 5-10-15-20
    structdict["proximity"] = getProximityFeatures(js, struct, contact, all_regions)
    return structdict


def CalcPhysicsFeatures(js, struct, structdir):
    # This function calculates the physical interaction features of the structure
    # Num_clashes: pairwise between regions and chains
    # Num_hbonds: pairwise between regions and chains (main chain-side chain, main-main, side-side?)
    # Num_hydrophobic
    # Num_ionic: positive and negative charged interaction (salt bridge)
    # Num aro-aro: Number of interactions between aromatic-aromatic groups
    # Num aro-sul: Number of interactions between aromatic and sulfur groups
    # Num catpi: Number of cation-pi interactions
    # Rosetta_Energy
    physdict = {}
    all_regions = getAllRegions(js)
    centroid_map = getCentroidMap(struct)
    # First save a hydrogenated version of the structure

    # Calculate the number of clashes of atoms in each region to any other atom
    physdict["clashes"] = getClashFeatures(js, struct, all_regions)
    # Calculate number of hydrogen bonds between regions
    physdict["hbonds"] = getHBondFeatures(js, struct, all_regions)
    # Calculate the number hydrophobic interactions between regions (hydrophobic residues within 4A of each other)
    physdict["hydrophobic"] = predictInteractions(
        js, struct, centroid_map, all_regions, interaction="hydrophobic"
    )
    # Calculate the number of ionic interactions between regions (positive and negative residues within 4A of each other)
    physdict["ionic"] = predictInteractions(
        js, struct, centroid_map, all_regions, interaction="ionic"
    )
    # Calculate the number of aromatic-aromatic interactions between regions (aromatic residues within 5A of each other)
    physdict["pi-pi"] = predictInteractions(
        js, struct, centroid_map, all_regions, interaction="pi-pi"
    )
    # Calculate the number of aromatic-sulfur interactions between regions (aromatic and sulfur residues within 5A of each other)
    physdict["aro-sul"] = predictInteractions(
        js, struct, centroid_map, all_regions, interaction="aro-sul"
    )
    # Calculate the number of cation-pi interactions between regions (positive and aromatic residues within 5A of each other)
    physdict["cat-pi"] = predictInteractions(
        js, struct, centroid_map, all_regions, interaction="cat-pi"
    )
    # Calculate the rosetta energy per residue-residue and get sum per region and chain
    physdict["energy"] = predictRosettaResidueEnergy(js, structdir, all_regions)
    
    return physdict

#########################################
# Physics helpers

def predictRosettaResidueEnergy(js, structdir, regions):
    """
    Predicts the Rosetta energy per residue-residue and returns a dictionary with the sum of energies for each region and chain.
    Parameters:
    js: dict
        The JSON data containing the structure information.
    sepdir: str
        The directory containing the separated structure files for Rosetta energy calculation.
    regions: dict
        A dictionary defining the regions and their corresponding chains and residue ranges.
    """
    # This function will calculate the Rosetta energy per residue-residue and get sum per region and chain
    # It will return a dict with three keys: 
    #   "total_energy": Will contain the terms used to calculate binding energy (fake FEP)  
    #   "region_energy": Will contain the region wise sum of regions (by self and with interaction)
    #   "interaction_energy": Will contain the comparables in the interaction with the sum and average energies within just the interaction
    # In order to calculate these, we first need to get the total energy, per-residue energy for the complex and individual 
    # From the complex we will return the residue-residue energy graph to calculate the interactions we are interested in
    rosetta_energy = {}
    seq_lens = {p: len(js["Sequence"][p]["Seq"]) for p in ["MHCa","MHCb","Peptide","TCRa","TCRb"]}
    s = 0
    cum_lens = {p: (s := s + seq_lens[p]) for p in ["MHCa","MHCb","Peptide","TCRa","TCRb"]}
    #Calculate the energies for the complex and each substructure
    # files = os.listdir(sepdir)
    pose_dict = {}
    poses = ["Complex", "TCR", "pMHC"]
    for p in poses:
        pdbfile = os.path.join(structdir, f"{p}.pdb")
        pose = pose_from_pdb(pdbfile)
        score = scorefxn(pose)
        residue_energies = [pose.energies().residue_total_energy(i) for i in range(1, pose.total_residue() + 1)]
        # print(len(residue_energies))
        # residue_energy_graph = np.array(pose.energies().residue_pair_energies_array()) This returns unweighted but we want the weighted per residue pair energies
        pose_dict[p] = {
            "score": score,
            "residue_energies": residue_energies,
        }
    #Calculate basic binding energies
    Gcom = pose_dict["Complex"]["score"]
    Gtcr = pose_dict["TCR"]["score"]
    Gpmhc = pose_dict["pMHC"]["score"]
    Gbind = Gcom - (Gtcr + Gpmhc)
    rosetta_energy["total_energy"] = {
        "Gcom": Gcom,
        "Gtcr": Gtcr,
        "Gpmhc": Gpmhc,
        "Gbind": Gbind,
    }

    #Calculated region wise energies. Each region will have the bound and unbound energies and the difference between them
    translate = {
        "Bound": {
            "Keys": {
                "MHCa": "Complex",
                "MHCb": "Complex",
                "Peptide": "Complex",
                "TCRa": "Complex",
                "TCRb": "Complex",
            },
            "Start": {
                "MHCa": 0,
                "MHCb": cum_lens["MHCa"],
                "Peptide": cum_lens["MHCb"],
                "TCRa": cum_lens["Peptide"],
                "TCRb": cum_lens["TCRa"],
            }
        },
        "Unbound": {
            "Keys": {
                "MHCa": "pMHC",
                "MHCb": "pMHC",
                "Peptide": "pMHC",
                "TCRa": "TCR",
                "TCRb": "TCR",
            },
            "Start": {
                "MHCa": 0,
                "MHCb": cum_lens["MHCa"],
                "Peptide": cum_lens["MHCb"],
                "TCRa": 0,
                "TCRb": cum_lens["TCRa"]-cum_lens["Peptide"],
            }
        }
    }
    
    region_energy = {}
    for region, chains in regions.items():
        try:
            region_energy[region] = {}
            chain = chains[0]
            res_range = chains[1]
            #First get the bound
            key = translate["Bound"]["Keys"][chain]
            start = translate["Bound"]["Start"][chain]
            bound = [pose_dict[key]["residue_energies"][i + start] for i in range(res_range[0], res_range[1])]
            bound_energy = sum(bound)
            #Then unbound
            key = translate["Unbound"]["Keys"][chain]
            start = translate["Unbound"]["Start"][chain]
            unbound = [pose_dict[key]["residue_energies"][i + start] for i in range(res_range[0], res_range[1])]
            unbound_energy = sum(unbound)
            region_energy[region] = {
                "bound": bound_energy,
                "unbound": unbound_energy,
                "difference": bound_energy - unbound_energy,
            }
        except Exception as e:
            print(f"Error calculating region energy for {region}: {e}")
            region_energy[region] = {
                "bound": np.nan,
                "unbound": np.nan,
                "difference": np.nan,
            }
    rosetta_energy["region_energy"] = region_energy
    #Calculate the interaction energies from the residue-residue energy graph (like with the distance matrix, we can use the residue ranges to get the submatrix of the energy graph and sum/average it)
    energies = calcWeightedEnergyGraph(os.path.join(structdir, f"Complex.pdb"))
    interaction_energy = {}
    comparables = [
        ["Peptide", "MHCha"],
        ["Peptide", "MHChb"],
        ["MHCha", "MHChb"],
        ["Peptide", "CDR3a"],
        ["Peptide", "CDR3b"],
        ["CDR3a", "CDR3b"],
    ]
    for comp in comparables:
        try:
            region1, region2 = comp
            r1 = regions[region1]
            r2 = regions[region2]
            chain1 = regions[region1][0]
            chain2 = regions[region2][0]
            res_range1 = r1[1]
            res_range2 = r2[1]
            start1 = translate["Bound"]["Start"][chain1]
            start2 = translate["Bound"]["Start"][chain2]
            submatrix = energies[res_range1[0] + start1: res_range1[1] + start1, res_range2[0] + start2: res_range2[1] + start2]
            # print(submatrix)
            interaction_energy[f"{region1}_{region2}"] = {
                "sum": np.sum(submatrix),
                "average": np.mean(submatrix),
                "median": np.median(submatrix),
                "max": np.max(submatrix),
                "min": np.min(submatrix)
            }
        except Exception as e:
            print(f"Error calculating interaction energy for {comp}: {e}")
            interaction_energy[f"{comp[0]}_{comp[1]}"] = {
                "sum": np.nan,
                "average": np.nan,
                "median": np.nan,
                "max": np.nan,
                "min": np.nan
            }
    rosetta_energy["interaction_energy"] = interaction_energy
    return rosetta_energy

def calcWeightedEnergyGraph(pdb_file):
    #This function will calculate the weighted residue-residue energy graph for a given PDB structure using Rosetta
    #It will return a numpy array representing the weighted energy graph
    pose = pose_from_pdb(pdb_file)
    score = scorefxn(pose)
    energies = pose.energies()
    energy_graph = energies.energy_graph()
    # (fa_atr 1) (fa_rep 0.55) (fa_sol 1) (fa_intra_rep 0.005)
    # (fa_intra_sol_xover4 1) (lk_ball_wtd 1) (fa_elec 1)
    # (pro_close 1.25) (hbond_sr_bb 1) (hbond_lr_bb 1)
    # (hbond_bb_sc 1) (hbond_sc 1) (dslf_fa13 1.25) (omega 0.4)
    # (fa_dun 0.7) (p_aa_pp 0.6) (yhh_planarity 0.625) (ref 1)
    # (rama_prepro 0.45)
    # res24 = pose.residue(24)
    # res20 = pose.residue(20)
    # energies = all_energies(res24,pose,scorefxn,"fa_atr",0)
    # # energies = pair_energy(scorefxn,"fa_atr",res20,res24,0)
    # print([i for i in energies])
    total = pose.total_residue()
    energies = np.zeros((total, total))
    for i in range(1, total + 1):
        for j in range(1, total + 1):
            edge = energy_graph.find_energy_edge(i, j)
            if edge:
                energies[i-1, j-1] = edge.dot(scorefxn.weights())
    return np.array(energies)

def predictInteractions(js, struct, dist, regions, interaction=None,cutoff = None):
    """
    Predicts the number of interactions between residues of specified types (e.g., hydrophobic, ionic, aromatic) within a certain distance threshold.
    Parameters:
    js: dict
        The JSON data containing the structure information.
    struct: Bio.PDB.Structure.Structure object
        The structure to analyze for interactions.
    dist: np.ndarray
        A pairwise distance matrix between residues side chain centroids.
    regions: dict
        A dictionary defining the regions and their corresponding chains and residue ranges.
    interaction: str
        The type of interaction to predict ("hydrophobic", "ionic", "aro-sul","cat-pi","pi-pi","amide-pi").
    """
    interaction_dict = {}
    chainseqs = [
            js["Sequence"][k]["Seq"] for k in ["MHCa", "MHCb", "Peptide", "TCRa", "TCRb"]
        ]
    seq = list("".join(chainseqs))
    allow_dict = {
        "hydrophobic": [["A", "V", "I", "L", "M", "F", "Y", "W"]],
        "ionic": [["R", "H", "K"], ["D", "E"]],
        "aro-sul": [["F", "Y", "W"], ["C", "M"]],
        "cat-pi": [["R", "H", "K"], ["F", "Y", "W"]],
        "pi-pi": [["F", "Y", "W"]],
        "amide-pi": [["N", "Q"], ["F", "Y", "W"]],
    }
    allow_cutoff = {
        "hydrophobic": 4.0,
        "ionic": 4.0,
        "aro-sul": 5.0,
        "cat-pi": 5.0,
        "pi-pi": 5.0,
        "amide-pi": 5.0,
    }
    comparables = [
        ["Peptide", "MHCha"],
        ["Peptide", "MHChb"],
        ["MHCha", "MHChb"],
        ["Peptide", "CDR3a"],
        ["Peptide", "CDR3b"],
        ["CDR3a", "CDR3b"],
        ["MHCha", "CDR1b"],
        ["MHCha", "CDR2b"],
        ["MHChb", "CDR1a"],
        ["MHChb", "CDR2a"],
    ]
    ##################################################
    
    if interaction is None:
        raise ValueError(
            "No interaction is requested, please include one of the interactions"
        )
    else:
        if cutoff is None:
            cutoff = allow_cutoff[interaction]
        allow = allow_dict[interaction]
    
    if len(allow) == 1:
        mask1 = np.array([True] * contact.shape[0], dtype=bool)
        mask2 = np.array([True] * contact.shape[1], dtype=bool)
        for i in range(len(seq)):
            if seq[i] not in allow[0]:
                mask1[i] = False
                mask2[i] = False
    elif len(allow) == 2:
        mask1 = np.array([True] * contact.shape[0], dtype=bool)
        mask2 = np.array([True] * contact.shape[1], dtype=bool)
        for i in range(len(seq)):
            if seq[i] not in allow[0]:
                mask1[i] = False
            if seq[i] not in allow[1]:
                mask2[i] = False
    ###################################
    # cutoff = allow_cutoff[interaction]
    # Now we can go ahead and use the
    for comp in comparables:
        region1 = comp[0]
        region2 = comp[1]
        if region1 not in regions or region2 not in regions:
            continue
        r1 = residue2contact(js, struct, regions[region1])
        r2 = residue2contact(js, struct, regions[region2])
        sub_contact = dist[r1[0] : r1[1], r2[0] : r2[1]]
        # Create the subcontact mask matching either direction
        sub_mask = (
            mask1[r1[0] : r1[1]][:, np.newaxis] & mask2[r2[0] : r2[1]][np.newaxis, :]
            | mask2[r1[0] : r1[1]][:, np.newaxis] & mask1[r2[0] : r2[1]][np.newaxis, :]
        )
        sub_contact_masked = sub_contact[sub_mask]
        inter = np.sum(sub_contact_masked < cutoff)
        interaction_dict[f"{region1}_{region2}"] = inter
    
    return interaction_dict


def getCentroidMap(struct):
    """
    Returns a np array of the pairwise distances between the centroids of each residue side chain.
    Parameters:
    struct: Bio.PDB.Structure.Structure object
        The structure to analyze for centroid distances.
    """

    # First get a list of all centroid coordinates for each residue side chain (excluding hydrogens)
    sc_coord = [
        getSCcentroid(res)
        for res in struct.get_residues()
        if res.get_resname() != "SOL"
    ]
    # Calculate pairwise distances between centroids
    n = len(sc_coord)
    centroid_map = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            centroid_map[i, j] = np.linalg.norm(sc_coord[i] - sc_coord[j])
            centroid_map[j, i] = centroid_map[i, j]
    return centroid_map


def getSCcentroid(residue):
    """
    Returns the centroid of the side chain of a residue (excluding hydrogens).
    Parameters:
    residue: Bio.PDB.Residue.Residue object
        The residue to analyze for side chain centroid.
    """
    sc_atoms = [
        atom
        for atom in residue.get_atoms()
        if atom.name not in ["N", "CA", "C", "O", "H"]
    ]
    sc_coords = np.array([atom.coord for atom in sc_atoms])
    centroid = np.mean(sc_coords, axis=0)
    return centroid


def getHBondFeatures(js, struct, regions, compare=None):
    """
    Returns a dictionary for each region wise interaction with the number of hydrogen bonds between residues of the region.
    Calculates Hbond using an angle of > 120 degrees and a distance of < 4A between donor and acceptor.
    The donor is defined as the atom bonded to the hydrogen, and the acceptor is the atom that is not bonded to the hydrogen.
    The function will calculate the number of hydrogen bonds between all region-wise residues for a defined set of residues (default all).
    If compare is None, it will use a default set of region pairs to compare. Otherwise, it will use the provided list of region pairs.
    Parameters:
    js: dict
        The JSON data containing the structure information.
    struct: Bio.PDB.Structure.Structure object
        The structure to analyze for hydrogen bonds.
    regions: dict
        A dictionary defining the regions and their corresponding chains and residue ranges.
    """
    # This function will calculate the number of hydrogen bonds between all region-wise residues for a defined set of residues (default all)
    hbondfeat = {}
    hbonds = calcHydrogenBonds(struct)
    seqlength = sum(
        [
            len(js["Sequence"][k]["Seq"])
            for k in ["MHCa", "MHCb", "Peptide", "TCRa", "TCRb"]
        ]
    )
    if hbonds.shape[0] != seqlength:
        raise ValueError(
            "Sum of chain lengths does not match hydrogen bond map shape, feature extraction failed. Sum of chain lengths: {}, hydrogen bond map shape: {}".format(
                seqlength, hbonds.shape
            )
        )
    if compare is None:
        comparables = [
            ["Peptide", "MHCha"],
            ["Peptide", "MHChb"],
            ["MHCha", "MHChb"],
            ["Peptide", "CDR3a"],
            ["Peptide", "CDR3b"],
            ["CDR3a", "CDR3b"],
            ["MHCha", "CDR1b"],
            ["MHCha", "CDR2b"],
            ["MHChb", "CDR1a"],
            ["MHChb", "CDR2a"],
        ]
    else:
        comparables = compare

    for comp in comparables:
        region1 = comp[0]
        region2 = comp[1]
        if region1 not in regions or region2 not in regions:
            continue
        r1 = residue2contact(js, struct, regions[region1])
        r2 = residue2contact(js, struct, regions[region2])
        sub_hbond = hbonds[r1[0] : r1[1], r2[0] : r2[1]]
        hbondfeat[f"{region1}_{region2}"] = sub_hbond.sum()
    return hbondfeat


def get_bonded_atom(atom):
    residue = atom.get_parent()
    close = 100
    for neighbor in residue.get_atoms():
        if "H" in neighbor.name:
            continue
        dist = atom - neighbor
        if neighbor != atom:  # and dist < 2.2:  # Check distance threshold
            if dist < close:
                bonded = neighbor
    return bonded


def is_hydrogen_bond(atom1, atom2):
    """
    Determines if two atoms form a hydrogen bond based on distance and angle criteria.
    Parameters:
    atom1: Bio.PDB.Atom.Atom object
        The first atom to analyze.
    atom2: Bio.PDB.Atom.Atom object
        The second atom to analyze.

    """
    distance_cutoff = 4  # 3.5  # Angstroms
    angle_cutoff = 120  # Degrees
    # Save the donor, acceptor and hydrogen atoms, we need to get the bonded atom to the hydrogen
    if "H" not in atom1.name and "H" in atom2.name:
        acceptor = atom1
        hydrogen = atom2
        donor = get_bonded_atom(atom2)

    elif "H" in atom1.name and "H" not in atom2.name:
        donor = get_bonded_atom(atom1)
        hydrogen = atom1
        acceptor = atom2
    else:
        return False, nan, nan, nan, nan, nan
    # if "C" in acceptor.name or "C" in donor.name: # in ['N','CA','CB']: #Maybe  carbon can act but need charges
    #      return False,nan,nan,nan,nan,nan
    if not donor:
        print(hydrogen, acceptor)
        print("no donor")
    # First test if within max distance
    distance = donor - acceptor
    if distance > distance_cutoff:
        return False
    # save the vectors making up the two bonds around the hydrogen ()
    vector_dh = hydrogen.coord - donor.coord
    vector_ha = acceptor.coord - hydrogen.coord
    # Get the angle of the vectors
    angle = np.degrees(
        np.arccos(
            np.dot(vector_dh, vector_ha)
            / (np.linalg.norm(vector_dh) * np.linalg.norm(vector_ha))
        )
    )

    # If angle is greater than the cutoff then we should save (don't know if goes over 180?)
    if angle > angle_cutoff:
        # We should return a bunch of True if angle and distance good
        return True
    else:
        return False


def calcHydrogenBonds(struct):
    """
    Calculates the number of hydrogen bonds in the structure between residues.
    Parameters:
    struct: Bio.PDB.Structure.Structure object
        The structure to analyze for hydrogen bonds.
    region1: tuple
        A tuple containing the chain ID and residue range for the first region (chain_id, start_residue, end_residue).
    region2: tuple
        A tuple containing the chain ID and residue range for the second region (chain_id, start_residue, end_residue).
    """

    # Return a np array of the number of hydrogen bonds pairwise between all residues (only calc if CA of residues are within 10A of each other)
    residues = [x for x in struct.get_residues() if x.get_resname() != "SOL"]
    hbonds = np.zeros((len(residues), len(residues)))
    for i in range(len(residues)):
        for j in range(len(residues)):
            numbonds = 0
            res1 = residues[i]
            res2 = residues[j]
            if res1["CA"] - res2["CA"] > 10:
                continue
            if res1.get_parent().get_id() == res2.get_parent().get_id():
                continue
            atoms1 = [p for p in res1.get_atoms()]
            atoms2 = [p for p in res2.get_atoms()]
            for atom1 in atoms1:
                for atom2 in atoms2:
                    if "H" not in atom1.name and "H" not in atom2.name:
                        continue
                    if "H" in atom1.name and "H" in atom2.name:
                        continue
                    is_bond = is_hydrogen_bond(
                        atom1, atom2
                    )  # Save the output of is_hydrogen_bond()
                    if is_bond:
                        numbonds += 1
            hbonds[i, j] = numbonds / 2

    return hbonds


def getClashFeatures(js, struct, regions):
    # This function will calculate the number of clashes between all region-wise residues for a defined set of residues (default all)
    clashfeat = {}
    translate = {"MHCa": "A", "MHCb": "B", "Peptide": "C", "TCRa": "D", "TCRb": "E"}
    # We will get a residue wise number of clashes so that we can extract the number from each region and chain.
    clashes = calcClashes(struct)
    for r in regions:
        chain = regions[r][0]
        ch = translate[chain]
        range = regions[r][1]
        clashfeat[r] = sum(clashes[ch][range[0] - 1 : range[1]])
    return clashfeat


def calcClashes(structure, clash_cutoff=0.63):
    """
    Calculates per chain and per residue the number of clashes based on the atomic radii.
    Parameters:
    struct: Bio.PDB.Structure.Structure object
        The structure to analyze for clashes.
    clash_cutoff: float
        The cutoff for determining clashes based on atomic radii.
    """

    # Setup initial dataframe
    # get chain id and lengths
    data = {chain.id: ([0] * len(chain)) for chain in structure.get_chains()}

    # Set what we count as a clash for each pair of atoms
    clash_cutoffs = {
        i + "_" + j: (clash_cutoff * (atom_radii[i] + atom_radii[j]))
        for i in atom_radii
        for j in atom_radii
    }
    # Extract atoms for which we have a radii
    atoms = [x for x in structure.get_atoms() if x.element in atom_radii]
    coords = np.array([a.coord for a in atoms], dtype="d")
    # Build a KDTree (speedy!!!)
    kdt = Bio.PDB.kdtrees.KDTree(coords)
    # Initialize a list to hold clashes
    clashes = []
    # Iterate through all atoms
    for atom_1 in atoms:
        # Get parent residue
        try:
            residue_1 = atom_1.get_parent()
            residue_1coord = residue_1["CA"].get_vector()
            residue_1pos = residue_1.get_full_id()[3][1]
            residue_1chain = residue_1.get_full_id()[2]
            # Find atoms that could be clashing
            kdt_search = kdt.search(
                np.array(atom_1.coord, dtype="d"), max(clash_cutoffs.values())
            )
            # Get index and distance of potential clashes
            potential_clash = [(a.index, a.radius) for a in kdt_search]
            for ix, atom_distance in potential_clash:
                atom_2 = atoms[ix]
                residue_2 = atom_2.get_parent()
                residue_2coord = residue_2["CA"].get_vector()
                residue_2pos = residue_2.get_full_id()[3][1]
                residue_2chain = residue_2.get_full_id()[2]
                res_distance = np.linalg.norm(residue_1coord - residue_2coord)
                # Exclude clashes from atoms in the same residue
                if residue_1pos == residue_2pos and residue_1chain == residue_2chain:
                    continue
                # Exclude clashes from peptide bonds
                elif (atom_2.name == "C" and atom_1.name == "N") or (
                    atom_2.name == "N" and atom_1.name == "C"
                ):
                    continue
                # Exclude clashes from disulphide bridges
                elif (
                    atom_2.name == "SG" and atom_1.name == "SG"
                ) and atom_distance > 1.88:
                    continue
                if atom_distance < clash_cutoffs[atom_2.element + "_" + atom_1.element]:
                    # Clashes should be a pd dataframe with the columns: Name, Residue 1, Residue 1 pos, Atom 1, Atom 1 XYZ, Residue 2, Residue 2 pos, Atom 2, Atom 2 XYZ, Residue distance, Atom distance
                    data[residue_1chain][residue_1pos - 1] += 1
        except:
            raise ValueError(
                "Error calculating clashes for atom: {} in residue: {} in chain: {}".format(
                    atom_1.name, residue_1pos, residue_1chain
                )
            )
            continue
    return data


########################################
# Sequence Helpers


def getSeqData(s):
    seqdict = {}
    # Num Residues
    seqdict["Num_Residues"] = calcNumResidues(s)
    # Num_hydrophobic/polar
    seqdict["Num_hydrophobic"] = calcNumResidues(s, type="Nonpolar")
    seqdict["Num_polar"] = calcNumResidues(s, type="Polar")
    # Num_acidic/base
    seqdict["Num_acidic"] = calcNumResidues(s, type="Acid")
    seqdict["Num_basic"] = calcNumResidues(s, type="Base")
    # Num_neg/pos
    seqdict["Num_negative"] = calcNumResidues(s, type="Negative")
    seqdict["Num_positive"] = calcNumResidues(s, type="Positive")
    # Num_aromatic/aliphatic
    seqdict["Num_aromatic"] = calcNumResidues(s, type="Aromatic")
    seqdict["Num_aliphatic"] = calcNumResidues(s, type="Aliphatic")
    # Mean_hydropathic
    seqdict["Sum_hydropathy"] = calcSum(s, type="Hydropathy")
    # Sum_mass
    seqdict["Sum_Mass"] = calcSum(s, type="Molar Mass")
    return seqdict


def calcNumResidues(seq, type="None"):
    master = pd.Series()
    if type == "None":
        num = len(seq)
        return num

    if type in ["Nonpolar", "Polar", "Acid", "Base"]:
        master = residue_prop["Polarity"]

    if type in ["Negative", "Positive"]:
        master = residue_prop["Charge"]
        if type == "Negative":
            type = -1
        else:
            type = 1

    if type in ["Aromatic", "Aliphatic", "Thiol", "Amide"]:
        master = residue_prop["Class"]

    master.index = residue_prop["Letter"].tolist()
    values = [master.loc[s] for s in list(seq)]
    # Find the number that match
    num = values.count(type)

    return num


def calcSum(seq, type="None"):
    master = pd.Series()
    if type == "Hydropathy":
        master = residue_prop.set_index("Letter")["Hydropathy"]

    if type == "Molar Mass":
        master = residue_prop.set_index("Letter")["Molarmass"]

    master.index = residue_prop["Letter"].tolist()
    values = np.array([master.loc[s] for s in list(seq)])
    sum = values.sum()
    return sum


######################################
# Structure Helpers
def getProximityFeatures(js, struct, contact, regions, compare=None):
    # This function will calculate the proportion of residues within 5, 10, 15, and 20 angstroms of each other for each region pair
    proximityfeat = {}
    if compare is None:
        comparables = [
            ["Peptide", "MHCha"],
            ["Peptide", "MHChb"],
            ["MHCha", "MHChb"],
            ["Peptide", "CDR3a"],
            ["Peptide", "CDR3b"],
            ["CDR3a", "CDR3b"],
            ["MHCha", "CDR1b"],
            ["MHCha", "CDR2b"],
            ["MHChb", "CDR1a"],
            ["MHChb", "CDR2a"],
        ]
    else:
        comparables = compare
    thresholds = [5, 10, 15, 20]
    for comp in comparables:
        try:
            region1, region2 = comp
            r1 = residue2contact(js, struct, regions[region1])
            r2 = residue2contact(js, struct, regions[region2])
            sub_contact = contact[r1[0] : r1[1], r2[0] : r2[1]]
            for t in thresholds:
                key = f"{region1}_{region2}_{t}"
                proximityfeat[key] = np.sum(sub_contact < t) / sub_contact.size
        except Exception as e:
            print(f"Error calculating proximity features for {comp}: {e}")
            for t in thresholds:
                key = f"{region1}_{region2}_{t}"
                proximityfeat[key] = np.nan

    return proximityfeat


def getDistanceFeatures(js, struct, contact, regions, allow=None, compare=None):
    # This function will calculate the min, max, and mean distance between all region-wise residues for a defined set of residues (default all)
    distfeat = {}
    chainseqs = [
        js["Sequence"][k]["Seq"] for k in ["MHCa", "MHCb", "Peptide", "TCRa", "TCRb"]
    ]
    seq = list("".join(chainseqs))
    # To be flexible, we want to define a set of residues to use in determining min, max, average distance
    # We will use allow to get a mask that is length of the contact map and True for residues we want to use in the distance calculations
    if allow is None:
        mask1 = np.ones(contact.shape[0], dtype=bool)
        mask2 = np.ones(contact.shape[1], dtype=bool)
    elif len(allow) == 1:
        mask1 = np.array([True] * contact.shape[0], dtype=bool)
        mask2 = np.array([True] * contact.shape[1], dtype=bool)
        for i in range(len(seq)):
            if seq[i] not in allow:
                mask1[i] = False
                mask2[i] = False
    elif len(allow) == 2:
        mask1 = np.array([True] * contact.shape[0], dtype=bool)
        mask2 = np.array([True] * contact.shape[1], dtype=bool)
        for i in range(len(seq)):
            if seq[i] not in allow[0]:
                mask1[i] = False
            if seq[i] not in allow[1]:
                mask2[i] = False
    # Since we likely don't want to do all the comparisons let's use the most important that lead to TCR recognition
    if compare is None:
        comparables = [
            ["Peptide", "MHCha"],
            ["Peptide", "MHChb"],
            ["MHCha", "MHChb"],
            ["Peptide", "CDR3a"],
            ["Peptide", "CDR3b"],
            ["CDR3a", "CDR3b"],
            ["MHCha", "CDR1b"],
            ["MHCha", "CDR2b"],
            ["MHChb", "CDR1a"],
            ["MHChb", "CDR2a"],
        ]
    else:
        comparables = compare
    # Now we can extract the distances from the contact map for each of the comparables and calculate the min, max, and mean distances
    for comp in comparables:
        region1 = comp[0]
        region2 = comp[1]
        if region1 not in regions or region2 not in regions:
            continue
        r1 = residue2contact(js, struct, regions[region1])
        r2 = residue2contact(js, struct, regions[region2])
        sub_contact = contact[r1[0] : r1[1], r2[0] : r2[1]]
        # Create the subcontact mask matching either direction
        sub_mask = (
            mask1[r1[0] : r1[1]][:, np.newaxis] & mask2[r2[0] : r2[1]][np.newaxis, :]
            | mask2[r1[0] : r1[1]][:, np.newaxis] & mask1[r2[0] : r2[1]][np.newaxis, :]
        )
        sub_contact_masked = sub_contact[sub_mask]
        if len(sub_contact_masked) == 0:
            distfeat[f"{region1}_{region2}"] = {
                "min": 1000,  # np.nan,
                "max": 1000,  # np.nan,
                "mean": 1000,  # np.nan,
            }
            continue
        min_dist = np.nanmin(sub_contact_masked)
        max_dist = np.nanmax(sub_contact_masked)
        mean_dist = np.nanmean(sub_contact_masked)
        distfeat[f"{region1}_{region2}"] = {
            "min": min_dist,
            "max": max_dist,
            "mean": mean_dist,
        }

    return distfeat


def getAngleFeatures(js, struct, regions, compare=None):
    # This function will calculate the angle between two regions based on the segment lines defined in the JSON
    anglefeat = {}
    if compare is None:
        comparables = [
            ["Peptide", "MHCha"],
            ["Peptide", "MHChb"],
            ["MHCha", "MHChb"],
            ["Peptide", "CDR3a"],
            ["Peptide", "CDR3b"],
            ["CDR3a", "CDR3b"],
            ["MHCha", "CDR1b"],
            ["MHCha", "CDR2b"],
            ["MHChb", "CDR1a"],
            ["MHChb", "CDR2a"],
        ]
    else:
        comparables = compare
    # Now we can extract the distances from the contact map for each of the comparables and calculate the min, max, and mean distances
    for comp in comparables:
        try:
            region1 = comp[0]
            region2 = comp[1]
            chain1 = regions[region1][0]
            chain2 = regions[region2][0]
            if region1 not in regions or region2 not in regions:
                continue
            seg1 = np.array(js["Components"]["Segments"][chain1][region1])
            seg2 = np.array(js["Components"]["Segments"][chain2][region2])
            # Normalize the segment vector directions
            v1 = seg1[1] - seg1[0]
            v2 = seg2[1] - seg2[0]
            angle = calc3DAngle(v1, v2)
        except Exception as e:
            print(f"Error calculating angle for {comp}: {e}")
            angle = np.nan
        anglefeat[f"{region1}_{region2}"] = {"angle": angle}
    return anglefeat


def getVolumeFeatures(js, struct, regions):
    # This function will calculate three values for each region and chain
    # bounding volume: The volume of each chain/region bounding sphere
    # structure volume: Volume taken up by the structure (convex hull)
    # fill: ratio bounding/structure volume
    volfeat = {}
    for r, v in regions.items():
        ch = v[0]
        range = v[1]
        # get the bounding volume from the JSON
        try:
            if r in js["Components"]["Bounding"][ch].keys():
                radius = js["Components"]["Bounding"][ch][r]["Sphere"][0]
            else:
                coords = getRegionCoords(struct, ch, range)
                center = coords.mean(axis=0)
                radius = np.max(np.linalg.norm(coords - center, axis=1))
            bounding_volume = (4 / 3) * np.pi * radius**3
            # Calculate the convex hull volume
            structure_volume = calcConvexHullVolume(struct, (ch, range))
            fill = structure_volume / bounding_volume
            volfeat[r] = {
                "bounding_volume": bounding_volume,
                "structure_volume": structure_volume,
                "fill": fill,
            }
        except Exception as e:
            print(f"Error calculating volume for {r}: {e}")
            volfeat[r] = {
                "bounding_volume": np.nan,
                "structure_volume": np.nan,
                "fill": np.nan,
            }
    return volfeat


def calcConvexHullVolume(struct, region):
    # This function will calculate the volume of the convex hull of the structure
    chain = region[0]
    range = region[1]
    coords = getRegionCoords(struct, chain, range)
    hull = scipy.spatial.ConvexHull(coords)
    return hull.volume


def getRegionCoords(struct, chain, range):
    # This function will return the coordinates of the region defined by the chain and range
    coords = []
    translate = {"MHCa": "A", "MHCb": "B", "Peptide": "C", "TCRa": "D", "TCRb": "E"}
    for model in struct:
        for c in model:
            if c.get_id() == translate[chain]:
                for r in c:
                    res_id = r.get_id()[1]
                    if res_id >= range[0] and res_id < range[1]:
                        try:
                            coords.append(r["CA"].get_coord())
                        except KeyError:
                            continue
    return np.array(coords)


def calc3DAngle(x, y):
    # This function will calculate the angle between two normalized vectors in 3D space
    x = np.array(x)
    x = x / np.linalg.norm(x)
    y = np.array(y)
    y = y / np.linalg.norm(y)
    cos_theta = np.dot(x, y)
    cos_theta = np.clip(cos_theta, -1.0, 1.0)
    return np.arccos(cos_theta)


def residue2contact(js, shape, region):
    # Take the range residue positions and return the contact map begin and end coordinates that correspond to the range
    chainnames = ["MHCa", "MHCb", "Peptide", "TCRa", "TCRb"]
    chainlength = [len(js["Sequence"][k]["Seq"]) for k in chainnames]
    select = region[0]
    range = region[1]
    index = chainnames.index(select)
    start = sum(chainlength[0:index])
    new_start = start + range[0]
    new_end = start + range[1]
    return (new_start, new_end)


def getContact(structure):

    residues = [r for r in structure.get_residues() if r.get_id()[0] == " "]
    distances = np.empty([len(residues), len(residues)])
    for x in range(len(residues)):
        for y in range(len(residues)):
            try:
                one = residues[x]["CA"].get_coord()
                two = residues[y]["CA"].get_coord()
                # distances[(x * len(residues)) + y] = np.linalg.norm(one-two)
                distances[x, y] = np.abs(np.linalg.norm(one - two))
            except:
                # print("No alpha carbon, distance nan")
                #  distances[(x * len(residues)) + y] = np.nan
                distances[x, y] = np.nan

    return distances


def calcCOM(js):
    # vtcr_com = np.concatenate(
    #     [
    #         np.array(js["Components"]["COM"]["TCRa"]["TCRva"]),
    #         np.array(js["Components"]["COM"]["TCRb"]["TCRvb"]),
    #     ]
    # ).mean(axis=0)
    # if js["Sequence"]["MHCb"]["Seq"] == "":
    #     mhc_com = np.array(js["Components"]["COM"]["MHCa"]["COM"])
    # else:
    #     mhc_com = np.concatenate(
    #         [
    #             np.array(js["Components"]["COM"]["MHCa"]["COM"]),
    #             np.array(js["Components"]["COM"]["MHCb"]["COM"]),
    #         ]
    #     ).mean(axis=0)
    vtcr_com = js["Components"]["COM"]["TCR"]["COM"]
    mhc_com = js["Components"]["COM"]["MHC"]["COM"]
    dx, dy, dz = np.subtract(vtcr_com, mhc_com)
    r = np.sqrt(np.sum(np.square(np.subtract(vtcr_com, mhc_com))))
    theta = np.degrees(np.arctan2(dy, dx))
    phi = np.degrees(np.arccos(dz / r))
    return r, theta, phi


def getSASA(struct, js, regions):
    # This function will calculate per chain and region the total surface area, the solvent accessible surface area, and the buried surface area, and percent of each
    # Total area = SASA per residue summed for each chain by itself (region sum)
    # SASA = SASA per residue summed for each chain in the complex (region sum)
    # BSA = Total area - SASA
    # percent buried = BSA / Total area
    chain_ids = [c.get_id() for c in struct.get_chains()]
    target = {"MHCa": "A", "MHCb": "B", "Peptide": "C", "TCRa": "D", "TCRb": "E"}
    # Now we want to subset target by which values are in chain_ids
    translate = {k: v for k, v in target.items() if v in chain_ids}

    sr = SASA.ShrakeRupley()  # ShrakeRupley()

    all = subset_chains(struct, chain_ids)
    sr.compute(all, level="R")

    all_sasa = {k: extractSASA(all, v) for k, v in translate.items()}
    sub_sasa = {}
    for k, v in translate.items():
        sub = subset_chains(struct, [v])
        sr.compute(sub, level="R")
        sub_sasa[k] = extractSASA(sub, v)

    sasadict = {}
    sa = np.array([item for sublist in sub_sasa.values() for item in sublist]).sum()
    sasa = np.array([item for sublist in all_sasa.values() for item in sublist]).sum()
    bsa = sa - sasa
    sasadict["Total"] = {
        "sa": sa,
        "sasa": sasa,
        "bsa": sa - sasa,
        "percent_buried": bsa / sa,
    }

    for r, v in regions.items():
        ch = v[0]
        range = v[1]
        sa = np.array(sub_sasa[ch][range[0] : range[1]]).sum()
        sasa = np.array(all_sasa[ch][range[0] : range[1]]).sum()
        bsa = sa - sasa
        percent_buried = bsa / sa
        sasadict[r] = {
            "sa": sa,
            "sasa": sasa,
            "bsa": bsa,
            "percent_buried": percent_buried,
        }
    return sasadict


def extractSASA(structure, ch):
    SASA = []
    for model in structure:
        for chain in model:
            if chain.get_id() == ch:
                for residue in chain:
                    SASA.append(residue.sasa)
    return SASA


def subset_chains(structure, chains):
    model = structure[0]
    struct = Structure.Structure("1")
    mod = Model.Model("1")
    struct.add(mod)
    for chain in chains:
        mod.add(model[chain])
    return struct


def getAllRegions(js):
    # This function will return a dict of all the regions and their corresponding sequences
    all_regions = {}
    for r in js["Region"]:
        all_regions[r] = (r, [0, len(js["Sequence"][r]["Seq"])])
        region = js["Region"][r]
        for k, v in region.items():
            if v == [0, 0] or v == [np.nan, np.nan]:
                continue
            all_regions[k] = (r, v)
    return all_regions


class ShrakeRupley:
    """Calculates SASAs using the Shrake-Rupley algorithm."""

    def __init__(self, probe_radius=1.40, n_points=100, radii_dict=None):
        """Initialize the class.

        :param probe_radius: radius of the probe in A. Default is 1.40, roughly
            the radius of a water molecule.
        :type probe_radius: float

        :param n_points: resolution of the surface of each atom. Default is 100.
            A higher number of points results in more precise measurements, but
            slows down the calculation.
        :type n_points: int

        :param radii_dict: user-provided dictionary of atomic radii to use in
            the calculation. Values will replace/complement those in the
            default ATOMIC_RADII dictionary.
        :type radii_dict: dict

        >>> sr = ShrakeRupley()
        >>> sr = ShrakeRupley(n_points=960)
        >>> sr = ShrakeRupley(radii_dict={"O": 3.1415})
        """
        if probe_radius <= 0.0:
            raise ValueError(
                f"Probe radius must be a positive number: {probe_radius} <= 0"
            )

        self.probe_radius = float(probe_radius)

        if n_points < 1:
            raise ValueError(
                f"Number of sphere points must be larger than 1: {n_points}"
            )
        self.n_points = n_points

        # Update radii list with user provided lists.
        self.radii_dict = ATOMIC_RADII.copy()
        if radii_dict is not None:
            self.radii_dict.update(radii_dict)

        # Pre-compute reference sphere
        self._sphere = self._compute_sphere()

    def _compute_sphere(self):
        """Return the 3D coordinates of n points on a sphere.

        Uses the golden spiral algorithm to place points 'evenly' on the sphere
        surface. We compute this once and then move the sphere to the centroid
        of each atom as we compute the ASAs.
        """
        n = self.n_points

        dl = np.pi * (3 - 5**0.5)
        dz = 2.0 / n

        longitude = 0
        z = 1 - dz / 2

        coords = np.zeros((n, 3), dtype=np.float32)
        for k in range(n):
            r = (1 - z * z) ** 0.5
            coords[k, 0] = math.cos(longitude) * r
            coords[k, 1] = math.sin(longitude) * r
            coords[k, 2] = z
            z -= dz
            longitude += dl

        return coords

    def compute(self, entity, level="A"):
        """Calculate surface accessibility surface area for an entity.

        The resulting atomic surface accessibility values are attached to the
        .sasa attribute of each entity (or atom), depending on the level. For
        example, if level="R", all residues will have a .sasa attribute. Atoms
        will always be assigned a .sasa attribute with their individual values.

        :param entity: input entity.
        :type entity: Bio.PDB.Entity, e.g. Residue, Chain, ...

        :param level: the level at which ASA values are assigned, which can be
            one of "A" (Atom), "R" (Residue), "C" (Chain), "M" (Model), or
            "S" (Structure). The ASA value of an entity is the sum of all ASA
            values of its children. Defaults to "A".
        :type entity: Bio.PDB.Entity

        >>> from Bio.PDB import PDBParser
        >>> from Bio.PDB.SASA import ShrakeRupley
        >>> p = PDBParser(QUIET=1)
        >>> # This assumes you have a local copy of 1LCD.pdb in a directory called "PDB"
        >>> struct = p.get_structure("1LCD", "PDB/1LCD.pdb")
        >>> sr = ShrakeRupley()
        >>> sr.compute(struct, level="S")
        >>> print(round(struct.sasa, 2))
        7053.43
        >>> print(round(struct[0]["A"][11]["OE1"].sasa, 2))
        9.64
        """
        is_valid = hasattr(entity, "level") and entity.level in {"R", "C", "M", "S"}
        if not is_valid:
            raise ValueError(
                f"Invalid entity type '{type(entity)}'. "
                "Must be Residue, Chain, Model, or Structure"
            )

        if level not in _ENTITY_HIERARCHY:
            raise ValueError(f"Invalid level '{level}'. Must be A, R, C, M, or S.")
        elif _ENTITY_HIERARCHY[level] > _ENTITY_HIERARCHY[entity.level]:
            raise ValueError(
                f"Level '{level}' must be equal or smaller than input entity: {entity.level}"
            )

        # Get atoms onto list for lookup
        atoms = list(entity.get_atoms())
        n_atoms = len(atoms)
        if not n_atoms:
            raise ValueError("Entity has no child atoms.")

        # Get coordinates as a numpy array
        # We trust DisorderedAtom and friends to pick representatives.
        coords = np.array([a.coord for a in atoms], dtype=np.float64)

        # Pre-compute atom neighbors using KDTree
        kdt = KDTree(coords, 10)

        # Pre-compute radius * probe table
        radii_dict = self.radii_dict
        radii = np.array([radii_dict[a.element] for a in atoms], dtype=np.float64)
        radii += self.probe_radius
        twice_maxradii = np.max(radii) * 2

        # Calculate ASAs
        asa_array = np.zeros((n_atoms, 1), dtype=np.int64)
        ptset = set(range(self.n_points))
        for i in range(n_atoms):
            r_i = radii[i]

            # Move sphere to atom
            s_on_i = (np.array(self._sphere, copy=True) * r_i) + coords[i]
            available_set = ptset.copy()

            # KDtree for sphere points
            kdt_sphere = KDTree(s_on_i, 10)

            # Iterate over neighbors of atom i
            for jj in kdt.search(coords[i], twice_maxradii):
                j = jj.index
                if i == j:
                    continue

                if jj.radius < (r_i + radii[j]):
                    # Remove overlapping points on sphere from available set
                    available_set -= {
                        pt.index for pt in kdt_sphere.search(coords[j], radii[j])
                    }

            asa_array[i] = len(available_set)  # update counts

        # Convert accessible point count to surface area in A**2
        f = radii * radii * (4 * np.pi / self.n_points)
        asa_array = asa_array * f[:, np.newaxis]

        # Set atom .sasa
        for i, atom in enumerate(atoms):
            atom.sasa = asa_array[i, 0]

        residue_sasa = []

        # Aggregate values per entity level if necessary
        if level != "A":
            entities = set(atoms)
            target = _ENTITY_HIERARCHY[level]
            for _ in range(target):
                entities = {e.parent for e in entities}

            atomdict = {a.full_id: idx for idx, a in enumerate(atoms)}
            for e in entities:
                e_atoms = [atomdict[a.full_id] for a in e.get_atoms()]
                e.sasa = asa_array[e_atoms].sum()
                residue_sasa.append(e.sasa)

        return residue_sasa

#########################################Uniform

def makeRow(name, seqdict, structdict, physdict):
    # This function will take the three dictionaries and combine them into a single row for a pandas dataframe
    row = {"name": name}
    seqrow = uniformDict(seqdict)
    structrow = uniformDict(structdict)
    physrow = uniformDict(physdict)
    row.update(seqrow)
    row.update(structrow)
    row.update(physrow)
    return row

def uniformDict(seqdict):
    return dict(flatdict.FlatDict(seqdict))


#################################################
######### Setup Variables ######################
#################################################
print("Setting up variables")

parser = argparse.ArgumentParser("simple_example")
parser.add_argument("-d","--Dir",help="Primary Directory in which the above structure will go/already exists", type=str)
parser.add_argument("-o","--Out",help="Output path for the results as a csv file", type=str)
parser.add_argument("--names-file", default=None,
                    help="Optional text file of structure names (one per line) to restrict "
                         "extraction to, instead of every subdirectory in --Dir - e.g. just the "
                         "structures a run of SetupData.py has just added")

args = parser.parse_args()

dir = args.Dir
output = args.Out
names_file = args.names_file


# reffile1 = "./data/ref/1ao7.pdb"

##############################################
########    Main function ####################
###########################################
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


# Test function
# name = "1ao7"
# name = "Fake_1"
# Dir = f"./dat/{name}/"

print("Running feature extraction on all subdirectories in: ", dir)
# Get the list of directories in Dir
subdirs = [d for d in os.listdir(dir) if os.path.isdir(os.path.join(dir, d))]
if names_file is not None:
    with open(names_file) as f:
        wanted = [line.strip() for line in f if line.strip()]
    missing = sorted(set(wanted) - set(subdirs))
    if missing:
        print(f"Warning: {len(missing)} names from {names_file} have no subdirectory in {dir}: {missing[:10]}")
    subdirs = [d for d in subdirs if d in wanted]
    print(f"Restricted to {len(subdirs)} of {len(wanted)} names from {names_file}")
if not subdirs:
    raise FileNotFoundError(f"No subdirectories found in: {dir}")
rows = []
#Now we will loop through the subdirectories and extract the features for each one with tqdm progress bar
for subdir in tqdm(subdirs):
    try:
        name = subdir
        Dir = os.path.join(dir, subdir)

        #Check if structure directory exists
        structdir = os.path.join(Dir, "Structure")
        if not os.path.exists(structdir):
            raise FileNotFoundError(f"Structure directory not found: {structdir}")
        #Make contact directory if it doesn't exist
        contact_dir = os.path.join(Dir, "Contact")
        if not os.path.exists(contact_dir):
            os.makedirs(contact_dir)

        # Now we will go ahead and test the feature extraction functions on the directory (The full version will be called from an sh script for generalizeability)
        # Read in the json file and the pdb file
        js = json.load(open(os.path.join(Dir, f"{name}.json")))
        structdir = os.path.join(Dir, "Structure")
        # struct = os.path.join(structdir, f"{name}.pdb")
        struct = os.path.join(structdir, f"Complex.pdb")

        parser = Bio.PDB.PDBParser(QUIET=False)
        structure = parser.get_structure(name, struct)
        # Sequence features
        seqdict = CalcSequenceFeatures(js)
        # Contacts
        contact = getContact(structure)
        #save the contact map to a csv file
        pd.DataFrame(contact).to_csv(os.path.join(contact_dir, f"contact.csv"), index=False)

        # Structure features
        structdict = CalcStructureFeatures(js, structure, contact)
        # Physics features
        physdict = CalcPhysicsFeatures(js, structure, structdir)

        #Save the features back to the json file
        js["Features"] = {"Sequence": seqdict, "Structure": structdict, "Physics": physdict}
        with open(os.path.join(Dir, f"{name}.json"), "w") as f:
            json.dump(js, f, indent=4, cls=NumpyEncoder)
        #Now we want to get the features into a single uniform dataframe row to append into a pandas dataframe
        row = makeRow(name,seqdict, structdict, physdict)
        rows.append(row)
    except Exception as e:
        print(f"Error processing subdirectory {subdir}: {e}")

df = pd.DataFrame(rows)
#Save the dataframe to a csv file with output
df.to_csv(os.path.join(output), index=False)
print(f"Extracted features for {len(df)} structure(s) -> {output}")

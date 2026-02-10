#!/usr/bin/env python3


#python ./scripts/ExtractFeatures.py --Dir ./ModelData/Immpred/Colabfold/ -ta 20 -ma 25
# python ./scripts/ExtractFeatures.py --Dir ./ModelData/KnownStructures/Colabfold/
# python ./scripts/ExtractFeatures.py --Dir ./ModelData/PredictedInteractions/Colabfold/
# python ./scripts/ExtractFeatures.py --Dir ./ModelData/Fake/Colabfold/

import os
import argparse 
from pathlib import Path
from struct import Struct
import subprocess
import sys
import urllib.request
import numpy as np
import pandas as pd
import json
import shutil
import time 
import scipy

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
from Bio.PDB import PDBParser, MMCIFIO, Superimposer, SASA, Structure, Model, Entity, NeighborSearch
from Bio.PDB.PDBIO import PDBIO
from Bio import SeqIO
from pymol import cmd

from skspatial.objects import Line, Points, Plane
from skspatial.plotting import plot_3d

from urllib.request import urlopen
import re
from sympy import false, true
from zmq import NULL

import MDAnalysis as mda
from MDAnalysis.analysis import align
# from MDAnalysis.tests.datafiles import CRD, PSF, DCD, DCD2
#import nglview as nv
from tqdm import tqdm

##TCR Analysis script
##This python script can be run on a standard directory to generate all of the analysis steps for TCR folding
##The script requires a file structure of:

##  /Project
##      /Structures
#           /Raw
#           /Aligned
#       /Contacts                       #Contains the 2d contact maps for each pdb
#       /Angles
#       /Physics
#           /Hbond                      #Contains per file list of hydrogen bonds
#       /Scores
#           /lDDT                       #Contains the linear lDDT for each pdb
#           /Deviation                  #Contains the linear deviations from ground truth
#           /bScores                    #Contains the b property extracted from each pdb, if alphafold, this is the plDDT
#           /DockQ
#               /json                   #Contains the json outputs of DockQ

## Then you must run this script from the Project directory




###############################################
##########  Functions   #######################
###############################################

#Run the functions that will be used in this script

def subset_chains(structure,chains):
    model = structure[0]
    struct = Structure.Structure('1')
    mod = Model.Model('1')
    struct.add(mod)
    for chain in chains:
        mod.add(model[chain])
    return struct

def download_read_pdb(pdbcode, datadir, keepfile=True):
    """
    Downloads a PDB file from the Internet and saves it in a data directory.
    Then it reads and returns the structure inside.
    :param pdbcode: The standard PDB ID e.g. '3ICB'
    :param datadir: The directory where the downloaded file will be saved
    :param keepfile: if False, then the downloaded file will be deleted (default: keep the downloaded file)
    :return: a Bio.PDB Structure object or None if something went wrong
    """
    pdbfilenm = download_pdb(pdbcode, datadir)
    if pdbfilenm is None:
        return None
    struct = read_pdb(pdbcode, pdbfilenm)
    if not keepfile:
        os.remove(pdbfilenm)
    return struct
        
def download_pdb(pdbcode, datadir, downloadurl="http://files.rcsb.org/download/"):
    """
    Downloads a PDB file from the Internet and saves it in a data directory.
    :param pdbcode: The standard PDB ID e.g. '3ICB' or '3icb'
    :param datadir: The directory where the downloaded file will be saved
    :param downloadurl: The base PDB download URL, cf.
        `https://www.rcsb.org/pages/download/http#structures` for details
        Note that the unencrypted HTTP protocol is used by default
        to avoid spurious OpenSSL errors...
    :return: the full path to the downloaded PDB file or None if something went wrong
    """
    pdbfn = pdbcode + ".pdb"
    url = downloadurl + pdbfn
    outfnm = os.path.join(datadir, pdbfn)
    try:
        urllib.request.urlretrieve(url, outfnm)
        return outfnm
    except Exception as err:
        # all sorts of things could have gone wrong...
        print(str(err), file=sys.stderr)
        return None

def read_pdb(pdbcode, pdbfilenm):
    """
    Read a PDB structure from a file.
    :param pdbcode: A PDB ID string
    :param pdbfilenm: The PDB file
    :return: a Bio.PDB.Structure object or None if something went wrong
    """
    try:
        pdbparser = Bio.PDB.PDBParser(QUIET=True)   # suppress PDBConstructionWarning
        struct = pdbparser.get_structure(pdbcode, pdbfilenm)
        return struct
    except Exception as err:
        print(str(err), file=sys.stderr)
        return None 

def extract_seqrecords(pdbcode, struct):
    """
    Extracts the sequence records from a Bio.PDB structure.
    :param pdbcode: the PDB ID of the structure, needed to add a sequence ID to the result
    :param struct: a Bio.PDB.Structure object
    :return: a list of Bio.SeqRecord objects
    """
    ppb = Bio.PDB.PPBuilder()
    seqrecords = []
    for i, chain in enumerate(struct.get_chains()):
        # extract and store sequences as list of SeqRecord objects
        pps = ppb.build_peptides(chain)    # polypeptides
        seq = pps[0].get_sequence() # just take the first, hope there's no chain break
        seqid = pdbcode + chain.id
        seqrec = Bio.SeqRecord.SeqRecord(seq, id=seqid, 
            description="Sequence #{}, {}".format(i+1, seqid))
        seqrecords.append(seqrec)
    return seqrecords

def get_calphas(struct):
    """
    Extracts the C-alpha atoms from a PDB structure.
    :param struct: A Bio.PDB.Structure object.
    :return: A list of Bio.PDB.Atom objects representing the C-alpha atoms in `struct`.
    """
    calphas = [ atom for atom in struct.get_atoms() if atom.get_fullname() == " CA " ]
    return calphas

def Norm(vec):
    mag = np.linalg.norm(vec) 
    norm = np.divide(vec,mag)
    return norm

def angle_between(v1, v2):
    # """ Returns the angle in radians between vectors 'v1' and 'v2'::

    #         >>> angle_between((1, 0, 0), (0, 1, 0))
    #         1.5707963267948966
    #         >>> angle_between((1, 0, 0), (1, 0, 0))
    #         0.0
    #         >>> angle_between((1, 0, 0), (-1, 0, 0))
    #         3.141592653589793
    # """
    return np.arccos(np.clip(np.dot(v1, v2), -1.0, 1.0)) #Just return in radians for now
    # return np.rad2deg(np.arccos(np.clip(np.dot(v1, v2), -1.0, 1.0)))

def pull_bonds(chain,dist = 5):
    disulfide_bonds = []
    for residue in chain:
            if residue.get_resname() == "CYS":
                for atom in residue:
                    if atom.get_name() == "SG":
                        for other_residue in chain:
                            if other_residue.get_resname() == "CYS" and other_residue != residue:
                                for other_atom in other_residue:
                                    if other_atom.get_name() == "SG":
                                        distance = atom - other_atom
                                        if distance < dist: # We will set to 5 becuase there seems to be slight discrepencies with actual disulfide bonds 2.5:  # Typical distance for a disulfide bond
                                            disulfide_bonds.append((residue, other_residue))
    return disulfide_bonds

def find_disulfide_bonds(pdb_file):
    parser = Bio.PDB.PDBParser()
    structure = parser.get_structure("structure", pdb_file)

    chains = structure.get_chains()
    chain_ids = [c.get_id() for c in chains]
    num_chains = len(chain_ids)#(len([c for c in chains]))
    TCRA_id = chain_ids[num_chains-2]
    TCRB_id = chain_ids[num_chains-1]
    TCRA = structure[0][TCRA_id]
    TCRB = structure[0][TCRB_id]
    TCRA_bonds = pull_bonds(TCRA)
    TCRB_bonds = pull_bonds(TCRB)

    return TCRA_bonds,TCRB_bonds

def find_disulfide_bonds_nolim(pdb_file):
    parser = Bio.PDB.PDBParser()
    structure = parser.get_structure("structure", pdb_file)

    chains = structure.get_chains()
    chain_ids = [c.get_id() for c in chains]
    num_chains = len(chain_ids)#(len([c for c in chains]))
    TCRA_id = chain_ids[num_chains-2]
    TCRB_id = chain_ids[num_chains-1]
    TCRA = structure[0][TCRA_id]
    TCRB = structure[0][TCRB_id]
    TCRA_bonds = pull_bonds(TCRA,100)
    TCRB_bonds = pull_bonds(TCRB,100)

    return TCRA_bonds,TCRB_bonds

def getAngles(name,file):
    #####This function calculates two angles associated with the TCR-pMHC interaction #########################
    #####Binding Angle: Lateral angle for the CDR1 and 2 interaction with the bun of the MHC hotdog############
    #####Incident Angle: Tilt of the TCR over the plane of the pMHC complex####################################
    #####BA formulae: dot product of the MHC cross angle and the angle formed between the centroids of the ####
    #####TCR heavy chain disulfide bridges                                                           ##########
    #####IA formulae: dot product of the MHC plane normal and the vertical angle of the TCR above it###########
    #Read in the PDB file
    parser = Bio.PDB.PDBParser(QUIET = False)
    structure = parser.get_structure(name, file)
    chains = structure.get_chains()
    classone = False
    num_chains = (len([c for c in chains]))
    # if num_chains < 4 or num_chains > 5:
    #     return
    chains = structure.get_chains()
    has_pep = False
    has_b = True
    for chain in chains:
        res_no = 0
        for r in chain.get_residues():
            res_no += 1
        if res_no < 25:
            has_pep = True


    chains = structure.get_chains()
    ids = []
    for chain in chains:
        ids.append(chain.get_id())

    if "B" in ids:
        has_b = True
        A = "A"
        B = "B"
        C = "C"
        D = "D"
        E = "E"
    else:
        A = "A"
        B = "C"
        C = "D"
        D = "E"

    chains = structure.get_chains()
    for chain in chains:
        tmp = [r for r in chain.get_residues() if r.get_id()[0] == " "]
        if num_chains == 4:
            if has_pep:
                classone = True
                if chain.get_id() == A: MHCA = tmp
                if chain.get_id() == C: TCRA = tmp
                if chain.get_id() == D: TCRB = tmp
            else:
                if chain.get_id() == A: MHCA = tmp
                if chain.get_id() == B: MHCB = tmp
                if chain.get_id() == C: TCRA = tmp
                if chain.get_id() == D: TCRB = tmp
        elif num_chains == 5:
            if chain.get_id() == A: MHCA = tmp
            if chain.get_id() == B: MHCB = tmp
            if chain.get_id() == D: TCRA = tmp
            if chain.get_id() == E: TCRB = tmp
        elif num_chains == 3:
            classone = True
            if chain.get_id() == A: MHCA = tmp
            if chain.get_id() == C: TCRA = tmp
            if chain.get_id() == D: TCRB = tmp
        else:
            print("Sequence has too many chains or too little chains")
            return [nan,nan]

    #Get useful vectors of 3d coord to construct angles
    if classone == False:
        vec = MHCA[45:77]+MHCB[53:63]+MHCB[66:90] #could be innacurate of start is greater than residue 1 but at least in our set it starts mhc at 1
        MHCvec = [r["CA"].get_coord() for r in vec]
    else:
        vec = MHCA[49:85]+MHCA[139:175]
        MHCvec = [r["CA"].get_coord() for r in vec]
    #Get the Cysteines from the alpha beta chain disulfide bridges 
    # There should only be 2 disulfide bridges in the heavy chain of the alpha and beta TCR and each should only have 2 cysteines
    # #May need to improve the method for getting disulfide bridge cysteines
    # TCRAcys = [r for r in TCRA if r.get_resname() == "CYS" and r.get_full_id()[3][1] > 18] 
    # TCRBcys = [r for r in TCRB if r.get_resname() == "CYS" and r.get_full_id()[3][1] > 18]#[r for r in TCRB if r.get_full_id()[3][1] in [23,106]] #Will need to limit light/heavy chain if not TCR3d by running on first 100 residues only
    # nullarray = [NULL,NULL]
    # if len(TCRAcys) != 2:
    #     return nullarray
    # if len(TCRBcys) != 2:
    #     return nullarray
    ###This is a lot better code to get those disulfide bonds because we have a function that grabs all of them and then we just take the first one
    Abonds,Bbonds = find_disulfide_bonds(file)
    ### check if there are bonds in there
    if len(Abonds) < 1:
        print(f"Could not find Disulfide bonds for {name}")
        return [nan,nan]

    if len(Bbonds) < 1:
        print(f"Could not find Disulfide bonds for {name}")
        return [nan,nan]

    TCRAsulfa = [a.get_coord() for a in Abonds[0][0].get_atoms() if a.get_name() == "SG"]
    TCRAsulfb = [a.get_coord() for a in Abonds[0][1].get_atoms() if a.get_name() == "SG"]
    TCRAcent = np.divide(np.add(TCRAsulfa,TCRAsulfb),2)
    TCRBsulfa = [a.get_coord() for a in Bbonds[0][0].get_atoms() if a.get_name() == "SG"]
    TCRBsulfb = [a.get_coord() for a in Bbonds[0][1].get_atoms() if a.get_name() == "SG"]
    TCRBcent = np.divide(np.add(TCRBsulfa,TCRBsulfb),2)

    #Extract the coordinates for the sulfur on the cysteines for the alpha and beta chains and calculated the centroids
    # TCRAsulfa = [a.get_coord() for a in TCRAcys[0].get_atoms() if a.get_name() == "SG"]
    # TCRAsulfb = [a.get_coord() for a in TCRAcys[1].get_atoms() if a.get_name() == "SG"]
    # TCRAcent = np.divide(np.add(TCRAsulfa,TCRAsulfb),2)
    # TCRBsulfa = [a.get_coord() for a in TCRBcys[0].get_atoms() if a.get_name() == "SG"]
    # TCRBsulfb = [a.get_coord() for a in TCRBcys[1].get_atoms() if a.get_name() == "SG"]
    # TCRBcent = np.divide(np.add(TCRBsulfa,TCRBsulfb),2)

    # print(name)

    #Calculated a line parallel to the alpha carbons of the MHC helices (buns)
    MHCline = Line.best_fit(np.array(MHCvec))
    MHCplane = Plane.best_fit(np.array(MHCvec)) #Calculate the plane of the helices
    MHCcross = np.cross(MHCline.direction,MHCplane.normal) #Calculate the horizontal cross line of the helices
    ###### Calculate the line vertically through the alpha beta chains
    TCRallcoord = np.array([r["CA"].get_coord() for r in TCRA] + [r["CA"].get_coord() for r in TCRB] )
    TCRline = Line.best_fit(TCRallcoord)
    TCRplane = Plane.best_fit(TCRallcoord)
    TCRcross = np.cross(TCRline.direction,TCRplane.normal) #Is this the best way to calculate the tilt vector of the TCR? It seems easy for the plane to be calculated along the z axis
    TCRlinenorm = Norm(np.subtract(TCRBcent[0],TCRAcent[0]))

    MHCcross2d = Norm(np.array([MHCcross[0],MHCcross[1]]))
    TCRline2d = Norm(np.array([TCRlinenorm[0],TCRlinenorm[1]]))
    MHCnormal2d = Norm(np.array([MHCplane.normal[0],MHCplane.normal[2]]))
    TCRnormal2d = Norm(np.array([TCRcross[0],TCRcross[2]]))

    BA = angle_between(MHCcross2d,TCRline2d) % np.pi #Calculates the dot product of the x and y 
    IA = angle_between(MHCnormal2d,TCRnormal2d) % np.pi
    array = [BA,IA]
    return array

def getContact(name,file):

    parser = Bio.PDB.PDBParser(QUIET = False)
    structure = parser.get_structure(name, file)

    residues = [r for r in structure.get_residues() if r.get_id()[0] == " "]
    distances = np.empty([len(residues),len(residues)])
    for x in range(len(residues)):
        for y in range(len(residues)):
          try: 
            one  = residues[x]["CA"].get_coord()
            two = residues[y]["CA"].get_coord()
            # distances[(x * len(residues)) + y] = np.linalg.norm(one-two)
            distances[x,y] = np.linalg.norm(one-two)
          except:
            # print("No alpha carbon, distance nan")
            #  distances[(x * len(residues)) + y] = np.nan
            distances[x,y] = np.nan

    return distances

#This only works with alphafold output (maybe esm fold)
def getlDDT(name,file):
  #lddt_filename = f"{jobname}/{jobname}{jobname_prefix}_scores_{tag}.json"
  f = open(file)#lddt_filename)
  data = json.load(f)
  return data["plddt"],data["ptm"],data["iptm"]

#Computes the score between the ground truth and aligned 
def computeRMSD(name,query,ref):
    ground = mda.Universe(ref)
    predicted = mda.Universe(query)
    merged = mda.Merge(ground.atoms, predicted.atoms)
    print(predicted)
    rmsds = align.alignto(predicted,  # mobile
                        ground,  # reference
                        select='name CA', # selection to operate on
                        match_atoms=True) # whether to match atoms
    align.alignto(predicted, ground, select='name CA')
    return rmsds,predicted

valid_amino_acids = {
    'LLP': 'K', 'TPO': 'T', 'CSS': 'C', 'OCS': 'C', 'CSO': 'C', 'PCA': 'E', 'KCX': 'K', \
    'CME': 'C', 'MLY': 'K', 'SEP': 'S', 'CSX': 'C', 'CSD': 'C', 'MSE': 'M', \
    'ALA': 'A', 'ASN': 'N', 'CYS': 'C', 'GLN': 'Q', 'HIS': 'H', 'LEU': 'L', \
    'MET': 'M', 'MHO': 'M', 'PRO': 'P', 'THR': 'T', 'TYR': 'Y', 'ARG': 'R', 'ASP': 'D', \
    'GLU': 'E', 'GLY': 'G', 'ILE': 'I', 'LYS': 'K', 'PHE': 'F', 'SER': 'S', \
    'TRP': 'W', 'VAL': 'V', 'SEC': 'U'
    }

def check_pdb_valid_row(valid_amino_acids, l):
    if (get_pdb_rname(l) in valid_amino_acids.keys()) and (l.startswith('ATOM') or l.startswith('HETA')):
        return True
    return False

def get_pdb_atom_name(l):
    return l[12: 16].strip()

def get_pdb_rnum(l):
    return int(l[22: 27].strip())

def get_pdb_rname(l):
    return l[17: 20].strip()

def get_pdb_xyz_cb(lines):
    xyz = {}
    for l in lines:
        if get_pdb_atom_name(l) == 'CB':
            xyz[get_pdb_rnum(l)] = (float(l[30:38].strip()), float(l[38:46].strip()), float(l[46:54].strip()))
    for l in lines:
        if (get_pdb_rnum(l) not in xyz) and get_pdb_atom_name(l) == 'CA':
            xyz[get_pdb_rnum(l)] = (float(l[30:38].strip()), float(l[38:46].strip()), float(l[46:54].strip()))
    return xyz

def get_pdb_xyz_ca(lines):
    xyz = {}
    for l in lines:
        if get_pdb_atom_name(l) == 'CA':
            xyz[get_pdb_rnum(l)] = (float(l[30:38].strip()), float(l[38:46].strip()), float(l[46:54].strip()))
    return xyz

def pdb2dmap(pdbfile):
    f = open(pdbfile, mode = 'r')
    flines = f.read()
    f.close()
    lines = flines.splitlines()
    templines = flines.splitlines()
    for l in templines:
        if not l.startswith('ATOM'): lines.remove(l)
    # We have filtered out all non ATOMs at this point
    rnum_rnames = {}
    for l in lines:
        atom = get_pdb_atom_name(l)
        if atom != 'CA': continue
        if not get_pdb_rname(l) in valid_amino_acids.keys():
            print ('' + get_pdb_rname(l) + ' is unknown amino acid in ' + l)
            return
        rnum_rnames[int(get_pdb_rnum(l))] = valid_amino_acids[get_pdb_rname(l)]
    seq = ""
    for i in range(max(rnum_rnames.keys())):
        if i+1 not in rnum_rnames:
            #print (rnum_rnames)
            #print ('Warning! residue not defined for rnum = ' + str(i+1))
            seq += '-'
        else:
            seq += rnum_rnames[i+1]
    L = len(seq)
    xyz_cb = get_pdb_xyz_cb(lines)
    total_valid_residues = len(xyz_cb)
    if len(xyz_cb) != L:
        print(rnum_rnames)
        for i in range(L):
            if i+1 not in xyz_cb: print('XYZ not defined for ' + str(i+1))
        print ('Warning! Something went wrong - len of cbxyz != seqlen!! ' + str(len(xyz_cb)) + ' ' +  str(L))
    cb_map = np.full((L, L), np.nan)
    for r1 in sorted(xyz_cb):
        (a, b, c) = xyz_cb[r1]
        for r2 in sorted(xyz_cb):
            (p, q, r) = xyz_cb[r2]
            cb_map[r1 - 1, r2 - 1] = sqrt((a-p)**2+(b-q)**2+(c-r)**2)
    return (total_valid_residues, cb_map, rnum_rnames)
# Helpers for metrics calculated using numpy scheme
def get_flattened(dmap):
  if dmap.ndim == 1:
    return dmap
  elif dmap.ndim == 2:
    return dmap[np.triu_indices_from(dmap, k=1)]
  else:
    assert False, "ERROR: the passes array has dimension not equal to 2 or 1!"

def get_separations(dmap):
  t_indices = np.triu_indices_from(dmap, k=1)
  separations = np.abs(t_indices[0] - t_indices[1])
  return separations

# return a 1D boolean array indicating where the sequence separation in the
# upper triangle meets the threshold comparison
def get_sep_thresh_b_indices(dmap, thresh, comparator):
  assert comparator in {'gt', 'lt', 'ge', 'le'}, "ERROR: Unknown comparator for thresholding!"
  dmap_flat = get_flattened(dmap)
  separations = get_separations(dmap)
  if comparator == 'gt':
    threshed = separations > thresh
  elif comparator == 'lt':
    threshed = separations < thresh
  elif comparator == 'ge':
    threshed = separations >= thresh
  elif comparator == 'le':
    threshed = separations <= thresh

  return threshed

# return a 1D boolean array indicating where the distance in the
# upper triangle meets the threshold comparison
def get_dist_thresh_b_indices(dmap, thresh, comparator):
  assert comparator in {'gt', 'lt', 'ge', 'le'}, "ERROR: Unknown comparator for thresholding!"
  dmap_flat = get_flattened(dmap)
  if comparator == 'gt':
    threshed = dmap_flat > thresh
  elif comparator == 'lt':
    threshed = dmap_flat < thresh
  elif comparator == 'ge':
    threshed = dmap_flat >= thresh
  elif comparator == 'le':
    threshed = dmap_flat <= thresh
  return threshed

# Calculate lDDT using numpy scheme
def get_LDDT(true_map, pred_map, R=15, sep_thresh=-1, T_set=[0.5, 1, 2, 4], precision=4):
    '''
    Mariani V, Biasini M, Barbato A, Schwede T.
    lDDT: a local superposition-free score for comparing protein structures and models using distance difference tests.
    Bioinformatics. 2013 Nov 1;29(21):2722-8.
    doi: 10.1093/bioinformatics/btt473.
    Epub 2013 Aug 27.
    PMID: 23986568; PMCID: PMC3799472.
    '''
    
    # Helper for number preserved in a threshold
    def get_n_preserved(ref_flat, mod_flat, thresh):
        err = np.abs(ref_flat - mod_flat)
        n_preserved = (err < thresh).sum()
        return n_preserved
    
    # flatten upper triangles
    true_flat_map = get_flattened(true_map)
    pred_flat_map = get_flattened(pred_map)
    
    # Find set L
    S_thresh_indices = get_sep_thresh_b_indices(true_map, sep_thresh, 'gt')
    R_thresh_indices = get_dist_thresh_b_indices(true_flat_map, R, 'lt')
    
    L_indices = S_thresh_indices & R_thresh_indices
    
    true_flat_in_L = true_flat_map[L_indices]
    pred_flat_in_L = pred_flat_map[L_indices]
    
    # Number of pairs in L
    L_n = L_indices.sum()
    
    # Calculated lDDT
    preserved_fractions = []
    for _thresh in T_set:
        _n_preserved = get_n_preserved(true_flat_in_L, pred_flat_in_L, _thresh)
        _f_preserved = _n_preserved / L_n
        preserved_fractions.append(_f_preserved)
    
    lDDT = np.mean(preserved_fractions)
    if precision > 0:
        lDDT = round(lDDT, precision)
    return lDDT

# Function to find distance
def distance(x,y): 
    
    x1 = x[0]
    x2 = y[0]
    y1 = x[1]
    y2 = y[1]
    z1 = x[2]
    z2 = y[2]

    d = sqrt(pow(x2 - x1, 2) +
                pow(y2 - y1, 2) +
                pow(z2 - z1, 2)* 1.0)
    return d


def smith_waterman(string1, string2, match_score=2, mismatch_penalty=-1, gap_penalty=-1):
    rows = len(string1) + 1
    cols = len(string2) + 1
    matrix = [[0 for _ in range(cols)] for _ in range(rows)]
    max_score = 0
    max_pos = None

    for i in range(1, rows):
        for j in range(1, cols):
            if string1[i - 1] == string2[j - 1]:
                diagonal_score = matrix[i - 1][j - 1] + match_score
            else:
                diagonal_score = matrix[i - 1][j - 1] + mismatch_penalty
            up_score = matrix[i - 1][j] + gap_penalty
            left_score = matrix[i][j - 1] + gap_penalty
            matrix[i][j] = max(0, diagonal_score, up_score, left_score)

            if matrix[i][j] > max_score:
                max_score = matrix[i][j]
                max_pos = (i, j)

    if max_pos is None:
        return None 

    align1 = []#""
    align2 = []#""
    current_i, current_j = max_pos
    start_index_string1 = current_i - 1

    while current_i > 0 and current_j > 0 and matrix[current_i][current_j] > 0:
        if string1[current_i - 1] == string2[current_j - 1]:
            align1.append(current_i-1) #= string1[current_i - 1] + align1
            align2.append(current_j - 1) #= string2[current_j - 1] + align2
            current_i -= 1
            current_j -= 1
            start_index_string1 = current_i 
        elif matrix[current_i][current_j] == matrix[current_i - 1][current_j] + gap_penalty:
            # align1 = string1[current_i - 1] + align1
            # align2 = "-" + align2
            current_i -= 1
            start_index_string1 = current_i 
        else:
            # align1 = "-" + align1
            # align2 = string2[current_j - 1] + align2
            current_j -= 1

    return start_index_string1, align1, align2

def get_common_index(pdb1,pdb2):
    parser = Bio.PDB.PDBParser()
    structure1 = parser.get_structure("ref",pdb1)
    structure2 = parser.get_structure("query",pdb2)
    seq1 = [x.get_resname() for x in structure1.get_residues()]
    seq2 = [x.get_resname() for x in structure2.get_residues()]
    return smith_waterman(seq1,seq2)

def get_alpha_carbon_coordinates(name,pdb_file):
    parser = Bio.PDB.PDBParser()
    structure = parser.get_structure(name,pdb_file)
    coordinates = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.has_id("CA"):
                    coordinates.append(residue["CA"].get_coord())

    return coordinates

#Calculates both the per residue lDDT and deviation which can be used to calculate global lDDT and global RMSD
def calc_scores(name,ref_file,pred_file,T_set=[0.5, 1, 2, 4]):
    ref_ca = get_alpha_carbon_coordinates(name,ref_file)
    pred_ca = get_alpha_carbon_coordinates(name,pred_file)
    if len(ref_ca) != len(pred_ca):
        print(f"Ref ca count: {len(ref_ca)}, Pred ca count: {len(pred_ca)}")
        print("Mismatching residues predicted, returning null")
        return "",""
    
    # get_common_index(ref_file,pred_file)
    length = len(ref_ca)

    deviation = []
    for p in range(len(ref_ca)):
        deviation.append(distance(ref_ca[p],pred_ca[p]))

    ref_ca_dist = []
    pred_ca_dist = []
    for L in range(len(ref_ca)):
        x = ref_ca[L]
        x_pred = pred_ca[L]
        for sub in range(len(ref_ca)):
            y = ref_ca[sub]
            y_pred = ref_ca[sub]
            ref_ca_dist.append(distance(x,y))
            pred_ca_dist.append(distance(x_pred,y_pred))
    
    ref_ca_sum = np.array(ref_ca_dist)
    pred_ca_sum = np.array(pred_ca_dist)
    ref_ca_sum = ref_ca_sum.reshape(length,length)
    pred_ca_sum = pred_ca_sum.reshape(length,length)

    #Now we have to calculate the rate of true/pred at each threshold for every residue
    lDDT = []
    for x in range(len(ref_ca)):
        Preserved = []
        for T in T_set:
            Frac = np.nan_to_num((ref_ca_sum[x] < T).sum() / (pred_ca_sum[x] < T).sum(),posinf=0)
            Preserved.append(Frac)
        lddt = np.mean(Preserved)
        lDDT.append(lddt)
    
    
    return deviation,lDDT

#Extracts the plDDT from an alphafold PDB file
def extractbscore(name,file):
    parser = Bio.PDB.PDBParser()
    structure = parser.get_structure(name, file)

    residues = [r for r in structure.get_residues() if r.get_id()[0] == " "]
    lddt = []
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:                                  
                    N = atom.get_name()
                    I = atom.get_id()
                    Y = atom.get_coord()                
                    V = atom.get_vector()
                    O = atom.get_occupancy()
                    B = atom.get_bfactor()     
                    if N == "CA":
                        lddt.append(B) 

    return lddt

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

def count_clashes(structure, clash_cutoff=0.63):
    #Setup initial dataframe
    data = []#pd.DataFrame()
    # Set what we count as a clash for each pair of atoms
    clash_cutoffs = {i + "_" + j: (clash_cutoff * (atom_radii[i] + atom_radii[j])) for i in atom_radii for j in atom_radii}
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
            residue_1coord = residue_1['CA'].get_vector()
            # Find atoms that could be clashing
            kdt_search = kdt.search(np.array(atom_1.coord, dtype="d"), max(clash_cutoffs.values()))
            # Get index and distance of potential clashes
            potential_clash = [(a.index, a.radius) for a in kdt_search]
            for ix, atom_distance in potential_clash:
                atom_2 = atoms[ix]
                residue_2 = atom_2.get_parent()
                residue_2coord = residue_2['CA'].get_vector()
                res_distance = np.linalg.norm(residue_1coord - residue_2coord)
                # Exclude clashes from atoms in the same residue
                if atom_1.parent.id == atom_2.parent.id:
                    continue
                # Exclude clashes from peptide bonds
                elif (atom_2.name == "C" and atom_1.name == "N") or (atom_2.name == "N" and atom_1.name == "C"):
                    continue
                # Exclude clashes from disulphide bridges
                elif (atom_2.name == "SG" and atom_1.name == "SG") and atom_distance > 1.88:
                    continue
                if atom_distance < clash_cutoffs[atom_2.element + "_" + atom_1.element]:
                    clashes.append((atom_1, atom_2))
                    #Clashes should be a pd dataframe with the columns: Name, Residue 1, Residue 1 pos, Atom 1, Atom 1 XYZ, Residue 2, Residue 2 pos, Atom 2, Atom 2 XYZ, Residue distance, Atom distance
                    df = pd.DataFrame({'Name': [name],
                                    'Res1': [residue_1.get_resname()],
                                    'Res1pos': [residue_1.get_full_id()[3][1]],
                                    'Atom1': [atom_1.get_name()],
                                    'Atom1_x': [atom_1.coord[0]],
                                    'Atom1_y': [atom_1.coord[1]],
                                    'Atom1_z': [atom_1.coord[2]],
                                    'Res2': [residue_2.get_resname()],
                                    'Res2pos': [residue_2.get_full_id()[3][1]],
                                    'Atom2': [atom_2.get_name()],
                                    'Atom2_x': [atom_2.coord[0]],
                                    'Atom2_y': [atom_2.coord[1]],
                                    'Atom2_z': [atom_2.coord[2]],
                                    'ResDistance': [res_distance],
                                    'AtomDistance': [atom_distance]})
                    data.append(df)
        except:
            continue
    if data:
        clash_df = pd.concat(data)
        return clash_df
    else:
        return []

__all__ = ["ShrakeRupley"]

_ENTITY_HIERARCHY = {
    "A": 0,
    "R": 1,
    "C": 2,
    "M": 3,
    "S": 4,
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

def extractDockQ(file):
    
    name = os.path.basename(file)[slice(0,3,1)]
    with open(file,'r') as f:
        data = json.load(f)

    map = data["best_mapping"]
    df = {'Name': name, 'GlobalDockQ': data['GlobalDockQ']}
    


    # mhca_pep = data["best_result"]["AC"]
    # if ("BC" in data["best_result"]):
    #     print("exists")
    #     mhcb_pep = data["best_result"]["BC"]
    #     mhcb_tcra = data["best_result"]["BD"]
    #     mhcb_tcra = data["best_result"]["BE"]
    # else:
    #     print("does not exist")
    #     mhcb_pep = {'DockQ':'NA','fnat': 'NA','clashes': 'NA'}
    #     mhcb_tcra = {'DockQ':'NA','fnat': 'NA','clashes': 'NA'}
    #     mhcb_tcrb = {'DockQ':'NA','fnat': 'NA','clashes': 'NA'}
    # # mhcb_pep = data["best_result"]["BC"]
    # mhca_tcra = data["best_result"]["AD"]
    # mhca_tcrb = data["best_result"]["AE"]
    # pep_tcra = data["best_result"]["CD"]
    # pep_tcrb = data["best_result"]["CE"]

    # df = {'Name': name, 
    #                         'GlobalDockQ': data['GlobalDockQ'],
    #                         'MHCA_Pep_DockQ': mhca_pep['DockQ'],
    #                         'MHCA_Pep_fnat': mhca_pep['fnat'],
    #                         'MHCA_Pep_clashes': mhca_pep['clashes'],
    #                         'MHCB_Pep_DockQ': mhcb_pep['DockQ'],
    #                         'MHCB_Pep_fnat': mhcb_pep['fnat'],
    #                         'MHCB_Pep_clashes': mhcb_pep['clashes'],
    #                         'MHCA_TCRA_DockQ': mhca_tcra['DockQ'],
    #                         'MHCA_TCRA_fnat': mhca_tcra['fnat'],
    #                         'MHCA_TCRA_clashes': mhca_tcra['clashes'],
    #                         'MHCB_TCRA_DockQ': mhca_pep['DockQ'],
    #                         'MHCB_TCRA_fnat': mhca_pep['fnat'],
    #                         'MHCB_TCRA_clashes': mhcb_tcra['clashes'],
    #                         'MHCA_TCRB_DockQ': mhca_tcrb['DockQ'],
    #                         'MHCA_TCRB_fnat': mhca_tcrb['fnat'],
    #                         'MHCA_TCRB_clashes': mhca_tcrb['clashes'],
    #                         'MHCB_TCRB_DockQ': mhcb_tcrb['DockQ'],
    #                         'MHCB_TCRB_fnat': mhcb_tcrb['fnat'],
    #                         'MHCB_TCRB_clashes': mhcb_tcrb['clashes'],
    #                         'Pep_TCRA_DockQ': pep_tcra['DockQ'],
    #                         'Pep_TCRA_fnat': pep_tcra['fnat'],
    #                         'Pep_TCRA_clashes': pep_tcra['clashes'],
    #                         'Pep_TCRB_DockQ': pep_tcrb['DockQ'],
    #                         'Pep_TCRB_fnat': pep_tcrb['fnat'],
    #                         'Pep_TCRB_clashes': pep_tcrb['clashes']}
    return df


def pdb_to_mmcif(pdb_file, mmcif_file):
    parser = PDBParser()
    structure = parser.get_structure("structure", pdb_file)

    io = MMCIFIO()
    io.set_structure(structure)
    io.save(mmcif_file)

def get_bonded_atom(atom):
    residue = atom.get_parent()
    close = 100
    for neighbor in residue.get_atoms():
        if 'H' in neighbor.name:
            continue
        dist = atom - neighbor
        if neighbor != atom :#and dist < 2.2:  # Check distance threshold
            if dist < close:
                 bonded = neighbor
    return bonded

def is_hydrogen_bond(atom1, atom2):
    #Use these as the maximum cutoffs just to quickly assess
    distance_cutoff = 4#3.5  # Angstroms
    angle_cutoff = 120  # Degrees
    #Save the donor, acceptor and hydrogen atoms, we need to get the bonded atom to the hydrogen
    if 'H' not in atom1.name and 'H' in atom2.name:
        acceptor = atom1
        hydrogen = atom2
        donor = get_bonded_atom(atom2)
        
    elif 'H' in atom1.name and 'H' not in atom2.name:
        donor = get_bonded_atom(atom1)
        hydrogen = atom1
        acceptor = atom2
    else:
        return False,nan,nan,nan,nan,nan
    # if "C" in acceptor.name or "C" in donor.name: # in ['N','CA','CB']: #Maybe  carbon can act but need charges 
    #      return False,nan,nan,nan,nan,nan
    if not donor:
            print(hydrogen, acceptor)
            print("no donor")
    ## Get rid on non-negative charge 
    if donor.get_bfactor() > 0 or acceptor.get_bfactor() >= 0:
        return False,nan,nan,nan,nan,nan
    #First test if within max distance
    distance = donor - acceptor
    if distance > distance_cutoff:
        return False,nan,nan,nan,nan,nan
    #save the vectors making up the two bonds around the hydrogen ()
    vector_dh = hydrogen.coord - donor.coord
    vector_ha = acceptor.coord - hydrogen.coord
    #Get the angle of the vectors
    angle = np.degrees(np.arccos(np.dot(vector_dh, vector_ha) / (np.linalg.norm(vector_dh) * np.linalg.norm(vector_ha))))
    
    #If angle is greater than the cutoff then we should save (don't know if goes over 180?)
    if angle > angle_cutoff:
        #We should return a bunch of values
        return True,donor,acceptor,hydrogen,distance,angle
    else:
        return False,nan,nan,nan,nan,nan

def get_Hydrogen_Bonds(struct):
    hydrogen_bonds = []

    residues = [x for x in struct.get_residues() if x.get_resname() != "SOL"]
    atoms = [x for x in struct.get_atoms() ]#if x.element in atom_radii]
    # residues = list(filter(lambda x: is_aa(x), model))
    max = 0
    total_res = 0
    #################Begin for loop############################
    for i in range(len(residues)):
        for j in range(len(residues)):
            # if (j > len(residues) - i - 1): #Do not repeat residue comparisons
            #             continue
            
            res1 = residues[i]
            res2 = residues[j]
            chain1 = res1.get_parent().get_id()
            chain2 = res2.get_parent().get_id()
            # if (res1.get_resname() == "SOL" or res2.get_resname() == "SOL"):
            #             continue
            if chain1 == chain2:
                continue
            if res1 == res2:
                continue
            if (res1['CA'] - res2['CA']) > 10:
                continue
            atoms1 = [p for p in res1.get_atoms()]
            atoms2 = [p for p in res2.get_atoms()]
            for atom1 in atoms1:
                for atom2 in atoms2:
                    if 'H' not in atom1.name and 'H' not in atom2.name:
                        continue
                    if 'H' in atom1.name and 'H' in atom2.name:
                        continue
                    is_bond,donor,acceptor,hydrogen,distance,angle = is_hydrogen_bond(atom1, atom2) #Save the output of is_hydrogen_bond()
                    if is_bond:
                        # print(res1,res2,donor, acceptor)
                        don = donor.get_parent()
                        if don == res1:
                            donor_res = 1
                            acceptor_res = 2
                        else:
                            donor_res = 2
                            acceptor_res = 1
                        dict = {'Res1': res1.get_resname(), 'Res1_num': i,'Res1_pos': res1.get_full_id()[3][1],'Chain1': res1.get_parent().get_id(),
                                'Res2': res2.get_resname(), 'Res2_num': j,'Res2_pos': res2.get_full_id()[3][1],'Chain2': res2.get_parent().get_id(),
                                'Donor_res': donor_res, 'Acceptor_res': acceptor_res,
                                'Donor': donor.get_name(),'Acceptor': acceptor.get_name(),'Hydrogen': hydrogen.get_name(),
                                'Bond_length': distance, 'Bond_angle': angle,
                                'Donor_Charge': donor.get_bfactor(),'Donor_Radius': donor.get_occupancy(),
                                'Hydrogen_Charge': hydrogen.get_bfactor(),'Hydrogen_Radius': hydrogen.get_occupancy(),
                                'Acceptor_Charge': acceptor.get_bfactor(),'Acceptor_Radius': acceptor.get_occupancy(),
                                'Donor_x': donor.get_coord()[0],'Donor_y': donor.get_coord()[1], 'Donor_z': donor.get_coord()[2],
                                'Acceptor_x': acceptor.get_coord()[0],'Acceptor_y': acceptor.get_coord()[1], 'Acceptor_z': acceptor.get_coord()[2],
                                'Hydrogen_x': hydrogen.get_coord()[0],'Hydrogen_y': hydrogen.get_coord()[1], 'Hydrogen_z': hydrogen.get_coord()[2]}
                        # print(dict)
                        hydrogen_bonds.append(dict)
    df = pd.DataFrame(hydrogen_bonds)
    return df

def center_of_mass(entity, geometric=False):
    # Copyright (C) 2010, Joao Rodrigues (anaryin@gmail.com)
    # This code is part of the Biopython distribution and governed by its
    # license.  Please see the LICENSE file that should have been included
    # as part of this package.
    """
    Returns gravitic [default] or geometric center of mass of an Entity.
    Geometric assumes all masses are equal (geometric=True)
    """

    # Structure, Model, Chain, Residue
    # if isinstance(entity, Entity.Entity):
    # atom_list = entity.get_atoms()
    # List of Atoms
    # elif hasattr(entity, "__iter__") and [x for x in entity if x.level == "A"]:
    atom_list = entity
    # else:  # Some other weirdo object
    #     raise ValueError(
    #         "Center of Mass can only be calculated from the following objects:\n"
    #         "Structure, Model, Chain, Residue, list of Atoms."
    #     )

    masses = []
    positions = [[], [], []]  # [ [X1, X2, ..] , [Y1, Y2, ...] , [Z1, Z2, ...] ]

    for atom in atom_list:
        masses.append(atom.mass)

        for i, coord in enumerate(atom.coord.tolist()):
            positions[i].append(coord)

    # If there is a single atom with undefined mass complain loudly.
    if "ukn" in set(masses) and not geometric:
        raise ValueError(
            "Some Atoms don't have an element assigned.\n"
            "Try adding them manually or calculate the geometrical center of mass instead"
        )
    if (sum(masses) == 0):
        print("Warning: All atoms have zero mass.")
        return False
    if geometric:
        return [sum(coord_list) / len(masses) for coord_list in positions]
    else:
        w_pos = [[], [], []]
        for atom_index, atom_mass in enumerate(masses):
            w_pos[0].append(positions[0][atom_index] * atom_mass)
            w_pos[1].append(positions[1][atom_index] * atom_mass)
            w_pos[2].append(positions[2][atom_index] * atom_mass)

        return [sum(coord_list) / sum(masses) for coord_list in w_pos]

def fetch_atoms(model, selection="A", atom_bounds=[1, 180]):
    """
    Function to fetch atoms from the defined "atom_bound" in "selection" of the "model"
    """
    selection = selection.upper()
    if not isinstance(atom_bounds, (list, tuple)) and len(atom_bounds) > 0:
        raise ValueError("expected non-empty list or tuple, got {}".format(atom_bounds))
    # making sure its a list of lists
    if not isinstance(atom_bounds[0], (list, tuple)):
        atom_bounds = [atom_bounds]
    if not all([len(b) == 2 for b in atom_bounds]):
        raise ValueError(
            "All bounds must be providing one upper "
            "and one lower bound, got {}".format(atom_bounds)
        )
    if not isinstance(selection, (tuple, list)):
        selection = [selection]
    result = []
    for sel in selection:
        for ref_res in model["%s" % sel]:
            resid = ref_res.get_id()[1]
            in_bounds = False
            for bounds in atom_bounds:
                in_bounds |= bounds[0] <= resid and resid <= bounds[1]
            if in_bounds:
                result.append(ref_res["CA"])
    return result

def fetch_entity(model, fetch_atoms=True, selection="A", res_ids=range(1, 180)):
    """
    Function to fetch atoms/resids from the defined "resid_bounds" in "selection" of the "model"
    """
    selection = selection.upper()
    # fetch atoms
    if fetch_atoms is True:
        result = []
        for sel in selection:
            for sample_res in model["%s" % sel]:
                resid = sample_res.get_id()[1]
                if resid in res_ids:
                    result.append(sample_res["CA"])
    # fetch_residues_indeces
    elif fetch_atoms is False:
        result = []
        for sel in selection:
            for sample_res in model["%s" % sel]:
                resid = sample_res.get_id()[1]
                if resid in res_ids:
                    result.append(resid)
    return result

def apply_transformation_to_atoms(model, rotmat, transvec):
    """
    Function to translate/rotate the model by the defined translation vector and rotation matrix
    """
    for chain in model:
        for res in chain:
            for atom in res:
                atom.transform(rotmat, transvec)

def is_chain_in_pdb(pdbid, input_chain_id):
    """
    Function to check if the input_chain_id is in pdb_chain_ids_list
    """

    pdb_chain_ids_list = []
    structure = PDBParser().get_structure("%s" % pdbid, "%s.pdb" % pdbid)
    model = structure[0]
    for chain in model:
        pdb_chain_ids_list.append(str(chain.get_id()))
    return input_chain_id in pdb_chain_ids_list

def number_of_chains(pdbid):

    pdb_chain_ids_list = []
    structure = PDBParser().get_structure("%s" % pdbid, "%s.pdb" % pdbid)
    model = structure[0]
    counter = 0
    for chain in model:
        counter += 1
    return counter

def str2bool(v):
    if isinstance(v, bool):
        return v
    if v is None:
        return True
    if v.lower() in ("TRUE", "True", "yes", "true", "t", "y", "1"):
        return True
    elif v.lower() in ("FALSE", "False", "no", "false", "f", "n", "0"):
        return False
    else:
        raise argparse.ArgumentTypeError("Boolean value expected.")

# reffile1 = "./FoldingModels/ref1.pdb"
# reffile2 = "./FoldingModels/ref2.pdb"

def tcr_mhci_geometrical_parameters(
    pdbid,
    mhc_a="A",
    tcr_a="D",
    tcr_b="E",
    persist_structure=True,
    mhc_a_init=1,
    mhc_a_final=179,
    tcr_a_init=1,
    tcr_a_final=109,
    tcr_b_init=1,
    tcr_b_final=116,
):
    """
    ARGUMENTS:
    pdbid = PDB-ID or name of the input PB file (str)
    mhc_a = ID of MHC alpha chain (str)
    tcr_a = ID of TCR alpha chain (str)
    tcr_b = ID of TCR beta chain (str)
    persist_structure = If true, save the processed structure as PDB-file (Boolean)
    """
    #################################################################
    # Define residues range to align and center of mass calculations#
    #################################################################
    mhc_a_resids = range(mhc_a_init, mhc_a_final + 1)
    tcr_a_resids = range(tcr_a_init, tcr_a_final + 1)
    tcr_b_resids = range(tcr_b_init, tcr_b_final + 1)

    ########################################################################################################
    # Import structure, align to reference, and calculate center of mass of CA atoms in MHCI binding groove#
    ########################################################################################################
    pdb_parser = Bio.PDB.PDBParser(QUIET=True)
    ref_structure = pdb_parser.get_structure("reference", reffile1)
    sample_structure = pdb_parser.get_structure("sample", pdbid)
    # Use the first model in the pdb-files for alignment
    ref_model = ref_structure[0]
    sample_model = sample_structure[0]
    # Iterate of all residues in each model in order to define proper atoms
    # Sample structure
    sample_resids = fetch_entity(
        sample_model, fetch_atoms=False, selection=mhc_a, res_ids=mhc_a_resids
    )
    sample_atoms = fetch_entity(
        sample_model, fetch_atoms=True, selection=mhc_a, res_ids=sample_resids
    )
    # Reference structure
    # ref_atoms = fetch_entity(
    #     ref_model, fetch_atoms=True, selection=mhc_a, res_ids=sample_resids
    # )

    # Initiate the superimposer:
    # super_imposer = Bio.PDB.Superimposer()
    # super_imposer.set_atoms(ref_atoms, sample_atoms)
    # super_imposer.apply(sample_model.get_atoms())
    # Calculate CoM of MHCI binding groove
    mhci_com = center_of_mass(sample_atoms, geometric=False)
    # Calculate CoM of vTCR
    tcr_atoms_for_com = fetch_entity(
        sample_model, fetch_atoms=True, selection=tcr_a, res_ids=tcr_a_resids
    )
    tcr_atoms_for_com += fetch_entity(
        sample_model, fetch_atoms=True, selection=tcr_b, res_ids=tcr_b_resids
    )
    vtcr_com = center_of_mass(tcr_atoms_for_com, geometric=True)
    if (not vtcr_com) or (not mhci_com):
        print("Warning: Center of mass calculation failed. "
              "Check if the input PDB file contains atoms with defined masses.")
        return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan


    ###################################
    # Calculate geometrical parameters#
    ###################################
    dx, dy, dz = np.subtract(vtcr_com, mhci_com)
    r = np.sqrt(np.sum(np.square(np.subtract(vtcr_com, mhci_com))))
    theta = np.degrees(np.arctan2(dy, dx))
    phi = np.degrees(np.arccos(dz / r))
    # print(
    #     "The Geomitrical parameters: r = {:.2f}, "
    #     "theta = {:.2f}, phi = {:.2f}".format(r, theta, phi)
    # )
    return r, theta, phi, dx, dy, dz

def tcr_mhcii_geometrical_parameters(
    pdbid,
    mhc_a="A",
    mhc_b="B",
    tcr_a="D",
    tcr_b="E",
    persist_structure=True,
    mhc_a_init=1,
    mhc_a_final=80,
    mhc_b_init=1,
    mhc_b_final=90,
    tcr_a_init=1,
    tcr_a_final=109,
    tcr_b_init=1,
    tcr_b_final=116,
):
    """
    ARGUMENTS:
    pdbid = PDB-ID or name of the input PB file (str)
    mhc_a = ID of MHC alpha chain (str)
    mhc_b = ID of MHC beta chain (str)
    tcr_a = ID of TCR alpha chain (str)
    tcr_b = ID of TCR beta chain (str)
    persist_structure = If true, save the processed structure as PDB-file (Boolean)
    """
    #################################################################
    # Define residues range to align and center of mass calculations#
    #################################################################
    mhc_a_resids = range(mhc_a_init, mhc_a_final + 1)
    mhc_b_resids = range(mhc_b_init, mhc_b_final + 1)
    tcr_a_resids = range(tcr_a_init, tcr_a_final + 1)
    tcr_b_resids = range(tcr_b_init, tcr_b_final + 1)

    #########################################################################################################
    # Import structure, align to reference, and calculate center of mass of CA atoms in MHCII binding groove#
    #########################################################################################################
    pdb_parser = Bio.PDB.PDBParser(QUIET=True)
    # ref_structure = pdb_parser.get_structure("reference", reffile2)
    sample_structure = pdb_parser.get_structure("sample", pdbid)
    # ref_model = ref_structure[0]
    sample_model = sample_structure[0]
    # Iterate of all residues in each model in order to define proper atoms
    # Sample structure
    sample_resids_a = fetch_entity(
        sample_model, fetch_atoms=False, selection=mhc_a, res_ids=mhc_a_resids
    )
    sample_resids_b = fetch_entity(
        sample_model, fetch_atoms=False, selection=mhc_b, res_ids=mhc_b_resids
    )

    sample_atoms = fetch_entity(sample_model, selection=mhc_a, res_ids=sample_resids_a)
    sample_atoms += fetch_entity(sample_model, selection=mhc_b, res_ids=sample_resids_b)
    # Reference structure
    # ref_atoms = fetch_entity(ref_model, selection=mhc_a, res_ids=sample_resids_a)
    # ref_atoms += fetch_entity(ref_model, selection=mhc_b, res_ids=sample_resids_b)

    # Initiate the superimposer:
    # super_imposer = Bio.PDB.Superimposer()
    # super_imposer.set_atoms(ref_atoms, sample_atoms)
    # super_imposer.apply(sample_model.get_atoms())
    # Calculate CoM of MHCII binding groove
    mhcii_com = center_of_mass(sample_atoms, geometric=True)
    # Calculate CoM of vTCR
    tcr_atoms_for_com = fetch_entity(
        sample_model, selection=tcr_a, res_ids=tcr_a_resids
    )
    tcr_atoms_for_com += fetch_entity(
        sample_model, selection=tcr_b, res_ids=tcr_b_resids
    )
    vtcr_com = center_of_mass(tcr_atoms_for_com, geometric=True)
    # print("MHC-CoM: ", [round(x, 2) for x in mhcii_com])
    # print("vTCR-CoM: ", [round(x, 2) for x in vtcr_com])
    #######################
    # Save final structure#
    #######################
    # if persist_structure:
    #     ali_structure = add_com_to_pdb(mhcii_com, vtcr_com, sample_model)
    #     io = Bio.PDB.PDBIO()
    #     io.set_structure(sample_structure)
    #     io.save("%s_aligned.pdb" % pdbid)

    ###################################
    # Calculate geometrical parameters#
    ###################################
    dx, dy, dz = np.subtract(vtcr_com, mhcii_com)
    r = np.sqrt(np.sum(np.square(np.subtract(vtcr_com, mhcii_com))))
    theta = np.degrees(np.arctan2(dy, dx))
    phi = np.degrees(np.arccos(dz / r))
    # print(
    #     "The Geomitrical parameters: r = {:.2f}, "
    #     "theta = {:.2f}, phi = {:.2f}".format(r, theta, phi)
    # )
    return r, theta, phi, dx, dy, dz

k_standard = 8.9875517873681764e9
angstrom_m = 10e10
electron = 1.60217646e-19
avogadro = 6.022e23
start_time = time.time()

def computeConstants():
    atoms = ['C','H','N','O','P','S']
    length = len(atoms)
    r_eqm_X = np.array([2.00, 1.00, 1.75, 1.60, 2.10, 2.00])
    eps_X = np.array([0.15, 0.02,0.16,0.20,0.20,0.20])
    r_eqm_X = np.repeat(r_eqm_X,length).reshape(length,length).swapaxes(0,1)
    eps_X = np.repeat(eps_X,length).reshape(length,length).swapaxes(0,1)
    eps_XY = np.sqrt(eps_X*eps_X.T)
    r_eqm_XY = (r_eqm_X + r_eqm_X.T)/2.0
    A = eps_XY*(r_eqm_XY**12)
    B = 2 * eps_XY * (r_eqm_XY**6)
    return A, B


def atom_atomPotential(atom1,atom2,coord1,coord2,A,B,charge1,charge2):
    atoms = ['C','H','N','O','P','S']
    for num, val in enumerate(atoms):
        if val == atom1:
            atom1 = num
        if val == atom2:
            atom2 = num
    dist = coord1 - coord2
    dist = np.dot(dist,dist)
    rootdist = pow(dist,0.5)
    k = k_standard*angstrom_m*electron**2 #adjust Coulomb constant for units
    LJ = A[atom1,atom2]/pow(dist,6) - B[atom1,atom2]/pow(dist,3) #Is this the right calculation? for LJ?
    Coulomb = k*(charge1*charge2)/rootdist
    Coulomb = (Coulomb/1000)*avogadro
    return LJ,Coulomb	

#For iterMol to calculated LJ and columb we need a list of all coordinates, a list of all atoms, a list of all bonded atoms, 
#       a list of charges (in the b value)
def getBondedAtoms(atoms):
    cutoff = 2.0  # Distance cutoff for bond
    ns = NeighborSearch(atoms)
    bonds = []
    for atom in atoms:
        bonded_atoms = []
        for neighbor in ns.search(atom.coord, cutoff, level='A'):
            if neighbor != atom:  # Exclude the atom itself
                bonded_atoms.append(neighbor.element)
        bonds.append(bonded_atoms)
    return bonds

def getLists(residue):
    atoms = [x for x in residue.get_atoms()]
    # atoms = [item for sublist in resatoms for item in sublist]
    #now we can grab things like the coordinates, b value, name, bonds
    coordinates = np.array([x.get_coord() for x in atoms])
    atom_names = [x.element for x in atoms]
    charges = [x.get_bfactor() for x in atoms]
    bonds = getBondedAtoms(atoms)

    return coordinates, atom_names, bonds, charges

def residue_Potential(res1, res2, A,B):
    res1_coords, res1_names, res1_bonds, res1_charges = getLists(res1)	
    res2_coords, res2_names, res2_bonds, res2_charges = getLists(res2)	
    coordlist = np.concatenate((res1_coords,res2_coords),axis=0)
    atomlist = res1_names + res2_names
    bonds = res1_bonds + res2_bonds
    charges = res1_charges + res2_charges
    totalLJ = float()
    totalCoul = float()
    for i in range(len(atomlist)):
        for j in [x for x in range(len(atomlist)-i) if x != i]:
            if (j in bonds[i]) or (i in bonds[j]):
                pass
            else:
                #Set cutoff of 1.4 for coloumb and 3 for LJ (2.5x sigma=3.4)
                x = coordlist[i,0] - coordlist[j,0]
                y = coordlist[i,1] - coordlist[j,1]
                z = coordlist[i,2] - coordlist[j,2]
                if x > 10 or y > 10 or z > 10:
                    continue
                temp,Coul = atom_atomPotential(atomlist[i],atomlist[j],coordlist[i,:],coordlist[j,:],A,B,charges[i],charges[j])
                totalLJ += temp
                totalCoul += Coul
    return totalLJ, totalCoul

def iterMol(name,pdbfile):
    energyVals = list()
    parser = Bio.PDB.PDBParser(QUIET = False)
    structure = parser.get_structure(name, pdbfile)
    residues = [x for x in structure.get_residues() if x.get_resname() != "SOL"]
    A, B = computeConstants()
    totalLJ = float()
    totalCoul = float()
    residue1 = []
    residue2 = []
    LJ = []
    Coul = []
    dists = []
    for i in range(len(residues)):
        counter = float()
        res1 = residues[i]
        for j in [x for x in range(len(residues)) if x != i]:
            res2 = residues[j]
            dist = np.linalg.norm(res1["CA"].get_coord() - res2["CA"].get_coord())
            if (dist) > 20: #Do not compare far away residues
                continue 
            lj, coul = residue_Potential(res1,res2,A,B)
            residue1.append(i)
            residue2.append(j)
            LJ.append(lj)
            Coul.append(coul)
            dists.append(dist)   
    df = pd.DataFrame({'Res1': residue1, 'Res2':residue2, 'LJ': LJ,'Coulomb': Coul, 'Distance': dists})
    return df

def renameChains(pdbfile,name,sequences):
    """
    Function to rename chains in a PDB file based on the sequence file.
    If no sequence file is provided, it will rename chains to A, B, C, etc.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", pdbfile)
    model = structure[0]
    chains = []
    for chain in model:
        chains.append(chain.get_id())
    
    row = sequences.loc[name]
    mhca = row.iloc[0]
    mhcb = row.iloc[1]
    pep = row.iloc[2]
    tcra = row.iloc[3]
    tcrb = row.iloc[4]
    #Set up the renamer
    if "P" in chains:
        # print("Peptide found")
        renames = {
        "M": "A",
        # "B": "B",
        "P": "C",
        "A": "D",
        "B": "E",
        }
    else:
        renames = {
        "A": "A",
        "B": "B",
        "C": "C",
        "D": "D",
        "E": "E",
        }
        if len(chains) == 5:
            mhc2 = True
        elif (len(chains) == 4):
            mhc2 = False
            renames[chains[0]] = "A"
            renames[chains[1]] = "C"
            renames[chains[2]] = "D"
            renames[chains[3]] = "E"
        else:
            renames[chains[0]] = "A"
            renames[chains[1]] = "D"
            renames[chains[2]] = "E"
    
    for model in structure:
        for chain in model:
            old_name = chain.get_id()
            new_name = renames.get(old_name)
            if new_name:
                # print(f"renaming chain {old_name} to {new_name}")
                chain.id = new_name
            N = 1
            for residue in chain:
                residue.id = (residue.id[0], N, residue.id[2])
                N += 1
            
            
    
    io = PDBIO()
    io.set_structure(structure)
    io.save(pdbfile)
##alpha = Residue 50 - 85 (alpha chain) #apparently i arbitrarily named these alpha and beta loop but the TCRB interacts with the "alpha" loop and vice versa
##beta = Residue 138 - 180 (alpha chain)
#####MHC Class 2
##alpha = Residue 46 - 77 (alpha chain)
##beta = Residue 50 - 91 (beta chain)

#Write a function that takes a series of points and calculates a vector as the line of best fit
def Line(atoms):
    if len(atoms) < 2:
        print("Not enough atoms to define a line")
        return np.array([nan,nan,nan])
    points = []
    for atom in atoms:
        points.append(atom.coord.tolist())
    points = np.array(points)
    x = points[:,0]
    y = points[:,1]
    z = points[:,2]
    A = np.vstack([x, np.ones(len(x))]).T
    m_x, c_x = np.linalg.lstsq(A, y, rcond=None)[0]
    m_y, c_y = np.linalg.lstsq(A, z, rcond=None)[0]
    direction = np.array([1, m_x, m_y])
    direction = direction / np.linalg.norm(direction)
    point_on_line = np.array([0, c_x, c_y])
    return direction


#Positions about
#Model  cdr3a     cdr3b     cdr2a     cdr2b    cdr1a     cdr1b    mhca    mhcb
#ground  109-136   110-134   58-70    58-69    25-32     25-40    57-86    139-176
#galpha  91-97     92-101    48-53    47-53    23-30     24-30    50-85    139-176
#exp     93-97     94-102    49-54    48-56    24-31     25-31    50-85    138-175
#fake    91-97     94-102    49-54    48-56    24-31     25-31    50-85    138-175
#imm     113-118   113-120   69-75    68-74    45-51     45-51    75-109   162-199



def calcCOMAngles(struct,ta,ma):
    tcr = struct[0]
    chains = [x.id for x in tcr.get_chains()]
    cdr3b_COM = np.array(center_of_mass(fetch_entity(tcr, selection="E", res_ids=range(92+ta,101+ta)), geometric=True))
    cdr3a_COM = np.array(center_of_mass(fetch_entity(tcr, selection="D", res_ids=range(91+ta,97+ta)), geometric=True))
    cdr2b_COM = np.array(center_of_mass(fetch_entity(tcr, selection="E", res_ids=range(48+ta,56+ta)), geometric=True))
    cdr2a_COM = np.array(center_of_mass(fetch_entity(tcr, selection="D", res_ids=range(49+ta,54+ta)), geometric=True))
    cdr1b_COM = np.array(center_of_mass(fetch_entity(tcr, selection="E", res_ids=range(25+ta,31+ta)), geometric=True))
    cdr1a_COM = np.array(center_of_mass(fetch_entity(tcr, selection="D", res_ids=range(24+ta,31+ta)), geometric=True))
    mhc_alpha_COM = np.array(center_of_mass(fetch_entity(tcr, selection="A",res_ids=range(50+ma,85+ma)), geometric=True))
    if len(chains) >4:
        mhc_beta_COM = np.array(center_of_mass(fetch_entity(tcr, selection="B",res_ids=range(50+ma,85+ma)), geometric=True)) 
    else:
        mhc_beta_COM = np.array(center_of_mass(fetch_entity(tcr, selection="A",res_ids=range(139+ma,176+ma)), geometric=True)) 
    pep_COM = np.array(center_of_mass(fetch_entity(tcr, selection="C"), geometric=True))

    #calc relevant distances

    cdr3a_pep_dist = np.linalg.norm(cdr3a_COM - pep_COM)
    cdr3b_pep_dist = np.linalg.norm(cdr3b_COM - pep_COM)
    cdr2a_mhcb_dist = np.linalg.norm(cdr2a_COM - mhc_beta_COM)
    cdr2b_mhca_dist = np.linalg.norm(cdr2b_COM - mhc_alpha_COM)
    cdr1a_mhcb_dist = np.linalg.norm(cdr1a_COM - mhc_beta_COM)
    cdr1b_mhca_dist = np.linalg.norm(cdr1b_COM - mhc_alpha_COM)

    #angles
    # cdr3a_start = fetch_entity(tcr, selection="D", res_ids=[93])[0].get_coord()
    # cdr3a_end = fetch_entity(tcr, selection="D", res_ids=[109])[0].get_coord()
    # cdr3b_start = fetch_entity(tcr, selection="E", res_ids=[95])[0].get_coord()
    # cdr3b_end = fetch_entity(tcr, selection="E", res_ids=[109])[0].get_coord()
    # pep_start = fetch_entity(tcr, selection="C", res_ids=[1])[0].get_coord()
    # pep_end = fetch_entity(tcr, selection="C", res_ids=[len(fetch_entity(tcr, selection="C"))])[0].get_coord()

    cdr3b_line = Line(fetch_entity(tcr, selection="E", res_ids=range(95,116)))
    cdr3a_line = Line(fetch_entity(tcr, selection="D", res_ids=range(93,109)))
    pep_line = Line(fetch_entity(tcr, selection="C"))
    #Now calcualte angles between lines
    cdr3a_pep_angle = angle_between(cdr3a_line, pep_line)
    cdr3b_pep_angle = angle_between(cdr3b_line, pep_line)

    #now return all coms and angles as a dataframe
    df = pd.DataFrame({'cdr3a_pep_dist':[cdr3a_pep_dist],'cdr3b_pep_dist':[cdr3b_pep_dist],
                       'cdr2a_mhcb_dist':[cdr2a_mhcb_dist],'cdr2b_mhca_dist':[cdr2b_mhca_dist],
                       'cdr1a_mhcb_dist':[cdr1a_mhcb_dist],'cdr1b_mhca_dist':[cdr1b_mhca_dist],
                       'cdr3a_pep_angle':[cdr3a_pep_angle],'cdr3b_pep_angle':[cdr3b_pep_angle]})
    return df

#Write a function that will take the pdb file and calculate all the hydrogen bonds between two ranges
def calculate_hydrogen_bonds_between_residues(res1, res2):
    hydrogen_bonds = []
    atoms1 = [atom for atom in res1.get_atoms() if atom.element in ['N', 'O', 'H']]
    atoms2 = [atom for atom in res2.get_atoms() if atom.element in ['N', 'O', 'H']]
    for i, donor in enumerate(atoms1):
        for j, acceptor in enumerate(atoms2):
            if 'H' in donor.element or 'H' in acceptor.element:
                continue
            # Identify hydrogen atoms bonded to the donor
            distance = donor - acceptor
            if distance > 3.5:
                    continue
            bonded_hydrogens = [atom for atom in res1.get_atoms() if 'H' in atom.element and donor - atom < 2.2]
            # print(bonded_hydrogens)
            for hydrogen in bonded_hydrogens:
                # Calculate angle D-H...A
                vec_DH = hydrogen.coord - donor.coord
                vec_HA = acceptor.coord - hydrogen.coord
                cos_angle = np.dot(vec_DH, vec_HA) / (np.linalg.norm(vec_DH) * np.linalg.norm(vec_HA))
                angle = np.degrees(np.arccos(cos_angle))
                is_bond = (distance <= 3.5) and (angle >= 120)
                if is_bond:
                    don = donor.get_parent()
                    if don == res1:
                        donor_res = 1
                        acceptor_res = 2
                    else:
                        donor_res = 2
                        acceptor_res = 1
                    dict = {'Res1': res1.get_resname(), 'Res1_num': res1.get_id()[1],'Res1_pos': res1.get_full_id()[3][1],'Chain1': res1.get_parent().get_id(),
                            'Res2': res2.get_resname(), 'Res2_num': res2.get_id()[1],'Res2_pos': res2.get_full_id()[3][1],'Chain2': res2.get_parent().get_id(),
                            'Donor_res': donor_res, 'Acceptor_res': acceptor_res,
                            'Donor': donor.get_name(),'Acceptor': acceptor.get_name(),'Hydrogen': hydrogen.get_name(),
                            'Bond_length': distance, 'Bond_angle': angle,
                            'Donor_Charge': donor.get_bfactor(),'Donor_Radius': donor.get_occupancy(),
                            'Acceptor_Charge': acceptor.get_bfactor(),'Acceptor_Radius': acceptor.get_occupancy()}
                    hydrogen_bonds.append(dict)
    hbonds_df = pd.DataFrame(hydrogen_bonds)
    return hbonds_df


def calcRegionHbond(tcr, chain1, chain2, region1, region2):
    hbond_list = list()
    model = tcr[0]
    chain1_obj = model[chain1]
    chain2_obj = model[chain2]
    #Go through each residue in region 1 and region 2 and calculate hbonds pairwise 
    for res1 in chain1_obj:
        resid1 = res1.get_id()[1]
        if resid1 not in region1:
            continue
        for res2 in chain2_obj:
            resid2 = res2.get_id()[1]
            if resid2 not in region2:
                continue
            hbonds_df = calculate_hydrogen_bonds_between_residues(res1, res2) #This assumes res1 are all donors and res2 are all acceptors
            hbonds_df = calculate_hydrogen_bonds_between_residues(res2, res1)
            if not hbonds_df.empty:
            #     hbonds_df['Chain1'] = chain1
            #     hbonds_df['Res1_num'] = resid1
            #     hbonds_df['Chain2'] = chain2
            #     hbonds_df['Res2_num'] = resid2
                hbond_list.append(hbonds_df)
    if hbond_list:
        all_hbonds_df = pd.concat(hbond_list, ignore_index=True)
    else:
        all_hbonds_df = pd.DataFrame()
    return all_hbonds_df.shape[0] #We divide by 2 because we counted each hbond twice (once for each direction)

def getHbondStructure(tcr,ta,ma):
    chains = [x.id for x in tcr.get_chains()]
    chains = chains[:-1] #Remove any hetatm chains
    mhcachain = "A"
    if len(chains) > 4:
        pepchain = "C"
        mhcbchain = "B"
        tcrachain = "D"
        tcrbchain = "E"
    else:
        pepchain = "B"
        mhcbchain = "A"
        tcrachain = "C"
        tcrbchain = "D" 

    pep = range(1, len(fetch_entity(tcr[0], selection=pepchain))+ 1)
    #Define regions
    cdr1a = range(25+ta,31+ta)
    cdr2a = range(49+ta,54+ta)
    cdr3a = range(91+ta,97+ta)
    cdr1b = range(24+ta,31+ta)
    cdr2b = range(48+ta,56+ta)
    cdr3b = range(92+ta,101+ta)
    mhc_a = range(50+ma,85+ma)
    if len(chains) > 4:
        hasb = True
        mhc_b = range(50+ma,85+ma)
    else:
        mhc_b = range(139+ma,176+ma)
        hasb = False
    

    
    cdr3a_pep = calcRegionHbond(tcr,tcrachain,pepchain,cdr3a,pep)
    cdr3b_pep = calcRegionHbond(tcr,tcrbchain,pepchain,cdr3b,pep)
    # if "B" in chains:
    #     cdr2a_mhcb = calcRegionHbond(tcr,"D","B",cdr2a,mhc_b)
    #     cdr1a_mhcb = calcRegionHbond(tcr,"D","B",cdr1a,mhc_b)
    # else:
    cdr2a_mhcb = calcRegionHbond(tcr,tcrachain,mhcbchain,cdr2a,mhc_b)
    cdr1a_mhcb = calcRegionHbond(tcr,tcrachain,mhcbchain,cdr1a,mhc_b)
    cdr2b_mhca = calcRegionHbond(tcr,tcrbchain,mhcachain,cdr2b,mhc_a)
    cdr1b_mhca = calcRegionHbond(tcr,tcrbchain,mhcachain,cdr1b,mhc_a)

    pep_mhca = calcRegionHbond(tcr,pepchain,mhcachain,pep,mhc_a)
    # if "B" in chains:
    #     pep_mhcb = calcRegionHbond(tcr,pepchain,mhcbchain,pep,mhc_b)
    # else:
    pep_mhcb = calcRegionHbond(tcr,pepchain,mhcbchain,pep,mhc_b)
    cdr3a_cdr3b = calcRegionHbond(tcr,tcrachain,tcrbchain,cdr3a,cdr3b)

    #Return the number of hydrogen bonds between each region as a dataframe
    all_hbonds_df = pd.DataFrame({'cdr3a_pep_hbonds':[cdr3a_pep],'cdr3b_pep_hbonds':[cdr3b_pep],
                       'cdr2a_mhcb_hbonds':[cdr2a_mhcb],'cdr2b_mhca_hbonds':[cdr2b_mhca],
                       'cdr1a_mhcb_hbonds':[cdr1a_mhcb],'cdr1b_mhca_hbonds':[cdr1b_mhca],
                       'pep_mhca_hbonds':[pep_mhca],'pep_mhcb_hbonds':[pep_mhcb],
                       'cdr3a_cdr3b_hbonds':[cdr3a_cdr3b]})
    return all_hbonds_df

#################################################
######### Setup Variables ######################
#################################################
print("Setting up variables")

parser = argparse.ArgumentParser("simple_example")
parser.add_argument("--Dir")#,"Primary Directory in which the above structure will go/already exists", type=str)
parser.add_argument("--Out")
parser.add_argument("-ta","--Tadjust",help="Number to add to TCRs default 0",default=0, type=int)
parser.add_argument("-ma","--Madjust",help="Number to add to MHC default 0",default=0, type=int)
# parser.add_argument("-t","--Template",help="Turns on if there is a reference template to check against", action="store_true")
# parser.add_argument("-m","--Modeller",help="Set to T if PDB output from TCRpMHCModeller to fix chains", action="store_true")
# parser.add_argument("-s","--Sequence_File",help="File containing the sequences of the TCR and MHC chains in the order they are in the PDB file", type=str, default=None)

args = parser.parse_args()

ta = args.Tadjust
ma = args.Madjust
#Set by options
# alphafold = false 
# template = false
# Modeller= false
# if args.Modeller:
#     Modeller = true
# if args.Alphafold:
#     alphafold = true
# # if args.Template:
# #     template = true
# if args.Ground:
#     GroundDir = args.Ground
#     template = true

# if args.Sequence_File:
#     SequenceFile = args.Sequence_File
#     Sequences = pd.read_csv(SequenceFile, index_col=0)


## Set the directory names
ProjectDir = args.Dir 
outputDir = args.Out 

StructDir = os.path.join(ProjectDir,"Structures/")


RawDir = os.path.join(ProjectDir,"Structures/Raw/")
#This will need to be set from script options and if not present need to create a variable

NormDir = os.path.join(ProjectDir,"Structures/Norm/")
# Path(NormDir).mkdir(parents=True, exist_ok=True)

AlignedDir = os.path.join(ProjectDir,"Structures/Aligned/")
# Path(AlignedDir).mkdir(parents=True, exist_ok=True)

CleanDir = os.path.join(ProjectDir,"Structures/Cleaned/")
# Path(CleanDir).mkdir(parents=True, exist_ok=True)

CifDir = os.path.join(ProjectDir,"Structures/CIF")     
# Path(CifDir).mkdir(parents=True, exist_ok=True)

BrokeDir = os.path.join(ProjectDir,"Structures/Sub/")     
# Path(BrokeDir).mkdir(parents=True, exist_ok=True)

PhysicsDir = os.path.join(ProjectDir,"Physics/")
# Path(PhysicsDir).mkdir(parents=True, exist_ok=True)

ScoreDir = os.path.join(ProjectDir,"Scores/")
# Path(ScoreDir).mkdir(parents=True, exist_ok=True)

ContactsDir = os.path.join(ProjectDir,"Contacts/")
# Path(ContactsDir).mkdir(parents=True, exist_ok=True)

PropDir = os.path.join(ProjectDir,"Properties/")
# Path(PropDir).mkdir(parents=True, exist_ok=True)

lDDTDir = os.path.join(ScoreDir,"lDDT/")                #"D:/Lab/TCRBinding/ModelData/KnownStructures/ColabFold/Scores/lDDT/"
# Path(lDDTDir).mkdir(parents=True, exist_ok=True)

COMDir = os.path.join(ScoreDir,"COM/")                #"D:/Lab/TCRBinding/ModelData/KnownStructures/ColabFold/Scores/lDDT/"
# Path(COMDir).mkdir(parents=True, exist_ok=True)

DockDir = os.path.join(ScoreDir,"DockQ/")                #"D:/Lab/TCRBinding/ModelData/KnownStructures/ColabFold/Scores/lDDT/"
# Path(DockDir).mkdir(parents=True, exist_ok=True)

JsonDir = os.path.join(ScoreDir,"DockQ/json/")                #"D:/Lab/TCRBinding/ModelData/KnownStructures/ColabFold/Scores/lDDT/"
# Path(JsonDir).mkdir(parents=True, exist_ok=True)

devDir = os.path.join(ScoreDir,"Deviation/")             #"D:/Lab/TCRBinding/ModelData/KnownStructures/ColabFold/Scores/Deviation/"
# Path(devDir).mkdir(parents=True, exist_ok=True)

bDir = os.path.join(ScoreDir,"Bscore/")             #This holds the pLDDT if alphafolded, if not it is whatever the b score is for that file, sometimes atomic constraints
# Path(bDir).mkdir(parents=True, exist_ok=True)

ssaDir = os.path.join(PropDir,"SASA/")     
# Path(ssaDir).mkdir(parents=True, exist_ok=True)

AngleDir = os.path.join(PropDir,"Angles/")     
# Path(AngleDir).mkdir(parents=True, exist_ok=True)

Hdir = os.path.join(PhysicsDir,"Hbond/")
# Path(Hdir).mkdir(parents=True, exist_ok=True)

Electrodir = os.path.join(PhysicsDir,"Electrostatics/")
# Path(Electrodir).mkdir(parents=True, exist_ok=True)

Sdir = os.path.join(PhysicsDir,"Disulfide/")
# Path(Sdir).mkdir(parents=True, exist_ok=True)

Edir = os.path.join(PhysicsDir,"Energy/")
# Path(Edir).mkdir(parents=True, exist_ok=True)

GroDir = os.path.join(ProjectDir,"Structures/Gromacs/")
# try:
#     shutil.rmtree(GroDir)
# except:
#     pass
# Path(GroDir).mkdir(parents=True, exist_ok=True)

ChargedDir = os.path.join(GroDir,"charged/")
# Path(os.path.join(GroDir,"raw/")).mkdir(parents=True, exist_ok=True)
# Path(os.path.join(GroDir,"index/")).mkdir(parents=True, exist_ok=True)
# Path(os.path.join(GroDir,"itp/")).mkdir(parents=True, exist_ok=True)
# Path(os.path.join(GroDir,"processed/")).mkdir(parents=True, exist_ok=True)
# Path(os.path.join(GroDir,"solvated/")).mkdir(parents=True, exist_ok=True)
# Path(os.path.join(GroDir,"topology/")).mkdir(parents=True, exist_ok=True)
# Path(os.path.join(GroDir,"charged/")).mkdir(parents=True, exist_ok=True)
# Path(os.path.join(GroDir,"hydrogen/")).mkdir(parents=True, exist_ok=True)
# Path(os.path.join(GroDir,"simulated/")).mkdir(parents=True, exist_ok=True)

reffile1 = "./data/ref/1ao7.pdb"

##############################################
########    Main function ####################
###########################################

print("Starting feature extraction")
#### Calculate the solvent accessible surface area of the pmHC and then in complex with TCR
parser = Bio.PDB.PDBParser(QUIET=True)

# hbond_df = pd.DataFrame()
# files = os.listdir(ChargedDir)
# for x in tqdm(files):
#     name = x.split(".",1)[0]
#     file = os.path.join(ChargedDir,x)
#     #Now run the hbond calculation on each file
#     structure = parser.get_structure("structure", file)
#     chains = [x.id for x in structure.get_chains()]
#     hbonds = getHbondStructure(structure)
#     hbonds['Name'] = name
#     hbond_df = pd.concat([hbond_df,hbonds], ignore_index=True, axis=0)




files = os.listdir(CleanDir)
total_ssa = []
pMHC_ssa = []
tcr_ssa = []
bsa = []
names = []
combined = pd.DataFrame()
for x in tqdm(files):
    name = x.split(".",1)[0]
    file = os.path.join(CleanDir,x)
    # we need to read in the structure twice and then chop one file down to just the pMHC and then calc SASA of each
    structure = parser.get_structure("structure", file)
    chains = [x.id for x in structure.get_chains()]
    # print(name,len(chains))
    if not "C" in chains:
        print(f"Skipping {name} as no peptide chain found")
        continue
    tcr_chains = chains[-2:]
    sub_chains = chains[:-2]
    all = parser.get_structure("structure", file)
    tcr = subset_chains(structure,tcr_chains)
    pMHC = subset_chains(structure,sub_chains)
    sr = SASA.ShrakeRupley() #ShrakeRupley()
    sr.compute(all, level="S")
    sr.compute(pMHC, level="S")
    sr.compute(tcr,level="S")

    # total_ssa.append(all.sasa)
    # pMHC_ssa.append(pMHC.sasa)
    # tcr_ssa.append(tcr.sasa) 
    # bsa.append((tcr.sasa + pMHC.sasa)-all.sasa)
    # names.append(name)
    ssa = pd.DataFrame({'Total_SSA':[all.sasa],'pMHC_SSA':[pMHC.sasa],'TCR_SSA':[tcr.sasa],
                        'BSA':[(tcr.sasa + pMHC.sasa)-all.sasa]})

    ####### Now calculate the COM for regions
    comdf = calcCOMAngles(structure,ta,ma)

    ###### Calculate hydrogen bonds
    chargefile = os.path.join(ChargedDir,x)
    if not os.path.exists(chargefile):
        print(f"Charged file for {name} not found, skipping hbonds")
        hbonddf = pd.DataFrame({'cdr3a_pep_hbonds':[nan],'cdr3b_pep_hbonds':[nan],
                       'cdr2a_mhcb_hbonds':[nan],'cdr2b_mhca_hbonds':[nan],
                       'cdr1a_mhcb_hbonds':[nan],'cdr1b_mhca_hbonds':[nan],
                       'pep_mhca_hbonds':[nan],'pep_mhcb_hbonds':[nan],
                       'cdr3a_cdr3b_hbonds':[nan]})
        continue
    else:
        structure = parser.get_structure("structure", chargefile)
        hbonddf = getHbondStructure(structure,ta,ma)

    #now combine the above columns into one row of a dataframe and update iteratively
    tmp = pd.concat([ssa,comdf,hbonddf],axis=1)
    tmp['Name'] = name
    combined = pd.concat([combined,tmp], ignore_index=True, axis=0)
    

# print(combined.head())
#Finally reread in the energy dataframe which should be in physicsDir under "Binding_Energies.csv" and merge with combined using the "Name" column on both df
if (not os.path.exists(os.path.join(PhysicsDir,"Binding_Energies.csv"))):
    print(f"Energy file not found at {os.path.join(PhysicsDir,'Binding_Energies.csv')}, using prodigy scores instead")
    #Read in and combine into one dataframe all of the individual files in Prodigy directory
    prodigyfiles = os.listdir(os.path.join(PhysicsDir,"Energy/"))
    energydf = pd.DataFrame()
    for x in prodigyfiles:
        tempdf = pd.read_csv(os.path.join(PhysicsDir,"Energy/",x),sep='\t',header=None)
        energydf = pd.concat([energydf,tempdf], ignore_index=True, axis=0)
    energydf.columns = ['Name', 'Prodigy_Energy']
else:
    energydf = pd.read_csv(os.path.join(PhysicsDir,"Binding_Energies.csv"))

finaldf = pd.merge(combined, energydf, on='Name')

print(finaldf.head())

finaldf.to_csv(os.path.join(outputDir,"features.csv"))

print(f"Saved features to {os.path.join(outputDir,'features.csv')}")

print("         finished")

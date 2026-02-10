#!/usr/bin/env python3

#Need to activate the conda environment in wsl before running
#conda activate wsl_tcr
#Run: python -Wi ./scripts/Calc_Energy.py ./ModelData/KnownStructures/GroundTruth/Structures/Norm/ ./ModelData/KnownStructures/GroundTruth/Structures/Sub/ ./ModelData/KnownStructures/GroundTruth/Physics/

#This script will calculate the deltaG for binding of the TCR-pMHC construct
#Energies are calculated from the scoring function in Rosetta, however, a more accurate score can be calculated through simulation https://ambermd.org/tutorials/advanced/tutorial3/section1.php
#We will actually calculate three different binding energies
#   - Gcomplex =  The free energy of the entire complex
#   - Gtcr = The energy score for the TCR
#   - Gpep = The energy score for just the peptide
#   - Gmhc = The energy score for just the mhc
#   - Gpep/mhc = The free energy of the peptide mhc complex
#   - Gtcr/mhc = The free energy of the tcr/mhc interaction minus the peptide
#   - Gbind = Gcomplex - (Gtcr/mhc + Gpep/mhc)
#   - Gbind2 = Gcomplex - (Gtcr + Gpmhc)

from math import nan
import pyrosetta
import pyrosetta.rosetta.core.pose as p
from pyrosetta import *
from rosetta.core.kinematics import MoveMap
from rosetta.core.kinematics import FoldTree
from rosetta.core.pack.task import TaskFactory
from rosetta.core.pack.task import operation
from rosetta.core.simple_metrics import metrics
from rosetta.core.select import residue_selector as selections
from rosetta.core import select
from rosetta.core.select.movemap import *
#Protocol Includes
from rosetta.protocols import minimization_packing as pack_min
from rosetta.protocols import relax as rel
from rosetta.protocols.antibody.residue_selector import CDRResidueSelector
from rosetta.protocols.antibody import *
from rosetta.protocols.loops import *
from rosetta.protocols.relax import FastRelax
init('-use_input_sc -input_ab_scheme AHo_Scheme -ignore_unrecognized_res \
     -ignore_zero_occupancy false -load_PDB_components false -relax:default_repeats 2 -no_fconfig')

from Bio.PDB import PDBParser
from pyrosetta.teaching import *

import sys
import os
import time
from tqdm import tqdm
import copy
import pickle
import pandas as pd

args = sys.argv
dir = args[1]
outdir = args[3]
chaindir = args[2]

# scorefxn = get_fa_scorefxn()
scorefxn = get_score_function(True)
data = []
############### Pyrosetta
# The task factory accepts all the task operations
# tf = TaskFactory()
# tf.push_back(operation.InitializeFromCommandline())
# tf.push_back(operation.RestrictToRepacking())
# packer = pack_min.PackRotamersMover()
# packer.task_factory(tf)
# #apply the mover
# packer.apply(pose)
fr = FastRelax()
fr.set_scorefxn(scorefxn)
fr.max_iter(100)

#Write a function to subset chains
from Bio.PDB import Select
from Bio.PDB import PDBIO

class SelectChains(Select):
    def __init__(self, chain_ids):
        self.chain_ids = chain_ids

    def accept_chain(self, chain):
        return chain.get_id() in self.chain_ids
#Write a function that will take in the pdb components and run the relaxation/scoring on them

def score_pdb(pdbfile):
    try:
        pose = pose_from_pdb(pdbfile)
        original_pose = pose.clone()
        fr.apply(pose)
        before = scorefxn(original_pose)
        after = scorefxn(pose)
        return before, after
    except:
        return nan, nan


for f in tqdm(os.listdir(dir)):
    name = f.split(".")[0]
    print("##############################################################################")
    print(f"Starting {name}")
    print("##############################################################################")
    #First we need to break the structure into it's components and save each one
    start = time.time()
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", os.path.join(dir,f))
    model = structure[0]
    chains = []
    for chain in model:
        chains.append(chain.get_id())
    if "C" in chains:
        has_pep = True
    else:
        has_pep = False
    if "B" in chains:
        has_mhcb = True
    else:
        has_mhcb = False
    #Save chains to individual files in output directory
    pmhc = ['A','B','C']
    tcrmhc = ['A','B','D','E']
    tcr = ['D','E']
    peptide = ['C']
    mhc = ['A','B']

    io = PDBIO()
    io.set_structure(structure)
    #Save the different components
    if not os.path.exists(chaindir):
        raise Exception(f"Directory {chaindir} does not exist, please create it")
    io.save(os.path.join(chaindir,name + "".join(pmhc) +".pdb"), SelectChains(pmhc))
    io.save(os.path.join(chaindir,name + "".join(tcrmhc) +".pdb"), SelectChains(tcrmhc))
    io.save(os.path.join(chaindir,name + "".join(tcr) +".pdb"), SelectChains(tcr))
    if has_pep:
        io.save(os.path.join(chaindir,name + "".join(peptide) +".pdb"), SelectChains(peptide))    
    io.save(os.path.join(chaindir,name + "".join(mhc) +".pdb"), SelectChains(mhc))

    #Now that the individual chains are saved, we can load them into pyrosetta and score them
    print("##########################  Starting Relax and Calc #####################################")
    allbefore, allafter = score_pdb(os.path.join(dir,f))
    print(f"Full complex scores: {allbefore} -> {allafter}")
    pmhcbefore, pmhcafter = score_pdb(os.path.join(chaindir,name + "".join(pmhc) +".pdb"))
    print(f"pmhc scores: {pmhcbefore} -> {pmhcafter}")
    tcrmhcbefore, tcrmhcafter = score_pdb(os.path.join(chaindir,name + "".join(tcrmhc) +".pdb"))
    print(f"tcrmhc scores: {tcrmhcbefore} -> {tcrmhcafter}")
    tcrbefore, tcrafter = score_pdb(os.path.join(chaindir,name + "".join(tcr) +".pdb"))
    print(f"tcr scores: {tcrbefore} -> {tcrafter}")
    if has_pep:
        peptidebefore, peptideafter = score_pdb(os.path.join(chaindir,name + "".join(peptide) +".pdb"))
        print(f"peptide scores: {peptidebefore} -> {peptideafter}")
    else:
        peptidebefore, peptideafter = nan, nan
        print(f"No peptide found")
    mhcbefore, mhcafter = score_pdb(os.path.join(chaindir,name + "".join(mhc) +".pdb"))
    print(f"mhc scores: {mhcbefore} -> {mhcafter}")


    ############################################################
    ######### Scoring#########################################
    ############################################################

    print("########################## Calculating scores #####################################")

    Gcom = allafter #Get the overall score of the complex
    
    Gtcr = tcrafter
    
    Gmhc = mhcafter
    Gtcrmhc = tcrmhcafter
    Gpepmhc = pmhcafter
    if has_pep:
        Gpep = peptideafter
        Gbind = Gcom - (Gtcrmhc + Gpep)
        Gbind2 = Gcom - (Gtcr + Gpepmhc)
    else:
        Gpep = nan
        Gbind = Gcom - (Gtcrmhc + 0) #If there is no peptide then Gbind1 = Gcom - Gtcrmhc
        Gbind2 = Gcom - (Gtcr + Gmhc)
    
    dat = {"Name": name, "MHCB": has_mhcb, "Peptide": has_pep,
           "Gcom": allbefore,"Gcom_relaxed": allafter, 
           "Gtcr": tcrbefore, "Gtcr_relaxed": tcrafter,
           "Gpep": peptidebefore,"Gpep_relaxed": peptideafter, 
           "Gmhc": mhcbefore,"Gmhc_relaxed": mhcafter,
           "Gpmhc": pmhcbefore,"Gpmhc_relaxed": pmhcafter, 
           "Gtcrmhc": tcrmhcbefore, "Gtcrmhc_relaxed": tcrmhcafter, 
           "Gbind1": Gbind, "Gbind2": Gbind2}
    data.append(dat)
    print(dat)
    
    end = time.time()
    print(f"Took {end - start} seconds")
    # except:
    #     print("Could not calculate scores")
    #     continue
        

    # print(f"Took {end - start} seconds")

df = pd.DataFrame(data)
print(df)
df.to_csv(os.path.join(outdir,"Binding_Energies.csv"),index = False)

print(f"Wrote dataframe to {outdir}")





#############Depraecated code##############
    # try:
    #     pose = pose_from_pdb(os.path.join(dir,f))
    #     original_pose = pose.clone()
    #     fr.apply(pose)
    #     before = scorefxn(original_pose)
    #     print("relaxed ",scorefxn(pose))
    #     print("original",before)
    # except:
    #     print("Could not read in model")
    #     continue
    # mhc = p.deep_copy(pose)
    # mhc_pep = p.deep_copy(pose)
    # pep = p.deep_copy(pose)
    # tcr = p.deep_copy(pose)
    # mhc_tcr = p.deep_copy(pose)
    # all = p.deep_copy(pose)

    # mhca_res = []
    # mhcb_res = []
    # pep_res = []
    # tcra_res = []
    # tcrb_res = []



    


    #Now we need to split the pose into just the tcr, the peptide, the mhc, the peptide mhc and the tcr mhc to calculate the above
    #Set up details
    
    # has_B = True
    # try:
    #     chainB_id = pyrosetta.rosetta.core.pose.get_chain_id_from_chain("B",pose)
    # except: 
    #     has_B = False
    
    # chain_length = []
    # has_pep = False
    # for n in range(pose.num_chains()):
    #     chain_length.append(len(p.get_resnums_for_chain_id(pose,n+1)))
    #     if len(p.get_resnums_for_chain_id(pose,n)) < 25 and len(p.get_resnums_for_chain_id(pose,n)) > 0:
    #         has_pep = True

    # if pose.num_chains() >= 5:
    #     has_mhcb = True
    # elif pose.num_chains() == 4:
    #     if has_pep:
    #         has_mhcb = False
    # else:
    #     has_mhcb = False
    #     has_pep = False
    
    # # print(chain_length)
    # # print(f"MHCB {has_mhcb}, Pep {has_pep}")
    # # print("############################### Making objects ####################################")

    # #Make objects
    # try:
    #     mhca_res = p.get_resnums_for_chain_id(pose,1)

    #     if has_pep: 
    #         if has_mhcb: #Has all 5 chains
    #             mhcb_res = p.get_resnums_for_chain_id(pose,2)
    #             pep_res = p.get_resnums_for_chain_id(pose,3)
    #             tcra_res = p.get_resnums_for_chain_id(pose,4)
    #             tcrb_res = p.get_resnums_for_chain_id(pose,5)

    #             mhc_array = p.get_resnums_for_chain_id(pose,1)
    #             mhc_array.extend(p.get_resnums_for_chain_id(pose,2))
    #             mhc_array2 = p.get_resnums_for_chain_id(pose,1)
    #             mhc_array2.extend(p.get_resnums_for_chain_id(pose,2))

    #             mhc_pep_array = p.get_resnums_for_chain_id(pose,1)
    #             mhc_pep_array.extend(p.get_resnums_for_chain_id(pose,2))
    #             mhc_pep_array.extend(p.get_resnums_for_chain_id(pose,3))

    #             tcr_array = p.get_resnums_for_chain_id(pose,4)
    #             tcr_array.extend(p.get_resnums_for_chain_id(pose,5))
    #             tcr_array2 = p.get_resnums_for_chain_id(pose,4)
    #             tcr_array2.extend(p.get_resnums_for_chain_id(pose,5))

    #             mhc_tcr_array = mhc_array2
    #             mhc_tcr_array.extend(tcr_array2)
    #         else: #Has 4 chains with one being a peptide
    #             pep_res = p.get_resnums_for_chain_id(pose,2)
    #             tcra_res = p.get_resnums_for_chain_id(pose,3)
    #             tcrb_res = p.get_resnums_for_chain_id(pose,4)

    #             mhc_array = p.get_resnums_for_chain_id(pose,1)
    #             mhc_array2 = p.get_resnums_for_chain_id(pose,1)

    #             mhc_pep_array = p.get_resnums_for_chain_id(pose,1)
    #             mhc_pep_array.extend(p.get_resnums_for_chain_id(pose,2))

    #             tcr_array = p.get_resnums_for_chain_id(pose,3)
    #             tcr_array.extend(p.get_resnums_for_chain_id(pose,4))

    #             mhc_tcr_array = p.get_resnums_for_chain_id(pose,1)
    #             mhc_tcr_array.extend(p.get_resnums_for_chain_id(pose,3))
    #             mhc_tcr_array.extend(p.get_resnums_for_chain_id(pose,4))
    #     else:
    #         if has_mhcb: #Has 4 chains with no pep
    #             mhcb_res = p.get_resnums_for_chain_id(pose,2)
    #             # pep_res = p.get_resnums_for_chain_id(pose,3)
    #             pep_res = []
    #             tcra_res = p.get_resnums_for_chain_id(pose,3)
    #             tcrb_res = p.get_resnums_for_chain_id(pose,5)

    #             mhc_array = p.get_resnums_for_chain_id(pose,1)
    #             mhc_array.extend(p.get_resnums_for_chain_id(pose,2))
    #             mhc_array2 = p.get_resnums_for_chain_id(pose,1)
    #             mhc_array2.extend(p.get_resnums_for_chain_id(pose,2))

    #             mhc_pep_array = p.get_resnums_for_chain_id(pose,1)
    #             mhc_pep_array.extend(p.get_resnums_for_chain_id(pose,2))
    #             # mhc_pep_array.extend(p.get_resnums_for_chain_id(pose,3))

    #             tcr_array = p.get_resnums_for_chain_id(pose,4)
    #             tcr_array.extend(p.get_resnums_for_chain_id(pose,5))
    #             tcr_array2 = p.get_resnums_for_chain_id(pose,4)
    #             tcr_array2.extend(p.get_resnums_for_chain_id(pose,5))

    #             mhc_tcr_array = mhc_array2
    #             mhc_tcr_array.extend(tcr_array2)
    #         else:
    #             if pose.num_chains() == 3:
    #                 mhc_array = p.get_resnums_for_chain_id(pose,1)
    #                 tcr_array = p.get_resnums_for_chain_id(pose,2)
    #                 tcr_array.extend(p.get_resnums_for_chain_id(pose,3))
    #                 mhc_tcr_array = p.get_resnums_for_chain_id(pose,1)
    #                 mhc_tcr_array.extend(p.get_resnums_for_chain_id(pose,2))
    #                 mhc_tcr_array.extend(p.get_resnums_for_chain_id(pose,3))
    #             else:
    #                 pep_res = [] #p.get_resnums_for_chain_id(pose,2)
    #                 tcra_res = p.get_resnums_for_chain_id(pose,2)
    #                 tcrb_res = p.get_resnums_for_chain_id(pose,3)

    #                 mhc_array = p.get_resnums_for_chain_id(pose,1)
    #                 mhc_array2 = p.get_resnums_for_chain_id(pose,1)

    #                 mhc_pep_array = p.get_resnums_for_chain_id(pose,1)
    #                 # mhc_pep_array.extend(p.get_resnums_for_chain_id(pose,2))

    #                 tcr_array = p.get_resnums_for_chain_id(pose,2)
    #                 tcr_array.extend(p.get_resnums_for_chain_id(pose,3))

    #                 mhc_tcr_array = p.get_resnums_for_chain_id(pose,1)
    #                 mhc_tcr_array.extend(p.get_resnums_for_chain_id(pose,2))
    #                 mhc_tcr_array.extend(p.get_resnums_for_chain_id(pose,3))
    # except:
    #     # print("Could not get model info")
    #     continue
    
    # print("########################## Building the poses #####################################")

    # try:
    #     # print(1)
    #     p.pdbslice(mhc,mhc_array)
    #     mhc.conformation().delete_residue_range_slow(1, len(mhc_array)-1);
    #     mhc.pdb_info().obsolete(False);
    #     # print(2)
    #     if has_pep:
    #         p.pdbslice(pep,pep_res)
    #         pep.conformation().delete_residue_range_slow(1, len(pep_res)-1);
    #         pep.pdb_info().obsolete(False);
        
    #         p.pdbslice(mhc_pep,mhc_pep_array)
    #         mhc_pep.conformation().delete_residue_range_slow(1, len(mhc_pep_array)-1);
    #         mhc_pep.pdb_info().obsolete(False);
    #     # print(3)
    #     p.pdbslice(tcr,tcr_array)
    #     tcr.conformation().delete_residue_range_slow(1, len(tcr_array)-1);
    #     tcr.pdb_info().obsolete(False);
    #     # print(4)
    #     p.pdbslice(mhc_tcr,mhc_tcr_array)
    #     mhc_tcr.conformation().delete_residue_range_slow(1, len(mhc_tcr_array)-1);
    #     mhc_tcr.pdb_info().obsolete(False)
    #     # print(5)
    # except:
    #     print("Could not make models")
    #     continue
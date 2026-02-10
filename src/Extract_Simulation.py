#!/usr/bin/env python3

#In this document we will use MDAnalysis to analyze and pull data from the molecular dynamics simulations

#run Syntax
# conda activate wsl_tcr 
# python3 ./scripts/Extract_Simulation.py ./Simulations/GroundTruth/


from pathlib import Path
import os
import sys
import time
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
# import nglview as nv
import MDAnalysis as mda
from MDAnalysis.analysis import rms, diffusionmap, align, gnm
from MDAnalysis.analysis.distances import dist, distance_array
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
from Bio.PDB import PDBList, PDBParser, Structure, Model, Chain, Residue, Atom

##########Pyrosetta
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

# scorefxn = get_score_function(True)
# fr = FastRelax() 
# fr.set_scorefxn(scorefxn)
# fr.max_iter(100)
# We probably shouldn't do relaxation because we already equilibrate right? We are more interested in the change in energy

#First set up the paths to load the data
mainDir = sys.argv[1]
dirs = [f.path for f in os.scandir(mainDir) if f.is_dir()]

resFile = mainDir + "regions.csv" #Read in a file that contains the regions of interest, we will these to set selection groups for the analysis functions
if os.path.isfile(resFile):
    regions_df = pd.read_csv(resFile)
    print(regions_df.head())
else:
    print("No regions file found, proceeding without it.")

############################
##### Functions ############
############################
#Write functions to extract data from the simulations

def atomgroup_from_resnum(universe, resnum_list):
    """Create an atom group from a list of residue numbers.

    Parameters
    ----------
    universe: MDAnalysis.core.universe.Universe
        MDAnalysis universe.
    resnum_list: list of int
        List of residue numbers.

    Returns
    -------
    atomgroup: MDAnalysis.core.groups.AtomGroup
        Atom group containing the selected residues.
    """
    selection_str = "resnum " + " ".join(map(str, resnum_list))
    atomgroup = universe.select_atoms(selection_str)
    return atomgroup

def rmsd_for_atomgroups(universe, selection1, selection2=None, cols = None):
    """Calulate the RMSD for selected atom groups.

    Parameters
    ----------
    universe: MDAnalysis.core.universe.Universe
        MDAnalysis universe.
    selection1: str
        Selection string for main atom group, also used during alignment.
    selection2: list of str, optional
        Selection strings for additional atom groups.

    Returns
    -------
    rmsd_df: pandas.core.frame.DataFrame
        DataFrame containing RMSD of the selected atom groups over time.
    """

    universe.trajectory[0]
    ref = universe
    rmsd_analysis = rms.RMSD(universe, ref, select=selection1, groupselections=selection2)
    rmsd_analysis.run()
    columns = [selection1, *cols] if selection2 else [selection1]
    rmsd_df = pd.DataFrame(np.round(rmsd_analysis.results.rmsd[:, 2:], 2), columns=columns)
    rmsd_df.index.name = "frame"
    return rmsd_df

# This function will extract frames from the trajectory and convert them into Biopython Structure objects for further analysis or visualization.
def extract_frames(universe, selection,step = 100, cols=None):
    """Extract frames for selected atom groups.

    Parameters
    ----------
    universe: MDAnalysis.core.universe.Universe
        MDAnalysis universe.
    selection: str
    """
    sele = u.select_atoms(selection)
    structure_list = []
    times = []
    frames = []
    i = -1
    for ts in universe.trajectory:
        i += 1
        if i % step != 0: # Only do every 100th frame which is about a nanosecond in time (10 ps per frame * 100 frames = 1 ns)
            continue
        # Create a new Biopython Structure object for the current frame
        times.append(ts.time)
        frames.append(ts.frame)
        bio_structure = Structure.Structure("frame_" + str(ts.frame))
        bio_model = Model.Model(0) # Biopython Structures typically have a single model
        
        # Iterate through residues and atoms in the MDAnalysis Universe
        # and add them to the Biopython Structure
        bio_chain = Chain.Chain("A") # Or use chain.chainid if available
        r = 0    
        for residue in sele.residues:
            bio_residue = Residue.Residue(
                r, 
                residue.resname,
                "A"
            )
            r += 1
            
            for atom in residue.atoms:
                bio_atom = Atom.Atom(
                    atom.name, 
                    atom.position, 
                    0,# atom.bfactor, 
                    1,# atom.occupancy, 
                    " ", #altloc
                    (" " + atom.name + " "), # fullname
                    atom.id,# atom.number, 
                    atom.type # element
                )
                bio_residue.add(bio_atom)
            bio_chain.add(bio_residue)
        bio_model.add(bio_chain)
        bio_structure.add(bio_model)
        structure_list.append(bio_structure)
    return structure_list, times, frames

#Take the output of the previous and do the pyrosetta calculation 
def calculate_energy(universe, selections, step = 100, cols=None):
    """Calculate energy for selected atom groups using PyRosetta.

    Parameters
    ----------
    universe: MDAnalysis.core.universe.Universe
        MDAnalysis universe.
    selection: str
        Selection string for main atom group.
    step: int, optional
        Step size for frame extraction.

    Returns
    -------
    energy_df: pandas.core.frame.DataFrame
        DataFrame containing energy of the selected atom groups over time.
    """
    structure_list, times, frames = extract_frames(universe, selections, step=step, cols=cols)
    
    energies = []
    for pose in structure_list:
        energy = scorefxn(pose)
        energies.append(energy)
    
    energy_df = pd.DataFrame({
        'frame': frames,
        'time': times,
        'energy': energies
    })
    
    return energy_df

def calc_dists(universe, selection1, selection2, step=100):
    #Calculate the distance between two selections over time. This will return a DataFrame with the distances for each frame.
    distances = []
    # times = []
    # frames = []
    for ts in universe.trajectory:
        if ts.frame % step != 0:
            continue
        # times.append(ts.time)
        # frames.append(ts.frame)
        sele1 = universe.select_atoms(selection1)
        sele2 = universe.select_atoms(selection2)
        dist_matrix = distance_array(sele1.positions, sele2.positions)
        # min_dist = np.min(dist_matrix)
        # Instead of just taking the minimum distance, we can calculate the average distance between the two selections for a more comprehensive analysis.
        avg_dist = np.mean(dist_matrix)
        distances.append(avg_dist)
    return distances

def calc_Hbonds(universe, selection1, selection2,name, step=10): #We use steps of 10 which are 100 ps apart
    #Calculate hydrogen bonds between two selections over time.
    hbond = HBA(universe, hydrogens_sel="name H*", acceptors_sel="name O* N*",d_a_cutoff=3, d_h_a_angle_cutoff=150,
            between=[[selection1,selection2],
                     [selection1,"resname TIP3 and around 10 not resname TIP3"],
                     [selection2,"resname TIP3 and around 10 not resname TIP3"]])
    hbond.run(step = step) 
    times = hbond.times
    frames = hbond.frames
    counts = hbond.count_by_time()
    hbond_results = pd.DataFrame({
        'Name': name,
        'frame': frames,
        'time': times,
        'hbond_count': counts
    })
    #Get individual hbonds as results
    df = pd.DataFrame(hbond.results.hbonds[:, :DISTANCE].astype(int),
                  columns=["Frame",
                           "Donor_ix",
                           "Hydrogen_ix",
                           "Acceptor_ix",])

    df["Distances"] = hbond.results.hbonds[:, DISTANCE]
    df["Angles"] = hbond.results.hbonds[:, ANGLE]

    df["Donor resname"] = u.atoms[df.Donor_ix].resnames
    df["Acceptor resname"] = u.atoms[df.Acceptor_ix].resnames
    df["Donor resid"] = u.atoms[df.Donor_ix].resids
    df["Acceptor resid"] = u.atoms[df.Acceptor_ix].resids
    df["Donor name"] = u.atoms[df.Donor_ix].names
    df["Acceptor name"] = u.atoms[df.Acceptor_ix].names
    return hbond, hbond_results, df

# def contacts_within_cutoff(u, group_a, group_b, radius=4.5):
#     timeseries = []
#     for i, n in zip(group_a, group_b):

#     for ts in u.trajectory:
#         # calculate distances between group_a and group_b
#         dist = contacts.distance_array(group_a.positions, group_b.positions)
#         # determine which distances <= radius
#         n_contacts = contacts.contact_matrix(dist, radius).sum()
#         timeseries.append([ts.frame, n_contacts])
#     return np.array(timeseries)

###########################
##### Constants ###########
###########################

cols = ['MHC', 'Peptide', 'TCRa', 'TCRb',
        'CDR3a', 'CDR3b', 'CDR2a', 'CDR2b', 'CDR1a', 'CDR1b']

FRAME = 0
DONOR = 1
HYDROGEN = 2
ACCEPTOR = 3
DISTANCE = 4
ANGLE = 5


###########################
##### Main Program ########
###########################

#Loop over the directories and extract the data
for d in dirs:
    elapsed = time.time()
    print("Processing Directory: ", d)
    #Check if md.xtc exists in the directory
    if not Path(d + "/md.xtc").is_file():
        print("No md.xtc file found in ", d)
        continue

    #Create analysis directory
    analysis_dir = Path(d + "/analysis")
    #If the analysis directory already exists skip this directory else create it
    if analysis_dir.is_dir():
        print("Analysis directory already exists, skipping this directory.")
        continue
    else:
        analysis_dir.mkdir(exist_ok=True)

    ###################################################################
    ##################Set the selections for the analysis##############
    ##################################################################
    sub = regions_df[regions_df['ID'] == Path(d).name]
    mhc_len = sub['MHCa_length'].values[0] + sub['MHCb_length'].values[0]
    tcra_len = sub['TCRa_length'].values[0]
    tcrb_len = sub['TCRb_length'].values[0]
    pep_len = sub['Peptide_length'].values[0]
    pep_start = mhc_len
    tcra_start = mhc_len + pep_len
    tcrb_start = mhc_len + pep_len + tcra_len
    
    mhc_sele = "resnum " + "0" + "-" + str(mhc_len - 1)
    mhca_sele = "resnum 50 - 85"  #Example MHCa selection for class 1
    mhcb_sele = "resnum 138 - 180"   #Example MHCb selection for class 1
    pep_sele = "resnum " + str(mhc_len) + "-" + str(mhc_len + pep_len - 1)
    tcra_sele = "resnum " + str(mhc_len + pep_len) + "-" + str(mhc_len + pep_len + tcra_len - 1)
    tcrb_sele = "resnum " + str(mhc_len + pep_len + tcra_len) + "-" + str(mhc_len + pep_len + tcra_len + tcrb_len - 1)

    cdr3a_sele = "resnum " + str(tcra_start + sub['CDR3a_start'].values[0]) + "-" + str(tcra_start + sub['CDR3a_end'].values[0])
    cdr3b_sele = "resnum " + str(tcrb_start + sub['CDR3b_start'].values[0]) + "-" + str(tcrb_start + sub['CDR3b_end'].values[0])
    cdr2a_sele = "resnum " + str(tcra_start + sub['CDR2a_start'].values[0]) + "-" + str(tcra_start + sub['CDR2a_end'].values[0])
    cdr2b_sele = "resnum " + str(tcrb_start + sub['CDR2b_start'].values[0]) + "-" + str(tcrb_start + sub['CDR2b_end'].values[0])
    cdr1a_sele = "resnum " + str(tcra_start + sub['CDR1a_start'].values[0]) + "-" + str(tcra_start + sub['CDR1a_end'].values[0])
    cdr1b_sele = "resnum " + str(tcrb_start + sub['CDR1b_start'].values[0]) + "-" + str(tcrb_start + sub['CDR1b_end'].values[0])
    
    ########################################
    #######Read in Trajectory###############
    ########################################
    gmxfile = d + "/md.xtc"
    topfile = d + "/npt.gro"
    if not Path(gmxfile).is_file() or not Path(topfile).is_file():
        print("No gro or xtc file found in ", d)
        continue
    #Load the MD simulation
    u = mda.Universe(d + "/npt.gro", d + "/md.xtc",guess_bonds=True)
    #Print the number of trajectory frames and the number of atoms in the system
    print(f"Number of frames: {len(u.trajectory)}")
    print(f"Number of atoms: {len(u.atoms)}")
    #Print the time step of the trajectory
    print(f"Time step: {u.trajectory.dt} ps")
    #Print the total simulation time
    total_time = len(u.trajectory) * u.trajectory.dt
    print(f"Total simulation time: {total_time/1000} ns")
    # u = mda.Universe(d + "/tcr.pdb", d + "/md.xtc",guess_bonds=True)
    # algin the trajectory to the first frame
    u.trajectory[0]
    alignment = align.AlignTraj(
        mobile=u, reference=u, select="protein", in_memory=True
    )
    alignment.run()

    ######################################################
    ###########      Extract RMSD Data      ##############
    ######################################################
    print("Calculating RMSD Data...")
    rmsd_df = rmsd_for_atomgroups(u, "protein", [mhc_sele, pep_sele, tcra_sele, tcrb_sele,
                                                 cdr3a_sele, cdr3b_sele, cdr2a_sele, cdr2b_sele, cdr1a_sele, cdr1b_sele],cols=cols)
    print(rmsd_df.head())
    
    rmsd_df.to_csv(analysis_dir / "rmsd.csv")

    ######################################################
    ###########   Extract Distance Data   ################
    ######################################################
    print("Calculating Distance Data...")
    #Calculate per frame distances between:
    # mhc and peptide
    # pep and CDR3
    # MHCa and CDR1b and CDR2b
    # MHCb and CDR1a and CDR2a
    cdr3a_pep_distances = calc_dists(u, pep_sele, cdr3a_sele)
    cdr3b_pep_distances = calc_dists(u, pep_sele, cdr3b_sele)
    mhca_cdr1b_distances = calc_dists(u, mhca_sele, cdr1b_sele)
    mhca_cdr2b_distances = calc_dists(u, mhca_sele, cdr2b_sele)
    mhcb_cdr1a_distances = calc_dists(u, mhcb_sele, cdr1a_sele)
    mhcb_cdr2a_distances = calc_dists(u, mhcb_sele, cdr2a_sele)
    distance_df = pd.DataFrame({
        'frame': range(0, len(cdr3a_pep_distances)*100, 100),
        'cdr3a_pep_dist': cdr3a_pep_distances,
        'cdr3b_pep_dist': cdr3b_pep_distances,
        'mhca_cdr1b_dist': mhca_cdr1b_distances,
        'mhca_cdr2b_dist': mhca_cdr2b_distances,
        'mhcb_cdr1a_dist': mhcb_cdr1a_distances,
        'mhcb_cdr2a_dist': mhcb_cdr2a_distances
    })
    distance_df.to_csv(analysis_dir / "distances.csv", index=False)

    ######################################################
    ###########   Calculate elastic network    ###########
    ######################################################

    print("Calculating Elastic Network Model...")
    # Calculate the elastic network model for the protein

    ######################################################
    ###########   Extract Hydrogen Bond Data   ###########
    ######################################################
    print("Calculating Hydrogen Bond Data...")
    #Extract all hydrogen bonds over all frames
    # hbond, mhca_pep_results, hbond_df = calc_Hbonds(u, mhca_sele, pep_sele,"mhca_pep", step=10)
    # hbond, mhcb_pep_results, hbond_df2 = calc_Hbonds(u, mhcb_sele, pep_sele,"mhcb_pep", step=10)
    hbond, tcra_pep_results, hbond_df3 = calc_Hbonds(u, cdr3a_sele, pep_sele,"cdr3a_pep", step=10)
    hbond, tcrb_pep_results, hbond_df4 = calc_Hbonds(u, cdr3b_sele, pep_sele,"cdr3b_pep", step=10)
    # hbond, cdr1b_mhca_results, hbond_df5 = calc_Hbonds(u, cdr1b_sele, mhca_sele,"cdr1b_mhca", step=10)
    # hbond, cdr2b_mhca_results, hbond_df6 = calc_Hbonds(u, cdr2b_sele, mhca_sele,"cdr2b_mhca", step=10)
    # hbond, cdr1a_mhcb_results, hbond_df7 = calc_Hbonds(u, cdr1a_sele, mhcb_sele,"cdr1a_mhcb", step=10)
    # hbond, cdr2a_mhcb_results, hbond_df8 = calc_Hbonds(u, cdr2a_sele, mhcb_sele,"cdr2a_mhcb", step=10)
    all_hbond_results = pd.concat([#mhca_pep_results, mhcb_pep_results, 
                                   tcra_pep_results, tcrb_pep_results],ignore_index=True)#,
                                   #cdr1b_mhca_results, cdr2b_mhca_results, cdr1a_mhcb_results, cdr2a_mhcb_results], ignore_index=True)
    all_hbond_results.to_csv(analysis_dir / "hbond_counts.csv", index=False)
    all_df = pd.concat([ hbond_df3, hbond_df4],ignore_index=True)#,hbond_df, hbond_df2, hbond_df5, hbond_df6, hbond_df7, hbond_df8], ignore_index=True)
    all_df.to_csv(analysis_dir / "hbond_details.csv", index=False)

    ######################################################
    ###########   Calculate per frame energy   ##############
    ######################################################
    #We will be using gmx mmpbsa instead so skip this for now
    #First grab slices of the trajectory for the thing
    # slices, times, frames = extract_frames(u, mhc_sele)
    # print(times)
    # print("Calculating Energy Data...")





    print("Finished Processing Directory: ", d)
    print("Time Elapsed: ", (time.time() - elapsed)/60, " minutes")
    print("Estimated time remaining: ", ((time.time() - elapsed)/60) * (len(dirs) - dirs.index(d) - 1), " minutes")
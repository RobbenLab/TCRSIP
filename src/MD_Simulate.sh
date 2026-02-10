#!/bin/bash
#This script will take a solvated structure and perform a molecular dynamics simulation and save the trajectories


#May have to run if can't find gromacs
source /usr/local/gromacs/bin/GMXRC
shopt -s extglob

# To run the script from base directory:
## ./scripts/MD_Simulate.sh ./Simulations/GroundTruth/1ao7
DIR=$1
# group="r ${2} \nname 17 mhc \nr ${2} \nname 18 tcr \nq"
# echo $group 
# TCR=$3
cd $DIR
rm !("tcr.pdb")
GMX_MAXCONSTRWARN=-1 #Suppress warnings about constraints

# simulate=1
# analyze=1
#Measure the time it takes to run the script
start_time=$(date +%s)

# if [ "$simulate" -eq 1] ; then
    
    #This script assumes that there is a folder structure for each pdb with a pdb file labeled tcr.pdb
    echo "processing the PDB file................."
    gmx pdb2gmx -f tcr.pdb -o processed.gro -water tip3p -ff "charmm27" -ignh #-chainsep id_or_ter #Last should keep chain info in the gro file Could throw off groupings
    #gmx editconf -f $groDir/raw/$name.gro -o $groDir/raw/$name.gro -c -d 1.2 
    gmx editconf -f processed.gro -o boxed.gro -c -box 13 #-bt dodecahedron #Put in a box #Use a larger box to avoid exploding?
    gmx solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top
    echo "Finished processing the PDB file........."

    #The first thing we need to do is add ions to the tpr
    echo "Adding Ions.............................."
    touch ions.mdp
    gmx grompp -f ions.mdp -c solvated.gro -p topol.top -o ions.tpr
    printf "SOL\n" | gmx genion -s ions.tpr -o solvated.gro -conc 0.15 -p topol.top -pname NA -nname CL -neutral
    echo "Finished adding Ions....................."

    #Next we need to do energy minimization
    echo "Starting Energy Minimization............."
    gmx grompp -f ../../Settings/emin-charmm.mdp -c solvated.gro -p topol.top -o em.tpr
    gmx mdrun -v -deffnm em -ntmpi 3 -ntomp 8 #use 24 threads
    echo "Finished Energy Minimization............."
    printf "Potential\n0\n" | gmx energy -f em.edr -o potential.xvg -xvg none #For downstream analysis in python

    #Next we can equalize the temperature
    echo "Equilibrating temperature................."
    gmx grompp -f ../../Settings/nvt-charmm.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
    gmx mdrun -ntmpi 1 -ntomp 24 -v -deffnm nvt
    echo "equilibrating temperature finished........"
    echo "Temperature" | gmx energy -f nvt.edr -o temperature.xvg -xvg none -b 20 #For analysis step downstream

    #That step equilibrates temperature so now we need to equilibrate pressure
    echo "Equilibrating pressure..................."
    gmx grompp -f ../../Settings/npt-charmm.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
    gmx mdrun -ntmpi 1 -ntomp 24 -v -deffnm npt
    echo "equilibrating pressure finished........"
    echo "Pressure" | gmx energy -f npt.edr -o pressure.xvg -xvg none -b 20 #For analysis step downstream

    #Now we can run the production simulation
    #Without any modifications, it takes ~40 minutes to run a 1 ns simulation at 2 fs steps
    echo "Running production simulation..........." 
    gmx grompp -f ../../Settings/md-mass-long-charmm.mdp -c npt.gro -r npt.gro -t npt.cpt -p topol.top -o md.tpr
    gmx mdrun -ntmpi 1 -ntomp 24 -v -deffnm md
    #Now rerun the simulation with the energy groups and -nb cpu to get the binding energy
    # printf $group | gmx make_ndx -f md.gro -o index.ndx
    # gmx grompp -f ../../Settings/md-mass-charmm-energy.mdp -c npt.gro -r npt.gro -t npt.cpt -p topol.top -o en.tpr -n index.ndx
    # gmx mdrun -ntmpi 1 -ntomp 24 -v -deffnm en -nb cpu -rerun md.xtc
    echo "Finished production simulation..........."

# fi

# if  [ "$analyze" -eq 1] ; then

    #Start the analysis of the simulation
    echo "Centering the simulation......."
    printf "1\n0\n" | gmx trjconv -s md.tpr -f md.xtc -o md_center.xtc -center -pbc mol #Output 0 is second option for whole system in centered trajectory
    echo "Finished analysis of the simulation......"

    #Now we can calculate the RMSD
    #Baseline, we are using 
    echo "Calculating RMSD........................"
    printf "4\n1\n" | gmx rms -s em.tpr -f md_center.xtc -o rmsd_xray.xvg -tu ns -xvg none
    echo "Finished calculating RMSD..............."

    #Now we can calculate the radius of gyration
    echo "Calculating Radius of Gyration..........."
    echo "1" | gmx gyrate -s em.tpr -f md_center.xtc -o radius_of_gyration.xvg -tu ns -xvg none
    echo "Finished calculating Radius of Gyration.."    

    # Now we can calculate the SASA
    echo "Calculating SASA........................"
    printf "4\n1\n" | gmx sasa -s em.tpr -f md_center.xtc -o sasa.xvg -tu ns -xvg none
    echo "Finished calculating SASA..............." 

    #Now we can calculate the hydrogen bonds
    echo "Calculating Hydrogen Bonds..............."
    printf "0\n0\n" | gmx hbond -s em.tpr -f md_center.xtc -num hbond.xvg -o hbond.ndx -dist hbond_dist.xvg -ang hbond_ang.xvg -tu ns -xvg none
    echo "Finished calculating Hydrogen Bonds......"

    #Calculate the energies
    echo "Calculating Energies....................."
    # printf "13" | gmx energy -f md.edr -o Energy.xvg -xvg none #Calculate total energy
    printf "9" | gmx energy -f md.edr -o Coloumb.xvg -xvg none -b 20 #Coloumbic energies
    printf "8" | gmx energy -f md.edr -o LJ.xvg -xvg none -b 20 #LJ energies
    echo "Finished calculating Energies............"

    #Now we can calculate the secondary structure
    # echo "Calculating Secondary Structure..........."   
    # printf "4\n1\n" | gmx do_dssp -s em.tpr -f md_center.xtc -o ss.xvg -tu ns -xvg none
    # echo "Finished calculating Secondary Structure.."

    #Now we can calculate the dihedral angles
    # echo "Calculating Dihedral Angles..............."
    # printf "4\n1\n" | gmx angle -s em.tpr -f md_center.xtc -o dihedral_angles.xvg -type dihedral -tu ns -xvg none
    # echo "Finished calculating Dihedral Angles......"

    # Now we can calculate the contacts (minimum distance for each atom pair)
    echo "Calculating Contacts......................"
    printf "4" | gmx mindist -s em.tpr -f md_center.xtc -o contacts.out -tu ns -xvg none
    echo "Finished calculating Contacts........."

# fi

#lastly report the methods:
echo "Simulation completed successfully!"
gmx report-methods -s md.tpr -o methods.txt



#Measure the time it took to run the script
end_time=$(date +%s)
execution_time=$((end_time - start_time))
echo "Total execution time: $execution_time seconds"



#!/bin/bash


#This script will take all directories in arg 1 and perform mmpbsa calculations on them

#Run using:
#scripts/Calc_mmpbsa.sh ./Simulations/GroundTruth/ 

#setup variables
DIR=$1
# cd $DIR
source /usr/local/gromacs/bin/GMXRC
shopt -s extglob


#Now write a for loop for each directory in DIR
for d in $DIR*/ ; do
    echo "Processing directory: $d"
    # dir="$d"
    # #check if dir exists
    if [ ! -d "$d" ]; then
        echo "Directory $d does not exist. Skipping..."
        continue
    fi
    mmpbsa_dir="$d/MMPBSA_Results/"
    #Make a new folder for the mmpbsa results
    mkdir -p $mmpbsa_dir
    #check if the mmpbsa results folder was created
    if [ ! -d "$mmpbsa_dir" ]; then
        echo "Failed to create directory $mmpbsa_dir. Skipping..."
        continue
    fi
    # head ${dir}md.xtc | echo
    # #Check if the md.gro files exist or use npt.gro instead
    groFile="$d/md.gro"
    if [ ! -f $groFile ]; then
        echo "File $groFile does not exist. using npt.gro instead."
        groFile="$d/npt.gro"
    fi
    #Check if md.xtc file exists
    xtcFile="$d/md.xtc"
    if [ ! -f $xtcFile ]; then
        echo "File $xtcFile does not exist. Skipping..."
        continue
    fi

    # #Now make the index file for the mmpbsa calculation
    # gmx make_ndx -f md.gro -o md_pbc.ndx
    # ri 0-189 &! r SOL &! r NA &! r CL
    # ri 191-450 &! r SOL &! r NA &! r CL
    # name 17 pmhc
    # name 18 tcr
    printf "ri 0-189 &! r SOL &! r NA &! r CL\nri 191-450 &! r SOL &! r NA &! r CL\nname 17 pmhc\nname 18 tcr\nq\n" | gmx make_ndx -f $groFile -o "$d/md_pbc.ndx"

    #Now redo the trajectory to remove the pbc
    # gmx trjconv -s md.tpr -f md.xtc -o md_pbc.xtc -center -n md_pbc.ndx -pbc noleap -ur compact -boxcenter rect
    printf "1\n0\n" | gmx trjconv -s "$d/md.tpr" -f "$d/md.xtc" -o "$d/md_pbc.xtc" -center -n "$d/md_pbc.ndx" -pbc nojump -ur compact -boxcenter rect

    # #Now run the mmpbsa calculation
    # #g_mmpbsa run -f md_pbc.xtc -s md.tpr -n md_pbc.ndx -unit1 pmhc -unit2 tcr -mme -mm energy_MM.xvg -decomp -mmcon contrib_MM.dat -dt 1 -tu ns -i ../../Settings/mmpbsa.mdp
    g_mmpbsa run -f "$d/md_pbc.xtc" -s "$d/md.tpr" -n "$d/md_pbc.ndx" -unit1 pmhc -unit2 tcr -mme -mm "$mmpbsa_dir/energy_MM.xvg" -decomp -mmcon "$mmpbsa_dir/contrib_MM.dat" -dt 1 -tu ns -i $DIR/../Settings/mmpbsa.mdp

done

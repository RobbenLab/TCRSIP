#!/bin/bash

#Tutorial site http://www.mdtutorials.com/gmx/lysozyme/05_EM.html
#The purpose of this script is to run gromacs on a single file to get the output hbonds
#There is alot more we can do with gromacs and simulations but that is not the objective at this time

#The script will take the name of the file pdb to be analyzed and the gromacs directory as input

#Must conda activate wsl_tcr

#Must run if first time:
#source /usr/local/gromacs/bin/GMXRC
source /usr/local/gromacs/bin/GMXRC
#Use: ./scripts/Get_HBonds.sh ./ModelData/KnownStructures/test/Structures/Norm/1bd2.pdb ./ModelData/KnownStructures/test/Structures/Gromacs
# ./scripts/Get_HBonds.sh ./ModelData/Immpred/Colabfold/Structures/Norm/Test_1 ./ModelData/Immpred/Colabfold/Structures/Gromacs

## ./scripts/Get_HBonds.sh "File" "Output Dir"

#First we need to get the file name for the pdb structure
file=$1
groDir=$2
# Hdir=$3

#Clear the temp files in gromacs
# rm -rf $groDir

# mkdir $groDir
# mkdir $groDir/raw
# mkdir $groDir/index
# mkdir $groDir/itp
# mkdir $groDir/processed
# mkdir $groDir/solvated
# mkdir $groDir/topology
# mkdir $groDir/charged
# mkdir $groDir/hydrogen


name=$(basename $file | cut -d"." -f1 ) #| cut -d "_" -f1)
echo $name

# rm -rf $Hdir/$name/

# mkdir $Hdir/$name/

rm -rf $groDir/hydrogen/$name
mkdir $groDir/hydrogen/$name

#Now we can first convert the pdb into a gromacs compatible file
echo "**************Beginning conversion******************"
# gmx pdb2gmx -f $file -o $groDir/raw/$name.gro -water spc -p $groDir/topology/$name.top -ff charmm27 -i $groDir/topology/$name.itp -n $groDir/index/$name.ndx
gmx pdb2gmx -f $file -o $groDir/raw/$name.gro -water tip3p -p $groDir/topology/$name.top -i $groDir/topology/$name.itp -ff charmm27
#So the problem is that the molecule is off center so when we place a box to solavte it does not cover use the following code to center the molecule
#https://gromacs.bioexcel.eu/t/molecules-slightly-out-of-solvent-box/408
#gmx editconf -f molecule.pdb -o boxed.pdb -c -d 1.2 -bt octahedron
gmx editconf -f $groDir/raw/$name.gro -o $groDir/raw/$name.gro -c -d 1.2 


#Now take the gromacs file and solvate it with water
echo "*******************Beginning solvation******************"
gmx solvate -cp $groDir/raw/$name.gro -o $groDir/solvated/$name.gro -p $groDir/topology/$name.top

#Combine the solvated and topology into one gromacs file
echo "*******************Beginning Processing******************"
# gmx pdb2gmx -f $groDir/solvated/$name.gro -o $groDir/processed/$name.processed.gro -water spc -p $groDir/topology/$name.top -ff charmm27 -i $groDir/topology/$name.itp -ignh -missing -ter
# gmx pdb2gmx -f $groDir/solvated/$name.gro -o $groDir/processed/$name.processed.gro -p $groDir/topology/$name.top -ff charmm27 

# gmx make_ndx -f $groDir/processed/$name.processed.gro -o $groDir/index/$name.ndx #Don't need index on these
# gmx make_ndx -f $groDir/raw/$name.gro -o $groDir/index/$name.ndx
gmx grompp -f ./settings/normal.mdp -c $groDir/solvated/$name.gro -p $groDir/topology/$name.top -o $groDir/processed/$name.tpr -maxwarn 500


#Now calculate all of the hydrogen bonds in the file
echo "*******************Calculating Hbonds******************"
# gmx hbond -f $groDir/processed/$name.processed.gro -n $groDir/index/$name.ndx -s $groDir/topology/$name.top -o $groDir/hbond/$name.hbond.ndx
# gmx hbond -f $groDir/solvated/$name.gro -n $groDir/index/$name.ndx -s $groDir/processed/$name.tpr -o $Hdir/$name/hbond.ndx -num $Hdir/$name/num.xvg -dist $Hdir/$name/dist.xvg -ang $Hdir/$name/ang.xvg
printf "1\n1\n" |gmx hbond-legacy -f $groDir/solvated/$name.gro -s $groDir/processed/$name.tpr -hbn $groDir/hydrogen/$name/hbond.ndx -hbm $groDir/hydrogen/$name/hbond.xpm -num $groDir/hydrogen/$name/num.xvg -dist $groDir/hydrogen/$name/dist.xvg -ang $groDir/hydrogen/$name/ang.xvg -merge no

#Output the gro file with charges
gmx editconf -f $groDir/processed/$name.tpr -o $groDir/charged/$name.pdb -grasp

#Output the gro with hydrogens
# gmx editconf -f $groDir/raw/$name.gro -o $groDir/hydrogen/$name/bonded.pdb

#Energy minimization function
# gmx mdrun -v -s $groDir/processed/$name.tpr                 #-deffnm em

#Calculate the coloumbic energies
# gmx energy -s $groDir/processed/$name.tpr -o $groDirenergy.xvg





###Base code used
# gmx pdb2gmx -f input.pdb -o output.gro -water spc
# gmx solvate -cp output.gro -cs spc216.gro -o output_solv.gro -p output.top
# gmx pdb2gmx -f output_solv.gro -o output_processed.gro -p output.top -ignh
# gmx make_ndx -f output_processed.gro -o index.ndx
# gmx hbond -f output_processed.gro -n index.ndx -s output.tpr -num hbond.xvg -hbn hbond.ndx

#Electrostatic interactions?
#http://www.mdtutorials.com/gmx/complex/09_analysis.html
#Coulombic interaction lennard jones energy
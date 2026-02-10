#!/bin/bash

#This script will take in the input of the aligned dir and the template dir and will run dockq on each protein 
#must be run in an environment where DockQ has been installed (pip install DockQ)

alnDir=$1
refDir=$2
jsonDir=$3

echo 'starting DockQ run on'
echo $alnDir

for f in $alnDir/*  
do
    name=$(basename $f |cut -d"." -f1 | cut -d "_" -f1)
    DockQ $alnDir/$name* $refDir/$name* --json $jsonDir/$name.json
done
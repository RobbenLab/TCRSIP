#!/bin/bash
##The following code will take the data in the input argument folder and run the full analysis and feature extraction pipeline.
##Results will be stored in the output argument folder

#Check if the correct number of positional arguments is provided
# if [ "$#" -ne 2 ]; then
#     echo "Usage: $0 <input_folder> <output_folder>"
#     exit 1
# fi

Rosetta=true
Gromacs=true

output_folder="./output"  #Default output folder
#Assign input and output folder variables
input_folder="$1"  
# output_folder="$2"

#get options for running no rosetta or gromacs
while [[ $# -gt 0 ]]; do
  case "$1" in
    -r|--no-rosetta)
      Rosetta=false #Will not run the rosetta energy calculations and will not extract rosetta features
      shift # shift once since flag has no argument
      ;;
    -g|--no-gromacs)
      Gromacs=false #Will not run the gromacs hbond calculations and will not extract gromacs features
      shift 2 # shift twice: once for flag, once for argument
      ;;
    -o|--output)
      output_folder="$2"
      shift 2 # shift twice: once for flag, once for argument
      ;;
    -*)
      echo "Unknown option $1"
      exit 1
      ;;
    *)
      # Handle positional arguments
      POSITIONAL_ARGS+=("$1")
      shift
      ;;
  esac
done


#Check if the input folder exists
if [ ! -d "$input_folder" ]; then
    echo "Error: Input folder '$input_folder' does not exist."
    exit 1
fi

#Detect if pyrosetta is installed
if ! python -c "import pyrosetta" &> /dev/null; then
    echo "Error: Pyrosetta is not installed. Running without."
    Rosetta=false
fi

#Detect if gromacs is installed
source /usr/local/gromacs/bin/GMXRC
if ! command -v gmx &> /dev/null; then
    echo "Error: Gromacs is not installed. Running without."
    Gromacs=false
fi

#Create the output folder if it does not exist and the file structure then transfer the contents of the input directory
mkdir -p "$output_folder"
mkdir -p "$output_folder/Structures/Raw"
cp -r "$input_folder/"* "$output_folder/Structures/Raw/"

#Run the full analysis pipeline with the Analyze_Folded.py script in the ./src folder
#Use the -a argument if the output is from alphafold to extract pLDDT from the bscore column of pdb
python -Wi ./src/Analyze_Folded.py $output_folder -r $Rosetta -gro $Gromacs 

#Now run the feature extraction pipeline with the Extract_Features.py script in the ./src folder
python ./src/ExtractFeatures.py --Dir $output_folder --Out $output_folder -ta 0 -ma 0

#Now run the provided R script to plot the features and do some basic analysis
#Rscript ./scripts/PlotFeatures.R $input_folder $output_folder
# TCRSIP
TCRSIP is a tool built upon AlphaFold2 based structural inference to predict TCR-pMHC specificity from inferred structures. The end-to-end model is available through a google colab notebook based webserver instance to be rapidly deployable by researchers. 

***Colab Notebook under construction***

The tool requires sequence information for the full length or truncated $\alpha$$\beta$TCR, the peptide sequence (valid amino acid only), and the full-length or truncated MHC molecule (Class 1 or 2). Sequences must include the variable region of both chains of the TCR and the interacting region of the MHC (Will not accept B2M). The model has only been validated on human TCR and not on gamma delta TCR or CD1 molecules. Putative TCR-pMHC can be run in batches and it is recommended to provide a google drive folder to save intermediate datasteps as each sample can take up to 20 min to run. 

The code hosted in this repository can be run to process and extract structural features in TCR-pMHC structures. Code to reproduce specific analyses from manuscript or simulation runs can also be found in the ./scripts directory.

## Installation

Requirements:
TCRSIP pipeline requires a version of the conda package manager to be installed.

To install the TCRSIP pipeline, clone the directory to a local directory:
```
git clone https://github.com/RobbenLab/TCRSIP
```
Create a new conda environment using the provided environment.yml file
```
conda env create -f environment.yml
```
You can then activate the environment with the following command:
```
conda activate new_env
```
Must also install PyMol:
```
conda install -c conda-forge pymol-open-source 
```
And install PyRosetta from website (requires license): 
https://www.pyrosetta.org/downloads

And install gromacs to /usr/local/ as specificied in the instructions:
https://manual.gromacs.org/current/install-guide/index.html 

## Running Processing

To run the processing pipeline, use the run.sh command from the main directory and point it towards a directory containing structures in the required format.
```
./{Folder}/Structures/Raw/....
```
Run.sh will run a pipeline that takes a folder of TCRpMHC pdb files as input and then extract features used for analysis and prediction of interaction.
```
bash run.sh ./data/example/ -o ./output
```
This code will run the 10 pdb's in the example folder and will output in the output directory specified. 
```
Run.sh
 
Usage: bash Run.sh /path/to/input/dir/ [options]
  -o | --output (./output) Output directory
  -r | --no-rosetta (off) Flag to turn off rosetta binding energy calculation (This saves time, ~7 min per file, but results in less predictive features) using Rosetta is on by default
  -g | --no-gromacs (off) Flag to turn off gromacs enabled Hbond calculations, using Gromacs is on by default

The command will run the **Analyze_Folded.py** and **ExtractFeatures.py** scripts to generate features from the input pdb files. Input directory is required to Temporary files and finished features will be output into the input directory. 
```

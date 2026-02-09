# TCRSIP
TCRSIP is a tool built upon AlphaFold2 based structural inference to predict TCR-pMHC specificity from inferred structures. The full model is available through a google colab notebook based webserver instance to be rapidly deployable by researchers. 

Initial commit

## Installation

Requirements:
TCRSIP pipeline requires a version of the conda package manager to be installed.

To install the TCRSIP pipeline, clone the directory to a local directory:


Create a new conda environment using the provided environment.yml file

conda env create -f environment.yml

You can then activate the environment with the following command:

conda activate new_env

Must also install PyMol:

conda install -c conda-forge pymol-open-source 

And install PyRosetta from website (requires license): 
https://www.pyrosetta.org/downloads


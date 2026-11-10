# TCRSIP
TCRSIP is a tool built upon AlphaFold2 based structural inference to predict TCR-pMHC specificity from inferred structures. The end-to-end model is available through a google colab notebook based webserver instance to be rapidly deployable by researchers. 

| | | Date Published | Change list | Status
| :--- | :--- | :--- | :--- | :--- |
| TCRSIP v1.2 | [link](https://colab.research.google.com/drive/1pQ95DXi0mVWrrmYC9HviP0dFkgwWS6gF?usp=sharing) | 5-26-26 | Added interactive dashboard, ranking of alphafold outputs, redid the outputs, added structural evaluations | Running |
| TCRSIP v1.1 | [link](https://colab.research.google.com/drive/1L8cF_c1xyjnIfIzjC_8hQLdj0JQe27vr?usp=sharing) | 3-29-26 | | Running |


The tool requires sequence information for the full length or truncated alpha beta TCR, the peptide sequence (valid amino acid only), and the full-length or truncated MHC molecule (Class 1 or 2). Sequences must include the variable region of both chains of the TCR and the interacting region of the MHC (Will not accept B2M). The model has only been validated on human TCR and not on gamma delta TCR or CD1 molecules. Putative TCR-pMHC can be run in batches and it is recommended to provide a google drive folder to save intermediate datasteps as each sample can take up to 20 min to run. 

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

Installing gromacs to /usr/local/ as specificied in the instructions is now deprecated but required for replicating dynamics simulations:
https://manual.gromacs.org/current/install-guide/index.html 

## Running Feature Extraction and Prediction Pipeline

To build the per-structure data folders the current pipeline runs on, use `Run.sh` from the main directory, pointing it at a directory of raw TCR-pMHC `.pdb` files and a CSV of reference chain sequences:
```
bash Run.sh ./data/example/ -s sequences.csv -o ./test
```
```
Run.sh

Usage: bash Run.sh [/path/to/raw/structures/ -s sequences.csv] [options]
  -s | --seqfile <file>    (required with a structure directory) CSV with reference chain
                           sequences per structure (columns: ID, MHCa, MHCb, Peptide, TCRa, TCRb)
  -o | --output <dir>      (.) Top-level output folder
  -p | --prefix <prefix>   Optional prefix prepended to every structure's folder/name
  --ml                     (off) Also train the full ML models (all of TCRclassifier.py's
                           DEFAULT_MODELS) and predict
  --dl                     (off) Also run the full DL (ResNet) model and predict
  --checkpoint <file>      (./model_checkpoints/full_v1.pt) DL checkpoint (--dl only)
```
`-o/--output` is a single top-level folder (default `./out`) - everything `Run.sh` writes goes into it or its subfolders: `<output>/results.csv`, `<output>/dat` (per-structure folders, features, predictions) and `<output>/input` (DL input maps, `--dl` only).

This runs `src/SetupData.py`, which matches each raw structure's chains to the reference sequences by similarity (so the original PDB chain ids/order don't matter), renames/renumbers them consistently, and writes a `<output>/dat/<name>/` folder per structure (`Structure/Complex.pdb` plus a `<name>.json` of extracted regions/annotations). It then runs `src/ExtractFeaturestoJSON.py` on just the structures it added (via the `new_structures.txt` manifest `SetupData.py` writes), producing `<output>/dat/new_features.csv`.

**The structure directory and `-s`/`--seqfile` may both be omitted** if `<output>/dat` already has every structure's features extracted (i.e. each `<output>/dat/<name>/<name>.json` already has a `Features` key from a prior run). This skips `SetupData.py` and `ExtractFeaturestoJSON.py` entirely and instead compiles the combined features CSV from what's already stored in each structure's json via `src/CompileFeaturesCSV.py` (no recomputation, no biopython/pyrosetta calls). `--ml`/`--dl`/the summary below then act on every structure in `<output>/dat`, instead of just a newly-added batch.

With `--ml`, `src/InferML.py` trains every model in `TCRclassifier.py`'s `DEFAULT_MODELS` (Logistic Regression, Random Forest, Gradient Boosting, Support Vector Machine, K-Nearest Neighbors, Multilayer Perceptron - each with its own tuned hyperparameters from `data/ML/Params/best_params.json`) on every currently labeled structure in `data/Features/features.csv`/`data/ML/Annotation.csv`, and predicts with each on `new_features.csv`, writing one row per structure with every model's `<Model>_Prob`/`<Model>_Pred` columns to `<output>/dat/ml_predictions.csv`.

With `--dl`, `src/Resnet/InferDL.py` builds the 2D input maps for the new structures (`src/Resnet/build_inputs.py`, into `<output>/input`) and runs the trained ResNet checkpoint (`--checkpoint`, default `./model_checkpoints/full_v1.pt` - expected alongside a `model_info.json` and `channel_stats.npz` from the same training run, same as `src/Resnet/run_models.py` writes) on them, writing `<output>/dat/dl_predictions.csv`.

Either flag also writes its results into each predicted structure's own `<output>/dat/<name>/<name>.json`, under a top-level `Inference` key (alongside `Scoring`/`Features`) - one sub-key per model (e.g. `Inference["Random Forest"]`, `Inference.DL`), each `{"Prob", "Pred", "Model"}`. Re-running a flag only overwrites its own models' keys, leaving every other model's (and the other flag's) results in place.

Both flags are off by default and assume this installation does not have NetTCR set up, so only our own ML and DL models are covered.

A summary is always printed at the end (`src/SummarizeStructures.py`, over the same batch - new or, in the bypass case above, every structure in `<output>/dat`): total structure count, the MHC class I/II mix, how many look potentially misfolded (`Scoring.Final` below a threshold), each model's binder call rate from whatever `Inference.<Model>` keys are present, and an ensemble call (mean probability across every model present) combining all of them. The same script also writes `<output>/results.csv`, one row per structure: `name`, `Peptide`, `CDR3a`/`CDR3b`, `MHC_class`, `Misfold_Score`/`Misfolded`, and every model's (and the ensemble's) `Prob`/`Pred`.

Every step (`SetupData.py`/`CompileFeaturesCSV.py`, `ExtractFeaturestoJSON.py`, `InferML.py`, `InferDL.py`, `SummarizeStructures.py`) logs a timestamp to stdout when it starts and when it finishes.

### Legacy pipeline (deprecated)

`Deprecated_Run.sh` is an older pipeline that extracts a flat feature CSV directly from raw structures (via `Analyze_Folded.py`/`ExtractFeatures.py`) without producing the `./dat/<name>/` folders above. It's incompatible with the current ML/DL/NetTCR pipeline and is kept only for backward compatibility with older feature.csv-based workflows.
```
bash Deprecated_Run.sh ./data/example/ -o ./output
```
```
Deprecated_Run.sh

Usage: bash Deprecated_Run.sh /path/to/input/dir/ [options]
  -o | --output (./output) Output directory
  -r | --no-rosetta (off) Flag to turn off rosetta binding energy calculation (This saves time, ~7 min per file, but results in less predictive features) using Rosetta is on by default
  -g | --no-gromacs (off) Flag to turn off gromacs enabled Hbond calculations, using Gromacs is on by default
```

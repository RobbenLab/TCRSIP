#!/bin/bash
## Builds the ./dat/<name>/ structure+JSON folders the current pipeline runs on, from a
## directory of raw TCR-pMHC .pdb structures plus a CSV of reference chain sequences, then
## extracts their features. This is the current (dat/JSON-based) pipeline - see
## src/SetupData.py and src/ExtractFeaturestoJSON.py. Deprecated_Run.sh is the older,
## feature.csv-based pipeline (Analyze_Folded.py/ExtractFeatures.py) and is deprecated but
## kept for reference/backward compatibility.
##
## ML/DL inference against the newly processed structures is optional (off by default,
## --ml/--dl to enable) and assumes this installation does not have NetTCR set up, so it
## only covers our own ML (src/InferML.py) and DL (src/Resnet/InferDL.py) models.

output_dir="./out"  # Top-level output folder; everything below goes into subfolders of it
prefix=""
run_ml=false
run_dl=false
checkpoint="./model_checkpoints/full_v1.pt"  # DL checkpoint (--dl only) - not an output, left as-is

step_start=0
log_start() {
    step_start=$SECONDS
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}
log_end() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1 ($((SECONDS - step_start))s)"
}

usage() {
  echo "Usage: bash Run.sh [/path/to/raw/structures/ -s sequences.csv] [options]"
  echo "  -o/--output is one top-level folder; everything this script writes goes into it"
  echo "  or its subfolders - <output>/results.csv, <output>/dat (per-structure folders,"
  echo "  features, predictions) and <output>/input (DL input maps)."
  echo "  The structure directory and -s/--seqfile may both be omitted if <output>/dat"
  echo "  already has every structure's features extracted - this skips SetupData.py and"
  echo "  feature extraction, and --ml/--dl/the summary act on everything already there"
  echo "  instead of just a newly-added batch."
  echo "  -s | --seqfile <file>    (required with a structure directory) CSV with reference"
  echo "                           chain sequences per structure (columns: ID, MHCa, MHCb,"
  echo "                           Peptide, TCRa, TCRb)"
  echo "  -o | --output <dir>      (./out) Top-level output folder"
  echo "  -p | --prefix <prefix>   Optional prefix prepended to every structure's folder/name"
  echo "  --ml                     (off) Also train every ML model and predict"
  echo "  --dl                     (off) Also run the full DL (ResNet) model and predict"
  echo "  --checkpoint <file>      (./model_checkpoints/full_v1.pt) DL checkpoint (--dl only)"
  echo "  NetTCR is not run here - this installation is assumed not to have it set up."
  echo "  A summary (structure count, MHC class mix, potential misfolds, model/ensemble"
  echo "  binder calls) is always printed at the end, and written per-structure to"
  echo "  <output>/results.csv (name, peptide, CDR3a/b, MHC class, misfold score, every"
  echo "  model's + the ensemble's Prob/Pred). Every step logs a timestamp when it"
  echo "  starts and finishes."
  exit 1
}

POSITIONAL_ARGS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    -s|--seqfile)
      seqfile="$2"
      shift 2
      ;;
    -o|--output)
      output_dir="$2"
      shift 2
      ;;
    -p|--prefix)
      prefix="$2"
      shift 2
      ;;
    --ml)
      run_ml=true
      shift
      ;;
    --dl)
      run_dl=true
      shift
      ;;
    --checkpoint)
      checkpoint="$2"
      shift 2
      ;;
    -h|--help)
      usage
      ;;
    -*)
      echo "Unknown option $1"
      usage
      ;;
    *)
      POSITIONAL_ARGS+=("$1")
      shift
      ;;
  esac
done
set -- "${POSITIONAL_ARGS[@]}"
input_structures_dir="$1"

dat_dir="$output_dir/dat"
input_dir="$output_dir/input"
features_csv="$dat_dir/new_features.csv"
if [ -n "$input_structures_dir" ]; then
    if [ ! -d "$input_structures_dir" ]; then
        echo "Error: Input folder '$input_structures_dir' does not exist."
        exit 1
    fi
    if [ -z "$seqfile" ]; then
        echo "Error: -s/--seqfile is required when a structure directory is given."
        usage
    fi
    if [ ! -f "$seqfile" ]; then
        echo "Error: Sequence file '$seqfile' does not exist."
        exit 1
    fi

    cmd=(python ./src/SetupData.py "$input_structures_dir" -s "$seqfile" -o "$dat_dir")
    if [ -n "$prefix" ]; then
        cmd+=(-p "$prefix")
    fi
    log_start "Matching raw structures in $input_structures_dir against $seqfile and writing per-structure folders to $dat_dir (SetupData.py)..."
    "${cmd[@]}"
    log_end "Done - new structures listed in $dat_dir/new_structures.txt (SetupData.py)"

    # Feature extraction, restricted to the structures SetupData.py just processed (its
    # new_structures.txt manifest) rather than every structure ever added to $dat_dir.
    names_file="$dat_dir/new_structures.txt"
    log_start "Extracting sequence/structure/physics features for the newly added structures (ExtractFeaturestoJSON.py)..."
    python ./src/ExtractFeaturestoJSON.py --Dir "$dat_dir" --Out "$features_csv" --names-file "$names_file"
    log_end "Done - features written to $features_csv (ExtractFeaturestoJSON.py)"
else
    echo "No structure directory given - skipping SetupData.py and feature extraction;"
    echo "using the features already extracted into $dat_dir."
    names_file=""
    log_start "Compiling already-extracted features for every structure in $dat_dir, without recomputing them (CompileFeaturesCSV.py)..."
    python ./src/CompileFeaturesCSV.py --Dir "$dat_dir" --Out "$features_csv"
    log_end "Done - features written to $features_csv (CompileFeaturesCSV.py)"
fi

names_file_arg=()
if [ -n "$names_file" ]; then
    names_file_arg=(--names-file "$names_file")
fi

if [ "$run_ml" = true ]; then
    log_start "Training every ML model on all labeled structures and predicting on the new structures (InferML.py)..."
    python ./src/InferML.py --new-features "$features_csv" --dat-dir "$dat_dir" \
        -o "$dat_dir/ml_predictions.csv"
    log_end "Done - predictions written to $dat_dir/ml_predictions.csv and each structure's Inference.<model> (InferML.py)"
fi

if [ "$run_dl" = true ]; then
    log_start "Building DL input maps and running the ResNet checkpoint ($checkpoint) on the new structures (InferDL.py)..."
    python ./src/Resnet/InferDL.py --dat-dir "$dat_dir" "${names_file_arg[@]}" \
        --input-dir "$input_dir" --checkpoint "$checkpoint" -o "$dat_dir/dl_predictions.csv"
    log_end "Done - predictions written to $dat_dir/dl_predictions.csv and each structure's Inference.DL (InferDL.py)"
fi

log_start "Summarizing all processed structures and compiling results.csv (SummarizeStructures.py)..."
python ./src/SummarizeStructures.py --dat-dir "$dat_dir" "${names_file_arg[@]}" -o "$output_dir/results.csv"
log_end "Done - results written to $output_dir/results.csv (SummarizeStructures.py)"

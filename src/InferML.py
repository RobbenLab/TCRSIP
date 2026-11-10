#!/usr/bin/env python3

# Trains every ML model in TCRclassifier.py's DEFAULT_MODELS (Logistic Regression, Random
# Forest, Gradient Boosting, Support Vector Machine, K-Nearest Neighbors, Multilayer
# Perceptron) on every engineered feature, with each model's tuned hyperparameters from
# data/ML/Params/best_params.json (the FullRun/All comparison, see
# data/ML/Results/FullRun/FullRun/All/test_summary_metrics.csv), then predicts a binder
# probability per model for newly processed structures that have no label yet.
#
# This is the deploy-time counterpart to TCRclassifier.py's evaluate_models/
# evaluate_fixed_split, which are for comparing models/feature sets against each other
# (train/test splits, CV, importance) rather than producing single models fit on
# everything to apply to new structures.
#
#   python src/InferML.py --new-features ./dat_out/new_features.csv -o ./ml_predictions.csv

import argparse
import json
import os

import pandas as pd
from sklearn.base import clone
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC

from TCRclassifier import DEFAULT_MODELS, getMisfolded, paramsToModels

ANNOTATION_CSV = "./data/ML/Annotation.csv"
FEATURES_CSV = "./data/Features/features.csv"
PARAMS_FILE = "./data/ML/Params/best_params.json"
DAT_DIR = "./dat"
MODEL_NAMES = list(DEFAULT_MODELS)
MISFOLD_THRESHOLD = 0.8
ID_COL, LABEL_COL = "name", "response"


def train_full_models(annotation_csv, features_csv, params_file, dat_dir, model_names, misfold_threshold):
    response_df = pd.read_csv(annotation_csv)
    feature_df = pd.read_csv(features_csv)

    misfolded = getMisfolded(dat_dir, threshold=misfold_threshold)
    feature_df = feature_df[~feature_df[ID_COL].isin(misfolded)]

    feature_columns = [c for c in feature_df.columns if c != ID_COL]
    merged = response_df.merge(feature_df, on=ID_COL, how="inner").dropna(subset=[LABEL_COL])
    X = merged[feature_columns].apply(pd.to_numeric, errors="coerce").reset_index(drop=True)
    y = merged[LABEL_COL].astype(int).reset_index(drop=True)
    print(f"Training {len(model_names)} model(s) on {len(X)} labeled structures "
         f"({int(y.sum())} binders), {len(feature_columns)} features")

    with open(params_file) as f:
        best_params = json.load(f)
    models = paramsToModels(best_params, DEFAULT_MODELS)

    pipes = {}
    for name in model_names:
        base = models[name]
        if isinstance(base, SVC):
            # Needed for predict_proba (Platt scaling via internal 5-fold CV); off by
            # default in DEFAULT_MODELS/best_params.json since the FullRun comparison
            # only needed decision_function's raw score for ROC_AUC, not a calibrated
            # probability. Doesn't change the predicted class, only adds this estimate.
            base = clone(base).set_params(probability=True)
        pipe = Pipeline([
            ("imputer", SimpleImputer(strategy="median", keep_empty_features=True)),
            ("scaler", StandardScaler()),
            ("model", clone(base)),
        ])
        pipe.fit(X, y)
        pipes[name] = pipe
        print(f"  trained {name!r}")
    return pipes, feature_columns


def write_inference_json(dat_dir, preds, model_name, key=None):
    """Record each prediction into its structure's <dat_dir>/<name>/<name>.json, under
    Inference.<key> - alongside Scoring/Features, so a structure's folder is the one place
    to look up everything known about it, not just this run's output CSV."""
    key = key or model_name
    missing = []
    for _, row in preds.iterrows():
        name = row[ID_COL]
        json_path = os.path.join(dat_dir, name, f"{name}.json")
        if not os.path.isfile(json_path):
            missing.append(name)
            continue
        with open(json_path) as f:
            js = json.load(f)
        js.setdefault("Inference", {})[key] = {
            "Prob": float(row["Prob"]),
            "Pred": int(row["Pred"]),
            "Model": model_name,
        }
        with open(json_path, "w") as f:
            json.dump(js, f, indent=4)
    if missing:
        print(f"Note: no {dat_dir}/<name>/<name>.json for {len(missing)} structure(s), "
             f"Inference.{key} not recorded: {missing[:10]}")


def predict_all(pipes, feature_columns, new_features_csv, threshold):
    new_df = pd.read_csv(new_features_csv)
    # A column can be legitimately absent for every row here (e.g. SASA:MHCb:* when this
    # batch is all MHC class I) without being a problem - reindex adds it back as NaN, and
    # each trained pipeline's imputer (fit on the full historical data) fills it with that
    # column's training-set median, same as it does for class-I rows during training.
    missing = set(feature_columns) - set(new_df.columns)
    if missing:
        print(f"Note: {new_features_csv} has no values for {len(missing)} feature(s) "
             f"(will be median-imputed): {sorted(missing)[:10]}{'...' if len(missing) > 10 else ''}")
    X_new = new_df.reindex(columns=feature_columns).apply(pd.to_numeric, errors="coerce")

    preds_by_model = {}
    for name, pipe in pipes.items():
        proba = pipe.predict_proba(X_new)[:, 1]
        preds_by_model[name] = pd.DataFrame({
            ID_COL: new_df[ID_COL],
            "Prob": proba,
            "Pred": (proba >= threshold).astype(int),
        })
    return preds_by_model


def wide_table(preds_by_model):
    wide = None
    for name, df in preds_by_model.items():
        cols = df.rename(columns={"Prob": f"{name}_Prob", "Pred": f"{name}_Pred"})
        wide = cols if wide is None else wide.merge(cols, on=ID_COL)
    return wide


def main():
    ap = argparse.ArgumentParser(
        description="Train every ML model on all labeled structures and predict on newly processed ones"
    )
    ap.add_argument("--new-features", required=True,
                    help="features CSV (from ExtractFeaturestoJSON.py) for the new, unlabeled structures")
    ap.add_argument("--annotation", default=ANNOTATION_CSV, help="labeled structures' response/source table")
    ap.add_argument("--features", default=FEATURES_CSV, help="labeled structures' feature matrix")
    ap.add_argument("--params", default=PARAMS_FILE, help="tuned hyperparameters (best_params.json)")
    ap.add_argument("--dat-dir", default=DAT_DIR, help="./dat/<name>/<name>.json folders, to exclude misfolded structures")
    ap.add_argument("--models", nargs="+", default=MODEL_NAMES, choices=MODEL_NAMES,
                    help="which classifiers to train and predict with (default: all)")
    ap.add_argument("--misfold-threshold", type=float, default=MISFOLD_THRESHOLD)
    ap.add_argument("--threshold", type=float, default=0.5, help="probability cutoff for the predicted class")
    ap.add_argument("-o", "--output", default="./ml_predictions.csv")
    args = ap.parse_args()

    pipes, feature_columns = train_full_models(
        args.annotation, args.features, args.params, args.dat_dir, args.models, args.misfold_threshold
    )
    preds_by_model = predict_all(pipes, feature_columns, args.new_features, args.threshold)

    wide = wide_table(preds_by_model)
    wide.to_csv(args.output, index=False)
    for name, preds in preds_by_model.items():
        write_inference_json(args.dat_dir, preds, name)

    print(f"Saved {len(wide)} structures' predictions from {len(preds_by_model)} model(s) to {args.output} "
         f"and to each structure's Inference.<model>")


if __name__ == "__main__":
    main()

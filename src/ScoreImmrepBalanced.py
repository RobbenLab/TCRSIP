#!/usr/bin/env python3

# Resamples each pipeline's existing "trained on core (or externally, for ERGO/NetTCR's
# paper-trained run), tested on immrep" per-sample predictions down to a ~50/50
# positive/negative balance - immrep is heavily skewed negative (~598 positive of 3484
# overall) - and recomputes the same Accuracy/ROC_AUC/PR_AUC/PPV/NPV metrics on each
# balanced resample. The point is to separate how much of the core->immrep
# generalization gap is simply class imbalance from how much is a genuine shift in
# peptide/MHC allele/TCR gene distribution between core and immrep.
#
# Each pipeline predicted on a slightly different immrep subset to begin with (ML/DL use
# the structured subset with a folded sample in ./dat, ERGO and NetTCR's paper-trained
# run cover the complete immrep set, NetTCR's own-trained run is filtered further by its
# own sequence-length limits) - this script doesn't try to align identical sample sets
# across pipelines, consistent with how every other cross-pipeline comparison in this
# project already tolerates that.
#
# For every (Method, Model, source fold) pool of predictions, N_RESAMPLES independent
# balanced subsamples are drawn (the majority class downsampled to match the minority
# class count) and scored separately - so even NetTCR/ERGO/each DL fold, which otherwise
# have no inherent repeated-measurement structure, end up with a real mean +- SD / point
# spread instead of a single fixed value.

import os

import numpy as np
import pandas as pd

from ScorePredictions import compute_binary_metrics

ANNOTATION_CSV = "./data/ML/Annotation.csv"
ML_PREDICTIONS = "./data/ML/Results/Base/non_immpred_to_immpred/All/test_predictions.csv"
DL_DIR = "./data/DL/Core_on_immrep"
DL_N_FOLDS = 5
NETTCR_PREDICTIONS = "./data/NetTCR/Comparisons/core_vs_immrep/predictions.csv"
NETTCR_PAPER_PREDICTIONS = "./data/NetTCR/Comparisons/paper_train_vs_immrep_full/predictions.csv"
ERGO_DIR = "./data/ERGO"
ERGO_RESULTS = {
    "ERGO-Auto": "ImmpredAutoresults.csv",
    "ERGO-LSTM": "ImmpredLMSTresults.csv",
}

OUT_FILE = "./data/Results/immrep_balanced_resamples.csv"
N_RESAMPLES = 50
RANDOM_STATE = 42


def balanced_resample_metrics(y_true, y_score, n_resamples, random_state):
    # Repeatedly downsamples the majority class to the minority class's count (without
    # replacement) and scores each resulting ~50/50 subset - the minority class itself is
    # also "resampled" every time (a full-size draw without replacement, i.e. just a
    # reshuffle) purely so both classes go through the same indexing path.
    y_true = np.asarray(y_true)
    y_score = np.asarray(y_score)
    pos_idx = np.flatnonzero(y_true == 1)
    neg_idx = np.flatnonzero(y_true == 0)
    n_minority = min(len(pos_idx), len(neg_idx))
    if n_minority == 0:
        return []

    rng = np.random.default_rng(random_state)
    rows = []
    for i in range(n_resamples):
        pos_sample = rng.choice(pos_idx, size=n_minority, replace=False)
        neg_sample = rng.choice(neg_idx, size=n_minority, replace=False)
        idx = np.concatenate([pos_sample, neg_sample])
        metrics = compute_binary_metrics(y_true[idx], y_score[idx])
        metrics["Repeat"] = i
        rows.append(metrics)
    return rows


def load_ml():
    df = pd.read_csv(ML_PREDICTIONS)
    pools = []
    for (model, fold), g in df.groupby(["Model", "Fold"]):
        pools.append(("ML", model, f"fold{fold}", g["response"].values, g["Probability"].values))
    return pools


def load_dl():
    pools = []
    for i in range(DL_N_FOLDS):
        path = os.path.join(DL_DIR, f"Fold{i}predictions.csv")
        df = pd.read_csv(path)
        pools.append(("DL", "ResNet", f"fold{i}", df["Label"].values, df["Prob"].values))
    return pools


def load_nettcr():
    # Two distinct NetTCR-2.0 runs, same Method (color) but explicitly different Model
    # labels so they aren't conflated: one trained on our own core set and tested on the
    # structured immrep subset (core_vs_immrep), the other trained purely on NetTCR-2.0's
    # own published training set and tested on the complete immrep test set
    # (paper_train_vs_immrep_full, see PrepNetTCR.py) - no exposure to our data at all.
    ours = pd.read_csv(NETTCR_PREDICTIONS)
    paper = pd.read_csv(NETTCR_PAPER_PREDICTIONS)
    return [
        ("NetTCR", "NetTCR-2.0 (ours)", "run0", ours["binder"].values, ours["prediction"].values),
        ("NetTCR", "NetTCR-2.0 (paper)", "run0", paper["binder"].values, paper["prediction"].values),
    ]


def load_ergo():
    # ERGO's own output carries a raw score but no label at all - attached positionally
    # from Annotation.csv's immrep rows, same approach as ScoreERGO.py (confirmed there:
    # Peptide matches 1:1 in row order).
    annotation = pd.read_csv(ANNOTATION_CSV)
    immrep = annotation[annotation["source"] == "immrep"].reset_index(drop=True)
    pools = []
    for model_name, filename in ERGO_RESULTS.items():
        preds = pd.read_csv(os.path.join(ERGO_DIR, filename))
        if len(preds) != len(immrep) or not (preds["Peptide"].values == immrep["Peptide"].values).all():
            raise ValueError(
                f"{filename} no longer lines up row-for-row with Annotation.csv's immrep rows."
            )
        pools.append(("ERGO", model_name, "run0", immrep["response"].values, preds["Score"].values))
    return pools


def main():
    os.makedirs(os.path.dirname(OUT_FILE), exist_ok=True)
    pools = load_ml() + load_dl() + load_nettcr() + load_ergo()

    rows = []
    for i, (method, model, fold, y_true, y_score) in enumerate(pools):
        for metrics in balanced_resample_metrics(y_true, y_score, N_RESAMPLES, random_state=RANDOM_STATE + i):
            metrics["Method"] = method
            metrics["Model"] = model
            metrics["Fold"] = fold
            rows.append(metrics)

    result = pd.DataFrame(rows)[
        ["Method", "Model", "Fold", "Repeat", "n", "n_positive", "Accuracy", "ROC_AUC", "PR_AUC", "PPV", "NPV"]
    ]
    result.to_csv(OUT_FILE, index=False)
    print(f"Wrote {len(result)} resample rows ({N_RESAMPLES} resamples x {len(pools)} (model, fold) pools) to {OUT_FILE}")
    print(
        result.groupby(["Method", "Model"])[["Accuracy", "ROC_AUC", "PR_AUC", "PPV", "NPV"]]
        .agg(["mean", "std"])
        .to_string()
    )


if __name__ == "__main__":
    main()

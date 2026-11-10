#!/usr/bin/env python3

# Builds a per-sample x per-model table of prediction outcomes (TP/TN/FP/FN) for the
# "Core+ImmRep (train/test)" FullRun comparison - the one fixed split (data/DL/FullRun/
# split.csv) all three pipelines were scored against, so unlike every other comparison in
# this project, the ML/DL/NetTCR held-out sets actually overlap sample-for-sample. Each
# pipeline additionally filtered that fixed split down in its own way though (ML drops
# misfolded structures, NetTCR drops non-class-I/too-long sequences), so this script
# restricts to the intersection of samples every pipeline actually produced a prediction
# for, rather than assuming the three sets are identical.
#
# A soft-voting "Ensemble" column is added alongside the individual models - the mean of
# the ML+DL models' predicted probability for that sample (NetTCR-2.0 excluded from the
# vote itself, though still covered by the same common-sample restriction and shown as
# its own column), thresholded at 0.5 and scored the same way as every other model, so
# its summary row is directly comparable.
#
# Each sample also gets a per-attribute "similarity to training" score (Allele, Peptide,
# TRBV - see compute_similarity) meant as a row annotation alongside the outcome
# heatmap: how common was this sample's MHC allele, peptide, and TRBV gene in what every
# model was actually trained on.

import os

import pandas as pd

from PrepNetTCR import ANNOTATION_CSV, FULLRUN_SPLIT_CSV, build_nettcr_table
from ScorePredictions import compute_binary_metrics

ML_PREDICTIONS = "./data/ML/Results/FullRun/FullRun/All/test_predictions.csv"
DL_PREDICTIONS = "./data/DL/FullRun/predictions.csv"
NETTCR_TEST = "./data/NetTCR/Comparisons/full_run/test.csv"
NETTCR_PREDICTIONS = "./data/NetTCR/Comparisons/full_run/predictions.csv"
SIMILARITY_ATTRS = ["Allele", "Peptide", "TRBV"]

OUT_DIR = "./data/Results"
OUTCOMES_FILE = os.path.join(OUT_DIR, "fullrun_sample_outcomes.csv")
SUMMARY_FILE = os.path.join(OUT_DIR, "fullrun_sample_outcomes_summary.csv")
SIMILARITY_FILE = os.path.join(OUT_DIR, "fullrun_sample_similarity.csv")


def classify(true, pred):
    if true == 1 and pred == 1:
        return "TP"
    if true == 0 and pred == 0:
        return "TN"
    if true == 0 and pred == 1:
        return "FP"
    return "FN"


def load_ml():
    df = pd.read_csv(ML_PREDICTIONS)
    return df.rename(columns={"response": "true", "Predicted": "pred", "Probability": "score"})[
        ["name", "Model", "true", "pred", "score"]
    ]


def load_dl():
    df = pd.read_csv(DL_PREDICTIONS)
    df = df[df["Set"] == "val"]
    df = df.rename(columns={"Name": "name", "response": "true", "Pred": "pred", "Prob": "score"})
    df["Model"] = "ResNet"
    return df[["name", "Model", "true", "pred", "score"]]


def load_nettcr():
    # predictions.csv carries no sample id at all - recover it by reconstructing the same
    # (annotation -> CDR3 extraction/filtering -> merge with the fixed split -> val rows)
    # pipeline PrepNetTCR.py used to build full_run/test.csv in the first place, then
    # confirming the row order still lines up before trusting the positional join.
    annotation = pd.read_csv(ANNOTATION_CSV)
    table = build_nettcr_table(annotation)
    core = table[table["source"].isin(["ground", "exp", "fake"])]
    immrep = table[table["source"] == "immrep"]
    core_immrep = pd.concat([core, immrep], ignore_index=True)
    fullrun_split = pd.read_csv(FULLRUN_SPLIT_CSV)
    fullrun_table = core_immrep.merge(fullrun_split[["name", "Set"]], on="name", how="inner")
    test = fullrun_table[fullrun_table["Set"] == "val"].reset_index(drop=True)

    preds = pd.read_csv(NETTCR_PREDICTIONS)
    if len(preds) != len(test) or not (
        test[["CDR3a", "CDR3b", "peptide", "binder"]].reset_index(drop=True)
        == preds[["CDR3a", "CDR3b", "peptide", "binder"]].reset_index(drop=True)
    ).all().all():
        raise ValueError(
            f"{NETTCR_PREDICTIONS} no longer lines up row-for-row with the reconstructed "
            "full_run test set - can't recover sample names this way anymore."
        )

    df = test[["name"]].copy()
    df["true"] = preds["binder"].values
    df["score"] = preds["prediction"].values
    df["pred"] = (df["score"] >= 0.5).astype(int)
    df["Model"] = "NetTCR-2.0 (ours)"
    return df[["name", "Model", "true", "pred", "score"]]


def compute_similarity(sample_names, annotation):
    # Per-attribute normalized similarity-to-training: for each of Allele/Peptide/TRBV,
    # count how many training samples shared that exact value, then normalize by the
    # largest such count seen among OUR samples specifically (not the training set's
    # overall modal value, which may belong to a value none of our samples even have) -
    # whichever of our samples has the most training support for that attribute scores
    # 1, a value never seen in training scores 0. Missing values are treated as their
    # own category rather than dropped, consistent with how bias_column CV grouping
    # handles them elsewhere in this project.
    fullrun_split = pd.read_csv(FULLRUN_SPLIT_CSV)
    train_names = fullrun_split[fullrun_split["Set"] == "train"]["name"]
    train_ann = annotation[annotation["name"].isin(train_names)].copy()
    test_ann = annotation[annotation["name"].isin(sample_names)].copy()

    result = test_ann[["name"]].copy()
    for attr in SIMILARITY_ATTRS:
        train_values = train_ann[attr].fillna("(missing)")
        test_values = test_ann[attr].fillna("(missing)")
        counts = train_values.value_counts()
        sample_counts = test_values.map(counts).fillna(0)
        result[attr] = (sample_counts / sample_counts.max()).values

    return result[["name"] + SIMILARITY_ATTRS]


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    model_predictions = pd.concat([load_ml(), load_dl(), load_nettcr()], ignore_index=True)

    # Each pipeline filtered the shared fixed split down differently (see module
    # docstring) - only samples every pipeline actually scored are usable here.
    names_per_model = model_predictions.groupby("Model")["name"].apply(set)
    common_names = set.intersection(*names_per_model)
    print(f"Samples per model: {dict(names_per_model.apply(len))}")
    print(f"Common to all {len(names_per_model)} models: {len(common_names)}")
    model_predictions = model_predictions[model_predictions["name"].isin(common_names)]

    # Soft-voting ensemble: mean predicted probability across the ML+DL models only
    # (NetTCR-2.0 is excluded - it's a far weaker, differently-trained model on a much
    # smaller/differently-filtered sample set, so folding it into the vote would mostly
    # just dilute the stronger models rather than add useful signal), thresholded at 0.5
    # and scored identically to every individual model below. It still only covers the
    # samples common to every model (NetTCR included), since it's compared against
    # NetTCR as its own column in the heatmap.
    ensemble = (
        model_predictions[model_predictions["Model"] != "NetTCR-2.0 (ours)"]
        .groupby("name")
        .agg(true=("true", "first"), score=("score", "mean"))
        .reset_index()
    )
    ensemble["pred"] = (ensemble["score"] >= 0.5).astype(int)
    ensemble["Model"] = "Ensemble"

    combined = pd.concat([model_predictions, ensemble[["name", "Model", "true", "pred", "score"]]], ignore_index=True)
    combined["Outcome"] = [classify(t, p) for t, p in zip(combined["true"], combined["pred"])]

    annotation = pd.read_csv(ANNOTATION_CSV)
    combined = combined.merge(annotation[["name", "source"]], on="name", how="left")

    similarity = compute_similarity(common_names, annotation)
    similarity.to_csv(SIMILARITY_FILE, index=False)
    combined = combined.merge(similarity, on="name", how="left")

    combined = combined[["name", "source", "Model", "true", "score", "pred", "Outcome"] + SIMILARITY_ATTRS]
    combined.to_csv(OUTCOMES_FILE, index=False)

    summary_rows = []
    for model_name, g in combined.groupby("Model"):
        metrics = compute_binary_metrics(g["true"], g["score"])
        metrics["Model"] = model_name
        summary_rows.append(metrics)
    summary = pd.DataFrame(summary_rows)[["Model", "n", "n_positive", "Accuracy", "ROC_AUC", "PR_AUC", "PPV", "NPV"]]
    summary = summary.sort_values("Accuracy", ascending=False)
    summary.to_csv(SUMMARY_FILE, index=False)

    print(f"Wrote {len(combined)} rows ({len(common_names)} samples x {combined['Model'].nunique()} models) to {OUTCOMES_FILE}")
    print(summary.to_string(index=False))
    print(f"Saved summary to {SUMMARY_FILE}")
    print(similarity[SIMILARITY_ATTRS].describe().loc[["mean", "min", "max"]].to_string())
    print(f"Saved similarity to {SIMILARITY_FILE}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3

# Checks whether the FullRun Ensemble's prediction confidence correlates with any of the
# engineered features in data/Features/features.csv. "Confidence" here is how far the
# Ensemble's predicted probability sits from 0.5 (scaled to 0-1), not the probability
# itself - a sample scored 0.95 or 0.05 should count the same way (very confident),
# which the raw probability alone doesn't capture. Correlated with every feature column
# across the 161 samples the Ensemble covers (see BuildFullRunSampleOutcomes.py), using
# both Pearson (linear) and Spearman (rank, more robust to outliers/non-linearity) so a
# feature showing up under one but not the other is still visible. Each feature column
# is median-imputed first, the same strategy evaluate_models()/PlotResults.R's PCA both
# use, so a column that's entirely NaN for these samples (and therefore constant after
# imputation) naturally drops out via a NaN correlation rather than crashing.

import os

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

OUTCOMES_CSV = "./data/Results/fullrun_sample_outcomes.csv"
FEATURES_CSV = "./data/Features/features.csv"

OUT_DIR = "./data/Results"
OUT_FILE = os.path.join(OUT_DIR, "ensemble_confidence_feature_correlation.csv")

TOP_N = 20


def main():
    outcomes = pd.read_csv(OUTCOMES_CSV)
    ensemble = outcomes[outcomes["Model"] == "Ensemble"][["name", "score"]].copy()
    ensemble["confidence"] = (ensemble["score"] - 0.5).abs() * 2

    features = pd.read_csv(FEATURES_CSV)
    merged = ensemble.merge(features, on="name", how="inner")
    print(f"Samples with both an Ensemble prediction and a feature row: {len(merged)} of {len(ensemble)}")

    feature_cols = [c for c in features.columns if c != "name"]
    confidence = merged["confidence"].values

    rows = []
    for col in feature_cols:
        values = pd.to_numeric(merged[col], errors="coerce")
        if values.isna().all():
            values = values.fillna(0)
        else:
            values = values.fillna(values.median())
        if values.std() == 0:
            continue

        pearson_r, pearson_p = pearsonr(values, confidence)
        spearman_r, spearman_p = spearmanr(values, confidence)
        rows.append(
            {
                "Feature": col,
                "Pearson_r": pearson_r,
                "Pearson_p": pearson_p,
                "Spearman_r": spearman_r,
                "Spearman_p": spearman_p,
            }
        )

    result = pd.DataFrame(rows)
    result["AbsPearson_r"] = result["Pearson_r"].abs()
    result = result.sort_values("AbsPearson_r", ascending=False).drop(columns="AbsPearson_r")
    result.to_csv(OUT_FILE, index=False)

    print(f"\nTop {TOP_N} features by |Pearson r| with Ensemble confidence:")
    print(result.head(TOP_N).to_string(index=False))
    print(f"\nSaved all {len(result)} feature correlations to {OUT_FILE}")


if __name__ == "__main__":
    main()

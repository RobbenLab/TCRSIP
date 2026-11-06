#!/usr/bin/env python3

# This script classifies TCR-pMHC interactions from the features extracted by
# ExtractFeaturestoJSON.py. `evaluate_models` runs a set of classifiers through
# stratified k-fold cross validation on a given subset of samples/features and
# writes metrics + comparison plots to a directory, so it can be called
# repeatedly on different subsets of the same underlying dataset (e.g. by
# MHC class, allele, or data source) without duplicating the CV/plotting logic.

import json
import os
import time
import warnings
from contextlib import contextmanager

import numpy as np
import pandas as pd
import matplotlib

# Force the non-interactive Agg backend before pyplot is used - this script only ever
# saves figures to disk, never displays them, and running with an interactive backend
# (e.g. TkAgg) alongside joblib's multiprocess workers (n_jobs=-1 in
# permutation_importance/RandomizedSearchCV) corrupts Tcl/Tk's global state across the
# fork/thread boundary, causing "Tcl_AsyncDelete: async handler deleted by the wrong
# thread", "RuntimeError: main thread is not in main loop", and leaked
# resource_tracker folder warnings at shutdown.
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.stats import randint, uniform, loguniform

# Some features (e.g. MHCb-chain measurements) are structurally undefined for
# class I MHC samples and are all-NaN within a fold; SimpleImputer just drops
# those columns for that fold, which is expected here rather than a problem.
warnings.filterwarnings("ignore", message="Skipping features without any observed values")

from sklearn.base import clone
from sklearn.model_selection import (
    StratifiedKFold,
    StratifiedGroupKFold,
    train_test_split,
    RandomizedSearchCV,
)
from sklearn.pipeline import Pipeline
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.inspection import permutation_importance
from sklearn.metrics import (
    accuracy_score,
    precision_score,
    recall_score,
    f1_score,
    roc_auc_score,
    average_precision_score,
    roc_curve,
    precision_recall_curve,
    confusion_matrix,
)

from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.svm import SVC
from sklearn.neighbors import KNeighborsClassifier
from sklearn.neural_network import MLPClassifier

# ==========================================
# DEFAULT CLASSIFIERS
# ==========================================
DEFAULT_MODELS = {
    "Logistic Regression": LogisticRegression(max_iter=1000, random_state=42),
    "Random Forest": RandomForestClassifier(random_state=42),
    "Gradient Boosting": GradientBoostingClassifier(random_state=42),
    "Support Vector Machine": SVC(random_state=42),
    "K-Nearest Neighbors": KNeighborsClassifier(),
    "Multilayer Perceptron": MLPClassifier(max_iter=500, random_state=42),
}

# Random-search hyperparameter spaces for evaulateParameters, keyed to match DEFAULT_MODELS.
# Each value is a sklearn param_distributions dict (list = discrete choices, distribution
# object = sampled continuously) for that model's own parameters - no pipeline "model__"
# prefix, that's added internally.
DEFAULT_PARAM_DISTRIBUTIONS = {
    "Logistic Regression": {
        "C": loguniform(1e-3, 1e3),
        "class_weight": [None, "balanced"],
    },
    "Random Forest": {
        "n_estimators": randint(100, 500),
        "max_depth": [None, 5, 10, 15, 20, 30],
        "min_samples_split": randint(2, 20),
        "min_samples_leaf": randint(1, 10),
        "max_features": ["sqrt", "log2", None],
        "class_weight": [None, "balanced"],
    },
    "Gradient Boosting": {
        "n_estimators": randint(50, 300),
        "learning_rate": loguniform(1e-3, 0.3),
        "max_depth": randint(2, 8),
        "subsample": uniform(0.6, 0.4),
        "min_samples_leaf": randint(1, 10),
    },
    "Support Vector Machine": {
        "C": loguniform(1e-2, 1e3),
        "kernel": ["rbf", "linear", "poly"],
        "gamma": loguniform(1e-4, 1e1),
        "class_weight": [None, "balanced"],
    },
    "K-Nearest Neighbors": {
        "n_neighbors": randint(3, 30),
        "weights": ["uniform", "distance"],
        "p": [1, 2],
    },
    "Multilayer Perceptron": {
        "hidden_layer_sizes": [(50,), (100,), (50, 50), (100, 50), (100, 100)],
        "alpha": loguniform(1e-5, 1e-1),
        "learning_rate_init": loguniform(1e-4, 1e-1),
        "activation": ["relu", "tanh"],
    },
}

METRIC_COLUMNS = ["Accuracy", "Precision", "Recall", "F1", "ROC_AUC", "PR_AUC", "FPR", "FNR", "PPV", "NPV"]
# Raw confusion-matrix counts. Kept out of METRIC_COLUMNS since that list also drives the
# metric-comparison bar plot, which assumes every value is a 0-1 rate.
COUNT_COLUMNS = ["TP", "TN", "FP", "FN"]


# ==========================================
# TIMING
# ==========================================
class _Timer:
    """Accumulates wall-clock time spent in named components across the whole run."""

    def __init__(self):
        self.totals = {}

    def add(self, key, seconds):
        self.totals[key] = self.totals.get(key, 0.0) + seconds

    @contextmanager
    def time(self, key):
        start = time.perf_counter()
        yield
        self.add(key, time.perf_counter() - start)

    def report(self, run_name, detail="full"):
        if detail not in ("none", "short", "full"):
            raise ValueError(f"timing_detail must be one of 'none', 'short', 'full'; got {detail!r}")
        if detail == "none":
            return
        total = sum(self.totals.values())
        if detail == "short":
            print(f"'{run_name}' total time: {total:.2f}s")
            return
        print(f"\nTiming breakdown for '{run_name}':")
        for key, seconds in sorted(self.totals.items(), key=lambda kv: -kv[1]):
            pct = (seconds / total * 100) if total else 0.0
            print(f"  {key:<26s} {seconds:9.2f}s  ({pct:5.1f}%)")
        print(f"  {'TOTAL':<26s} {total:9.2f}s")


# ==========================================
# MAIN ENTRY POINT
# ==========================================
def evaluate_models(
    name,
    validation_ids,
    feature_columns,
    response_df,
    feature_df,
    output_dir,
    test_ids=None,
    test_size=0.2,
    bias_column=None,
    models=None,
    id_col="name",
    label_col="response",
    n_splits=5,
    random_state=42,
    importance_method="both",
    perm_n_repeats=5,
    perm_scoring="roc_auc",
    perm_n_jobs=-1,
    perm_max_samples=1.0,
    timing_detail="full",
):
    """
    Run stratified k-fold cross validation for a dict of classifiers on
    `validation_ids`, and write metrics/predictions/feature importance/plots
    to `output_dir`, prefixed "validation_". A test set is always evaluated
    too: each fold's fitted pipeline is additionally scored against a fixed,
    held-out set of ids never touched by the CV split, with its own set of
    "test_"-prefixed outputs - useful for a true generalization check
    alongside the CV estimate.

    validation_ids: list of sample identifiers (values in `id_col`) to run
        stratified k-fold CV over. If `test_ids` is None, a `test_size`
        fraction of this set is first split off (stratified by label) to use
        as the test set, and the CV only runs over the remainder.
    test_ids: optional list of sample identifiers held out from the CV split
        entirely; each fold's fitted pipeline is scored against this same
        fixed set. If None (default), the test set is instead a random
        `test_size` fraction of `validation_ids`. Feature importance is only
        computed for the validation folds, not for the test set.
    test_size: fraction of `validation_ids` to hold out as the test set when
        `test_ids` is not given (default 0.2). Ignored if `test_ids` is given.
    bias_column: optional column in `response_df`. The train/test split (explicit
        test_ids, or the random test_size split) happens first either way and does not
        consider bias_column at all - a group can end up on both sides. bias_column
        only affects how the *CV folds within the train side* are built: if None
        (default), plain StratifiedKFold with exactly `n_splits` folds (random,
        label-balanced). If given, every sample sharing a value of this column is kept
        in a single fold via StratifiedGroupKFold - and since a fixed fold count can
        land on very unevenly sized folds depending on how the group sizes happen to
        fall, `n_splits` instead becomes an upper bound: fold counts from 2 up to
        `n_splits` (capped at the number of unique bias_column values) are tried, and
        whichever minimizes the gap between the largest and smallest fold is used. The
        number of unique values and the chosen fold sizes are printed, and a
        "validation_<bias_column>_fold_composition.csv" matrix (fold x bias_column
        value -> sample count) is written to `output_dir`.
    feature_columns: list of columns in `feature_df` to use as predictors.
    response_df: dataframe containing `id_col` and `label_col`.
    feature_df: dataframe containing `id_col` and `feature_columns`.
    models: dict of {name: unfitted sklearn estimator}; defaults to DEFAULT_MODELS.
    output_dir: directory under which a `name` subdirectory is created for output.
    importance_method: one of "native", "permutation", "both", or "none".
        "native" uses coef_/feature_importances_ (fast, only Logistic
        Regression/Random Forest/Gradient Boosting support it). "permutation"
        is model-agnostic (works for all models) but is by far the most
        expensive component - see the timing report printed at the end.
        "both" (default) computes whichever apply; "none" skips importance
        entirely.
    perm_n_repeats/perm_scoring/perm_n_jobs/perm_max_samples: passed through to
        sklearn.inspection.permutation_importance when permutation importance
        is computed.
    timing_detail: one of "full" (default, per-component breakdown), "short"
        (just the total run time), or "none" (no timing report printed).

    Returns (validation_fold_results, validation_summary) dataframes.
    """
    if importance_method not in ("native", "permutation", "both", "none"):
        raise ValueError(
            f"importance_method must be one of 'native', 'permutation', 'both', 'none'; got {importance_method!r}"
        )
    if timing_detail not in ("none", "short", "full"):
        raise ValueError(
            f"timing_detail must be one of 'none', 'short', 'full'; got {timing_detail!r}"
        )
    timer = _Timer()
    ##################################################
    print(f"Evaluating models for {name}")
    # out_dir = os.path.join(output_dir, name)
    out_dir = output_dir
    os.makedirs(out_dir, exist_ok=True)
    ##################################################
    models = models if models is not None else DEFAULT_MODELS
    with timer.time("data_prep"):
        # 1. Train/test split, same either way regardless of bias_column - explicit
        # test_ids are used as-is, otherwise a plain random test_size split off
        # validation_ids. No group-awareness here; a bias_column group can end up on
        # both sides of the train/test boundary.
        if test_ids is not None:
            X, y, sample_ids, bias_values = _build_dataset(
                validation_ids, feature_columns, response_df, feature_df, id_col, label_col, bias_column
            )
            X_test, y_test, ids_test, _ = _build_dataset(
                test_ids, feature_columns, response_df, feature_df, id_col, label_col
            )
        else:
            X_pool, y_pool, ids_pool, bias_pool = _build_dataset(
                validation_ids, feature_columns, response_df, feature_df, id_col, label_col, bias_column
            )
            split_arrays = [X_pool, y_pool, ids_pool] if bias_column is None else [X_pool, y_pool, ids_pool, bias_pool]
            split_result = train_test_split(
                *split_arrays, test_size=test_size, stratify=y_pool, random_state=random_state
            )
            if bias_column is None:
                X, X_test, y, y_test, sample_ids, ids_test = split_result
                bias_values = None
            else:
                X, X_test, y, y_test, sample_ids, ids_test, bias_values, _ = split_result

            X, y, sample_ids = X.reset_index(drop=True), y.reset_index(drop=True), sample_ids.reset_index(drop=True)
            X_test, y_test, ids_test = X_test.reset_index(drop=True), y_test.reset_index(drop=True), ids_test.reset_index(drop=True)
            if bias_values is not None:
                bias_values = bias_values.reset_index(drop=True)

        # 2. CV folds within the train side only. Plain StratifiedKFold, unless
        # bias_column is given, in which case StratifiedGroupKFold keeps every sample
        # sharing a bias_column value in the same fold.
        if bias_column is None:
            folds = list(StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=random_state).split(X, y))
        else:
            # A NaN bias_column value isn't "the same group" as any other NaN row, and
            # StratifiedGroupKFold can't sort NaN against real values anyway - give each
            # one its own unique placeholder so it's treated as its own group of one.
            missing_placeholders = pd.Series(
                [f"__missing_{bias_column}_{i}__" for i in bias_values.index], index=bias_values.index
            )
            bias_values = bias_values.fillna(missing_placeholders)
            folds = _choose_group_folds(X, y, bias_values, bias_column, n_splits, random_state)

    if bias_column is not None:
        with timer.time("bias_column_composition"):
            fold_of_sample = np.empty(len(X), dtype=int)
            for fold_idx, (_, val_idx) in enumerate(folds, start=1):
                fold_of_sample[val_idx] = fold_idx
            composition = pd.crosstab(pd.Series(fold_of_sample, name="Fold"), bias_values.rename(bias_column))
            composition.to_csv(os.path.join(out_dir, f"validation_{bias_column}_fold_composition.csv"))
    ##################################################
    val_state = _init_run_state(models)
    test_state = _init_run_state(models)
    importance_rows = []
    ##################################################
    print("Training models...")
    for fold_idx, (train_idx, val_idx) in enumerate(tqdm(folds), start=1):
        X_train, X_val = X.iloc[train_idx], X.iloc[val_idx]
        y_train, y_val = y.iloc[train_idx], y.iloc[val_idx]
        ids_val = sample_ids.iloc[val_idx]

        for model_name, estimator in models.items():
            pipe = Pipeline(
                [
                    ("imputer", SimpleImputer(strategy="median", keep_empty_features=True)),
                    ("scaler", StandardScaler()),
                    ("model", clone(estimator)),
                ]
            )
            with timer.time("model_fit_predict"):
                pipe.fit(X_train, y_train)
                y_pred = pipe.predict(X_val)
                y_score = _decision_scores(pipe, X_val)

            with timer.time("validation_score_record"):
                _record_predictions(
                    val_state, model_name, fold_idx, ids_val, y_val, y_pred, y_score, id_col, label_col
                )

            with timer.time("test_score_record"):
                y_pred_test = pipe.predict(X_test)
                y_score_test = _decision_scores(pipe, X_test)
                _record_predictions(
                    test_state, model_name, fold_idx, ids_test, y_test, y_pred_test, y_score_test, id_col, label_col
                )

            if importance_method in ("native", "both"):
                with timer.time("native_importance"):
                    importance_rows.extend(
                        _native_importance_rows(pipe, X_val.columns, model_name, fold_idx)
                    )

            if importance_method in ("permutation", "both"):
                with timer.time("permutation_importance"):
                    importance_rows.extend(
                        _permutation_importance_rows(
                            pipe,
                            X_val,
                            y_val,
                            model_name,
                            fold_idx,
                            n_repeats=perm_n_repeats,
                            scoring=perm_scoring,
                            n_jobs=perm_n_jobs,
                            max_samples=perm_max_samples,
                            random_state=random_state,
                        )
                    )
    ####################################################
    print("Summarizing folds...")
    validation_fold_results, validation_summary = _finalize_run(val_state, out_dir, "validation", timer)
    _finalize_run(test_state, out_dir, "test", timer)

    feature_importance = pd.DataFrame(importance_rows)
    with timer.time("summarize_importance"):
        importance_summary = _summarize_feature_importance(feature_importance)
    with timer.time("validation_write_csvs"):
        feature_importance.to_csv(os.path.join(out_dir, "validation_feature_importance.csv"), index=False)
        importance_summary.to_csv(os.path.join(out_dir, "validation_feature_importance_summary.csv"), index=False)

    timer.report(name, detail=timing_detail)

    return validation_fold_results, validation_summary


# ==========================================
# FIXED TRAIN/TEST SPLIT (no CV)
# ==========================================
def evaluate_fixed_split(
    name,
    train_ids,
    test_ids,
    feature_columns,
    response_df,
    feature_df,
    output_dir,
    models=None,
    id_col="name",
    label_col="response",
    importance_method="native",
    perm_n_repeats=5,
    perm_scoring="roc_auc",
    perm_n_jobs=-1,
    perm_max_samples=1.0,
    timing_detail="short",
):
    """
    Fit each model once on `train_ids` and score it once on `test_ids` - no inner k-fold
    CV, unlike evaluate_models. Meant for comparisons against an externally fixed
    train/test assignment (e.g. the DL pipeline's data/DL/FullRun/split.csv) where the
    split itself must match exactly rather than be rebuilt by evaluate_models' own CV/
    random-split logic.

    Writes the same "test_"-prefixed fold_metrics/summary_metrics/predictions/plot
    outputs evaluate_models writes for its fixed test set (fold index is always 1, since
    there's only one fit per model), so createSummary reads it identically. Also writes
    "feature_importance.csv"/"feature_importance_summary.csv" when importance_method
    requests importance.

    Returns (test_fold_results, test_summary) dataframes.
    """
    if importance_method not in ("native", "permutation", "both", "none"):
        raise ValueError(
            f"importance_method must be one of 'native', 'permutation', 'both', 'none'; got {importance_method!r}"
        )
    timer = _Timer()
    print(f"Evaluating fixed split for {name}")
    out_dir = output_dir
    os.makedirs(out_dir, exist_ok=True)
    models = models if models is not None else DEFAULT_MODELS

    with timer.time("data_prep"):
        X_train, y_train, _, _ = _build_dataset(
            train_ids, feature_columns, response_df, feature_df, id_col, label_col
        )
        X_test, y_test, ids_test, _ = _build_dataset(
            test_ids, feature_columns, response_df, feature_df, id_col, label_col
        )

    test_state = _init_run_state(models)
    importance_rows = []

    print("Training models...")
    for model_name, estimator in tqdm(models.items()):
        pipe = Pipeline(
            [
                ("imputer", SimpleImputer(strategy="median", keep_empty_features=True)),
                ("scaler", StandardScaler()),
                ("model", clone(estimator)),
            ]
        )
        with timer.time("model_fit_predict"):
            pipe.fit(X_train, y_train)
            y_pred = pipe.predict(X_test)
            y_score = _decision_scores(pipe, X_test)

        with timer.time("test_score_record"):
            _record_predictions(
                test_state, model_name, 1, ids_test, y_test, y_pred, y_score, id_col, label_col
            )

        if importance_method in ("native", "both"):
            with timer.time("native_importance"):
                importance_rows.extend(
                    _native_importance_rows(pipe, X_train.columns, model_name, 1)
                )

        if importance_method in ("permutation", "both"):
            with timer.time("permutation_importance"):
                importance_rows.extend(
                    _permutation_importance_rows(
                        pipe,
                        X_test,
                        y_test,
                        model_name,
                        1,
                        n_repeats=perm_n_repeats,
                        scoring=perm_scoring,
                        n_jobs=perm_n_jobs,
                        max_samples=perm_max_samples,
                        random_state=42,
                    )
                )

    print("Summarizing...")
    test_fold_results, test_summary = _finalize_run(test_state, out_dir, "test", timer)

    if importance_rows:
        feature_importance = pd.DataFrame(importance_rows)
        with timer.time("summarize_importance"):
            importance_summary = _summarize_feature_importance(feature_importance)
        with timer.time("write_importance_csvs"):
            feature_importance.to_csv(os.path.join(out_dir, "feature_importance.csv"), index=False)
            importance_summary.to_csv(os.path.join(out_dir, "feature_importance_summary.csv"), index=False)

    timer.report(name, detail=timing_detail)
    return test_fold_results, test_summary


# ==========================================
# HYPERPARAMETER SEARCH
# ==========================================
def evaluateParameters(
    ids,
    feature_columns,
    response_df,
    feature_df,
    output_dir,
    models=None,
    param_distributions=None,
    id_col="name",
    label_col="response",
    n_splits=5,
    n_iter=25,
    scoring="roc_auc",
    random_state=42,
    n_jobs=-1,
    verbose=0,
):
    """
    Runs a randomized hyperparameter search (RandomizedSearchCV) per model over `ids`,
    inside the same impute->scale->model pipeline evaluate_models uses, and returns a
    dict of {model_name: best_params} - each model's best hyperparameters as plain
    estimator.set_params()-style kwargs (no pipeline "model__" prefix), which can be
    turned into a `models` dict for evaluate_models via paramsToModels().

    Writes to output_dir: best_params.json (the returned dict), plus one
    "<model_name>_search_results.csv" per searched model recording every candidate
    tried and its CV score (sklearn's search.cv_results_), for later reference.

    ids: sample identifiers to search over (same convention as evaluate_models).
    feature_columns/response_df/feature_df/id_col/label_col: same as evaluate_models.
    output_dir: directory to write best_params.json and the per-model search-results
        CSVs to; created if it doesn't exist.
    models: dict of {name: unfitted estimator}; defaults to DEFAULT_MODELS.
    param_distributions: dict of {name: sklearn param_distributions dict}; defaults to
        DEFAULT_PARAM_DISTRIBUTIONS. A model name with no entry (or an empty one) is
        skipped - its default parameters are kept as-is and it's left out of the result.
    n_splits: number of CV folds used to score each candidate.
    n_iter: number of random parameter combinations to try per model.
    scoring: sklearn scoring string used to rank candidates.
    """
    os.makedirs(output_dir, exist_ok=True)
    models = models if models is not None else DEFAULT_MODELS
    param_distributions = (
        param_distributions if param_distributions is not None else DEFAULT_PARAM_DISTRIBUTIONS
    )

    X, y, _, _ = _build_dataset(ids, feature_columns, response_df, feature_df, id_col, label_col)
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=random_state)

    best_params = {}
    for model_name, estimator in models.items():
        distributions = param_distributions.get(model_name)
        if not distributions:
            print(f"No param_distributions for '{model_name}', skipping search (keeping defaults).")
            continue

        pipe = Pipeline(
            [
                ("imputer", SimpleImputer(strategy="median", keep_empty_features=True)),
                ("scaler", StandardScaler()),
                ("model", clone(estimator)),
            ]
        )
        prefixed_distributions = {f"model__{k}": v for k, v in distributions.items()}

        search = RandomizedSearchCV(
            pipe,
            param_distributions=prefixed_distributions,
            n_iter=n_iter,
            scoring=scoring,
            cv=skf,
            random_state=random_state,
            n_jobs=n_jobs,
            verbose=verbose,
        )
        print(f"Searching hyperparameters for {model_name} ({n_iter} candidates x {n_splits} folds)...")
        search.fit(X, y)

        best_params[model_name] = {
            k.removeprefix("model__"): v for k, v in search.best_params_.items()
        }
        print(f"  best {scoring}={search.best_score_:.4f}: {best_params[model_name]}")

        # Every candidate tried for this model, with its per-fold/mean/std CV scores -
        # kept so the full search (not just the winner) can be inspected later.
        results = pd.DataFrame(search.cv_results_)
        results.columns = [c.removeprefix("param_model__") for c in results.columns]
        safe_name = model_name.replace(" ", "_")
        results.to_csv(os.path.join(output_dir, f"{safe_name}_search_results.csv"), index=False)

    with open(os.path.join(output_dir, "best_params.json"), "w") as f:
        json.dump(best_params, f, indent=2, default=_json_default)

    return best_params


def _json_default(obj):
    # json.dump's fallback for values RandomizedSearchCV hands back that plain json
    # doesn't know how to serialize - mainly numpy scalar types from the scipy.stats
    # distributions (e.g. randint/loguniform sample as np.int64/np.float64). Converts
    # to native int/float rather than stringifying, so the round-tripped value is still
    # directly usable by estimator.set_params() after paramsToModels() reloads it.
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    raise TypeError(f"Object of type {type(obj)} is not JSON serializable")


def paramsToModels(best_params, models=None):
    # Builds a `models` dict ready to pass into evaluate_models from the output of
    # evaluateParameters, by cloning each base estimator and applying its best
    # hyperparameters. Model names in best_params with no matching entry in `models`
    # (e.g. after a round-trip through JSON with a typo) are skipped.
    models = models if models is not None else DEFAULT_MODELS
    return {
        name: clone(models[name]).set_params(**params)
        for name, params in best_params.items()
        if name in models
    }


# ==========================================
# GROUPED CV FOLD SELECTION
# ==========================================
def _choose_group_folds(X, y, bias_values, bias_column, max_splits, random_state):
    # StratifiedGroupKFold can't split a group across folds, so a fixed fold count can
    # land on very unevenly sized folds depending on how the group sizes happen to fall
    # (e.g. one big group forced into a single fold). Searching fold counts from 2 up to
    # max_splits (capped at the number of groups available) and keeping whichever
    # minimizes the gap between the largest and smallest fold gives the split more
    # chances to land on a balanced grouping. Returns the winning candidate's folds
    # directly, so they aren't computed twice.
    n_groups = bias_values.nunique()
    if n_groups < 2:
        raise ValueError(
            f"bias_column {bias_column!r} needs at least 2 unique values to build grouped CV folds, found {n_groups}"
        )
    max_candidate = min(max_splits, n_groups)

    best_n_splits, best_folds, best_gap = None, None, None
    for candidate in range(2, max_candidate + 1):
        sgkf = StratifiedGroupKFold(n_splits=candidate, shuffle=True, random_state=random_state)
        candidate_folds = list(sgkf.split(X, y, groups=bias_values))
        fold_sizes = [len(val_idx) for _, val_idx in candidate_folds]
        gap = max(fold_sizes) - min(fold_sizes)
        if best_gap is None or gap < best_gap or (gap == best_gap and candidate > best_n_splits):
            best_n_splits, best_folds, best_gap = candidate, candidate_folds, gap

    fold_sizes = [len(val_idx) for _, val_idx in best_folds]
    print(
        f"bias_column {bias_column!r}: {n_groups} unique variables, "
        f"chose n_splits={best_n_splits} (searched 2-{max_candidate}), fold sizes: {fold_sizes}"
    )
    return best_folds


# ==========================================
# DATA PREP
# ==========================================
def _build_dataset(ids, feature_columns, response_df, feature_df, id_col, label_col, bias_column=None):
    missing_cols = set(feature_columns) - set(feature_df.columns)
    if missing_cols:
        raise ValueError(f"feature_columns not found in feature_df: {missing_cols}")
    if id_col in feature_columns:
        raise ValueError(f"feature_columns must not include the id column {id_col!r}")
    duplicate_cols = {c for c in feature_columns if list(feature_columns).count(c) > 1}
    if duplicate_cols:
        raise ValueError(f"feature_columns contains duplicate entries: {duplicate_cols}")
    if bias_column is not None and bias_column not in response_df.columns:
        raise ValueError(f"bias_column {bias_column!r} not found in response_df")

    ids = pd.Index(ids).unique()
    resp_cols = [id_col, label_col] if bias_column is None else [id_col, label_col, bias_column]
    resp = response_df[response_df[id_col].isin(ids)][resp_cols]
    feats = feature_df[feature_df[id_col].isin(ids)][[id_col, *feature_columns]]

    merged = resp.merge(feats, on=id_col, how="inner").dropna(subset=[label_col])
    if merged.empty:
        raise ValueError("No samples remain after filtering by ids/matching response and feature rows.")

    X = merged[feature_columns].apply(pd.to_numeric, errors="coerce").reset_index(drop=True)
    y = merged[label_col].astype(int).reset_index(drop=True)
    sample_ids = merged[id_col].reset_index(drop=True)
    bias_values = merged[bias_column].reset_index(drop=True) if bias_column is not None else None
    return X, y, sample_ids, bias_values


def _decision_scores(pipe, X_test):
    model = pipe.named_steps["model"]
    if hasattr(model, "predict_proba"):
        return pipe.predict_proba(X_test)[:, 1]
    if hasattr(model, "decision_function"):
        return pipe.decision_function(X_test)
    return None


def _prediction_rows(ids_test, model_name, fold_idx, y_test, y_pred, y_score, id_col, label_col):
    rows = []
    for i, sample_id in enumerate(ids_test):
        rows.append(
            {
                id_col: sample_id,
                "Model": model_name,
                "Fold": fold_idx,
                label_col: int(y_test.iloc[i]),
                "Predicted": int(y_pred[i]),
                "Probability": float(y_score[i]) if y_score is not None else np.nan,
            }
        )
    return rows

# def identifyMisfolded

# ==========================================
# FEATURE IMPORTANCE
# ==========================================
def _native_importance_rows(pipe, feature_columns, model_name, fold_idx):
    # Only some estimators expose a built-in importance (coef_ for linear
    # models, feature_importances_ for tree ensembles); others (SVC, KNN, MLP)
    # rely solely on permutation importance for a comparable measure.
    model = pipe.named_steps["model"]
    if hasattr(model, "coef_"):
        values = np.abs(np.ravel(model.coef_))
    elif hasattr(model, "feature_importances_"):
        values = np.ravel(model.feature_importances_)
    else:
        return []
    return [
        {
            "Feature": feature,
            "Model": model_name,
            "Fold": fold_idx,
            "Method": "native",
            "Importance": value,
        }
        for feature, value in zip(feature_columns, values)
    ]


def _permutation_importance_rows(
    pipe, X_test, y_test, model_name, fold_idx, n_repeats, scoring, n_jobs, max_samples, random_state
):
    result = permutation_importance(
        pipe,
        X_test,
        y_test,
        n_repeats=n_repeats,
        scoring=scoring,
        n_jobs=n_jobs,
        max_samples=max_samples,
        random_state=random_state,
    )
    return [
        {
            "Feature": feature,
            "Model": model_name,
            "Fold": fold_idx,
            "Method": "permutation",
            "Importance": mean,
        }
        for feature, mean in zip(X_test.columns, result.importances_mean)
    ]


def _summarize_feature_importance(feature_importance):
    if feature_importance.empty:
        return feature_importance
    summary = (
        feature_importance.groupby(["Feature", "Model", "Method"])["Importance"]
        .agg(["mean", "std"])
        .rename(columns={"mean": "Importance_Mean", "std": "Importance_Std"})
        .reset_index()
        .sort_values(["Model", "Method", "Importance_Mean"], ascending=[True, True, False])
    )
    return summary


# ==========================================
# RUN STATE (shared bookkeeping for validation folds and the fixed test set)
# ==========================================
def _init_run_state(models):
    return {
        "fold_rows": [],
        "pred_rows": [],
        "curves": {model_name: {"roc": [], "pr": []} for model_name in models},
        "conf_mats": {model_name: np.zeros((2, 2), dtype=int) for model_name in models},
    }


def _record_predictions(state, model_name, fold_idx, ids_true, y_true, y_pred, y_score, id_col, label_col):
    state["fold_rows"].append(_score_fold(model_name, fold_idx, y_true, y_pred, y_score))
    state["conf_mats"][model_name] += confusion_matrix(y_true, y_pred, labels=[0, 1])
    if y_score is not None:
        fpr, tpr, _ = roc_curve(y_true, y_score)
        state["curves"][model_name]["roc"].append((fpr, tpr, roc_auc_score(y_true, y_score)))
        precision, recall, _ = precision_recall_curve(y_true, y_score)
        state["curves"][model_name]["pr"].append((recall, precision))
    state["pred_rows"].extend(
        _prediction_rows(ids_true, model_name, fold_idx, y_true, y_pred, y_score, id_col, label_col)
    )


def _finalize_run(state, out_dir, prefix, timer):
    with timer.time(f"{prefix}_summarize"):
        fold_results = pd.DataFrame(state["fold_rows"])
        summary = _summarize_folds(fold_results)
        predictions = pd.DataFrame(state["pred_rows"])

    with timer.time(f"{prefix}_write_csvs"):
        fold_results.to_csv(os.path.join(out_dir, f"{prefix}_fold_metrics.csv"), index=False)
        summary.to_csv(os.path.join(out_dir, f"{prefix}_summary_metrics.csv"), index=False)
        predictions.to_csv(os.path.join(out_dir, f"{prefix}_predictions.csv"), index=False)

    with timer.time(f"{prefix}_plot_metric_comparison"):
        _plot_metric_comparison(summary, out_dir, prefix)
    with timer.time(f"{prefix}_plot_roc_curves"):
        _plot_roc_curves(state["curves"], out_dir, prefix)
    with timer.time(f"{prefix}_plot_pr_curves"):
        _plot_pr_curves(state["curves"], out_dir, prefix)
    with timer.time(f"{prefix}_plot_confusion_matrices"):
        _plot_confusion_matrices(state["conf_mats"], out_dir, prefix)

    return fold_results, summary


# ==========================================
# METRICS
# ==========================================
def _score_fold(name, fold_idx, y_test, y_pred, y_score):
    tn, fp, fn, tp = confusion_matrix(y_test, y_pred, labels=[0, 1]).ravel()
    row = {
        "Model": name,
        "Fold": fold_idx,
        "Accuracy": accuracy_score(y_test, y_pred),
        "Precision": precision_score(y_test, y_pred, zero_division=0),
        "Recall": recall_score(y_test, y_pred, zero_division=0),
        "F1": f1_score(y_test, y_pred, zero_division=0),
        "ROC_AUC": np.nan,
        "PR_AUC": np.nan,
        "TP": tp,
        "TN": tn,
        "FP": fp,
        "FN": fn,
        "FPR": fp / (fp + tn) if (fp + tn) > 0 else np.nan,
        "FNR": fn / (fn + tp) if (fn + tp) > 0 else np.nan,
        # PPV/"accuracy of positive examples": TP/(TP+FP) - equivalent to Precision above,
        # kept as its own named column to pair with NPV below.
        "PPV": tp / (tp + fp) if (tp + fp) > 0 else np.nan,
        # NPV/"accuracy of negative examples": TN/(TN+FN)
        "NPV": tn / (tn + fn) if (tn + fn) > 0 else np.nan,
    }
    if y_score is not None:
        row["ROC_AUC"] = roc_auc_score(y_test, y_score)
        row["PR_AUC"] = average_precision_score(y_test, y_score)
    return row


def _summarize_folds(fold_results):
    summary = fold_results.groupby("Model")[METRIC_COLUMNS + COUNT_COLUMNS].agg(["mean", "std"])
    summary.columns = [f"{metric}_{stat}" for metric, stat in summary.columns]
    summary = summary.sort_values("Accuracy_mean", ascending=False).reset_index()
    return summary


# ==========================================
# PLOTS
# ==========================================
def _plot_metric_comparison(summary, output_dir, prefix):
    fig, ax = plt.subplots(figsize=(13, 6))
    n_metrics = len(METRIC_COLUMNS)
    x = np.arange(len(summary))
    width = 0.8 / n_metrics

    for i, metric in enumerate(METRIC_COLUMNS):
        ax.bar(
            x + i * width,
            summary[f"{metric}_mean"],
            width,
            yerr=summary[f"{metric}_std"],
            capsize=3,
            label=metric,
        )

    ax.set_xticks(x + width * (n_metrics - 1) / 2)
    ax.set_xticklabels(summary["Model"], rotation=30, ha="right")
    ax.set_ylabel("Score")
    ax.set_ylim(0, 1.05)
    ax.set_title(f"Model comparison (mean ± std across folds) - {prefix}")
    ax.legend(loc="lower right", ncol=5, fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, f"{prefix}_metric_comparison.png"), dpi=150)
    plt.close(fig)


def _plot_roc_curves(curves, output_dir, prefix):
    fig, ax = plt.subplots(figsize=(7, 7))
    mean_fpr = np.linspace(0, 1, 100)

    for name, data in curves.items():
        if not data["roc"]:
            continue
        tprs = []
        aucs = []
        for fpr, tpr, auc in data["roc"]:
            tprs.append(np.interp(mean_fpr, fpr, tpr))
            tprs[-1][0] = 0.0
            aucs.append(auc)
        mean_tpr = np.mean(tprs, axis=0)
        mean_tpr[-1] = 1.0
        ax.plot(mean_fpr, mean_tpr, label=f"{name} (AUC={np.mean(aucs):.2f})")

    ax.plot([0, 1], [0, 1], linestyle="--", color="grey", label="Chance")
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.set_title(f"Mean ROC curves - {prefix}")
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, f"{prefix}_roc_curves.png"), dpi=150)
    plt.close(fig)


def _plot_pr_curves(curves, output_dir, prefix):
    fig, ax = plt.subplots(figsize=(7, 7))
    mean_recall = np.linspace(0, 1, 100)

    for name, data in curves.items():
        if not data["pr"]:
            continue
        precisions = []
        for recall, precision in data["pr"]:
            order = np.argsort(recall)
            precisions.append(np.interp(mean_recall, recall[order], precision[order]))
        mean_precision = np.mean(precisions, axis=0)
        ax.plot(mean_recall, mean_precision, label=name)

    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.set_title(f"Mean Precision-Recall curves - {prefix}")
    ax.legend(loc="lower left", fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, f"{prefix}_pr_curves.png"), dpi=150)
    plt.close(fig)


def _plot_confusion_matrices(conf_mats, output_dir, prefix):
    names = list(conf_mats.keys())
    fig, axes = plt.subplots(1, len(names), figsize=(4 * len(names), 4))
    if len(names) == 1:
        axes = [axes]

    for ax, name in zip(axes, names):
        cm = conf_mats[name]
        im = ax.imshow(cm, cmap="Blues")
        ax.set_title(name, fontsize=9)
        ax.set_xticks([0, 1])
        ax.set_yticks([0, 1])
        ax.set_xticklabels(["0", "1"])
        ax.set_yticklabels(["0", "1"])
        ax.set_xlabel("Predicted")
        ax.set_ylabel("True")
        for i in range(2):
            for j in range(2):
                ax.text(j, i, cm[i, j], ha="center", va="center", color="black")

    fig.suptitle(f"Aggregated confusion matrices (all folds) - {prefix}")
    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, f"{prefix}_confusion_matrices.png"), dpi=150)
    plt.close(fig)

############################################
#############Run function###################
############################################

def runTraining(splits, features, response_df, feature_df, models=DEFAULT_MODELS, skip_existing=True):
    for split_name, split_data in splits.items():
        train_ids = split_data["train"]
        test_ids = split_data["test"]
        output_dir = split_data["Dir"]
        os.makedirs(output_dir, exist_ok=True)
        for feature_name, feature_columns in features.items():
            dir = os.path.join(output_dir, split_name, feature_name)
            os.makedirs(dir, exist_ok=True)
            # Call the training function here with the appropriate parameters
            test_name = split_name + "-" + feature_name
            metrics_file = os.path.join(dir, "test_summary_metrics.csv")
            if skip_existing and os.path.exists(metrics_file):
                print(f"Skipping split: {split_name}, feature set: {feature_name} - {metrics_file} already exists")
                reportTestMetrics(metrics_file)
                print("#############################################################")
                continue
            print(f"Running evaluation for split: {split_name}, feature set: {feature_name}")
            evaluate_models(
                    test_name,
                    train_ids,
                    feature_columns,
                    response_df,
                    feature_df,
                    test_ids=test_ids,
                    output_dir=dir,
                    importance_method="native",
                    timing_detail = "short",
                    models=models
                )
            #Print the test performance metrics summary
            print(f"Completed evaluation for split: {split_name}, feature set: {feature_name}")
            reportTestMetrics(os.path.join(dir, f"test_summary_metrics.csv"))
            print("#############################################################")

def runFixedSplit(splits, features, response_df, feature_df, models=DEFAULT_MODELS, skip_existing=True):
    # Like runTraining, but for a pre-defined, non-CV train/test split (e.g. the DL
    # pipeline's fixed FullRun split) - one fit per model, scored once against the fixed
    # test set, no inner k-fold CV. Keeps the same <Dir>/<split_name>/<feature_name>/
    # directory layout runTraining uses, so createSummary reads both identically.
    for split_name, split_data in splits.items():
        train_ids = split_data["train"]
        test_ids = split_data["test"]
        output_dir = split_data["Dir"]
        os.makedirs(output_dir, exist_ok=True)
        for feature_name, feature_columns in features.items():
            dir = os.path.join(output_dir, split_name, feature_name)
            os.makedirs(dir, exist_ok=True)
            run_name = split_name + "-" + feature_name
            metrics_file = os.path.join(dir, "test_summary_metrics.csv")
            if skip_existing and os.path.exists(metrics_file):
                print(f"Skipping split: {split_name}, feature set: {feature_name} - {metrics_file} already exists")
                reportTestMetrics(metrics_file)
                print("#############################################################")
                continue
            print(f"Running fixed-split evaluation for split: {split_name}, feature set: {feature_name}")
            evaluate_fixed_split(
                run_name,
                train_ids,
                test_ids,
                feature_columns,
                response_df,
                feature_df,
                output_dir=dir,
                importance_method="native",
                models=models,
            )
            print(f"Completed fixed-split evaluation for split: {split_name}, feature set: {feature_name}")
            reportTestMetrics(os.path.join(dir, "test_summary_metrics.csv"))
            print("#############################################################")


def reportTestMetrics(file):
    if os.path.exists(file):
        df = pd.read_csv(file)
        print("Accuracy - " + ", ".join(f"{col1}: {col2:.3f} ± {col3:.3f}" for col1, col2, col3 in zip(df['Model'], df['Accuracy_mean'], df['Accuracy_std'])))
        print("ROC_AUC - " + ", ".join(f"{col1}: {col2:.3f} ± {col3:.3f}" for col1, col2, col3 in zip(df['Model'], df['ROC_AUC_mean'], df['ROC_AUC_std'])))
    else:
        print(f"Test metrics file not found: {file}")


def createSummary(summary_dirs, output_file):
    """
    Compiles the results scattered across one or more runTraining output trees into a
    single dataframe, and saves it to output_file.

    Each directory in summary_dirs is expected to look like:
        <summary_dir>/<split>/<feature_set>/test_summary_metrics.csv
    (i.e. what runTraining's `dir = os.path.join(output_dir, split_name, feature_name)`
    writes). Every test_summary_metrics.csv found is read in and tagged with which
    Dir/Split/FeatureSet it came from before being concatenated together, so each row of
    the result is one (dir, split, feature set, model) combination.
    """
    rows = []
    for summary_dir in summary_dirs:
        for split_name in sorted(os.listdir(summary_dir)):
            split_path = os.path.join(summary_dir, split_name)
            if not os.path.isdir(split_path):
                continue
            for feature_name in sorted(os.listdir(split_path)):
                feature_path = os.path.join(split_path, feature_name)
                metrics_file = os.path.join(feature_path, "test_summary_metrics.csv")
                if not os.path.exists(metrics_file):
                    continue
                df = pd.read_csv(metrics_file)
                df.insert(0, "Dir", summary_dir)
                df.insert(1, "Split", split_name)
                df.insert(2, "FeatureSet", feature_name)
                rows.append(df)

    if not rows:
        print(f"No test_summary_metrics.csv files found under {summary_dirs}")
        return pd.DataFrame()

    combined = pd.concat(rows, ignore_index=True)
    combined.to_csv(output_file, index=False)
    return combined


def getMisfolded(json_dir, threshold=0.5):
    dirs = os.listdir(json_dir)
    json_scores = []
    for d in dirs:
        json_file = os.path.join(json_dir, d, f"{d}.json")
        if not os.path.exists(json_file):
            continue
        with open(json_file, "r") as f:
            data = json.load(f)
        score = {"name": d}
        score.update(data.get("Scoring", {}))
        json_scores.append(score)
    json_scores_df = pd.DataFrame(json_scores)
    high_score_names = json_scores_df.loc[json_scores_df["Final"] < threshold]["name"].tolist()
    return high_score_names

############################################################
################# Main #####################################
############################################################

if __name__ == "__main__":

    ##########################################
    #############Setup########################
    #########################################

    # cls_data: sample id -> response label (and metadata used to define subsets)
    cls_data = pd.read_csv("./data/ML/Annotation.csv")
    # feature_matrix: sample id -> structural/sequence features
    feature_matrix = pd.read_csv("./data/Features/features.csv")
    #First, subset the feature matrix by the misfolded samples in the ./dat directory
    misfolded_samples = getMisfolded("./dat", threshold=0.8) #Using a 0.8 threshold should get rid of 60 poor folded or non-useful samples
    feature_matrix = feature_matrix[~feature_matrix["name"].isin(misfolded_samples)]
    #Subset the cls_data by the name column in the feature_matrix
    cls_data = cls_data[cls_data["name"].isin(feature_matrix["name"])]
    all_ids = cls_data["name"]
    feature_columns = [c for c in feature_matrix.columns if c != "name"]
    base_columns = feature_columns[0:1]
    sequence_columns = feature_columns[1:164] #columns 2-165
    structural_columns = feature_columns[164:364] #columns 166-365
    physics_columns = feature_columns[364:] #columns 366-end
    features = {'All': feature_columns,
                'Seq': base_columns + sequence_columns,
                'Str': base_columns + structural_columns, 
                'Phys': base_columns + physics_columns} #list of lists with the different feature subsets to run as models
    print("Found ", len(all_ids), "samples and ", len(feature_columns), "features in dataset.")
    #Break it down by source
    for source, source_df in cls_data.groupby("source"):
        print(f"Source: {source}, Number of samples: {len(source_df)}")


    ############################################
    ##########Set the model runs################
    ############################################
    #Save the test set which is all features and every sample not from the "immrep" source
    test = cls_data[cls_data["source"] != "immrep"]["name"]
    cls_scan = cls_data.copy()
    cls_data = cls_data[cls_data["source"] != "scan"]
    cls_core = cls_data[cls_data["source"].isin(["exp","fake","ground"])]
    # For each set of features ['all','sequence','structural','physics'] we want to create one model trained and tested with subset
    # We will want to make the following splits of the train and test data
    # 1. All tested on all (CV), non-immpred tested on immpred, immpred trained on immpred
    basic = {
        'test': {
                    'train': cls_data[cls_data["source"].isin(["exp","fake","ground"])]["name"],
                    'test': None,
                    'Dir': "./data/ML/Results/Base/"
                },
        # 'exp_to_ground':  {
        #     'train': cls_data[cls_data["source"].isin(["exp","fake"])]["name"],
        #     'test': cls_data[cls_data["source"].isin(["ground"])]["name"],
        #     'Dir': "./data/ML/Results/Base/"
        # },
        'all_to_all': {
            'train':  cls_data[cls_data["source"] != "scan"]["name"],
            'test': None,
            'Dir': "./data/ML/Results/Base/"
        },
        'non_immpred_to_immpred': {
            'train': cls_data[~cls_data["source"].isin(["immrep","scan"])]["name"],
            'test': cls_data[cls_data["source"] == "immrep"]["name"],
            'Dir': "./data/ML/Results/Base/"
        },
        'immpred_to_immpred': {
            'train': cls_data[cls_data["source"] == "immrep"]["name"],
            'test': cls_data[cls_data["source"] == "immrep"]["name"],
            'Dir': "./data/ML/Results/Base/"
        },
        'all_to_scan': {
            'train': cls_scan[cls_scan["source"] != "scan"]["name"],
            'test': cls_scan[cls_scan["source"] == "scan"]["name"],
            'Dir': "./data/ML/Results/Base/"
        }
    }
    # 2. Each batch (source) train and test on rest
    batch = {
        'Exp_Fake': {
            'train': cls_data[cls_data["source"].isin(["exp","fake"])]["name"],
            'test': cls_data[cls_data["source"].isin(["exp","fake"])]["name"],
            'Dir': "./data/ML/Results/Batch/"
        },
        'Exp_Fake_to_ground': {
            'train': cls_data[~cls_data["source"].isin(["exp","fake"])]["name"],
            'test': cls_data[cls_data["source"].isin(["ground"])]["name"],
            'Dir': "./data/ML/Results/Batch/"
        },
        'all_to_fake': {
                    'train': cls_data[~cls_data["source"].isin(["fake"])]["name"],
                    'test': cls_data[cls_data["source"].isin(["fake"])]["name"],
                    'Dir': "./data/ML/Results/Batch/"
                }
    }
    # 3. MHC class 1 test on class 2 (vice versa), each MHC class (> 100 rows) tested on all other MHC classes, each peptide (> 100 rows) tested on all other peptides
    
    #mhc type greater than 100
    mhc_types = cls_data["Type"].value_counts()[cls_data["Type"].value_counts() > 100].index
    print("MHC types with more than 100 entries:", len(mhc_types))
    #peptide greater than 100
    peptides = cls_data["Peptide"].value_counts()[cls_data["Peptide"].value_counts() > 100].index
    print("Peptides with more than 100 entries:", len(peptides))


    Bias = {
        'MHC_class_1_to_class_1': {
            'train': feature_matrix[feature_matrix["MHC_class"] == 1]["name"],
            'test': feature_matrix[feature_matrix["MHC_class"] == 1]["name"],
            'Dir': "./data/ML/Results/Bias/"
        },
        'MHC_class_1_to_class_2': {
            'train': feature_matrix[feature_matrix["MHC_class"] == 1]["name"],
            'test': feature_matrix[feature_matrix["MHC_class"] == 2]["name"],
            'Dir': "./data/ML/Results/Bias/"
        },
        'MHC_class_2_to_class_1': {
            'train': feature_matrix[feature_matrix["MHC_class"] == 2]["name"],
            'test': feature_matrix[feature_matrix["MHC_class"] == 1]["name"],
            'Dir': "./data/ML/Results/Bias/"
        },
        'MHC_class_2_to_class_2': {
            'train': feature_matrix[feature_matrix["MHC_class"] == 2]["name"],
            'test': feature_matrix[feature_matrix["MHC_class"] == 2]["name"],
            'Dir': "./data/ML/Results/Bias/"
        }
    }
    Bias.update(
        {
            f'MHC_class_{i}': {
                'train': cls_core[cls_core["Type"] == i]["name"],
                'test': cls_core[cls_core["Type"] != i]["name"],
                'Dir': "./data/ML/Results/Bias/"
            }
            for i in mhc_types
        }
    )
    # Bias.update(
    #     {
    #         f'peptide_{i}': {
    #             'train': cls_core[cls_core["Peptide"] == i]["name"],
    #             'test': cls_core[cls_core["Peptide"] != i]["name"],
    #             'Dir': "./data/ML/Results/Bias/"
    #         }
    #         for i in peptides
    #     }
    # )

    ############################################
    ##############Run Models####################
    ############################################

    # Run once on every sample shared between the two tables

    param_dir = "./data/ML/Params"
    param_file = os.path.join(param_dir, "best_params.json")
    if not os.path.exists(param_file):
        print("No existing parameter file found. Evaluating parameters...")
        best_params = evaluateParameters(
            test,
            feature_columns,
            cls_data,
            feature_matrix,
            param_dir,
            n_splits=5,
            n_iter=25,
        )
    else:
        print("Found existing parameter file.")
        with open(param_file, "r") as f:
            best_params = json.load(f)
    models = paramsToModels(best_params)

    #############################################
    ########Test#################################
    #############################################
    # print(test)
    # evaluate_models(
    #     "Test",
    #     test,
    #     feature_columns,
    #     cls_data,
    #     feature_matrix,
    #     models = models,
    #     output_dir="./data/ML/Results/Test",
    #     importance_method="both",
    #     timing_detail = "short"
    # )
    #############################################
    #############Calculate Effects################
    #############################################
    # evaluate_models(
    #         "Test",
    #         test,
    #         feature_columns,
    #         cls_data,
    #         feature_matrix,
    #         models = models,
    #         bias_column="Type",
    #         output_dir="./data/ML/Results/Effects/MHCtype",
    #         importance_method="both",
    #         timing_detail = "short"
    #     )
    # evaluate_models(
    #             "Test",
    #             test,
    #             feature_columns,
    #             cls_data,
    #             feature_matrix,
    #             models = models,
    #             bias_column="Peptide",
    #             output_dir="./data/ML/Results/Effects/Peptide",
    #             importance_method="both",
    #             timing_detail = "short"
    #         )
    # evaluate_models(
            #     "Test",
            #     test,
            #     feature_columns,
            #     cls_data,
            #     feature_matrix,
            #     models = models,
            #     bias_column="TRBV",
            #     output_dir="./data/ML/Results/Effects/TRBV",
            #     importance_method="both",
            #     timing_detail = "short"
            # )
    # 4. Fixed external split - the exact train/val assignment the DL pipeline used for
    # its FullRun comparison (./data/DL/FullRun/split.csv), so the ML and NetTCR results
    # are scored against the identical holdout set. No CV here: one fit per model,
    # scored once against the fixed val set (see runFixedSplit/evaluate_fixed_split).
    fullrun_split = pd.read_csv("./data/DL/FullRun/split.csv")
    fullrun_split = fullrun_split[fullrun_split["name"].isin(cls_data["name"])]
    FullRun = {
        'FullRun': {
            'train': fullrun_split[fullrun_split["Set"] == "train"]["name"],
            'test': fullrun_split[fullrun_split["Set"] == "val"]["name"],
            'Dir': "./data/ML/Results/FullRun/"
        }
    }

    ################################################
    ############Calculate all#######################
    ################################################
    runTraining(basic, features, cls_scan, feature_matrix, models=models)
    runTraining(batch, features, cls_scan, feature_matrix, models=models)
    runTraining(Bias, features, cls_scan, feature_matrix, models=models)
    runFixedSplit(FullRun, features, cls_scan, feature_matrix, models=models)
    createSummary(["./data/ML/Results/Base", "./data/ML/Results/Bias", "./data/ML/Results/Batch", "./data/ML/Results/FullRun"], output_file="./data/ML/Results/results.csv")

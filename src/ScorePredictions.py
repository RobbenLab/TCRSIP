#!/usr/bin/env python3

# Shared metric computation for any model that outputs raw per-sample binder
# probabilities rather than already-aggregated metrics (NetTCR-2.0, ERGO-II, and our own
# ML/DL prediction files all do this) - one place to compute the same Accuracy/ROC_AUC/
# PR_AUC/PPV/NPV set everywhere instead of reimplementing it per source. Attaching the
# true label to a particular raw-output file (e.g. ERGO's, which doesn't carry one at
# all) is still source-specific and stays in that source's own script.

import pandas as pd
from sklearn.metrics import accuracy_score, average_precision_score, confusion_matrix, roc_auc_score


def compute_binary_metrics(y_true, y_score, threshold=0.5):
    """
    Accuracy/ROC_AUC/PR_AUC/PPV/NPV (plus n/n_positive) for a set of raw binder
    probabilities scored against known binary labels.

    y_true: array-like of 0/1 binder labels.
    y_score: array-like of predicted binder probabilities, same order/length as y_true.
    threshold: probability cutoff used to derive the predicted class (default 0.5).
    """
    y_true = pd.Series(y_true).reset_index(drop=True).astype(int)
    y_score = pd.Series(y_score).reset_index(drop=True).astype(float)
    y_pred = (y_score >= threshold).astype(int)

    tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=[0, 1]).ravel()
    return {
        "n": len(y_true),
        "n_positive": int(y_true.sum()),
        "Accuracy": accuracy_score(y_true, y_pred),
        "ROC_AUC": roc_auc_score(y_true, y_score),
        "PR_AUC": average_precision_score(y_true, y_score),
        # PPV ("positive set accuracy"): of samples predicted positive, how many truly are.
        "PPV": tp / (tp + fp) if (tp + fp) > 0 else float("nan"),
        # NPV ("negative set accuracy"): of samples predicted negative, how many truly are.
        "NPV": tn / (tn + fn) if (tn + fn) > 0 else float("nan"),
    }


def score_predictions_file(predictions_csv, label_col, score_col, threshold=0.5, filter_fn=None):
    """
    Convenience wrapper for a predictions CSV that already carries both the true label
    and the predicted score as columns (e.g. NetTCR's own --outfile, or our DL
    pipeline's predictions.csv) - unlike ERGO's raw results, which need the label
    attached separately first since they don't include one at all (see ScoreERGO.py).

    filter_fn, if given, is applied to the loaded dataframe before scoring (e.g. to
    select one Set/Fold out of a file covering more than one).
    """
    df = pd.read_csv(predictions_csv)
    if filter_fn is not None:
        df = filter_fn(df)
    return compute_binary_metrics(df[label_col], df[score_col], threshold=threshold)

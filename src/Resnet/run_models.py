#!/usr/bin/env python3

# Model runs for the TCR ResNet. Each entry in RUNS is one model (or one set of CV folds):
# which samples it uses, which features, and whether to cross-validate. Everything else -
# datasets, loaders, normalization, training, checkpoints, predictions, TensorBoard - is
# built from that entry at runtime, so a run is one line:
#
#   run_model("curated_all")      # one run
#   run_all()                     # every run in RUNS
#
# Samples come from Annotation.csv (name, response, source, Allele, species, Peptide, TRAV,
# TRBV, Type); their input maps are looked up by name in Input_structures.csv (written by
# build_inputs.py). Features are the always-included position channels plus any of the
# sequence, structure and physics groups. Standalone apart from tcr_dataset.py and
# tcr_model.py, so these three files plus data/input/ and Annotation.csv run in Colab:
#
#   from run_models import CONFIG, run_all
#   CONFIG.update(inputs="/content/drive/MyDrive/TCR_DL/input/Input_structures.csv",
#                 labels="/content/drive/MyDrive/TCR_DL/Annotation.csv",
#                 out_dir="/content/drive/MyDrive/TCR_DL/results")   # results + tensorboard
#   run_all()
#
# Output per model under results/<run>/<fold0.. or single>/:
#   split.csv, checkpoint.pt, best.pt, history.csv, val_predictions.csv,
#   test_predictions.csv, metrics.json, train.log, and model_info.json +
#   channel_stats.npz (everything needed to apply the model to new structures).
# Plus results/results.csv (one row per model) and results/results_summary.csv (mean/std
# over folds). TensorBoard logs go to results/tensorboard/<run>/<fold>:
#   tensorboard --logdir results/tensorboard        (Colab: %tensorboard --logdir ...)
# A finished model is skipped and an interrupted one resumes from its checkpoint, so the
# script can just be re-run after a disconnect.

import argparse
import json
import os
import random
import shutil
import time

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.optim as optim
from sklearn.model_selection import StratifiedGroupKFold, StratifiedKFold, train_test_split
from torch.utils.tensorboard import SummaryWriter

from tcr_dataset import (FEATURE_GROUPS, TCRMapDataset, channel_stats, check_canvas_fit,
                         load_inputs, make_loader, preload_maps)
from tcr_model import build_resnet, compute_metrics, predict, train_model

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

##############################################
###       Settings                         ###
##############################################
# Model/optimizer defaults are the notebook's: ResNet [3, 4, 6, 3], Adam, lr 1e-5, weight
# decay 0.1, batch 16, 100 epochs with early stopping after 10% of them without an
# improvement in validation loss. Any of these can be overridden per run.
CONFIG = {
    "inputs": os.path.join(PROJECT_DIR, "data", "input", "Input_structures.csv"),
    "labels": os.path.join(PROJECT_DIR, "Annotation.csv"),
    "out_dir": os.path.join(PROJECT_DIR, "results"),  # <out_dir>/<run>/<fold>/ and <out_dir>/tensorboard/
    "name_col": "name",          # sample name column in Annotation.csv
    "label_col": "response",     # label column in Annotation.csv
    "exclude_flags": ["short_tcr_a", "short_tcr_b"],  # drop fragment TCR chains (Flags column)
    "epochs": 100,
    "batch_size": 8,
    "lr": 1e-5,
    "weight_decay": 0.005,
    "layers": [3, 4, 6],      # blocks per residual stage; the length is the depth
    "base_filters": 64,          # channels in the first stage, doubling each stage after
    "max_filters": None,         # ceiling on that doubling (None = uncapped)
    "pool_size": 4,              # adaptive pool before the classifier (None = notebook's AvgPool2d(7))
    "hidden": 1024,              # width of the fully connected layer
    "dropout": 0.2,              # dropout on the hidden layer (0 = the notebook's behaviour)
    "conv_dropout": 0,         # Dropout2d on the conv stack before the flatten
    "track_test": True,          # also score the test set every epoch (curve only, never
                                 # used for selection); costs one extra pass per epoch
    "patience": 10,              # stop after this many epochs with no better validation loss
    "min_delta": 0.0,            # improvement smaller than this does not count as better
    "checkpoint_every": 5,
    "amp": False,                # mixed precision; not needed at 512 (fp32 fits in ~6 GB)
    "test_size": 0.2,            # held-out test fraction when a run defines no test set
    "val_size": 0.15,            # validation fraction (early stopping) when not cross-validating
    "val_folds": 6,              # grouped runs without cv: upper bound on the fold search
    "normalize": True,           # per-channel z-score, stats from the training samples
    "threshold": 0.5,            # probability above which a sample is called positive
    "pretrain_epochs": 0,        # >0: warm up on a class-balanced split before the real one
    "pretrain_balance": "undersample",   # or "oversample"
    "save_train_predictions": True,  # also score the training split with the best weights
    "preload": False,            # cache maps in RAM (~17 MB each, ~41 GB for all)
    "workers": min(8, os.cpu_count() or 1),
    "device": "cuda:0" if torch.cuda.is_available() else "cpu",
    "seed": 42,
    "params_file": os.path.join(PROJECT_DIR, "results", "best_params.json"),  # written by hypersearch.py
    "use_params": True,          # use the tuned parameters from params_file when it exists
    "params_key": "default",     # which entry of params_file to use
    "tb_dir": None,              # TensorBoard root; None = <out_dir>/tensorboard
    "tb_tag": None,              # tag prefix = the collapsible section; None = the run name
}

##############################################
###       Runs                             ###
##############################################
# Run name -> settings. All optional:
#   select    {column: value or [values]} picking the samples from the annotation columns
#             (name, response, source, Allele, species, Peptide, TRAV, TRBV, Type,
#             MHC_class, CDR3a, CDR3b, Len_*, Flags); default every sample
#   query     extra pandas query on the same table, e.g. "MHC_class == 1 and Len_Peptide <= 10"
#   limit     keep only this many samples, drawn label-stratified (for quick test runs)
#   test      {column: ...} and/or "test_query": held-out test samples, removed from the
#             training pool. Without it, a stratified test_size fraction is split off;
#             with "test_size": 0 there is no test set (train on everything).
#   features  any of "sequence", "structure", "physics" (position channels are always
#             included); default all three
#   cv        number of stratified CV folds -> that many models, each validated on its own
#             fold and scored on the test set (test performance). 0/None (default) trains
#             one model with a val_size validation split (inference model).
#   use_params  False ignores the tuned parameters from hypersearch.py for this run; a
#               string picks a different entry of params_file than "default". Tuned
#               parameters are used automatically, and anything set explicitly in a run
#               still wins over them.
#   rerun       True deletes this run's previous output and trains it again (otherwise a
#               finished model is skipped and an interrupted one resumes from its checkpoint)
#   importance  True runs the gradient feature importance (importance.py) after training,
#               or a dict of its options, e.g. {"n": 16, "split": "val", "method": "grad"}
#   pretrain_epochs  warm up for this many epochs on a class-balanced version of the fold's
#               own train/val sets, then fine-tune on the real split. Useful where the
#               classes are skewed (immrep is 17% positive), since a model trained straight
#               on a skewed split can reach a low loss by predicting the majority class.
#   pretrain_balance  "undersample" (default, drops majority samples) or "oversample"
#   tb_group    folder the run's TensorBoard logs go in, so the view can be scoped to a
#               group of runs (tensorboard --logdir results/tensorboard/<tb_group>)
#   group       an annotation column (Peptide, TRAV, TRBV, Allele, Type, ...). Every sample
#               sharing a value stays in one fold and the held-out test set takes whole
#               groups, so the model is scored on values it never trained on. Fold sizes
#               come out uneven - the fold count is chosen to even them out as far as the
#               group sizes allow.
#   plus any CONFIG key to override for this run, e.g. "epochs": 50, "lr": 1e-4
#
# Examples:
#   RUNS["curated_cv"] = {"select": {"source": ["exp", "fake", "ground"]}, "cv": 5}
#   RUNS["curated_to_immrep_structure"] = {
#       "select": {"source": ["exp", "fake", "ground"]},
#       "test": {"source": "immrep"},
#       "features": ["structure"],
#   }
#   RUNS["inference_all"] = {"features": ["sequence", "structure", "physics"], "test_size": 0}
#   RUNS["quick"] = {"limit": 200, "epochs": 4, "rerun": True, "importance": {"n": 4}}
# Shared settings for the model runs below
CORE = ["exp", "fake", "ground"]   # the curated dataset: experimental, decoys, crystal
CV = 5                             # folds per run; raise for tighter error bars, at 1 model per fold

RUNS = {
    # Quick pipeline test: a small stratified subset, every feature, few epochs. Checks that
    # the whole path runs end to end (~2 min) before starting real runs.
    # "test_all": {"tb_group": "test", "select": {"source": CORE}, "features": FEATURE_GROUPS,
    #              "limit": 200, "epochs": 100, "cv": 3, "rerun": True, "importance": {"n": 4}},

    ##########################################################################
    # 1. Baseline: the core dataset with every feature group
    ##########################################################################
    "core_all": {"tb_group": "core", "select": {"source": CORE}, "features": FEATURE_GROUPS,
                 "cv": CV, "importance": True},

    ##########################################################################
    # 2. Feature ablation: the same samples and the same folds as core_all, so the
    #    differences are the features and nothing else. Position channels (chain, region,
    #    sequence separation) are always present, so each run is position + one group.
    ##########################################################################
    "core_seq": {"tb_group": "core", "select": {"source": CORE}, "features": ["sequence"],
                 "cv": CV, "importance": True},
    "core_str": {"tb_group": "core", "select": {"source": CORE}, "features": ["structure"],
                 "cv": CV, "importance": True},
    "core_phys": {"tb_group": "core", "select": {"source": CORE}, "features": ["physics"],
                  "cv": CV, "importance": True},
    # Control: position channels only. Nothing about the residues or their contacts, so a
    # good score here means the model is reading chain lengths or padding geometry rather
    # than biology - the number the three ablations above have to beat. Delete if not wanted.
    "core_pos": {"tb_group": "core", "select": {"source": CORE}, "features": [],
                 "cv": CV, "importance": True},

    ##########################################################################
    # 3. Generalization: train on the predicted structures (exp) and decoys (fake), test on
    #    the crystal structures (ground), which the model never sees in training.
    #    Note: ground is 100% positive, so the test set has one class - ROC AUC and PR AUC
    #    come out NaN and the meaningful numbers are accuracy and recall (the share of real
    #    complexes called positive). Add negatives to the test set if you want ranking
    #    metrics here.
    ##########################################################################
    "expfake_to_ground": {"tb_group": "transfer", "select": {"source": ["exp", "fake"]}, "test": {"source": "ground"},
                          "features": FEATURE_GROUPS, "cv": CV, "importance": True,"pretrain_epochs": 20, "pretrain_balance": "oversample"},

    ##########################################################################
    # 4. MHC class transfer: train on one class, test on the other. Both directions use the
    #    core sources, since every class II sample comes from them (123 of them: 69 negative,
    #    54 positive, against 1050 class I).
    #    Note the asymmetry: class II puts its beta chain on MHCb, which is empty for every
    #    class I complex, so a class I-trained model has never seen that chain occupied -
    #    and the class II direction trains on ~80 samples per fold, so expect wide folds.
    ##########################################################################
    "mhc1_to_mhc2": {"tb_group": "transfer", "select": {"source": CORE, "MHC_class": 1},
                     "test": {"source": CORE, "MHC_class": 2},
                     "features": FEATURE_GROUPS, "cv": CV, "importance": True,
                     "pretrain_epochs": 20, "pretrain_balance": "oversample"},
    "mhc2_to_mhc1": {"tb_group": "transfer", "select": {"source": CORE, "MHC_class": 2},
                     "test": {"source": CORE, "MHC_class": 1},
                     "features": FEATURE_GROUPS, "cv": CV, "importance": True,
                     "pretrain_epochs": 20, "pretrain_balance": "oversample"},

    ##########################################################################
    # 5. The immrep set (1022 samples, 17% positive against core's 55%). Two questions:
    #    does a core-trained model transfer to it, and does training on both do better.
    ##########################################################################
    # Core only, immrep held out entirely - the honest generalization test
    "core_to_immrep": {"tb_group": "transfer", "select": {"source": CORE}, "test": {"source": "immrep"},
                       "features": FEATURE_GROUPS, "cv": CV, "importance": True},
    # Both pooled, tested on a held-out fraction of the pool
    "core_immrep": {"tb_group": "transfer", "select": {"source": CORE + ["immrep"]}, "features": FEATURE_GROUPS,
                    "cv": CV, "importance": True,"pretrain_epochs": 20, "pretrain_balance": "oversample"},

    ##########################################################################
    # 6. The scan set: 26 point mutants of a few complexes (62% positive), held out from a
    #    model trained on everything else. A small, hard test - a handful of mutations away
    #    from structures the model has seen, so it asks whether the model reads the
    #    substitution rather than the overall complex.
    ##########################################################################
    "all_to_scan": {"tb_group": "transfer", "select": {"source": CORE + ["immrep"]}, "test": {"source": "scan"},
                    "features": FEATURE_GROUPS, "cv": CV, "importance": True, "pretrain_epochs": 20, "pretrain_balance": "oversample"},

    ##########################################################################
    # 7. The model to keep: every sample, every feature, one model, nothing held out for
    #    testing (test_size 0). A validation split is still taken for early stopping, and
    #    the importance run reports on it since there is no test set. Use this one for
    #    inference on new structures; the runs above are what the numbers come from.
    ##########################################################################
    "all_final": {"tb_group": "final", "select": {"source": CORE + ["immrep"]}, "features": FEATURE_GROUPS, "test_size": 0, 
                  "importance": True,"pretrain_epochs": 20, "pretrain_balance": "oversample"},

    ##########################################################################
    # 8. Bias: one run per factor, with every sample sharing a value of that factor kept in
    #    a single fold and the test set made of whole groups. So the model is scored on
    #    peptides / V genes / alleles it has never trained on, and the drop against
    #    core_all (same samples, random folds) is how much of its performance came from
    #    recognising the background rather than the interaction.
    #    Folds come out uneven by construction - the fold count is chosen to even them out
    #    as far as the group sizes allow.
    ##########################################################################
    "bias_peptide": {"tb_group": "bias", "select": {"source": CORE}, "group": "Peptide",
                     "features": FEATURE_GROUPS, "cv": CV, "importance": True},
    "bias_trbv": {"tb_group": "bias", "select": {"source": CORE}, "group": "TRBV",
                  "features": FEATURE_GROUPS, "cv": CV, "importance": True},
    "bias_trav": {"tb_group": "bias", "select": {"source": CORE}, "group": "TRAV",
                  "features": FEATURE_GROUPS, "cv": CV, "importance": True},
    "bias_allele": {"tb_group": "bias", "select": {"source": CORE}, "group": "Type",
                    "features": FEATURE_GROUPS, "cv": CV, "importance": True},
}


##############################################
###       Samples                          ###
##############################################

def load_samples(cfg=None):
    """Every sample with both a label in Annotation.csv and an input map in
    Input_structures.csv: the Annotation.csv columns plus the input columns (File, Path,
    MHC_class, CDR3a, CDR3b, Len_*, Anchor_index, TCR_side, Flags, ...)."""
    cfg = {**CONFIG, **(cfg or {})}
    name = cfg["name_col"]
    labels = pd.read_csv(cfg["labels"], index_col=0)  # first column is an R row index
    inputs = load_inputs(cfg["inputs"]).rename(columns={"Name": name})
    overlap = [c for c in inputs.columns if c in labels.columns and c != name]
    samples = labels.merge(inputs.drop(columns=overlap), on=name, how="inner")
    print(f"{len(samples)} samples with a label and an input map ({len(labels) - len(samples)} "
          f"labelled without a map, {len(inputs) - len(samples)} maps without a label)")
    if cfg["exclude_flags"]:
        flags = samples["Flags"].fillna("").str.split(";")
        drop = flags.apply(lambda f: any(x in cfg["exclude_flags"] for x in f))
        if drop.any():
            print(f"Excluding {drop.sum()} samples flagged {cfg['exclude_flags']}")
        samples = samples[~drop].reset_index(drop=True)
    check_canvas_fit(samples.rename(columns={name: "Name"}))
    return samples


def select_samples(samples, select=None, query=None):
    """Subset samples by annotation columns ({column: value or [values]}) and/or a pandas
    query string."""
    keep = pd.Series(True, index=samples.index)
    for col, values in (select or {}).items():
        if col not in samples:
            raise KeyError(f"No column {col!r}; available: {list(samples.columns)}")
        values = list(values) if isinstance(values, (list, tuple, set, pd.Series)) else [values]
        keep &= samples[col].isin(values)
    selected = samples[keep]
    if query:
        selected = selected.query(query)
    return selected.reset_index(drop=True)


def describe(samples, label_col, prefix=""):
    by_source = (", " + ", ".join(f"{s} {int(r['size'])}" for s, r in
                                  samples.groupby("source")[label_col].agg(["size", "mean"]).iterrows())
                 if "source" in samples else "")
    print(f"{prefix}{len(samples)} samples, {samples[label_col].mean():.0%} positive{by_source}")


##############################################
###       Splits                           ###
##############################################

def group_values(samples, column):
    """Group labels for grouped CV. A missing value is not evidence that two samples belong
    together, so each one becomes its own group and they spread across the folds."""
    values = samples[column]
    missing = values.isna()
    out = values.astype("object").copy()
    out[missing] = [f"__missing_{i}" for i in range(int(missing.sum()))]
    return out.to_numpy()


def group_folds(y, groups, n_splits, seed, log=print):
    """Stratified folds that never split a group across folds.

    The requested fold count is used when it works. Because a group cannot be divided, a
    few large groups can force a fold to come out empty or single-class - then the fold
    count steps down until every fold is usable. Of the random seeds tried at the winning
    count, the one whose folds are closest to equal size is kept, which is the most even
    spread the group sizes allow."""
    n_groups = len(set(groups))
    if n_groups < 2:
        raise ValueError(f"grouped CV needs at least 2 groups, found {n_groups}")
    for n in range(min(int(n_splits), n_groups), 1, -1):
        best = None
        for attempt in range(8):      # same fold count, different shuffles
            folds = list(StratifiedGroupKFold(n_splits=n, shuffle=True, random_state=seed + attempt)
                         .split(np.zeros(len(y)), y, groups))
            sizes = [len(va) for _, va in folds]
            if min(sizes) == 0 or any(len(set(y[va])) < 2 for _, va in folds):
                continue              # an empty or single-class fold cannot be scored
            # Spread relative to an equal split, so the measure is comparable across counts
            imbalance = (max(sizes) - min(sizes)) / (len(y) / n)
            if best is None or imbalance < best[0]:
                best = (imbalance, folds, sizes)
        if best is not None:
            imbalance, folds, sizes = best
            note = "" if n == n_splits else f" (asked for {n_splits}; a fold was unusable above {n})"
            log(f"  grouped folds: {n_groups} groups -> {n} folds of {sizes}, "
                f"largest group {max(np.bincount(pd.factorize(groups)[0])) / len(y):.0%} of the "
                f"samples, imbalance {imbalance:.0%} of an equal fold{note}")
            return folds
    raise ValueError(f"cannot build grouped folds from {n_groups} groups without an empty or "
                     f"single-class fold")


def make_splits(samples, run, cfg):
    """Resolve a run's sample settings into (train_df, [(train_idx, val_idx) per model],
    test_df). Each model trains on train_df.iloc[train_idx] and validates on the val_idx
    rows; test_df is shared and may be empty."""
    label = cfg["label_col"]
    pool = select_samples(samples, run.get("select"), run.get("query"))
    if pool.empty:
        raise ValueError(f"No samples match select={run.get('select')} query={run.get('query')!r}")
    if run.get("limit") and run["limit"] < len(pool):
        pool, _ = train_test_split(pool, train_size=int(run["limit"]), random_state=cfg["seed"],
                                   stratify=pool[label])
        pool = pool.reset_index(drop=True)

    group = run.get("group")
    if run.get("test") or run.get("test_query"):
        test = select_samples(samples, run.get("test"), run.get("test_query"))
        train = pool[~pool[cfg["name_col"]].isin(test[cfg["name_col"]])].reset_index(drop=True)
    elif cfg["test_size"] and group:
        # Hold out whole groups, so the test set shares no group value with training
        n = max(2, int(round(1 / cfg["test_size"])))
        splitter = StratifiedGroupKFold(n_splits=n, shuffle=True, random_state=cfg["seed"])
        tr_i, te_i = next(iter(splitter.split(np.zeros(len(pool)), pool[label],
                                             group_values(pool, group))))
        train, test = pool.iloc[tr_i].reset_index(drop=True), pool.iloc[te_i].reset_index(drop=True)
    elif cfg["test_size"]:
        train, test = train_test_split(pool, test_size=cfg["test_size"], random_state=cfg["seed"],
                                       stratify=pool[label])
        train, test = train.reset_index(drop=True), test.reset_index(drop=True)
    else:
        train, test = pool, pool.iloc[:0]

    y = train[label].to_numpy()
    if group:
        if group not in train:
            raise KeyError(f"No column {group!r} to group by; available: {list(train.columns)}")
        folds = group_folds(y, group_values(train, group), run.get("cv") or cfg["val_folds"],
                            cfg["seed"])
        if not run.get("cv"):
            folds = folds[:1]          # one model, validated on a whole held-out group set
    elif run.get("cv"):
        skf = StratifiedKFold(int(run["cv"]), shuffle=True, random_state=cfg["seed"])
        folds = list(skf.split(np.zeros(len(train)), y))
    else:
        tr, va = train_test_split(np.arange(len(train)), test_size=cfg["val_size"],
                                  random_state=cfg["seed"], stratify=y)
        folds = [(tr, va)]
    return train, folds, test


##############################################
###       Training one model               ###
##############################################

def set_seed(seed):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)


def balance_classes(df, label_col, how="undersample", seed=42):
    """Equal numbers of each class. "undersample" drops majority samples (fast, but throws
    data away); "oversample" repeats minority samples (keeps everything, at the cost of
    seeing duplicates within an epoch)."""
    counts = df[label_col].value_counts()
    if len(counts) < 2:
        return df
    n = counts.min() if how == "undersample" else counts.max()
    parts = [g.sample(n, replace=len(g) < n, random_state=seed) for _, g in df.groupby(label_col)]
    return pd.concat(parts).sample(frac=1, random_state=seed).reset_index(drop=True)


def make_logger(path):
    def log(msg):
        print(msg, flush=True)
        with open(path, "a") as f:
            f.write(msg + "\n")
    return log


def train_one(run_name, run, cfg, train_df, val_df, test_df, fold_name, fold_seed):
    """Build the loaders for one model from the run settings, train it, and write every
    output into out_dir/<run>/<fold_name>/. Returns the metrics row."""
    run_dir = os.path.join(cfg["out_dir"], run_name, fold_name)
    metrics_file = os.path.join(run_dir, "metrics.json")
    if run.get("rerun") and os.path.isdir(run_dir):
        # Start clean: drop the checkpoint, weights, predictions and any importance output
        shutil.rmtree(run_dir)
        shutil.rmtree(os.path.join(cfg["out_dir"], "tensorboard", run.get("tb_group", ""),
                                   run_name, fold_name), ignore_errors=True)
        print(f"Rerunning {run_name}/{fold_name}: removed the previous output")
    elif os.path.exists(metrics_file):
        print(f"Skipping {run_name}/{fold_name}: already finished")
        with open(metrics_file) as f:
            result = json.load(f)
        add_importance(run_name, run, cfg, run_dir, result)
        return result
    os.makedirs(run_dir, exist_ok=True)
    log = make_logger(os.path.join(run_dir, "train.log"))
    set_seed(fold_seed)
    device = torch.device(cfg["device"])
    features = run.get("features", FEATURE_GROUPS)
    name_col, label_col = cfg["name_col"], cfg["label_col"]

    log(f"=== {run_name}/{fold_name}: train {len(train_df)} ({train_df[label_col].mean():.0%} pos), "
        f"val {len(val_df)} ({val_df[label_col].mean():.0%} pos), test {len(test_df)}")
    pd.concat([df.assign(Set=s)[[name_col, label_col, "Set"]] for df, s in
               [(train_df, "train"), (val_df, "val"), (test_df, "test")] if len(df)]
              ).to_csv(os.path.join(run_dir, "split.csv"), index=False)

    # Normalization from this model's training samples only (over all stored channels, so
    # the same stats serve any feature selection)
    kw = dict(features=features, name_col=name_col, label_col=label_col)
    stats = None
    if cfg["normalize"]:
        stats_file = os.path.join(run_dir, "channel_stats.npz")
        if os.path.exists(stats_file):
            with np.load(stats_file) as d:
                stats = (d["mean"], d["std"])
        else:
            t = time.time()
            stats = channel_stats(TCRMapDataset(train_df, **kw))
            np.savez(stats_file, mean=stats[0], std=stats[1])
            log(f"Channel stats from {len(train_df)} training maps in {time.time() - t:.0f}s")

    if cfg["preload"]:
        t = time.time()
        paths = pd.concat([df["Path"] for df in (train_df, val_df, test_df) if len(df)])
        preload_maps(paths)
        log(f"Preloaded {len(paths)} maps in {time.time() - t:.0f}s")

    train_ds = TCRMapDataset(train_df, stats=stats, **kw)
    val_ds = TCRMapDataset(val_df, stats=stats, **kw)
    train_loader = make_loader(train_ds, cfg["batch_size"], shuffle=True, workers=cfg["workers"])
    val_loader = make_loader(val_ds, cfg["batch_size"], workers=cfg["workers"])
    test_loader = None
    if cfg["track_test"] and len(test_df):
        test_loader = make_loader(TCRMapDataset(test_df, stats=stats, **kw), cfg["batch_size"],
                                  workers=cfg["workers"])

    channels = train_ds.getChannels()
    model = build_resnet(channels, img_size=train_ds.canvas_size, layers=cfg["layers"],
                         dropout=cfg["dropout"], conv_dropout=cfg["conv_dropout"],
                         hidden=cfg["hidden"], base_filters=cfg["base_filters"],
                         pool_size=cfg["pool_size"], max_filters=cfg["max_filters"]).to(device)
    optimizer = optim.Adam(model.parameters(), lr=cfg["lr"], weight_decay=cfg["weight_decay"])
    criterion = nn.CrossEntropyLoss()
    log(f"Features {train_ds.groups} -> {channels} channels, "
        f"{sum(p.numel() for p in model.parameters()) / 1e6:.1f}M parameters, device {device}")

    # Everything needed to apply this model to new structures
    with open(os.path.join(run_dir, "model_info.json"), "w") as f:
        json.dump({"run": run_name, "fold": fold_name, "features": train_ds.features,
                   "channels": train_ds.channel_names, "stored_channels": train_ds.stored_channels,
                   "layers": cfg["layers"], "hidden": cfg["hidden"], "dropout": cfg["dropout"],
                   "base_filters": cfg["base_filters"], "pool_size": cfg["pool_size"],
                   "max_filters": cfg["max_filters"],
                   "conv_dropout": cfg["conv_dropout"], "canvas_size": train_ds.canvas_size,
                   "anchor": train_ds.anchor, "pad_value": train_ds.pad_value,
                   "normalized": stats is not None}, f, indent=2)

    # One writer per split in sibling folders: TensorBoard overlays them on one card, and
    # the run name in front of every tag makes each run its own collapsible section
    # TensorBoard's run sidebar is a flat, regex-filterable list, so the folder a run sits
    # in is what lets you scope the view: point --logdir at one tb_group to see only those
    # runs, or type the group into the filter box when looking at everything.
    tb_root = cfg["tb_dir"] or os.path.join(cfg["out_dir"], "tensorboard")
    tb_tag = cfg["tb_tag"] or run_name
    tb_dir = os.path.join(tb_root, run.get("tb_group", ""), run_name, fold_name)
    writers = {s: SummaryWriter(os.path.join(tb_dir, s)) for s in ("train", "val", "test")}
    writers["train"].add_text(f"{tb_tag}/settings",
                              "  \n".join(f"{k}: {v}" for k, v in
                                          {"run": run_name, **run, "channels": channels, **cfg}.items()))
    # Phase 1: pre-train on a class-balanced version of this fold's own train and validation
    # sets, so the model cannot reach a low loss by learning the prior. Balancing happens
    # inside the fold rather than from a separate split, which keeps phase 1 from training
    # on anything phase 2 validates against. Skipped when phase 2 already has a checkpoint,
    # since that checkpoint already holds the pre-trained weights.
    pre_epochs = int(run.get("pretrain_epochs", cfg["pretrain_epochs"]) or 0)
    if pre_epochs and not os.path.exists(os.path.join(run_dir, "checkpoint.pt")):
        how = run.get("pretrain_balance", cfg["pretrain_balance"])
        bal_train = balance_classes(train_df, label_col, how, cfg["seed"])
        bal_val = balance_classes(val_df, label_col, how, cfg["seed"])
        log(f"Pre-training {pre_epochs} epochs on a {how}d split: train {len(bal_train)} "
            f"({bal_train[label_col].mean():.0%} pos, was {len(train_df)} at "
            f"{train_df[label_col].mean():.0%}), val {len(bal_val)} ({bal_val[label_col].mean():.0%} pos)")
        pre_dir = os.path.join(run_dir, "pretrain")
        pre_writers = {s: SummaryWriter(os.path.join(tb_dir, f"pretrain_{s}")) for s in ("train", "val")}
        pre_optimizer = optim.Adam(model.parameters(), lr=cfg["lr"], weight_decay=cfg["weight_decay"])
        train_model(model, pre_optimizer, criterion, pre_epochs,
                    make_loader(TCRMapDataset(bal_train, stats=stats, **kw), cfg["batch_size"],
                                shuffle=True, workers=cfg["workers"]),
                    make_loader(TCRMapDataset(bal_val, stats=stats, **kw), cfg["batch_size"],
                                workers=cfg["workers"]),
                    device, pre_dir, writers=pre_writers, tag_prefix=f"{tb_tag}_pretrain",
                    amp=cfg["amp"], patience=pre_epochs,   # a fixed warm-up, no early stopping
                    min_delta=cfg["min_delta"], checkpoint_every=cfg["checkpoint_every"], log=log)
        best_pre = os.path.join(pre_dir, "best.pt")
        if os.path.exists(best_pre):
            model.load_state_dict(torch.load(best_pre, map_location=device, weights_only=True))
        for w in pre_writers.values():
            w.close()
        result_pretrain = {"Pretrain_epochs": pre_epochs, "Pretrain_N": len(bal_train)}
        log(f"Fine-tuning on the full split: train {len(train_df)} "
            f"({train_df[label_col].mean():.0%} pos)")
    else:
        result_pretrain = {}

    train_model(model, optimizer, criterion, cfg["epochs"], train_loader, val_loader, device,
                run_dir, test_loader=test_loader, writers=writers, tag_prefix=tb_tag,
                amp=cfg["amp"], patience=cfg["patience"],
                min_delta=cfg["min_delta"], checkpoint_every=cfg["checkpoint_every"], log=log)

    # Score the best (lowest validation loss) weights
    model.load_state_dict(torch.load(os.path.join(run_dir, "best.pt"), map_location=device,
                                     weights_only=True))
    hist = pd.read_csv(os.path.join(run_dir, "history.csv"))
    result = {"Run": run_name, "Fold": fold_name, "Features": "+".join(train_ds.features) or "none",
              **result_pretrain,
              "Channels": channels, "N_train": len(train_df), "N_val": len(val_df), "N_test": len(test_df),
              "Best_epoch": int(hist["Val_loss"].idxmin()) + 1, "Epochs_run": len(hist)}
    # Predictions for every split from the best weights, with the called class next to the
    # probability and the annotation columns carried along
    splits = [("Train", train_df)] if cfg["save_train_predictions"] else []
    splits += [("Val", val_df), ("Test", test_df)]
    meta_cols = [c for c in train_df.columns if c not in ("Path", "File")]
    frames = []
    for split, df in splits:
        if not len(df):
            continue
        ds = TCRMapDataset(df, stats=stats, **kw)
        preds, loss = predict(model, make_loader(ds, cfg["batch_size"], workers=cfg["workers"]),
                              device, criterion, cfg["amp"])
        preds["Pred"] = (preds["Prob"] >= cfg["threshold"]).astype(int)
        preds["Correct"] = (preds["Pred"] == preds["Label"]).astype(int)
        preds.insert(1, "Set", split.lower())
        frames.append(preds.merge(df[meta_cols].rename(columns={name_col: "Name"}),
                                  on="Name", how="left"))
        result[f"{split}_loss"] = loss
        result.update({f"{split}_{k}": v for k, v in
                       compute_metrics(preds["Label"], preds["Prob"], cfg["threshold"]).items()})
        if split == "Test" and test_loader is None:
            # Not tracked per epoch, so mark the single best-weights score on the cards
            for k in ("Accuracy", "ROC_AUC", "PR_AUC"):
                writers["test"].add_scalar(f"{tb_tag}/{k}", result[f"Test_{k}"], result["Best_epoch"])
            writers["test"].add_scalar(f"{tb_tag}/Loss", loss, result["Best_epoch"])
    predictions = pd.concat(frames, ignore_index=True)
    predictions.to_csv(os.path.join(run_dir, "predictions.csv"), index=False)
    for w in writers.values():
        w.close()
    with open(metrics_file, "w") as f:
        json.dump(result, f, indent=2, default=float)
    if "Test_Accuracy" in result:
        log("Test: acc {:.3f}, ROC AUC {:.3f}, PR AUC {:.3f}, F1 {:.3f}".format(
            result["Test_Accuracy"], result["Test_ROC_AUC"], result["Test_PR_AUC"], result["Test_F1"]))

    del model, optimizer, train_loader, val_loader
    if device.type == "cuda":
        torch.cuda.empty_cache()
    add_importance(run_name, run, cfg, run_dir, result, log=log)
    return result


def add_importance(run_name, run, cfg, run_dir, result, log=print):
    """Run the gradient feature importance for a finished model when the run asks for it
    (importance.py is imported here so matplotlib is only needed for runs that use it)."""
    spec = run.get("importance")
    if not spec:
        return
    done = os.path.join(run_dir, "importance", "channel_importance.csv")
    if os.path.exists(done):
        log(f"Importance for {run_name}/{os.path.basename(run_dir)}: already done")
        return
    opts = dict(spec) if isinstance(spec, dict) else {}
    # Default to the test predictions, but fall back to validation when a run has no test set
    opts.setdefault("split", "test" if result.get("N_test") else "val")
    from importance import run_importance
    try:
        run_importance(run_dir, inputs=cfg["inputs"], device=cfg["device"], **opts)
    except Exception as e:
        log(f"Importance failed for {run_name}/{os.path.basename(run_dir)}: {e}")


##############################################
###       Running                          ###
##############################################

def run_model(run_name, run=None, cfg=None, samples=None):
    """Train one run: every CV fold, or a single model when the run sets no cv. Returns a
    list of metrics rows."""
    run = RUNS[run_name] if run is None else run
    base = {**CONFIG, **(cfg or {})}
    tuned = {}
    use = run.get("use_params", base["use_params"])
    if use:
        from hypersearch import load_params
        key = use if isinstance(use, str) else base["params_key"]
        tuned = load_params(key, base)
        if tuned:
            print(f"{run_name}: tuned parameters {key!r} - "
                  + ", ".join(f"{k}={v:.3g}" if isinstance(v, float) else f"{k}={v}"
                              for k, v in tuned.items()))
    # Run settings win over tuned parameters, which win over the config defaults
    cfg = {**CONFIG, **(cfg or {}), **tuned, **{k: v for k, v in run.items() if k in CONFIG}}
    samples = load_samples(cfg) if samples is None else samples
    train, folds, test = make_splits(samples, run, cfg)
    cv = bool(run.get("cv"))
    print(f"\n### {run_name}: {len(folds)} model{'s' if cv else ''} "
          f"({'%d-fold CV' % len(folds) if cv else 'single split'}), "
          f"features {run.get('features', FEATURE_GROUPS)}")
    describe(train, cfg["label_col"], prefix="  train pool: ")
    if len(test):
        describe(test, cfg["label_col"], prefix="  test set:   ")

    results = []
    for k, (tr, va) in enumerate(folds):
        results.append(train_one(run_name, run, cfg, train.iloc[tr], train.iloc[va], test,
                                 f"fold{k}" if cv else "single", cfg["seed"] + k))
    return results


def run_all(runs=None, cfg=None):
    """Train every run in RUNS (or just the named ones). Returns the results DataFrame."""
    if not RUNS:
        raise ValueError("No runs defined - add them to RUNS in run_models.py")
    selected = runs or list(RUNS)
    unknown = [r for r in selected if r not in RUNS]
    if unknown:
        raise ValueError(f"Unknown runs {unknown}; defined: {list(RUNS)}")
    cfg = {**CONFIG, **(cfg or {})}
    os.makedirs(cfg["out_dir"], exist_ok=True)
    samples = load_samples(cfg)
    results = []
    for run_name in selected:
        results += run_model(run_name, RUNS[run_name], cfg, samples)
        write_results(results, cfg["out_dir"])
    return pd.DataFrame(results)


def write_results(results, out_dir):
    df = pd.DataFrame(results)
    df.to_csv(os.path.join(out_dir, "results.csv"), index=False)
    metrics = [c for c in df.columns if c.startswith(("Val_", "Test_")) and df[c].dtype != object]
    summary = df.groupby(["Run", "Features"], sort=False)[metrics].agg(["mean", "std"])
    summary.columns = [f"{a}_{b}" for a, b in summary.columns]
    summary.reset_index().to_csv(os.path.join(out_dir, "results_summary.csv"), index=False)
    return df


def list_runs(runs=None, cfg=None):
    """Print each run's splits and channel count without training anything."""
    cfg = {**CONFIG, **(cfg or {})}
    if not RUNS:
        print("No runs defined - add them to RUNS in run_models.py")
        return
    selected = {k: RUNS[k] for k in (runs or list(RUNS))}
    samples = load_samples(cfg)
    for run_name, run in selected.items():
        rcfg = {**cfg, **{k: v for k, v in run.items() if k in CONFIG}}
        train, folds, test = make_splits(samples, run, rcfg)
        ds = TCRMapDataset(train.head(1), features=run.get("features", FEATURE_GROUPS),
                           name_col=rcfg["name_col"], label_col=rcfg["label_col"])
        tr, va = folds[0]
        print(f"  {run_name:28s} {len(folds)} model(s)  train {len(tr):5d}  val {len(va):4d}  "
              f"test {len(test):5d}  {ds.getChannels():2d} channels {ds.groups}")


def main():
    p = argparse.ArgumentParser(description="Train the models defined in RUNS")
    p.add_argument("runs", nargs="*",
                   help="run names to train; given, they override the default of training every run in RUNS")
    p.add_argument("--list", action="store_true", help="print the runs and their splits, then exit")
    args = p.parse_args()
    if args.list:
        list_runs(runs=args.runs or None)
        return
    run_all(args.runs or None)   # no names: train every run in RUNS


if __name__ == "__main__":
    main()

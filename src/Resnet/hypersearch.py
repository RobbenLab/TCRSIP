#!/usr/bin/env python3

# Hyperparameter search with cross validation for the TCR ResNet, coarse to fine.
#
# Each trial is a sampled parameter set trained through the ordinary run logic in
# run_models.py (same splits, loaders, checkpointing and TensorBoard) - this module only
# samples the parameters, scores the trials and records the winner. Trials are scored on
# generalization (validation loss and how far it separates from the training loss), across
# the CV folds; the test set is never involved.
#
# It is a separate step from run_models.py: search once, then every run picks the winning
# parameters up automatically from results/best_params.json (CONFIG["use_params"]).
#
#   python src/hypersearch.py                    # search on SEARCH_RUN below, save as "default"
#   python src/hypersearch.py --configs 24 8 3 --round-epochs 10 30 60
#   python src/hypersearch.py --metric Performance            # strongest model, no bias terms
#   python src/hypersearch.py --strategy halving --rounds 8 24 60   # cheaper, survivors continue
#   python src/hypersearch.py --strategy random --epochs 30         # one flat round
#   python src/hypersearch.py --run test_all     # use that run's samples/features instead
#   python src/hypersearch.py --show             # print the saved parameters
#   python src/run_models.py                     # the runs now use them
#
#   # in a notebook
#   from hypersearch import search, get_params
#   params = get_params()                  # saved parameters, or run the search if there are none
#
# Trials go to results/hyperparameters/<key>/trial00.../fold0..., with trials.csv next to
# them. A finished trial is skipped on a rerun, so an interrupted search resumes where it
# stopped.

import argparse
import json
import os
import shutil

import numpy as np
import pandas as pd
import torch

import run_models
from run_models import CONFIG, RUNS, load_samples
from tcr_dataset import FEATURE_GROUPS

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Samples and features the search trains on - one search serves every run. The curated
# sources are used rather than everything: they are close to balanced (55% positive against
# 37% for the full set, which is dragged down by immrep at 17%), and on a skewed pool the
# cheapest way to a low loss is to predict the majority class, which is not the behaviour
# worth tuning for. "limit" keeps it affordable; raise it (or drop it) for a fuller search.
SEARCH_RUN = {"select": {"source": ["exp", "fake", "ground"]}, "features": FEATURE_GROUPS,
              "limit": 600}

# Key the winning parameters are stored under in best_params.json, and the one every run
# reads by default
PARAMS_KEY = "default"

# Sampled per trial: ("loguniform", lo, hi) | ("uniform", lo, hi) | ("choice", [values])
SPACE = {
    # Architecture: one entry per residual stage, so the list length is the depth. The
    # original [3, 4, 6, 3] was picked to fit a 1024 input on the GPU, not for accuracy;
    # at 512 the shallower nets fit easily and carry far fewer parameters (2.9M for
    # [2, 2] against 29.9M for [3, 4, 6, 3]), which matters when the problem is overfitting.
    # Six stages is the ceiling: the input is 512, so each stage after the first halves the
    # map down to 4x4 at six, 2x2 at seven and a single cell at eight - past six there is
    # no spatial structure left for the pooling to summarize.
    "layers": ("choice", [[2, 2], [2, 2, 2], [3, 4, 6], [2, 2, 2, 2], [3, 4, 6, 3],
                          [3, 4, 6, 3, 3], [2, 2, 2, 2, 2, 2], [3, 4, 6, 3, 3, 3]]),
    "base_filters": ("choice", [32, 64]),        # channels in the first stage
    # Cap on the channel doubling: without one, a 7-stage net reaches 4096 channels and
    # 783M parameters, so depth and size cannot be told apart
    "max_filters": ("choice", [256, 512, None]),
    "hidden": ("choice", [256, 512, 1024]),      # width of the classifier layer
    # Discrete half-decade grids rather than continuous draws: the same value recurs across
    # trials, so a difference between two trials is the other settings and not a slightly
    # different learning rate. Refinement moves to a neighbouring grid value.
    "lr": ("choice", [1e-6,  1e-5, 1e-4, 1e-3, 1e-2]),        # notebook default 1e-5
    "weight_decay": ("choice", [1.0, 0.5, 0.05, 0.005]),       # notebook default 0.1
    "batch_size": ("choice", [8, 16, 32]),       # notebook default 16
    "dropout": ("choice", [0.0, 0.2, 0.4, 0.6]),  # hidden layer; 0 = the notebook's behaviour
    "conv_dropout": ("choice", [0.0, 0.1, 0.2]),  # Dropout2d on the conv stack
}

# Trial settings: shorter than a real run, CV so the score is not one lucky split.
# The search is about generalization, so trials are scored on validation loss, not on a
# discrimination metric - a model that overfits still has a high training AUC and often a
# decent validation AUC while its validation loss climbs. Alternatives for `metric`:
# All of the composite metrics are built from terms that are 1.0 at chance, so no term wins
# on scale alone. Performance = Val_loss/ln(2) + (1 - Val_ROC_AUC)/0.5 - calibration and
# ranking together, which is the "best model" term the others build on.
#   Performance    quality only. Overfitting is already penalised through the validation
#                  loss, so this is the honest default when you want the strongest model
#   Balanced       (default) Performance + divergence_weight * Divergence. Same, with an
#                  explicit penalty for the training and validation curves separating.
#                  Keep divergence_weight well below 1: at 1.0 a model that learns nothing
#                  has zero divergence and beats a better model that generalised fine
#   Generalization Balanced plus a penalty for models still improving when training ended
#                  (underfit), so neither failure mode can win
#   AUCDiv         (1 - Val_ROC_AUC)/0.5 + w * Divergence - ranking only, ignores calibration
#   LossDiv        Val_loss/ln(2) + w * Divergence - calibration only, ignores ranking
#   Divergence     area between the validation and training curves over the run, divided by
#                  the area under the validation curve. Pure overfitting, no quality term:
#                  a model that learns nothing scores 0, so read it, do not optimise it
#   Underfit       best epoch / epochs run; 1.0 = still improving when training stopped
#   Val_loss / Gap / Rise / Overfit / Val_ROC_AUC / Val_PR_AUC / Val_Accuracy
# n_iter is the number of parameter sets tried. With params_per_fold (the default) each
# one is trained on a different CV fold, taken round-robin, so the whole budget goes into
# sampling the space instead of re-testing the same parameters on every fold. The winner is
# then re-trained across all cv folds (confirm_best) to check it was not just a lucky fold.
#
# Strategy: "halving" (default) spends the budget coarse to fine. Every config is trained
# for rounds[0] epochs, the worst are dropped, and the survivors *continue from their
# checkpoints* for the next round rather than restarting - so the promising settings end up
# with more training than a flat search would give them, and the hopeless ones are cut
# after a few epochs. After the first round it also resamples around the survivors
# (narrowing the continuous ranges), which is the fine-grained half of the search.
# "random" is the previous behaviour: every config trained for `epochs` and ranked once.
SEARCH_DEFAULTS = {"n_iter": 30, "cv": 5, "epochs": 30, "patience": 8,
                   "metric": "Balanced", "higher_is_better": False,
                   "divergence_weight": 0.25, "underfit_weight": 0.5,
                   "gap_weight": 1.0, "rise_weight": 1.0,
                   "params_per_fold": True, "confirm_best": True,
                   "strategy": "resample",
                   # resample strategy: (configs, epochs) per round. Many short runs explore
                   # the space, then fewer longer ones around whatever that found - so the
                   # budget concentrates on the promising region, but every model inside a
                   # round is trained for the same epochs and so is directly comparable.
                   "schedule": [[24, 10], [8, 30], [3, 60]],
                   "carry_forward": 1,      # best configs retrained at the next round's budget
                   "rounds": [8, 24, 60],   # cumulative epoch budget per round
                   "factor": 3,             # keep 1/factor of the configs each round
                   "refine": True,          # resample around the first round's survivors
                   "refine_n": None,        # how many (default: half of n_iter)
                   "refine_width": 3.0}     # continuous values resampled within x/÷ this

# Diagnostics computed per fold from history.csv, averaged over the folds of a trial
FOLD_DIAGNOSTICS = ["Val_loss", "Loss_ratio", "Train_loss", "Gap", "Divergence", "Rise", "Underfit",
                    "Val_ROC_AUC", "Val_PR_AUC", "Val_Accuracy", "Best_epoch", "Epochs_run"]


def sample_params(space, rng):
    params = {}
    for name, spec in space.items():
        kind = spec[0]
        if kind == "loguniform":
            params[name] = float(np.exp(rng.uniform(np.log(spec[1]), np.log(spec[2]))))
        elif kind == "uniform":
            params[name] = float(rng.uniform(spec[1], spec[2]))
        elif kind == "choice":
            value = spec[1][rng.integers(len(spec[1]))]
            params[name] = value.item() if hasattr(value, "item") else value
        else:
            raise ValueError(f"Unknown spec {spec!r} for {name}")
    return params


def sample_near(space, survivors, rng, width=3.0):
    """Sample a parameter set near the surviving ones: continuous values are drawn from a
    band around a randomly chosen survivor (clipped to the original range), categorical
    values from those the survivors actually used."""
    anchor = survivors[int(rng.integers(len(survivors)))]
    params = {}
    for name, spec in space.items():
        kind, value = spec[0], anchor[name]
        if kind == "loguniform":
            lo = max(spec[1], value / width)
            hi = min(spec[2], value * width)
            params[name] = float(np.exp(rng.uniform(np.log(lo), np.log(hi))))
        elif kind == "uniform":
            span = (spec[2] - spec[1]) * 0.15
            params[name] = float(np.clip(rng.normal(value, span), spec[1], spec[2]))
        elif kind == "choice":
            options = list(spec[1])
            seen = [s[name] for s in survivors]
            unique = [v for i, v in enumerate(seen) if v not in seen[:i]]
            numeric = all(isinstance(v, (int, float)) and not isinstance(v, bool) for v in options)
            if numeric:
                # An ordered grid: step to a neighbouring value so refinement can move off
                # the surviving values instead of only reshuffling them
                anchor_i = options.index(value)
                step = int(rng.integers(-1, 2))
                params[name] = options[int(np.clip(anchor_i + step, 0, len(options) - 1))]
            else:
                params[name] = unique[int(rng.integers(len(unique)))]
        else:
            raise ValueError(f"Unknown spec {spec!r} for {name}")
    return params


def _chance_loss(pos_rate):
    """Cross entropy of always predicting the class prior."""
    if pos_rate is None or not 0 < pos_rate < 1:
        return float(np.log(2))
    p = float(pos_rate)
    return float(-(p * np.log(p) + (1 - p) * np.log(1 - p)))


def fold_diagnostics(history, pos_rate=None):
    """Overfitting diagnostics for one fold's history.csv. The best epoch is the lowest
    validation loss, matching how training picks its weights.
    pos_rate: the validation set's positive rate. The loss of always predicting that rate is
    what "no better than chance" costs, and it is only ln(2) when the classes are balanced -
    at 17% positive it is 0.45, so without this a trivial model looks strong on loss."""
    h = history.reset_index(drop=True)
    best = int(h["Val_loss"].idxmin())
    val_min = float(h.loc[best, "Val_loss"])
    # How far validation loss climbed after its minimum - the overfitting signal itself
    after = h.loc[best:, "Val_loss"]
    rise = float(after.max() - val_min) if len(after) > 1 else 0.0
    # Area between the curves over the whole run, divided by the area under the validation
    # curve: a dimensionless measure of how far the two separated, rather than their
    # distance at one epoch. Only the part where validation sits above training counts -
    # training loss is measured with dropout on and averaged over the epoch, so it sits
    # above validation early in a run without that meaning anything.
    above = np.maximum(h["Val_loss"].to_numpy() - h["Train_loss"].to_numpy(), 0.0)
    if len(h) > 1:
        divergence = float(np.trapezoid(above) / max(np.trapezoid(h["Val_loss"].to_numpy()), 1e-8))
    else:
        divergence = float(above[0] / max(h.loc[0, "Val_loss"], 1e-8))
    return {"Val_loss": val_min,
            "Train_loss": float(h.loc[best, "Train_loss"]),
            "Gap": float(h.loc[best, "Val_loss"] - h.loc[best, "Train_loss"]),
            "Divergence": divergence,
            "Rise": rise,
            "Val_ROC_AUC": float(h.loc[best, "Val_ROC_AUC"]),
            "Val_PR_AUC": float(h.loc[best, "Val_PR_AUC"]),
            "Val_Accuracy": float(h.loc[best, "Val_Accuracy"]),
            "Loss_ratio": val_min / _chance_loss(pos_rate),
            "Best_epoch": best + 1,
            "Epochs_run": len(h),
            # 1.0 means the best epoch was the last one trained: validation loss was still
            # falling when training ended, so the model is underfit rather than converged
            "Underfit": (best + 1) / len(h)}


def score_trial(diagnostics, opts):
    """Collapse a trial's mean diagnostics into the single number trials are ranked by.
    The composite metrics pair a quality term with a divergence term, both dimensionless
    and both 1.0 at chance, so neither dominates by scale:
      Balanced = (1 - Val_ROC_AUC) / 0.5  +  w * Divergence     (discrimination + overfitting)
      LossDiv  = Val_loss / ln(2)         +  w * Divergence     (calibration + overfitting)
    A divergence term alone is not enough: a model that learns nothing has identical train
    and validation curves and would score a perfect 0, so the quality term has to be there."""
    d, w = diagnostics, opts["divergence_weight"]
    # Quality: calibration (loss) and ranking (AUC), each 1.0 at chance so neither dominates
    # Loss_ratio is the validation loss over what always predicting the class prior costs,
    # so 1.0 is "no better than the prior" whatever the balance (older trial tables without
    # the column fall back to the balanced assumption)
    loss_term = d["Loss_ratio"] if "Loss_ratio" in d else d["Val_loss"] / np.log(2)
    performance = loss_term + (1.0 - d["Val_ROC_AUC"]) / 0.5
    if opts["metric"] == "Performance":
        return performance
    if opts["metric"] == "Balanced":
        return performance + w * d["Divergence"]
    if opts["metric"] == "Generalization":
        # Adds an underfitting penalty: a model still improving when training stopped has
        # not shown what it can do, and its curves have not had a chance to separate
        return (performance + w * d["Divergence"]
                + opts["underfit_weight"] * max(0.0, d["Underfit"] - 0.9) / 0.1)
    if opts["metric"] == "AUCDiv":
        return (1.0 - d["Val_ROC_AUC"]) / 0.5 + w * d["Divergence"]
    if opts["metric"] == "LossDiv":
        return loss_term + w * d["Divergence"]
    if opts["metric"] == "Overfit":
        return (d["Val_loss"] + opts["gap_weight"] * max(0.0, d["Gap"])
                + opts["rise_weight"] * d["Rise"])
    if opts["metric"] not in d:
        raise ValueError(f"Unknown metric {opts['metric']!r}; available: {sorted(d)}, "
                         f"'Performance', 'Balanced', 'Generalization', 'AUCDiv', 'LossDiv' "
                         f"or 'Overfit'")
    return d[opts["metric"]]


def params_path(cfg):
    return cfg.get("params_file") or os.path.join(cfg["out_dir"], "best_params.json")


def save_params(run_name, params, cfg, extra=None):
    """Merge this run's winning parameters into the shared best_params.json."""
    path = params_path(cfg)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    saved = {}
    if os.path.exists(path):
        with open(path) as f:
            saved = json.load(f)
    saved[run_name] = {**params, **(extra or {})}
    with open(path, "w") as f:
        json.dump(saved, f, indent=2, default=float)
    return path


def load_params(run_name, cfg=None):
    """The saved parameters for a run name, or {} if the search has not been run."""
    cfg = {**CONFIG, **(cfg or {})}
    path = params_path(cfg)
    if not os.path.exists(path):
        return {}
    with open(path) as f:
        saved = json.load(f)
    entry = saved.get(run_name, {})
    return {k: v for k, v in entry.items() if k in CONFIG}


def _setup(key, base_run, space, cfg, opts, samples):
    """Shared context: where trials are written, the CV split they share, and the config
    every trial is trained with."""
    base_run = dict(SEARCH_RUN if base_run is None else base_run)
    search_dir = os.path.join(cfg["out_dir"], "hyperparameters", key)
    os.makedirs(search_dir, exist_ok=True)
    # Trials share the run's samples and features but never its flags: no importance, no
    # rerun (so an interrupted search resumes), no tuned parameters bleeding in
    base_run = {k: v for k, v in base_run.items()
                if k not in ("importance", "rerun", "cv", "epochs", "patience", *space)}
    # Trials log into the main TensorBoard tree under one collapsible section per search,
    # so every trial's train/val curves land on the same cards and can be compared there
    trial_cfg = {**cfg, "out_dir": search_dir, "track_test": False, "use_params": False,
                 "tb_dir": os.path.join(cfg["out_dir"], "tensorboard", "hyperparameters"),
                 "tb_tag": f"hp_{key}"}
    cv = int(opts["cv"])
    # One split of the samples, reused by every trial so parameter sets are comparable
    train_df, folds, test_df = run_models.make_splits(samples, {**base_run, "cv": cv},
                                                      {**trial_cfg, **base_run})
    return {"search_dir": search_dir, "base_run": base_run, "trial_cfg": trial_cfg,
            "train": train_df, "folds": folds, "test": test_df, "cv": cv, "samples": samples,
            "key": key, "cfg": cfg}


def _show(params):
    return ", ".join(f"{k}={v:.3g}" if isinstance(v, float) else f"{k}={v}"
                     for k, v in params.items())


def _reopen_for_more_epochs(run_dir, budget):
    """Prepare a finished trial to continue training to a larger budget. Returns "done" if
    it has nothing left to give (it early-stopped, or already reached the budget)."""
    checkfile = os.path.join(run_dir, "checkpoint.pt")
    if not os.path.exists(checkfile):
        return "fresh"
    ck = torch.load(checkfile, map_location="cpu", weights_only=False)
    if ck.get("stopped"):
        return "done"                       # early stopped: more epochs would not help
    if len(ck.get("history", [])) >= budget:
        return "done"                       # already trained this far
    ck["finished"] = False                  # let train_model pick up where it left off
    torch.save(ck, checkfile)
    metrics = os.path.join(run_dir, "metrics.json")
    if os.path.exists(metrics):
        os.remove(metrics)                  # so train_one does not skip it as finished
    return "continue"


def _clear_if_stale(run_dir, params):
    """Trial names are reused between searches, so a folder holding a different parameter
    set belongs to an older run - its checkpoint would not load into this architecture."""
    saved = os.path.join(run_dir, "params.json")
    if not os.path.exists(saved):
        return
    with open(saved) as f:
        old = json.load(f)
    if any(old.get(k) != v for k, v in params.items()):
        print(f"  {os.path.basename(os.path.dirname(run_dir))}: previous trial used different "
              f"parameters, retraining from scratch")
        shutil.rmtree(run_dir)


def _save_params(run_dir, name, params, fold_k, budget, stage):
    _clear_if_stale(run_dir, params)
    os.makedirs(run_dir, exist_ok=True)
    with open(os.path.join(run_dir, "params.json"), "w") as f:
        json.dump({"trial": name, "fold": fold_k, "epochs": budget, "stage": stage,
                   **params}, f, indent=2, default=float)


def _run_trial(ctx, name, params, fold_k, budget, opts, seed):
    """Train a trial to `budget` epochs (continuing from its checkpoint if it has one) and
    return its diagnostics."""
    fold_name = f"fold{fold_k}"
    run_dir = os.path.join(ctx["search_dir"], name, fold_name)
    state = _reopen_for_more_epochs(run_dir, budget)
    if state != "done":
        run = {**ctx["base_run"], **params, "epochs": budget, "patience": opts["patience"],
               "importance": False}
        cfg_i = {**ctx["trial_cfg"], **{k: v for k, v in run.items() if k in CONFIG}}
        tr, va = ctx["folds"][fold_k]
        run_models.train_one(name, run, cfg_i, ctx["train"].iloc[tr], ctx["train"].iloc[va],
                             ctx["test"], fold_name, seed)
    hist = os.path.join(run_dir, "history.csv")
    if not os.path.exists(hist):
        raise RuntimeError(f"{name}: no history.csv under {run_dir}")
    return fold_diagnostics(pd.read_csv(hist), _fold_pos_rate(run_dir)), state


def _fold_pos_rate(run_dir):
    """Validation positive rate for a trained fold, from the metrics it wrote."""
    metrics = os.path.join(run_dir, "metrics.json")
    if not os.path.exists(metrics):
        return None
    with open(metrics) as f:
        m = json.load(f)
    return m["Val_N_pos"] / m["Val_N"] if m.get("Val_N") else None


def _record(name, params, fold_k, budget, diag, opts, extra=None):
    row = {"Trial": name, "Fold": f"fold{fold_k}", "Budget": budget, **params,
           "Score": score_trial(diag, opts),
           **{k: diag[k] for k in FOLD_DIAGNOSTICS}, **(extra or {})}
    return row


def _report(name, budget, row, opts, note=""):
    print(f"--- {name} @{budget}ep: {opts['metric']} {row['Score']:.3f}  "
          f"(val loss {row['Val_loss']:.3f} = {row['Loss_ratio']:.2f}x prior, "
          f"train {row['Train_loss']:.3f}, "
          f"div {row['Divergence']:.3f}, rise {row['Rise']:.3f}, val auc {row['Val_ROC_AUC']:.3f}, "
          f"best epoch {row['Best_epoch']:.0f}/{row['Epochs_run']:.0f}){note}")


def _rank(rows, opts):
    return sorted(rows, key=lambda r: r["Score"], reverse=bool(opts["higher_is_better"]))


def _search_halving(ctx, space, opts, rng, rows):
    """Coarse to fine: train everything briefly, keep the best, let the survivors continue
    from their checkpoints for longer, and resample around them in between."""
    rounds = list(opts["rounds"])
    n = int(opts["n_iter"])
    cv = ctx["cv"]
    seed = ctx["cfg"]["seed"]

    def evaluate(name, params, fold_k, budget, stage, index):
        _save_params(os.path.join(ctx["search_dir"], name, f"fold{fold_k}"), name, params,
                     fold_k, budget, stage)
        diag, state = _run_trial(ctx, name, params, fold_k, budget, opts, seed + index)
        row = _record(name, params, fold_k, budget, diag, opts, {"Stage": stage})
        rows.append(row)
        _report(name, budget, row, opts, "  (converged)" if state == "done" else "")
        pd.DataFrame(rows).to_csv(os.path.join(ctx["search_dir"], "trials.csv"), index=False)
        return row

    # Round 1, coarse: every config gets the same short budget
    print(f"\n=== Round 1/{len(rounds)}: {n} configs x {rounds[0]} epochs (coarse)")
    live = []
    for i in range(n):
        params = sample_params(space, rng)
        name = f"trial{i:03d}"
        fold_k = i % cv
        print(f"\n--- {name} (fold{fold_k}): {_show(params)}")
        row = evaluate(name, params, fold_k, rounds[0], "coarse", i)
        live.append({"name": name, "params": params, "fold": fold_k, "row": row, "index": i})

    # Refinement: new configs drawn around the survivors, judged at the same budget
    if opts["refine"] and len(rounds) > 1:
        keep = max(2, len(live) // int(opts["factor"]))
        top_names = [row["Trial"] for row in _rank([t["row"] for t in live], opts)[:keep]]
        by_name = {t["name"]: t["params"] for t in live}
        survivors = [by_name[n] for n in top_names]
        n_refine = int(opts["refine_n"] or max(1, n // 2))
        print(f"\n=== Refinement: {n_refine} configs resampled around the top {keep}, "
              f"{rounds[0]} epochs")
        for j in range(n_refine):
            params = sample_near(space, survivors, rng, opts["refine_width"])
            name = f"trial{n + j:03d}"
            fold_k = (n + j) % cv
            print(f"\n--- {name} (fold{fold_k}, refined): {_show(params)}")
            row = evaluate(name, params, fold_k, rounds[0], "refine", n + j)
            live.append({"name": name, "params": params, "fold": fold_k, "row": row,
                         "index": n + j})

    # Halving rounds: the survivors continue from their checkpoints
    for r, budget in enumerate(rounds[1:], start=2):
        keep = max(1, len(live) // int(opts["factor"]))
        order = _rank([t["row"] for t in live], opts)
        kept_names = {row["Trial"] for row in order[:keep]}
        live = [t for t in live if t["name"] in kept_names]
        print(f"\n=== Round {r}/{len(rounds)}: top {len(live)} continue to {budget} epochs")
        for t in live:
            t["row"] = evaluate(t["name"], t["params"], t["fold"], budget, f"round{r}", t["index"])

    return _rank([t["row"] for t in live], opts)[0]


def _search_resample(ctx, space, opts, rng, rows):
    """Rounds of fresh models on a growing budget. Round 1 samples the whole space with many
    short runs; each later round trains fewer, longer models drawn around the best of the
    previous round, and retrains the best configs themselves at the new budget. Models are
    never continued from a checkpoint, so within a round every score is on equal terms, and
    the extra epochs go to the region that looked good rather than to whichever individual
    model happened to survive."""
    schedule = [(int(n), int(e)) for n, e in opts["schedule"]]
    cv, seed = ctx["cv"], ctx["cfg"]["seed"]
    index = 0
    previous = []          # [{params, row}] from the round just finished

    for r, (n_configs, budget) in enumerate(schedule, start=1):
        if r == 1:
            configs = [(sample_params(space, rng), "new") for _ in range(n_configs)]
            print(f"\n=== Round {r}/{len(schedule)}: {n_configs} configs x {budget} epochs "
                  f"(sampling the whole space)")
        else:
            keep = max(2, len(previous) // int(opts["factor"]))
            top = _rank([p["row"] for p in previous], opts)[:keep]
            names = [row["Trial"] for row in top]
            by_name = {p["row"]["Trial"]: p["params"] for p in previous}
            survivors = [by_name[nm] for nm in names]
            # The best configs earn more epochs: retrained from scratch at the new budget
            carried = [(by_name[nm], "carried")
                       for nm in names[:min(int(opts["carry_forward"]), n_configs)]]
            configs = carried + [(sample_near(space, survivors, rng, opts["refine_width"]), "new")
                                 for _ in range(n_configs - len(carried))]
            print(f"\n=== Round {r}/{len(schedule)}: {n_configs} configs x {budget} epochs "
                  f"({len(carried)} carried forward, {n_configs - len(carried)} resampled "
                  f"around the best {len(survivors)} of {len(previous)})")

        current = []
        for params, origin in configs:
            name = f"trial{index:03d}"
            fold_k = index % cv
            print(f"\n--- {name} (fold{fold_k}, round {r}, {origin}): {_show(params)}")
            _save_params(os.path.join(ctx["search_dir"], name, f"fold{fold_k}"), name, params,
                         fold_k, budget, f"round{r}")
            diag, _ = _run_trial(ctx, name, params, fold_k, budget, opts, seed + index)
            row = _record(name, params, fold_k, budget, diag, opts,
                          {"Stage": f"round{r}", "Origin": origin})
            rows.append(row)
            current.append({"params": params, "row": row})
            _report(name, budget, row, opts)
            pd.DataFrame(rows).to_csv(os.path.join(ctx["search_dir"], "trials.csv"), index=False)
            index += 1
        previous = current

    # The winner comes from the last round, where every model had the full budget
    return _rank([p["row"] for p in previous], opts)[0]


def _search_random(ctx, space, opts, rng, rows):
    """Every config trained for the same budget, ranked once (the original behaviour)."""
    cv, seed = ctx["cv"], ctx["cfg"]["seed"]
    for i in range(int(opts["n_iter"])):
        params = sample_params(space, rng)
        name = f"trial{i:03d}"
        if opts["params_per_fold"]:
            fold_k = i % cv
            print(f"\n--- {name} (fold{fold_k}): {_show(params)}")
            diag, _ = _run_trial(ctx, name, params, fold_k, opts["epochs"], opts, seed + i)
            row = _record(name, params, fold_k, opts["epochs"], diag, opts, {"Stage": "flat"})
        else:
            print(f"\n--- {name}: {_show(params)}")
            run = {**ctx["base_run"], **params, "cv": cv, "epochs": opts["epochs"],
                   "patience": opts["patience"], "importance": False}
            results = run_models.run_model(name, run, ctx["trial_cfg"], ctx["samples"])
            per_fold = [fold_diagnostics(
                pd.read_csv(os.path.join(ctx["search_dir"], name, str(x["Fold"]), "history.csv")),
                _fold_pos_rate(os.path.join(ctx["search_dir"], name, str(x["Fold"]))))
                for x in results]
            mean = {k: float(np.nanmean([d[k] for d in per_fold])) for k in FOLD_DIAGNOSTICS}
            row = {"Trial": name, "Fold": "all", "Budget": opts["epochs"], **params,
                   "Score": float(np.nanmean([score_trial(d, opts) for d in per_fold])),
                   **mean, "Stage": "flat"}
        rows.append(row)
        _report(name, row["Budget"], row, opts)
        pd.DataFrame(rows).to_csv(os.path.join(ctx["search_dir"], "trials.csv"), index=False)
    return _rank(rows, opts)[0]


def search(key=PARAMS_KEY, base_run=None, space=None, cfg=None, samples=None, force=False,
           **options):
    """Hyperparameter search, saved under `key`. base_run gives the samples and features to
    search on (default SEARCH_RUN). force deletes the previous search for this key first -
    needed when the search space changed, since trial names are reused and an old
    checkpoint cannot be loaded into a different architecture.
    Returns (best_params, trials DataFrame)."""
    opts = {**SEARCH_DEFAULTS, **options}
    space = space or SPACE
    cfg = {**CONFIG, **(cfg or {})}
    samples = load_samples(cfg) if samples is None else samples
    if force:
        for stale in (os.path.join(cfg["out_dir"], "hyperparameters", key),
                      os.path.join(cfg["out_dir"], "tensorboard", "hyperparameters")):
            if os.path.isdir(stale):
                print(f"--force: removing the previous search at {stale}")
                shutil.rmtree(stale)
    ctx = _setup(key, base_run, space, cfg, opts, samples)
    rng = np.random.default_rng(cfg["seed"])
    direction = "lowest" if not opts["higher_is_better"] else "highest"

    if opts["strategy"] == "resample":
        sched = [(int(n), int(e)) for n, e in opts["schedule"]]
        total = sum(n * e for n, e in sched)
        shape = " -> ".join(f"{n}x{e}ep" for n, e in sched)
        print(f"Hyperparameter search ({key}): {shape} ({sum(n for n, _ in sched)} models, "
              f"~{total} epochs), fresh models each round drawn around the best of the last; "
              f"one CV fold each (of {opts['cv']}); keeping the {direction} {opts['metric']}")
    elif opts["strategy"] == "halving":
        rounds = list(opts["rounds"])
        total = opts["n_iter"] * rounds[0]
        if opts["refine"] and len(rounds) > 1:
            total += int(opts["refine_n"] or max(1, opts["n_iter"] // 2)) * rounds[0]
        print(f"Hyperparameter search ({key}): halving over {rounds} epochs, "
              f"keeping 1/{opts['factor']} each round, from {opts['n_iter']} configs"
              + (f" plus {opts['refine_n'] or max(1, opts['n_iter'] // 2)} refined"
                 if opts["refine"] else "")
              + f"; one CV fold each (of {opts['cv']}); keeping the {direction} {opts['metric']}")
    else:
        shape = (f"{opts['n_iter']} parameter sets, one CV fold each (of {opts['cv']})"
                 if opts["params_per_fold"] else
                 f"{opts['n_iter']} parameter sets x {opts['cv']} folds")
        print(f"Hyperparameter search ({key}): {shape}, {opts['epochs']} epochs each; "
              f"keeping the {direction} {opts['metric']}")

    rows = []
    searcher = {"halving": _search_halving, "resample": _search_resample,
                "random": _search_random}[opts["strategy"]]
    best = searcher(ctx, space, opts, rng, rows)

    trials = pd.DataFrame(rows)
    trials = trials.sort_values(["Budget", "Score"], ascending=[False, not opts["higher_is_better"]])
    trials.to_csv(os.path.join(ctx["search_dir"], "trials.csv"), index=False)
    best_params = {k: (list(v) if isinstance(v, (list, tuple))
                       else int(v) if k in ("batch_size", "base_filters", "hidden")
                       else None if v is None
                       else float(v)) for k, v in ((k, best[k]) for k in space)}

    # One fold is a noisy read, so re-train the winner across every fold before saving it
    confirmed = {}
    final_budget = int(best["Budget"])
    if opts["confirm_best"] and ctx["cv"] > 1:
        print(f"\nConfirming {best['Trial']} across all {ctx['cv']} folds at {final_budget} epochs")
        run = {**ctx["base_run"], **best_params, "cv": ctx["cv"], "epochs": final_budget,
               "patience": opts["patience"], "importance": False}
        results = run_models.run_model("best_confirmed", run, ctx["trial_cfg"], samples)
        per_fold = [fold_diagnostics(
            pd.read_csv(os.path.join(ctx["search_dir"], "best_confirmed", str(r["Fold"]),
                                     "history.csv")),
            _fold_pos_rate(os.path.join(ctx["search_dir"], "best_confirmed", str(r["Fold"]))))
            for r in results]
        confirmed = {f"_confirmed_{k}": float(np.nanmean([d[k] for d in per_fold]))
                     for k in ("Val_loss", "Divergence", "Rise", "Val_ROC_AUC")}
        confirmed["_confirmed_score"] = float(np.nanmean([score_trial(d, opts) for d in per_fold]))
        print(f"Confirmed across {len(per_fold)} folds: {opts['metric']} "
              f"{confirmed['_confirmed_score']:.3f} (val loss {confirmed['_confirmed_Val_loss']:.3f}, "
              f"div {confirmed['_confirmed_Divergence']:.3f}, rise {confirmed['_confirmed_Rise']:.3f}, "
              f"val auc {confirmed['_confirmed_Val_ROC_AUC']:.3f})")

    path = save_params(key, best_params, cfg,
                       extra={"_metric": opts["metric"], "_score": float(best["Score"]),
                              "_val_loss": float(best["Val_loss"]),
                              "_divergence": float(best["Divergence"]),
                              "_rise": float(best["Rise"]),
                              "_val_roc_auc": float(best["Val_ROC_AUC"]),
                              "_strategy": opts["strategy"], "_epochs": final_budget,
                              "_trials": int(trials["Trial"].nunique()), "_cv": opts["cv"],
                              **confirmed})
    print(f"\nBest trial {best['Trial']} ({final_budget} epochs): {opts['metric']} "
          f"{best['Score']:.3f} (val loss {best['Val_loss']:.3f}, train {best['Train_loss']:.3f}, "
          f"div {best['Divergence']:.3f}, rise {best['Rise']:.3f}, val auc {best['Val_ROC_AUC']:.3f})")
    print("  " + _show(best_params))
    print(f"Saved to {path} as {key!r} - run_models.py picks these up automatically")
    return best_params, trials


def get_params(key=PARAMS_KEY, base_run=None, cfg=None, force=False, **options):
    """The saved parameters, running the search only if there are none (or force)."""
    cfg = {**CONFIG, **(cfg or {})}
    if not force:
        saved = load_params(key, cfg)
        if saved:
            print(f"Using saved parameters {key!r} from {params_path(cfg)}: " + ", ".join(
                f"{k}={v:.3g}" if isinstance(v, float) else f"{k}={v}" for k, v in saved.items()))
            return saved
    return search(key, base_run, cfg=cfg, force=force, **options)[0]


def main():
    p = argparse.ArgumentParser(description="Random hyperparameter search with CV")
    p.add_argument("--run", help="take the samples/features from this run in RUNS instead of SEARCH_RUN")
    p.add_argument("--key", default=PARAMS_KEY, help="name to save the parameters under")
    p.add_argument("--n-iter", type=int, default=SEARCH_DEFAULTS["n_iter"])
    p.add_argument("--cv", type=int, default=SEARCH_DEFAULTS["cv"])
    p.add_argument("--epochs", type=int, default=SEARCH_DEFAULTS["epochs"])
    p.add_argument("--patience", type=int, default=SEARCH_DEFAULTS["patience"])
    p.add_argument("--metric", default=SEARCH_DEFAULTS["metric"],
                   help="Balanced (default), Performance, Generalization, AUCDiv, LossDiv, "
                        "Overfit, Divergence, Val_loss, Val_ROC_AUC, ...")
    p.add_argument("--maximize", action="store_true",
                   help="higher is better (needed for the AUC/accuracy metrics)")
    p.add_argument("--divergence-weight", type=float, default=SEARCH_DEFAULTS["divergence_weight"],
                   help="weight on the divergence term (keep well below 1)")
    p.add_argument("--underfit-weight", type=float, default=SEARCH_DEFAULTS["underfit_weight"],
                   help="weight on the underfitting penalty in Generalization")
    p.add_argument("--gap-weight", type=float, default=SEARCH_DEFAULTS["gap_weight"])
    p.add_argument("--rise-weight", type=float, default=SEARCH_DEFAULTS["rise_weight"])
    p.add_argument("--strategy", default=SEARCH_DEFAULTS["strategy"],
                   choices=["resample", "halving", "random"],
                   help="resample (default): rounds of fresh models, each round drawn around the "
                        "best so far - every model trained for the same epochs. halving: survivors "
                        "continue from their checkpoints (cheaper, but budget favours slow learners). "
                        "random: one flat round")
    p.add_argument("--configs", nargs="+", type=int, default=None,
                   help="resample strategy: configs per round, e.g. --configs 24 8 3")
    p.add_argument("--round-epochs", nargs="+", type=int, default=None,
                   help="resample strategy: epochs per round, e.g. --round-epochs 10 30 60")
    p.add_argument("--carry-forward", type=int, default=SEARCH_DEFAULTS["carry_forward"],
                   help="best configs retrained at the next round's larger budget")
    p.add_argument("--rounds", nargs="+", type=int, default=SEARCH_DEFAULTS["rounds"],
                   help="cumulative epoch budget per halving round, e.g. --rounds 8 24 60")
    p.add_argument("--factor", type=int, default=SEARCH_DEFAULTS["factor"],
                   help="keep 1/factor of the configs each round")
    p.add_argument("--no-refine", action="store_true",
                   help="skip resampling around the first round's survivors")
    p.add_argument("--refine-n", type=int, default=None, help="how many refined configs")
    p.add_argument("--same-params-all-folds", action="store_true",
                   help="random strategy only: train each set on every fold instead of one")
    p.add_argument("--no-confirm", action="store_true",
                   help="skip re-training the winner across all folds")
    p.add_argument("--out-dir", default=CONFIG["out_dir"])
    p.add_argument("--show", action="store_true", help="print the saved parameters and exit")
    p.add_argument("--force", action="store_true",
                   help="delete the previous search for this key and start over (needed after "
                        "changing the search space)")
    args = p.parse_args()
    cfg = {"out_dir": args.out_dir}
    if args.show:
        path = params_path({**CONFIG, **cfg})
        if not os.path.exists(path):
            print(f"No saved parameters yet ({path}) - run: python src/hypersearch.py")
            return
        with open(path) as f:
            print(json.dumps(json.load(f), indent=2))
        return
    base_run = RUNS[args.run] if args.run else None
    get_params(args.key, base_run=base_run, cfg=cfg, force=args.force, n_iter=args.n_iter,
               cv=args.cv, epochs=args.epochs, patience=args.patience, metric=args.metric,
               higher_is_better=args.maximize, divergence_weight=args.divergence_weight,
               gap_weight=args.gap_weight,
               rise_weight=args.rise_weight, params_per_fold=not args.same_params_all_folds,
               confirm_best=not args.no_confirm, strategy=args.strategy, rounds=args.rounds,
               factor=args.factor, refine=not args.no_refine, refine_n=args.refine_n,
               underfit_weight=args.underfit_weight, carry_forward=args.carry_forward,
               **({"schedule": [[n, e] for n, e in zip(args.configs, args.round_epochs)]}
                  if args.configs and args.round_epochs else {}))


if __name__ == "__main__":
    main()

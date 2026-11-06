#!/usr/bin/env python3

# Feature importance for a trained TCR ResNet, by input gradients.
#
# For a model trained by run_models.py, samples are taken from its predictions and split
# into true positives, true negatives, false positives and false negatives; the gradient
# of the positive-class margin with respect to the input map is averaged over the samples
# in each category. Importance is reported per channel and per region pair (peptide, MHC
# helices, CDR1-3 of each chain), and drawn as 3D stacked bars (region x region, stacked by
# feature group), a saliency map and a per-channel bar chart.
#
#   python src/importance.py results/test_all/single
#   python src/importance.py results/test_all/single --split val --n 16 --method grad
#
#   # in a notebook
#   from importance import run_importance
#   tables = run_importance("results/test_all/single")
#
# Output in <model_dir>/importance/:
#   channel_importance.csv, region_importance.csv, samples.csv   (the table view)
#   importance_3d_<category>.png, saliency_map_<category>.png, channel_importance.png

import argparse
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registers the 3d projection)

from tcr_dataset import FEATURE_GROUPS, TCRMapDataset, load_inputs
from tcr_model import build_resnet

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Regions reported on both axes, in structural order (pMHC side first, then TCR). Each is
# a position channel, optionally intersected with a chain channel; "Other" is every
# remaining residue of the complex (MHC floor, TCR framework).
REGIONS = [
    ("Peptide", "pos_chain_Peptide", None),
    ("MHCha", "pos_region_MHCha", None),
    ("MHChb", "pos_region_MHChb", None),
    ("CDR1a", "pos_region_CDR1", "pos_chain_TCRa"),
    ("CDR2a", "pos_region_CDR2", "pos_chain_TCRa"),
    ("CDR3a", "pos_region_CDR3", "pos_chain_TCRa"),
    ("CDR1b", "pos_region_CDR1", "pos_chain_TCRb"),
    ("CDR2b", "pos_region_CDR2", "pos_chain_TCRb"),
    ("CDR3b", "pos_region_CDR3", "pos_chain_TCRb"),
]
OTHER = "Other"
CATEGORIES = ["TP", "FN", "TN", "FP"]   # positives first, each followed by its error type

# Categorical slots 1-4 of the reference palette, in order (validated for adjacent pairs);
# every chart also ships the CSV table view, which covers the low-contrast slots.
GROUP_COLORS = {"position": "#2a78d6", "sequence": "#eb6834",
                "structure": "#1baf7a", "physics": "#eda100"}
INK = {"primary": "#0b0b0b", "secondary": "#52514e", "muted": "#8a8880", "surface": "#fcfcfb"}
# Sequential: one hue, light -> dark (never a rainbow)
BLUE_RAMP = LinearSegmentedColormap.from_list("blue_ramp", ["#fcfcfb", "#9ec4ef", "#2a78d6", "#123a6b"])


##############################
###### Model and samples
##############################

def load_trained(model_dir, device):
    """Load a model written by run_models.py: (model, info, stats)."""
    with open(os.path.join(model_dir, "model_info.json")) as f:
        info = json.load(f)
    stats = None
    if info.get("normalized"):
        with np.load(os.path.join(model_dir, "channel_stats.npz")) as d:
            stats = (d["mean"], d["std"])
    model = build_resnet(len(info["channels"]), img_size=info["canvas_size"], layers=info["layers"],
                         dropout=info.get("dropout", 0.0),
                         conv_dropout=info.get("conv_dropout", 0.0),
                         hidden=info.get("hidden", 1024),
                         base_filters=info.get("base_filters", 64),
                         pool_size=info.get("pool_size", 4),
                         max_filters=info.get("max_filters"))
    model.load_state_dict(torch.load(os.path.join(model_dir, "best.pt"), map_location=device,
                                     weights_only=True))
    return model.to(device).eval(), info, stats


def pick_samples(model_dir, split="test", n=8, threshold=0.5, seed=0):
    """Highest-confidence n samples per category (TP, FN, TN, FP) from the model's
    predictions. Returns a DataFrame with Name, Label, Prob, Category."""
    combined = os.path.join(model_dir, "predictions.csv")
    if os.path.exists(combined):
        preds = pd.read_csv(combined)
        preds = preds[preds["Set"] == split].reset_index(drop=True)
        if preds.empty:
            raise ValueError(f"{combined} has no {split!r} rows")
    else:   # models written before predictions.csv existed
        preds = pd.read_csv(os.path.join(model_dir, f"{split}_predictions.csv"))
    pred_pos = preds["Prob"] >= threshold
    cat = pd.Series(np.where(preds["Label"] == 1,
                             np.where(pred_pos, "TP", "FN"),
                             np.where(pred_pos, "FP", "TN")), index=preds.index)
    preds = preds.assign(Category=cat)
    # Rank by how strongly the model committed to its call, so each category shows what
    # drove that outcome rather than borderline cases
    conf = np.where(pred_pos, preds["Prob"], 1 - preds["Prob"])
    picked = (preds.assign(_conf=conf).sort_values("_conf", ascending=False)
              .groupby("Category").head(n).drop(columns="_conf"))
    order = {c: i for i, c in enumerate(CATEGORIES)}
    return picked.sort_values(["Category", "Prob"], key=lambda s: s.map(order) if s.name == "Category" else s)


##############################
###### Saliency and regions
##############################

def saliency(model, dataset, idx, device, method="grad_input"):
    """|gradient| (or |gradient x input|) of the positive-class margin for one sample,
    with the padding zeroed. Returns (C, S, S) float32 and the predicted probability."""
    image, _, _ = dataset[idx]
    x = image.unsqueeze(0).to(device).requires_grad_(True)
    out = model(x)
    margin = out[0, 1] - out[0, 0]
    model.zero_grad(set_to_none=True)
    margin.backward()
    grad = x.grad[0].detach()
    sal = (grad * x[0].detach()).abs() if method == "grad_input" else grad.abs()
    sal = sal.cpu().numpy().astype(np.float32)
    # Only the complex itself carries signal; padding attribution is an artifact
    mask = complex_mask(dataset, idx)
    sal[:, ~mask, :] = 0
    sal[:, :, ~mask] = 0
    prob = torch.softmax(out[0].detach(), 0)[1].item()
    return sal, prob


def complex_mask(dataset, idx):
    """Boolean (canvas_size,) mask of the rows/cols holding the complex (not padding)."""
    row = dataset.samples.iloc[idx]
    L = int(row["Length"])
    start = dataset.anchor - int(row["Anchor_index"])
    mask = np.zeros(dataset.canvas_size, dtype=bool)
    mask[start:start + L] = True
    return mask


def region_masks(dataset, idx):
    """Boolean (canvas_size,) mask per region, from the position channels of the stored
    map (read unnormalized, so the chain/region channels keep their 1/2 coding)."""
    x, anchor_index = dataset.load_map(idx, all_channels=True)
    offset = dataset.anchor - anchor_index
    chan = {name: k for k, name in enumerate(dataset.stored_channels)}

    def diag(name, value):
        d = np.diag(x[chan[name]]) == value
        m = np.zeros(dataset.canvas_size, dtype=bool)
        m[offset:offset + len(d)] = d
        return m

    masks = {}
    for name, region_ch, chain_ch in REGIONS:
        m = diag(region_ch, 1 if region_ch.startswith("pos_chain") else 2)
        if chain_ch:
            m &= diag(chain_ch, 1)
        masks[name] = m
    covered = np.logical_or.reduce(list(masks.values()))
    masks[OTHER] = complex_mask(dataset, idx) & ~covered
    return masks


def aggregate(sal, dataset, idx, groups):
    """Mean importance per channel, and per (region, region) block per feature group.
    Block values are means over the cells in the block, so block size does not inflate
    them. Returns (channel_mean (C,), region_table rows)."""
    mask = complex_mask(dataset, idx)
    channel_mean = sal[:, mask, :][:, :, mask].mean(axis=(1, 2))
    masks = region_masks(dataset, idx)
    names = [r[0] for r in REGIONS] + [OTHER]
    group_index = {g: np.array([k for k, gg in enumerate(groups) if gg == g]) for g in set(groups)}
    rows = []
    for r1 in names:
        m1 = masks[r1]
        if not m1.any():
            continue
        block_rows = sal[:, m1, :]
        for r2 in names:
            m2 = masks[r2]
            if not m2.any():
                continue
            block = block_rows[:, :, m2]                  # (C, n1, n2)
            per_channel = block.mean(axis=(1, 2))
            row = {"Region_i": r1, "Region_j": r2, "N_pairs": int(m1.sum() * m2.sum())}
            for g, cols in group_index.items():
                row[g] = float(per_channel[cols].mean()) if len(cols) else 0.0
            rows.append(row)
    return channel_mean, rows


##############################
###### Plots
##############################

def style_axes(ax):
    ax.set_facecolor(INK["surface"])
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(INK["muted"])
        ax.spines[spine].set_linewidth(0.8)
    ax.tick_params(colors=INK["secondary"], labelsize=8, length=3, width=0.8)
    ax.grid(True, axis="y", color=INK["muted"], alpha=0.25, linewidth=0.6)
    ax.set_axisbelow(True)


def plot_region_3d(region_df, groups_present, category, out_file, n_samples):
    """3D bars: region x region, stacked by feature group."""
    names = [r[0] for r in REGIONS] + [OTHER]
    names = [n for n in names if n in set(region_df.Region_i)]
    pos = {n: i for i, n in enumerate(names)}
    fig = plt.figure(figsize=(11, 8.5), facecolor=INK["surface"])
    ax = fig.add_subplot(111, projection="3d")
    ax.set_facecolor(INK["surface"])

    for _, r in region_df.iterrows():
        bottom = 0.0
        for g in groups_present:                      # fixed order -> consistent stacking
            h = float(r[g])
            if h <= 0:
                continue
            ax.bar3d(pos[r.Region_i] + 0.14, pos[r.Region_j] + 0.14, bottom, 0.72, 0.72, h,
                     color=GROUP_COLORS[g], shade=False, edgecolor=INK["surface"], linewidth=0.3)
            bottom += h

    ax.set_xticks([i + 0.5 for i in range(len(names))])
    ax.set_yticks([i + 0.5 for i in range(len(names))])
    ax.set_xticklabels(names, rotation=48, ha="right", fontsize=8, color=INK["secondary"])
    ax.set_yticklabels(names, rotation=-18, ha="left", fontsize=8, color=INK["secondary"])
    ax.set_zlabel("mean |gradient x input|", fontsize=9, color=INK["secondary"], labelpad=8)
    ax.tick_params(colors=INK["secondary"])
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.set_facecolor(INK["surface"])
        axis.pane.set_edgecolor(INK["muted"])
        axis.pane.set_alpha(0.35)
        axis._axinfo["grid"].update(color=INK["muted"], linewidth=0.5, alpha=0.3)
    ax.view_init(elev=26, azim=-121)
    ax.set_title(f"Region-pair importance - {category} ({n_samples} samples)",
                 fontsize=12, color=INK["primary"], pad=16, loc="left")
    handles = [plt.Rectangle((0, 0), 1, 1, color=GROUP_COLORS[g]) for g in groups_present]
    ax.legend(handles, groups_present, loc="upper right", frameon=False, fontsize=9,
              labelcolor=INK["secondary"], title="feature group",
              title_fontproperties={"size": 9}, bbox_to_anchor=(1.02, 0.92))
    fig.tight_layout()
    fig.savefig(out_file, dpi=200, facecolor=INK["surface"])
    plt.close(fig)


def plot_saliency_map(sal_mean, masks, category, out_file, n_samples):
    """Importance summed over channels on the canvas, with the regions marked."""
    total = sal_mean.sum(axis=0)
    rows = np.where(total.sum(axis=1) > 0)[0]
    lo, hi = rows.min(), rows.max() + 1
    fig, ax = plt.subplots(figsize=(7.2, 6.4), facecolor=INK["surface"])
    vmax = np.percentile(total[lo:hi, lo:hi], 99.5) or total.max()
    im = ax.imshow(total[lo:hi, lo:hi], cmap=BLUE_RAMP, vmin=0, vmax=vmax, interpolation="nearest")
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cb.set_label("importance (sum over channels)", fontsize=9, color=INK["secondary"])
    cb.ax.tick_params(colors=INK["secondary"], labelsize=8)
    cb.outline.set_visible(False)
    centers, labels = [], []
    for name, _, _ in REGIONS:
        idx = np.where(masks[name])[0]
        if not len(idx):
            continue
        a, b = idx.min() - lo, idx.max() + 1 - lo
        for line in (a, b):
            ax.axhline(line, color=INK["secondary"], linewidth=0.5, alpha=0.5)
            ax.axvline(line, color=INK["secondary"], linewidth=0.5, alpha=0.5)
        centers.append((a + b) / 2)
        labels.append(name)
    # Ticks rather than text, so the labels get their own margin instead of being clipped
    ax.set_yticks(centers)
    ax.set_yticklabels(labels, fontsize=7, color=INK["secondary"])
    ax.set_xticks(centers)
    ax.set_xticklabels(labels, fontsize=7, color=INK["secondary"], rotation=45, ha="right")
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title(f"Saliency map - {category} ({n_samples} samples)", fontsize=12,
                 color=INK["primary"], pad=12, loc="left")
    fig.tight_layout()
    fig.savefig(out_file, dpi=200, facecolor=INK["surface"])
    plt.close(fig)


def plot_channel_importance(channel_df, out_file, top=25):
    """Per-channel importance, one panel per category, bars colored by feature group."""
    cats = [c for c in CATEGORIES if c in set(channel_df.Category)]
    order = (channel_df.groupby("Channel")["Importance"].mean().sort_values(ascending=False)
             .head(top).index[::-1])
    fig, axes = plt.subplots(1, len(cats), figsize=(3.6 * len(cats), 0.26 * len(order) + 2.2),
                             sharey=True, facecolor=INK["surface"])
    axes = np.atleast_1d(axes)
    for ax, cat in zip(axes, cats):
        sub = channel_df[channel_df.Category == cat].set_index("Channel").loc[order]
        ax.barh(range(len(order)), sub["Importance"],
                color=[GROUP_COLORS[g] for g in sub["Group"]], height=0.72)
        style_axes(ax)
        ax.grid(True, axis="x", color=INK["muted"], alpha=0.25, linewidth=0.6)
        ax.grid(False, axis="y")
        ax.set_yticks(range(len(order)))
        ax.set_yticklabels(order, fontsize=7)
        ax.set_title(cat, fontsize=10, color=INK["primary"], loc="left")
        ax.set_xlabel("mean importance", fontsize=8, color=INK["secondary"])
    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for g, c in GROUP_COLORS.items()
               if g in set(channel_df.Group)]
    labels = [g for g in GROUP_COLORS if g in set(channel_df.Group)]
    axes[-1].legend(handles, labels, loc="lower right", frameon=False, fontsize=8,
                    labelcolor=INK["secondary"])
    fig.suptitle(f"Top {len(order)} channels by importance", fontsize=12, color=INK["primary"],
                 x=0.01, ha="left")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(out_file, dpi=200, facecolor=INK["surface"])
    plt.close(fig)


##############################
###### Run
##############################

def run_importance(model_dir, inputs=None, split="test", n=8, method="grad_input",
                   threshold=0.5, device=None, plots=True):
    """Gradient importance for one trained model. Returns
    {"channels": df, "regions": df, "samples": df}."""
    device = torch.device(device or ("cuda:0" if torch.cuda.is_available() else "cpu"))
    model, info, stats = load_trained(model_dir, device)
    out_dir = os.path.join(model_dir, "importance")
    os.makedirs(out_dir, exist_ok=True)

    picked = pick_samples(model_dir, split, n, threshold)
    inputs = inputs or os.path.join(PROJECT_DIR, "data", "input", "Input_structures.csv")
    all_inputs = load_inputs(inputs)
    # Keep only the prediction columns: predictions.csv also carries the annotation columns,
    # and letting those collide with Input_structures.csv suffixes them to _x/_y, leaving the
    # dataset without the Anchor_index / TCR_side it needs
    keep = [c for c in ("Name", "Label", "Prob", "Pred", "Category") if c in picked]
    samples = picked[keep].merge(all_inputs, on="Name", how="inner")
    missing = set(picked.Name) - set(samples.Name)
    if missing:
        print(f"Warning: {len(missing)} sampled structures are not in {inputs}")
    print(f"{model_dir}: {len(samples)} samples from {split} "
          + ", ".join(f"{c} {int((samples.Category == c).sum())}" for c in CATEGORIES))

    ds = TCRMapDataset(samples, features=info["features"], name_col="Name", label_col="Label",
                       stats=stats, canvas_size=info["canvas_size"], anchor=info["anchor"],
                       pad_value=info["pad_value"])
    if ds.channel_names != info["channels"]:
        raise ValueError("Inputs were built with different channels than the model was trained on")
    groups = [ds.stored_groups[ds.stored_channels.index(c)] for c in ds.channel_names]
    groups_present = [g for g in ["position"] + FEATURE_GROUPS if g in set(groups)]

    channel_rows, region_rows, sal_by_cat, masks_by_cat = [], [], {}, {}
    for cat in CATEGORIES:
        idxs = list(np.where(samples.Category.to_numpy() == cat)[0])
        if not idxs:
            print(f"  {cat}: no samples")
            continue
        sal_sum, chan_sum, blocks = None, None, []
        for i in idxs:
            sal, prob = saliency(model, ds, i, device, method)
            chan_mean, rows = aggregate(sal, ds, i, groups)
            sal_sum = sal if sal_sum is None else sal_sum + sal
            chan_sum = chan_mean if chan_sum is None else chan_sum + chan_mean
            blocks.append(pd.DataFrame(rows))
        sal_by_cat[cat] = sal_sum / len(idxs)
        masks_by_cat[cat] = region_masks(ds, idxs[0])
        channel_rows += [{"Category": cat, "Channel": c, "Group": g, "Importance": v, "N": len(idxs)}
                         for c, g, v in zip(ds.channel_names, groups, chan_sum / len(idxs))]
        region_df = (pd.concat(blocks).groupby(["Region_i", "Region_j"], as_index=False)
                     .mean(numeric_only=True).assign(Category=cat, N=len(idxs)))
        region_rows.append(region_df)
        print(f"  {cat}: {len(idxs)} samples")

    channels = pd.DataFrame(channel_rows)
    regions = pd.concat(region_rows, ignore_index=True) if region_rows else pd.DataFrame()
    channels.to_csv(os.path.join(out_dir, "channel_importance.csv"), index=False)
    regions.to_csv(os.path.join(out_dir, "region_importance.csv"), index=False)
    samples.drop(columns=[c for c in ("Path", "File") if c in samples]).to_csv(
        os.path.join(out_dir, "samples.csv"), index=False)

    if plots and len(channels):
        plot_channel_importance(channels, os.path.join(out_dir, "channel_importance.png"))
        for cat in sal_by_cat:
            sub = regions[regions.Category == cat]
            plot_region_3d(sub, groups_present, cat,
                           os.path.join(out_dir, f"importance_3d_{cat}.png"), int(sub.N.iloc[0]))
            plot_saliency_map(sal_by_cat[cat], masks_by_cat[cat], cat,
                              os.path.join(out_dir, f"saliency_map_{cat}.png"),
                              int(sub.N.iloc[0]))
        print(f"Wrote tables and plots to {out_dir}")
    return {"channels": channels, "regions": regions, "samples": samples}


def main():
    p = argparse.ArgumentParser(description="Gradient feature importance for a trained model")
    p.add_argument("model_dir", help="a run_models.py model folder, e.g. results/test_all/single")
    p.add_argument("--inputs", default=None, help="Input_structures.csv (default: data/input/)")
    p.add_argument("--split", default="test", choices=["test", "val", "train"],
                   help="predictions to sample from")
    p.add_argument("--n", type=int, default=8, help="samples per category")
    p.add_argument("--method", default="grad_input", choices=["grad_input", "grad"])
    p.add_argument("--threshold", type=float, default=0.5)
    p.add_argument("--device", default=None)
    args = p.parse_args()
    run_importance(args.model_dir, args.inputs, args.split, args.n, args.method,
                   args.threshold, args.device)


if __name__ == "__main__":
    main()

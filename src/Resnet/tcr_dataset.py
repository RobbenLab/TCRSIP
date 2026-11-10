#!/usr/bin/env python3

# Dataset and DataLoader for the 2D TCR-pMHC input maps written by build_inputs.py, used by
# every model. Standalone (numpy, pandas, torch only) so it can be dropped into Colab.
#
# Each .npz holds an unpadded (C, L, L) float16 map of the trimmed complex. At load time
# the map is placed on a canvas_size x canvas_size canvas so that the last peptide residue
# lands on the anchor on both axes (the pMHC side fills the canvas up to the anchor, the
# TCR side starts right after it) and the rest is filled with pad_value.
#
# Features: the position channels (chain, region, sequence separation) are always used;
# any combination of the sequence, structure and physics channel groups can be added.
#
#   inputs = load_inputs("data/input/Input_structures.csv")
#   ds = TCRMapDataset(inputs, features=["sequence", "structure"])
#   loader = make_loader(ds, batch_size=16, shuffle=True)
#   images, labels, names = next(iter(loader))    # (image, label, name) as in the notebook

import os

import numpy as np
import pandas as pd
import torch
from torch.utils.data import DataLoader, Dataset

# Canvas: last peptide residue at "anchor" on both axes. With the trimming rules in
# build_inputs.py the pMHC side is <= 204 residues and the TCR side <= 247, so 512 with the
# anchor at 255 fits everything. Padding value as in the notebook.
CANVAS = {"canvas_size": 512, "anchor": 255, "pad_value": -1.0}

# Channel groups: position is always included, the others are selectable
BASE_GROUP = "position"
FEATURE_GROUPS = ["sequence", "structure", "physics"]

# Shared RAM cache: file path -> (all channels float16 (C, L, L), anchor_index)
_CACHE = {}


def load_inputs(path):
    """Read Input_structures.csv (written by build_inputs.py) and resolve each File
    relative to the CSV's folder, so the data folder can be moved (e.g. to Google Drive)
    without editing paths."""
    df = pd.read_csv(path)
    root = os.path.dirname(os.path.abspath(path))
    df["Path"] = [os.path.join(root, f) for f in df["File"]]
    return df


def _read_npz(path):
    with np.load(path) as d:
        return d["x"], int(d["anchor_index"])


def preload_maps(paths):
    """Decompress maps into the shared RAM cache (~17 MB each; skips cached ones)."""
    for p in paths:
        if p not in _CACHE:
            _CACHE[p] = _read_npz(p)


def check_canvas_fit(samples, canvas_size=CANVAS["canvas_size"], anchor=CANVAS["anchor"]):
    """Raise if any sample does not fit the canvas, using the Anchor_index / TCR_side
    columns of Input_structures.csv (no maps are loaded)."""
    over = samples[(samples["Anchor_index"] > anchor) | (samples["TCR_side"] > canvas_size - anchor - 1)]
    if len(over):
        raise ValueError(f"{len(over)} samples do not fit the {canvas_size} canvas anchored at {anchor} "
                         f"(pMHC side <= {anchor + 1}, TCR side <= {canvas_size - anchor - 1}): "
                         f"{list(over.iloc[:, 0][:10])}")


class TCRMapDataset(Dataset):
    """samples: DataFrame with Path (from load_inputs) or File (+ data_dir), Anchor_index
    and TCR_side - i.e. rows of Input_structures.csv, optionally joined with labels.
    features: optional channel groups to add to the always-included position channels,
        any of "sequence", "structure", "physics" (default: all three).
    name_col / label_col: sample name and label columns; without a label column every
        label is -1 (inference).
    stats: optional (mean, std) from channel_stats(), over all stored channels or just the
        selected ones; applied before padding so the padding value is not shifted.
    preload: decompress every map into the shared RAM cache up front."""

    def __init__(self, samples, features=FEATURE_GROUPS, name_col="Name", label_col=None, stats=None,
                 data_dir=None, transform=None, canvas_size=None, anchor=None, pad_value=None,
                 preload=False):
        self.samples = samples.reset_index(drop=True)
        if len(self.samples) == 0:
            raise ValueError("No samples")
        if "Path" in self.samples:
            self.paths = list(self.samples["Path"])
        else:
            self.paths = [os.path.join(data_dir or "", f) for f in self.samples["File"]]
        self.names = list(self.samples[name_col])
        self.targets = (self.samples[label_col].astype(int).to_numpy() if label_col
                        else np.full(len(self.samples), -1))
        self.transform = transform
        self.canvas_size = canvas_size or CANVAS["canvas_size"]
        self.anchor = anchor if anchor is not None else CANVAS["anchor"]
        self.pad_value = pad_value if pad_value is not None else CANVAS["pad_value"]
        check_canvas_fit(self.samples.assign(_name=self.names)[["_name", "Anchor_index", "TCR_side"]],
                         self.canvas_size, self.anchor)

        unknown = [f for f in features if f not in FEATURE_GROUPS]
        if unknown:
            raise ValueError(f"Unknown features {unknown}; choose from {FEATURE_GROUPS} "
                             f"({BASE_GROUP} channels are always included)")
        self.features = [f for f in FEATURE_GROUPS if f in features]
        self.groups = [BASE_GROUP] + self.features

        # Resolve the channels from the group labels stored in the files
        with np.load(self.paths[0]) as d:
            self.stored_channels = list(d["channels"])
            self.stored_groups = list(d["groups"])
        self.index = [k for k, g in enumerate(self.stored_groups) if g in self.groups]
        self.channel_names = [self.stored_channels[k] for k in self.index]

        self.stats = None
        if stats is not None:
            mean, std = (np.asarray(s, dtype=np.float32) for s in stats)
            if len(mean) == len(self.stored_channels):
                mean, std = mean[self.index], std[self.index]
            if len(mean) != len(self.index):
                raise ValueError(f"stats have {len(mean)} channels, expected "
                                 f"{len(self.stored_channels)} or {len(self.index)}")
            std = np.where(std > 0, std, 1.0)
            self.stats = (mean.reshape(-1, 1, 1), std.reshape(-1, 1, 1).astype(np.float32))
        if preload:
            preload_maps(self.paths)

    def __len__(self):
        return len(self.samples)

    def getChannels(self):
        return len(self.index)

    def load_map(self, idx, all_channels=False):
        """Unpadded map as float32 (C, L, L) - selected channels unless all_channels - and
        the anchor_index."""
        p = self.paths[idx]
        x, anchor_index = _CACHE[p] if p in _CACHE else _read_npz(p)
        if len(x) != len(self.stored_channels):
            raise ValueError(f"{p} was built with a different channel list")
        if not all_channels:
            x = x[self.index]
        return x.astype(np.float32), anchor_index

    def place(self, x, anchor_index):
        """Put an unpadded (C, L, L) map on the canvas with the last peptide residue at the
        anchor."""
        C, L, _ = x.shape
        offset = self.anchor - anchor_index
        canvas = np.full((C, self.canvas_size, self.canvas_size), self.pad_value, dtype=np.float32)
        canvas[:, offset:offset + L, offset:offset + L] = x
        return canvas

    def __getitem__(self, idx):
        x, anchor_index = self.load_map(idx)
        if self.stats is not None:
            x = (x - self.stats[0]) / self.stats[1]
        image = torch.from_numpy(self.place(x, anchor_index))
        if self.transform:
            image = self.transform(image)
        return image, int(self.targets[idx]), self.names[idx]


def make_loader(dataset, batch_size=16, shuffle=False, workers=None, pin_memory=None):
    """Standard DataLoader for a TCRMapDataset (map loading/padding runs in worker
    processes)."""
    workers = min(8, os.cpu_count() or 1) if workers is None else workers
    pin_memory = torch.cuda.is_available() if pin_memory is None else pin_memory
    return DataLoader(dataset, batch_size=batch_size, shuffle=shuffle, num_workers=workers,
                      pin_memory=pin_memory, persistent_workers=workers > 0)


def channel_stats(dataset, all_channels=True):
    """Per-channel mean and std over the real (unpadded) residue pairs of a dataset.
    Compute on the training set only and pass as stats= to every dataset. With
    all_channels the stats cover every stored channel, so one set serves any features."""
    n = 0
    s = ss = None
    for i in range(len(dataset)):
        x, _ = dataset.load_map(i, all_channels=all_channels)
        x = x.reshape(x.shape[0], -1)
        if s is None:
            s = np.zeros(len(x))
            ss = np.zeros(len(x))
        n += x.shape[1]
        s += x.sum(1, dtype=np.float64)
        ss += np.square(x, dtype=np.float64).sum(1)
    mean = s / n
    std = np.sqrt(np.maximum(ss / n - mean ** 2, 0))
    return mean.astype(np.float32), std.astype(np.float32)

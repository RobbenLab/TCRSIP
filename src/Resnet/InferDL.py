#!/usr/bin/env python3

# Deploy-time DL inference: builds the 2D input maps (build_inputs.py) for newly
# processed ./dat/<name>/ structures and runs the trained full ResNet checkpoint on them.
#
#   python src/Resnet/InferDL.py --dat-dir ./dat --names-file ./dat/new_structures.txt \
#       --input-dir ./input --checkpoint ./model_checkpoints/full_v1.pt -o ./dl_predictions.csv
#
# --checkpoint is a state_dict (torch.save(model.state_dict(), ...), same as run_models.py's
# best.pt) - model_info.json and channel_stats.npz, written alongside it by the same
# training run (see run_models.py), are expected in the same directory. Same convention as
# importance.py's load_trained, just under whatever filename the checkpoint is shipped as.

import argparse
import json
import os
import subprocess
import sys

import numpy as np
import pandas as pd
import torch

from tcr_dataset import TCRMapDataset, load_inputs, make_loader
from tcr_model import build_resnet

HERE = os.path.dirname(os.path.abspath(__file__))

DAT_DIR = "./dat"
INPUT_DIR = "./input"
CHECKPOINT = "./model_checkpoints/full_v1.pt"


def build_inputs_for(dat_dir, input_dir, names_file, workers):
    cmd = [sys.executable, "build_inputs.py", "-i", os.path.abspath(dat_dir),
          "-o", os.path.abspath(input_dir), "--workers", str(workers)]
    if names_file:
        cmd += ["--names-file", os.path.abspath(names_file)]
    subprocess.run(cmd, cwd=HERE, check=True)


def load_trained(checkpoint_path, device):
    ckpt_dir = os.path.dirname(os.path.abspath(checkpoint_path))
    with open(os.path.join(ckpt_dir, "model_info.json")) as f:
        info = json.load(f)
    stats = None
    if info.get("normalized"):
        with np.load(os.path.join(ckpt_dir, "channel_stats.npz")) as d:
            stats = (d["mean"], d["std"])
    model = build_resnet(len(info["channels"]), img_size=info["canvas_size"], layers=info["layers"],
                         dropout=info.get("dropout", 0.0), conv_dropout=info.get("conv_dropout", 0.0),
                         hidden=info.get("hidden", 1024), base_filters=info.get("base_filters", 64),
                         pool_size=info.get("pool_size", 4), max_filters=info.get("max_filters"))
    model.load_state_dict(torch.load(checkpoint_path, map_location=device, weights_only=True))
    return model.to(device).eval(), info, stats


def predict(checkpoint, input_dir, names_file, batch_size, threshold, device=None):
    device = torch.device(device or ("cuda:0" if torch.cuda.is_available() else "cpu"))
    model, info, stats = load_trained(checkpoint, device)

    samples = load_inputs(os.path.join(input_dir, "Input_structures.csv"))
    if names_file:
        with open(names_file) as f:
            wanted = {line.strip() for line in f if line.strip()}
        samples = samples[samples["Name"].isin(wanted)]
    if samples.empty:
        raise ValueError("No structures to predict on")

    ds = TCRMapDataset(samples, features=info["features"], name_col="Name", label_col=None,
                       stats=stats, canvas_size=info["canvas_size"], anchor=info["anchor"],
                       pad_value=info["pad_value"])
    if ds.channel_names != info["channels"]:
        raise ValueError("Input maps were built with a different channel set than the model was trained on")
    loader = make_loader(ds, batch_size=batch_size, shuffle=False)

    rows = []
    with torch.no_grad():
        for images, _, names in loader:
            probs = torch.softmax(model(images.to(device)), dim=1)[:, 1].cpu().numpy()
            rows += [{"Name": n, "Prob": float(p), "Pred": int(p >= threshold)}
                     for n, p in zip(names, probs)]
    return pd.DataFrame(rows), info


def write_inference_json(dat_dir, preds, model_name, key="DL"):
    """Record each prediction into its structure's <dat_dir>/<name>/<name>.json, under
    Inference.<key> - alongside Scoring/Features, so a structure's folder is the one place
    to look up everything known about it, not just this run's output CSV."""
    missing = []
    for _, row in preds.iterrows():
        name = row["Name"]
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


def main():
    ap = argparse.ArgumentParser(
        description="Build DL input maps for newly processed structures and run the full ResNet checkpoint on them"
    )
    ap.add_argument("--dat-dir", default=DAT_DIR, help="./dat/<name>/ structure+JSON folders (from SetupData.py)")
    ap.add_argument("--names-file", default=None,
                    help="Restrict to these structure names (one per line) - e.g. SetupData.py's "
                         "new_structures.txt - instead of every structure in --dat-dir")
    ap.add_argument("--input-dir", default=INPUT_DIR, help="Where to write the built input maps (default: ./input)")
    ap.add_argument("--checkpoint", default=CHECKPOINT,
                    help="Trained ResNet state_dict, with model_info.json/channel_stats.npz "
                         "alongside it (default: ./model_checkpoints/full_v1.pt)")
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--batch-size", type=int, default=16)
    ap.add_argument("--threshold", type=float, default=0.5)
    ap.add_argument("--device", default=None)
    ap.add_argument("-o", "--output", default="./dl_predictions.csv")
    args = ap.parse_args()

    if not os.path.isfile(args.checkpoint):
        raise FileNotFoundError(
            f"Checkpoint {args.checkpoint} not found. It should be a ResNet state_dict (as "
            "run_models.py writes to best.pt) with model_info.json and channel_stats.npz "
            "alongside it in the same directory."
        )

    os.makedirs(args.input_dir, exist_ok=True)
    build_inputs_for(args.dat_dir, args.input_dir, args.names_file, args.workers)

    preds, info = predict(args.checkpoint, args.input_dir, args.names_file, args.batch_size, args.threshold, args.device)
    preds.to_csv(args.output, index=False)
    model_name = info.get("run") or os.path.basename(args.checkpoint)
    write_inference_json(args.dat_dir, preds, model_name)
    print(f"Saved {len(preds)} predictions to {args.output} and to each structure's Inference.DL")


if __name__ == "__main__":
    main()

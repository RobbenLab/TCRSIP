#!/usr/bin/env python3

# ResNet model, training loop and evaluation for the TCR-pMHC input maps. The model and
# training logic are ported from TCR_Resnet_Full.ipynb. Standalone (torch, numpy, pandas,
# scikit-learn) so it can be dropped into Colab.
#
# train_model() writes into one run directory:
#   checkpoint.pt          model + optimizer + epoch, saved every checkpoint_every epochs;
#                          an interrupted run (e.g. a Colab disconnect) resumes from it
#   best.pt                model weights with the lowest validation loss
#   val_predictions.csv    validation predictions from the best epoch
#   history.csv            per-epoch losses/metrics
# and logs Loss/Accuracy/ROC_AUC per epoch to TensorBoard via the writer passed in.

import os
import time

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from sklearn.metrics import (accuracy_score, average_precision_score, f1_score,
                             precision_score, recall_score, roc_auc_score)


##############################
###### Model (from notebook)
##############################

class ResidualBlock(nn.Module):
    def __init__(self, in_channels, out_channels, stride=1, downsample=None):
        super(ResidualBlock, self).__init__()
        self.conv1 = nn.Sequential(
            nn.Conv2d(in_channels, out_channels, kernel_size=3, stride=stride, padding=1),
            nn.BatchNorm2d(out_channels),
            nn.ReLU())
        self.bn1 = nn.BatchNorm2d(out_channels)
        self.conv2 = nn.Sequential(
            nn.Conv2d(out_channels, out_channels, kernel_size=3, stride=1, padding=1),
            nn.BatchNorm2d(out_channels))
        self.bn2 = nn.BatchNorm2d(out_channels)
        self.downsample = downsample
        self.relu = nn.ReLU()
        self.out_channels = out_channels

    def forward(self, x):
        residual = x
        out = self.conv1(x)
        out = self.bn1(out)
        out = self.conv2(out)
        out = self.bn2(out)
        if self.downsample:
            residual = self.downsample(x)
        out += residual
        out = self.relu(out)
        return out


class ResNet(nn.Module):
    def __init__(self, block, layers, channels, num_classes=2, img_size=512, hidden=1024,
                 dropout=0.0, conv_dropout=0.0, base_filters=64, pool_size=4, max_filters=None):
        """layers: blocks per residual stage, one entry per stage - its length sets the
        depth. Stage i has base_filters * 2**i channels and halves the map (the first keeps
        its size), so the flattened width is derived from the stage count rather than fixed.
        pool_size: adaptive average pool to pool_size x pool_size before the classifier, so
        the fully connected layer stays the same size whatever the depth or input size. It
        is capped at the feature map size, since pooling to more cells than the map has
        would just duplicate them. pool_size=None instead uses the notebook's
        AvgPool2d(7, stride=1).
        max_filters: ceiling on the channel doubling. Without it a deep net is mostly a
        parameter-count experiment (7 stages reaches 4096 channels and 783M parameters);
        with a cap, depth can be varied on its own."""
        super(ResNet, self).__init__()
        self.hidden = hidden
        self.in_filter = base_filters
        self.conv1 = nn.Sequential(
            nn.Conv2d(channels, base_filters, kernel_size=7, stride=2, padding=3),
            nn.BatchNorm2d(base_filters),
            nn.ReLU())
        self.max_pool = nn.MaxPool2d(kernel_size=3, stride=2, padding=1)
        # One stage per entry of `layers`; the first keeps the resolution, the rest halve it
        def stage_filters(i):
            f = base_filters * 2 ** i
            return min(f, max_filters) if max_filters else f

        self.stages = nn.Sequential(*[
            self.make_layer(block, stage_filters(i), n, stride=1 if i == 0 else 2)
            for i, n in enumerate(layers)])
        out_filters = stage_filters(len(layers) - 1)
        # conv1 and max_pool each halve the map, then every stage after the first
        feat = img_size // (4 * 2 ** (len(layers) - 1))
        if pool_size:
            # Pooling to more cells than the map holds would just repeat them
            side = min(int(pool_size), max(feat, 1))
            self.avg_pool = nn.AdaptiveAvgPool2d(side)
        else:
            self.avg_pool = nn.AvgPool2d(7, stride=1)   # the notebook's pooling
            side = feat - 6
        # Below a 4x4 map the pooling has no spatial structure left to summarize
        if feat < 4:
            raise ValueError(f"{len(layers)} stages leave a {feat}x{feat} map for a {img_size} "
                             f"input; each stage after the first halves it, so a {img_size} "
                             f"input allows {int(np.log2(img_size // 4 // 4)) + 1} stages at most. "
                             f"Use fewer stages or a larger input.")
        if side < 1:
            raise ValueError(f"pooling leaves a {side}x{side} map; use a smaller pool_size")
        self.fc1 = nn.Linear(out_filters * side * side, self.hidden)
        # The notebook built this layer but never called it in forward; it is applied here
        # on the hidden layer, where almost all the parameters are (512*10*10 -> 1024)
        self.drop = nn.Dropout(dropout)
        # Optional 2D dropout on the convolutional stack, before the flatten
        self.conv_drop = nn.Dropout2d(conv_dropout)
        self.act1 = nn.ReLU()
        self.fc2 = nn.Linear(self.hidden, num_classes)

    def make_layer(self, block, filters, blocks, stride=1):
        downsample = None
        if stride != 1 or self.in_filter != filters:
            downsample = nn.Sequential(
                nn.Conv2d(self.in_filter, filters, kernel_size=1, stride=stride),
                nn.BatchNorm2d(filters),
            )
        layers = []
        layers.append(block(self.in_filter, filters, stride, downsample))
        self.in_filter = filters
        for i in range(1, blocks):
            layers.append(block(self.in_filter, filters))
        return nn.Sequential(*layers)

    def forward(self, x):
        x = self.conv1(x)
        x = self.max_pool(x)
        x = self.stages(x)
        x = self.avg_pool(x)
        x = self.conv_drop(x)
        x = x.view(x.size(0), -1)  # Flatten layer along batches
        x = self.fc1(x)
        x = self.act1(x)
        x = self.drop(x)
        x = self.fc2(x)
        return x


def build_resnet(channels, img_size=512, layers=(3, 4, 6, 3), dropout=0.0, conv_dropout=0.0,
                 hidden=1024, base_filters=64, pool_size=4, max_filters=None):
    return ResNet(ResidualBlock, list(layers), channels, num_classes=2, img_size=img_size,
                  hidden=hidden, dropout=dropout, conv_dropout=conv_dropout,
                  base_filters=base_filters, pool_size=pool_size, max_filters=max_filters)


##############################
###### Evaluation
##############################

def compute_metrics(labels, probs, threshold=0.5):
    labels = np.asarray(labels).astype(int)
    probs = np.asarray(probs, dtype=float)
    pred = (probs >= threshold).astype(int)
    both = len(np.unique(labels)) == 2
    return {
        "N": len(labels),
        "N_pos": int(labels.sum()),
        "Accuracy": accuracy_score(labels, pred),
        "ROC_AUC": roc_auc_score(labels, probs) if both else np.nan,
        "PR_AUC": average_precision_score(labels, probs) if both else np.nan,
        "F1": f1_score(labels, pred, zero_division=0),
        "Precision": precision_score(labels, pred, zero_division=0),
        "Recall": recall_score(labels, pred, zero_division=0),
    }


@torch.no_grad()
def predict(model, loader, device, criterion=None, amp=False):
    """Run the model over a loader. Returns (DataFrame Name/Label/Prob, mean loss)."""
    model.eval()
    names, labels, probs, losses = [], [], [], []
    for images, y, name in loader:
        images = images.to(device, non_blocking=True)
        y = y.to(device, non_blocking=True)
        with torch.autocast(device.type, dtype=torch.float16, enabled=amp):
            out = model(images)
        out = out.float()
        if criterion is not None:
            losses.append(criterion(out, y).item() * len(y))
        probs += torch.softmax(out, dim=1)[:, 1].cpu().tolist()
        labels += y.cpu().tolist()
        names += list(name)
    df = pd.DataFrame({"Name": names, "Label": labels, "Prob": probs})
    return df, (sum(losses) / len(df) if losses else np.nan)


##############################
###### Training (from notebook)
##############################

def train_model(model, optimizer, criterion, num_epochs, train_loader, valid_loader, device,
                run_dir, test_loader=None, writers=None, tag_prefix="", amp=False, patience=20,
                min_delta=0.0, checkpoint_every=5, log=print):
    """Train with early stopping on validation loss, checkpointing and resume. Training
    stops after `patience` epochs with no improvement of more than `min_delta`; the best
    weights are in run_dir/best.pt. Returns the history DataFrame.
    writers: {"train": SummaryWriter, "val": SummaryWriter} writing to sibling folders, so
    the two curves land on one TensorBoard card; tag_prefix (the run name) becomes the
    collapsible section the cards sit in.
    test_loader: optional, scored every epoch purely for the curve - model selection and
    early stopping always use the validation loss, never the test set."""
    os.makedirs(run_dir, exist_ok=True)
    checkfile = os.path.join(run_dir, "checkpoint.pt")
    bestfile = os.path.join(run_dir, "best.pt")
    scaler = torch.amp.GradScaler("cuda", enabled=amp)
    history = []
    start_epoch = 0
    best = {"epoch": -1, "loss": float("inf")}

    if os.path.exists(checkfile):
        ck = torch.load(checkfile, map_location=device, weights_only=False)
        model.load_state_dict(ck["model_state_dict"])
        optimizer.load_state_dict(ck["optimizer_state_dict"])
        start_epoch = ck["epoch"] + 1
        best = ck["best"]
        history = ck["history"]
        log(f"Resuming from checkpoint at epoch {start_epoch}")
        if ck.get("finished"):
            log("Run already finished")
            return pd.DataFrame(history)

    stopped = False
    for epoch in range(start_epoch, num_epochs):
        model.train()
        start_time = time.time()
        losses, names, labels, probs = [], [], [], []
        for images, y, name in train_loader:
            images = images.to(device, non_blocking=True)
            y = y.to(device, non_blocking=True)
            with torch.autocast(device.type, dtype=torch.float16, enabled=amp):
                outputs = model(images)
                loss = criterion(outputs.float(), y)
            optimizer.zero_grad()
            scaler.scale(loss).backward()
            scaler.step(optimizer)
            scaler.update()
            losses.append(loss.item() * len(y))
            probs += torch.softmax(outputs.detach().float(), dim=1)[:, 1].cpu().tolist()
            labels += y.cpu().tolist()
        train_loss = sum(losses) / len(labels)
        train_m = compute_metrics(labels, probs)
        t = time.time() - start_time

        # Validation (the notebook never switched to eval mode; predict() does)
        val_df, val_loss = predict(model, valid_loader, device, criterion, amp)
        val_m = compute_metrics(val_df["Label"], val_df["Prob"])

        # Scored for the curve only; never used for selection or stopping
        test_loss, test_m = np.nan, {}
        if test_loader is not None:
            test_df, test_loss = predict(model, test_loader, device, criterion, amp)
            test_m = compute_metrics(test_df["Label"], test_df["Prob"])

        row = {"Epoch": epoch, "Time": t, "Train_loss": train_loss, "Val_loss": val_loss,
               **{f"Train_{k}": v for k, v in train_m.items() if k not in ("N", "N_pos")},
               **{f"Val_{k}": v for k, v in val_m.items() if k not in ("N", "N_pos")}}
        if test_m:
            row["Test_loss"] = test_loss
            row.update({f"Test_{k}": v for k, v in test_m.items() if k not in ("N", "N_pos")})
        history.append(row)
        if writers:
            pre = f"{tag_prefix}/" if tag_prefix else ""
            splits = [("train", train_loss, train_m), ("val", val_loss, val_m)]
            if test_m:
                splits.append(("test", test_loss, test_m))
            for split, loss, m in splits:
                w = writers.get(split)
                if w is None:
                    continue
                w.add_scalar(f"{pre}Loss", loss, epoch)
                w.add_scalar(f"{pre}Accuracy", m["Accuracy"], epoch)
                w.add_scalar(f"{pre}ROC_AUC", m["ROC_AUC"], epoch)
                w.add_scalar(f"{pre}PR_AUC", m["PR_AUC"], epoch)
            writers["train"].add_scalar(f"{pre}Epoch_time_s", t, epoch)

        # Save best (lowest validation loss)
        if val_loss < best["loss"] - min_delta:
            best = {"epoch": epoch, "loss": val_loss, **{k: val_m[k] for k in ("Accuracy", "ROC_AUC")}}
            torch.save(model.state_dict(), bestfile)
            # Snapshot of the best epoch; train_one later writes the full predictions.csv
            snapshot = val_df.assign(Pred=(val_df["Prob"] >= 0.5).astype(int))
            snapshot.to_csv(os.path.join(run_dir, "val_predictions.csv"), index=False)
        elif (epoch - best["epoch"]) >= patience:
            stopped = True

        log("Epoch [{}/{}] {:.0f}s  train loss {:.4f} acc {:.3f} auc {:.3f} | val loss {:.4f} "
            "acc {:.3f} auc {:.3f}{}{}".format(
                epoch + 1, num_epochs, t, train_loss, train_m["Accuracy"], train_m["ROC_AUC"],
                val_loss, val_m["Accuracy"], val_m["ROC_AUC"],
                "" if not test_m else " | test auc {:.3f}".format(test_m["ROC_AUC"]),
                "  *" if best["epoch"] == epoch else ""))
        if stopped:
            log(f"Early stopping: no validation loss improvement for {patience} epochs")

        finished = stopped or epoch == num_epochs - 1
        if finished or (epoch + 1) % checkpoint_every == 0:
            torch.save({"epoch": epoch, "model_state_dict": model.state_dict(),
                        "optimizer_state_dict": optimizer.state_dict(), "best": best,
                        "history": history, "finished": finished,
                        # stopped = early stopping, i.e. training more epochs would not help;
                        # finished without stopped just means the epoch budget ran out
                        "stopped": stopped}, checkfile)
        if stopped:
            break

    pd.DataFrame(history).to_csv(os.path.join(run_dir, "history.csv"), index=False)
    log("Best epoch {}: val loss {:.4f}, acc {:.3f}, auc {:.3f}".format(
        best["epoch"] + 1, best["loss"], best.get("Accuracy", np.nan), best.get("ROC_AUC", np.nan)))
    return pd.DataFrame(history)

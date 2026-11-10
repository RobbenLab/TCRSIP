#!/usr/bin/env python3

# Toy run of the ResNet from example/TCR_Resnet_Full.ipynb, used to check that the
# model trains on this server with the notebook's default settings before we swap in
# the restructured TCR inputs. The model and training loop are ported from the
# notebook; the data is a binary image task (CIFAR-10 cat vs dog, upsampled to the
# notebook's 1024x1024 input size) or a synthetic N-channel task that mimics the
# encoded contact-map tensors (BLOSUM62 encoding = 20 + 5 chain + 1 distance = 26).
#
#   conda activate tcr_dl
#   python src/toy_resnet.py                              # defaults from the notebook
#   python src/toy_resnet.py --data synthetic --channels 26
#   python src/toy_resnet.py --img-size 224 --epochs 3    # quick smoke test

import argparse
import os
import time

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, Dataset, Subset
from torch.utils.tensorboard import SummaryWriter
from torchvision import datasets, transforms

from tcr_model import ResidualBlock, ResNet

PROJECT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


##############################
###### Data
##############################

class BinaryCIFAR(Dataset):
    """CIFAR-10 restricted to two classes, returning (image, label, name) like the
    notebook's CustomImageDataset so the training loop is unchanged."""

    def __init__(self, root, train, img_size, n, classes=(3, 5), seed=42):
        tf = transforms.Compose([
            transforms.Resize(img_size),
            transforms.ToTensor(),
            transforms.Normalize((0.4914, 0.4822, 0.4465), (0.247, 0.243, 0.261)),
        ])
        base = datasets.CIFAR10(root, train=train, download=True, transform=tf)
        targets = np.array(base.targets)
        idx = np.where(np.isin(targets, classes))[0]
        rng = np.random.default_rng(seed)
        idx = rng.permutation(idx)[:n]
        self.data = Subset(base, idx)
        self.idx = idx
        self.pos_class = classes[1]
        self.channels = 3

    def __len__(self):
        return len(self.data)

    def getChannels(self):
        return self.channels

    def __getitem__(self, i):
        image, target = self.data[i]
        label = int(target == self.pos_class)
        return image, label, f"cifar_{self.idx[i]}"


class SyntheticMaps(Dataset):
    """Random N-channel square 'contact maps'; positives get a brighter block in a
    random location so there is something learnable. Padded with -1 like makeEncoded."""

    def __init__(self, n, channels, img_size, seed=42):
        self.n = n
        self.channels = channels
        self.img_size = img_size
        self.seed = seed
        self.labels = np.random.default_rng(seed).integers(0, 2, n)

    def __len__(self):
        return self.n

    def getChannels(self):
        return self.channels

    def __getitem__(self, i):
        rng = np.random.default_rng(self.seed + i + 1)
        s = self.img_size
        L = int(rng.integers(s // 2, s))  # "sequence length" before padding
        image = np.full((self.channels, s, s), -1, dtype=np.float32)
        image[:, :L, :L] = rng.normal(0, 1, (self.channels, L, L))
        label = int(self.labels[i])
        if label:
            b = max(4, s // 16)
            x, y = rng.integers(0, L - b, 2)
            image[:, x:x + b, y:y + b] += 1.0
        return torch.from_numpy(image), label, f"synthetic_{i}"


##############################
###### Training (from notebook)
##############################

def trainModel(model, optimizer, criterion, num_epochs, train_loader, valid_loader, writer,
               device, saveDir=None, checkDir=None, amp=False):
    m = nn.Softmax(dim=1)
    best_epoch = 0
    best_test_acc = 0
    best_test_loss = 100000
    start_epoch = 0
    early_stop_thresh = num_epochs * 0.1  # Stop if val loss does not improve after 10% of epochs
    checkpoint_every = max(1, num_epochs // 10)
    scaler = torch.amp.GradScaler("cuda", enabled=amp)
    print(f"Running model for {num_epochs} epochs")

    if checkDir is not None:
        checkfile = os.path.join(checkDir, "checkpoint.pt")
        if os.path.exists(checkfile):
            checkpoint = torch.load(checkfile, weights_only=True)
            model.load_state_dict(checkpoint['model_state_dict'])
            start_epoch = checkpoint['epoch']
            print("Model Checkpointed!")

    print(f"Starting at Epoch: {start_epoch}")
    for epoch in range(start_epoch, num_epochs):
        model.train()
        if device.type == "cuda":
            torch.cuda.reset_peak_memory_stats(device)
        start_time = time.time()
        correct = 0
        total = 0
        losses = []
        for i, (images, labels, name) in enumerate(train_loader):
            images = images.to(device, non_blocking=True)
            labels = labels.to(device, non_blocking=True)

            with torch.autocast(device.type, dtype=torch.float16, enabled=amp):
                outputs = model(images)
                loss = criterion(outputs, labels)
            losses.append(loss.item())
            _, predicted = torch.max(outputs.data, 1)
            total += labels.size(0)
            correct += (predicted == labels).sum().item()

            optimizer.zero_grad()
            scaler.scale(loss).backward()
            scaler.step(optimizer)
            scaler.update()

        trloss = np.mean(losses)
        tracc = 100 * correct / total
        writer.add_scalar("Loss/train", trloss, epoch)
        writer.add_scalar("Accuracy/train", tracc, epoch)
        t = time.time() - start_time

        ############### Checkpointing
        if checkDir is not None and epoch > 0 and epoch % checkpoint_every == 0:
            torch.save({
                'epoch': epoch,
                'model_state_dict': model.state_dict(),
                'optimizer_state_dict': optimizer.state_dict(),
                'loss': trloss,
                'acc': tracc,
            }, checkfile)

        ################### Validation
        # The notebook never switched to eval mode, so BatchNorm used batch statistics
        # during validation; eval() here gives the standard behaviour.
        model.eval()
        with torch.no_grad(), torch.autocast(device.type, dtype=torch.float16, enabled=amp):
            correct = 0
            total = 0
            losses = []
            names = []
            probs = []
            labs = []
            for images, labels, name in valid_loader:
                images = images.to(device, non_blocking=True)
                labels = labels.to(device, non_blocking=True)
                outputs = model(images)
                probability = m(outputs.float())
                loss = criterion(outputs, labels)
                losses.append(loss.item())
                _, predicted = torch.max(outputs.data, 1)
                total += labels.size(0)
                correct += (predicted == labels).sum().item()
                probs += probability[:, 1].cpu().tolist()
                labs += labels.cpu().tolist()
                names += list(name)

            tsloss = np.mean(losses)
            tsacc = 100 * correct / total
            writer.add_scalar("Loss/test", tsloss, epoch)
            writer.add_scalar("Accuracy/test", tsacc, epoch)

        ########### Save best
        if tsloss < best_test_loss:
            best_test_loss = tsloss
            best_test_acc = tsacc
            best_epoch = epoch
            if saveDir is not None:
                torch.save(model.state_dict(), os.path.join(saveDir, "best.pt"))
                results = pd.DataFrame({'Names': names, 'Predictions': probs, 'Labels': labs})
                results.to_csv(os.path.join(saveDir, "best_results.csv"), index=False)
        elif (epoch - best_epoch) > early_stop_thresh:
            print("Early Stopping at epoch {}".format(epoch))
            break

        mem = ""
        if device.type == "cuda":
            mem = ", Peak GPU mem: {:.1f} GB".format(torch.cuda.max_memory_allocated(device) / 1e9)
        print('Epoch [{}/{}], Time: {:.2f} s, Training Loss: {:.4f}, Training Accuracy: {:.2f}, '
              'Validation Loss: {:.4f}, Validation Accuracy: {:.2f} %{}'.format(
                  epoch + 1, num_epochs, t, trloss, tracc, tsloss, tsacc, mem))

    print('Finished model, Best Epoch: {}, Best Validation Loss: {:.4f}, Best Validation Accuracy: {:.2f}'.format(
        best_epoch, best_test_loss, best_test_acc))


##############################
###### Main
##############################

def main():
    p = argparse.ArgumentParser(description="Toy run of the notebook ResNet on a binary image task")
    p.add_argument("--data", choices=["cifar", "synthetic"], default="cifar")
    p.add_argument("--channels", type=int, default=26, help="input channels (synthetic only)")
    p.add_argument("--img-size", type=int, default=1024)
    p.add_argument("--n-train", type=int, default=512)
    p.add_argument("--n-val", type=int, default=128)
    # Defaults below match the notebook's parameter cell
    p.add_argument("--epochs", type=int, default=100)
    p.add_argument("--batch-size", type=int, default=16)
    p.add_argument("--lr", type=float, default=1e-5)
    p.add_argument("--decay", type=float, default=0.1)
    p.add_argument("--layers", type=int, nargs=4, default=[3, 4, 6, 3])
    p.add_argument("--amp", action="store_true", help="mixed precision (fp16) training")
    p.add_argument("--device", default="cuda:0" if torch.cuda.is_available() else "cpu")
    p.add_argument("--workers", type=int, default=8)
    p.add_argument("--data-dir", default=os.path.join(PROJECT_DIR, "data", "toy"))
    p.add_argument("--out-dir", default=os.path.join(PROJECT_DIR, "model_checkpoints", "toy_resnet"))
    p.add_argument("--seed", type=int, default=42)
    args = p.parse_args()

    torch.manual_seed(args.seed)
    device = torch.device(args.device)
    print(f"torch {torch.__version__}, device: {device}"
          + (f" ({torch.cuda.get_device_name(device)})" if device.type == "cuda" else ""))

    if args.data == "cifar":
        training_data = BinaryCIFAR(args.data_dir, True, args.img_size, args.n_train, seed=args.seed)
        test_data = BinaryCIFAR(args.data_dir, False, args.img_size, args.n_val, seed=args.seed)
    else:
        training_data = SyntheticMaps(args.n_train, args.channels, args.img_size, seed=args.seed)
        test_data = SyntheticMaps(args.n_val, args.channels, args.img_size, seed=args.seed + 10**6)

    loader_kw = dict(batch_size=args.batch_size, num_workers=args.workers,
                     pin_memory=device.type == "cuda")
    train_loader = DataLoader(training_data, shuffle=True, **loader_kw)
    valid_loader = DataLoader(test_data, shuffle=False, **loader_kw)
    channel_num = training_data.getChannels()

    model = ResNet(ResidualBlock, args.layers, channel_num, num_classes=2,
                   img_size=args.img_size).to(device)
    print(f"Input: {channel_num} x {args.img_size} x {args.img_size}, "
          f"parameters: {sum(p.numel() for p in model.parameters()) / 1e6:.1f}M")

    criterion = nn.CrossEntropyLoss()
    optimizer = optim.Adam(model.parameters(), lr=args.lr, weight_decay=args.decay)

    run_name = f"{args.data}_{args.img_size}px_b{args.batch_size}"
    save_dir = os.path.join(args.out_dir, run_name)
    os.makedirs(save_dir, exist_ok=True)
    writer = SummaryWriter(os.path.join(args.out_dir, "runs", run_name))
    trainModel(model, optimizer, criterion, args.epochs, train_loader, valid_loader, writer,
               device, saveDir=save_dir, checkDir=None, amp=args.amp)
    writer.close()


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
02_run_g4mer.py — score mRNA region sequences for rG4 formation with G4mer.

G4mer (Biociphers/g4mer) is a 6-mer tokenized RNA language model fine-tuned from
mRNAbert. It accepts a MAXIMUM OF 70 nt per input. The authors' recommendation for
longer sequences is to scan with a 70-nt sliding window and take the maximum rG4
score across windows, which is what this script does.

Per-window scores are written in full alongside the collapsed per-transcript
summary, so alternative aggregations (mean, fraction above threshold, positional)
can be explored later without re-running inference.

Same script runs locally (5'UTR pilot) and on the cluster (CDS / 3'UTR) — only the
--region and --device flags change.

Usage:
    python g4mer/02_run_g4mer.py --region utr5 --selftest
    python g4mer/02_run_g4mer.py --region cds --device cuda --batch-size 512
"""

import argparse
import gzip
import os
import sys
import time

import torch
from transformers import AutoTokenizer, AutoModelForSequenceClassification

WINDOW_DEFAULT = 70  # G4mer hard limit: max 70 nt per input sequence
KMER = 6            # G4mer tokenization is 6-mer

# Model card example — a known G-rich rG4-forming sequence. Used by --selftest.
SELFTEST_SEQ = "GGGAGGGCGCGTGTGGTGAGAGGAGGGAGGGAAGGAAGGCGGAGGAAGGA"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--region", required=True, choices=["utr5", "cds", "utr3"],
                   help="mRNA region to score; selects input FASTA and output names")
    p.add_argument("--fasta", default=None,
                   help="Input FASTA (default: <basedir>/fasta/<region>_sequences.fa)")
    p.add_argument("--outdir", default=None,
                   help="Output directory (default: <repo>/output/g4mer)")
    p.add_argument("--basedir", default=None,
                   help="Base dir holding fasta/ (default: directory containing this script)")
    p.add_argument("--model", default="Biociphers/g4mer",
                   help="HuggingFace model id or local path")
    p.add_argument("--window", type=int, default=WINDOW_DEFAULT,
                   help="Sliding window size in nt (G4mer max is 70)")
    p.add_argument("--stride", type=int, default=10,
                   help="Sliding window stride in nt")
    p.add_argument("--batch-size", type=int, default=256)
    p.add_argument("--threshold", type=float, default=0.5,
                   help="Probability cutoff for the frac_above summary statistic")
    p.add_argument("--device", default="auto", choices=["auto", "cuda", "mps", "cpu"])
    p.add_argument("--limit", type=int, default=None,
                   help="Score only the first N transcripts (pilot runs)")
    p.add_argument("--selftest", action="store_true",
                   help="Score the model card example sequence and assert it reads as rG4")
    return p.parse_args()


def pick_device(requested):
    if requested != "auto":
        return torch.device(requested)
    if torch.cuda.is_available():
        return torch.device("cuda")
    if torch.backends.mps.is_available():
        return torch.device("mps")
    return torch.device("cpu")


def read_fasta(path):
    """Yield (name, sequence). Handles the line-wrapped FASTA that writeXStringSet emits."""
    name, chunks = None, []
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.upper())
    if name is not None:
        yield name, "".join(chunks)


def windows(seq, window, stride):
    """Tile a sequence into (start, end, subseq), 0-based half-open starts.

    Sequences at or below the window length are scored whole in a single pass.
    The final window is always anchored to the 3' end so the tail is never dropped
    when (len - window) is not a multiple of stride.
    """
    n = len(seq)
    if n <= window:
        return [(0, n, seq)]
    starts = list(range(0, n - window + 1, stride))
    if starts[-1] != n - window:
        starts.append(n - window)
    return [(s, s + window, seq[s:s + window]) for s in starts]


def to_kmers(seq, k=KMER):
    return " ".join(seq[i:i + k] for i in range(len(seq) - k + 1))


class Scorer:
    def __init__(self, model_id, device, batch_size):
        self.device = device
        self.batch_size = batch_size
        self.tokenizer = AutoTokenizer.from_pretrained(model_id)
        self.model = AutoModelForSequenceClassification.from_pretrained(model_id)
        self.model.to(device)
        self.model.eval()

    @torch.no_grad()
    def score(self, seqs):
        """Return P(rG4) for a list of raw nucleotide strings."""
        out = []
        for i in range(0, len(seqs), self.batch_size):
            batch = seqs[i:i + self.batch_size]
            enc = self.tokenizer([to_kmers(s) for s in batch],
                                 return_tensors="pt", padding=True, truncation=True)
            enc = {k: v.to(self.device) for k, v in enc.items()}
            logits = self.model(**enc).logits
            probs = torch.softmax(logits.float(), dim=1)[:, 1]
            out.extend(probs.detach().cpu().tolist())
        return out


def main():
    args = parse_args()

    basedir = args.basedir or os.path.dirname(os.path.abspath(__file__))
    repo_root = os.path.dirname(basedir)
    fasta_path = args.fasta or os.path.join(basedir, "fasta", f"{args.region}_sequences.fa")
    outdir = args.outdir or os.path.join(repo_root, "output", "g4mer")
    os.makedirs(outdir, exist_ok=True)

    if args.window > WINDOW_DEFAULT:
        sys.exit(f"ERROR: --window {args.window} exceeds the G4mer 70 nt input limit.")
    if not os.path.exists(fasta_path):
        sys.exit(f"ERROR: FASTA not found: {fasta_path}\nRun g4mer/01_generate_fasta.R first.")

    device = pick_device(args.device)
    print(f"Device: {device}", flush=True)
    print(f"Model:  {args.model}", flush=True)
    print(f"Input:  {fasta_path}", flush=True)

    scorer = Scorer(args.model, device, args.batch_size)

    if args.selftest:
        p = scorer.score([SELFTEST_SEQ])[0]
        print(f"Self-test (model card rG4 example): P(rG4) = {p:.4f}", flush=True)
        assert 0.0 <= p <= 1.0, "probability out of range"
        assert p > 0.5, (
            f"Known rG4-forming control scored {p:.4f} (expected > 0.5). "
            "Check that the model id and 6-mer tokenization are correct."
        )
        print("Self-test PASSED", flush=True)

    win_path = os.path.join(outdir, f"g4mer_windows_{args.region}.tsv.gz")
    sum_path = os.path.join(outdir, f"g4mer_summary_{args.region}.tsv")

    n_tx = 0
    n_win_total = 0
    t0 = time.time()

    with gzip.open(win_path, "wt") as win_fh, open(sum_path, "w") as sum_fh:
        win_fh.write("transcript_id_clean\tregion\twindow_start\twindow_end\tg4_prob\n")
        sum_fh.write(
            "transcript_id_clean\tregion\tregion_length\tn_windows\t"
            "g4mer_max\tg4mer_mean\tg4mer_frac_above\tg4mer_argmax_start\n"
        )

        for name, seq in read_fasta(fasta_path):
            if args.limit is not None and n_tx >= args.limit:
                break
            if len(seq) < KMER:
                continue

            wins = windows(seq, args.window, args.stride)
            probs = scorer.score([w[2] for w in wins])

            for (start, end, _), p in zip(wins, probs):
                win_fh.write(f"{name}\t{args.region}\t{start}\t{end}\t{p:.6f}\n")

            n = len(probs)
            g_max = max(probs)
            g_mean = sum(probs) / n
            frac_above = sum(1 for p in probs if p >= args.threshold) / n
            argmax_start = wins[probs.index(g_max)][0]

            sum_fh.write(
                f"{name}\t{args.region}\t{len(seq)}\t{n}\t"
                f"{g_max:.6f}\t{g_mean:.6f}\t{frac_above:.6f}\t{argmax_start}\n"
            )

            n_tx += 1
            n_win_total += n
            if n_tx % 1000 == 0:
                rate = n_win_total / max(time.time() - t0, 1e-9)
                print(f"  {n_tx} transcripts, {n_win_total} windows "
                      f"({rate:.0f} win/s)", flush=True)

    elapsed = time.time() - t0
    print(f"\nDone: {n_tx} transcripts, {n_win_total} windows in {elapsed/60:.1f} min",
          flush=True)
    print(f"Per-window scores: {win_path}", flush=True)
    print(f"Summary:           {sum_path}", flush=True)


if __name__ == "__main__":
    main()

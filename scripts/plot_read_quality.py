"""
Plot the distribution of per-read average Phred quality scores for a FASTQ file.

Accepts gzipped or plain FASTQ. Marks the Q30 threshold, a common baseline
for high-confidence base calls.

Usage:
    python plot_read_quality.py <reads.fastq[.gz]> <output_plot.png>
"""

import sys
import gzip
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

Q30 = 30


def open_fastq(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def main():
    fastq_path = sys.argv[1]
    out_path = sys.argv[2]

    avg_quals = []
    with open_fastq(fastq_path) as f:
        for i, line in enumerate(f):
            if i % 4 == 3:
                quals = [ord(c) - 33 for c in line.strip()]
                avg_quals.append(sum(quals) / len(quals))

    n = len(avg_quals)
    frac_q30 = sum(q >= Q30 for q in avg_quals) / n

    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.hist(avg_quals, bins=60, color="#6699CC", edgecolor="white")
    ax.axvline(Q30, linestyle="--", color="#EC5F67")
    ax.set_xlabel("Average Phred quality per read")
    ax.set_ylabel("Number of reads")
    ax.set_title("Read quality distribution")
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)

    print(f"reads: {n}, mean avg quality: {sum(avg_quals)/n:.1f}, fraction >= Q30: {frac_q30:.3f}")


if __name__ == "__main__":
    main()

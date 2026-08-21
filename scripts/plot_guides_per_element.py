"""
Plot the distribution of gRNAs per targeted element (pRE) in a gRNA library FASTA.

Guide names of the form "<element_coordinates>.<guide_index>" (e.g.
chr9:93276907-93277330.0) are grouped by element coordinates; headers not
starting with "chr" (e.g. non-targeting controls) are excluded, since they
don't represent a tiled element.

Usage:
    python plot_guides_per_element.py <gRNA_library.fasta> <output_plot.png>
"""

import sys
from collections import Counter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def main():
    fasta_path = sys.argv[1]
    out_path = sys.argv[2]

    element_counts = Counter()
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">chr"):
                guide_id = line[1:].strip()
                element = guide_id.rsplit(".", 1)[0]
                element_counts[element] += 1

    counts = list(element_counts.values())
    max_count = max(counts)

    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.hist(counts, bins=range(1, max_count + 2), color="#6699CC", edgecolor="white", align="left")
    ax.set_xlabel("gRNAs per pRE")
    ax.set_ylabel("Number of pREs")
    ax.set_title("gRNA tiling density per pRE")
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)

    print(f"pREs: {len(counts)}, guides: {sum(counts)}, mean/pRE: {sum(counts)/len(counts):.1f}, max: {max_count}")


if __name__ == "__main__":
    main()

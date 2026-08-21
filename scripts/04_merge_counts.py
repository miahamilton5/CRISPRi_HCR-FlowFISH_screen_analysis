"""
Merge per-sample .counts files (from 03_align_and_count.sh) into a single
raw count table, ready for QC and RPM normalization (05_qc_and_normalize_counts.R).

Guide names in the gRNA library FASTAs are of the form
"<element_coordinates>.<guide_index>" (e.g. chr9:93276907-93277330.0); the
element coordinates (everything before the last ".N") are used as the "id"
(gene/element) column that MAGeCK RRA aggregates guides by.

Usage:
    python 04_merge_counts.py <output_table.csv> <sample1>=<sample1.counts> [<sample2>=<sample2.counts> ...]

Example:
    python 04_merge_counts.py GD_FAM120A_raw_counts.csv \
        GD_FAM120A_b15_rep1=GD_FAM120A_b15_rep1.counts \
        GD_FAM120A_t15_rep1=GD_FAM120A_t15_rep1.counts \
        GD_FAM120A_bulk_rep1=GD_FAM120A_bulk_rep1.counts
"""

import sys
import pandas as pd


def load_counts(path):
    return pd.read_csv(path, sep="\t", header=None, names=["grna", "count"])


def main():
    out_path = sys.argv[1]
    sample_args = sys.argv[2:]

    table = None
    for arg in sample_args:
        sample_name, counts_path = arg.split("=", 1)
        df = load_counts(counts_path).rename(columns={"count": sample_name})
        table = df if table is None else table.merge(df, on="grna", how="outer")

    table = table.fillna(0)
    table["id"] = table["grna"].str.rsplit(".", n=1).str[0]

    sample_cols = [c for c in table.columns if c not in ("grna", "id")]
    table = table[["grna", "id"] + sample_cols]
    table.to_csv(out_path, index=False)


if __name__ == "__main__":
    main()

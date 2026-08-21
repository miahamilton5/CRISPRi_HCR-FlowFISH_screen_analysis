#!/usr/bin/env bash
# Trim a demultiplexed Read 1 FASTQ to the 20-bp protospacer window.
#
# Each sample's inline barcode is a different length (3-9 nt, depending on the
# staggered barcode used for that sample), so the start position of the
# protospacer within Read 1 differs by sample. This trims to the 20-bp window
# immediately following that sample's barcode + the 18-bp U6 handle
# (GGAAAGGACGAAACACCG).
#
# Usage:
#   ./02_trim_reads.sh <barcode_length> <read1.fastq> <read1_trimmed.fastq>
#
# Example (5-nt barcode):
#   ./02_trim_reads.sh 5 sample_read1.fastq sample_read1_trimmed.fastq

set -euo pipefail

BARCODE_LENGTH="$1"
INFILE="$2"
OUTFILE="$3"

U6_HANDLE_LEN=18
PROTOSPACER_LEN=20

FIRST_BP=$((BARCODE_LENGTH + U6_HANDLE_LEN + 1))
LAST_BP=$((FIRST_BP + PROTOSPACER_LEN - 1))

fastx_trimmer -Q33 -f "$FIRST_BP" -l "$LAST_BP" -i "$INFILE" -o "$OUTFILE"

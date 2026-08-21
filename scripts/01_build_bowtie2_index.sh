#!/usr/bin/env bash
# Build a Bowtie2 index from a locus gRNA library FASTA.
#
# Usage:
#   ./01_build_bowtie2_index.sh <gRNA_library.fasta> <index_prefix>
#
# Example:
#   ./01_build_bowtie2_index.sh ../gRNA_libraries/FAM120A_gRNA_lib.fasta FAM120A_index

set -euo pipefail

FASTA="$1"
INDEX_PREFIX="$2"

bowtie2-build -f "$FASTA" "$INDEX_PREFIX"

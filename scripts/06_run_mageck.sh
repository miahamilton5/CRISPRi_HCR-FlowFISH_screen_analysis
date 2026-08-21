#!/usr/bin/env bash
# Run MAGeCK RRA comparing bottom (treatment) vs. top (control) sorted bins.
#
# Usage:
#   ./06_run_mageck.sh <count_table.csv> <output_prefix> <treatment_samples> <control_samples> <control_gRNAs.txt>
#
# Example (4 replicates):
#   ./06_run_mageck.sh GD_FAM120A_RPM.csv GD_FAM120A \
#     GD_FAM120A_b15_rep1,GD_FAM120A_b15_rep2,GD_FAM120A_b15_rep3,GD_FAM120A_b15_rep4 \
#     GD_FAM120A_t15_rep1,GD_FAM120A_t15_rep2,GD_FAM120A_t15_rep3,GD_FAM120A_t15_rep4 \
#     ../gRNA_libraries/control_gRNAs.txt

set -euo pipefail

COUNT_TABLE="$1"
OUTPUT_PREFIX="$2"
TREATMENT_SAMPLES="$3"
CONTROL_SAMPLES="$4"
CONTROL_GRNAS="$5"

mageck test \
  -k "$COUNT_TABLE" \
  -t "$TREATMENT_SAMPLES" \
  -c "$CONTROL_SAMPLES" \
  -n "$OUTPUT_PREFIX" \
  --paired \
  --norm-method control \
  --control-sgrna "$CONTROL_GRNAS" \
  --pdf-report

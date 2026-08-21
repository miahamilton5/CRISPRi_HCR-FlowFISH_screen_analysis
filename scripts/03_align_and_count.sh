#!/usr/bin/env bash
# Align trimmed Read 1 to the gRNA library index and generate per-guide read counts.
#
# Usage:
#   ./03_align_and_count.sh <index_prefix> <read1_trimmed.fastq> <sample_name>
#
# Produces:
#   <sample_name>.sam / .bam / .sorted.bam / .sorted.bam.bai
#   <sample_name>.counts               (guide name <tab> read count)
#   <sample_name>.alignment_summary.txt (bowtie2's read count and % alignment stats)

set -euo pipefail

INDEX_PREFIX="$1"
TRIMMED_FASTQ="$2"
SAMPLE="$3"

bowtie2 --end-to-end -a --norc -x "$INDEX_PREFIX" -U "$TRIMMED_FASTQ" -S "${SAMPLE}.sam" 2> "${SAMPLE}.alignment_summary.txt"
cat "${SAMPLE}.alignment_summary.txt"

samtools view -bS "${SAMPLE}.sam" > "${SAMPLE}.bam"
samtools sort "${SAMPLE}.bam" -o "${SAMPLE}.sorted.bam"
samtools index "${SAMPLE}.sorted.bam"
samtools idxstats "${SAMPLE}.sorted.bam" | awk '{print $1"\t"$3}' | grep -v "^\*" > "${SAMPLE}.counts"

#!/bin/bash
# =============================================================================
# Script:      consensus_peaks.sh
# Author:      Christopher M. Rogers (https://github.com/RogersChrisM/)
# Description:
#   Generates high-confidence consensus peaksets from per-replicate MACS2
#   narrowPeak files. Supports ATAC-seq and occupancy (ChIP-seq/CUT&RUN)
#   peaksets. Automatically infers data type and merge distance from input
#   filename conventions. Requires all input files to be the same data type.
#   Outputs union, counted, and consensus BED files.
#
# Usage:
#   ./consensus_peaks.sh <peak_list.txt> [min_overlap]
#
#   peak_list.txt : single-column file listing per-replicate narrowPeak files
#   min_overlap   : optional, minimum number of files required to support a
#                   peak (default: majority, i.e. floor(N/2) + 1)
#
# Input Conventions:
#   ATAC-seq    : *.markdup.rmchrM_peaks.rmblck.narrowPeak  (merge dist: 50 bp)
#   Occupancy   : *.markdup.uq_peaks.rmblck.narrowPeak      (merge dist: 0 bp)
#
# Outputs:
#   <prefix>_union_peaks.bed            : all merged peaks across replicates
#   <prefix>_union_peaks_with_counts.bed: union peaks with per-file overlap counts
#   <prefix>_consensus_peaks.bed        : high-confidence consensus peaks (BED3)
#
# Examples:
#   # ATAC consensus, default majority overlap
#   ./consensus_peaks.sh atac_peak_list.txt
#
#   # Occupancy consensus, explicit min_overlap
#   ./consensus_peaks.sh bcl11a_peak_list.txt 2
#
# Dependencies:
#   bedtools (PATH or module system)
#
# Associated Package:
#   HATKit
#
# Creation Date: 2025-04-12T18:05:00
#   Host: L241568
#   OS: Darwin 25.4.0
#   Bash: 3.2.57(1)-release
#   User: crogers
# =============================================================================

source "$(dirname "$0")/utils.sh"

set -euo pipefail

if ! load_dependency "bedtools"; then
        exit 1
fi

PEAK_LIST=$1
MIN_OVERLAP_ARG=${2:-""}
if [[ ! -f "$PEAK_LIST" ]]; then
    echo "Error: Peak list file '$PEAK_LIST' not found."
    exit 1
fi

PREFIX=$(basename "$PEAK_LIST" | sed 's/\..*$//')
mapfile -t PEAKS < "$PEAK_LIST"

N_FILES=${#PEAKS[@]}
if (( N_FILES == 0 )); then
        echo "ERROR: no peak files in input: '$PEAK_LIST'."
        exit 1
fi

MIN_OVERLAP=${MIN_OVERLAP_ARG:-$(( (N_FILES / 2) + 1 ))}
if (( MIN_OVERLAP > N_FILES )); then
        echo "ERROR: min_overlap ($MIN_OVERLAP) exceeds number of peak files ($N_FILES)."
        exit 1
fi

FIRST_FILE=$(head -1 "$PEAK_LIST")
if [[ "$FIRST_FILE" == *"rmchrM"* ]]; then
        MERGE_DIST=50
        EXPECTED_PATTERN="rmchrM"
        echo "Detected ATAC-seq peaks, using merge distance: ${MERGE_DIST} bp"
elif [[ "$FIRST_FILE" == *".uq_"* ]]; then
        MERGE_DIST=0
        EXPECTED_PATTERN="uq_"
        echo "Detected Occupancy peaks, using merge distance: ${MERGE_DIST} bp"
else
        echo "ERROR: could not infer data type from filename '${FIRST_FILE}'."
        exit 1
fi
for FILE in "${PEAKS[@]}"; do
        if [["$FILE" != *"$EXPECTED_PATTERN"* ]]; then
                echo "ERROR: file '${FILE}' does not match the expected pattern '${EXPECTED_PATTERN}'."
                echo "Please ensure all input files must be the same data type (ATAC or occupancy)."
                exit 1
        fi
done

cat "${PEAKS[@]}" | sort -k1,1 -k2,2n | bedtools merge -d "$MERGE_DIST" -i - > "${PREFIX}_union_peaks.bed"
bedtools intersect -a "${PREFIX}_union_peaks.bed" -b "${PEAKS[@]}" -c > "${PREFIX}_union_peaks_with_counts.bed"
awk -v min="$MIN_OVERLAP" '$4 >= min {print $1"\t"$2"\t"$3}' "${PREFIX}_union_peaks_with_counts.bed" > "${PREFIX}_consensus_peaks.bed"

echo "Done! Outputs:"
echo " - ${PREFIX}_union_peaks.bed"
echo " - ${PREFIX}_union_peaks_with_counts.bed"
echo " - ${PREFIX}_consensus_peaks.bed (min_overlap=$MIN_OVERLAP)"
# --- Signature ---
# Author: CM Rogers (https://github.com/RogersChrisM/)
# Date: 2026-05-07
# SHA256: e1ef7db72bc3abd0db01a10dff83d48778540cd78bfc035368e28b37f51241f2

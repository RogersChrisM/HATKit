#!/bin/bash
# =============================================================================
# Script:      consensus_peaks.sh
# Author:      Christopher M. Rogers (https://github.com/RogersChrisM/)
# Description:
#   Generates high-confidence consensus peaksets from per-replicate MACS2
#   narrowPeak or BED files. Supports ATAC-seq and occupancy (ChIP-seq/CUT&RUN)
#   peaksets. Assay type is explicitly declared in the input peak list and used
#   to determine the appropriate merge distance. Requires all input files to be
#   the same assay type. Outputs union, counted, and consensus BED files.
#
# Usage:
#   ./consensus_peaks.sh <peak_list.txt> [min_overlap]
#
#   peak_list.txt : tab-delimited, two-column file:
#                     col 1 — path to per-replicate .narrowPeak or .bed file
#                     col 2 — assay type (e.g. ATAC, ATAC-seq, ChIP, ChIP-seq,
#                             CUT&RUN, CUTNRUN, CNR)
#   min_overlap   : optional, minimum number of files required to support a
#                   peak (default: majority, i.e. floor(N/2) + 1)
#
# Input Conventions:
#   ATAC-seq          : merge dist 50 bp
#   ChIP-seq, CUT&RUN : merge dist 0 bp
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
# Peak List Format:
#   rep1.markdup.rmchrM_peaks.rmblck.narrowPeak    ATAC-seq
#   rep2.markdup.rmchrM_peaks.rmblck.narrowPeak    ATAC-seq
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

set -euo pipefail

source "$(dirname "$0")/utils.sh"

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

declare -a PEAKS
declare -a ASSAYS
while IFS=$'\t' read -r FILE ASSAY; do
    PEAKS+=("$FILE")
    ASSAYS+=("$ASSAY")
done < "$PEAK_LIST"

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

NORMALIZED_ASSAY=$(echo "${ASSAYS[0],,}" | sed 's/-.*$//')
case "$NORMALIZED_ASSAY" in
    atac)
        MERGE_DIST=50
        DATA_TYPE="ATAC-seq"
        ;;
    chip)
        MERGE_DIST=0
        DATA_TYPE="ChIP-seq"
        ;;
    cut\&run|cutnrun|cnr)
        MERGE_DIST=0
        DATA_TYPE="CUT&RUN"
        ;;
    *)
        echo "ERROR: unrecognized assay type '${ASSAYS[0]}'."
        echo "Accepted types: ATAC, ATAC-seq, ChIP, ChIP-seq, CUT&RUN, CUTNRUN, CNR"
        exit 1
        ;;
esac
echo "Detected ${DATA_TYPE}, using merge distance: ${MERGE_DIST} bp"

for i in "${!PEAKS[@]}"; do
    FILE="${PEAKS[$i]}"
    ASSAY=$(echo "${ASSAYS[$i],,}" | sed 's/-.*$//')

    if [[ "$FILE" != *.narrowPeak && "$FILE" != *.bed ]]; then
        echo "ERROR: file '${FILE}' is not a .narrowPeak or .bed file."
        echo "Accepted formats: .narrowPeak, .bed"
        exit 1
    fi

    if [[ "$ASSAY" != "$NORMALIZED_ASSAY" ]]; then
        echo "ERROR: file '${FILE}' has assay type '${ASSAYS[$i]}' but expected '${ASSAYS[0]}'."
        echo "All input files must be the same assay type."
        exit 1
    fi
done
echo "All input files confirmed as ${DATA_TYPE}."

cat "${PEAKS[@]}" | sort -k1,1 -k2,2n | bedtools merge -d "$MERGE_DIST" -i - > "${PREFIX}_union_peaks.bed"
bedtools intersect -a "${PREFIX}_union_peaks.bed" -b "${PEAKS[@]}" -c > "${PREFIX}_union_peaks_with_counts.bed"
awk -v min="$MIN_OVERLAP" '$4 >= min {print $1"\t"$2"\t"$3}' "${PREFIX}_union_peaks_with_counts.bed" > "${PREFIX}_consensus_peaks.bed"

echo "Done! Outputs:"
echo " - ${PREFIX}_union_peaks.bed"
echo " - ${PREFIX}_union_peaks_with_counts.bed"
echo " - ${PREFIX}_consensus_peaks.bed (assay=${DATA_TYPE}, n_files=${N_FILES}, min_overlap=${MIN_OVERLAP}, merge_dist=${MERGE_DIST}bp)"
# --- Signature ---
# Author: CM Rogers (https://github.com/RogersChrisM/)
# Date: 2026-05-07
# SHA256: 8963ff153c7e468a655ed0ef4951b014f9f043a58a7a487f21b10db997e92be0

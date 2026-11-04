#!/bin/bash
#
# Pulls primary metabolite KEGG pathways that differ significantly in 
# percent completeness between HA and FL for further investigation.
#
# For each pathway group, filters for host_association_sig == "yes",
# sorts by host_association_mw_stat descending, and takes the top 10.
#
# Optionally takes a list of pathway IDs (col 1) to include
# regardless of significance.
# These are appended after the significant hits and deduplicated.
#
# Column reference (pathway_sumstats.tsv):
#   1=Pathway  2=Group  3=Supergroup  4=Annot  5=Orthcounts
#   12=host_association_mw_stat  15=host_association_sig
#   22=aerotolerance_mw_stat     25=aerotolerance_sig
#
# To look at differences in aerotolerance instead, change: 
#    COL_MW_STAT=12
#    COL_SIG=15
# To: 
#    COL_MW_STAT=22
#    COL_SIG=25
#
#
# Usage:
#   bash pull_pathways.sh \
#       pathway_sumstats.tsv \
#       key_pathways.tsv \
#       [optional_pathway_list.txt]
#
set -euo pipefail

KOFAM_SUM=$1
OUT_TSV=$2
OPTIONAL_LIST=${3:-""}   # path to optional pathway ID list, or empty

# Column indices
COL_PATHWAY=1
COL_GROUP=2
COL_MW_STAT=12
COL_SIG=15

path_groups=(
    "Metabolism of cofactors and vitamins"
    "Amino acid metabolism"
    "Lipid metabolism"
)

# Write header once
head -n 1 "${KOFAM_SUM}" > "${OUT_TSV}"

# Temp file to collect all significant rows before deduplication
TMP_SIG=$(mktemp)
trap "rm -f $TMP_SIG" EXIT

# Significant rows: filter for sig, sort by MW stat desc, top 10 per group
for pathway in "${path_groups[@]}"; do
    awk -F'\t' \
        -v group="$pathway" \
        -v sig_col="$COL_SIG" \
        -v group_col="$COL_GROUP" \
        '$group_col == group && $sig_col == "yes"' \
        "${KOFAM_SUM}" \
    | sort -t$'\t' -k"${COL_MW_STAT}","${COL_MW_STAT}" -rn \
    | head -n 10 \
    >> "$TMP_SIG"
done

# Optional forced-inclusion pathways
TMP_FORCED=$(mktemp)
trap "rm -f $TMP_SIG $TMP_FORCED" EXIT

if [ -n "$OPTIONAL_LIST" ] && [ -f "$OPTIONAL_LIST" ]; then
    echo "Including forced pathways from ${OPTIONAL_LIST}..." >&2

    while IFS= read -r pid || [ -n "$pid" ]; do
        pid=$(echo "$pid" | tr -d '[:space:]')
        [ -z "$pid" ] && continue

        # Pull the row for this pathway ID from the full summary file
        awk -F'\t' \
            -v pid="$pid" \
            -v id_col="$COL_PATHWAY" \
            'NR > 1 && $id_col == pid' \
            "${KOFAM_SUM}" \
        >> "$TMP_FORCED"
    done < "$OPTIONAL_LIST"

    echo "  $(wc -l < "$TMP_FORCED") forced pathway rows found" >&2
fi

# Combine significant rows with forced rows
# Deduplicate on pathway ID by keeping first occurrence
cat "$TMP_SIG" "$TMP_FORCED" \
    | awk -F'\t' -v id_col="$COL_PATHWAY" \
        '!seen[$id_col]++' \
    >> "${OUT_TSV}"

# Summary
n_sig=$(wc -l < "$TMP_SIG")
n_forced=$(wc -l < "$TMP_FORCED")
n_total=$(( $(wc -l < "${OUT_TSV}") - 1 ))   # subtract header

echo "Done." >&2
echo "  Significant rows : ${n_sig}" >&2
echo "  Forced rows      : ${n_forced}" >&2
echo "  Total in output  : ${n_total} (after deduplication)" >&2
echo "  Output           : ${OUT_TSV}" >&2

#!/bin/bash

#Runs Diamond BLAST on a single FAA file.

#Usage: diamond.sh \
#    path/to/proteome.faa \
#    path/to/output/directory \
#    path/to/diamond/database \
#    path/to/temporary/directory

set -euo pipefail
#need diamond env

#set variables
PROTEOME=$1
OUTDIR=$2
DIAMOND_DB=$3
TMP_DIR=$4

BASENAME=$(basename "${PROTEOME}" .faa)
SAMPLE_OUTDIR="${OUTDIR}/${BASENAME}"

mkdir -p "${OUTDIR}"
mkdir -p "${SAMPLE_OUTDIR}"
mkdir -p "${TMP_DIR}"

# for each faa run blast against itself to get -bs
diamond blastp \
    --query "${PROTEOME}" \
    --db "${PROTEOME}" \
    --threads 12 \
    --out "${SAMPLE_OUTDIR}/${BASENAME}"_BS.out

# run diamond blastp to get -bn
diamond blastp \
    --query "${PROTEOME}" \
    --db "${DIAMOND_DB}" \
    --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen slen qlen nident positive staxids qcovhsp \
    --max-target-seqs 500 \
    --threads 12 \
    -b 2 \
    -c 1 \
    --tmpdir "${TMP_DIR}" \
    --out "${SAMPLE_OUTDIR}/${BASENAME}"_BN.out

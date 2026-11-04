#!/bin/bash
#
#Runs metaeuk easy-predict for gene prediction on a single genome.
#
#Usage: metaeuk.sh \
#    path/to/genome.fna \
#    path/to/output/directory \
#    path/to/metaeuk.sif \
#    path/to/protein_db

set -euo pipefail

# Genome file
GENOME=$1

# Output directory
OUTDIR=$2
mkdir -p "${OUTDIR}"

# MetaEuk apptainer image
METAEUK=$3

# Protein db
PROT=$4

BASENAME=$(basename "${GENOME}" .fna)

apptainer exec "${METAEUK}" \
    metaeuk easy-predict \
    "${GENOME}" \
    "${PROT}" \
    "${OUTDIR}/metaeuk_out" \
    "${OUTDIR}/tmp"

echo "MetaEuk completed for ${BASENAME}"

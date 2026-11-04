#!/bin/bash

#Runs BUSCO for completeness on a single FAA file in protein mode.

#Usage: busco.sh \
#    path/to/proteome.faa \
#    path/to/output/directory \
#    path/to/bind/directory \
#    path/to/busco.sif \
#    path/to/lineage

set -euo pipefail

#set variables
PROTEOME=$1   # Input FAA file
OUTDIR=$2     # Output directory
BIND=$3       # Host path to bind-mount into the container
BUSCO_SIF=$4  # Apptainer image
LINEAGE=$5    # BUSCO lineage dataset path

BASENAME=$(basename "${PROTEOME}" .faa)

mkdir -p "${OUTDIR}"

#protein flag + utilizing eukaryota_odb12 database
apptainer exec --bind "${BIND}:${BIND}" "${BUSCO_SIF}" \
    busco \
    --in "${PROTEOME}" \
    --out "${BASENAME}" \
    --out_path "${OUTDIR}" \
    --lineage_dataset "${LINEAGE}" \
    --mode proteins \
    --cpu "${SLURM_CPUS_PER_TASK:-1}" \
    --force

echo "BUSCO completed for ${BASENAME}"

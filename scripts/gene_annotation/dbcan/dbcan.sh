#!/bin/bash

#Runs dbcan for gene annotation on a single FAA file.

#Usage: dbcan.sh \
#    path/to/proteome.faa \
#    path/to/output/directory \
#    path/to/bind/directory \
#    path/to/dbcan.sif \
#    path/to/dbcan/protein/directory


#set variables
PROTEOME=$1   # Input FAA file
OUTDIR=$2     # Output directory
BIND=$3       # Host path to bind-mount into the container
DBCAN_SIF=$4  # Apptainer image
DBCAN_DIR=$5    # Path to dbcan protein database

BASENAME=$(basename "${PROTEOME}" .faa)

mkdir -p "${OUTDIR}"

# Run dbcan
apptainer exec --bind "${BIND}:${BIND}" "${DBCAN_SIF}" \
    run_dbcan \
    --db_dir "${DBCAN_DIR}" \
    --hmm_cpu 12 \
    --dia_cpu 12 \
    --out_dir "${OUTDIR}" \
    --out_pre "${BASENAME}_" \
    "${PROTEOME}" protein

echo "Dbcan completed for ${BASENAME}"

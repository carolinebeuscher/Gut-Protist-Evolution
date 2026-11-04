#!/bin/bash

#Runs KofamScan for gene annotation on a single FAA file.

#Usage: kofam.sh \
#    path/to/proteome.faa \
#    path/to/output/directory \
#    path/to/kofam.sif \
#    path/to/kofam/profiles \
#    path/to/list/of/kos
set -euo pipefail

#set variables
PROTEOME=$1   # Input FAA file
OUTDIR=$2     # Output directory
KOFAM_SIF=$3  # Apptainer image
KOFAM_PROFILES=$4   #kofam profiles
KO_LIST=$5          #list of KOs


BASENAME=$(basename "${PROTEOME}" .faa)
mkdir -p "${OUTDIR}"


apptainer exec "${KOFAM_SIF}" exec_annotation \
    -p "${KOFAM_PROFILES}" \
    -k "${KO_LIST}" \
    --cpu "${SLURM_CPUS_PER_TASK:-1}" \
    -f detail-tsv \
    --e-value 1e-5 \
    --tmp-dir "${OUTDIR}/${BASENAME}_tmp" \
    -o "${OUTDIR}/${BASENAME}_ko_annotations.tsv" \
    "${PROTEOME}"

echo "KofamScan completed for ${BASENAME}"

#!/bin/bash

#Runs Interproscan for gene annotation on a single FAA file.

#Usage: interproscan.sh \
#    path/to/proteome.faa \
#    path/to/output/directory \
#    path/to/interproscan/script \
#    number_of_cpus

#~/shared/apps/interproscan-5.66-98.0/interproscan.sh

set -euo pipefail

#set variables
PROTEOME=$1   # Input FAA file
OUTDIR=$2     # Output directory
INTERPRO=$3   # location of interproscan.sh
CPU=$4        # number of cpus used to run 

BASENAME=$(basename "${PROTEOME}" .faa)
mkdir -p "${OUTDIR}"

# run InterProScan
"${INTERPRO}" \
    -exclappl ProSiteProfiles,ProSitePatterns \
    -cpu "${CPU}" \
    -iprlookup \
    -goterms \
    -pa \
    -d "${OUTDIR}" \
    -f TSV \
    -i "${PROTEOME}"

# Rename to the expected output filename
mv "${OUTDIR}/${BASENAME}.faa.tsv" "${OUTDIR}/${BASENAME}_interproscan.tsv"

echo "InterProScan completed for ${BASENAME}"



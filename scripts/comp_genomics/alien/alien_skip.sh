#!/bin/bash

#Runs Alien Index for HGT analysis on a single sample.

#Usage: alien_skip.sh \
#    path/to/output/directory \
#    path/to/alien/script \
#    path/to/blast/self/output \
#    path/to/blast/others/output \
#    recipient/id \
#    list/of/skip/lineages \
#    path/to/sqlite/taxa/db

set -euo pipefail

#need env with ete3, python, and diamond

#set variables
OUTDIR=$1       # Output directory
ALIEN=$2        # location of alien_index_multiskip.py
BS=$3           # blast output against self 
BN=$4           # blast output against others
RECIPIENT=$5    #NCBI taxonomy ID of the recipient/focal group (-r); 2759 for eukaryotes
SKIP=$6         #skipped taxon #5741,42740,163343,555406,2611341,5794
SQLITE=$7       #path to sqlite taxa db 

BASENAME=$(basename "${BN}" _BN.out)
mkdir -p "${OUTDIR}"

#run alien
python "${ALIEN}" \
    -bs "${BS}" \
    -bn "${BN}" \
    -r "${RECIPIENT}" \
    -s "${SKIP}" \
    -t "${SQLITE}" \
    -o "${OUTDIR}/${BASENAME}/${BASENAME}"_alien.tsv


echo "Alien completed for ${BASENAME}"

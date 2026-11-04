#!/bin/bash

#Runs Braker3 for gene prediction on a single genome.
#
#Usage: braker.sh \
#    path/to/genome.fna \
#    path/to/output/directory \
#    path/to/braker.sif \
#    path/to/protein_db \
#    name/of/species



GENOME=$1           # Genome file
OUTDIR=$2           # Output directory
BRAKER=$3           # Braker apptainer image
PROT=$4             # Protein db
SPECIES=${5:-""}    # Optional existing species 
    #--species=tetrahymena
BASENAME=$(basename "${GENOME}" .fna)   #basename of sample

mkdir -p "${OUTDIR}"

# Augustus Config directory — must be writable and unique per job
export AUGUSTUS_CONFIG_PATH="${OUTDIR}/augustus_config"
mkdir -p "${AUGUSTUS_CONFIG_PATH}"

# Execute braker
[ -z "$SPECIES" ] && SPECIES="Sp_${BASENAME}"

apptainer exec "${BRAKER}" braker.pl \
    --genome="${GENOME}" \
    --prot_seq="${PROT}" \
    --workingdir="${OUTDIR}" \
    --threads="${SLURM_CPUS_PER_TASK:-1}" \
    --min_contig=10000 \
    --species="${SPECIES}" \
    --useexisting \
    --AUGUSTUS_CONFIG_PATH="${AUGUSTUS_CONFIG_PATH}"


echo "BRAKER completed for ${BASENAME}"

#!/bin/bash

#Runs Orthofinder for comparative genomics on a directory of fasta files.
#Flags: 
#    -S: use diamond blast 
#    -og: stop the analysis immediately after inferring orthogroups 

#Usage: orthofinder.sh \
#    path/to/proteome.faa/dir \
#    path/to/outdir \
#    number/of/threads \
#    path/to/previous/results/dir

#need orthofinder in env

set -euo pipefail

# -a (OrthoFinder algorithm thread count) is kept low deliberately: at
# higher values (was 6) the "Initial processing of each species" step has
# been observed to stall indefinitely ("ERROR: Stalled for 200.0s") partway
# through, wasting a full all-vs-all diamond run (~1.5 days for 204
# species). Lower -a is slower per-species but has not reproduced the
# stall. Revisit if a future OrthoFinder version fixes the underlying
# multiprocessing issue.

#set variables
PROTEOME_DIR=$1   # path to faa directory
OUTDIR=$2
THREADS=${3:-${SLURM_CPUS_PER_TASK:-3}}
PREV_RESULTS=${4:-""}   # optional: path to a previous Results_* directory

# Execute program
if [ -z "$PREV_RESULTS" ]; then
    orthofinder \
        -S diamond \
        -og \
        -f "${PROTEOME_DIR}" \
        -o "${OUTDIR}" \
	-n "here" \
        -t "${THREADS}" \
        -a 1
else
    # -b: pass in previous results dir so OrthoFinder only recomputes
    # what's needed for the new samples.
    #
    # OrthoFinder's add-species workflow requires -f to contain ONLY
    # species that are NOT already in the previous results -- passing
    # already-present species again is unsupported and can error or
    # produce duplicate/incorrect entries. PROTEOME_DIR is the full
    # current protein_dir though, which may substantially overlap with
    # PREV_RESULTS (e.g. a baseline covering only a subset of clades).
    # So build a staging dir containing symlinks to just the species in
    # PROTEOME_DIR that are absent from PREV_RESULTS/SpeciesIDs.txt, and
    # pass that as -f instead.
    SPECIES_IDS_FILE="${PREV_RESULTS}/SpeciesIDs.txt"
    if [ ! -f "$SPECIES_IDS_FILE" ]; then
        echo "ERROR: ${SPECIES_IDS_FILE} not found."
        echo "PREV_RESULTS should be a previous run's WorkingDirectory"
        echo "(contains SpeciesIDs.txt, SequenceIDs.txt, etc.), not the"
        echo "top-level Results_* directory."
        exit 1
    fi

    NEW_SPECIES_DIR=$(mktemp -d)
    trap 'rm -rf "$NEW_SPECIES_DIR"' EXIT

    N_NEW=0
    echo "[$(date)] Diffing ${PROTEOME_DIR} against ${SPECIES_IDS_FILE}:"
    for f in "${PROTEOME_DIR}"/*.faa; do
        [ -e "$f" ] || continue
        base=$(basename "$f")
        if grep -qF "$base" "$SPECIES_IDS_FILE"; then
            echo "  already present, skipping: $base"
        else
            echo "  new, adding: $base"
            ln -s "$(realpath "$f")" "${NEW_SPECIES_DIR}/${base}"
            N_NEW=$((N_NEW + 1))
        fi
    done

    if [ "$N_NEW" -eq 0 ]; then
        echo "ERROR: no species in ${PROTEOME_DIR} are new relative to"
        echo "${SPECIES_IDS_FILE} -- nothing to add."
        exit 1
    fi
    echo "[$(date)] Adding ${N_NEW} new species to previous results at ${PREV_RESULTS}"

    # OrthoFinder rejects -o outright whenever -b is used (even together
    # with -f): "ERROR: Incompatible arguments, -o (non-default output
    # directory) can only be used with a new OrthoFinder run using option
    # -f". Without -o, add-mode results land under
    # <NEW_SPECIES_DIR>/OrthoFinder/Results_<name> — let it write there,
    # then move that into the stable OUTDIR/Results_here location after.
    RESULTS_PARENT="${NEW_SPECIES_DIR}/OrthoFinder"
    rm -rf "${RESULTS_PARENT}/Results_here"

    orthofinder \
        -f "${NEW_SPECIES_DIR}" \
        -b "${PREV_RESULTS}" \
	-n "here" \
        -t "${THREADS}" \
        -a 1

    if [ ! -d "${RESULTS_PARENT}/Results_here" ]; then
        echo "ERROR: expected OrthoFinder to write ${RESULTS_PARENT}/Results_here"
        echo "Contents of ${RESULTS_PARENT}:"
        ls -la "${RESULTS_PARENT}" 2>&1 || echo "  (directory does not exist)"
        exit 1
    fi

    rm -rf "${OUTDIR}/Results_here"
    mv "${RESULTS_PARENT}/Results_here" "${OUTDIR}/Results_here"
    rm -rf "${RESULTS_PARENT}"
fi

echo "Orthofinder completed."






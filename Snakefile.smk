"""
Snakemake workflow for host-associated protist convergent evolution

Requires: .fna genome files already downloaded.

Run (dry-run first):
    snakemake -s Snakefile.smk --profile profiles/slurm -n

Run:
    snakemake -s Snakefile.smk --profile profiles/slurm
"""

import os
import glob

configfile: "config/config.yaml"

#############################################################################
# Define directories and collect samples from the genome directory
#############################################################################
# Resolve {root} and {project} placeholders in all config string values.
# To switch users, change only config["root"] and config["project"].
_ROOT    = config["root"]
_PROJECT = config["project"]   # relative to root; full path = root/project

def p(key):
    """Return config[key] with {root} and {project} substituted."""
    val = config[key]
    if isinstance(val, str):
        return val.format(root=_ROOT, project=_PROJECT)
    return val

def rp(val):
    """Resolve a raw config string value (for use inside rules)."""
    if isinstance(val, str):
        return val.format(root=_ROOT, project=_PROJECT)
    return val

GENOME_DIR      = p("genome_dir")
GP_OUTDIR       = p("gp_outdir")
OUTDIR          = p("outdir")
PROTEIN_DIR     = p("protein_dir")   # existing FAA directory; predicted FAAs land here too
KOFAM_OUTDIR    = p("kofam_outdir")
BUSCO_OUTDIR    = p("busco_outdir")
INTERPRO_OUTDIR = p("interproscan_outdir")
DBCAN_OUTDIR    = p("dbcan_outdir")
ORTHO_OUTDIR    = p("orthofinder_outdir")
DIAMOND_OUTDIR  = p("diamond_outdir")
ALIEN_OUTDIR    = p("alien_outdir")
FILTERED_FAA_DIR = p("filtered_faa_dir")

# Samples that need gene prediction (have .fna but no .faa yet)
FNA_SAMPLES = [
    os.path.basename(f).replace(".fna", "")
    for f in sorted(glob.glob(os.path.join(GENOME_DIR, "*.fna")))
]

# Samples with existing .faa files (already downloaded/predicted)
FAA_SAMPLES = [
    os.path.basename(f).replace(".faa", "")
    for f in sorted(glob.glob(os.path.join(PROTEIN_DIR, "*.faa")))
]

# Full sample set: union of both, preserving order (fna first, then faa-only)
SAMPLES = FNA_SAMPLES + [s for s in FAA_SAMPLES if s not in FNA_SAMPLES]

if not SAMPLES:
    raise ValueError(
        f"No samples found.\n"
        f"  genome_dir ({GENOME_DIR}): {len(FNA_SAMPLES)} .fna files\n"
        f"  protein_dir ({PROTEIN_DIR}): {len(FAA_SAMPLES)} .faa files"
    )

print(f"Samples: {len(FNA_SAMPLES)} from .fna + "
      f"{len(FAA_SAMPLES)} from .faa = {len(SAMPLES)} total")

### HELPER FUNCTIONS ###
#############################################################################
# Grouping helper functions
# resolve per-group column from config group lists 
#       - config["dbcan_groups"]
#       - config["interpro_groups"] 
# They are lists of {group: "Metamonada", col: "clade"} dicts, allowing 
# each group to be looked up in a different column.
# 
# ie metamonada and amoebozoa are in a different clade than stramenopiles, 
# apicomplexans, and ciliates in ncbi so when grouping_key.tsv is created and
# used, we can't just sort on one column
#############################################################################
def _group_col_map(group_list):
    """Return {group_value: col_name} from a list of {group, col} dicts."""
    return {entry["group"]: entry["col"] for entry in (group_list or [])}

def _group_names(group_list):
    """Return just the group name strings from a {group, col} list."""
    return [entry["group"] for entry in (group_list or [])]

DBCAN_GROUP_COL_MAP   = _group_col_map(config.get("dbcan_groups",   []))
INTERPRO_GROUP_COL_MAP = _group_col_map(config.get("interpro_groups", []))
DBCAN_GROUPS          = _group_names(config.get("dbcan_groups",   []))
INTERPRO_GROUPS       = _group_names(config.get("interpro_groups", []))


# module helper function:
# read checkpoint output to get the list of modules to assess
############################################################################
def get_modules(wildcards):
    """
    Reads key_pathways.tsv (produced by the pull_kofam_pathways checkpoint) to get
    the list of significant pathway IDs (column "Pathway"), then maps each pathway
    to its KEGG module ID(s) via pathway_module_ko.tsv ("pathway_id" -> "module_id"),
    excluding NA. key_pathways.tsv itself is pathway-keyed, not module-keyed, so this
    join is required — it is not a plain column read.
    Returns a sorted list of unique module IDs.
    """
    import csv
    tsv = checkpoints.pull_kofam_pathways.get(**wildcards).output.pathway_list
    pathways = set()
    with open(tsv) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            pid = row.get("Pathway", "").strip()
            if pid and pid != "NA":
                pathways.add(pid)

    modules = set()
    with open(rp(config["pathway_module_ko"])) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            pid = row.get("pathway_id", "").strip()
            mid = row.get("module_id", "").strip()
            if pid in pathways and mid and mid != "NA":
                modules.add(mid)
    return sorted(modules)


def get_filtered_samples(wildcards):
    """
    Reads kept_samples.txt produced by the filter_proteomes checkpoint.
    Returns the list of samples that passed BUSCO quality filtering.
    All downstream rules that operate on .faa files use this instead of SAMPLES,
    so quarantined samples are automatically excluded from the rest of the DAG.
    """
    kept_file = checkpoints.filter_proteomes.get(**wildcards).output.kept
    with open(kept_file) as fh:
        # kept_samples.txt contains short_summary filenames; extract sample name
        samples = []
        for line in fh:
            line = line.strip()
            if line:
                # e.g. "short_summary.specific.Gint.txt" -> "Gint"
                name = line.replace("short_summary.specific.", "").replace(".txt", "")
                samples.append(name)
    return samples


############################################################################
# Target rule!!
#      Collects full DAG
############################################################################
rule all:
    input:
        # ── Gene prediction (only for samples with .fna files) ──────────────
        expand(
            os.path.join(GP_OUTDIR, "{sample}", "gene_prediction.done"),
            sample=FNA_SAMPLES,
        ),

        # ── BUSCO completeness ────────────────────────────────────────────────
        os.path.join(BUSCO_OUTDIR, "busco_summary.tsv"),

        # ── Grouping key ──────────────────────────────────────────────────────
        "grouping_key.tsv",

        # ── KofamScan (filtered samples only) ────────────────────────────────
        os.path.join(KOFAM_OUTDIR, "kofam_summary.tsv"),
        os.path.join(KOFAM_OUTDIR, "pathway_sumstats.tsv"),

        # ── Module assessment (checkpoint-gated) ─────────────────────────────
        lambda wc: expand(
            os.path.join(KOFAM_OUTDIR, "module_results", "{module}_functionality.tsv"),
            module=get_modules(wc),
        ),
        os.path.join(KOFAM_OUTDIR, "module_stats.tsv"),

        # ── InterProScan ──────────────────────────────────────────────────────
        os.path.join(INTERPRO_OUTDIR, "all_interpro_summary.tsv"),
        expand(
            os.path.join(INTERPRO_OUTDIR, "stats", "{group}_mannwhitney.tsv"),
            group=INTERPRO_GROUPS,
        ),
        os.path.join(INTERPRO_OUTDIR, "stats", "all_lineages_"
                     + config.get("interpro_stat_type", "mannwhitney") + ".tsv"),

        # ── dbCAN ─────────────────────────────────────────────────────────────
        os.path.join(DBCAN_OUTDIR, "dbcan_combined.tsv"),
        expand(
            os.path.join(DBCAN_OUTDIR, "stats", "{group}_dbcan_stats.tsv"),
            group=DBCAN_GROUPS,
        ),
        os.path.join(DBCAN_OUTDIR, "stats", "all_lineages_dbcan_stats.tsv"),

        # ── OrthoFinder ───────────────────────────────────────────────────────
        os.path.join(ORTHO_OUTDIR, "orthogroup_stats.done"),
        os.path.join(ORTHO_OUTDIR, "sig_orthogroup_annotations.tsv"),

        # ── AlienIndex ────────────────────────────────────────────────────────
        os.path.join(ALIEN_OUTDIR, "alien_summary.tsv"),
        os.path.join(ALIEN_OUTDIR, "alien_ipr_counts.tsv"),



############################################################################
# Rule: gene_prediction
# Scripts: 
#        - scripts/gene_prediction/braker.sh
#        - scripts/gene_prediction/metaeuk.sh
#        - scripts/gene_prediction/reduce_metaeuk_faa_name.py
# Goal: to convert FNA fasta files to FAA files if needed and move those new
#    FAA files to PROTEIN_DIR
# Produces: 
#    - {sample}.faa in PROTEIN_DIR
# Notes: 
#        - If scripts/braker.sh exits non-zero, calls scripts/metaeuk.sh
#        - If metaeuk.sh is called so is reduce_metaeuk_faa_name.py
#                - Raw MetaEuk headers are long pipe-delimited strings; 
#                this script keeps only the first 4 fields 
#                (e.g. contig|start|end|strand) which downstream tools 
#                (InterProScan, kofamscan) handle correctly.
############################################################################
rule gene_prediction:
    input:
        genome = os.path.join(GENOME_DIR, "{sample}.fna"),
        bprot   = rp(config["braker_protein_db"]),
        mprot = rp(config["meta_protein_db"])
    output:
        done = os.path.join(GP_OUTDIR, "{sample}", "gene_prediction.done"),
        merged_faa = os.path.join(PROTEIN_DIR, "{sample}.faa"), 
    params:
        braker_outdir  = os.path.join(GP_OUTDIR, "{sample}", "braker"),
        metaeuk_outdir = os.path.join(GP_OUTDIR, "{sample}", "metaeuk"),
        braker_sif     = rp(config["braker_sif"]),
        metaeuk_sif    = rp(config["metaeuk_sif"]),
        braker_species = config.get("species", ""),
    threads: config.get("threads", 4)

    log:
        os.path.join(GP_OUTDIR, "{sample}", "logs", "gene_prediction.log"),
    conda:
        "envs/biopython.yaml"
    shell:
        r"""
        set -euo pipefail
        exec > >(tee -a {log}) 2>&1
        mkdir -p "$(dirname {log})"

        BRAKER_OK=false
        BRAKER_FAA="{params.braker_outdir}/augustus.hints.aa"

        echo "[$(date)] Attempting BRAKER3 for {wildcards.sample}"

        bash scripts/gene_prediction/braker.sh \
            "{input.genome}" \
            "{params.braker_outdir}" \
            "{params.braker_sif}" \
            "{input.bprot}" \
            "{params.braker_species}" \
            && BRAKER_OK=true \
            || echo "[$(date)] BRAKER3 exited non-zero — will check output before deciding"

        # braker.sh may swallow braker.pl's internal exit
        # code (as seen with the --useexisting / species-dir-exists error).
        # Even if the exit code looked clean, confirm the expected output exists.
        if $BRAKER_OK; then
            if [ ! -s "${{BRAKER_FAA}}" ]; then
                echo "[$(date)] BRAKER3 exited 0 but augustus.hints.aa is missing or empty — treating as failure"
                BRAKER_OK=false
            fi
        else
            # Non-zero exit: also check whether output somehow exists anyway
            # (e.g. a partial run that still produced usable proteins)
            if [ -s "${{BRAKER_FAA}}" ]; then
                echo "[$(date)] BRAKER3 exited non-zero but augustus.hints.aa exists — salvaging output"
                BRAKER_OK=true
            fi
        fi

        if $BRAKER_OK; then
            echo "[$(date)] BRAKER3 succeeded for {wildcards.sample}"
            FAA_SRC="${{BRAKER_FAA}}"
            echo "tool=braker" > {output.done}
        else
            echo "[$(date)] BRAKER3 failed for {wildcards.sample} — falling back to MetaEuk"
            echo "[$(date)] Attempting MetaEuk for {wildcards.sample}"

            bash scripts/gene_prediction/metaeuk.sh \
                "{input.genome}" \
                "{params.metaeuk_outdir}" \
                "{params.metaeuk_sif}" \
                "{input.mprot}"

            echo "[$(date)] MetaEuk succeeded for {wildcards.sample}"

            METAEUK_RAW="{params.metaeuk_outdir}/metaeuk_out.fas"
            METAEUK_REDUCED="{params.metaeuk_outdir}/metaeuk_out_reduced.fas"

            echo "[$(date)] Reducing MetaEuk header names for {wildcards.sample}"
            python3 scripts/gene_prediction/reduce_metaeuk_faa_name.py "${{METAEUK_RAW}}" \
                > "${{METAEUK_REDUCED}}"

            FAA_SRC="${{METAEUK_REDUCED}}"
            echo "tool=metaeuk" > {output.done}
        fi

        echo "[$(date)] Copying ${{FAA_SRC}} → {output.merged_faa}"
        mkdir -p "{PROTEIN_DIR}"
        cp "${{FAA_SRC}}" "{output.merged_faa}"
        echo "[$(date)] Done — {wildcards.sample}.faa is in {PROTEIN_DIR}"
        exec 1>&- 2>&-
        wait
        """


############################################################################
# Rule: proteome_completeness
# Script: 
#    - scripts/data_cleaning/busco.sh
# Goal: to filter out low quality proteomes
# Produces: 
#   - short_summary.specific.{sample}.txt file for each sample
############################################################################
rule proteome_completeness:
    input:
        proteome = os.path.join(PROTEIN_DIR, "{sample}.faa"),
        lineage  = rp(config["busco_lineage"])
    output:
        summary = os.path.join(rp(config["busco_outdir"]), "{sample}", "short_summary.specific.{sample}.txt"),
        done    = os.path.join(BUSCO_OUTDIR, "done", "{sample}_proteome_completeness.done"),
    params:
        busco_outdir = os.path.join(rp(config["busco_outdir"]), "{sample}"),
        busco_sif    = rp(config["busco_sif"]),
        busco_bind   = rp(config["busco_bind"]),
    threads: config.get("threads", 4)
    log:
        os.path.join(BUSCO_OUTDIR, "{sample}", "logs", "proteome_completeness.log"),
    shell:
        r"""
        set -euo pipefail
        exec > >(tee -a {log}) 2>&1
        mkdir -p "$(dirname {log})"
        echo "[$(date)] Attempting BUSCO for {wildcards.sample}"

        bash scripts/data_cleaning/busco.sh \
            "{input.proteome}" \
            "{params.busco_outdir}" \
            "{params.busco_bind}" \
            "{params.busco_sif}" \
            "{input.lineage}"

        cp {params.busco_outdir}/{wildcards.sample}/short_summary.specific.*.{wildcards.sample}.txt \
        {output.summary}

        echo "[$(date)] BUSCO completed for {wildcards.sample}"
        echo "done" > {output.done}
        exec 1>&- 2>&-
        wait

        """


############################################################################
# Rule: filter_proteomes  (aggregate — runs once, after all BUSCOs finish)
# Scripts:
#   - scripts/data_cleaning/sum_busco.sh
#   - scripts/data_cleaning/filter_busco.py
# Goal: combines all short_summary.specific.{sample}.txt files for each 
# sample into a combined tsv, then filter proteomes in study based on 
# relative completeness 
# Produces: 
#   - busco_summary.tsv
# Notes: 
#   - Filters out samples whose %_complete is an outlier on the low end,
#    defined as below Q1 - 1.5 * IQR (the standard Tukey fence).
#   - Removed FAAs are moved to a quarantine folder rather than deleted so
#    you can inspect them.
############################################################################
checkpoint filter_proteomes:
    input:
        summaries = expand(
            os.path.join(BUSCO_OUTDIR, "{sample}", "short_summary.specific.{sample}.txt"),
            sample=SAMPLES,
        ),
    output:
        tsv     = os.path.join(BUSCO_OUTDIR, "busco_summary.tsv"),
        kept    = os.path.join(BUSCO_OUTDIR, "kept_samples.txt"),
        removed = os.path.join(BUSCO_OUTDIR, "removed_samples.txt"),
    params:
        protein_dir    = PROTEIN_DIR,
        quarantine_dir = os.path.join(PROTEIN_DIR, "low_completeness"),
        min_complete   = config.get("busco_min_complete", 10),
    threads: 1
    log:
        os.path.join(BUSCO_OUTDIR, "logs", "filter_proteomes.log"),
    shell:
        r"""
        set -euo pipefail
        exec > >(tee -a {log}) 2>&1
        mkdir -p "$(dirname {log})"

        echo "[$(date)] Summarising BUSCO results"
        bash scripts/data_cleaning/sum_busco.sh {input.summaries} > {output.tsv}

        echo "[$(date)] Filtering proteomes by completeness (floor + IQR method)"
        python3 scripts/data_cleaning/filter_busco.py \
            "{output.tsv}" \
            "{params.protein_dir}" \
            "{params.quarantine_dir}" \
            "{output.kept}" \
            "{output.removed}" \
            "{params.min_complete}"

        echo "[$(date)] filter_proteomes done"
        exec 1>&- 2>&-
        wait
        """
############################################################################
# Rule: filter_short_seqs
# Goal: Remove sequences shorter than min_seq_len (default 50 aa) from each
#   filtered proteome before annotation. Short peptides cause spurious hits
#   in InterProScan, KofamScan, dbCAN, and DIAMOND.
# Produces:
#   - {filtered_faa_dir}/{sample}.faa
############################################################################
rule filter_short_seqs:
    input:
        proteome = os.path.join(PROTEIN_DIR, "{sample}.faa"),
        kept     = os.path.join(BUSCO_OUTDIR, "kept_samples.txt"),  # ensures checkpoint has run
    output:
        proteome = os.path.join(FILTERED_FAA_DIR, "{sample}.faa"),
    params:
        min_len = config.get("min_seq_len", 50),
    threads: 1
    log:
        os.path.join(OUTDIR, "logs", "filter_short_seqs", "{sample}.log"),
    conda:
        "envs/biopython.yaml"
    shell:
        r"""
        set -euo pipefail
        exec > >(tee -a {log}) 2>&1
        mkdir -p "$(dirname {log})"
        mkdir -p "{FILTERED_FAA_DIR}"

        echo "[$(date)] Filtering sequences shorter than {params.min_len} aa for {wildcards.sample}"

        seqkit seq -m {params.min_len} "{input.proteome}" > "{output.proteome}"

        BEFORE=$(grep -c "^>" "{input.proteome}" || true)
        AFTER=$(grep -c "^>" "{output.proteome}" || true)
        echo "[$(date)] {wildcards.sample}: $BEFORE sequences in, $AFTER sequences out"
        """

############################################################################
# Grouping key: get_taxonomy + add_clustering rules
# Scripts: 
#   - scripts/data_clustering/ncbi_taxonomy.py
#   - scripts/data_clustering/add_clustering.py
# Goal: create a clustering key that links sample id, species taxonomy, and
# aerotolerance and host association labels
# Produces: 
#   - nothing, if a grouping_key.tsv is already provided; copies provided
#   file into grouping_key.tsv
#   - if no file is provided; grouping_key.tsv
#
#   provided_grouping_key: ""            → run get_taxonomy + add_clustering
#   provided_grouping_key: "path/to/file" → copy that file to grouping_key.tsv
#
############################################################################
if config.get("provided_grouping_key"):

    rule use_provided_grouping_key:
        """Copy a pre-built grouping key, skipping get_taxonomy + add_clustering."""
        input:
            rp(config["provided_grouping_key"]),
        output:
            "grouping_key.tsv",
        log:
            os.path.join(OUTDIR, "logs", "use_provided_grouping_key.log"),
        shell:
            r"""
            set -euo pipefail
            exec > >(tee -a {log}) 2>&1
            echo "[$(date)] Using provided grouping key: {input}"
            cp "{input}" "{output}"
            echo "[$(date)] grouping_key.tsv ready"
            exec 1>&- 2>&-
            wait
            """

else:

    # -----------------------------------------------------------------
    # Rule: get_taxonomy
    # Script: scripts/data_clustering/ncbi_taxonomy.py
    # Goal: Get lineage information (genus → clade) from given species
    # list
    # Produces: taxonomy_raw.tsv
    # -----------------------------------------------------------------
    rule get_taxonomy:
        input:
            species = rp(config["species_list"]),
        output:
            tsv = "taxonomy_raw.tsv",
        log:
            os.path.join(OUTDIR, "logs", "get_taxonomy.log"),
        conda:
            "envs/ete3.yaml"
        shell:
            r"""
            set -euo pipefail
            exec > >(tee -a {log}) 2>&1

            echo "[$(date)] Fetching NCBI taxonomy for {input.species}"

            python3 scripts/data_clustering/ncbi_taxonomy.py \
                "{input.species}" \
                > "{output.tsv}"

            echo "[$(date)] get_taxonomy done — $(wc -l < {output.tsv}) rows written"
            exec 1>&- 2>&-
            wait
            """

    # -----------------------------------------------------------------
    # Rule: add_clustering
    # Script: scripts/data_clustering/add_clustering.py
    # Goal: annotate raw taxonomy TSV with: 
    #   host_association : HA (host-associated) or FL (free-living)
    #   aerotolerance    : aerobic or anaerobic
    #   sample_id        : matched from sample_key.tsv
    # Produces: final grouping_key.tsv used for downstream analyses.
    # -----------------------------------------------------------------
    rule add_clustering:
        input:
            tsv        = "taxonomy_raw.tsv",
            ha         = rp(config["HA_list"]),
            aerobic    = rp(config["aerobic_list"]),
            sample_key = rp(config["sample_key"]),
        output:
            tsv = "grouping_key.tsv",
        log:
            os.path.join(OUTDIR, "logs", "add_clustering.log"),
        shell:
            r"""
            set -euo pipefail
            exec > >(tee -a {log}) 2>&1

            echo "[$(date)] Adding host_association, aerotolerance, and sample_id columns"

            python3 scripts/data_clustering/add_clustering.py \
                "{input.tsv}" \
                "{input.ha}" \
                "{input.aerobic}" \
                "{input.sample_key}" \
                "{output.tsv}"

            echo "[$(date)] add_clustering done"
            exec 1>&- 2>&-
            wait
            """

############################################################################
# Rule: run_kofamscan 
# Script: 
#   - scripts/gene_annotation/kofam/kofam.sh
# Goal: run kofamscan on each proteome
# Produces: 
#    - "{sample}_ko_annotations.tsv"
############################################################################
rule run_kofamscan:
    input:
        proteome = os.path.join(FILTERED_FAA_DIR, "{sample}.faa"),

    output:
        annotations = os.path.join(rp(config["kofam_outdir"]), "{sample}", "{sample}_ko_annotations.tsv"),
        done        = os.path.join(KOFAM_OUTDIR, "{sample}", "run_kofamscan.done"),
    params:
        kofam_outdir = os.path.join(rp(config["kofam_outdir"]), "{sample}"),
        kofam_list  = rp(config["kofam_list"]),
        kofam_sif    = rp(config["kofam_sif"]),
        kofam_profiles   = rp(config["kofam_profiles"]),
    threads: config.get("threads", 4)

    log:
        os.path.join(KOFAM_OUTDIR, "{sample}", "logs", "run_kofamscan.log"),
    shell:
        r"""
        set -euo pipefail
        exec > >(tee -a {log}) 2>&1
        mkdir -p "$(dirname {log})"

        echo "[$(date)] Running KofamScan for {wildcards.sample}"

        bash scripts/gene_annotation/kofam/kofam.sh \
            "{input.proteome}" \
            "{params.kofam_outdir}" \
            "{params.kofam_sif}" \
            "{params.kofam_profiles}" \
            "{params.kofam_list}"
        
        echo "[$(date)] KofamScan completed for {wildcards.sample}"
        echo "done" > {output.done}
        exec 1>&- 2>&-
        wait
        """


############################################################################
# Rule: summarize_kofamscan 
# Script: 
#   - scripts/gene_annotation/kofam/kofam_summary.py
# Goal: Summarize all samples kofamscan output into one summary table
# Produces: 
#   - kofam_summary.tsv
# Note: 
#   - aggregate — runs once, after all samples finish
############################################################################
rule summarize_kofamscan:
    input:
        annotations = lambda wc: expand(
            os.path.join(KOFAM_OUTDIR, "{sample}", "{sample}_ko_annotations.tsv"),
            sample=get_filtered_samples(wc),
        ),
        kegg_pathway_ids           = rp(config["kegg_pathway_ids"]),
        kegg_orth_and_pathway_corr = rp(config["kegg_orth_and_pathway_corr"]),
    output:
        summary = os.path.join(KOFAM_OUTDIR, "kofam_summary.tsv"),
        done    = os.path.join(KOFAM_OUTDIR, "summarize_kofamscan.done"),
    threads: 1
    log:
        os.path.join(KOFAM_OUTDIR, "logs", "summarize_kofamscan.log"),
    shell:
        r"""
        set -euo pipefail
        exec > >(tee -a {log}) 2>&1
        mkdir -p "$(dirname {log})"

        echo "[$(date)] Summarizing KofamScan outputs"

        # Write one annotation file path per line for kofam_summary.py
        FILE_LIST=$(mktemp)
        trap "rm -f $FILE_LIST" EXIT
        printf '%s\n' {input.annotations} > "$FILE_LIST"

        python3 scripts/gene_annotation/kofam/kofam_summary.py \
            "{input.kegg_pathway_ids}" \
            "{input.kegg_orth_and_pathway_corr}" \
            "$FILE_LIST" \
            > "{output.summary}"

        echo "[$(date)] summarize_kofamscan done"
        echo "done" > {output.done}
        exec 1>&- 2>&-
        wait
        """




############################################################################
# Rule: pathway_summary_stats
# Script: 
#   - scripts/gene_annotation/kofam/pathways_summary_stats.py
# Goal: Joins kofam_summary.tsv with kegg_pathway_hierarchy.tsv and 
# grouping_key.tsv to compute per-pathway min/max/mean completion rates 
# for HA vs FL and aerobic vs anaerobic, with Mann-Whitney U test + BH 
# correction.
# Produces: 
#   - pathway_sumstats.tsv
# Note: 
#   - Runs after summarize_kofamscan; its output feeds pull_kofam_pathways.
############################################################################
rule pathway_summary_stats:
    input:
        kofam_summary = os.path.join(KOFAM_OUTDIR, "kofam_summary.tsv"),
        pathway_hier  = rp(config["pathway_hierarchy"]),
        grouping_key  = ancient("grouping_key.tsv"),
    output:
        tsv  = os.path.join(KOFAM_OUTDIR, "pathway_sumstats.tsv"),
        done = os.path.join(KOFAM_OUTDIR, "pathway_summary_stats.done"),
    params:
        alpha = config.get("stats_alpha", 0.05),
    threads: 2
    log:
        os.path.join(KOFAM_OUTDIR, "logs", "pathway_summary_stats.log"),
    conda:
        "envs/scipy.yaml"
    shell:
        r"""
        set -euo pipefail
        exec > >(tee -a {log}) 2>&1
        mkdir -p "$(dirname {log})"

        echo "[$(date)] Computing pathway group statistics"

        python3 scripts/gene_annotation/kofam/pathways_summary_stats.py \
            --kofam-summary "{input.kofam_summary}" \
            --pathway-hier  "{input.pathway_hier}" \
            --grouping-key  "{input.grouping_key}" \
            --output        "{output.tsv}" \
            --alpha         {params.alpha}

        echo "[$(date)] pathway_summary_stats done"
        echo "done" > {output.done}
        exec 1>&- 2>&-
        wait
        """

############################################################################
# Rule: pull_kofam_pathways
# Script: 
#   - scripts/gene_annotation/kofam/pull_pathways.sh
# Goal: Filters pathway_sumstats.tsv for three primary metabolite 
# pathway groups:
#   - Metabolism of cofactors and vitamins
#   - Amino acid metabolism
#   - Lipid metabolism
# keeping only rows where host_association_sig == "yes", sorted by MW stat 
# descending, top 10 per group. 
# Produces: 
#   - key_pathways.tsv
# Note: 
#   - An optional pathway list forces additional pathways in
#   regardless of significance (set forced_pathways: "" in config to skip).
#   - Can change to look at aerotolerance split instead of HA by changing
#   column numbers in pull_pathways.sh
############################################################################
checkpoint pull_kofam_pathways:
    input:
        summary = os.path.join(KOFAM_OUTDIR, "pathway_sumstats.tsv"),
        **({
            "forced_pathways": rp(config["forced_pathways"])
        } if config.get("forced_pathways") else {}),
    output:
        pathway_list = os.path.join(KOFAM_OUTDIR, "key_pathways.tsv"),
        done         = os.path.join(KOFAM_OUTDIR, "pull_kofam_pathways.done"),
    params:
        forced = config.get("forced_pathways", ""),
    threads: 1
    log:
        os.path.join(KOFAM_OUTDIR, "logs", "pull_kofam_pathways.log"),
    shell:
        r"""
        set -euo pipefail
        exec > >(tee -a {log}) 2>&1
        mkdir -p "$(dirname {log})"

        echo "[$(date)] Pulling significant KofamScan pathways"

        FORCED="{params.forced}"

        bash scripts/gene_annotation/kofam/pull_pathways.sh \
            "{input.summary}" \
            "{output.pathway_list}" \
            $([ -n "$FORCED" ] && echo "$FORCED" || true)

        echo "[$(date)] pull_kofam_pathways done"
        echo "done" > {output.done}
        exec 1>&- 2>&-
        wait
        """

############################################################################
# Rule: assess_modules
# Scripts: 
#   - scripts/gene_annotation/kofam/KEGG_modules/module_parser.py
#   - scripts/gene_annotation/kofam/KEGG_modules/join_tsv.py
# Goal: For each module in key_pathways.tsv, pull its boolean KO 
# definition from pathway_module_ko.tsv, write its definition to a temp file, 
# run module_parser.py for every sample annotation to determine functionality 
# of that module in a given sample, then join the functionality results with 
# grouping_key.tsv
# Produces:    
#   - {module}_functionality.tsv
# Note: 
#   - per module, aggregate over all samples
############################################################################
rule assess_modules:
    input:
        pathway_module_ko = rp(config["pathway_module_ko"]),
        key_pathways      = os.path.join(KOFAM_OUTDIR, "key_pathways.tsv"),
        grouping_key      = ancient("grouping_key.tsv"),
        annotations       = lambda wc: expand(
            os.path.join(rp(config["kofam_outdir"]), "{sample}", "{sample}_ko_annotations.tsv"),
            sample=get_filtered_samples(wc),
        ),
    output:
        tsv  = os.path.join(KOFAM_OUTDIR, "module_results", "{module}_functionality.tsv"),
        done = os.path.join(KOFAM_OUTDIR, "{module}_assess_modules.done"),
    threads: 1
    log:
        os.path.join(KOFAM_OUTDIR, "logs", "{module}_assess_modules.log"),
    shell:
        r"""
        set -euo pipefail
        exec > >(tee -a {log}) 2>&1
        mkdir -p "$(dirname {log})"

        echo "[$(date)] Assessing module {wildcards.module} for all samples"

        mkdir -p "$(dirname {output.tsv})"

        # Extract the boolean definition for this module from col 5 of pathway_module_ko.tsv
        # (col 3 = module_id, col 5 = module_definition)
        MODULE_DEF=$(awk -F'\t' -v mod="{wildcards.module}" '$3 == mod {{print $5; exit}}' "{input.pathway_module_ko}")

        if [ -z "$MODULE_DEF" ]; then
            echo "ERROR: {wildcards.module} not found in {input.pathway_module_ko}"
            exit 1
        fi

        echo "Definition: $MODULE_DEF"

        # Write definition to a temp file for module_parser.py
        TMP_DEF=$(mktemp)
        TMP_FUNC=$(mktemp)
        trap "rm -f $TMP_DEF $TMP_FUNC" EXIT
        echo "$MODULE_DEF" > "$TMP_DEF"

        # Write functionality header to temp file
        echo -e "sample_id\tfunctional\tmissing_kos\tpresent_kos" > "$TMP_FUNC"

        # Run module_parser.py for each sample and append one row per sample
        for annotation in {input.annotations}; do
            sample=$(basename "$annotation" _ko_annotations.tsv)
            echo "  Processing $sample..."
            result=$(python3 scripts/gene_annotation/kofam/KEGG_modules/module_parser.py -1 "$annotation" -2 "$TMP_DEF")
            echo -e "$sample\t$result" >> "$TMP_FUNC"
        done

        # Left-join functionality results with grouping_key.tsv on sample_id
        echo "[$(date)] Joining with grouping_key.tsv..."
        python3 scripts/gene_annotation/kofam/KEGG_modules/join_tsv.py \
            "$TMP_FUNC" \
            "{input.grouping_key}" \
            "sample_id" \
            "sample_id" \
            "{output.tsv}"

        echo "[$(date)] Done — $(wc -l < {output.tsv}) rows in {output.tsv}"
        echo "done" > {output.done}
        exec 1>&- 2>&-
        wait
        """

############################################################################
# Rule: module_stats  (aggregate — runs once, after all assess_modules finish)
# Script: 
#   - scripts/gene_annotation/kofam/KEGG_modules/module_stats.py
# Goal: 
#   - For each module x grouping column (host_association, aerotolerance, phylum,
#   etc.), test whether functional rates differ significantly across groups:
#       - 2 groups (Fisher's exact test (odds ratio))
#       - >2 groups (Chi-squared test)
# Produces:
#   - module_stats.tsv
#       - where one row per module x significant grouping column combination
# Note: 
#   - Alpha (significance cutoff) is set in config and passed as a parameter.
############################################################################
rule module_stats:
    input:
        module_files      = lambda wc: expand(
            os.path.join(KOFAM_OUTDIR, "module_results", "{module}_functionality.tsv"),
            module=get_modules(wc),
        ),
        pathway_module_ko = rp(config["pathway_module_ko"]),
    output:
        tsv  = os.path.join(KOFAM_OUTDIR, "module_stats.tsv"),
        done = os.path.join(KOFAM_OUTDIR, "module_stats.done"),
    params:
        alpha        = config.get("stats_alpha", 0.05),
        top_kos      = config.get("top_kos", 5),
        grouping_cols = " ".join(config.get("grouping_cols",
            ["host_association", "aerotolerance", "phylum", "class", "order", "family"])),
    threads: 1
    log:
        os.path.join(KOFAM_OUTDIR, "logs", "module_stats.log"),
    conda:
        "envs/scipy.yaml"
    shell:
        r"""
        set -euo pipefail
        exec > >(tee -a {log}) 2>&1
        mkdir -p "$(dirname {log})"

        echo "[$(date)] Running module statistics (alpha={params.alpha})"

        # If no modules passed significance thresholds, key_pathways.tsv is empty
        # and input.module_files will be empty — guard against calling the script
        # with no arguments, which causes an argparse error.
        MODULE_FILES="{input.module_files}"
        if [ -z "$MODULE_FILES" ]; then
            echo "[$(date)] No significant modules found — writing empty module_stats.tsv"
            echo -e "module_id\tmodule_name\ttest\tstatistic\tpvalue\tpadj" > "{output.tsv}"
            echo "done" > "{output.done}"
            exec 1>&- 2>&-
            wait
            exit 0
        fi

        python3 scripts/gene_annotation/kofam/KEGG_modules/module_stats.py \
            --module-files {input.module_files} \
            --pathway-module-ko "{input.pathway_module_ko}" \
            --output "{output.tsv}" \
            --alpha {params.alpha} \
            --top-kos {params.top_kos} \
            --grouping-cols {params.grouping_cols}

        echo "[$(date)] module_stats done"
        echo "done" > {output.done}
        exec 1>&- 2>&-
        wait
        """

############################################################################
# Rule: run_interproscan 
# Script:   
#   - scripts/gene_annotation/interproscan/interproscan.sh
# Goal: 
#   - Run interproscan for gene annotation on each sample
# Produces: 
#   - {sample}_interproscan.tsv
############################################################################
rule run_interproscan:
    input:
        proteome = os.path.join(FILTERED_FAA_DIR, "{sample}.faa"),
    output:
        annotations = os.path.join(rp(config["interproscan_outdir"]), "{sample}", "{sample}_interproscan.tsv"),
        done        = os.path.join(INTERPRO_OUTDIR, "{sample}", "run_interproscan.done"),
    params:
        interproscan_outdir = os.path.join(rp(config["interproscan_outdir"]), "{sample}"),
        interpro    = rp(config["interpro"]),
    threads: 8   # matches cpus_per_task in profiles/slurm/config.yaml run_interproscan
    log:
        os.path.join(INTERPRO_OUTDIR, "{sample}", "logs", "run_interproscan.log"),
    conda:
        "envs/interproscan.yaml"
    shell:
        r"""
        set -euo pipefail
        exec > >(tee -a {log}) 2>&1
        mkdir -p "$(dirname {log})"

        echo "[$(date)] Running Interproscan for {wildcards.sample}"

        bash scripts/gene_annotation/interproscan/interproscan.sh \
            "{input.proteome}" \
            "{params.interproscan_outdir}" \
            "{params.interpro}" \
            "{threads}"

        echo "[$(date)] Interproscan completed for {wildcards.sample}"
        echo "done" > {output.done}
        exec 1>&- 2>&-
        wait
        """
    
############################################################################
# Rule: clean_interproscan 
# Script: 
#   - scripts/gene_annotation/interproscan/clean_interpro.py
# Goal: To reduce redundancy in interproscan outputs 
# Produces: 
#   - {sample}_interprosum.tsv
############################################################################
rule clean_interproscan:
    input:
        proteome = os.path.join(FILTERED_FAA_DIR, "{sample}.faa"),
        annotations = os.path.join(rp(config["interproscan_outdir"]), "{sample}", "{sample}_interproscan.tsv"),
    output:
        annotations = os.path.join(rp(config["interproscan_outdir"]), "{sample}", "{sample}_interprosum.tsv"),
        done        = os.path.join(INTERPRO_OUTDIR, "{sample}", "clean_interproscan.done"),
    params:
        interproscan_outdir = os.path.join(rp(config["interproscan_outdir"]), "{sample}"),
        interpro    = rp(config["interpro"]),
    threads: 1
    log:
        os.path.join(INTERPRO_OUTDIR, "{sample}", "logs", "clean_interproscan.log"),
    conda:
        "envs/biopython.yaml"
    shell:
        r"""
        set -euo pipefail
        exec > >(tee -a {log}) 2>&1
        mkdir -p "$(dirname {log})"

        echo "[$(date)] Cleaning Interproscan for {wildcards.sample}"

        python3 scripts/gene_annotation/interproscan/clean_interpro.py \
            "{input.proteome}" \
            "{input.annotations}" \
            > "{output.annotations}" 
        
        echo "[$(date)] Interproscan Cleaning completed for {wildcards.sample}"
        echo "done" > {output.done}
        exec 1>&- 2>&-
        wait
        """

############################################################################
# Rule: summarize_interproscan
# Script: 
#    - scripts/gene_annotation/interproscan/interpro_summary.py
# Goal: 
#   - Combine per-species interproscan summary tsvs into one
#   multi-species table, annotated with grouping_key.tsv
# Produces: 
#   - all_interpro_summary.tsv
# Note:
#   - aggregate — runs once all samples are done
############################################################################
rule summarize_interproscan:
    input:
        grouping_key = ancient("grouping_key.tsv"),
        ipr          = rp(config["ipr_entry_list"]),
        tsv_files    = lambda wc: expand(
            os.path.join(INTERPRO_OUTDIR, "{sample}", "{sample}_interprosum.tsv"),
            sample=get_filtered_samples(wc),
        ),
        faa_files    = lambda wc: expand(
            os.path.join(FILTERED_FAA_DIR, "{sample}.faa"),
            sample=get_filtered_samples(wc),
        ),
    output:
        summary = os.path.join(INTERPRO_OUTDIR, "all_interpro_summary.tsv"),
        done    = os.path.join(INTERPRO_OUTDIR, "summarize_interproscan.done"),
    threads: 1
    log:
        os.path.join(INTERPRO_OUTDIR, "logs", "summarize_interproscan.log"),
    shell:
        r"""
        set -euo pipefail
        exec > >(tee -a {log}) 2>&1
        mkdir -p "$(dirname {log})"

        echo "[$(date)] Summarizing InterProScan results"

        # Write file lists to temp files for the Python script
        TSV_LIST=$(mktemp)
        FAA_LIST=$(mktemp)
        trap "rm -f $TSV_LIST $FAA_LIST" EXIT

        printf '%s\n' {input.tsv_files} > "$TSV_LIST"
        printf '%s\n' {input.faa_files} > "$FAA_LIST"

        python3 scripts/gene_annotation/interproscan/interpro_summary.py \
            --tsv-list     "$TSV_LIST" \
            --faa-list     "$FAA_LIST" \
            --grouping-key "{input.grouping_key}" \
            --ipr          "{input.ipr}" \
            --output       "{output.summary}"

        echo "[$(date)] summarize_interproscan done"
        echo "done" > {output.done}
        exec 1>&- 2>&-
        wait
        """

############################################################################
# Rule: interpro_sum_stats  (per taxonomic group)
# Script: 
#   - scripts/gene_annotation/interproscan/interpro_sum_stats.py
# Goal: 
#   - to calculate summary statistics for each taxonomic group 
#   (config["interpro_groups"]) on host_association or aerotolerance using
# ttests and mann whitney 
# Produces: 
#   - {group}_ttest.tsv
#   - {group}_mannwhitney.tsv.
############################################################################
INTERPRO_STATS_OUTDIR = os.path.join(INTERPRO_OUTDIR, "stats")

rule interpro_sum_stats:
    input:
        summary = os.path.join(INTERPRO_OUTDIR, "all_interpro_summary.tsv"),
    output:
        ttest  = os.path.join(INTERPRO_STATS_OUTDIR, "{group}_ttest.tsv"),
        mw     = os.path.join(INTERPRO_STATS_OUTDIR, "{group}_mannwhitney.tsv"),
        done   = os.path.join(INTERPRO_OUTDIR, "{group}_interpro_sum_stats.done"),
    params:
        group_col   = lambda wc: INTERPRO_GROUP_COL_MAP.get(wc.group, "phylum"),
        compare_col = config.get("interpro_compare_col", "host_association"),
        out_prefix  = os.path.join(INTERPRO_STATS_OUTDIR, "{group}"),
    threads: 1
    log:
        os.path.join(INTERPRO_OUTDIR, "logs", "{group}_interpro_sum_stats.log"),
    conda:
        "envs/scipy.yaml"
    shell:
        r"""
        set -euo pipefail
        exec > >(tee -a {log}) 2>&1
        mkdir -p "$(dirname {log})"

        echo "[$(date)] Running InterPro stats for group '{wildcards.group}'"

        mkdir -p "{INTERPRO_STATS_OUTDIR}"

        python3 scripts/gene_annotation/interproscan/interpro_sum_stats.py \
            --input       "{input.summary}" \
            --group       "{wildcards.group}" \
            --group-col   "{params.group_col}" \
            --compare-col "{params.compare_col}" \
            --output      "{params.out_prefix}"

        echo "[$(date)] interpro_sum_stats done for '{wildcards.group}'"
        echo "done" > {output.done}
        exec 1>&- 2>&-
        wait
        """


############################################################################
# Rule: merge_interpro_stats  (aggregate — runs once all interpro_sum_stats finish)
# Script: 
#   - scripts/gene_annotation/interproscan/all_interpro_sum_stats.py
# Goal: 
#   -  Merge per-lineage stat TSVs into one wide-format table for comparison
#   across lineages
# Produces: 
#   - all_lineages_mannwhitney.tsv
# Note: 
#   - config["interpro_stat_type"] : "mannwhitney" or "ttest" (default: mannwhitney)
#   - rp(config["ipr_compare_list"])   : text file with IPR IDs to include (one per line)
############################################################################
rule merge_interpro_stats:
    input:
        stat_files = expand(
            os.path.join(INTERPRO_OUTDIR, "stats",
                         "{group}_" + config.get("interpro_stat_type", "mannwhitney") + ".tsv"),
            group=INTERPRO_GROUPS,
        ),
        ipr_list = rp(config["ipr_compare_list"]),
    output:
        summary = os.path.join(INTERPRO_OUTDIR, "stats", "all_lineages_"
                               + config.get("interpro_stat_type", "mannwhitney") + ".tsv"),
        done    = os.path.join(INTERPRO_OUTDIR, "merge_interpro_stats.done"),
    params:
        labels = ",".join(INTERPRO_GROUPS),
    threads: 1
    log:
        os.path.join(INTERPRO_OUTDIR, "logs", "merge_interpro_stats.log"),
    shell:
        r"""
        set -euo pipefail
        exec > >(tee -a {log}) 2>&1
        mkdir -p "$(dirname {log})"

        echo "[$(date)] Merging per-lineage InterPro stats"

        TSV_LIST=$(mktemp)
        trap "rm -f $TSV_LIST" EXIT

        printf '%s\n' {input.stat_files} > "$TSV_LIST"

        python3 scripts/gene_annotation/interproscan/all_interpro_sum_stats.py \
            --tsv_list "$TSV_LIST" \
            --ipr_list "{input.ipr_list}" \
            --labels   "{params.labels}" \
            --output   "{output.summary}"

        echo "[$(date)] merge_interpro_stats done"
        echo "done" > {output.done}
        exec 1>&- 2>&-
        wait
        """
############################################################################
# Rule: run_dbcan
# Script: 
#   - scripts/gene_annotation/dbcan/dbcan.sh
# Goal: to annoate CAZyme families, and integrates Cazyme Gene Clusters (CGCs) 
#  and substrate predictions.
# Produces: 
#   - {sample}_overview.txt
############################################################################
rule run_dbcan:
    input:
        proteome = os.path.join(FILTERED_FAA_DIR, "{sample}.faa"),
    output:
        summary = os.path.join(rp(config["dbcan_outdir"]), "{sample}", "{sample}_overview.txt"),
        done    = os.path.join(DBCAN_OUTDIR, "{sample}", "run_dbcan.done"),
    params:
        dbcan_outdir = os.path.join(rp(config["dbcan_outdir"]), "{sample}"),
        dbcan_sif    = rp(config["dbcan_sif"]),
        dbcan_bind   = rp(config["dbcan_bind"]),
        dbcan_db     = rp(config["dbcan_db"]),
    threads: config.get("threads", 4)

    log:
        os.path.join(DBCAN_OUTDIR, "{sample}", "logs", "run_dbcan.log"),
    shell:
        r"""
        set -euo pipefail
        exec > >(tee -a {log}) 2>&1
        mkdir -p "$(dirname {log})"
        echo "[$(date)] Attempting dbCAN for {wildcards.sample}"

        bash scripts/gene_annotation/dbcan/dbcan.sh \
            "{input.proteome}" \
            "{params.dbcan_outdir}" \
            "{params.dbcan_bind}" \
            "{params.dbcan_sif}" \
            "{params.dbcan_db}"

        echo "[$(date)] dbCAN completed for {wildcards.sample}"
        echo "done" > {output.done}
        exec 1>&- 2>&-
        wait

        """

############################################################################
# Rule: clean_dbcan 
# Script: 
#   - scripts/gene_annotation/dbcan/clean_dbcan.py
# Goal: processes dbCAN overview files and collapses predictions 
# from 3 tools to one summary prediction
# Produces:
#   - {sample}_dbcan_summary.txt
# Note: 
#   - Any error disagreements are deleted
############################################################################
rule clean_dbcan:
    input:
        annotations = os.path.join(rp(config["dbcan_outdir"]), "{sample}", "{sample}_overview.txt"),
    output:
        summary = os.path.join(rp(config["dbcan_outdir"]), "{sample}", "{sample}_dbcan_summary.txt"),
        done        = os.path.join(DBCAN_OUTDIR, "{sample}", "clean_dbcan.done"),
    threads: 1
    log:
        os.path.join(DBCAN_OUTDIR, "{sample}", "logs", "clean_dbcan.log"),
    shell:
        r"""
        set -euo pipefail
        exec > >(tee -a {log}) 2>&1
        mkdir -p "$(dirname {log})"

        echo "[$(date)] Cleaning dbCAN for {wildcards.sample}"

        python3 scripts/gene_annotation/dbcan/clean_dbcan.py \
            "{input.annotations}" \
            > "{output.summary}"
        
        # remove error disagreements
        sed -i '/ERROR/d' "{output.summary}"
        
        echo "[$(date)] dbCAN cleaning completed for {wildcards.sample}"
        echo "done" > {output.done}
        exec 1>&- 2>&-
        wait
        """


############################################################################
# Rule: summarize_dbcan  
# Script: 
#   - scripts/gene_annotation/dbcan/combine_dbcan.py
# Goal: Combines per-species dbCAN CAZyme summary TXTs into one multi-species
# table, annotated with grouping_key.tsv.
# Produces: 
#   - dbcan_combined.tsv
# Note: 
#   - aggregate — runs once all dbCAN jobs finish
############################################################################
rule summarize_dbcan:
    input:
        grouping_key = ancient("grouping_key.tsv"),
        txt_files    = lambda wc: expand(
            os.path.join(DBCAN_OUTDIR, "{sample}", "{sample}_dbcan_summary.txt"),
            sample=get_filtered_samples(wc),
        ),
        faa_files    = lambda wc: expand(
            os.path.join(FILTERED_FAA_DIR, "{sample}.faa"),
            sample=get_filtered_samples(wc),
        ),
    output:
        summary = os.path.join(DBCAN_OUTDIR, "dbcan_combined.tsv"),
        done    = os.path.join(DBCAN_OUTDIR, "summarize_dbcan.done"),
    threads: 1
    log:
        os.path.join(DBCAN_OUTDIR, "logs", "summarize_dbcan.log"),
    shell:
        r"""
        set -euo pipefail
        exec > >(tee -a {log}) 2>&1
        mkdir -p "$(dirname {log})"

        echo "[$(date)] Combining dbCAN CAZyme results"

        TXT_LIST=$(mktemp)
        FAA_LIST=$(mktemp)
        trap "rm -f $TXT_LIST $FAA_LIST" EXIT

        printf '%s\n' {input.txt_files} > "$TXT_LIST"
        printf '%s\n' {input.faa_files} > "$FAA_LIST"

        python3 scripts/gene_annotation/dbcan/combine_dbcan.py \
            --txt-list     "$TXT_LIST" \
            --faa-list     "$FAA_LIST" \
            --grouping-key "{input.grouping_key}" \
            --output       "{output.summary}"

        echo "[$(date)] summarize_cazyme done"
        echo "done" > {output.done}
        exec 1>&- 2>&-
        wait
        """


############################################################################
# Rule: dbcan_stats
# Script: 
#   - scripts/gene_annotation/dbcan/dbcan_stats.py
# Goal: to calculate summary statistics for each taxonomic group 
#   (config["dbcan_groups"]) on host_association or aerotolerance using
#   ttests and mann whitney 
# Produces: 
#   - {group}_dbcan_stats.tsv
#           - per taxonomic group — one job per group
############################################################################
DBCAN_STATS_OUTDIR = os.path.join(DBCAN_OUTDIR, "stats")

rule dbcan_stats:
    input:
        combined = os.path.join(DBCAN_OUTDIR, "dbcan_combined.tsv"),
    output:
        tsv  = os.path.join(DBCAN_STATS_OUTDIR, "{group}_dbcan_stats.tsv"),
        done = os.path.join(DBCAN_OUTDIR, "{group}_dbcan_stats.done"),
    params:
        group_col   = lambda wc: DBCAN_GROUP_COL_MAP.get(wc.group, "phylum"),
        compare_col = config.get("dbcan_compare_col", "host_association"),
        metric      = config.get("dbcan_metric",      "pct_genes"),
        alpha       = config.get("stats_alpha",        0.05),
    threads: 1
    log:
        os.path.join(DBCAN_OUTDIR, "logs", "{group}_dbcan_stats.log"),
    conda:
        "envs/scipy_numpy.yaml"
    shell:
        r"""
        set -euo pipefail
        exec > >(tee -a {log}) 2>&1
        mkdir -p "$(dirname {log})"

        echo "[$(date)] Running dbCAN stats for group '{wildcards.group}'"

        mkdir -p "{DBCAN_STATS_OUTDIR}"

        python3 scripts/gene_annotation/dbcan/dbcan_stats.py \
            --input       "{input.combined}" \
            --output      "{output.tsv}" \
            --group       "{wildcards.group}" \
            --group-col   "{params.group_col}" \
            --compare-col "{params.compare_col}" \
            --metric      "{params.metric}" \
            --alpha       {params.alpha}

        echo "[$(date)] dbcan_stats done for '{wildcards.group}'"
        echo "done" > {output.done}
        exec 1>&- 2>&-
        wait
        """

############################################################################
# Rule: merge_dbcan_stats  (aggregate — runs once all dbcan_stats finish)
# Script: 
#   - scripts/gene_annotation/dbcan/all_dbcan_stats.py
# Goal: 
#   - Merges per-lineage dbCAN stats TSVs into one wide-format table for
#   comparison across lineages. (Mirrors merge_interpro_stats)
# Produces: 
#   - all_lineages_dbcan_stats.tsv
# Note: 
#   = config["dbcan_cazyme_list"] : optional text file with CAZyme IDs to 
#     include; leave empty ("") to include all families
############################################################################
rule merge_dbcan_stats:
    input:
        stat_files = expand(
            os.path.join(DBCAN_OUTDIR, "stats", "{group}_dbcan_stats.tsv"),
            group=DBCAN_GROUPS,
        ),
        **({
            "cazyme_list": config["dbcan_cazyme_list"]
        } if config.get("dbcan_cazyme_list") else {}),
    output:
        summary = os.path.join(DBCAN_OUTDIR, "stats", "all_lineages_dbcan_stats.tsv"),
        done    = os.path.join(DBCAN_OUTDIR, "merge_dbcan_stats.done"),
    params:
        labels      = ",".join(DBCAN_GROUPS),
        cazyme_list = config.get("dbcan_cazyme_list", ""),
    threads: 1
    log:
        os.path.join(DBCAN_OUTDIR, "logs", "merge_dbcan_stats.log"),
    shell:
        r"""
        set -euo pipefail
        exec > >(tee -a {log}) 2>&1
        mkdir -p "$(dirname {log})"

        echo "[$(date)] Merging per-lineage dbCAN stats"

        TSV_LIST=$(mktemp)
        trap "rm -f $TSV_LIST" EXIT

        printf '%s\n' {input.stat_files} > "$TSV_LIST"

        CAZYME_ARG=""
        if [ -n "{params.cazyme_list}" ]; then
            CAZYME_ARG="--cazyme-list {params.cazyme_list}"
        fi

        python3 scripts/gene_annotation/dbcan/all_dbcan_stats.py \
            --tsv-list "$TSV_LIST" \
            --labels   "{params.labels}" \
            --output   "{output.summary}" \
            $CAZYME_ARG

        echo "[$(date)] merge_dbcan_stats done"
        echo "done" > {output.done}
        exec 1>&- 2>&-
        wait
        """

############################################################################
# Rule: run_orthofinder
# Script: 
#   -  scripts/comp_genomics/orthofinder/orthofinder.sh
# Goal: To run orthofinder for orthogroup formation on whole dataset
# Modes (set via config orthofinder_mode):
#   - fresh: wipes existing output and runs from scratch
#   - add:   adds new samples to a previous run via -b flag
#            (requires orthofinder_prev_results in config)
# Produces (alot; these are key products): 
#   - Orthogroups.GeneCount.tsv
#   - Orthogroups.tsv
############################################################################
rule run_orthofinder:
    input:
        proteomes = lambda wc: expand(
            os.path.join(FILTERED_FAA_DIR, "{sample}.faa"),
            sample=get_filtered_samples(wc),
        ),
    output:
        genecounts  = os.path.join(rp(config["orthofinder_outdir"]), "Results_here", "Orthogroups", "Orthogroups.GeneCount.tsv"),
        orthogroups = os.path.join(rp(config["orthofinder_outdir"]), "Results_here", "Orthogroups", "Orthogroups.tsv"),
        done        = os.path.join(ORTHO_OUTDIR, "run_orthofinder.done"),
    threads: config.get("threads", 24),
    params:
        orthofinder_outdir = ORTHO_OUTDIR,
        protein_dir        = PROTEIN_DIR,
        mode               = config.get("orthofinder_mode", "fresh"),
        prev_results       = rp(config.get("orthofinder_prev_results", "")),
    log:
        os.path.join(p("outdir"), "comparative_genomics", "logs", "run_orthofinder.log"),
    conda:
        "envs/orthofinder.yaml"
    shell:
        r"""
        set -euo pipefail

        MODE="{params.mode}"
        PREV="{params.prev_results}"

        # Fresh mode: wipe output BEFORE the tee redirect is open,
        # so rm -rf doesn't delete the log directory we're writing to.
        if [ "$MODE" = "fresh" ]; then
            REAL_OUTDIR=$(realpath -m "{params.orthofinder_outdir}" 2>/dev/null \
                          || echo "{params.orthofinder_outdir}")
            # A previous run's Results_here (the baseline "add" mode builds
            # from) lives inside this same directory. Back it up before
            # wiping so a fresh run -- deliberate or from a stale config on
            # a retry -- can never destroy the only copy of prior results.
            if [ -d "$REAL_OUTDIR/Results_here" ]; then
                BACKUP_DIR="${{REAL_OUTDIR}}_backup_$(date +%Y%m%d_%H%M%S)"
                mkdir -p "$BACKUP_DIR"
                mv "$REAL_OUTDIR/Results_here" "$BACKUP_DIR/Results_here"
            fi
            rm -rf "$REAL_OUTDIR"
        fi

        # (Re)create the log directory, then open the tee redirect.
        mkdir -p "$(dirname {log})"
        exec > >(tee -a {log}) 2>&1

        if [ "$MODE" = "fresh" ]; then
            echo "[$(date)] Fresh run — previous output removed"
            bash scripts/comp_genomics/orthofinder/orthofinder.sh \
                "{params.protein_dir}" \
                "{params.orthofinder_outdir}" \
                "{threads}" \
                ""

        elif [ "$MODE" = "add" ]; then
            if [ -z "$PREV" ]; then
                echo "ERROR: orthofinder_mode=add requires orthofinder_prev_results to be set in config"
                exit 1
            fi
            echo "[$(date)] Adding new samples using previous results at $PREV"
            bash scripts/comp_genomics/orthofinder/orthofinder.sh \
                "{params.protein_dir}" \
                "{params.orthofinder_outdir}" \
                "{threads}" \
                "$PREV"

        else
            echo "ERROR: Unknown orthofinder_mode '$MODE'. Must be 'fresh' or 'add'."
            exit 1
        fi

        echo "[$(date)] OrthoFinder completed."
        echo "done" > {output.done}
        exec 1>&- 2>&-
        wait
        """


############################################################################
# Rule: orthogroup_stats
# Script: 
#   - scripts/comp_genomics/orthofinder/orthogroup_stats.py
# Goal:
#   - Compares orthogroup gene counts across groupings using:
#       -  Mann-Whitney with BH correction (2 groups)
#       - or Kruskal-Wallis with BH correction (>2 groups)
# Produces: 
#   - orthogroups_{col}.tsv
# Note: 
#   - Optionally takes a list of genes (ie like cadherin.txt) 
#   and will pull the orthogroups associated with those genes and 
#   will calculate stats for those orthogroups 
############################################################################
############################################################################
# Rule: extract_domain_genes  (optional — runs before orthogroup_stats)
# Script:
#   - scripts/gene_annotation/interproscan/extract_domain_genes.py
# Goal: Find genes carrying InterPro domains that were significant in
# merge_interpro_stats' output (all_lineages_{stat}.tsv), by scanning
# per-sample {sample}_interprosum.tsv files for those domain IDs.
# Produces:
#   - a gene list file (one gene ID per line) — the same format
#   orthogroup_stats.py's --genes-of-interest expects, so its output feeds
#   directly into that rule alongside (or instead of) a manually curated list.
# Note:
#   - Optional: set config["domain_significant_genes_outfile"] (a path) to
#   enable; leave "" (default) to skip.
#   - config["domain_significant_in"]: "any" (default, significant in at
#   least one lineage) or "all" (significant in every lineage).
############################################################################
if config.get("domain_significant_genes_outfile"):

    rule extract_domain_genes:
        input:
            stats_tsv = os.path.join(INTERPRO_OUTDIR, "stats", "all_lineages_"
                         + config.get("interpro_stat_type", "mannwhitney") + ".tsv"),
            interpro_files = lambda wc: expand(
                os.path.join(INTERPRO_OUTDIR, "{sample}", "{sample}_interprosum.tsv"),
                sample=get_filtered_samples(wc),
            ),
        output:
            genes = rp(config["domain_significant_genes_outfile"]),
            done  = os.path.join(ORTHO_OUTDIR, "extract_domain_genes.done"),
        params:
            significant_in = config.get("domain_significant_in", "any"),
        threads: 1
        log:
            os.path.join(ORTHO_OUTDIR, "logs", "extract_domain_genes.log"),
        shell:
            r"""
            set -euo pipefail
            exec > >(tee -a {log}) 2>&1
            mkdir -p "$(dirname {log})"
            mkdir -p "$(dirname {output.genes})"

            echo "[$(date)] Extracting genes for significant InterPro domains"

            TSV_LIST=$(mktemp)
            trap "rm -f $TSV_LIST" EXIT
            printf '%s\n' {input.interpro_files} > "$TSV_LIST"

            python3 scripts/gene_annotation/interproscan/extract_domain_genes.py \
                --stats-tsv         "{input.stats_tsv}" \
                --interpro-tsv-list "$TSV_LIST" \
                --output            "{output.genes}" \
                --significant-in    {params.significant_in}

            echo "[$(date)] extract_domain_genes done"
            echo "done" > {output.done}
            exec 1>&- 2>&-
            wait
            """

############################################################################
# Rule: orthogroup_stats
# Script:
#   - scripts/comp_genomics/orthofinder/orthogroup_stats.py
# Goal:
#   - Compares orthogroup gene counts across groupings using:
#       -  Mann-Whitney with BH correction (2 groups)
#       - or Kruskal-Wallis with BH correction (>2 groups)
# Produces:
#   - orthogroups_{col}.tsv
# Note:
#   - Optionally takes one or more gene-of-interest text files (one gene ID
#   per line) — a manually curated list (config["genes_of_interest"]) and/or
#   the domain-derived list from extract_domain_genes above — and flags
#   which orthogroup(s) each gene falls in, tagged by which list matched.
############################################################################
rule orthogroup_stats:
    input:
        orthogroups      = os.path.join(rp(config["orthofinder_outdir"]), "Results_here", "Orthogroups", "Orthogroups.GeneCount.tsv"),
        orthofinder_done = os.path.join(ORTHO_OUTDIR, "run_orthofinder.done"),
        grouping_key = ancient("grouping_key.tsv"),
        genes_of_interest_files = (
            ([config["genes_of_interest"]] if config.get("genes_of_interest") else [])
            + ([rp(config["domain_significant_genes_outfile"])] if config.get("domain_significant_genes_outfile") else [])
        ),
        **({
            "orthogroups_tsv": os.path.join(rp(config["orthofinder_outdir"]), "Results_here", "Orthogroups", "Orthogroups.tsv"),
        } if (config.get("genes_of_interest") or config.get("domain_significant_genes_outfile")) else {}),
    output:
        done = os.path.join(ORTHO_OUTDIR, "orthogroup_stats.done"),
        tsvs = expand(
            os.path.join(ORTHO_OUTDIR, "stats", "orthogroups_{col}.tsv"),
            col=config.get("grouping_cols",
                ["host_association", "aerotolerance", "phylum", "class", "order"]),
        ),
    params:
        outdir             = os.path.join(ORTHO_OUTDIR, "stats"),
        grouping_cols      = " ".join(config.get("grouping_cols",
            ["host_association", "aerotolerance", "phylum", "class", "order"])),
        min_samples        = config.get("ortho_min_samples", 2),
        min_nonzero        = config.get("ortho_min_nonzero", 2),
        alpha              = config.get("stats_alpha", 0.05),
        orthogroups_tsv    = os.path.join(rp(config["orthofinder_outdir"]), "Results_here", "Orthogroups", "Orthogroups.tsv") if (config.get("genes_of_interest") or config.get("domain_significant_genes_outfile")) else "",
        genes_of_interest_files = " ".join(
            ([config["genes_of_interest"]] if config.get("genes_of_interest") else [])
            + ([rp(config["domain_significant_genes_outfile"])] if config.get("domain_significant_genes_outfile") else [])
        ),
    threads: 1
    log:
        os.path.join(ORTHO_OUTDIR, "logs", "orthogroup_stats.log"),
    conda:
        "envs/scipy.yaml"
    shell:
        r"""
        set -euo pipefail
        exec > >(tee -a {log}) 2>&1
        mkdir -p "$(dirname {log})"

        echo "[$(date)] Running orthogroup statistics"

        mkdir -p "{params.outdir}"

        python3 scripts/comp_genomics/orthofinder/orthogroup_stats.py \
            --orthogroups   "{input.orthogroups}" \
            --grouping-key  "{input.grouping_key}" \
            --outdir        "{params.outdir}" \
            --grouping-cols {params.grouping_cols} \
            --min-samples   {params.min_samples} \
            --min-nonzero   {params.min_nonzero} \
            --alpha         {params.alpha} \
            $([ -n "{params.orthogroups_tsv}" ] && echo "--orthogroups-tsv {params.orthogroups_tsv}") \
            $([ -n "{params.genes_of_interest_files}" ] && echo "--genes-of-interest {params.genes_of_interest_files}")

        echo "[$(date)] orthogroup_stats done"
        echo "done" > {output.done}
        exec 1>&- 2>&-
        wait
        """

############################################################################
# Rule: annotate_sig_orthogroups
# Script:
#   - scripts/comp_genomics/orthofinder/collect_sig_orthogroup_annots.py
#       - scripts/comp_genomics/orthofinder/get_annots_for_ogroup.py
# Goal: to identify all significant orthogroups from orthogroup_stats output TSVs,
# and fetch their per-gene annotations from the pipeline's own InterProScan
# output ({sample}_interprosum.tsv, from clean_interproscan) using
# get_annots_for_ogroup.py, and combine into one TSV.
# Produces:
#   - sig_orthogroup_annotations.tsv
# Note:
#   - Depends on orthogroup_stats and clean_interproscan being complete.
############################################################################
rule annotate_sig_orthogroups:
    input:
        stats_done    = os.path.join(ORTHO_OUTDIR, "orthogroup_stats.done"),
        orthogroups   = os.path.join(rp(config["orthofinder_outdir"]), "Results_here", "Orthogroups", "Orthogroups.tsv"),
        orthofinder_done = os.path.join(ORTHO_OUTDIR, "run_orthofinder.done"),
        interpro_files = lambda wc: expand(
            os.path.join(INTERPRO_OUTDIR, "{sample}", "{sample}_interprosum.tsv"),
            sample=get_filtered_samples(wc),
        ),
    output:
        tsv  = os.path.join(ORTHO_OUTDIR, "sig_orthogroup_annotations.tsv"),
        done = os.path.join(ORTHO_OUTDIR, "annotate_sig_orthogroups.done"),
    params:
        stats_dir       = os.path.join(ORTHO_OUTDIR, "stats"),
        interpro_outdir = INTERPRO_OUTDIR,
        alpha           = config.get("stats_alpha", 0.05),
    threads: 1
    log:
        os.path.join(ORTHO_OUTDIR, "logs", "annotate_sig_orthogroups.log"),
    shell:
        r"""
        set -euo pipefail
        exec > >(tee -a {log}) 2>&1
        mkdir -p "$(dirname {log})"

        echo "[$(date)] Annotating significant orthogroups"

        python3 scripts/comp_genomics/orthofinder/collect_sig_orthogroup_annots.py \
            --stats-dir       "{params.stats_dir}" \
            --orthogroups     "{input.orthogroups}" \
            --interpro-outdir "{params.interpro_outdir}" \
            --output          "{output.tsv}" \
            --alpha           {params.alpha}

        echo "[$(date)] annotate_sig_orthogroups done"
        echo "done" > {output.done}
        exec 1>&- 2>&-
        wait
        """


############################################################################
# Rule: run_diamond
# Script: 
#   - scripts/comp_genomics/alien/diamond.sh
# Goal: Run diamond (self vs self and self vs db) on all proteomes
# Produces: 
#   - {sample}_BS.out (Self output)
#   - {sample}_BN.out (DB output)
############################################################################
rule run_diamond:
    input:
        proteome = os.path.join(FILTERED_FAA_DIR, "{sample}.faa"),
    output:
        bs   = os.path.join(DIAMOND_OUTDIR, "{sample}", "{sample}_BS.out"),
        bn   = os.path.join(DIAMOND_OUTDIR, "{sample}", "{sample}_BN.out"),
        done = os.path.join(DIAMOND_OUTDIR, "{sample}", "run_diamond.done"),
    params:
        diamond_outdir = DIAMOND_OUTDIR,
        diamond_db     = rp(config["diamond_db"]),
        tmp_dir        = rp(config["tmp_dir"]),
    threads: config.get("threads", 4)

    log:
        os.path.join(ALIEN_OUTDIR, "{sample}", "logs", "run_diamond.log"),
    conda:
        "envs/diamond.yaml"
    shell:
        r"""
        set -euo pipefail
        exec > >(tee -a {log}) 2>&1
        mkdir -p "$(dirname {log})"
        echo "[$(date)] Running Diamond for {wildcards.sample}"

        bash scripts/comp_genomics/alien/diamond.sh \
            "{input.proteome}" \
            "{params.diamond_outdir}" \
            "{params.diamond_db}" \
            "{params.tmp_dir}"

        echo "[$(date)] Diamond completed for {wildcards.sample}"
        echo "done" > {output.done}
        exec 1>&- 2>&-
        wait

        """


############################################################################
# Rule: run_alien
# Script: 
#   - scripts/comp_genomics/alien/alien_skip.sh
# Goal: Identify alien index scores for each sample
# Produces: 
#   - {sample}_alien.tsv
############################################################################
rule run_alien:
    input:
        bs = os.path.join(DIAMOND_OUTDIR, "{sample}", "{sample}_BS.out"),
        bn = os.path.join(DIAMOND_OUTDIR, "{sample}", "{sample}_BN.out"),
    output:
        summary = os.path.join(ALIEN_OUTDIR, "{sample}", "{sample}_alien.tsv"),
        done    = os.path.join(ALIEN_OUTDIR, "{sample}", "run_alien.done"),
    params:
        alien_outdir    = ALIEN_OUTDIR,
        alien_script    = rp(config["alien_script"]),
        recipient     = config["recipient"],
        skip_lineages = config["skip_lineages"],
        sqlite_taxa   = rp(config["sqlite_taxa"]),

    threads: config.get("threads", 4)

    log:
        os.path.join(ALIEN_OUTDIR, "{sample}", "logs", "run_alien.log"),
    conda:
        "envs/alien.yaml"
    shell:
        r"""
        set -euo pipefail
        exec > >(tee -a {log}) 2>&1
        mkdir -p "$(dirname {log})"
        echo "[$(date)] Running Alien for {wildcards.sample}"

        bash scripts/comp_genomics/alien/alien_skip.sh \
            "{params.alien_outdir}" \
            "{params.alien_script}" \
            "{input.bs}" \
            "{input.bn}" \
            "{params.recipient}" \
            "{params.skip_lineages}" \
            "{params.sqlite_taxa}"

        echo "[$(date)] Alien completed for {wildcards.sample}"
        echo "done" > {output.done}
        exec 1>&- 2>&-
        wait

        """

############################################################################
# Rule: alien_add_annotations
# Script: 
#   - scripts/comp_genomics/alien/add_ipr_caz.py
# Goal: add dbcan and interproscan annotations to alien output 
# Produces: 
#   - {sample}_alien_ipr.tsv
############################################################################
rule alien_add_annotations:
    input:
        alien = os.path.join(ALIEN_OUTDIR,        "{sample}", "{sample}_alien.tsv"),
        ipr   = os.path.join(INTERPRO_OUTDIR, "{sample}", "{sample}_interprosum.tsv"),
        dbcan = os.path.join(DBCAN_OUTDIR,    "{sample}", "{sample}_dbcan_summary.txt"),
    output:
        tsv  = os.path.join(ALIEN_OUTDIR, "{sample}", "{sample}_alien_ipr.tsv"),
        done = os.path.join(ALIEN_OUTDIR, "{sample}", "alien_add_annotations.done"),
    params:
        outdir = os.path.join(ALIEN_OUTDIR, "{sample}"),
    threads: 1
    log:
        os.path.join(ALIEN_OUTDIR, "{sample}", "logs", "alien_add_annotations.log"),
    shell:
        r"""
        set -euo pipefail
        exec > >(tee -a {log}) 2>&1
        mkdir -p "$(dirname {log})"

        echo "[$(date)] Annotating alien output for {wildcards.sample}"

        python3 scripts/comp_genomics/alien/add_ipr_caz.py \
            --alien  "{input.alien}" \
            --ipr    "{input.ipr}" \
            --dbcan  "{input.dbcan}" \
            --outdir "{params.outdir}"

        echo "[$(date)] Done annotating alien output for for {wildcards.sample}"
        echo "done" > {output.done}
        exec 1>&- 2>&-
        wait
        """


############################################################################
# Rule: summarize_alien  (aggregate — runs once all per-sample alien jobs done)
# Script: 
#   - scripts/comp_genomics/alien/run_alien_summary.py
#       - scripts/comp_genomics/alien/summarize_alien.py
# Goal: 
#   - For each sample, looks up its ingroup term from the column specified by
#   config["alien_ingroup_col"] in grouping_key.tsv (e.g. order, phylum),
#   runs summarize_alien.py with that term as --ingroup
# Produces:
#   - alien_summary.tsv          — one row per sample, summary stats
#   - alien_outside_ingroup.tsv  — all outside-ingroup detail rows, all samples
# Note: 
#   - Both tables include all grouping columns for downstream filtering by clade.
############################################################################
rule summarize_alien:
    input:
        grouping_key = ancient("grouping_key.tsv"),
        alien_files  = lambda wc: expand(
            os.path.join(ALIEN_OUTDIR, "{sample}", "{sample}_alien.tsv"),
            sample=get_filtered_samples(wc),
        ),
    output:
        summary = os.path.join(ALIEN_OUTDIR, "alien_summary.tsv"),
        detail  = os.path.join(ALIEN_OUTDIR, "alien_outside_ingroup.tsv"),
        done    = os.path.join(ALIEN_OUTDIR, "summarize_alien.done"),
    params:
        ingroup_col     = config.get("alien_ingroup_col", "order"),
        threshold       = config.get("alien_threshold", 0.2),
        ingroup_species = config.get("alien_ingroup_species", ""),
    threads: 1
    log:
        os.path.join(ALIEN_OUTDIR, "logs", "summarize_alien.log"),
    shell:
        r"""
        set -euo pipefail
        exec > >(tee -a {log}) 2>&1
        mkdir -p "$(dirname {log})"

        echo "[$(date)] Summarizing AlienIndex results"

        ALIEN_LIST=$(mktemp)
        trap "rm -f $ALIEN_LIST" EXIT
        printf '%s\n' {input.alien_files} > "$ALIEN_LIST"

        python3 scripts/comp_genomics/alien/run_alien_summary.py \
            --alien-files     "$ALIEN_LIST" \
            --grouping-key    "{input.grouping_key}" \
            --ingroup-col     "{params.ingroup_col}" \
            --summary-out     "{output.summary}" \
            --detail-out      "{output.detail}" \
            --threshold       {params.threshold} \
            $([ -n "{params.ingroup_species}" ] && \
                echo "--ingroup-species {params.ingroup_species}" || true)

        echo "[$(date)] summarize_alien done"
        echo "done" > {output.done}
        exec 1>&- 2>&-
        wait
        """


############################################################################
# Rule: filter_alien_tsv  (aggregate — runs after all per-sample alien jobs)
# Script: 
#   - scripts/comp_genomics/alien/filter_alien_tsv.py
# Goal: Counts IPR domain hits per species across alien output TSVs, filtered by
# alien score cutoff
# Produces: 
#   - alien_ipr_counts.tsv
############################################################################
rule filter_alien_tsv:
    input:
        # Needs the IPR-annotated alien output (15 columns: the original 14
        # plus Interpro), not the plain 14-column {sample}_alien.tsv from
        # run_alien — filter_alien_tsv.py reads column 15 for domain IDs.
        alien_files = lambda wc: expand(
            os.path.join(ALIEN_OUTDIR, "{sample}", "{sample}_alien_ipr.tsv"),
            sample=get_filtered_samples(wc),
        ),
        domain_file = rp(config["ipr_entry_list"]),
    output:
        counts    = os.path.join(ALIEN_OUTDIR, "alien_ipr_counts.tsv"),
        done      = os.path.join(ALIEN_OUTDIR, "filter_alien_tsv.done"),
    params:
        cutoff = config.get("alien_threshold", 0.2),
    threads: 1
    log:
        os.path.join(ALIEN_OUTDIR, "logs", "filter_alien_tsv.log"),
    shell:
        r"""
        set -euo pipefail
        exec > >(tee -a {log}) 2>&1
        mkdir -p "$(dirname {log})"

        echo "[$(date)] Counting IPR domains in alien TSVs (cutoff={params.cutoff})"

        ALIEN_LIST=$(mktemp)
        trap "rm -f $ALIEN_LIST" EXIT
        printf '%s\n' {input.alien_files} > "$ALIEN_LIST"


        python3 scripts/comp_genomics/alien/filter_alien_tsv.py \
            --path-file   "$ALIEN_LIST" \
            --domain-file "{input.domain_file}" \
            --cutoff      {params.cutoff} \
            --output      "{output.counts}"

        echo "[$(date)] filter_alien_tsv done"
        echo "done" > {output.done}
        exec 1>&- 2>&-
        wait
        """


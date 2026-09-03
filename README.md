# Convergent Evolution in Gut Eukaryotes: Snakemake Workflow

Snakemake pipeline for comparative genomic analysis of host-associated protists. Starting from genome assemblies (`.fna`) and/or pre-existing proteomes (`.faa`), the workflow performs gene prediction, quality filtering, functional annotation, statistical clustering, and evolutionary analysis to identify convergent gene families and pathways associated with host-association and aerotolerance.

---

## Table of Contents

1. [Requirements](#requirements)
2. [Directory Structure](#directory-structure)
3. [Configuration](#configuration)
4. [How to Run](#how-to-run)
5. [Running via SLURM with `run_snakemake.sh`](#running-via-slurm-with-run_snakemakesh)
6. [Adding New Genomes](#adding-new-genomes)
7. [Workflow Overview](#workflow-overview)
8. [Rules and Scripts](#rules-and-scripts)

---

## Requirements

- Snakemake with SLURM profile (`profiles/slurm/`)
- Conda (environments defined in `envs/`)
- Singularity/Apptainer (for BRAKER3, MetaEuk, BUSCO, KofamScan, dbCAN .sif containers)
- Python 3 with: `biopython`, `ete3`, `scipy`, `numpy`

---

## Directory Structure

```
snakemake/
├── Snakefile.smk                  # Main workflow
├── config/
│   └── config.yaml                # All configurable paths and parameters
├── profiles/
│   └── slurm/                     # Snakemake SLURM profile
├── envs/                          # Conda environment YAML files
│   ├── biopython.yaml
│   ├── ete3.yaml
│   ├── scipy.yaml
│   ├── scipy_numpy.yaml
│   ├── orthofinder.yaml
│   ├── diamond.yaml
│   ├── interproscan.yaml
│   └── alien.yaml
├── scripts/
│   ├── gene_prediction/
│   │   ├── braker.sh
│   │   ├── metaeuk.sh
│   │   └── reduce_metaeuk_faa_name.py
│   ├── data_cleaning/
│   │   ├── busco.sh
│   │   ├── sum_busco.sh
│   │   └── filter_busco.py
│   ├── data_clustering/
│   │   ├── ncbi_taxonomy.py
│   │   ├── add_clustering.py
│   │   ├── filter_representative_genomes.py    # step 0 of adding new genomes — see below
│   │   ├── generate_sample_ids.py
│   │   ├── rename_xlsx_genomes.py
│   │   ├── parse_genome_xlsx.py
│   │   ├── combine_taxonomy.py
│   │   └── rename_and_filter_faa.py
│   ├── create_input_files/
│   │   ├── build_kegg_pathway_hierarchy.py
│   │   ├── build_kegg_pathway_table.py
│   │   ├── get_ciliates.sh
│   │   └── fix_ko_labels.sh
│   ├── gene_annotation/
│   │   ├── kofam/
│   │   │   ├── kofam.sh
│   │   │   ├── kofam_summary.py
│   │   │   ├── pathways_summary_stats.py
│   │   │   ├── pull_pathways.sh
│   │   │   └── KEGG_modules/
│   │   │       ├── module_parser.py
│   │   │       ├── join_tsv.py
│   │   │       └── module_stats.py
│   │   ├── interproscan/
│   │   │   ├── interproscan.sh
│   │   │   ├── clean_interpro.py
│   │   │   ├── interpro_summary.py
│   │   │   ├── interpro_sum_stats.py
│   │   │   ├── all_interpro_sum_stats.py
│   │   │   └── extract_domain_genes.py
│   │   └── dbcan/
│   │       ├── dbcan.sh
│   │       ├── clean_dbcan.py
│   │       ├── combine_dbcan.py
│   │       ├── dbcan_stats.py
│   │       └── all_dbcan_stats.py
│   └── comp_genomics/
│       ├── orthofinder/
│       │   ├── orthofinder.sh
│       │   ├── orthogroup_stats.py
│       │   ├── collect_sig_orthogroup_annots.py
│       │   └── get_annots_for_ogroup.py
│       └── alien/
│           ├── diamond.sh
│           ├── alien_skip.sh
│           ├── add_ipr_caz.py
│           ├── run_alien_summary.py
│           ├── summarize_alien.py
│           └── filter_alien_tsv.py
├── input_files/                   # Static reference files (paths set in config.yaml)
│   ├── aerobic.txt                # One species name per line; used to label samples as aerobic in grouping_key
│   ├── HA.txt                     # One species name per line; used to label samples as host-associated in grouping_key
│   ├── species_list.txt           # One species name per line; input to NCBI taxonomy lookup (skipped if provided_grouping_key is set)
│   ├── sample_key.tsv             # Maps species names to sample IDs (short names matching .faa filenames)
│   ├── grouping_key.tsv           # Pre-built grouping key; set provided_grouping_key in config to use this and skip taxonomy generation
│   ├── pathway_hierarchy.tsv      # KEGG pathway hierarchy (pathway ID → category → name); used in pathway_summary_stats
│   ├── pathway_module_ko.tsv      # Maps KEGG module IDs to their boolean KO definitions; used in assess_modules and module_stats
│   ├── kegg_pathways.txt          # List of KEGG pathway IDs to include in KofamScan summarization (config[kegg_pathway_ids])
│   ├── kegg_orthology_pathways.txt# Two-column tab-separated file mapping KEGG pathway IDs to KO terms (e.g. ko00010<tab>K00001); used in summarize_kofamscan
│   ├── entry.list                 # InterPro entry list (IPR accessions + metadata); used in summarize_interproscan and filter_alien_tsv
│   └── domain_ids.txt             # List of IPR accessions to include in cross-lineage InterPro stat comparisons (config[ipr_compare_list])
├── fastas/
│   ├── {genome_dir}/                  # Input: .fna genome files
│   └── {protein_dir}/                 # Input/Output: .faa proteome files
│      └── low_completeness/          # Quarantine for failed BUSCO samples
├── filtered_faa/            # Length-filtered proteomes (≥ min_seq_len aa); input to all annotation rules
├── outdir/
│   ├── {gp_outdir}/                   # Gene prediction outputs (per sample)
│   ├── {busco_outdir}/                # BUSCO outputs (per sample + summaries)
│   ├── {kofam_outdir}/                # KofamScan outputs
│      └── module_results/            # Per-module functionality tables
│   ├── {interproscan_outdir}/         # InterProScan outputs
│      └── stats/                     # Per-group and merged stat tables
│   ├── {dbcan_outdir}/                # dbCAN outputs
│      └── stats/                     # Per-group and merged stat tables
│   ├── {orthofinder_outdir}/          # OrthoFinder outputs
│      └── stats/                     # Per-grouping-column orthogroup stats
│   ├── {alien_outdir}/                # AlienIndex outputs
├── grouping_key.tsv               # Sample taxonomy + phenotype labels
├── taxonomy_raw.tsv               # Intermediate NCBI taxonomy output
└── logs/                          # SLURM head-job logs
```

Paths in `{}` are resolved from `config/config.yaml` at runtime.

---

## Configuration

All paths and parameters are set in `config/config.yaml`:

| Key | Description |
|---|---|
| `root` / `project` | Base paths; `{root}` and `{project}` are substituted into all other path values |
| `genome_dir` | Directory of input `.fna` genome assemblies |
| `protein_dir` | Directory of input/output `.faa` proteomes |
| `provided_grouping_key` | Path to a pre-built grouping key TSV; leave empty to generate from NCBI taxonomy |
| `species_list` | Input species list for taxonomy lookup (used when `provided_grouping_key` is empty) |
| `braker_sif` / `metaeuk_sif` | Singularity images for gene predictors |
| `busco_sif` / `busco_lineage` | BUSCO container and lineage dataset |
| `busco_min_complete` | Absolute floor (%_complete, default 10) for `filter_proteomes` — a sample is removed if it's below this floor *or* the dynamic Tukey (IQR) fence, whichever is stricter |
| `kofam_sif` / `kofam_profiles` / `kofam_list` | KofamScan container and HMM database |
| `interpro` | InterProScan installation path |
| `dbcan_sif` / `dbcan_db` | dbCAN container and database |
| `diamond_db` | DIAMOND database for AlienIndex |
| `stats_alpha` | Significance threshold (default: 0.05) |
| `grouping_cols` | Columns used for statistical comparisons |
| `dbcan_groups` / `interpro_groups` | Taxonomic groups for per-group stats, each with a `group` name and `col` (the taxonomy column to look it up in) |
| `orthofinder_mode` | `fresh` (default) or `add` (add samples to a prior run); `orthofinder_prev_results` gives the path to the prior `Results_*` dir and supports `{root}`/`{project}` substitution |
| `genes_of_interest` | Optional path to a text file of gene IDs (one per line) to flag in `orthogroup_stats` output |
| `domain_significant_genes_outfile` | Optional path — if set, runs `extract_domain_genes` before `orthogroup_stats` to build a second gene-of-interest list from genes carrying InterPro domains significant in `merge_interpro_stats`' output. Combines with `genes_of_interest` above (both get passed to `orthogroup_stats.py`); leave empty to skip |
| `domain_significant_in` | `any` (default — significant in at least one lineage) or `all` (significant in every lineage), used by `extract_domain_genes` |
| `forced_pathways` | Optional file of pathway IDs to force into module assessment |
| `filtered_faa_dir` | Output directory for length-filtered proteomes produced by `filter_short_seqs` |
| `min_seq_len` | Minimum protein sequence length in aa to retain (default: 50); used by `filter_short_seqs` |
| `sqlite_taxa` | Full absolute path to ete3 taxa SQLite database (`~/.etetoolkit/taxa.sqlite`); but do not use `~` — use the full path |
| `alien_script` | Path to `alien_index_multiskip.py` |
| `recipient` | NCBI taxonomy ID of the recipient/focal group (e.g. 2759 for Eukaryota) |
| `skip_lineages` | Comma-separated NCBI taxids to skip in AlienIndex scoring |
| `alien_threshold` | AlienIndex score cutoff for summarize_alien and filter_alien_tsv (default: 0.2) |
| `alien_ingroup_col` | Grouping key column used to define the ingroup for summarize_alien (e.g. order) |

---

## How to Run

**Before running, edit the following:**

| File/Directory | Description |
|---|---|
| `fastas/fna` | Add new genome files to this directory |
| `fastas/faa` | Add new proteome files to this directory |
| `species_list.txt` | Add full species name |
| `sample_key.tsv` | Add sample ID and full species name |
| `HA.txt` | If species is HA, add full species name |
| `aerobic.txt` | If species is aerobic, add full species name |
| `grouping_key.tsv` | if you want to keep the existing grouping key, add full lineage of species. If not, remove current grouping key. |

**Before running, edit the following in `config/config.yaml`:**

| Field | Description |
|---|---|
| `root` | Your home root |
| `project` | path to the snakemake project directory from root |
| `braker species` | Opt: specify a species database you want braker to train with; example, tetrahymena |
| `forced_pathways` | Opt: KEGG pathway IDs to always include in key_pathways.tsv regardless of significance |
| `interpro_groups and dbcan_groups` | Confirm all your species fall into one of the listed taxonomic groups and add additional groups as needed |
| `interpro_compare_col and dbcan_compare_col` | column with two labels to compared; ex are host_association or aerotolerance |
| `dbcan_cazyme_list` | Opt: CAZyme family IDs to include in merged table, leave empty to include all families found in any lineage TSV |
| `genes_of_interest` | Opt: path to text file with gene IDs of interest for orthogroup analysis |
| `domain_significant_genes_outfile` | Opt: path to write a second, domain-derived gene-of-interest list (see `extract_domain_genes`); leave empty to skip |
| `orthofinder_mode` | Fresh = run orthofinder from scratch; Add = add to a previous run and provide path in `orthofinder_prev_results` |
| `busco_min_complete` | Absolute BUSCO completeness floor (default 10); raise/lower if `filter_proteomes` is removing too much or too little |
| `recipient` | Recipient ID for alien index |
| `skip_lineages` | Lineage IDs to skip in alien index |
| `sqlite_taxa` | Path to an updated SQLite Taxa database |
| `alien_ingroup_species` | Opt: extra species to treat as ingroup regardless of lineage |
| `etc!` | Any other parameters for specific tools you may want to change like significance cutoffs |



**Dry run (always recommended first):**
```bash
snakemake -s Snakefile.smk --profile profiles/slurm -n
```

**Run with config overrides:**
```bash
snakemake -s Snakefile.smk \
    --profile profiles/slurm \
    --config genome_dir="path/to/fna" protein_dir="path/to/faa" \
    --conda-prefix /path/to/conda/envs \
    --rerun-incomplete
```

**Resume an interrupted run:**

Add `--rerun-incomplete` (already included in `run_snakemake.sh`) to pick up from where the workflow left off without re-running completed jobs.

---

## Running via SLURM with `run_snakemake.sh`

`run_snakemake.sh` submits the Snakemake head process as a SLURM job, which then dispatches individual rule jobs to the cluster via the SLURM profile. Before submitting, please make the same necessary changes to the files listed above as you would when running interactively. 

**In addition, please edit the following in `run_snakemake.sh` before submitting:**

| Field | Description |
|---|---|
| `--account` | Your SLURM account |
| `--mail-user` | Notification email |
| `--time` | Wall time for the head job (default 7 days; extend to `21-00:00:00` for large runs) |
| `--conda-prefix` | Path where conda environments will be built/cached |
| `cd /path/to/snakemake` | Working directory containing `Snakefile.smk` |


> **Important:** Since `config/config.yaml` already defines `genome_dir` and `protein_dir` via `{root}/{project}`, you can omit the `--config genome_dir=... protein_dir=...` overrides entirely and let the config file handle them.

**Submit:**
```bash
sbatch run_snakemake.sh
```

Head job logs are written to `logs/snakemake_head_<JOBID>.log`. Per-rule logs land in `<rule_outdir>/<sample>/logs/`.

---

## Adding New Genomes

If you have a batch of new genomes summarized in an xlsx file (columns: genome/sample name, genome assembly type, GTDB-style `Taxonomy` string, and optionally a `Representative genome for each species` flag), the following scripts in `scripts/data_clustering/` prepare them for the pipeline. **These are standalone scripts, not Snakemake rules** — run them manually, in order, before `snakemake` itself runs; nothing here is picked up automatically by a pipeline rerun.

1. **`filter_representative_genomes.py`** — (optional) keep only one genome per species, filtering the xlsx to rows where `Representative genome for each species == "Yes"`. Errors out if any species has zero or more than one such row.
   ```bash
   python3 scripts/data_clustering/filter_representative_genomes.py \
       genome_summary.xlsx representative_genomes.xlsx
   ```

2. **`generate_sample_ids.py`** — reads the xlsx directly and assigns each genome a short `sample_key.tsv`-style code (genus[:2] + species[:3], e.g. `Dasytricha ruminantium` → `Darum`; `Genus sp. TAG` → `GeTAG`). Duplicate-species genomes get `-1`, `-2`, ... suffixes; genuine cross-species code collisions get an `x` (`xx`, ...) appended rather than colliding. Writes two files but modifies nothing:
   ```bash
   python3 scripts/data_clustering/generate_sample_ids.py \
       representative_genomes.xlsx input_files/sample_key.tsv \
       new_sample_key_rows.tsv new_genome_rename_map.tsv
   ```
   Review `new_sample_key_rows.tsv`, then append it: `cat new_sample_key_rows.tsv >> input_files/sample_key.tsv`.

3. **`rename_xlsx_genomes.py`** — applies `new_genome_rename_map.tsv` to the xlsx's genome-name column, writing a **new** xlsx (the original is left untouched). Errors out if any row is missing a mapping rather than writing a partially-renamed file.
   ```bash
   python3 scripts/data_clustering/rename_xlsx_genomes.py \
       representative_genomes.xlsx new_genome_rename_map.tsv representative_genomes_renamed.xlsx
   ```

4. **`parse_genome_xlsx.py`** — parses the renamed xlsx into a TSV with `grouping_key.tsv`'s rank columns plus a `genome_type` column (`monoculture` / `MAG` / `SAG`, mapped from the xlsx's assembly-type column).
   ```bash
   python3 scripts/data_clustering/parse_genome_xlsx.py \
       representative_genomes_renamed.xlsx new_genomes_taxonomy.tsv
   ```

5. **`combine_taxonomy.py`** — merges the new taxonomy TSV with the existing `grouping_key.tsv`, computing `host_association`/`aerotolerance` for the new rows directly from `HA.txt`/`aerobic.txt` (the same rule `add_clustering` uses). If a new genome's species already has a row in the existing `grouping_key.tsv`, the new row is **dropped** and the existing one is kept.
   ```bash
   python3 scripts/data_clustering/combine_taxonomy.py \
       input_files/grouping_key.tsv new_genomes_taxonomy.tsv \
       input_files/HA.txt input_files/aerobic.txt \
       combined_grouping_key.tsv
   ```
   Review the output, then replace `input_files/grouping_key.tsv` with it (make sure `provided_grouping_key` in config still points here).

6. **`rename_and_filter_faa.py`** — copies each new genome's `.faa` file into a directory, renamed to its short sample_id. Only processes what's in the rename map (so pass the representative-filtered map to keep this a filter too), and skips any genome whose species already has a `.faa` in the output directory — so it never creates a same-species duplicate. Copies rather than moving/renaming in place; never deletes anything.
   ```bash
   python3 scripts/data_clustering/rename_and_filter_faa.py \
       new_genome_rename_map.tsv new_sample_key_rows.tsv input_files/sample_key.tsv \
       /path/to/source_faa_dir /path/to/protein_dir
   ```

---

## Workflow Overview

The pipeline runs in five broad stages. `grouping_key.tsv` is a shared dependency consumed by multiple rules across stages.

**Stage 1 — Gene prediction and QC.** `.fna` assemblies are converted to `.faa` proteomes via `gene_prediction` (BRAKER3, falling back to MetaEuk). Both newly predicted and pre-existing `.faa` files feed into `proteome_completeness` (BUSCO), then `filter_proteomes` (checkpoint), which quarantines low-completeness samples — anything below an absolute floor (`busco_min_complete`, catches categorical failures) *or* the dynamic Tukey IQR fence (catches moderate outliers relative to the rest of the batch), whichever removes more. All downstream rules operate only on the samples that pass this filter.

**Stage 1.5 — Sequence length filtering.** After `filter_proteomes`, `filter_short_seqs` removes sequences shorter than `min_seq_len` aa (default 50) from each passing proteome, writing length-filtered copies to `filtered_faa_dir`. All annotation rules (KofamScan, InterProScan, dbCAN, DIAMOND, OrthoFinder) use these filtered proteomes as input. BUSCO runs on the original unfiltered proteomes so completeness scores are not affected.

**Stage 2 — Grouping key.** Built in parallel with annotation: either generated from NCBI taxonomy via `get_taxonomy` → `add_clustering`, or copied from a pre-built file. Feeds as a shared input into `pathway_summary_stats`, `assess_modules`, `summarize_interproscan`, `summarize_dbcan`, `orthogroup_stats`, and `summarize_alien`.

**Stage 3 — Functional annotation (three parallel branches).**
- KofamScan: `run_kofamscan` → `summarize_kofamscan` → `pathway_summary_stats` → `pull_kofam_pathways` (checkpoint) → `assess_modules` → `module_stats`
- InterProScan: `run_interproscan` → `clean_interproscan` → `summarize_interproscan` → `interpro_sum_stats` → `merge_interpro_stats`
- dbCAN: `run_dbcan` → `clean_dbcan` → `summarize_dbcan` → `dbcan_stats` → `merge_dbcan_stats`

**Stage 4 — Comparative genomics (two parallel branches).**
- OrthoFinder: `run_orthofinder` → *(optional)* `extract_domain_genes` → `orthogroup_stats` → `annotate_sig_orthogroups`. `extract_domain_genes` (gated by `domain_significant_genes_outfile`) builds a gene-of-interest list from genes carrying InterPro domains significant in `merge_interpro_stats`' output; `orthogroup_stats` accepts this alongside (or instead of) a manually curated `genes_of_interest` list, tagging matched orthogroups by which list(s) they came from. `annotate_sig_orthogroups` pulls per-gene domain evidence from the pipeline's own InterProScan output (`clean_interproscan`).
- AlienIndex: `run_diamond` → `run_alien` → `alien_add_annotations`. This step joins each sample's alien scores with its `clean_interproscan` and `clean_dbcan` outputs to annotate horizontally transferred genes with domain and CAZyme calls. Results then feed into `summarize_alien` and `filter_alien_tsv`.

**Key cross-branch dependencies.**
- `clean_interproscan` output feeds into both `alien_add_annotations` and `annotate_sig_orthogroups`
- `clean_dbcan` output feeds into `alien_add_annotations`
---

## Rules and Scripts

### Gene Prediction

**`gene_prediction`** — Converts `.fna` assemblies to `.faa` proteomes. Runs BRAKER3 first; falls back to MetaEuk if BRAKER3 fails or produces no output.

| | |
|---|---|
| Input | `{genome_dir}/{sample}.fna`, BRAKER protein DB, MetaEuk protein DB |
| Output | `{protein_dir}/{sample}.faa`, `{gp_outdir}/{sample}/gene_prediction.done` |
| Scripts | `scripts/gene_prediction/braker.sh` — runs BRAKER3 via Singularity |
| | `scripts/gene_prediction/metaeuk.sh` — runs MetaEuk via Singularity |
| | `scripts/gene_prediction/reduce_metaeuk_faa_name.py` — truncates long MetaEuk pipe-delimited headers to the first 4 fields for compatibility with InterProScan and KofamScan |

---

### Quality Control

**`proteome_completeness`** — Runs BUSCO in protein mode on each proteome to assess completeness.

| | |
|---|---|
| Input | `{protein_dir}/{sample}.faa`, BUSCO lineage dataset |
| Output | `{busco_outdir}/{sample}/short_summary.specific.{sample}.txt` |
| Script | `scripts/data_cleaning/busco.sh` |

**`filter_proteomes`** *(checkpoint)* — Aggregates all BUSCO summaries and removes outlier samples: below the absolute `busco_min_complete` floor (catches categorical failures like a near-empty proteome, which a purely relative test can miss once the batch mixes genome types with different expected completeness ranges) *or* below the dynamic Tukey fence (Q1 − 1.5×IQR, catches moderate outliers relative to the rest of the batch) — whichever removes more. Removed `.faa` files are moved to a quarantine folder rather than deleted. All downstream rules use `get_filtered_samples()` to operate only on samples that pass this filter.

| | |
|---|---|
| Input | All per-sample `short_summary.specific.{sample}.txt` files |
| Output | `{busco_outdir}/busco_summary.tsv`, `kept_samples.txt`, `removed_samples.txt` |
| Scripts | `scripts/data_cleaning/sum_busco.sh` — concatenates BUSCO summary files |
| | `scripts/data_cleaning/filter_busco.py` — applies IQR filter and writes kept/removed lists |

**`filter_short_seqs`** — Removes sequences shorter than `min_seq_len` amino acids (default 50) from each proteome that passed BUSCO filtering. Short peptides cause spurious hits in InterProScan, KofamScan, dbCAN, and DIAMOND. Runs after `filter_proteomes`; all annotation rules use the filtered proteomes from `filtered_faa_dir` as their input.

| | |
|---|---|
| Input | `{protein_dir}/{sample}.faa` (post-BUSCO-filter samples only) |
| Output | `{filtered_faa_dir}/{sample}.faa` |
| Script | none — direct `seqkit seq -m {min_seq_len}` call |

---

### Grouping Key

The grouping key links each sample to its taxonomy (phylum through species) and phenotype labels (`host_association`, `aerotolerance`). It is used by all downstream statistical rules.

**Option A — Auto-generate:**

**`get_taxonomy`** — Queries NCBI taxonomy for each species in `species_list`.

| | |
|---|---|
| Input | `config[species_list]` |
| Output | `taxonomy_raw.tsv` |
| Script | `scripts/data_clustering/ncbi_taxonomy.py` |

**`add_clustering`** — Annotates raw taxonomy with host-association and aerotolerance labels, and maps species names to sample IDs.

| | |
|---|---|
| Input | `taxonomy_raw.tsv`, `HA_list`, `aerobic_list`, `sample_key` |
| Output | `grouping_key.tsv` |
| Script | `scripts/data_clustering/add_clustering.py` |

**Option B — Provide your own:** Set `provided_grouping_key` in config to a pre-built TSV path. The workflow copies it to `grouping_key.tsv` and skips taxonomy lookup entirely.

---

### KofamScan (KEGG Ortholog Annotation)

**`run_kofamscan`** — Annotates each proteome with KEGG Ortholog (KO) terms using KofamScan via Singularity.

| | |
|---|---|
| Input | `{filtered_faa_dir}/{sample}.faa` |
| Output | `{kofam_outdir}/{sample}/{sample}_ko_annotations.tsv` |
| Script | `scripts/gene_annotation/kofam/kofam.sh` |

**`summarize_kofamscan`** — Merges per-sample KO annotations into a single pathway-level summary table.

| | |
|---|---|
| Input | All `{sample}_ko_annotations.tsv`, KEGG pathway IDs, KO-pathway correspondence |
| Output | `{kofam_outdir}/kofam_summary.tsv` |
| Script | `scripts/gene_annotation/kofam/kofam_summary.py` |

**`pathway_summary_stats`** — Joins the KofamScan summary with pathway hierarchy and grouping key; computes min/max/mean pathway completion rates for HA vs. FL and aerobic vs. anaerobic, with Mann-Whitney U + Benjamini-Hochberg correction.

| | |
|---|---|
| Input | `kofam_summary.tsv`, pathway hierarchy, `grouping_key.tsv` |
| Output | `{kofam_outdir}/pathway_sumstats.tsv` |
| Script | `scripts/gene_annotation/kofam/pathways_summary_stats.py` |

**`pull_kofam_pathways`** *(checkpoint)* — Filters `pathway_sumstats.tsv` for three metabolite groups (cofactors/vitamins, amino acid metabolism, lipid metabolism), keeping the top 10 significant pathways per group. Optional forced pathways can be added via config.

| | |
|---|---|
| Input | `pathway_sumstats.tsv`, optional `forced_pathways` file |
| Output | `{kofam_outdir}/key_pathways.tsv` |
| Script | `scripts/gene_annotation/kofam/pull_pathways.sh` |

**`assess_modules`** — For each KEGG module in `key_pathways.tsv`, evaluates module functionality (complete/incomplete) per sample using the module's boolean KO definition. Runs one job per module over all filtered samples.

| | |
|---|---|
| Input | `pathway_module_ko.tsv`, `key_pathways.tsv`, all `{sample}_ko_annotations.tsv`, `grouping_key.tsv` |
| Output | `{kofam_outdir}/module_results/{module}_functionality.tsv` |
| Scripts | `scripts/gene_annotation/kofam/KEGG_modules/module_parser.py` — evaluates boolean module completeness per sample |
| | `scripts/gene_annotation/kofam/KEGG_modules/join_tsv.py` — left-joins functionality results with grouping key |

**`module_stats`** — Tests whether module functionality rates differ across grouping columns (Fisher's exact test for 2 groups, chi-squared for >2 groups) with BH correction.

| | |
|---|---|
| Input | All `{module}_functionality.tsv`, `pathway_module_ko.tsv` |
| Output | `{kofam_outdir}/module_stats.tsv` |
| Script | `scripts/gene_annotation/kofam/KEGG_modules/module_stats.py` |

---

### InterProScan (Domain / Family Annotation)

**`run_interproscan`** — Annotates each proteome with protein domains and families using InterProScan.

| | |
|---|---|
| Input | `{filtered_faa_dir}/{sample}.faa` |
| Output | `{interproscan_outdir}/{sample}/{sample}_interproscan.tsv` |
| Script | `scripts/gene_annotation/interproscan/interproscan.sh` |

**`clean_interproscan`** — Reduces redundancy in InterProScan output, collapsing to a per-gene summary.

| | |
|---|---|
| Input | `{filtered_faa_dir}/{sample}.faa`, `{sample}_interproscan.tsv` |
| Output | `{interproscan_outdir}/{sample}/{sample}_interprosum.tsv` |
| Script | `scripts/gene_annotation/interproscan/clean_interpro.py` |

**`summarize_interproscan`** — Combines per-sample InterPro summaries into one multi-species table annotated with the grouping key.

| | |
|---|---|
| Input | All `{sample}_interprosum.tsv`, all `{filtered_faa_dir}/{sample}.faa` files, `grouping_key.tsv`, IPR entry list |
| Output | `{interproscan_outdir}/all_interpro_summary.tsv` |
| Script | `scripts/gene_annotation/interproscan/interpro_summary.py` |

**`interpro_sum_stats`** — Computes per-group Mann-Whitney U and t-test statistics for domain enrichment across host-association or aerotolerance. Runs one job per taxonomic group defined in `config[interpro_groups]`. If a group has only one value for the compare column, empty TSVs are written for both `_ttest.tsv` and `_mannwhitney.tsv` instead of erroring.

| | |
|---|---|
| Input | `all_interpro_summary.tsv` |
| Output | `stats/{group}_mannwhitney.tsv`, `stats/{group}_ttest.tsv` |
| Script | `scripts/gene_annotation/interproscan/interpro_sum_stats.py` |

**`merge_interpro_stats`** — Merges all per-group stat TSVs into one wide-format cross-lineage table.

| | |
|---|---|
| Input | All `{group}_mannwhitney.tsv` (or ttest), IPR compare list |
| Output | `stats/all_lineages_mannwhitney.tsv` |
| Script | `scripts/gene_annotation/interproscan/all_interpro_sum_stats.py` |

---

### dbCAN (CAZyme Annotation)

**`run_dbcan`** — Annotates each proteome with CAZyme families, CGCs, and substrate predictions using dbCAN via Singularity.

| | |
|---|---|
| Input | `{filtered_faa_dir}/{sample}.faa` |
| Output | `{dbcan_outdir}/{sample}/{sample}_overview.txt` |
| Script | `scripts/gene_annotation/dbcan/dbcan.sh` |

**`clean_dbcan`** — Collapses dbCAN's three-tool predictions into a single consensus summary; removes any rows flagged as ERROR.

| | |
|---|---|
| Input | `{sample}_overview.txt` |
| Output | `{dbcan_outdir}/{sample}/{sample}_dbcan_summary.txt` |
| Script | `scripts/gene_annotation/dbcan/clean_dbcan.py` |

**`summarize_dbcan`** — Combines per-sample dbCAN summaries into one multi-species table annotated with the grouping key.

| | |
|---|---|
| Input | All `{sample}_dbcan_summary.txt`, all `{filtered_faa_dir}/{sample}.faa` files, `grouping_key.tsv` |
| Output | `{dbcan_outdir}/dbcan_combined.tsv` |
| Script | `scripts/gene_annotation/dbcan/combine_dbcan.py` |

**`dbcan_stats`** — Computes per-group Mann-Whitney U statistics for CAZyme family enrichment. Runs one job per taxonomic group in `config[dbcan_groups]`. If a group has only one value for the compare column (e.g. only HA Apicomplexa), an empty TSV is written instead of erroring.

| | |
|---|---|
| Input | `dbcan_combined.tsv` |
| Output | `stats/{group}_dbcan_stats.tsv` |
| Script | `scripts/gene_annotation/dbcan/dbcan_stats.py` |

**`merge_dbcan_stats`** — Merges all per-group dbCAN stat TSVs into one wide-format cross-lineage table.

| | |
|---|---|
| Input | All `{group}_dbcan_stats.tsv`, optional CAZyme ID filter list |
| Output | `stats/all_lineages_dbcan_stats.tsv` |
| Script | `scripts/gene_annotation/dbcan/all_dbcan_stats.py` |

---

### OrthoFinder (Orthogroup Analysis)

**`run_orthofinder`** — Clusters all proteomes into orthogroups. Supports two modes: `fresh` (start from scratch) or `add` (extend a previous run with new samples using OrthoFinder's `-b` flag).

| | |
|---|---|
| Input | All `{filtered_faa_dir}/{sample}.faa` files (filtered samples only) |
| Output | `Orthogroups/Orthogroups.tsv`, `Orthogroups/Orthogroups.GeneCount.tsv` |
| Script | `scripts/comp_genomics/orthofinder/orthofinder.sh` |

**`extract_domain_genes`** *(optional — runs before `orthogroup_stats`)* — Finds genes carrying InterPro domains that were significant in `merge_interpro_stats`' output (a domain counts if `{lineage}_significant_bh == "yes"` in `domain_significant_in` lineage(s), `any` by default), by scanning per-sample `{sample}_interprosum.tsv` files for those domain IDs. Writes a gene list in the same format `orthogroup_stats.py`'s `--genes-of-interest` expects. Enabled by setting `domain_significant_genes_outfile` in config; skipped (default) if left empty.

| | |
|---|---|
| Input | `all_lineages_{mannwhitney,ttest}.tsv`, all `{sample}_interprosum.tsv` |
| Output | Gene list at `domain_significant_genes_outfile` |
| Script | `scripts/gene_annotation/interproscan/extract_domain_genes.py` |

**`orthogroup_stats`** — Tests whether orthogroup gene counts differ across each grouping column using Mann-Whitney U (2 groups) or Kruskal-Wallis (>2 groups) with BH correction. Optionally filters to orthogroups containing genes of interest, accepting multiple gene-list files at once (a manually curated `genes_of_interest` list and/or `extract_domain_genes`' output) — matched orthogroups are tagged with which list(s) they came from via a `gene_of_interest_source` column. If a grouping column has fewer than 2 testable groups (e.g. only HA Apicomplexa with no FL counterpart), an empty TSV with just the header is written so Snakemake does not fail.

| | |
|---|---|
| Input | `Orthogroups.GeneCount.tsv`, `grouping_key.tsv`, optional `genes_of_interest` and/or `extract_domain_genes` gene lists |
| Output | `stats/orthogroups_{col}.tsv` for each grouping column (columns include `gene_of_interest`, `gene_of_interest_source`) |
| Script | `scripts/comp_genomics/orthofinder/orthogroup_stats.py` |

**`annotate_sig_orthogroups`** — Identifies significant orthogroups from stat TSVs and fetches per-gene annotations from the pipeline's own InterProScan output (`{sample}_interprosum.tsv`, from `clean_interproscan`).

| | |
|---|---|
| Input | Orthogroup stats TSVs, `Orthogroups.tsv`, all `{sample}_interprosum.tsv` |
| Output | `{orthofinder_outdir}/sig_orthogroup_annotations.tsv` |
| Scripts | `scripts/comp_genomics/orthofinder/collect_sig_orthogroup_annots.py` |
| | `scripts/comp_genomics/orthofinder/get_annots_for_ogroup.py` |

---

### AlienIndex (Horizontal Gene Transfer Detection)

**`run_diamond`** — Runs DIAMOND BLASTP for each sample in two modes: self-vs-self (`_BS.out`) and self-vs-reference-database (`_BN.out`). These outputs are the inputs to the AlienIndex score calculation.

| | |
|---|---|
| Input | `{filtered_faa_dir}/{sample}.faa` |
| Output | `{diamond_outdir}/{sample}/{sample}_BS.out`, `{diamond_outdir}/{sample}/{sample}_BN.out` |
| Script | `scripts/comp_genomics/alien/diamond.sh` |

**`run_alien`** — Calculates AlienIndex scores for each sample using the DIAMOND outputs. Scores reflect the relative similarity of each gene to recipient vs. donor lineages.

| | |
|---|---|
| Input | `{diamond_outdir}/{sample}/{sample}_BS.out`, `{diamond_outdir}/{sample}/{sample}_BN.out` |
| Output | `{alien_outdir}/{sample}/{sample}_alien.tsv` |
| Script | `scripts/comp_genomics/alien/alien_skip.sh` |

**`alien_add_annotations`** — Annotates AlienIndex output with InterProScan domain hits and dbCAN CAZyme family calls.

| | |
|---|---|
| Input | `{sample}_alien.tsv`, `{sample}_interprosum.tsv`, `{sample}_dbcan_summary.txt` |
| Output | `{alien_outdir}/{sample}/{sample}_alien_ipr.tsv` |
| Script | `scripts/comp_genomics/alien/add_ipr_caz.py` |

**`summarize_alien`** — Aggregates per-sample AlienIndex results across all filtered samples, using the ingroup lineage from `grouping_key.tsv` to distinguish outside-ingroup genes.

| | |
|---|---|
| Input | All `{sample}_alien.tsv`, `grouping_key.tsv` |
| Output | `{alien_outdir}/alien_summary.tsv`, `alien_outside_ingroup.tsv` |
| Scripts | `scripts/comp_genomics/alien/run_alien_summary.py` |
| | `scripts/comp_genomics/alien/summarize_alien.py` |

**`filter_alien_tsv`** — Counts InterPro domain hits per species across AlienIndex outputs, filtered by an alien score cutoff.

| | |
|---|---|
| Input | All `{sample}_alien.tsv`, IPR entry list |
| Output | `{alien_outdir}/alien_ipr_counts.tsv` |
| Script | `scripts/comp_genomics/alien/filter_alien_tsv.py` |

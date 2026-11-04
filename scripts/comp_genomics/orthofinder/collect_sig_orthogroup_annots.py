#!/usr/bin/env python3
"""
Collects significant orthogroups from orthogroup_stats output TSVs,
runs get_annots_for_ogroup.py on each, and combines all annotations
into one output TSV.

Usage:
    python3 collect_sig_orthogroup_annots.py \\
        --stats-dir       orthogroup_stats/ \\
        --orthogroups     Orthogroups.tsv \\
        --interpro-outdir path/to/gene_annotation/interproscan \\
        --output          sig_orthogroup_annotations.tsv \\
        --alpha           0.05
"""

import argparse
import csv
import os
import subprocess
import sys
import tempfile

# orthogroups_{col}.tsv's gene_of_interest column can exceed the csv module's
# default 131072-byte field limit when an orthogroup accumulates a lot of
# matched gene IDs (e.g. a common domain in --genes-of-interest hitting many
# genes in one orthogroup) — raise it as high as this platform allows.
_maxint = sys.maxsize
while True:
    try:
        csv.field_size_limit(_maxint)
        break
    except OverflowError:
        _maxint = int(_maxint / 10)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--stats-dir",       required=True,
                        help="Directory containing orthogroups_{col}.tsv files "
                             "from orthogroup_stats.py")
    parser.add_argument("--orthogroups",     required=True,
                        help="Orthogroups.tsv from OrthoFinder")
    parser.add_argument("--interpro-outdir", required=True,
                        help="InterProScan outdir; expects "
                             "{interpro_outdir}/{sample}/{sample}_interprosum.tsv")
    parser.add_argument("--output",          required=True,
                        help="Output TSV path")
    parser.add_argument("--alpha",           type=float, default=0.05,
                        help="BH-corrected p-value cutoff for significance (default: 0.05)")
    return parser.parse_args()


def collect_significant_ogs(stats_dir, alpha):
    """
    Read all orthogroups_{col}.tsv files in stats_dir and return a dict:
        {orthogroup_id: {grouping_col, p_value_bh, lowest_mean_group,
                         groups_compared, gene_of_interest, gene_of_interest_source}}
    If the same OG is significant in multiple columns, all are kept.
    Returns a list of (og, col, row) tuples.
    """
    sig_ogs = []   # (og_id, col_name, row_dict)

    for fname in sorted(os.listdir(stats_dir)):
        # Only process the main per-column files, not the _genes_of_interest ones
        if not (fname.startswith("orthogroups_") and
                fname.endswith(".tsv") and
                "genes_of_interest" not in fname):
            continue

        col = fname.replace("orthogroups_", "").replace(".tsv", "")
        fpath = os.path.join(stats_dir, fname)

        with open(fpath) as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                if row.get("significant_bh", "no") == "yes":
                    try:
                        p = float(row.get("p_value_bh", 1.0))
                    except ValueError:
                        continue
                    if p <= alpha:
                        sig_ogs.append((row["orthogroup"], col, row))

    return sig_ogs


def get_annots_for_og(og_id, orthogroups_tsv, interpro_outdir):
    """
    Call get_annots_for_ogroup.py for one orthogroup.
    Returns list of output lines (excluding the header).
    """
    script = os.path.join(os.path.dirname(__file__), "get_annots_for_ogroup.py")
    result = subprocess.run(
        [sys.executable, script,
         "--ogroup",          og_id,
         "--orthogroups",     orthogroups_tsv,
         "--interpro-outdir", interpro_outdir],
        capture_output=True, text=True
    )

    if result.returncode != 0:
        print(f"  WARNING: get_annots_for_ogroup.py failed for {og_id}: "
              f"{result.stderr.strip()}", file=sys.stderr)
        return [], []

    lines = result.stdout.strip().split("\n")
    if not lines or not lines[0]:
        return [], []

    header = lines[0].split("\t")
    data   = [l.split("\t") for l in lines[1:] if l]
    return header, data


def main():
    args = parse_args()

    print("Collecting significant orthogroups...", file=sys.stderr)
    sig_ogs = collect_significant_ogs(args.stats_dir, args.alpha)

    if not sig_ogs:
        print(f"No significant orthogroups found at alpha={args.alpha}",
              file=sys.stderr)
        # Write empty output with header so Snakemake output file still exists
        with open(args.output, "w") as fh:
            fh.write("orthogroup\tsig_grouping_col\tp_value_bh\tlowest_mean_group\t"
                     "groups_compared\tgene_of_interest\tgene_of_interest_source\t"
                     "Sample\tGene\tCDD-IPR\tCoils\tFunFam\tGene3D\tMobiDBLite\t"
                     "NCBIfam\tPANTHER\tPfam\tPIRSF\tPRINTS\tSFLD\tSMART\t"
                     "SUPERFAMILY\tInterpro\tGO\n")
        return

    # Deduplicate OGs but keep all significant col entries for the metadata
    unique_ogs   = sorted({og for og, col, row in sig_ogs})
    og_meta      = {}   # og -> list of (col, row) — may be multiple columns
    for og, col, row in sig_ogs:
        og_meta.setdefault(og, []).append((col, row))

    print(f"  {len(unique_ogs)} unique significant orthogroups across "
          f"{len(sig_ogs)} grouping-column combinations", file=sys.stderr)

    # Get annotation header from first successful OG call
    annot_header = None
    rows_written = 0

    with open(args.output, "w") as fh:
        for og in unique_ogs:
            print(f"  Fetching annotations for {og}...", file=sys.stderr)
            header, annot_rows = get_annots_for_og(
                og, args.orthogroups, args.interpro_outdir)

            if not header:
                print(f"  WARNING: no annotations found for {og}", file=sys.stderr)
                continue

            # Write header once on first successful result
            if annot_header is None:
                # Prepend metadata columns; skip 'Ogroup' from annot header (col 0)
                meta_cols   = ["orthogroup", "sig_grouping_col", "p_value_bh",
                               "lowest_mean_group", "groups_compared",
                               "gene_of_interest", "gene_of_interest_source"]
                annot_cols  = header[1:]   # drop 'Ogroup' — we use orthogroup from meta
                annot_header = meta_cols + annot_cols
                fh.write("\t".join(annot_header) + "\n")

            # Write one row per annotation × per significant grouping column
            for col, row in og_meta[og]:
                meta_vals = [
                    og,
                    col,
                    row.get("p_value_bh",               "NA"),
                    row.get("lowest_mean_group",         "NA"),
                    row.get("groups_compared",           "NA"),
                    row.get("gene_of_interest",          "NA"),
                    row.get("gene_of_interest_source",   "NA"),
                ]
                for annot_row in annot_rows:
                    # annot_row[0] is Ogroup (redundant), skip it
                    fh.write("\t".join(meta_vals + annot_row[1:]) + "\n")
                    rows_written += 1

    print(f"Done. {rows_written} annotation rows written to {args.output}",
          file=sys.stderr)


if __name__ == "__main__":
    main()

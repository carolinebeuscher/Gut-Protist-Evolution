#!/usr/bin/env python3
"""
Combine per-species dbCAN CAZyme summary files into one multi-species table.

Usage:
    python3 combine_dbcan.py \\
        --txt-list      list_of_dbcan_txts.txt \\
        --faa-list      list_of_faa_files.txt \\
        --grouping-key  grouping_key.tsv \\
        --output        cazyme_combined.tsv
"""

import argparse
import csv
import os
import sys
from collections import defaultdict

# Can edit to be any column in grouping_key.tsv: 
GROUPING_COLS = ["host_association", "aerotolerance", "phylum", "class", "order", "clade_1", "clade_2"]


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--txt-list",     required=True,
                        help="Text file with paths to per-species dbCAN summary TXTs")
    parser.add_argument("--faa-list",     required=True,
                        help="Text file with paths to .faa files (one per line)")
    parser.add_argument("--grouping-key", required=True,
                        help="grouping_key.tsv produced by add_clustering rule")
    parser.add_argument("--output",       required=True,
                        help="Output TSV path")
    return parser.parse_args()


def load_grouping_key(path):
    """
    Load grouping key — indexed by sample_id (filename stem) so it matches
    the species names derived from TXT/FAA filenames. Also indexes by species
    name as a fallback for cases where sample_id is absent or doesn't match.
    """
    grouping = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        available = reader.fieldnames or []
        missing   = [c for c in GROUPING_COLS if c not in available]
        if missing:
            print(f"WARNING: grouping columns not in grouping key "
                  f"(will be NA): {missing}")
        for row in reader:
            sid = row.get("sample_id", "").strip()
            sp  = row.get("species",   "").strip()
            entry = {c: row.get(c, "NA").strip()
                     for c in GROUPING_COLS if c in available}
            if sid:
                grouping[sid] = entry
            if sp and sp not in grouping:
                grouping[sp] = entry   # fallback; don't overwrite sample_id entry
    return grouping



def count_faa_sequences(faa_path):
    """Count sequences in a FASTA file."""
    count = 0
    with open(faa_path) as f:
        for line in f:
            if line.startswith(">"):
                count += 1
    return count


def compute_scale_factor(total_genes_map):
    """Normalise genome sizes to [0, 1] across all species."""
    sizes    = list(total_genes_map.values())
    lo, hi   = min(sizes), max(sizes)
    if hi == lo:
        return {sp: 1.0 for sp in total_genes_map}
    return {sp: round((n - lo) / (hi - lo), 4) for sp, n in total_genes_map.items()}


def process_cazyme_txt(txt_path):
    """
    Parse a per-species dbCAN summary file.
    Returns {cazyme_family: {"genes": set(), "occurrence_sum": int}}.
    """
    domain_data = defaultdict(lambda: {"genes": set(), "occurrence_sum": 0})

    with open(txt_path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) < 2:
                continue
            gene_id   = cols[0].strip()
            cazyme_id = cols[1].strip()

            if gene_id.lower() in ("gene_id", "geneid", "gene", "#gene_id"):
                continue
            if not cazyme_id or cazyme_id == "-":
                continue

            domain_data[cazyme_id]["genes"].add(gene_id)
            domain_data[cazyme_id]["occurrence_sum"] += 1

    return domain_data


def main():
    args = parse_args()

    # Load file lists
    with open(args.txt_list) as f:
        txt_paths = [l.strip() for l in f if l.strip()]
    with open(args.faa_list) as f:
        faa_paths = [l.strip() for l in f if l.strip()]

    if not txt_paths:
        sys.exit("ERROR: txt_list is empty.")
    if not faa_paths:
        sys.exit("ERROR: faa_list is empty.")

    # Load grouping key
    grouping = load_grouping_key(args.grouping_key)

    # Build species -> total_genes from faa files
    faa_map = {os.path.basename(p).replace(".faa", ""): count_faa_sequences(p)
               for p in faa_paths}

    scale_factors = compute_scale_factor(faa_map)

    out_cols = (["species"] + GROUPING_COLS +
                ["cazyme_family_id", "gene_count", "total_genes",
                 "pct_genes", "occurrence_sum", "pct_occurrences", "scale_factor"])

    rows      = []
    unmatched = []

    for txt_path in txt_paths:
        txt_basename = os.path.basename(txt_path)

        # Match txt file to a species in faa_map by prefix
        species = next((sp for sp in faa_map if txt_basename.startswith(sp)), None)
        if species is None:
            unmatched.append(txt_basename)
            print(f"  WARNING: no .faa match for '{txt_basename}' — skipping")
            continue

        total_genes = faa_map[species]
        print(f"  Processing {species} ({total_genes} genes)...")

        sp_grouping = grouping.get(species, {})
        if not sp_grouping:
            print(f"  WARNING: '{species}' not found in grouping_key.tsv — "
                  f"grouping columns will be NA")

        try:
            domain_data = process_cazyme_txt(txt_path)
        except Exception as e:
            print(f"  ERROR processing {txt_path}: {e} — skipping")
            continue

        for cazyme_id, data in domain_data.items():
            gene_count     = len(data["genes"])
            occurrence_sum = data["occurrence_sum"]
            pct_genes      = round(gene_count     / total_genes * 100, 4) if total_genes else 0
            pct_occ        = round(occurrence_sum / total_genes * 100, 4) if total_genes else 0

            row = {"species":          species,
                   "cazyme_family_id": cazyme_id,
                   "gene_count":       gene_count,
                   "total_genes":      total_genes,
                   "pct_genes":        pct_genes,
                   "occurrence_sum":   occurrence_sum,
                   "pct_occurrences":  pct_occ,
                   "scale_factor":     scale_factors[species]}

            for col in GROUPING_COLS:
                row[col] = sp_grouping.get(col, "NA")

            rows.append(row)

    with open(args.output, "w") as fh:
        fh.write("\t".join(out_cols) + "\n")
        for row in rows:
            fh.write("\t".join(str(row.get(c, "NA")) for c in out_cols) + "\n")

    print(f"\nDone. {len(rows)} rows written to {args.output}")
    if unmatched:
        print(f"Unmatched files (no .faa counterpart): {unmatched}")


if __name__ == "__main__":
    main()

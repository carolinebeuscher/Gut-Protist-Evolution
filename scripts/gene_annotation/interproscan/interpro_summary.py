#!/usr/bin/env python3
"""
Combines per-species InterProScan summary TSVs into one multi-species table.

Usage:
    python3 interpro_summary.py \\
        --tsv-list      list_of_interproscan_tsvs.txt \\
        --faa-list      list_of_faa_files.txt \\
        --grouping-key  grouping_key.tsv \\
        --ipr           entry.list \\
        --output        all_interpro_sum.tsv
"""

import argparse
import re
import os
import csv
from collections import defaultdict


# Columns to pull from grouping_key.tsv and carry into the output.
# Can edit to include any column present in grouping_key.tsv:
GROUPING_COLS = ["host_association", "aerotolerance", "phylum", "class", "order", "clade_1", "clade_2"]


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--tsv-list",     required=True,
                        help="Text file with paths to per-species InterProScan summary TSVs")
    parser.add_argument("--faa-list",     required=True,
                        help="Text file with paths to .faa files (one per line)")
    parser.add_argument("--grouping-key", required=True,
                        help="grouping_key.tsv produced by add_clustering rule")
    parser.add_argument("--ipr",          required=True,
                        help="IPR entry list: col 1 = IPR ID, col 3 = domain name")
    parser.add_argument("--output",       required=True,
                        help="Output TSV path")
    return parser.parse_args()



def count_faa_sequences(faa_path):
    """Count number of sequences in a FASTA file."""
    count = 0
    with open(faa_path) as f:
        for line in f:
            if line.startswith(">"):
                count += 1
    return count


def parse_interpro_field(interpro_str):
    """
    Parse the InterPro column (col 15, index 14).
    Returns a list of (ipr_id, domain_name, count) tuples.
    """
    interpro_str = interpro_str.strip()
    if not interpro_str or interpro_str == "-":
        return []

    results = []
    entries = re.split(r',(?=IPR\d+)', interpro_str)

    for entry in entries:
        entry = entry.strip()
        if not entry:
            continue
        count = 1
        count_match = re.search(r'\s*\((\d+)\)\s*$', entry)
        if count_match:
            count = int(count_match.group(1))
            entry = entry[:count_match.start()].strip()
        parts = entry.split("-", 1)
        if len(parts) != 2:
            continue
        ipr_id, domain_name = parts[0].strip(), parts[1].strip()
        if not re.match(r'^IPR\d+$', ipr_id):
            continue
        results.append((ipr_id, domain_name, count))

    return results


def process_tsv(tsv_path, ipr_targets):
    """
    Parse one per-species InterProScan summary TSV.
    Returns {(ipr_id, domain_name): {genes, occurrence_sum, go_terms, superfamily}}.
    """
    domain_data = {
        (ipr_id, name): {"genes": set(), "occurrence_sum": 0}
        for ipr_id, name in ipr_targets.items()
    }
    domain_go = {k: set() for k in domain_data}
    domain_sf = {k: set() for k in domain_data}

    with open(tsv_path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.startswith("Gene\t"):
                continue
            cols = line.split("\t")
            if len(cols) < 15:
                continue

            gene_id      = cols[0].strip()
            interpro_str = cols[14].strip()
            go_str       = cols[15].strip() if len(cols) > 15 else "-"
            sf_str       = cols[13].strip() if len(cols) > 13 else "-"

            for ipr_id, domain_name, count in parse_interpro_field(interpro_str):
                if ipr_id not in ipr_targets:
                    continue
                canonical = ipr_targets[ipr_id]
                key = (ipr_id, canonical)
                domain_data[key]["genes"].add(gene_id)
                domain_data[key]["occurrence_sum"] += count
                if go_str and go_str != "-":
                    for t in go_str.split(","):
                        t = t.strip()
                        if t and t != "-":
                            domain_go[key].add(t)
                if sf_str and sf_str != "-":
                    for s in sf_str.split(","):
                        s = s.strip()
                        if s and s != "-":
                            domain_sf[key].add(s)

    for key in domain_data:
        domain_data[key]["go_terms"]   = ",".join(sorted(domain_go[key])) or "NA"
        domain_data[key]["superfamily"] = ",".join(sorted(domain_sf[key])) or "NA"

    return domain_data


def compute_scale_factor(total_genes_map):
    """Normalise genome sizes to [0, 1] across all species."""
    sizes = list(total_genes_map.values())
    lo, hi = min(sizes), max(sizes)
    if hi == lo:
        return {sp: 1.0 for sp in total_genes_map}
    return {sp: round((n - lo) / (hi - lo), 4) for sp, n in total_genes_map.items()}



def main():
    args = parse_args()

    # Load file lists
    with open(args.tsv_list) as f:
        tsv_paths = [l.strip() for l in f if l.strip()]
    with open(args.faa_list) as f:
        faa_paths = [l.strip() for l in f if l.strip()]

    # Load grouping key
    # index by sample_id (filename stem) so it matches the species names 
    # derived from TSV/FAA filenames. 
    # Also index by species name as a fallback for cases where sample_id 
    # is absent or doesn't match.
    
    grouping = {}   # sample_id (or species) -> {col: value}
    with open(args.grouping_key) as f:
        reader = csv.DictReader(f, delimiter="\t")
        available_cols = reader.fieldnames or []
        for row in reader:
            sid = row.get("sample_id", "").strip()
            sp  = row.get("species",   "").strip()
            entry = {c: row.get(c, "NA").strip()
                     for c in GROUPING_COLS if c in available_cols}
            if sid:
                grouping[sid] = entry
            if sp and sp not in grouping:
                grouping[sp] = entry   # fallback; don't overwrite sample_id entry

    missing_grouping_cols = [c for c in GROUPING_COLS if c not in available_cols]
    if missing_grouping_cols:
        print(f"WARNING: grouping columns not found in grouping key "
              f"and will be NA: {missing_grouping_cols}")


    # Load IPR targets: col 1 = IPR ID, col 3 = domain name
    ipr_targets = {}
    with open(args.ipr) as f:
        for line in f:
            cols = line.strip().split("\t")
            if len(cols) >= 3 and cols[0].strip():
                ipr_targets[cols[0].strip()] = cols[2].strip()

    # Build species -> total_genes from faa files
    faa_map = {os.path.basename(p).replace(".faa", ""): count_faa_sequences(p)
               for p in faa_paths}

    scale_factors = compute_scale_factor(faa_map)

    # Output columns
    out_cols = (["species"] + GROUPING_COLS +
                ["domain_id", "domain_name", "gene_count", "occurrence_sum",
                 "total_genes", "pct_genes", "pct_occurrences",
                 "scale_factor", "go_terms", "superfamily"])

    rows      = []
    unmatched = []

    for tsv_path in tsv_paths:
        tsv_basename = os.path.basename(tsv_path)

        # Match TSV to a species in faa_map by prefix
        species = next((sp for sp in faa_map if tsv_basename.startswith(sp)), None)
        if species is None:
            unmatched.append(tsv_basename)
            print(f"  WARNING: no .faa match for '{tsv_basename}' — skipping")
            continue

        total_genes = faa_map[species]
        print(f"  Processing {species} ({total_genes} genes)...")

        sp_grouping = grouping.get(species, {})
        if not sp_grouping:
            print(f"  WARNING: '{species}' not found in grouping_key.tsv — "
                  f"grouping columns will be NA")

        try:
            domain_data = process_tsv(tsv_path, ipr_targets)
        except Exception as e:
            print(f"  ERROR processing {tsv_path}: {e} — skipping")
            continue

        for (ipr_id, domain_name), data in domain_data.items():
            gene_count     = len(data["genes"])
            occurrence_sum = data["occurrence_sum"]
            pct_genes      = round(gene_count     / total_genes * 100, 4) if total_genes else 0
            pct_occ        = round(occurrence_sum / total_genes * 100, 4) if total_genes else 0

            row = {"species":          species,
                   "domain_id":        ipr_id,
                   "domain_name":      domain_name,
                   "gene_count":       gene_count,
                   "occurrence_sum":   occurrence_sum,
                   "total_genes":      total_genes,
                   "pct_genes":        pct_genes,
                   "pct_occurrences":  pct_occ,
                   "scale_factor":     scale_factors[species],
                   "go_terms":         data["go_terms"],
                   "superfamily":      data["superfamily"]}

            # Add grouping columns
            for col in GROUPING_COLS:
                row[col] = sp_grouping.get(col, "NA")

            rows.append(row)

    with open(args.output, "w") as fh:
        fh.write("\t".join(out_cols) + "\n")
        for row in rows:
            fh.write("\t".join(str(row.get(c, "NA")) for c in out_cols) + "\n")

    print(f"\nDone. {len(rows)} rows written to {args.output}")
    if unmatched:
        print(f"Unmatched TSV files: {unmatched}")


if __name__ == "__main__":
    main()

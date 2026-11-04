#!/usr/bin/env python3
"""
Merge per-lineage dbCAN CAZyme stats TSVs into one wide-format table.

For each CAZyme family ID in a provided list, finds matching rows across all
input TSVs and merges them into a single output table with one column set per
lineage. Families with all-zero means across every lineage are dropped.


Usage:
    python3 all_dbcan_stats.py \\
        --tsv-list   paths/to/dbcan_stats TSV.txt \\
        --cazyme-list optional/list/of/cazyme_ids.txt \\
        --output     combined_dbcan_stats.tsv \\
        [--labels    Metamonada,Amoebozoa,Apicomplexa]

if cazyme_ids.txt is not included, all families present in any TSV are included

"""

import argparse
import csv
import os
import re
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--tsv-list",    required=True,
                        help="Text file with paths to lineage dbcan_stats TSVs (one per line)")
    parser.add_argument("--cazyme-list", default=None,
                        help="Optional text file with CAZyme family IDs to include "
                             "(one per line). If omitted, all families are included.")
    parser.add_argument("--output",      required=True,
                        help="Output TSV path")
    parser.add_argument("--labels",      default=None,
                        help="Optional comma-separated labels for each TSV in the same "
                             "order as tsv-list (e.g. 'Metamonada,Amoebozoa'). "
                             "If not provided, labels are inferred from filenames.")
    return parser.parse_args()


def infer_label(filepath):
    """
    Infer a short lineage label from the filename.
    e.g. Metamonada_dbcan_stats.tsv → 'Metamonada'
         metamonad_stats.tsv        → 'metamonad'
    Strips common suffixes, returns first underscore-delimited token.
    """
    base = os.path.basename(filepath)
    base = re.sub(r'\.(tsv|txt|csv)$', '', base, flags=re.IGNORECASE)
    for suffix in ['_dbcan_stats', '_dbcan', '_stats', '_cazyme']:
        base = re.sub(suffix + '$', '', base, flags=re.IGNORECASE)
    return base.split('_')[0]


def load_tsv(filepath):
    """
    Load a dbcan_stats TSV into {cazyme_family_id: row_dict}.
    Returns (header_cols, family_dict).
    """
    family_dict = {}
    with open(filepath) as f:
        reader = csv.DictReader(f, delimiter="\t")
        header_cols = list(reader.fieldnames or [])
        for row in reader:
            fid = row.get("cazyme_family_id", "").strip()
            if fid:
                family_dict[fid] = dict(row)
    return header_cols, family_dict


def is_zero(val):
    """Return True if val is 0.0 or non-numeric (treated as missing)."""
    try:
        return float(val) == 0.0
    except (TypeError, ValueError):
        return True


def main():
    args = parse_args()

    # Load TSV paths
    with open(args.tsv_list) as f:
        tsv_paths = [l.strip() for l in f if l.strip()]
    if not tsv_paths:
        sys.exit("ERROR: --tsv-list is empty.")

    # Load optional CAZyme ID filter list
    cazyme_filter = None
    if args.cazyme_list:
        with open(args.cazyme_list) as f:
            cazyme_filter = [l.strip() for l in f if l.strip()]
        if not cazyme_filter:
            sys.exit("ERROR: --cazyme-list is empty.")
        print(f"Filtering to {len(cazyme_filter)} CAZyme families from {args.cazyme_list}")

    # Determine labels
    if args.labels:
        labels = [l.strip() for l in args.labels.split(",")]
        if len(labels) != len(tsv_paths):
            sys.exit(f"ERROR: --labels has {len(labels)} entries "
                     f"but tsv-list has {len(tsv_paths)}.")
    else:
        labels = [infer_label(p) for p in tsv_paths]
        print(f"Inferred lineage labels: {labels}")
        print("  (Override with --labels L1,L2,... if these are wrong)")

    # Load all TSVs
    lineage_data = []   # list of (label, header_cols, family_dict)
    all_families = []   # ordered union of all family IDs seen

    for path, label in zip(tsv_paths, labels):
        if not os.path.exists(path):
            sys.exit(f"ERROR: File not found: {path}")
        header_cols, family_dict = load_tsv(path)
        lineage_data.append((label, header_cols, family_dict))
        print(f"  Loaded '{label}': {len(family_dict)} families "
              f"from {os.path.basename(path)}")
        for fid in family_dict:
            if fid not in all_families:
                all_families.append(fid)

    # Apply ID filter or use full union
    target_families = cazyme_filter if cazyme_filter else all_families

    # Shared columns (not repeated per lineage)
    skip_cols = {"cazyme_family_id"}

    # Build output header
    out_header = ["cazyme_family_id"]
    for label, header_cols, _ in lineage_data:
        for col in header_cols:
            if col not in skip_cols:
                out_header.append(f"{label}_{col}")

    # Build one row per CAZyme family
    rows = []
    not_found = []
    all_zero_dropped = []

    for fid in target_families:
        fid = fid.strip()
        if not fid:
            continue

        found_in_any = any(fid in fd for _, _, fd in lineage_data)
        if not found_in_any:
            not_found.append(fid)
            continue

        row = {"cazyme_family_id": fid}

        for label, header_cols, family_dict in lineage_data:
            stat_cols = [c for c in header_cols if c not in skip_cols]
            if fid in family_dict:
                for col in stat_cols:
                    row[f"{label}_{col}"] = family_dict[fid].get(col, "NA")
            else:
                for col in stat_cols:
                    row[f"{label}_{col}"] = "NA"

        # Drop families where all mean_* values are zero/NA across every lineage
        all_zero = all(
            all(is_zero(row.get(f"{label}_{col}"))
                for col in header_cols
                if col.startswith("mean_"))
            for label, header_cols, _ in lineage_data
        )
        if all_zero:
            all_zero_dropped.append(fid)
            continue

        rows.append(row)

    # Write output
    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=out_header, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    # Summary
    n_found = sum(1 for fid in target_families
                  if any(fid in fd for _, _, fd in lineage_data))
    n_all   = sum(1 for fid in target_families
                  if all(fid in fd for _, _, fd in lineage_data))

    print(f"\nDone. {len(rows)} CAZyme family rows written to {args.output}")
    print(f"  Present in at least one lineage : {n_found}/{len(target_families)}")
    print(f"  Present in ALL {len(lineage_data)} lineages    : {n_all}/{len(target_families)}")
    if all_zero_dropped:
        print(f"  Dropped (all-zero means)        : {len(all_zero_dropped)}")
    if not_found:
        preview = ", ".join(not_found[:10])
        ellipsis = "..." if len(not_found) > 10 else ""
        print(f"  Not found in any lineage ({len(not_found)}): {preview}{ellipsis}")


if __name__ == "__main__":
    main()

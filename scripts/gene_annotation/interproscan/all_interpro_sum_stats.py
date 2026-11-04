#!/usr/bin/env python3

"""
Compare domain statistics across multiple lineage TSV files.
For each domain_id in a provided list, finds matching rows across all input TSVs
and merges them into one output table.

Usage:
  python3 all_interpro_sum_stats.py \
      --tsv_list paths/to/each/lineage/Mann-Whitney or t-test/tsv.txt \
      --ipr_list list/of/ipr/accessions/one/per/line.txt \
      --output combined_lineages.tsv

Note: 
    - Lineage label is inferred from the filename, or set with --labels
"""

import argparse
import csv
import os
import re
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description="Merge per-lineage domain stats into one wide-format table.")
    parser.add_argument("--tsv_list", required=True,
                        help="Text file with paths to lineage TSV files (one per line)")
    parser.add_argument("--ipr_list", required=True,
                        help="Text file with IPR domain IDs to include (one per line)")
    parser.add_argument("--output", required=True,
                        help="Output TSV path")
    parser.add_argument("--labels", default=None,
                        help="Optional comma-separated labels for each TSV in the same order "
                             "as tsv_list (e.g. 'S,M,A'). If not provided, labels are "
                             "inferred from filenames.")
    return parser.parse_args()


def infer_label(filepath):
    """
    Infer a short lineage label from the filename.
    e.g. metamonad_mannwhitney.tsv  -> 'metamonad'
         amoebozoan_fl_ha_mw.tsv   -> 'amoebozoan'
    Strips common suffixes and returns the first underscore-delimited token.
    """
    base = os.path.basename(filepath)
    base = re.sub(r'\.(tsv|txt|csv)$', '', base, flags=re.IGNORECASE)
    for suffix in ['_mannwhitney', '_mw', '_ttest', '_stats', '_fl_ha', '_results']:
        base = re.sub(suffix + '$', '', base, flags=re.IGNORECASE)
    return base.split('_')[0]


def load_tsv(filepath):
    """
    Load a lineage TSV into a dict: {domain_id: row_dict}.
    Returns (header_cols, domain_dict).
    """
    domain_dict = {}
    header_cols = []
    with open(filepath) as f:
        reader = csv.DictReader(f, delimiter="\t")
        header_cols = list(reader.fieldnames or [])
        for row in reader:
            domain_id = row.get("domain_id", "").strip()
            if domain_id:
                domain_dict[domain_id] = dict(row)
    return header_cols, domain_dict


def main():
    args = parse_args()

    # Load TSV paths
    with open(args.tsv_list) as f:
        tsv_paths = [l.strip() for l in f if l.strip()]
    if not tsv_paths:
        sys.exit("ERROR: tsv_list is empty.")

    # Load IPR domain IDs
    with open(args.ipr_list) as f:
        ipr_ids = [l.strip() for l in f if l.strip()]
    if not ipr_ids:
        sys.exit("ERROR: ipr_list is empty.")

    # Determine labels
    if args.labels:
        labels = [l.strip() for l in args.labels.split(",")]
        if len(labels) != len(tsv_paths):
            sys.exit(f"ERROR: --labels has {len(labels)} entries but tsv_list has {len(tsv_paths)}.")
    else:
        labels = [infer_label(p) for p in tsv_paths]
        print(f"Inferred lineage labels: {labels}")
        print("  (Override with --labels L1,L2,L3 if these are wrong)")

    # Load all TSVs
    lineage_data = []   # list of (label, header_cols, domain_dict)
    for path, label in zip(tsv_paths, labels):
        if not os.path.exists(path):
            sys.exit(f"ERROR: File not found: {path}")
        header_cols, domain_dict = load_tsv(path)
        lineage_data.append((label, header_cols, domain_dict))
        print(f"  Loaded '{label}': {len(domain_dict)} domains from {os.path.basename(path)}")

    # Columns shared across all lineages (not repeated per lineage)
    skip_cols = {"domain_id", "domain_name"}

    # Build output header:
    # domain_id | domain_name | <label>_col ... for each lineage
    out_header = ["domain_id", "domain_name"]
    for label, header_cols, _ in lineage_data:
        for col in header_cols:
            if col not in skip_cols:
                out_header.append(f"{label}_{col}")

    # Build one output row per IPR domain
    rows = []
    not_found_in_any = []

    for ipr_id in ipr_ids:
        ipr_id = ipr_id.strip()
        if not ipr_id:
            continue

        # Get domain_name from whichever lineage has it
        domain_name = "NA"
        found_in_any = False
        for _, _, domain_dict in lineage_data:
            if ipr_id in domain_dict:
                domain_name = domain_dict[ipr_id].get("domain_name", "NA")
                found_in_any = True
                break

        if not found_in_any:
            not_found_in_any.append(ipr_id)
            continue  # skip — not present in any lineage

        row = {"domain_id": ipr_id, "domain_name": domain_name}

        for label, header_cols, domain_dict in lineage_data:
            stat_cols = [c for c in header_cols if c not in skip_cols]
            if ipr_id in domain_dict:
                for col in stat_cols:
                    row[f"{label}_{col}"] = domain_dict[ipr_id].get(col, "NA")
            else:
                for col in stat_cols:
                    row[f"{label}_{col}"] = "NA"

        # Exclude domains where mean_FL and mean_HA are both 0 (or NA)
        # across every lineage — these have no signal anywhere
        def is_zero(val):
            try:
                return float(val) == 0.0
            except (TypeError, ValueError):
                return True  # treat NA/missing as zero for this check

        # Check if all mean values are zero/NA across every lineage.
        # Detect column names dynamically rather than assuming mean_FL/mean_HA.
        all_zero = all(
            all(is_zero(row.get(f"{label}_{col}"))
                for col in header_cols
                if col.startswith("mean_") and col not in {"domain_id", "domain_name"})
            for label, header_cols, _ in lineage_data
        )
        if all_zero:
            not_found_in_any.append(ipr_id)
            continue

        rows.append(row)

    # Write output
    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=out_header, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    # Summary stats
    n_found = sum(
        1 for ipr_id in ipr_ids
        if any(ipr_id in d for _, _, d in lineage_data)
    )
    n_all = sum(
        1 for ipr_id in ipr_ids
        if all(ipr_id in d for _, _, d in lineage_data)
    )

    print(f"\nDone! {len(rows)} domain rows written to {args.output}")
    print(f"  Present in at least one lineage : {n_found}/{len(ipr_ids)}")
    print(f"  Present in ALL {len(lineage_data)} lineages    : {n_all}/{len(ipr_ids)}")
    if not_found_in_any:
        preview = ', '.join(not_found_in_any[:10])
        ellipsis = '...' if len(not_found_in_any) > 10 else ''
        print(f"  Not found in any lineage ({len(not_found_in_any)}): {preview}{ellipsis}")


if __name__ == "__main__":
    main()

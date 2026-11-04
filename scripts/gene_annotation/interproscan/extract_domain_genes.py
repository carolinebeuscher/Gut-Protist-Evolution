#!/usr/bin/env python3
"""
Finds genes carrying InterPro domains that were significant in
merge_interpro_stats' output (e.g. all_lineages_mannwhitney.tsv), by scanning
per-sample {sample}_interprosum.tsv files (from clean_interproscan) for those
domain IDs.

A domain counts as significant if its "{lineage}_significant_bh" column is
"yes" for at least one lineage (--significant-in any, the default) or for
every lineage present (--significant-in all). Restrict which lineages count
with --lineages; default is every lineage found in --stats-tsv.

Output is one gene ID per line — the same format orthogroup_stats.py's
--genes-of-interest expects, so this feeds directly into it.

Usage:
    python3 extract_domain_genes.py \\
        --stats-tsv         all_lineages_mannwhitney.tsv \\
        --interpro-tsv-list list_of_interprosum_tsv_paths.txt \\
        --output            domain_significant_genes.txt \\
        [--significant-in   any|all] \\
        [--lineages         label1 label2 ...]
"""

import argparse
import csv
import os
import re
import sys

# interprosum.tsv's Interpro column can exceed the csv module's default
# 131072-byte field limit for a heavily-annotated gene — raise it as high
# as this platform allows.
_maxint = sys.maxsize
while True:
    try:
        csv.field_size_limit(_maxint)
        break
    except OverflowError:
        _maxint = int(_maxint / 10)


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--stats-tsv", required=True,
                        help="all_lineages_mannwhitney.tsv (or _ttest.tsv) "
                             "from merge_interpro_stats")
    parser.add_argument("--interpro-tsv-list", required=True,
                        help="Text file with paths to {sample}_interprosum.tsv "
                             "files (one per line)")
    parser.add_argument("--output", required=True,
                        help="Output gene list path (one gene ID per line)")
    parser.add_argument("--significant-in", choices=["any", "all"], default="any",
                        help="Require significance in 'any' or 'all' lineage "
                             "columns (default: any)")
    parser.add_argument("--lineages", nargs="+", default=None,
                        help="Restrict to specific lineage labels (default: "
                             "every lineage found in --stats-tsv)")
    return parser.parse_args()


def find_significant_domains(stats_tsv, mode, requested_lineages):
    with open(stats_tsv) as f:
        reader = csv.DictReader(f, delimiter="\t")
        fieldnames = reader.fieldnames or []
        sig_cols_all = [c for c in fieldnames if c.endswith("_significant_bh")]
        if not sig_cols_all:
            sys.exit(f"ERROR: no '*_significant_bh' columns found in {stats_tsv}")

        found_lineages = [c[:-len("_significant_bh")] for c in sig_cols_all]
        if requested_lineages:
            missing = set(requested_lineages) - set(found_lineages)
            if missing:
                sys.exit(f"ERROR: lineage(s) {sorted(missing)} not found in "
                         f"{stats_tsv}. Available: {found_lineages}")
            use_cols = [f"{l}_significant_bh" for l in requested_lineages]
        else:
            use_cols = sig_cols_all

        print(f"Checking significance across lineage(s) "
              f"{[c[:-len('_significant_bh')] for c in use_cols]} "
              f"(mode={mode})", file=sys.stderr)

        significant = set()
        for row in reader:
            flags = [row.get(c, "NA") == "yes" for c in use_cols]
            is_sig = any(flags) if mode == "any" else all(flags)
            if is_sig:
                significant.add(row["domain_id"].strip())

    return significant


def parse_interpro_field(interpro_str):
    """
    Extract IPR IDs from a clean_interpro.py 'Interpro' cell, e.g.
    "IPR002126-Cadherin-like_dom,IPR015919-Cadherin-like_sf (2)".
    Mirrors interpro_summary.py's parse_interpro_field().
    """
    interpro_str = interpro_str.strip()
    if not interpro_str or interpro_str == "-":
        return []
    ids = []
    for entry in re.split(r',(?=IPR\d+)', interpro_str):
        entry = entry.strip()
        if not entry:
            continue
        entry = re.sub(r'\s*\(\d+\)\s*$', '', entry)
        ipr_id = entry.split("-", 1)[0].strip()
        if re.match(r'^IPR\d+$', ipr_id):
            ids.append(ipr_id)
    return ids


def main():
    args = parse_args()

    significant_domains = find_significant_domains(
        args.stats_tsv, args.significant_in, args.lineages)
    print(f"{len(significant_domains)} significant domain(s) found",
          file=sys.stderr)

    if not significant_domains:
        print("No significant domains — writing empty gene list", file=sys.stderr)
        open(args.output, "w").close()
        return

    with open(args.interpro_tsv_list) as f:
        tsv_paths = [l.strip() for l in f if l.strip()]

    genes = set()
    for tsv_path in tsv_paths:
        if not os.path.exists(tsv_path):
            print(f"WARNING: {tsv_path} not found — skipping", file=sys.stderr)
            continue
        with open(tsv_path) as f:
            reader = csv.DictReader(f, delimiter="\t")
            if "Interpro" not in (reader.fieldnames or []):
                print(f"WARNING: no 'Interpro' column in {tsv_path} — skipping",
                      file=sys.stderr)
                continue
            for row in reader:
                ipr_ids = parse_interpro_field(row.get("Interpro", ""))
                if significant_domains.intersection(ipr_ids):
                    genes.add(row["Gene"].strip())

    with open(args.output, "w") as fh:
        for gene in sorted(genes):
            fh.write(gene + "\n")

    print(f"Done. {len(genes)} gene(s) written to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()

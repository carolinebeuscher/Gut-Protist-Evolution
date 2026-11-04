#!/usr/bin/env python3

"""
Add InterPro domain annotations and dbCAN CAZyme annotations to an alien
index TSV by matching the alien 'query' column to:
  - 'Gene' column in an InterProScan summary TSV
  - 'gene' column in a dbCAN overview txt file

Note:
    --dbcan is optional; omit it to skip dbCAN annotation
    --cols option lets you add other IPR columns (ex:  --cols Interpro Pfam)

Output is named like: Gint_alien_ipr.tsv

Usage:
  python3 add_ipr_to_alien.py \
      --alien Gint_alien.tsv \
      --ipr Gint_prot.fasta_sum.tsv \
      --dbcan Gint_prot_overview_dbcan.txt
"""

import argparse
import csv
import os
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description="Add InterProScan and dbCAN domain columns to an alien index TSV.")
    parser.add_argument("--alien", nargs="+", required=True,
                        help="One or more alien index TSV files")
    parser.add_argument("--ipr", nargs="+", required=True,
                        help="One or more InterProScan summary TSV files, "
                             "in the same order as --alien")
    parser.add_argument("--dbcan", nargs="+", default=None,
                        help="One or more dbCAN overview txt files (gene<tab>dbcan_id), "
                             "in the same order as --alien. Optional.")
    parser.add_argument("--cols", nargs="+",
                        default=["Interpro"],
                        help="Column(s) from the InterProScan TSV to add "
                             "(default: Interpro). Options include: CDD-IPR, "
                             "Gene3D, Hamap, PANTHER, Pfam, PIRSF, PRINTS, "
                             "ProSitePatterns, ProSiteProfiles, SFLD, SMART, "
                             "SUPERFAMILY, TIGRFAM, Interpro")
    parser.add_argument("--outdir", default=None,
                        help="Output directory (default: same dir as each alien file)")
    return parser.parse_args()


def load_ipr(ipr_path, cols):
    """
    Load InterProScan summary TSV into a dict: {gene_id: {col: value}}.
    Only loads the requested columns.
    """
    gene_map = {}
    with open(ipr_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        available = set(reader.fieldnames or [])
        missing = [c for c in cols if c not in available]
        if missing:
            print(f"  WARNING: Column(s) not found in {os.path.basename(ipr_path)}: "
                  f"{missing}. Available: {sorted(available)}")
        cols_to_use = [c for c in cols if c in available]

        for row in reader:
            gene = row.get("Gene", "").strip()
            if not gene:
                continue
            gene_map[gene] = {c: row.get(c, "-").strip() for c in cols_to_use}

    return gene_map, cols_to_use


def load_dbcan(dbcan_path):
    """
    Load a dbCAN overview txt file (gene<tab>dbcan_id) into a dict:
    {gene_id: dbcan_id}. 
    A gene with multiple rows (multiple CAZyme families) are joined with commas.
    """
    gene_map = {}
    with open(dbcan_path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) < 2:
                continue
            gene = cols[0].strip()
            cazyme = cols[1].strip()

            # Skip header if present
            if gene.lower() in ("gene", "gene_id", "#gene"):
                continue
            if not cazyme or cazyme == "-":
                continue

            if gene in gene_map:
                gene_map[gene] = gene_map[gene] + "," + cazyme
            else:
                gene_map[gene] = cazyme

    return gene_map


def process_trio(alien_path, ipr_path, dbcan_path, cols, outdir):
    print(f"\nProcessing: {os.path.basename(alien_path)}")

    # Load IPR annotations
    print(f"  Loading IPR file: {os.path.basename(ipr_path)}")
    ipr_map, cols_found = load_ipr(ipr_path, cols)
    print(f"  {len(ipr_map)} genes loaded from IPR file.")

    if not cols_found:
        print(f"  ERROR: None of the requested IPR columns found. Skipping.")
        return

    # Load dbCAN annotations if provided
    dbcan_map = {}
    if dbcan_path:
        print(f"  Loading dbCAN file: {os.path.basename(dbcan_path)}")
        dbcan_map = load_dbcan(dbcan_path)
        print(f"  {len(dbcan_map)} genes with dbCAN annotations loaded.")

    # Build output path
    alien_stem = os.path.splitext(os.path.basename(alien_path))[0]  # e.g. Gint_alien
    out_name   = f"{alien_stem}_ipr.tsv"
    out_dir    = outdir if outdir else os.path.dirname(os.path.abspath(alien_path))
    out_path   = os.path.join(out_dir, out_name)

    ipr_matched     = 0
    ipr_unmatched   = 0
    dbcan_matched   = 0
    dbcan_unmatched = 0

    with open(alien_path) as infile, open(out_path, "w") as outfile:
        reader = csv.DictReader(infile, delimiter="\t")
        if not reader.fieldnames:
            print(f"  ERROR: Could not read header from {alien_path}. Skipping.")
            return

        # Always add dbcan column if a dbcan file was provided
        extra_cols = cols_found + (["dbcan"] if dbcan_path else [])
        out_fields = list(reader.fieldnames) + extra_cols

        writer = csv.DictWriter(outfile, fieldnames=out_fields, delimiter="\t",
                                extrasaction="ignore", lineterminator="\n")
        writer.writeheader()

        for row in reader:
            query = row.get("query", "").strip()

            # Add IPR columns
            ipr_data = ipr_map.get(query)
            if ipr_data:
                row.update(ipr_data)
                ipr_matched += 1
            else:
                for c in cols_found:
                    row[c] = "-"
                ipr_unmatched += 1

            # Add dbCAN column
            if dbcan_path:
                dbcan_val = dbcan_map.get(query)
                if dbcan_val:
                    row["dbcan"] = dbcan_val
                    dbcan_matched += 1
                else:
                    row["dbcan"] = "-"
                    dbcan_unmatched += 1

            writer.writerow(row)

    print(f"  IPR   — matched: {ipr_matched} | unmatched: {ipr_unmatched}")
    if dbcan_path:
        print(f"  dbCAN — matched: {dbcan_matched} | unmatched: {dbcan_unmatched}")
    print(f"  Written to: {out_path}")


def main():
    args = parse_args()

    if len(args.alien) != len(args.ipr):
        sys.exit(f"ERROR: --alien and --ipr must have the same number of files "
                 f"({len(args.alien)} alien vs {len(args.ipr)} ipr provided).")

    if args.dbcan and len(args.dbcan) != len(args.alien):
        sys.exit(f"ERROR: --dbcan must have the same number of files as --alien "
                 f"({len(args.alien)} alien vs {len(args.dbcan)} dbcan provided).")

    # Validate all input files exist
    all_files = args.alien + args.ipr + (args.dbcan or [])
    for path in all_files:
        if not os.path.exists(path):
            sys.exit(f"ERROR: File not found: {path}")

    if args.outdir:
        os.makedirs(args.outdir, exist_ok=True)

    dbcan_files = args.dbcan if args.dbcan else [None] * len(args.alien)

    for alien_path, ipr_path, dbcan_path in zip(args.alien, args.ipr, dbcan_files):
        process_trio(alien_path, ipr_path, dbcan_path, cols=args.cols, outdir=args.outdir)

    print("\nDone!")


if __name__ == "__main__":
    main()

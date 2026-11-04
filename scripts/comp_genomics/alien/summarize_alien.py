#!/usr/bin/env python3
"""
Summarizes AlienIndex output TSV files.

For queries with alien_score > THRESHOLD (default 0.2):
  - Count queries whose skipped_hits contain species OUTSIDE the ingroup
  - Count queries whose skipped_hits contain only ingroup species
  - Count queries with NO entries in skipped_hits
  - Write a new TSV of the "outside ingroup" queries

Skipped-hits cell format (pipe-delimited fields, comma-separated hits):
  hitID|species name|class|bitscore,hitID|species name|class|bitscore,...

Usage
-----
  python summarize_alien_output.py \\
      --input  results.alien.tsv \\
      --ingroup Streptophyta \\
      [--threshold 0.2] \\
      [--output outside_ingroup.tsv]

The --ingroup string is case-insensitive and can have multiple values
    --ingroup "Streptophyta,Embryophyta"
    - matched as a substring against the CLASS field of ecah skipped hit

Can use --ingroup-species to skip specific species that lack a proper NCBI
class (e.g. "No phylum or order in NCBI taxid") but are still ingroup members
"""

import argparse
import csv
import sys
from pathlib import Path


EMPTY_MARKERS = {"", "-", "--", "N/A", "NA", "none", "None"}


def parse_skipped_hits(cell: str) -> list[dict]:
    """
    Parse the skipped_hits cell into a list of dicts with keys:
      hit_id, species, tax_class, bitscore
    Returns an empty list when the cell is blank / '--'.
    """
    cell = cell.strip()
    if cell in EMPTY_MARKERS:
        return []

    hits = []
    for entry in cell.split(","):
        entry = entry.strip()
        if not entry:
            continue
        parts = entry.split("|")
        if len(parts) < 3:
            # malformed – keep it but mark class as unknown
            hits.append({"hit_id": parts[0] if parts else entry,
                         "species": parts[1] if len(parts) > 1 else "",
                         "tax_class": "UNKNOWN",
                         "bitscore": ""})
        else:
            hits.append({
                "hit_id":    parts[0],
                "species":   parts[1],
                "tax_class": parts[2],
                "bitscore":  parts[3] if len(parts) > 3 else "",
            })
    return hits


def is_ingroup(tax_class: str, species: str,
               ingroup_terms: list[str], ingroup_species: list[str]) -> bool:
    """
    Return True if the hit belongs to the ingroup by either:
      - tax_class matching any ingroup class term (case-insensitive substring), OR
      - species matching any ingroup_species term (case-insensitive substring)
    The species fallback catches entries whose class is 'No phylum or order in
    NCBI taxid' or similarly uninformative.
    """
    tc = tax_class.lower()
    sp = species.lower()
    return (any(term in tc for term in ingroup_terms) or
            any(term in sp for term in ingroup_species))



def main():
    parser = argparse.ArgumentParser(
        description="Summarise alien output TSV and categorise skipped hits."
    )
    parser.add_argument("--input",  "-i", required=True,
                        help="Path to the alien output .tsv file")
    parser.add_argument("--ingroup", required=True,
                        help="Ingroup taxon string(s) to match against the class "
                             "field of skipped hits. Comma-separated for multiple "
                             "terms, e.g. 'Streptophyta,Embryophyta'")
    parser.add_argument("--ingroup-species", default="",
                        help="Additional species name(s) to treat as ingroup, "
                             "matched against the species field. Useful for taxa "
                             "that lack a proper NCBI class. Comma-separated, "
                             "e.g. 'Blastocystis hominis,Blastocystis sp. NMH'")
    parser.add_argument("--threshold", "-t", type=float, default=0.2,
                        help="Alien score threshold (default: 0.2)")
    parser.add_argument("--output", "-o", default=None,
                        help="Output TSV for queries hitting outside the ingroup "
                             "(default: <input>_outside_ingroup.tsv)")
    args = parser.parse_args()

    ingroup_terms = [t.strip().lower() for t in args.ingroup.split(",") if t.strip()]
    ingroup_species = [t.strip().lower() for t in args.ingroup_species.split(",") if t.strip()]
    input_path = Path(args.input)
    if not input_path.exists():
        sys.exit(f"ERROR: input file not found: {input_path}")

    output_path = Path(args.output) if args.output else \
        input_path.with_name(input_path.stem + "_outside_ingroup.tsv")

    #read & classify
    total_rows          = 0
    above_threshold     = 0
    outside_ingroup_rows = []
    ingroup_only_count  = 0
    no_hits_count       = 0
    outside_count       = 0

    fieldnames = None

    with open(input_path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = reader.fieldnames

        if fieldnames is None:
            sys.exit("ERROR: could not read header from input file.")

        # Normalise column names (strip whitespace)
        fieldnames = [f.strip() for f in fieldnames]

        # Find the alien_score and skipped_hits columns (case-insensitive)
        col_map = {f.lower(): f for f in fieldnames}
        score_col   = col_map.get("alien_score")
        skipped_col = col_map.get("skipped_hits")

        if score_col is None:
            sys.exit("ERROR: could not find 'alien_score' column in the header.")
        if skipped_col is None:
            sys.exit("ERROR: could not find 'skipped_hits' column in the header.")

        for raw_row in reader:
            # Normalise keys
            row = {k.strip(): v for k, v in raw_row.items() if k}
            total_rows += 1

            score_val = row.get(score_col, "").strip()
            try:
                score = float(score_val)
            except ValueError:
                continue  # skip rows with non-numeric scores

            if score <= args.threshold:
                continue

            above_threshold += 1
            skipped_cell = row.get(skipped_col, "").strip()
            hits = parse_skipped_hits(skipped_cell)

            if not hits:
                no_hits_count += 1
            else:
                has_outside = any(
                    not is_ingroup(h["tax_class"], h["species"],
                                   ingroup_terms, ingroup_species)
                    for h in hits
                )
                if has_outside:
                    outside_count += 1
                    outside_ingroup_rows.append(row)
                else:
                    ingroup_only_count += 1

    # write outside-ingroup table
    with open(output_path, "w", newline="", encoding="utf-8") as out_fh:
        writer = csv.DictWriter(out_fh, fieldnames=fieldnames, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(outside_ingroup_rows)

    # print summary
    #printing input_file	ingroup_term	alien_score	total_queries	queries_above_threshold	no_skipped_hits	skipped_ingroup_hits	skipped_outgroup_hits
    print(f"{input_path}\t{', '.join(ingroup_terms)}\t{args.threshold}\t{total_rows}\t{above_threshold}\t{no_hits_count}\t{ingroup_only_count}\t{outside_count}")


if __name__ == "__main__":
    main()

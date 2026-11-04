#!/usr/bin/env python3
"""
Parses a genome summary xlsx (sample name, genome assembly type, GTDB-style
taxonomy string) into a TSV with the same rank columns as grouping_key.tsv,
plus a genome_type column. Does NOT compute host_association/aerotolerance —
those are added later by add_clustering.py once this is combined with the
existing taxonomy data.

Expects columns (by header name, any order): "Genome name", "Genome assembly",
"Taxonomy" (GTDB format: d__Domain;p__Phylum;c__Class;o__Order;f__Family;
g__Genus;s__Genus species).

Sample names are written with "_" replaced by "-" to match genomes already
renamed for hyphen-only naming.

Output columns:
    sample_id, species, genus, family, order, class, phylum, domain,
    clade_1, clade_2, clade_3, genome_type

Usage:
    python3 parse_genome_xlsx.py \\
        RCG_450_summary.xlsx \\
        new_genomes_taxonomy.tsv

Install dependencies:
    pip install pandas openpyxl
"""

import sys

try:
    import pandas as pd
except ImportError:
    sys.exit("ERROR: run: pip install pandas openpyxl")


ASSEMBLY_TYPE_MAP = {
    "Monoculture genome":            "monoculture",
    "Metagenome-assembled genome":   "MAG",
    "Single-cell amplified genome":  "SAG",
}

GTDB_RANK_PREFIXES = {
    "d": "domain",
    "p": "phylum",
    "c": "class",
    "o": "order",
    "f": "family",
    "g": "genus",
    "s": "species",
}

OUT_COLS = [
    "sample_id", "species", "genus", "family", "order", "class",
    "phylum", "domain", "clade_1", "clade_2", "clade_3", "genome_type",
]


def parse_gtdb_taxonomy(tax_string, row_label):
    """
    Parse "d__X;p__X;c__X;o__X;f__X;g__X;s__X" into a rank dict.
    Raises on missing/unparseable ranks rather than silently writing NA,
    since a bad parse here would silently mis-taxonomize a sample.
    """
    ranks = {}
    for token in tax_string.split(";"):
        token = token.strip()
        if "__" not in token:
            sys.exit(f"ERROR: unparseable taxonomy token '{token}' for "
                      f"sample '{row_label}' (full string: {tax_string!r})")
        prefix, _, value = token.partition("__")
        rank_name = GTDB_RANK_PREFIXES.get(prefix)
        if rank_name is None:
            sys.exit(f"ERROR: unrecognized rank prefix '{prefix}__' for "
                      f"sample '{row_label}' (full string: {tax_string!r})")
        ranks[rank_name] = value.strip() or "NA"

    missing = [r for r in GTDB_RANK_PREFIXES.values() if r not in ranks]
    if missing:
        sys.exit(f"ERROR: taxonomy string for sample '{row_label}' is "
                  f"missing rank(s) {missing} (full string: {tax_string!r})")
    return ranks


def main():
    if len(sys.argv) < 3:
        sys.exit("Usage: parse_genome_xlsx.py <input.xlsx> <output.tsv>")

    xlsx_path  = sys.argv[1]
    output_path = sys.argv[2]

    df = pd.read_excel(xlsx_path, sheet_name=0, engine="openpyxl")

    required_cols = ["Genome name", "Genome assembly", "Taxonomy"]
    missing_cols = [c for c in required_cols if c not in df.columns]
    if missing_cols:
        sys.exit(f"ERROR: missing expected column(s) {missing_cols} in "
                  f"{xlsx_path}. Found columns: {list(df.columns)}")

    rows_out = []
    seen_ids = {}

    for _, row in df.iterrows():
        raw_name = str(row["Genome name"]).strip()
        if not raw_name or raw_name.lower() == "nan":
            sys.exit(f"ERROR: row with missing 'Genome name' in {xlsx_path}")
        sample_id = raw_name.replace("_", "-")

        if sample_id in seen_ids:
            sys.exit(f"ERROR: duplicate sample_id '{sample_id}' derived from "
                      f"'Genome name' values '{seen_ids[sample_id]}' and "
                      f"'{raw_name}' — check for a naming collision after "
                      f"underscore-to-hyphen conversion")
        seen_ids[sample_id] = raw_name

        assembly_raw = str(row["Genome assembly"]).strip()
        genome_type = ASSEMBLY_TYPE_MAP.get(assembly_raw)
        if genome_type is None:
            sys.exit(f"ERROR: unrecognized 'Genome assembly' value "
                      f"'{assembly_raw}' for sample '{raw_name}'. "
                      f"Known values: {list(ASSEMBLY_TYPE_MAP)}")

        tax_string = str(row["Taxonomy"]).strip()
        ranks = parse_gtdb_taxonomy(tax_string, raw_name)

        rows_out.append({
            "sample_id": sample_id,
            "species":   ranks["species"],
            "genus":     ranks["genus"],
            "family":    ranks["family"],
            "order":     ranks["order"],
            "class":     ranks["class"],
            "phylum":    ranks["phylum"],
            "domain":    ranks["domain"],
            "clade_1":   "NA",
            "clade_2":   "NA",
            "clade_3":   "NA",
            "genome_type": genome_type,
        })

    with open(output_path, "w") as fh:
        fh.write("\t".join(OUT_COLS) + "\n")
        for r in rows_out:
            fh.write("\t".join(str(r[c]) for c in OUT_COLS) + "\n")

    type_counts = {}
    for r in rows_out:
        type_counts[r["genome_type"]] = type_counts.get(r["genome_type"], 0) + 1

    print(f"Wrote {len(rows_out)} rows to {output_path}", file=sys.stderr)
    print(f"  genome_type counts: {type_counts}", file=sys.stderr)


if __name__ == "__main__":
    main()

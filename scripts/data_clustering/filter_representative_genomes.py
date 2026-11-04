#!/usr/bin/env python3
"""
Filters the genome summary xlsx down to one representative genome per species
— rows where "Representative genome for each species" == "Yes" — before any
of the other new-genome scripts run. This is now step 0 of the workflow:

    filter_representative_genomes.py  (this script)
    generate_sample_ids.py
    rename_xlsx_genomes.py
    parse_genome_xlsx.py
    combine_taxonomy.py

Errors out (rather than silently proceeding) if any species has zero or more
than one "Yes" row, since either would mean the filtered set no longer has
exactly one genome per species.

Usage:
    python3 filter_representative_genomes.py \\
        RCG_450_summary.xlsx \\
        RCG_450_summary_representative.xlsx

Install dependencies:
    pip install pandas openpyxl
"""

import sys

try:
    import pandas as pd
except ImportError:
    sys.exit("ERROR: run: pip install pandas openpyxl")

REP_COL = "Representative genome for each species"


def main():
    if len(sys.argv) != 3:
        sys.exit("Usage: filter_representative_genomes.py <input.xlsx> "
                  "<output.xlsx>")

    xlsx_path   = sys.argv[1]
    output_path = sys.argv[2]

    df = pd.read_excel(xlsx_path, sheet_name=0, engine="openpyxl")

    for col in (REP_COL, "Taxonomy", "Genome name"):
        if col not in df.columns:
            sys.exit(f"ERROR: missing expected column '{col}' in {xlsx_path}. "
                      f"Found columns: {list(df.columns)}")

    filtered = df[df[REP_COL] == "Yes"].copy()
    if filtered.empty:
        sys.exit(f"ERROR: no rows with {REP_COL!r} == 'Yes' found in {xlsx_path}")

    species = filtered["Taxonomy"].str.extract(r"s__([^;]+)")[0].str.strip()
    counts = species.value_counts()
    dupes = counts[counts > 1]
    if not dupes.empty:
        sys.exit(f"ERROR: {len(dupes)} species have more than one "
                  f"'{REP_COL}' == 'Yes' row — expected exactly one per "
                  f"species: {dupes.to_dict()}")

    filtered.to_excel(output_path, sheet_name="Sheet1", index=False)

    print(f"Kept {len(filtered)}/{len(df)} rows (one representative genome "
          f"per species) -> {output_path}", file=sys.stderr)


if __name__ == "__main__":
    main()

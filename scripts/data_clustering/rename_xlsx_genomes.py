#!/usr/bin/env python3
"""
Replaces the "Genome name" column of the genome summary xlsx with the short
sample_ids from generate_sample_ids.py's rename map, so every later step
(parse_genome_xlsx.py, combine_taxonomy.py, and sample_key.tsv/add_clustering.py)
already sees the short codes as the sample name.

Writes a NEW xlsx rather than overwriting the input, so you can review before
replacing the original. Every "Genome name" value in the input must have an
entry in the rename map — this errors out (rather than leaving a row
unrenamed) if any are missing, which would otherwise silently produce a mixed
long-name/short-code file.

Usage:
    python3 rename_xlsx_genomes.py \\
        RCG_450_summary.xlsx \\
        new_genome_rename_map.tsv \\
        RCG_450_summary_renamed.xlsx

Install dependencies:
    pip install pandas openpyxl
"""

import sys

try:
    import pandas as pd
except ImportError:
    sys.exit("ERROR: run: pip install pandas openpyxl")


def read_rename_map(path):
    mapping = {}
    with open(path) as fh:
        header = fh.readline()  # "original_name\tnew_sample_id"
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            original, new_id = line.split("\t")
            mapping[original] = new_id
    return mapping


def main():
    if len(sys.argv) < 4:
        sys.exit("Usage: rename_xlsx_genomes.py <input.xlsx> "
                  "<rename_map.tsv> <output.xlsx>")

    xlsx_path   = sys.argv[1]
    map_path    = sys.argv[2]
    output_path = sys.argv[3]

    mapping = read_rename_map(map_path)

    df = pd.read_excel(xlsx_path, sheet_name=0, engine="openpyxl")
    if "Genome name" not in df.columns:
        sys.exit(f"ERROR: 'Genome name' column not found in {xlsx_path}. "
                  f"Found columns: {list(df.columns)}")

    original_names = df["Genome name"].astype(str).str.strip()
    missing = sorted(set(original_names) - set(mapping))
    if missing:
        sys.exit(f"ERROR: {len(missing)} 'Genome name' value(s) in {xlsx_path} "
                  f"have no entry in {map_path} — refusing to write a "
                  f"partially-renamed file. Missing: {missing[:10]}"
                  f"{' ...' if len(missing) > 10 else ''}")

    df["Genome name"] = original_names.map(mapping)
    df.to_excel(output_path, sheet_name="Sheet1", index=False)

    print(f"Renamed {len(mapping)} genome names -> {output_path}", file=sys.stderr)
    print(f"  (original {xlsx_path} left untouched)", file=sys.stderr)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Left-joins two TSV files on specified key columns.

All rows from the left file are kept. Matching rows from the right file are
appended. Non-matching rows get NA for all right-file columns.
Duplicate key columns from the right file are dropped.

Usage:
    python3 join_tsv.py \\
        left.tsv \\
        right.tsv \\
        left_key_col \\
        right_key_col \\
        output.tsv
"""

import sys


def load_tsv(path):
    """Return (header_list, [{col: val, ...}, ...])."""
    with open(path) as fh:
        lines = [l.rstrip("\n") for l in fh if l.strip()]
    if not lines:
        sys.exit(f"ERROR: {path} is empty")
    header = lines[0].split("\t")
    rows   = [dict(zip(header, l.split("\t"))) for l in lines[1:]]
    return header, rows


def main():
    left_path      = sys.argv[1]
    right_path     = sys.argv[2]
    left_key_col   = sys.argv[3]
    right_key_col  = sys.argv[4]
    output_path    = sys.argv[5]

    left_header,  left_rows  = load_tsv(left_path)
    right_header, right_rows = load_tsv(right_path)

    if left_key_col not in left_header:
        sys.exit(f"ERROR: '{left_key_col}' not found in {left_path}\n"
                 f"  Columns: {left_header}")
    if right_key_col not in right_header:
        sys.exit(f"ERROR: '{right_key_col}' not found in {right_path}\n"
                 f"  Columns: {right_header}")

    # Build lookup: right_key_value -> row dict
    right_lookup = {}
    for row in right_rows:
        key = row[right_key_col]
        right_lookup[key] = row

    # Right-file columns to append (skip the key col to avoid duplication)
    right_extra_cols = [c for c in right_header if c != right_key_col]

    # Output header: all left cols + right extra cols
    out_header = left_header + right_extra_cols

    unmatched = []

    with open(output_path, "w") as fh:
        fh.write("\t".join(out_header) + "\n")

        for row in left_rows:
            key        = row.get(left_key_col, "")
            right_row  = right_lookup.get(key, {})

            if not right_row:
                unmatched.append(key)

            out_vals = [row.get(c, "NA") for c in left_header]
            out_vals += [right_row.get(c, "NA") for c in right_extra_cols]
            fh.write("\t".join(out_vals) + "\n")

    if unmatched:
        print(f"WARNING: {len(unmatched)} sample(s) in left file had no match "
              f"in {right_path}:", file=sys.stderr)
        for s in unmatched:
            print(f"  {s}", file=sys.stderr)

    print(f"Wrote {output_path} ({len(left_rows)} rows, "
          f"{len(out_header)} columns)", file=sys.stderr)


if __name__ == "__main__":
    main()

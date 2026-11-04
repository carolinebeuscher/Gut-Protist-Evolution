#!/usr/bin/env python3
"""
summarize_alien.py wrapper that:
  1. Reads grouping_key.tsv to find each sample's ingroup term from a
     configurable taxonomic column (e.g. order, phylum, class)
  2. Runs summarize_alien.py for each sample, using that term as --ingroup
  3. Combines per-sample summary stats into one TSV
  4. Combines per-sample outside-ingroup detail rows into one TSV
     with sample_id and all grouping columns prepended

Both output tables include all grouping columns from grouping_key.tsv!

Usage:
    python3 run_alien_summary.py \\
        --alien-files    list_of_alien_tsvs.txt \\
        --grouping-key   grouping_key.tsv \\
        --ingroup-col    order \\
        --summary-out    alien_summary.tsv \\
        --detail-out     alien_outside_ingroup.tsv \\
        [--threshold     0.2] \\
        [--ingroup-species "Blastocystis hominis,Blastocystis sp. NMH"]
"""

import argparse
import csv
import os
import subprocess
import sys
import tempfile


GROUPING_COLS = ["host_association", "aerotolerance", "phylum", "class", "order", "family"]


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--alien-files",     required=True,
                        help="Text file with paths to per-sample alien output TSVs")
    parser.add_argument("--grouping-key",    required=True,
                        help="grouping_key.tsv from add_clustering rule")
    parser.add_argument("--ingroup-col",     required=True,
                        help="Column in grouping_key.tsv to use as ingroup term "
                             "(e.g. order, phylum, class)")
    parser.add_argument("--summary-out",     required=True,
                        help="Output TSV for per-sample summary statistics")
    parser.add_argument("--detail-out",      required=True,
                        help="Output TSV for combined outside-ingroup detail rows")
    parser.add_argument("--threshold",       type=float, default=0.2,
                        help="Alien score threshold passed to summarize_alien.py "
                             "(default: 0.2)")
    parser.add_argument("--ingroup-species", default="",
                        help="Additional species to treat as ingroup regardless of "
                             "class (passed to summarize_alien.py --ingroup-species)")
    return parser.parse_args()


def load_grouping_key(path, ingroup_col):
    """
    Load grouping_key.tsv into {sample_id: {col: value}}.
    Returns (grouping dict, list of available grouping cols in file).
    """
    grouping = {}
    available_cols = []

    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        available_cols = [c for c in (reader.fieldnames or [])
                          if c in GROUPING_COLS]
        if ingroup_col not in (reader.fieldnames or []):
            sys.exit(f"ERROR: --ingroup-col '{ingroup_col}' not found in "
                     f"{path}.\n  Available columns: {reader.fieldnames}")
        for row in reader:
            sid = row.get("sample_id", "").strip() or row.get("species", "").strip()
            if sid:
                grouping[sid] = {c: row.get(c, "NA").strip()
                                 for c in GROUPING_COLS
                                 if c in (reader.fieldnames or [])}

    return grouping, available_cols


def infer_sample_id(alien_path):
    """
    Infer sample ID from alien TSV filename.
    """
    base = os.path.basename(alien_path)
    for suffix in (".alien.tsv", "_alien.tsv", ".tsv"):
        if base.endswith(suffix):
            return base[: -len(suffix)]
    return base


def run_summarize_alien(alien_path, ingroup_term, threshold,
                        ingroup_species, detail_tmp):
    """
    Run summarize_alien.py for one sample. Returns the summary line (str)
    or None on failure.
    """
    script = os.path.join(os.path.dirname(__file__), "summarize_alien.py")
    cmd = [
        sys.executable, script,
        "--input",     alien_path,
        "--ingroup",   ingroup_term,
        "--threshold", str(threshold),
        "--output",    detail_tmp,
    ]
    if ingroup_species:
        cmd += ["--ingroup-species", ingroup_species]

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"  WARNING: summarize_alien.py failed for {alien_path}:\n"
              f"    {result.stderr.strip()}", file=sys.stderr)
        return None

    # summarize_alien.py prints one summary line to stdout
    line = result.stdout.strip()
    return line if line else None


def main():
    args = parse_args()

    # Load file list
    with open(args.alien_files) as f:
        alien_paths = [l.strip() for l in f if l.strip()]
    if not alien_paths:
        sys.exit("ERROR: --alien-files is empty")

    # Load grouping key
    print("Loading grouping key...", file=sys.stderr)
    grouping, available_cols = load_grouping_key(args.grouping_key, args.ingroup_col)
    print(f"  {len(grouping)} samples loaded, ingroup column: '{args.ingroup_col}'",
          file=sys.stderr)

    # Summary table header
    summary_cols = (["sample_id"] + available_cols +
                    ["ingroup_term", "alien_threshold",
                     "total_queries", "queries_above_threshold",
                     "no_skipped_hits", "skipped_ingroup_hits",
                     "skipped_outgroup_hits"])

    # Detail table header is built from the first alien TSV + prepended cols
    detail_header_written = False
    detail_prepend_cols   = ["sample_id"] + available_cols

    summary_rows  = []
    unmatched     = []

    with open(args.summary_out, "w") as sum_fh, \
         open(args.detail_out,  "w") as det_fh:

        sum_writer = csv.DictWriter(sum_fh, fieldnames=summary_cols,
                                    delimiter="\t", extrasaction="ignore")
        sum_writer.writeheader()

        for alien_path in alien_paths:
            sample_id = infer_sample_id(alien_path)
            print(f"  Processing {sample_id}...", file=sys.stderr)

            sp_grouping = grouping.get(sample_id, {})
            if not sp_grouping:
                # Try matching by prefix in case names don't align exactly
                match = next((k for k in grouping if sample_id.startswith(k)
                               or k.startswith(sample_id)), None)
                if match:
                    sp_grouping = grouping[match]
                    sample_id   = match
                else:
                    unmatched.append(sample_id)
                    print(f"  WARNING: '{sample_id}' not found in grouping_key.tsv "
                          f"— grouping columns will be NA", file=sys.stderr)

            ingroup_term = sp_grouping.get(args.ingroup_col, "")
            if not ingroup_term or ingroup_term == "NA":
                print(f"  WARNING: no ingroup term for '{sample_id}' in column "
                      f"'{args.ingroup_col}' — skipping", file=sys.stderr)
                continue

            # Run summarize_alien.py with a temp file for the detail output
            with tempfile.NamedTemporaryFile(suffix=".tsv", delete=False,
                                             mode="w") as tmp:
                tmp_path = tmp.name

            try:
                summary_line = run_summarize_alien(
                    alien_path, ingroup_term, args.threshold,
                    args.ingroup_species, tmp_path)

                if summary_line is None:
                    continue

                # Parse the summary line:
                # input_file\tingroup_term\talien_score\ttotal\tabove\tno_hits\tingroup\toutside
                parts = summary_line.split("\t")
                if len(parts) < 8:
                    print(f"  WARNING: unexpected summary format for {sample_id}",
                          file=sys.stderr)
                    continue

                sum_row = {"sample_id": sample_id,
                           "ingroup_term":              parts[1],
                           "alien_threshold":           parts[2],
                           "total_queries":             parts[3],
                           "queries_above_threshold":   parts[4],
                           "no_skipped_hits":           parts[5],
                           "skipped_ingroup_hits":      parts[6],
                           "skipped_outgroup_hits":     parts[7]}
                for col in available_cols:
                    sum_row[col] = sp_grouping.get(col, "NA")
                sum_writer.writerow(sum_row)

                # Append detail rows with prepended sample + grouping columns
                if os.path.exists(tmp_path) and os.path.getsize(tmp_path) > 0:
                    with open(tmp_path) as tmp_in:
                        det_reader = csv.DictReader(tmp_in, delimiter="\t")
                        det_fieldnames = det_reader.fieldnames or []

                        # Write detail header once
                        if not detail_header_written:
                            full_header = detail_prepend_cols + det_fieldnames
                            det_fh.write("\t".join(full_header) + "\n")
                            detail_header_written = True

                        for det_row in det_reader:
                            prepend_vals = [sample_id] + [
                                sp_grouping.get(c, "NA") for c in available_cols]
                            det_vals = [det_row.get(c, "NA")
                                        for c in det_fieldnames]
                            det_fh.write(
                                "\t".join(prepend_vals + det_vals) + "\n")

            finally:
                if os.path.exists(tmp_path):
                    os.unlink(tmp_path)

    # Write an empty detail header if no samples produced detail rows
    if not detail_header_written:
        with open(args.detail_out, "w") as det_fh:
            det_fh.write("\t".join(detail_prepend_cols) + "\n")

    if unmatched:
        print(f"\nWARNING: {len(unmatched)} samples not found in grouping_key.tsv: "
              f"{unmatched}", file=sys.stderr)

    print(f"\nDone.\n"
          f"  Summary : {args.summary_out}\n"
          f"  Detail  : {args.detail_out}", file=sys.stderr)


if __name__ == "__main__":
    main()

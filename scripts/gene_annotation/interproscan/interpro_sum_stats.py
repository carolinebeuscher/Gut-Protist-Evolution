#!/usr/bin/env python3
"""
Compares gene_count for InterPro domains between two groups (e.g. FL vs HA)
within a taxonomic subset of the InterProScan summary table.

Runs both Welch's t-test and Mann-Whitney U test, with Benjamini-Hochberg
FDR correction applied across all domains.

Outputs two TSVs:
    {output_prefix}_ttest.tsv
    {output_prefix}_mannwhitney.tsv

Usage:
    python3 interpro_sum_stats.py \\
        --input        all_interpro_summary.tsv \\
        --group        Metamonada \\
        --group-col    phylum \\
        --compare-col  host_association \\
        --output       metamonad_fl_ha

  --group-col   : column to subset on          (default: phylum)
  --group       : value inside group-col that you want to sort by 
  --compare-col : column with two group labels (default: host_association)
"""

import argparse
import csv
import math
import statistics
import sys
from collections import defaultdict

try:
    from scipy import stats as scipy_stats
except ImportError:
    sys.exit("ERROR: scipy is required. Run: pip install scipy")


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input",       required=True,
                        help="all_interpro_summary.tsv from summarize_interproscan")
    parser.add_argument("--group",       required=True,
                        help="Value of --group-col to subset on (e.g. Metamonada)")
    parser.add_argument("--group-col",   default="phylum",
                        help="Column to subset by (default: phylum)")
    parser.add_argument("--compare-col", default="host_association",
                        help="Column whose two values are compared (default: host_association)")
    parser.add_argument("--output",      required=True,
                        help="Output prefix; produces {prefix}_ttest.tsv and {prefix}_mannwhitney.tsv")
    return parser.parse_args()



def benjamini_hochberg(p_values):
    """BH FDR correction. None values (untestable) are passed through."""
    indexed = [(i, p) for i, p in enumerate(p_values) if p is not None]
    n = len(indexed)
    adjusted = [None] * len(p_values)
    if n == 0:
        return adjusted
    indexed.sort(key=lambda x: x[1])
    bh_values = [p * n / (rank + 1) for rank, (_, p) in enumerate(indexed)]
    running_min = 1.0
    for rank in range(n - 1, -1, -1):
        running_min = min(running_min, bh_values[rank])
        bh_values[rank] = running_min
    for rank, (orig_i, _) in enumerate(indexed):
        adjusted[orig_i] = round(min(bh_values[rank], 1.0), 6)
    return adjusted


def summary_stats(values):
    """Return (n, mean, median, min, max)."""
    n = len(values)
    if n == 0:
        return 0, None, None, None, None
    return (n,
            round(statistics.mean(values), 4),
            round(statistics.median(values), 4),
            min(values),
            max(values))



def main():
    args = parse_args()

    # Load and subset to the requested group
    domain_data  = defaultdict(lambda: defaultdict(list))  # domain_id -> {label: [counts]}
    domain_names = {}
    domain_go    = {}
    domain_sf    = {}
    compare_labels = set()   # the two values of --compare-col we see

    with open(args.input) as f:
        reader = csv.DictReader(f, delimiter="\t")

        # Validate columns exist
        fieldnames = reader.fieldnames or []
        for col in (args.group_col, args.compare_col):
            if col not in fieldnames:
                sys.exit(f"ERROR: column '{col}' not found in {args.input}.\n"
                         f"  Available columns: {fieldnames}")

        for row in reader:
            if row[args.group_col].strip().lower() != args.group.strip().lower():
                continue

            domain_id  = row["domain_id"].strip()
            label      = row[args.compare_col].strip()
            gene_count = int(row["gene_count"])

            compare_labels.add(label)
            domain_data[domain_id][label].append(gene_count)
            domain_names[domain_id] = row["domain_name"]

            if domain_id not in domain_go:
                domain_go[domain_id] = row.get("go_terms",   "NA") or "NA"
            if domain_id not in domain_sf:
                domain_sf[domain_id] = row.get("superfamily", "NA") or "NA"

    if not domain_data:
        sys.exit(f"ERROR: No rows found for {args.group_col}='{args.group}'.\n"
                 f"  Check the value matches your column exactly (case-insensitive).")

    if len(compare_labels) != 2:
        print(f"WARNING: --compare-col '{args.compare_col}' must have exactly 2 unique "
              f"values in the subset for '{args.group}', found: {compare_labels}. "
              f"Writing empty outputs.", file=sys.stderr)
        # Use whatever label(s) exist so the header is still meaningful
        found = sorted(compare_labels)
        la = found[0] if len(found) > 0 else "A"
        lb = found[1] if len(found) > 1 else "B"
        shared_cols = [
            "domain_id", "domain_name",
            f"n_{la}", f"n_{lb}",
            f"mean_{la}", f"mean_{lb}",
            f"median_{la}", f"median_{lb}",
            f"min_{la}", f"min_{lb}",
            f"max_{la}", f"max_{lb}",
            "pi_enrichment", "log2fc", "go_terms", "superfamily",
        ]
        ttest_path = f"{args.output}_ttest.tsv"
        mw_path    = f"{args.output}_mannwhitney.tsv"
        with open(ttest_path, "w") as f:
            f.write("\t".join(shared_cols + ["t_statistic", "t_pvalue", "t_pvalue_bh", "significant_bh", "note"]) + "\n")
        with open(mw_path, "w") as f:
            f.write("\t".join(shared_cols + ["mw_statistic", "mw_pvalue", "mw_pvalue_bh", "significant_bh", "note"]) + "\n")
        sys.exit(0)


    label_a, label_b = sorted(compare_labels)
    print(f"Group       : {args.group_col} = '{args.group}'")
    print(f"Comparing   : '{label_a}' vs '{label_b}' ({args.compare_col})")
    print(f"Domains     : {len(domain_data)}")

    # Run tests per domain
    results = []

    for domain_id, counts in sorted(domain_data.items()):
        vals_a = counts.get(label_a, [])
        vals_b = counts.get(label_b, [])

        n_a, mean_a, median_a, min_a, max_a = summary_stats(vals_a)
        n_b, mean_b, median_b, min_b, max_b = summary_stats(vals_b)

        all_identical = len(set(vals_a + vals_b)) == 1

        # Welch's t-test
        if n_a >= 2 and n_b >= 2 and not all_identical:
            t_stat, t_p = scipy_stats.ttest_ind(vals_a, vals_b, equal_var=False)
            t_stat, t_p, t_note = round(float(t_stat), 4), round(float(t_p), 6), ""
        elif all_identical:
            t_stat, t_p, t_note = None, None, "all values identical"
        else:
            t_stat, t_p, t_note = None, None, \
                f"insufficient data (n_{label_a}={n_a}, n_{label_b}={n_b}; need >=2 per group)"

        # Mann-Whitney U
        if n_a >= 1 and n_b >= 1 and not all_identical:
            try:
                mw_stat, mw_p = scipy_stats.mannwhitneyu(vals_a, vals_b, alternative="two-sided")
                mw_stat, mw_p, mw_note = round(float(mw_stat), 4), round(float(mw_p), 6), ""
            except ValueError as e:
                mw_stat, mw_p, mw_note = None, None, str(e)
        elif all_identical:
            mw_stat, mw_p, mw_note = None, None, "all values identical"
        else:
            mw_stat, mw_p, mw_note = None, None, \
                f"insufficient data (n_{label_a}={n_a}, n_{label_b}={n_b})"

        # Enrichment (b relative to a)
        if mean_a is None or mean_b is None:
            pi_enrichment = log2fc = None
        elif median_a == 0:
            pi_enrichment = round(mean_b, 4)
            log2fc        = round(math.log2((mean_b + 1) / (mean_a + 1)), 4)
        else:
            pi_enrichment = round(mean_b / mean_a, 4)
            log2fc        = round(math.log2((mean_b + 1) / (mean_a + 1)), 4)

        results.append({
            "domain_id":             domain_id,
            "domain_name":           domain_names[domain_id],
            f"n_{label_a}":          n_a,
            f"n_{label_b}":          n_b,
            f"mean_{label_a}":       mean_a,
            f"mean_{label_b}":       mean_b,
            f"median_{label_a}":     median_a,
            f"median_{label_b}":     median_b,
            f"min_{label_a}":        min_a,
            f"min_{label_b}":        min_b,
            f"max_{label_a}":        max_a,
            f"max_{label_b}":        max_b,
            "pi_enrichment":         pi_enrichment,
            "log2fc":                log2fc,
            "go_terms":              domain_go.get(domain_id, "NA"),
            "superfamily":           domain_sf.get(domain_id, "NA"),
            "t_stat":                t_stat,
            "t_pvalue":              t_p,
            "t_note":                t_note,
            "mw_stat":               mw_stat,
            "mw_pvalue":             mw_p,
            "mw_note":               mw_note,
        })

    # BH correction
    t_bh  = benjamini_hochberg([r["t_pvalue"]  for r in results])
    mw_bh = benjamini_hochberg([r["mw_pvalue"] for r in results])
    for r, tb, mb in zip(results, t_bh, mw_bh):
        r["t_pvalue_bh"]  = tb
        r["mw_pvalue_bh"] = mb

    def fmt(v):
        return "NA" if v is None else v

    def sig_flag(bh_p):
        if bh_p is None: return "NA"
        return "yes" if bh_p < 0.05 else "no"

    shared_cols = [
        "domain_id", "domain_name",
        f"n_{label_a}",      f"n_{label_b}",
        f"mean_{label_a}",   f"mean_{label_b}",
        f"median_{label_a}", f"median_{label_b}",
        f"min_{label_a}",    f"min_{label_b}",
        f"max_{label_a}",    f"max_{label_b}",
        "pi_enrichment", "log2fc",
        "go_terms", "superfamily",
    ]

    # Write t-test output
    ttest_path = f"{args.output}_ttest.tsv"
    with open(ttest_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(shared_cols + ["t_statistic", "t_pvalue", "t_pvalue_bh", "significant_bh", "note"])
        for r in results:
            writer.writerow([fmt(r[c]) for c in shared_cols] + [
                fmt(r["t_stat"]), fmt(r["t_pvalue"]), fmt(r["t_pvalue_bh"]),
                sig_flag(r["t_pvalue_bh"]), r["t_note"],
            ])

    # Write Mann-Whitney output
    mw_path = f"{args.output}_mannwhitney.tsv"
    with open(mw_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(shared_cols + ["mw_statistic", "mw_pvalue", "mw_pvalue_bh", "significant_bh", "note"])
        for r in results:
            writer.writerow([fmt(r[c]) for c in shared_cols] + [
                fmt(r["mw_stat"]), fmt(r["mw_pvalue"]), fmt(r["mw_pvalue_bh"]),
                sig_flag(r["mw_pvalue_bh"]), r["mw_note"],
            ])

    t_sig  = sum(1 for r in results if r["t_pvalue_bh"]  is not None and r["t_pvalue_bh"]  < 0.05)
    mw_sig = sum(1 for r in results if r["mw_pvalue_bh"] is not None and r["mw_pvalue_bh"] < 0.05)
    print(f"T-test results       → {ttest_path}  ({t_sig} significant BH<0.05)")
    print(f"Mann-Whitney results → {mw_path}  ({mw_sig} significant BH<0.05)")


if __name__ == "__main__":
    main()

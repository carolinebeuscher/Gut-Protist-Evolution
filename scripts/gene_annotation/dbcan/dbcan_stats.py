#!/usr/bin/env python3
"""
Compute FL vs HA (or any two-group comparison) summary statistics for each
CAZyme family within a taxonomic subset, using pct_genes as the metric.
Outputs a wide-format table with stat columns per comparison group.

Requires the output of combine_dbcan.py as input.

Usage:
  python3 dbcan_stats.py \\
      --input      cazyme_combined.tsv \\
      --output     cazyme_stats.tsv \\
      --group      Metamonada \\
      --group-col  clade_1 \\
      --compare-col host_association

  --group-col   : column to subset rows on           (default: phylum)
  --group : the value fo group-col you want to filter by 
  --compare-col : column with two labels to compare  (default: host_association)

"""

import argparse
import csv
import math
import sys
from collections import defaultdict

try:
    import numpy as np
    from scipy import stats
except ImportError:
    sys.exit("ERROR: scipy and numpy are required: pip install scipy numpy")


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input",        required=True,
                        help="Combined CAZyme TSV from combine_dbcan.py")
    parser.add_argument("--output",       required=True,
                        help="Output wide-format stats TSV")
    parser.add_argument("--group",        required=True,
                        help="Value of --group-col to subset on (e.g. Metamonada)")
    parser.add_argument("--group-col",    default="phylum",
                        help="Column to subset by (default: phylum)")
    parser.add_argument("--compare-col",  default="host_association",
                        help="Column whose two values are compared (default: host_association)")
    parser.add_argument("--metric",       default="pct_genes",
                        choices=["pct_genes", "pct_occurrences"],
                        help="Metric to compare (default: pct_genes)")
    parser.add_argument("--alpha",        type=float, default=0.05,
                        help="Significance threshold after BH correction (default: 0.05)")
    return parser.parse_args()


def bh_correction(pvalues):
    """Benjamini-Hochberg FDR correction. Returns adjusted p-values in input order."""
    n = len(pvalues)
    if n == 0:
        return []
    indexed = sorted(enumerate(pvalues), key=lambda x: x[1])
    adjusted = [0.0] * n
    prev = 1.0
    for rank, (orig_idx, pval) in enumerate(reversed(indexed), start=1):
        adj = min(pval * n / (n - rank + 1), prev)
        prev = adj
        adjusted[orig_idx] = adj
    return adjusted


def safe_log2fc(mean_a, mean_b):
    """log2((mean_a + 1e-6) / (mean_b + 1e-6))."""
    return round(math.log2((mean_a + 1e-6) / (mean_b + 1e-6)), 6)


def compute_group_stats(vals_a, vals_b, label_a, label_b):
    """
    Summary stats and Mann-Whitney U for two groups.
    Returns a dict keyed with the actual label names so output columns
    are e.g. mean_HA / mean_FL rather than mean_A / mean_B.
    Raw p-value is stored; BH correction is applied later.
    """
    a = np.array(vals_a, dtype=float)
    b = np.array(vals_b, dtype=float)

    def _stats(arr, label):
        n = len(arr)
        return {
            f"n_{label}":      n,
            f"mean_{label}":   round(float(np.mean(arr)),   6) if n > 0 else float("nan"),
            f"median_{label}": round(float(np.median(arr)), 6) if n > 0 else float("nan"),
            f"min_{label}":    round(float(np.min(arr)),    6) if n > 0 else float("nan"),
            f"max_{label}":    round(float(np.max(arr)),    6) if n > 0 else float("nan"),
        }

    result = {}
    result.update(_stats(a, label_a))
    result.update(_stats(b, label_b))

    mean_a = result[f"mean_{label_a}"]
    mean_b = result[f"mean_{label_b}"]
    result["log2fc"] = (safe_log2fc(mean_a, mean_b)
                        if not (math.isnan(mean_a) or math.isnan(mean_b))
                        else float("nan"))

    if len(a) >= 1 and len(b) >= 1:
        try:
            mw_stat, mw_pval = stats.mannwhitneyu(a, b, alternative="two-sided")
        except Exception:
            mw_stat, mw_pval = float("nan"), float("nan")
    else:
        mw_stat, mw_pval = float("nan"), float("nan")

    result["mw_statistic"] = round(float(mw_stat), 4) if not math.isnan(mw_stat) else "NA"
    result["mw_pvalue"]    = mw_pval   # raw — BH applied after all families processed

    return result


def fmt(val, decimals=6):
    if isinstance(val, float) and math.isnan(val):
        return "NA"
    if isinstance(val, float):
        return f"{val:.{decimals}f}"
    return str(val)


def main():
    args = parse_args()

    # Load combined CAZyme table 
    # data[cazyme_family_id][label] = [metric_value, ...]
    data          = defaultdict(lambda: defaultdict(list))
    cazymes_seen  = set()
    compare_labels = set()

    with open(args.input) as f:
        reader = csv.DictReader(f, delimiter="\t")
        fieldnames = reader.fieldnames or []

        for col in (args.group_col, args.compare_col, args.metric):
            if col not in fieldnames:
                sys.exit(f"ERROR: column '{col}' not found in {args.input}.\n"
                         f"  Available columns: {fieldnames}")

        for row in reader:
            if row[args.group_col].strip().lower() != args.group.strip().lower():
                continue

            cazyme_id = row["cazyme_family_id"].strip()
            label     = row[args.compare_col].strip()

            try:
                value = float(row[args.metric])
            except (ValueError, KeyError):
                continue

            data[cazyme_id][label].append(value)
            cazymes_seen.add(cazyme_id)
            compare_labels.add(label)

    if not cazymes_seen:
        sys.exit(f"ERROR: no rows found for {args.group_col}='{args.group}'.\n"
                 f"  Check the value matches your column exactly (case-insensitive).")

    if len(compare_labels) != 2:
        print(f"WARNING: --compare-col '{args.compare_col}' must have exactly 2 unique "
              f"values in the subset for '{args.group}', found: {compare_labels}. "
              f"Writing empty output.", file=sys.stderr)
        # Write empty output with header so Snakemake finds the expected file
        # Use whatever label(s) exist so the header is still meaningful
        found = sorted(compare_labels)
        la = found[0] if len(found) > 0 else "A"
        lb = found[1] if len(found) > 1 else "B"
        empty_header = ["cazyme_family_id",
                        f"n_{la}", f"n_{lb}",
                        f"mean_{la}", f"mean_{lb}",
                        f"median_{la}", f"median_{lb}",
                        f"min_{la}", f"min_{lb}",
                        f"max_{la}", f"max_{lb}",
                        "log2fc", "mw_statistic", "mw_pvalue", "mw_pvalue_bh", "significant_bh"]
        with open(args.output, "w") as f:
            f.write("\t".join(empty_header) + "\n")
        sys.exit(0)

    label_a, label_b = sorted(compare_labels)
    cazymes = sorted(cazymes_seen)

    print(f"Group       : {args.group_col} = '{args.group}'")
    print(f"Comparing   : '{label_a}' vs '{label_b}' ({args.compare_col})")
    print(f"CAZyme families: {len(cazymes)}")

    # Compute raw stats
    raw_stats = {}
    for cazyme_id in cazymes:
        vals_a = data[cazyme_id].get(label_a, [])
        vals_b = data[cazyme_id].get(label_b, [])
        raw_stats[cazyme_id] = compute_group_stats(vals_a, vals_b, label_a, label_b)

    # BH correction across all CAZyme families 
    valid = [(caz, raw_stats[caz]["mw_pvalue"])
             for caz in cazymes
             if isinstance(raw_stats[caz]["mw_pvalue"], float)
             and not math.isnan(raw_stats[caz]["mw_pvalue"])]

    if valid:
        adj = bh_correction([p for _, p in valid])
        for (caz, _), adj_p in zip(valid, adj):
            raw_stats[caz]["mw_pvalue_bh"]   = adj_p
            raw_stats[caz]["significant_bh"] = "yes" if adj_p < args.alpha else "no"

    for caz in cazymes:
        if "mw_pvalue_bh" not in raw_stats[caz]:
            raw_stats[caz]["mw_pvalue_bh"]   = "NA"
            raw_stats[caz]["significant_bh"] = "NA"
        # Format raw p-value now that BH is done
        pv = raw_stats[caz]["mw_pvalue"]
        raw_stats[caz]["mw_pvalue"] = fmt(pv) if not (isinstance(pv, float) and math.isnan(pv)) else "NA"

    # Build output
    stat_suffixes = (
        [f"n_{label_a}", f"n_{label_b}",
         f"mean_{label_a}", f"mean_{label_b}",
         f"median_{label_a}", f"median_{label_b}",
         f"min_{label_a}", f"min_{label_b}",
         f"max_{label_a}", f"max_{label_b}",
         "log2fc", "mw_statistic", "mw_pvalue", "mw_pvalue_bh", "significant_bh"]
    )

    out_header = ["cazyme_family_id"] + stat_suffixes

    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=out_header, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        for cazyme_id in cazymes:
            s   = raw_stats[cazyme_id]
            row = {"cazyme_family_id": cazyme_id}
            for suffix in stat_suffixes:
                val = s.get(suffix, "NA")
                row[suffix] = (fmt(val) if isinstance(val, float)
                               and not math.isnan(val) else str(val)
                               if not (isinstance(val, float) and math.isnan(val))
                               else "NA")
            writer.writerow(row)

    n_sig = sum(1 for caz in cazymes
                if raw_stats[caz].get("significant_bh") == "yes")
    print(f"\nDone. {len(cazymes)} CAZyme families written to {args.output}")
    print(f"  Significant (BH < {args.alpha}): {n_sig}")


if __name__ == "__main__":
    main()

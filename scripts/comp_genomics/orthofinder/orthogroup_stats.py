#!/usr/bin/env python3
"""
Compares orthogroup gene counts across groupings using:
    - Mann-Whitney U for 2 groups
    - Kruskal-Wallis for > 2 groups
with Benjamini-Hochberg FDR correction across all orthogroups.

Groupings can be any column in grouping_key.tsv 

Output: one TSV per grouping column, written to --outdir. The column set
depends on how many unique labels that grouping column has:

Usage:
    python3 orthogroup_stats.py \\
        --orthogroups   Orthogroups.GeneCount.tsv \\
        --grouping-key  grouping_key.tsv \\
        --outdir        orthogroup_stats/ \\
        --grouping-cols taxonomic_group host_association aerotolerance \\
        --min-samples   2 \\
        --min-nonzero   2

Install dependencies:
    pip install scipy
"""

import argparse
import csv
import math
import os
import statistics
import sys
from collections import defaultdict

# Orthogroups.tsv's per-species gene-list cells can exceed the csv module's
# default 131072-byte field limit for large multi-copy orthogroups — raise
# it as high as this platform allows.
_maxint = sys.maxsize
while True:
    try:
        csv.field_size_limit(_maxint)
        break
    except OverflowError:
        _maxint = int(_maxint / 10)

try:
    from scipy.stats import mannwhitneyu, kruskal
except ImportError:
    sys.exit("ERROR: scipy is required. Run: pip install scipy")


DEFAULT_GROUPING_COLS = ["taxonomic_group", "host_association", "aerotolerance",
                          "phylum", "class", "order"]


def benjamini_hochberg(p_values):
    """BH FDR correction. None values are passed through unchanged."""
    indexed = [(i, p) for i, p in enumerate(p_values) if p is not None]
    n = len(indexed)
    adjusted = [None] * len(p_values)
    if n == 0:
        return adjusted
    indexed.sort(key=lambda x: x[1])
    bh = [p * n / (rank + 1) for rank, (_, p) in enumerate(indexed)]
    running_min = 1.0
    for rank in range(n - 1, -1, -1):
        running_min = min(running_min, bh[rank])
        bh[rank] = running_min
    for rank, (orig_i, _) in enumerate(indexed):
        adjusted[orig_i] = round(min(bh[rank], 1.0), 6)
    return adjusted


### data loading ###

def compute_taxonomic_group(row):
    """
    Buckets a grouping_key.tsv row into Ciliates/Apicomplexa/Stramenopiles/
    Amoebozoa/Metamonada/Other via the same phylum/clade_1/clade_2 priority
    chain used in alien_sum_stats.py's load_sample_groups_and_ha(). A sample
    gets exactly one label from this mutually-exclusive if/elif chain.
    """
    phylum  = row.get("phylum",  "").strip()
    clade_1 = row.get("clade_1", "").strip()
    clade_2 = row.get("clade_2", "").strip()

    if phylum == "Apicomplexa":
        return "Apicomplexa"
    elif clade_1 == "Amoebozoa":
        return "Amoebozoa"
    elif clade_1 == "Metamonada":
        return "Metamonada"
    elif phylum == "Ciliophora":
        return "Ciliates"
    elif clade_2 == "Stramenopiles":
        return "Stramenopiles"
    else:
        return "Other"


def load_grouping_key(path, grouping_cols):
    """
    Load grouping_key.tsv into {sample_id: {col: value}}.
    Also accepts species as a fallback key if sample_id is absent.

    "taxonomic_group" is handled specially: it's never a literal column in
    grouping_key.tsv, so it's computed per row via compute_taxonomic_group()
    instead of read directly.
    """
    grouping = {}
    want_taxonomic_group = "taxonomic_group" in grouping_cols
    real_cols = [c for c in grouping_cols if c != "taxonomic_group"]

    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        available = set(reader.fieldnames or [])
        missing = [c for c in real_cols if c not in available]
        if missing:
            print(f"WARNING: grouping columns not found in grouping_key.tsv "
                  f"(will be NA): {missing}", file=sys.stderr)
        if want_taxonomic_group and not ({"phylum", "clade_1", "clade_2"} & available):
            print("WARNING: none of phylum/clade_1/clade_2 found in "
                  "grouping_key.tsv — 'taxonomic_group' will be 'Other' for "
                  "every sample", file=sys.stderr)

        for row in reader:
            # Prefer sample_id; fall back to species for matching column headers
            key = row.get("sample_id", "").strip() or row.get("species", "").strip()
            if key:
                entry = {c: row.get(c, "NA").strip()
                         for c in real_cols if c in available}
                if want_taxonomic_group:
                    entry["taxonomic_group"] = compute_taxonomic_group(row)
                grouping[key] = entry
    return grouping


def load_orthogroups(path):
    """
    Load Orthogroups.GeneCount.tsv into:
        species_ids : list of species column names (excluding 'Total' if present)
        og_counts   : {orthogroup_id: {species_id: count}}
    """
    species_ids = []
    og_counts   = {}

    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        # First column is 'Orthogroup'; last column is often 'Total' — drop it
        species_ids = [c for c in (reader.fieldnames or [])
                       if c not in ("Orthogroup", "Total")]
        for row in reader:
            og = row["Orthogroup"].strip()
            og_counts[og] = {sp: int(row.get(sp, 0) or 0)
                             for sp in species_ids}

    return species_ids, og_counts


### Statistics ##

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


def run_primary_test(valid_groups):
    """
    valid_groups: {label: [count, count, ...]}, >=2 labels, each already
    filtered to >= --min-samples.
    Mann-Whitney U for 2 groups, Kruskal-Wallis for >2.
    Returns (test_name, stat, p_value) — stat/p_value are None if untestable
    (all values identical, or a scipy ValueError).
    """
    labels = sorted(valid_groups.keys())
    arrays = [valid_groups[l] for l in labels]
    all_vals = [v for arr in arrays for v in arr]
    test_name = "mannwhitney" if len(labels) == 2 else "kruskal"

    if len(set(all_vals)) == 1:
        return test_name, None, None

    try:
        if len(labels) == 2:
            stat, p = mannwhitneyu(arrays[0], arrays[1], alternative="two-sided")
        else:
            stat, p = kruskal(*arrays)
        return test_name, round(float(stat), 4), round(float(p), 6)
    except ValueError:
        return test_name, None, None


def order_two_labels(labels):
    """
    Decide display/log2fc order for a 2-label grouping column. log2fc
    direction always follows alphabetical order (mean_later+1)/(mean_
    earlier+1) — same convention as interpro_sum_stats.py — which happens
    to put HA in the numerator for host_association, matching
    alien_sum_stats.py. Column display order puts "HA" first if present
    (matching that same HA-first convention used elsewhere in this
    pipeline), else alphabetical.
    """
    label_a, label_b = sorted(labels)  # alphabetical: used for log2fc direction
    if "HA" in (label_a, label_b):
        display_first  = "HA"
        display_second = label_a if label_b == "HA" else label_b
    else:
        display_first, display_second = label_a, label_b
    return label_a, label_b, display_first, display_second


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--orthogroups",   required=True,
                        help="Orthogroups.GeneCount.tsv from OrthoFinder")
    parser.add_argument("--grouping-key",  required=True,
                        help="grouping_key.tsv from add_clustering rule")
    parser.add_argument("--outdir",        required=True,
                        help="Directory to write one TSV per grouping column")
    parser.add_argument("--grouping-cols", nargs="+", default=DEFAULT_GROUPING_COLS,
                        help="Columns to test, or the special pseudo-column "
                             "'taxonomic_group' (default: {})".format(
                             " ".join(DEFAULT_GROUPING_COLS)))
    parser.add_argument("--min-samples",   type=int, default=2,
                        help="Min samples per group to run a test (default: 2)")
    parser.add_argument("--min-nonzero",   type=int, default=2,
                        help="Min samples with count > 0 across all groups "
                             "to include an orthogroup (default: 2)")
    parser.add_argument("--alpha",         type=float, default=0.05,
                        help="BH-corrected p-value cutoff for the "
                             "'significant_bh' column (default: 0.05)")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # Load data
    print("Loading grouping key...", file=sys.stderr)
    grouping = load_grouping_key(args.grouping_key, args.grouping_cols)
    print(f"  {len(grouping)} samples loaded", file=sys.stderr)

    print("Loading orthogroups...", file=sys.stderr)
    species_ids, og_counts = load_orthogroups(args.orthogroups)
    print(f"  {len(og_counts)} orthogroups × {len(species_ids)} species", file=sys.stderr)

    # Warn about species in OrthoFinder output not in grouping key
    unmatched = [sp for sp in species_ids if sp not in grouping]
    if unmatched:
        print(f"WARNING: {len(unmatched)} species in orthogroups not found in "
              f"grouping_key.tsv (will be excluded from tests): {unmatched}",
              file=sys.stderr)

    matched_species = [sp for sp in species_ids if sp in grouping]

    base_cols = ["orthogroup", "test", "groups_compared",
                 "stat_value", "p_value", "p_value_bh", "significant_bh"]

    for col in args.grouping_cols:
        print(f"Testing grouping column: {col}", file=sys.stderr)

        # Assign each matched species to a group label for this column
        species_label = {}
        for sp in matched_species:
            label = grouping[sp].get(col, "NA")
            if label and label != "NA":
                species_label[sp] = label

        label_counts = defaultdict(int)
        for label in species_label.values():
            label_counts[label] += 1
        qualifying_labels = sorted(l for l, n in label_counts.items()
                                   if n >= args.min_samples)

        outpath = os.path.join(args.outdir, f"orthogroups_{col}.tsv")

        if len(qualifying_labels) < 2:
            print(f"  Skipping '{col}' — fewer than 2 groups with >= "
                  f"--min-samples ({args.min_samples}) — writing empty output",
                  file=sys.stderr)
            with open(outpath, "w") as fh:
                fh.write("\t".join(base_cols) + "\n")
            continue

        two_group = len(qualifying_labels) == 2
        if two_group:
            label_a, label_b, disp1, disp2 = order_two_labels(qualifying_labels)
            out_cols = base_cols + [
                f"{disp1}_mean", f"{disp1}_median", f"{disp1}_min", f"{disp1}_max",
                f"{disp2}_mean", f"{disp2}_median", f"{disp2}_min", f"{disp2}_max",
                "log2fc",
            ]
        else:
            out_cols = base_cols + ["significant_group", "direction",
                                     "group_means", "group_medians",
                                     "group_mins", "group_maxs", "group_ns"]

        print(f"  Groups: {qualifying_labels}", file=sys.stderr)

        # Build results for all orthogroups
        results = []

        for og, counts in og_counts.items():
            groups = defaultdict(list)
            for sp, label in species_label.items():
                if label in qualifying_labels:
                    groups[label].append(counts.get(sp, 0))

            if len(groups) < 2:
                continue

            all_counts = [v for vals in groups.values() for v in vals]
            if sum(1 for v in all_counts if v > 0) < args.min_nonzero:
                continue

            test_name, stat, p = run_primary_test(groups)
            stats_by_label = {l: summary_stats(v) for l, v in groups.items()}

            row = {
                "orthogroup":      og,
                "test":            test_name,
                "groups_compared": ",".join(sorted(groups.keys())),
                "stat_value":      stat,
                "p_value":         p,
            }

            if two_group:
                _, mean_a, median_a, min_a, max_a = stats_by_label[label_a]
                _, mean_b, median_b, min_b, max_b = stats_by_label[label_b]
                vals_for_disp = {label_a: (mean_a, median_a, min_a, max_a),
                                 label_b: (mean_b, median_b, min_b, max_b)}
                mean1, median1, min1, max1 = vals_for_disp[disp1]
                mean2, median2, min2, max2 = vals_for_disp[disp2]
                row[f"{disp1}_mean"], row[f"{disp1}_median"] = mean1, median1
                row[f"{disp1}_min"],  row[f"{disp1}_max"]    = min1, max1
                row[f"{disp2}_mean"], row[f"{disp2}_median"] = mean2, median2
                row[f"{disp2}_min"],  row[f"{disp2}_max"]    = min2, max2
                row["log2fc"] = round(math.log2((mean_b + 1) / (mean_a + 1)), 4)
            else:
                means   = {l: s[1] for l, s in stats_by_label.items()}
                medians = {l: s[2] for l, s in stats_by_label.items()}
                mins    = {l: s[3] for l, s in stats_by_label.items()}
                maxs    = {l: s[4] for l, s in stats_by_label.items()}
                ns      = {l: s[0] for l, s in stats_by_label.items()}

                # Kruskal-Wallis is omnibus — it says whether the groups
                # differ at all, not which one. Identify the single group
                # whose mean deviates furthest from everyone else's pooled
                # mean (one-vs-rest deviation, not a separate hypothesis
                # test — significant_bh above is still the test to trust;
                # this just points at the group likely driving it).
                best_label, best_direction, best_abs_dev = "NA", "NA", -1
                for label in groups:
                    rest_vals = [v for l, vs in groups.items() if l != label for v in vs]
                    if not rest_vals:
                        continue
                    deviation = means[label] - statistics.mean(rest_vals)
                    if abs(deviation) > best_abs_dev:
                        best_abs_dev = abs(deviation)
                        best_label = label
                        best_direction = ("higher" if deviation > 0
                                          else ("lower" if deviation < 0 else "equal"))
                row["significant_group"] = best_label
                row["direction"]         = best_direction

                row["group_means"]   = ";".join(f"{l}={means[l]}"   for l in sorted(means))
                row["group_medians"] = ";".join(f"{l}={medians[l]}" for l in sorted(medians))
                row["group_mins"]    = ";".join(f"{l}={mins[l]}"    for l in sorted(mins))
                row["group_maxs"]    = ";".join(f"{l}={maxs[l]}"    for l in sorted(maxs))
                row["group_ns"]      = ";".join(f"{l}={ns[l]}"      for l in sorted(ns))

            results.append(row)

        if not results:
            print(f"  No testable orthogroups for '{col}' — writing empty output",
                  file=sys.stderr)
            with open(outpath, "w") as fh:
                fh.write("\t".join(out_cols) + "\n")
            continue

        bh = benjamini_hochberg([r["p_value"] for r in results])
        for r, bh_p in zip(results, bh):
            r["p_value_bh"] = bh_p
            r["significant_bh"] = ("yes" if bh_p is not None and bh_p < args.alpha
                                   else ("NA" if bh_p is None else "no"))

        def fmt(v):
            return "NA" if v is None else str(v)

        with open(outpath, "w") as fh:
            fh.write("\t".join(out_cols) + "\n")
            for r in results:
                fh.write("\t".join(fmt(r.get(c, "")) for c in out_cols) + "\n")

        n_sig = sum(1 for r in results if r["significant_bh"] == "yes")
        print(f"  {len(results)} orthogroups tested, {n_sig} significant "
              f"(BH<{args.alpha}) → {outpath}", file=sys.stderr)

    print("Done.", file=sys.stderr)


if __name__ == "__main__":
    main()

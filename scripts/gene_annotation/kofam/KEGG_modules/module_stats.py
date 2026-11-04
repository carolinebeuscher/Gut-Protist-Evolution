#!/usr/bin/env python3
"""
Computes per-module statistics comparing functional vs non-functional counts
across grouping variables

For each module x grouping column combination runs:
  - Fisher's exact test (odds ratio + p-value) if there's 2 groups
  - >Chi-squared test (statistic + p-value) of there's >2 groups

A row is written only when at least one group has >=2 samples with non-NaN
functional values. If a module is significant on more than one column,
multiple rows are written.

Usage:
    python3 module_stats.py \\
        --module-files M00001_functionality.tsv M00002_functionality.tsv ... \\
        --pathway-module-ko pathway_module_ko.tsv \\
        --output module_stats.tsv \\
        --alpha 0.05 \\
        --top-kos 5 \\
        --grouping-cols host_association aerotolerance phylum class order

"""

import sys
import os
import argparse
import csv
from collections import defaultdict, Counter
from scipy.stats import fisher_exact, chi2_contingency


# Default grouping columns — any column present in the functionality TSV is valid
DEFAULT_GROUPING_COLS = [
    "host_association",
    "aerotolerance",
    "phylum",
    "class",
    "order",
    "family",
]


### Data loading ##

def load_functionality_tsv(path):
    """
    Read a {module}_functionality.tsv into a list of row dicts.
    Converts the 'functional' column to int (1=True, 0=False, None=unparseable).
    """
    rows = []
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            val = row.get("functional", "").strip()
            if val == "True":
                row["_functional"] = 1
            elif val == "False":
                row["_functional"] = 0
            else:
                row["_functional"] = None   # skip unparseable rows in tests
            rows.append(row)
    return rows


def load_pathway_info(pathway_module_ko):
    """
    Build {module_id: {"pathway_ids": [...], "pathway_names": [...],
                        "module_name": str}} from pathway_module_ko.tsv.
    """
    info = defaultdict(lambda: {"pathway_ids": [], "pathway_names": [],
                                "module_name": "NA"})
    with open(pathway_module_ko) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            mid  = row.get("module_id", "").strip()
            if not mid or mid == "NA":
                continue
            pid  = row.get("pathway_id",   "NA").strip()
            pnm  = row.get("pathway_name", "NA").strip()
            mnm  = row.get("module_name",  "NA").strip()
            if pid not in info[mid]["pathway_ids"]:
                info[mid]["pathway_ids"].append(pid)
            if pnm not in info[mid]["pathway_names"]:
                info[mid]["pathway_names"].append(pnm)
            info[mid]["module_name"] = mnm
    return info


### Statistical tests ###

def fisher_test(group_a, group_b):
    """
    Fisher's exact test on a 2x2 table:
        [[a_true, a_false], [b_true, b_false]]
    Returns (odds_ratio, p_value).
    """
    a_true  = sum(group_a)
    a_false = len(group_a) - a_true
    b_true  = sum(group_b)
    b_false = len(group_b) - b_true

    table = [[a_true, a_false],
             [b_true, b_false]]
    odds_ratio, p_value = fisher_exact(table, alternative="two-sided")
    return odds_ratio, p_value


def chi2_test(groups):
    """
    Chi-squared test across >2 groups.
    groups: {group_label: [0/1, ...]}
    Returns (chi2_stat, p_value), or (None, None) if the contingency table is
    degenerate — e.g. every group is 100% functional or 0% functional, so an
    entire column (true_count or false_count) sums to zero across all groups.
    chi2_contingency can't compute expected frequencies in that case (they'd
    all be zero), and there's no variation to test anyway.
    """
    # Build contingency table: rows = groups, cols = [true_count, false_count]
    table = []
    for vals in groups.values():
        true_count  = sum(vals)
        false_count = len(vals) - true_count
        table.append([true_count, false_count])

    col_totals = [sum(col) for col in zip(*table)]
    if any(total == 0 for total in col_totals):
        return None, None

    chi2, p, dof, expected = chi2_contingency(table)
    return chi2, p


### Missing KO analysis ###

def top_missing_kos(rows, group_col, group_label, n_top):
    """
    For samples in the specified group, count how often each KO appears
    in the missing_kos column. Returns the top n_top as a comma-separated string.
    """
    ko_counter = Counter()
    for row in rows:
        if row.get(group_col, "").strip() != group_label:
            continue
        missing = row.get("missing_kos", "").strip()
        if missing and missing != "NA":
            for ko in missing.split():
                if ko.startswith("K") and ko[1:].isdigit():
                    ko_counter[ko] += 1

    if not ko_counter:
        return "NA"
    top = [ko for ko, _ in ko_counter.most_common(n_top)]
    return ",".join(top)


### Core analysis per module ###

def analyse_module(rows, grouping_cols, alpha, n_top):
    """
    For each grouping column, test whether functional rates differ across groups.
    Returns a list of result dicts (one per significant grouping column).
    """
    results = []

    for col in grouping_cols:
        # Check the column exists in this TSV
        if rows and col not in rows[0]:
            continue

        # Group rows by the value of this column; skip NA/empty
        groups = defaultdict(list)
        for row in rows:
            label = row.get(col, "").strip()
            func  = row["_functional"]
            if label and label != "NA" and func is not None:
                groups[label].append(func)

        # Need at least 2 groups each with >=2 samples to run a test
        valid_groups = {k: v for k, v in groups.items() if len(v) >= 2}
        if len(valid_groups) < 2:
            continue

        group_labels = sorted(valid_groups.keys())

        # Choose test based on number of groups
        if len(valid_groups) == 2:
            a_label, b_label = group_labels
            stat, p = fisher_test(valid_groups[a_label], valid_groups[b_label])
            stat_name      = "odds_ratio"
            groups_compared = "{} vs {}".format(a_label, b_label)
        else:
            stat, p = chi2_test(valid_groups)
            stat_name      = "chi2"
            groups_compared = ",".join(group_labels)

        if p is None:
            # Degenerate table — no variation to test (e.g. every group is
            # 100% or 0% functional)
            continue

        sig = "yes" if p < alpha else "no"

        # Identify the more dysfunctional group (lowest proportion of functional=True)
        rates = {label: sum(vals) / len(vals)
                 for label, vals in valid_groups.items()}
        dysfunctional_group = min(rates, key=rates.get)

        # Most common missing KOs in the dysfunctional group
        missing_kos = top_missing_kos(rows, col, dysfunctional_group, n_top)

        results.append({
            "grouping_col":       col,
            "groups_compared":    groups_compared,
            "sig":                sig,
            "stat_name":          stat_name,
            "stat_value":         round(stat, 6),
            "p_value":            round(p, 6),
            "dysfunctional_group": dysfunctional_group,
            "sig_KOs":            missing_kos,
        })

    return results


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("--module-files",       nargs="+", required=True,
                        help="{module}_functionality.tsv files")
    parser.add_argument("--pathway-module-ko",  required=True,
                        help="pathway_module_ko.tsv from build_kegg_pathway_table.py")
    parser.add_argument("--output",             required=True,
                        help="Output TSV path")
    parser.add_argument("--alpha",              type=float, default=0.05,
                        help="Significance cutoff (default: 0.05)")
    parser.add_argument("--top-kos",            type=int,   default=5,
                        help="Number of top missing KOs to report (default: 5)")
    parser.add_argument("--grouping-cols",      nargs="+",
                        default=DEFAULT_GROUPING_COLS,
                        help="Columns to test (default: {})".format(
                             " ".join(DEFAULT_GROUPING_COLS)))
    args = parser.parse_args()

    print("Loading pathway info...", file=sys.stderr)
    pathway_info = load_pathway_info(args.pathway_module_ko)

    out_cols = [
        "module_id", "module_name", "pathway_ids", "pathway_names",
        "grouping_col", "groups_compared",
        "sig", "stat_name", "stat_value", "p_value",
        "dysfunctional_group", "sig_KOs",
    ]

    rows_written = 0

    with open(args.output, "w") as fh:
        fh.write("\t".join(out_cols) + "\n")

        for module_file in sorted(args.module_files):
            # Derive module_id from filename: M00001_functionality.tsv -> M00001
            module_id = os.path.basename(module_file).replace("_functionality.tsv", "")
            print("  Processing {}...".format(module_id), file=sys.stderr)

            rows = load_functionality_tsv(module_file)
            if not rows:
                print("    WARNING: empty file, skipping", file=sys.stderr)
                continue

            info        = pathway_info.get(module_id, {})
            module_name = info.get("module_name", "NA")
            pathway_ids = ",".join(info.get("pathway_ids",   ["NA"]))
            pathway_nms = ",".join(info.get("pathway_names", ["NA"]))

            results = analyse_module(rows, args.grouping_cols, args.alpha, args.top_kos)

            if not results:
                print("    No testable grouping columns found", file=sys.stderr)
                continue

            for res in results:
                row = {
                    "module_id":           module_id,
                    "module_name":         module_name,
                    "pathway_ids":         pathway_ids,
                    "pathway_names":       pathway_nms,
                    "grouping_col":        res["grouping_col"],
                    "groups_compared":     res["groups_compared"],
                    "sig":                 res["sig"],
                    "stat_name":           res["stat_name"],
                    "stat_value":          res["stat_value"],
                    "p_value":             res["p_value"],
                    "dysfunctional_group": res["dysfunctional_group"],
                    "sig_KOs":             res["sig_KOs"],
                }
                fh.write("\t".join(str(row[c]) for c in out_cols) + "\n")
                rows_written += 1

    sig_count = sum(1 for r in open(args.output)
                    if "\tyes\t" in r)
    print("Done. {} rows written ({} significant at alpha={}).".format(
          rows_written, sig_count, args.alpha), file=sys.stderr)


if __name__ == "__main__":
    main()

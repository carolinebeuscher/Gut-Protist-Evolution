#!/usr/bin/env python3
"""
Joins kofam_summary.tsv with kegg_pathway_hierarchy.tsv and grouping_key.tsv
to produce a per-pathway statistics table comparing HA vs FL and aerobic vs anaerobic.

For each comparison (host_association, aerotolerance) per pathway:
  - Computes min / max / mean percent completion per group
  - Runs Mann-Whitney U test (two-sided)
  - Applies Benjamini-Hochberg FDR correction across all pathways

Usage:
    python3 pathways_summary_stats.py \\
        --kofam-summary   kofam_summary.tsv \\
        --pathway-hier    kegg_pathway_hierarchy.tsv \\
        --grouping-key    grouping_key.tsv \\
        --output          lifestyle_sumstats.tsv \\
        [--alpha          0.05]

"""

import argparse
import sys
from collections import defaultdict

try:
    import pandas as pd
    from scipy.stats import mannwhitneyu
except ImportError:
    sys.exit("ERROR: run: pip install pandas scipy")


def benjamini_hochberg(p_values):
    """BH FDR correction. None values are passed through."""
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


def run_mw(group_a, group_b):
    """
    Mann-Whitney U test between two lists of floats.
    Returns (stat, p) or (None, None) if insufficient data.
    """
    a = [v for v in group_a if v is not None]
    b = [v for v in group_b if v is not None]
    if len(a) < 2 or len(b) < 2:
        return None, None
    if len(set(a + b)) == 1:
        return None, None   # all identical — test undefined
    try:
        stat, p = mannwhitneyu(a, b, alternative="two-sided")
        return round(float(stat), 4), round(float(p), 6)
    except ValueError:
        return None, None


def group_stats(values):
    """Return (min, max, mean) for a list, or (None, None, None) if empty."""
    v = [x for x in values if x is not None]
    if not v:
        return None, None, None
    return round(min(v), 4), round(max(v), 4), round(sum(v) / len(v), 4)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--kofam-summary", required=True,
                        help="kofam_summary.tsv from summarize_kofamscan")
    parser.add_argument("--pathway-hier",  required=True,
                        help="kegg_pathway_hierarchy.tsv (Pathway col in ko00010 format)")
    parser.add_argument("--grouping-key",  required=True,
                        help="grouping_key.tsv from add_clustering rule")
    parser.add_argument("--output",        required=True,
                        help="Output TSV path")
    parser.add_argument("--alpha",         type=float, default=0.05,
                        help="Significance cutoff for BH-corrected p-value (default 0.05)")
    args = parser.parse_args()

    # Load pathway hierarchy 
    print("Loading pathway hierarchy...", file=sys.stderr)
    hier = pd.read_csv(args.pathway_hier, sep="\t", dtype=str)
    # Expect columns: Pathway, Group, Supergroup, Annot
    hier["Pathway"] = hier["Pathway"].str.strip()
    hier_dict = hier.set_index("Pathway").to_dict("index")
    print(f"  {len(hier_dict)} pathways loaded", file=sys.stderr)

    # Load grouping key
    print("Loading grouping key...", file=sys.stderr)
    gk = pd.read_csv(args.grouping_key, sep="\t", dtype=str)

    # Build {sample_id: host_association} and {sample_id: aerotolerance}
    # Also try matching on species name in case sample_id doesn't match column headers
    ha_map   = {}   # sample_id -> "HA" or "FL"
    aero_map = {}   # sample_id -> "aerobic" or "anaerobic"

    for _, row in gk.iterrows():
        sid = str(row.get("sample_id", "")).strip()
        sp  = str(row.get("species",   "")).strip()
        ha  = str(row.get("host_association", "")).strip()
        ae  = str(row.get("aerotolerance",    "")).strip()
        for key in (sid, sp):
            if key and key != "nan":
                ha_map[key]   = ha
                aero_map[key] = ae

    print(f"  {len(ha_map)} sample mappings loaded", file=sys.stderr)

    # Load kofam summary
    print("Loading kofam summary...", file=sys.stderr)
    kofam = pd.read_csv(args.kofam_summary, sep="\t", dtype=str)

    # Identify percent columns: {sample}_percent
    # kofam_summary.py header template: sample + "_percent"
    percent_cols = [c for c in kofam.columns if c.endswith("_percent")]
    sample_names = [c.replace("_percent", "") for c in percent_cols]
    print(f"  {len(kofam)} pathways × {len(percent_cols)} samples", file=sys.stderr)

    # Map sample name → HA/FL and aerobic/anaerobic
    # Try exact match first, then prefix match for cases where kofam
    # truncated the name (kofam_summary.py takes basename split on '_')[0])
    def lookup(name, mapping):
        if name in mapping:
            return mapping[name]
        for key in mapping:
            if key.startswith(name) or name.startswith(key):
                return mapping[key]
        return None

    sample_ha   = {s: lookup(s, ha_map)   for s in sample_names}
    sample_aero = {s: lookup(s, aero_map) for s in sample_names}

    unmatched_ha   = [s for s, v in sample_ha.items()   if v is None]
    unmatched_aero = [s for s, v in sample_aero.items() if v is None]
    if unmatched_ha:
        print(f"  WARNING: {len(unmatched_ha)} samples not matched in "
              f"grouping_key host_association: {unmatched_ha[:5]}",
              file=sys.stderr)

    # Build per-pathway stats 
    print("Computing per-pathway group statistics...", file=sys.stderr)

    # Convert percent columns to float once
    for col in percent_cols:
        kofam[col] = pd.to_numeric(kofam[col], errors="coerce")

    rows = []

    for _, pathway_row in kofam.iterrows():
        pid    = str(pathway_row.get("Pathway", "")).strip()
        annot  = str(pathway_row.get("Annot",   "")).strip()
        orthc  = str(pathway_row.get("Orthcounts", "")).strip()

        # Look up hierarchy info
        hier_info = hier_dict.get(pid, {})
        group      = hier_info.get("Group",      "NA")
        supergroup = hier_info.get("Supergroup", "NA")
        hier_annot = hier_info.get("Annot",      annot)   # prefer hierarchy name

        # Collect per-sample percent values split by group
        fl_vals      = []
        ha_vals      = []
        aerobic_vals = []
        anaerobic_vals = []

        for sample, pcol in zip(sample_names, percent_cols):
            val = pathway_row[pcol]
            if pd.isna(val):
                val = None

            ha_label   = sample_ha.get(sample)
            aero_label = sample_aero.get(sample)

            if ha_label == "FL":
                fl_vals.append(val)
            elif ha_label == "HA":
                ha_vals.append(val)

            if aero_label == "aerobic":
                aerobic_vals.append(val)
            elif aero_label == "anaerobic":
                anaerobic_vals.append(val)

        fl_min,  fl_max,  fl_mean  = group_stats(fl_vals)
        ha_min,  ha_max,  ha_mean  = group_stats(ha_vals)
        ae_min,  ae_max,  ae_mean  = group_stats(aerobic_vals)
        an_min,  an_max,  an_mean  = group_stats(anaerobic_vals)

        ha_stat,   ha_p   = run_mw(fl_vals,      ha_vals)
        aero_stat, aero_p = run_mw(aerobic_vals, anaerobic_vals)

        rows.append({
            "Pathway":    pid,
            "Group":      group,
            "Supergroup": supergroup,
            "Annot":      hier_annot,
            "Orthcounts": orthc,
            # Host association
            "FL_min":    fl_min,   "FL_max":    fl_max,   "FL_mean":    fl_mean,
            "HA_min":    ha_min,   "HA_max":    ha_max,   "HA_mean":    ha_mean,
            "host_association_mw_stat": ha_stat,
            "host_association_p":       ha_p,
            "host_association_p_bh":    None,   # filled below
            "host_association_sig":     None,
            # Aerotolerance
            "aerobic_min":    ae_min,  "aerobic_max":    ae_max,  "aerobic_mean":    ae_mean,
            "anaerobic_min":  an_min,  "anaerobic_max":  an_max,  "anaerobic_mean":  an_mean,
            "aerotolerance_mw_stat": aero_stat,
            "aerotolerance_p":       aero_p,
            "aerotolerance_p_bh":    None,
            "aerotolerance_sig":     None,
        })

    # BH correction 
    ha_bh   = benjamini_hochberg([r["host_association_p"] for r in rows])
    aero_bh = benjamini_hochberg([r["aerotolerance_p"]    for r in rows])

    for r, hbh, abh in zip(rows, ha_bh, aero_bh):
        r["host_association_p_bh"] = hbh
        r["aerotolerance_p_bh"]    = abh
        r["host_association_sig"]  = ("yes" if hbh is not None and hbh < args.alpha
                                      else ("NA" if hbh is None else "no"))
        r["aerotolerance_sig"]     = ("yes" if abh is not None and abh < args.alpha
                                      else ("NA" if abh is None else "no"))

    # Write output 
    out_cols = [
        "Pathway", "Group", "Supergroup", "Annot", "Orthcounts",
        "FL_min",  "FL_max",  "FL_mean",
        "HA_min",  "HA_max",  "HA_mean",
        "host_association_mw_stat", "host_association_p",
        "host_association_p_bh",    "host_association_sig",
        "aerobic_min",   "aerobic_max",   "aerobic_mean",
        "anaerobic_min", "anaerobic_max", "anaerobic_mean",
        "aerotolerance_mw_stat", "aerotolerance_p",
        "aerotolerance_p_bh",    "aerotolerance_sig",
    ]

    result_df = pd.DataFrame(rows, columns=out_cols)
    result_df.to_csv(args.output, sep="\t", index=False)

    n_ha_sig   = (result_df["host_association_sig"] == "yes").sum()
    n_aero_sig = (result_df["aerotolerance_sig"]    == "yes").sum()
    print(f"Done. {len(rows)} pathways written to {args.output}", file=sys.stderr)
    print(f"  Significant (BH<{args.alpha}): "
          f"host_association={n_ha_sig}, aerotolerance={n_aero_sig}",
          file=sys.stderr)


if __name__ == "__main__":
    main()

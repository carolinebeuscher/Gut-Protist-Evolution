#!/usr/bin/env python3
"""
Filters proteomes based on BUSCO completeness using an absolute floor
combined with the Tukey IQR fence.

A sample is removed if its %_complete falls below EITHER:
  - min_complete (an absolute floor, e.g. 10%), or
  - Q1 - 1.5*IQR (the dynamic Tukey fence computed from the batch)

The floor catches categorical failures (e.g. a near-empty proteome from a
failed gene prediction) that the IQR fence alone can miss once the batch
includes genome types with a wide, legitimately-lower expected completeness
range (MAGs/SAGs) — a heterogeneous batch inflates Q1/IQR enough that the
fence can end up well below even genuinely broken samples. The IQR fence
still catches moderate outliers relative to the rest of the batch that a
fixed floor alone wouldn't.

Removed FAA files are moved to a quarantine directory rather than deleted.

Usage:
    python3 filter_busco.py \
        busco_summary.tsv \
        path/to/protein_dir \
        path/to/quarantine_dir \
        kept_samples.txt \
        removed_samples.txt \
        [min_complete]        # default: 10
"""

import sys
import os
import shutil

def main():
    tsv_path       = sys.argv[1]
    protein_dir    = sys.argv[2]
    quarantine_dir = sys.argv[3]
    kept_out       = sys.argv[4]
    removed_out    = sys.argv[5]
    min_complete   = float(sys.argv[6]) if len(sys.argv) > 6 else 10.0

    # Load TSV written by sum_busco.sh
    header = None
    rows   = []
    with open(tsv_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if header is None:
                header = line.split("\t")
                continue
            rows.append(line.split("\t"))

    try:
        name_col     = header.index("inputfile_name")
        complete_col = header.index("%_complete")
    except ValueError as e:
        sys.exit(f"ERROR: expected column not found in {tsv_path}: {e}")

    samples = []
    for row in rows:
        try:
            pct = float(row[complete_col])
        except (ValueError, IndexError):
            print(f"WARNING: could not parse %_complete for {row[name_col]}, skipping",
                  file=sys.stderr)
            continue
        samples.append((row[name_col], pct))

    if not samples:
        sys.exit("ERROR: no valid rows found in BUSCO summary TSV")

    completeness = [pct for _, pct in samples]

    # IQR filter 
    sorted_vals = sorted(completeness)
    n = len(sorted_vals)

    def quartile(data, q):
        """Simple quartile without numpy."""
        pos = (len(data) - 1) * q
        lo, hi = int(pos), int(pos) + 1
        frac = pos - lo
        if hi >= len(data):
            return data[lo]
        return data[lo] + frac * (data[hi] - data[lo])

    q1  = quartile(sorted_vals, 0.25)
    q3  = quartile(sorted_vals, 0.75)
    iqr = q3 - q1
    lower_fence = q1 - 1.5 * iqr

    print(f"Samples      : {n}")
    print(f"Q1           : {q1:.2f}%")
    print(f"Q3           : {q3:.2f}%")
    print(f"IQR          : {iqr:.2f}")
    print(f"IQR fence    : {lower_fence:.2f}%  (Q1 - 1.5*IQR)")
    print(f"Absolute floor: {min_complete:.2f}%")
    print(f"Effective cutoff: {max(lower_fence, min_complete):.2f}%  "
          f"(whichever is higher removes more)")
    print()

    # Quarantine low-completeness FAAs 
    os.makedirs(quarantine_dir, exist_ok=True)

    kept, removed = [], []

    for filename, pct in samples:
        # filename is the short_summary filename; extract the sample name from it
        # e.g. "short_summary.specific.mysample.txt" → "mysample"
        sample_name = filename.replace("short_summary.specific.", "").replace(".txt", "")
        faa_src  = os.path.join(protein_dir,    f"{sample_name}.faa")
        faa_dest = os.path.join(quarantine_dir, f"{sample_name}.faa")

        below_floor = pct < min_complete
        below_fence = pct < lower_fence
        if below_floor or below_fence:
            reason = ("below absolute floor" if below_floor and not below_fence
                       else "below IQR fence" if below_fence and not below_floor
                       else "below both floor and IQR fence")
            if os.path.exists(faa_src):
                shutil.move(faa_src, faa_dest)
                print(f"REMOVED  {sample_name:40s}  {pct:.2f}%  ({reason})")
            else:
                print(f"WARNING: {faa_src} not found — skipping move", file=sys.stderr)
            removed.append(sample_name)
        else:
            print(f"KEPT     {sample_name:40s}  {pct:.2f}%")
            kept.append(sample_name)

    #Write kept / removed lists 
    with open(kept_out, "w") as fh:
        fh.write("\n".join(kept) + "\n")
    with open(removed_out, "w") as fh:
        fh.write("\n".join(removed) + "\n")

    print(f"\nKept {len(kept)} samples, removed {len(removed)} samples.")
    if removed:
        print(f"Removed FAAs moved to: {quarantine_dir}")

if __name__ == "__main__":
    main()

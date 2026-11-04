#!/usr/bin/env python3
"""
Combines the existing grouping_key.tsv with the output of parse_genome_xlsx.py
(new genomes) into one final grouping_key.tsv, keeping every column from both
(including the new genome_type column).

Unlike the original version of this script, this does NOT route through
add_clustering.py — that would drop sample_id and re-derive it from
sample_key.tsv's per-species lookup, which collapses every genome of the same
species onto one (comma-joined) sample_id. Since generate_sample_ids.py +
rename_xlsx_genomes.py already gave every new genome its own correct,
per-row sample_id, this script keeps that as-is and instead computes
host_association/aerotolerance directly — the same HA_list/aerobic_list
species lookup add_clustering.py uses.

If a new genome's species already has a row in the old grouping_key.tsv, the
new row is DROPPED and the old row is kept — the old, curated entry wins.
(Pair this with filter_representative_genomes.py upstream, which already
narrows the new genomes down to one per species, so at most one new row can
ever be dropped per species.)

Old rows get a constant genome_type (see [old_genome_type]) since the old
grouping_key.tsv predates that column.

Output columns:
    sample_id, species, genus, family, order, class, phylum, domain,
    clade_1, clade_2, clade_3, host_association, aerotolerance, genome_type

Usage:
    python3 combine_taxonomy.py \\
        input_files/grouping_key.tsv \\
        new_genomes_taxonomy.tsv \\
        input_files/HA.txt \\
        input_files/aerobic.txt \\
        combined_grouping_key.tsv \\
        [old_genome_type]        # default: monoculture
"""

import sys

OLD_HEADER = ["sample_id", "species", "genus", "family", "order", "class",
              "phylum", "domain", "clade_1", "clade_2", "clade_3",
              "host_association", "aerotolerance"]

NEW_HEADER = ["sample_id", "species", "genus", "family", "order", "class",
              "phylum", "domain", "clade_1", "clade_2", "clade_3",
              "genome_type"]

OUT_HEADER = ["sample_id", "species", "genus", "family", "order", "class",
              "phylum", "domain", "clade_1", "clade_2", "clade_3",
              "host_association", "aerotolerance", "genome_type"]


def load_set(path):
    """Return a set of stripped, non-empty lines from a file."""
    with open(path) as fh:
        return {line.strip() for line in fh if line.strip()}


def read_tsv_rows(path, expected_header):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        if header != expected_header:
            sys.exit(f"ERROR: unexpected header in {path}.\n"
                      f"  expected: {expected_header}\n"
                      f"  found:    {header}")
        rows = []
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            rows.append(line.split("\t"))
    return rows


def main():
    if len(sys.argv) < 6:
        sys.exit("Usage: combine_taxonomy.py <old_grouping_key.tsv> "
                  "<new_taxonomy.tsv> <HA_list.txt> <aerobic_list.txt> "
                  "<output.tsv> [old_genome_type]")

    old_path        = sys.argv[1]
    new_path        = sys.argv[2]
    ha_path         = sys.argv[3]
    aerobic_path    = sys.argv[4]
    output_path     = sys.argv[5]
    old_genome_type = sys.argv[6] if len(sys.argv) > 6 else "monoculture"

    old_rows = read_tsv_rows(old_path, OLD_HEADER)
    new_rows = read_tsv_rows(new_path, NEW_HEADER)
    ha_species      = load_set(ha_path)
    aerobic_species = load_set(aerobic_path)

    # old: sample_id..aerotolerance kept as-is; append constant genome_type
    old_out = [row + [old_genome_type] for row in old_rows]
    old_species = {row[1] for row in old_out}

    # new: sample_id..clade_3 kept as-is (cols 0:11), genome_type at 11;
    # host_association/aerotolerance computed fresh via species lookup,
    # same rule as add_clustering.py: absence from the list means FL/anaerobic.
    # A species already present in the old grouping_key.tsv is dropped here —
    # the existing curated row wins over the new genome for that species.
    new_out = []
    dropped = []
    unmatched_ha, unmatched_aero = [], []
    for row in new_rows:
        species = row[1]
        if species in old_species:
            dropped.append((row[0], species))
            continue
        genome_type = row[11]
        host_association = "HA" if species in ha_species else "FL"
        aerotolerance = "aerobic" if species in aerobic_species else "anaerobic"
        if species not in ha_species:
            unmatched_ha.append(species)
        if species not in aerobic_species:
            unmatched_aero.append(species)
        new_out.append(row[0:11] + [host_association, aerotolerance, genome_type])

    # Sanity check: sample_id must be unique across the combined table —
    # a collision here would silently merge two distinct genomes downstream.
    all_ids = [row[0] for row in old_out] + [row[0] for row in new_out]
    dupes = sorted({sid for sid in all_ids if all_ids.count(sid) > 1})
    if dupes:
        sys.exit(f"ERROR: sample_id collision(s) between old and new rows: "
                  f"{dupes}")

    with open(output_path, "w") as fh:
        fh.write("\t".join(OUT_HEADER) + "\n")
        for row in old_out + new_out:
            fh.write("\t".join(row) + "\n")

    print(f"Wrote {len(old_out)} old + {len(new_out)} new = "
          f"{len(old_out) + len(new_out)} rows to {output_path}",
          file=sys.stderr)
    print(f"  old rows tagged genome_type='{old_genome_type}'", file=sys.stderr)
    print(f"  new rows: {len(new_out) - len(unmatched_ha)} HA / "
          f"{len(unmatched_ha)} FL (not in {ha_path})", file=sys.stderr)
    print(f"  new rows: {len(new_out) - len(unmatched_aero)} aerobic / "
          f"{len(unmatched_aero)} anaerobic (not in {aerobic_path})",
          file=sys.stderr)
    if dropped:
        print(f"\nDropped {len(dropped)} new genome(s) — species already in "
              f"the old grouping_key.tsv, old row kept:", file=sys.stderr)
        for sample_id, species in dropped:
            print(f"  {sample_id}\t{species}", file=sys.stderr)


if __name__ == "__main__":
    main()

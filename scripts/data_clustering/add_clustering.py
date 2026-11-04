#!/usr/bin/env python3
"""
Adds host_association, aerotolerance, and sample_id columns to a taxonomy TSV using:
  - taxonomy.tsv      : tsv output from get_taxonomy
  - HA_list.txt       : list of host-associated species (one species per line)
        - if a species isn't in here, labeled as FL
  - aerobic_list.txt  : list of aerobic species (one species per line)
        - if a species isn't in here, labeled as anaerobic
  - sample_key.tsv    : maps basenames to species name (sample_id<TAB>species)
        - if a species isn't in here, labeled as NA

Usage:
    python3 add_clustering.py \\
        taxonomy_raw.tsv \\
        HA_list.txt \\
        aerobic_list.txt \\
        sample_key.tsv \\
        grouping_key.tsv
"""

import sys


def load_set(path):
    """Return a set of stripped, non-empty lines from a file."""
    with open(path) as fh:
        return {line.strip() for line in fh if line.strip()}


def load_sample_key(path):
    """
    Load a two-column TSV (sample_id<TAB>species) into {species: sample_id}.
    If a species maps to multiple sample IDs, they are joined with a comma.
    """
    species_to_sample = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                print(f"WARNING: skipping malformed line in sample key: {line!r}",
                      file=sys.stderr)
                continue
            sample_id = parts[0].strip()
            species   = parts[1].strip()
            if species in species_to_sample:
                species_to_sample[species] = species_to_sample[species] + "," + sample_id
            else:
                species_to_sample[species] = sample_id
    return species_to_sample


def main():
    taxonomy_path    = sys.argv[1]
    ha_path          = sys.argv[2]
    aerobic_path     = sys.argv[3]
    sample_key_path  = sys.argv[4]
    output_path      = sys.argv[5]

    ha_species       = load_set(ha_path)
    aerobic_species  = load_set(aerobic_path)
    species_to_sample = load_sample_key(sample_key_path)

    unmatched = []

    with open(taxonomy_path) as fin, open(output_path, "w") as fout:

        # Write header with the three new columns prepended for easy joining
        header = fin.readline().rstrip("\n")
        fout.write("sample_id\t" + header + "\thost_association\taerotolerance\n")

        for line in fin:
            line = line.rstrip("\n")
            if not line:
                continue

            species = line.split("\t")[0]

            host_association = "HA"      if species in ha_species     else "FL"
            aerotolerance    = "aerobic" if species in aerobic_species else "anaerobic"
            sample_id        = species_to_sample.get(species, "NA")

            if sample_id == "NA":
                unmatched.append(species)

            fout.write(f"{sample_id}\t{line}\t{host_association}\t{aerotolerance}\n")

    if unmatched:
        print(f"WARNING: {len(unmatched)} species had no matching sample_id "
              f"in {sample_key_path}:", file=sys.stderr)
        for sp in unmatched:
            print(f"  {sp}", file=sys.stderr)

    print(f"Wrote grouping key to {output_path}", file=sys.stderr)


if __name__ == "__main__":
    main()

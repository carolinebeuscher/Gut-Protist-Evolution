#!/usr/bin/env python3
"""
Copies each new genome's .faa file into a directory, renamed to its short
sample_id — and by construction only the representative genomes, since it
only processes what's listed in the rename map. Pass it the rename map
produced by running generate_sample_ids.py on the REPRESENTATIVE-filtered
xlsx (i.e. after filter_representative_genomes.py), not the full 450-genome
map — that's what makes this a filter as well as a rename.

Before copying, checks whether a .faa for that SAME SPECIES already exists in
the output directory — by scanning every "{code}.faa" already there and
looking up its species via old_sample_key.tsv (existing genomes) or
new_sample_key_rows.tsv (genomes this script already copied, e.g. from a
previous run). If a match is found, the new file is skipped and whatever is
already in the directory is kept as-is — this is what makes duplicates like
Encau.faa + Encau-1.faa impossible: only the first .faa for a species ever
lands in the output directory, nothing gets overwritten.

Source .faa files are looked up by the rename map's original_name first, and
by that name with "_" replaced by "-" as a fallback (in case the file on disk
was already renamed for hyphen-only naming but the xlsx wasn't).

Copies rather than moves/renames in place, so nothing in the source directory
is touched. This script also never deletes anything — if you already have
duplicate .faa files in the output directory from before, this won't clean
those up; it only prevents new ones going forward.

Usage:
    python3 rename_and_filter_faa.py \\
        new_genome_rename_map.tsv \\
        new_sample_key_rows.tsv \\
        input_files/sample_key.tsv \\
        path/to/source_faa_dir \\
        path/to/output_faa_dir
"""

import os
import shutil
import sys


def read_rename_map(path):
    mapping = []
    with open(path) as fh:
        fh.readline()  # header: original_name\tnew_sample_id
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            original, new_id = line.split("\t")
            mapping.append((original, new_id))
    return mapping


def read_code_to_species(path):
    code_to_species = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            code_to_species[parts[0].strip()] = parts[1].strip()
    return code_to_species


def find_source_faa(source_dir, original_name):
    for candidate in (original_name, original_name.replace("_", "-")):
        path = os.path.join(source_dir, candidate + ".faa")
        if os.path.isfile(path):
            return path
    return None


def species_already_in_output(output_dir, code_to_species):
    """
    Scan existing "{code}.faa" files in output_dir and return the set of
    species already represented there, using code_to_species (old + new
    combined) to translate filenames to species.
    """
    present = set()
    unrecognized = []
    if not os.path.isdir(output_dir):
        return present, unrecognized
    for fname in os.listdir(output_dir):
        if not fname.endswith(".faa"):
            continue
        code = fname[:-len(".faa")]
        species = code_to_species.get(code)
        if species:
            present.add(species)
        else:
            unrecognized.append(fname)
    return present, unrecognized


def main():
    if len(sys.argv) != 6:
        sys.exit("Usage: rename_and_filter_faa.py <rename_map.tsv> "
                  "<new_sample_key_rows.tsv> <old_sample_key.tsv> "
                  "<source_faa_dir> <output_faa_dir>")

    map_path          = sys.argv[1]
    new_sample_key    = sys.argv[2]
    old_sample_key    = sys.argv[3]
    source_dir        = sys.argv[4]
    output_dir        = sys.argv[5]

    mapping = read_rename_map(map_path)
    new_code_to_species = read_code_to_species(new_sample_key)
    old_code_to_species = read_code_to_species(old_sample_key)
    code_to_species = {**old_code_to_species, **new_code_to_species}

    os.makedirs(output_dir, exist_ok=True)
    species_present, unrecognized = species_already_in_output(output_dir, code_to_species)

    copied, missing, skipped_dupe, skipped_no_species = [], [], [], []
    for original_name, new_id in mapping:
        species = new_code_to_species.get(new_id)
        if species is None:
            skipped_no_species.append(new_id)
            continue
        if species in species_present:
            skipped_dupe.append((original_name, new_id, species))
            continue

        src = find_source_faa(source_dir, original_name)
        if src is None:
            missing.append(original_name)
            continue

        dst = os.path.join(output_dir, f"{new_id}.faa")
        shutil.copy2(src, dst)
        species_present.add(species)  # so later entries for the same species also skip
        copied.append((original_name, new_id))

    print(f"Copied {len(copied)}/{len(mapping)} .faa files -> {output_dir}",
          file=sys.stderr)

    if skipped_dupe:
        print(f"\nSkipped {len(skipped_dupe)} genome(s) — a .faa for that "
              f"species was already in {output_dir}:", file=sys.stderr)
        for original_name, new_id, species in skipped_dupe:
            print(f"  {new_id} ({species}) — kept existing file", file=sys.stderr)

    if missing:
        print(f"\nWARNING: {len(missing)} genome(s) had no matching .faa file "
              f"in {source_dir} (tried both the original name and its "
              f"underscore-to-hyphen form):", file=sys.stderr)
        for name in missing:
            print(f"  {name}", file=sys.stderr)

    if skipped_no_species:
        print(f"\nWARNING: {len(skipped_no_species)} sample_id(s) from "
              f"{map_path} were not found in {new_sample_key} — skipped "
              f"since their species is unknown: {skipped_no_species}",
              file=sys.stderr)

    if unrecognized:
        print(f"\nNOTE: {len(unrecognized)} .faa file(s) already in "
              f"{output_dir} don't match any code in {old_sample_key} or "
              f"{new_sample_key}, so they weren't counted as covering a "
              f"species: {unrecognized[:10]}"
              f"{' ...' if len(unrecognized) > 10 else ''}", file=sys.stderr)


if __name__ == "__main__":
    main()

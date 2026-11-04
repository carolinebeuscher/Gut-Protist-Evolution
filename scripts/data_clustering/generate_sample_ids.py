#!/usr/bin/env python3
"""
Generates short sample_key.tsv-style codes for new genomes, following the
existing convention (verified against input_files/sample_key.tsv):
    genus[:2] + species_epithet[:3]           e.g. Dasytricha ruminantium -> Darum
For informal "Genus sp. TAG" species:
    genus[:2] + TAG (verbatim)                e.g. Dasytricha sp. RCG4    -> DaRCG4
When multiple new genomes share the exact same species (or the computed code
already exists in sample_key.tsv for that same species), every one of them is
suffixed -1, -2, ... in sample_id order — none are left unsuffixed once there's
more than one, so an existing bare code in sample_key.tsv is never collided with.

When two DIFFERENT species compute to the same base code (e.g. "Entodinium
simplex" and "Entodinium simulans" both reducing to "Ensim"), the later one
(by sort order) gets an "x" appended to its base — "Ensimx" — before any -1/-2
suffixing is applied, so the two species' codes never overlap. "x" is repeated
("xx", "xxx", ...) in the rare case of a three-or-more-way collision.

This is step 1 of the new-genome workflow — it runs directly off the raw
summary xlsx (before parse_genome_xlsx.py), so its rename map can be applied
to the xlsx's "Genome name" column (see rename_xlsx_genomes.py) before
parse_genome_xlsx.py ever runs. That way every later step (parse_genome_xlsx.py,
combine_taxonomy.py, and eventually add_clustering.py via sample_key.tsv)
already sees the short codes as the sample name.

Does NOT modify sample_key.tsv — writes the new rows to a separate file you
review and append yourself (`cat <output> >> sample_key.tsv`).

Usage:
    python3 generate_sample_ids.py \\
        RCG_450_summary.xlsx \\
        input_files/sample_key.tsv \\
        new_sample_key_rows.tsv \\
        new_genome_rename_map.tsv

Install dependencies:
    pip install pandas openpyxl
"""

import re
import sys

try:
    import pandas as pd
except ImportError:
    sys.exit("ERROR: run: pip install pandas openpyxl")

from collections import defaultdict

SP_PATTERN = re.compile(r"\bsp\.\s*(\S+)")

GTDB_RANK_PREFIXES = {
    "d": "domain", "p": "phylum", "c": "class", "o": "order",
    "f": "family", "g": "genus", "s": "species",
}


def compute_base_code(genus, species):
    m = SP_PATTERN.search(species)
    if m:
        tag = m.group(1)
        return genus[:2] + tag
    parts = species.split()
    epithet = parts[1] if len(parts) > 1 else parts[0]
    return genus[:2] + epithet[:3]


def parse_gtdb_taxonomy(tax_string, row_label):
    """
    Parse "d__X;p__X;c__X;o__X;f__X;g__X;s__X" into a rank dict.
    Mirrors parse_genome_xlsx.py's parser — kept duplicated rather than
    imported, matching this codebase's convention of standalone scripts.
    """
    ranks = {}
    for token in tax_string.split(";"):
        token = token.strip()
        if "__" not in token:
            sys.exit(f"ERROR: unparseable taxonomy token '{token}' for "
                      f"sample '{row_label}' (full string: {tax_string!r})")
        prefix, _, value = token.partition("__")
        rank_name = GTDB_RANK_PREFIXES.get(prefix)
        if rank_name is None:
            sys.exit(f"ERROR: unrecognized rank prefix '{prefix}__' for "
                      f"sample '{row_label}' (full string: {tax_string!r})")
        ranks[rank_name] = value.strip() or "NA"

    missing = [r for r in ("genus", "species") if r not in ranks]
    if missing:
        sys.exit(f"ERROR: taxonomy string for sample '{row_label}' is "
                  f"missing rank(s) {missing} (full string: {tax_string!r})")
    return ranks


def read_xlsx_genomes(xlsx_path):
    df = pd.read_excel(xlsx_path, sheet_name=0, engine="openpyxl")

    required_cols = ["Genome name", "Taxonomy"]
    missing_cols = [c for c in required_cols if c not in df.columns]
    if missing_cols:
        sys.exit(f"ERROR: missing expected column(s) {missing_cols} in "
                  f"{xlsx_path}. Found columns: {list(df.columns)}")

    rows = []
    seen = set()
    for _, row in df.iterrows():
        raw_name = str(row["Genome name"]).strip()
        if not raw_name or raw_name.lower() == "nan":
            sys.exit(f"ERROR: row with missing 'Genome name' in {xlsx_path}")
        if raw_name in seen:
            sys.exit(f"ERROR: duplicate 'Genome name' value '{raw_name}' in "
                      f"{xlsx_path}")
        seen.add(raw_name)

        tax_string = str(row["Taxonomy"]).strip()
        ranks = parse_gtdb_taxonomy(tax_string, raw_name)

        rows.append({
            "sample_id": raw_name,
            "species":   ranks["species"],
            "genus":     ranks["genus"],
        })
    return rows


def read_old_sample_key(path):
    """Returns (used_codes set, {code: species})."""
    used_codes = set()
    code_to_species = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            code, species = parts[0].strip(), parts[1].strip()
            used_codes.add(code)
            code_to_species[code] = species
    return used_codes, code_to_species


def main():
    if len(sys.argv) < 5:
        sys.exit("Usage: generate_sample_ids.py <RCG_450_summary.xlsx> "
                  "<old_sample_key.tsv> <new_sample_key_rows.tsv> "
                  "<new_genome_rename_map.tsv>")

    xlsx_path = sys.argv[1]
    old_sample_key_path = sys.argv[2]
    out_rows_path = sys.argv[3]
    out_rename_path = sys.argv[4]

    new_rows = read_xlsx_genomes(xlsx_path)
    used_codes, code_to_species = read_old_sample_key(old_sample_key_path)

    # Group new genomes by species so duplicates get consistent suffixing
    by_species = defaultdict(list)
    for row in new_rows:
        by_species[row["species"]].append(row)

    assigned = []       # (sample_id, code, species)
    disambiguated = []  # (species, original_base, final_base) — for the summary
    unresolved = []     # (species, base_code, conflicting_species, sample_ids) — should stay empty

    base_owner = {}  # base_code (pre-suffix) -> species that claimed it

    for species in sorted(by_species):
        members = sorted(by_species[species], key=lambda r: r["sample_id"])
        genus = members[0]["genus"]
        base = compute_base_code(genus, species)

        # Disambiguate the base itself against anything already claiming it
        # (an old exact code, or a base already claimed by a different new
        # species this run) before any -1/-2 suffixing is applied.
        final_base = base
        while (code_to_species.get(final_base) not in (None, species)
               or base_owner.get(final_base) not in (None, species)):
            final_base += "x"
        if final_base != base:
            disambiguated.append((species, base, final_base))
        base_owner[final_base] = species

        need_suffix = len(members) > 1 or final_base in used_codes
        candidates = ([final_base] if not need_suffix
                      else [f"{final_base}-{i}" for i in range(1, len(members) + 1)])

        # Defensive final check — should be unreachable given the base
        # disambiguation above, but never silently overwrite a collision.
        conflicts = sorted({code_to_species[c] for c in candidates
                             if c in code_to_species and code_to_species[c] != species})
        if conflicts:
            unresolved.append((species, final_base, conflicts,
                                [m["sample_id"] for m in members]))
            continue

        for member, code in zip(members, candidates):
            used_codes.add(code)
            code_to_species[code] = species
            assigned.append((member["sample_id"], code, species))

    with open(out_rows_path, "w") as fh:
        for sample_id, code, species in assigned:
            fh.write(f"{code}\t{species}\n")

    with open(out_rename_path, "w") as fh:
        fh.write("original_name\tnew_sample_id\n")
        for sample_id, code, species in assigned:
            fh.write(f"{sample_id}\t{code}\n")

    print(f"Assigned {len(assigned)} codes -> {out_rows_path}", file=sys.stderr)
    print(f"Wrote rename map -> {out_rename_path}", file=sys.stderr)
    print(f"  (nothing was modified in {old_sample_key_path} or any genome "
          f"files — review, then: cat {out_rows_path} >> {old_sample_key_path})",
          file=sys.stderr)

    if disambiguated:
        print(f"\n{len(disambiguated)} species had a base code collision with "
              f"a different species, resolved with 'x':", file=sys.stderr)
        for species, base, final_base in disambiguated:
            print(f"  '{species}': '{base}' -> '{final_base}'", file=sys.stderr)

    if unresolved:
        print(f"\nWARNING: {len(unresolved)} species left UNASSIGNED despite "
              f"'x' disambiguation — this shouldn't normally happen; "
              f"investigate manually:", file=sys.stderr)
        for species, base, conflicts, sample_ids in unresolved:
            print(f"  '{species}' (code '{base}') collides with "
                  f"{conflicts} — affected genomes: {sample_ids}",
                  file=sys.stderr)


if __name__ == "__main__":
    main()

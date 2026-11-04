#!/usr/bin/env python3
"""
Builds a TSV of KEGG pathway -> module -> KO relationships.

Output columns:
    pathway_id, pathway_name, module_id, module_name, ko_id

Data sources:
    --pathways   local file (ko_id<TAB>pathway_name, one pathway per line)
    --ko-pathway local file (pathway_id<TAB>ko_id, one KO per line)
    module names and KO->module links fetched from KEGG REST API

The KO->module join is done locally: for every KO already in your
pathway->KO file, we look up which module(s) it belongs to.
Rows without a module assignment are still written with module_id = NA.
** this is okay because some kos are used in a pathway but not in any of the modules 

Usage:
    python3 build_kegg_pathway_table.py \\
        --pathways       kegg_pathways.txt \\
        --ko-pathway     kegg_orthology_pathways.txt \\
        --output         pathway_module_ko.tsv

Install dependency:
    pip install requests
"""

import sys
import argparse
import time
import requests
from collections import defaultdict


BASE = "https://rest.kegg.jp"


def kegg_get(endpoint, retries=3, wait=5):
    url = "{}/{}".format(BASE, endpoint)
    for attempt in range(1, retries + 1):
        try:
            r = requests.get(url, timeout=60)
            r.raise_for_status()
            return r.text
        except requests.RequestException as e:
            print("  Attempt {}/{} failed for {}: {}".format(
                  attempt, retries, url, e), file=sys.stderr)
            if attempt < retries:
                time.sleep(wait)
    sys.exit("ERROR: could not fetch {} after {} attempts".format(url, retries))


def load_tsv_dict(path, key_col=0, val_col=1):
    """Load a two-column TSV into a dict {key: name}."""
    d = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) > max(key_col, val_col):
                d[parts[key_col].strip()] = parts[val_col].strip()
    return d


def load_tsv_multidict(path, key_col=0, val_col=1):
    """Load a two-column TSV into a dict {key: [val, ...]}."""
    d = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) > max(key_col, val_col):
                d[parts[key_col].strip()].append(parts[val_col].strip())
    return d


def parse_link(text):
    """Parse a KEGG link response into {source: [target, ...]}."""
    mapping = defaultdict(list)
    for line in text.strip().split("\n"):
        parts = line.split("\t")
        if len(parts) == 2:
            mapping[parts[0].strip()].append(parts[1].strip())
    return mapping


def parse_list(text):
    """Parse a KEGG list response into {id: name}."""
    names = {}
    for line in text.strip().split("\n"):
        parts = line.split("\t", 1)
        if len(parts) == 2:
            names[parts[0].strip()] = parts[1].strip()
    return names


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--pathways",   required=True,
                        help="kegg_pathways.txt: pathway_id<TAB>pathway_name")
    parser.add_argument("--ko-pathway", required=True,
                        help="kegg_orthology_pathways.txt: pathway_id<TAB>ko_id")
    parser.add_argument("--output",     required=True,
                        help="Output TSV path")
    args = parser.parse_args()

    # 1. Load local files (reliable)
    print("Loading pathway names from {}...".format(args.pathways), file=sys.stderr)
    pathway_names = load_tsv_dict(args.pathways, key_col=0, val_col=1)
    print("  {} pathways loaded".format(len(pathway_names)), file=sys.stderr)

    print("Loading pathway->KO from {}...".format(args.ko_pathway), file=sys.stderr)
    pathway_to_kos = load_tsv_multidict(args.ko_pathway, key_col=0, val_col=1)
    print("  {} pathway->KO pairs loaded".format(
          sum(len(v) for v in pathway_to_kos.values())), file=sys.stderr)

    # Collect all unique KOs so we only look up modules for KOs we actually have
    all_kos = {ko for kos in pathway_to_kos.values() for ko in kos}
    print("  {} unique KOs".format(len(all_kos)), file=sys.stderr)

    # 2. Fetch module names and KO->module links from API
    print("Fetching module names from KEGG API...", file=sys.stderr)
    module_names_raw = parse_list(kegg_get("list/module"))
    # Normalise: strip "md:" prefix for consistent keys
    module_names = {
        mid.replace("md:", ""): name
        for mid, name in module_names_raw.items()
    }
    print("  {} module names fetched".format(len(module_names)), file=sys.stderr)

    print("Fetching KO -> module links from KEGG API...", file=sys.stderr)
    ko_to_modules_raw = parse_link(kegg_get("link/module/ko"))
    # Normalise: strip "ko:" from keys and "md:" from values
    ko_to_modules = {
        k.replace("ko:", ""): [v.replace("md:", "") for v in vals]
        for k, vals in ko_to_modules_raw.items()
    }
    print("  {} KO->module links fetched".format(len(ko_to_modules)), file=sys.stderr)

    # 3. Write output TSV
    print("Writing to {}...".format(args.output), file=sys.stderr)

    rows_written   = 0
    no_module_kos  = 0

    with open(args.output, "w") as fh:
        fh.write("pathway_id\tpathway_name\tmodule_id\tmodule_name\tko_id\n")

        for pathway_id in sorted(pathway_to_kos.keys()):
            pathway_name = pathway_names.get(pathway_id, "NA")
            kos = pathway_to_kos[pathway_id]

            for ko in sorted(kos):
                modules = ko_to_modules.get(ko, [])

                if not modules:
                    # Still write the row; module columns are NA
                    fh.write("{}\t{}\tNA\tNA\t{}\n".format(
                        pathway_id, pathway_name, ko))
                    no_module_kos += 1
                    rows_written += 1
                else:
                    for mod_id in sorted(modules):
                        module_name = module_names.get(mod_id, "NA")
                        fh.write("{}\t{}\t{}\t{}\t{}\n".format(
                            pathway_id, pathway_name,
                            mod_id, module_name,
                            ko))
                        rows_written += 1

    print("Done. {} rows written.".format(rows_written), file=sys.stderr)
    if no_module_kos:
        print("  {} KOs had no module assignment (written with module=NA).".format(
              no_module_kos), file=sys.stderr)


if __name__ == "__main__":
    main()

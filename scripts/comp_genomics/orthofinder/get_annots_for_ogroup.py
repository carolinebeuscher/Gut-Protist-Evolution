#! /usr/bin/env python3
"""
Gets annotations for a given orthogroup from per-sample InterProScan output
({sample}_interprosum.tsv files from clean_interproscan), one row per gene.

Usage:
    python3 get_annots_for_ogroup.py \\
        --ogroup    OG0000170 \\
        --orthogroups Orthogroups.tsv \\
        --interpro-outdir path/to/gene_annotation/interproscan
"""

import sys
import os
import argparse


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--ogroup",          required=True,
                        help="Orthogroup ID to look up (e.g. OG0000170)")
    parser.add_argument("--orthogroups",     required=True,
                        help="Orthogroups.tsv from OrthoFinder")
    parser.add_argument("--interpro-outdir", required=True,
                        help="InterProScan outdir; expects "
                             "{interpro_outdir}/{sample}/{sample}_interprosum.tsv "
                             "(clean_interproscan's output)")
    return parser.parse_args()


def main():
    args = parse_args()

    with open(args.orthogroups) as f:
        sns = f.readline().strip("\n").split("\t")[1:]

        ogroup_genes = None
        for line in f:
            line = line.strip("\n")
            o = line.split("\t")[0]
            if o == args.ogroup:
                ogroup_genes = line.split("\t")[1:]
                break

    if ogroup_genes is None:
        print(f"WARNING: orthogroup '{args.ogroup}' not found in {args.orthogroups}",
              file=sys.stderr)
        return

    print("Ogroup\tSample\tGene\tCDD-IPR\tCoils\tFunFam\tGene3D\tMobiDBLite\tNCBIfam\t"
          "PANTHER\tPfam\tPIRSF\tPRINTS\tSFLD\tSMART\tSUPERFAMILY\tInterpro\tGO")

    for i, sample in enumerate(sns):
        genes = [e.strip() for e in ogroup_genes[i].split(",") if e.strip()]
        if not genes:
            continue

        interpro_path = os.path.join(args.interpro_outdir, sample,
                                      f"{sample}_interprosum.tsv")
        if not os.path.exists(interpro_path):
            print(f"  WARNING: no InterProScan summary for '{sample}' at "
                  f"{interpro_path}", file=sys.stderr)
            continue

        with open(interpro_path) as ef:
            ef.readline()  # skip header
            for line in ef:
                line = line.strip("\n")
                if not line:
                    continue
                gene = line.split("\t")[0]
                if gene in genes:
                    print(f"{args.ogroup}\t{sample}\t{line}")


if __name__ == "__main__":
    main()

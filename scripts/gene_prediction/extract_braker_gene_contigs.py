#!/usr/bin/env python3
"""
Extracts a gene_id -> contig mapping from BRAKER3's augustus.hints.gtf files
(one per sample)

Reads each GTF's transcript feature lines specifically 
    <contig>\tAUGUSTUS\ttranscript\t<start>\t<end>\t.\t<strand>\t.\t<transcript_id>

Usage:
    python3 extract_braker_gene_contigs.py \\
        --gtf-list  gtf_paths.txt \\
        --output    gene_contig_map.tsv

--gtf-list should have one path per line. Each GTF's sample name is inferred
from its path as the directory N levels up
"""

import argparse
import os
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--gtf-list", required=True,
                        help="Text file with paths to augustus.hints.gtf files "
                             "(one per line, one per sample)")
    parser.add_argument("--output", required=True,
                        help="Output TSV: sample\\tgene_id\\tcontig")
    parser.add_argument("--sample-dir-level", type=int, default=3,
                        help="How many directories up from the gtf file the "
                             "sample-named directory is (default: 3, matching "
                             "{gp_outdir}/{sample}/braker/Augustus/augustus.hints.gtf)")
    return parser.parse_args()


def sample_from_path(gtf_path, level):
    path = os.path.abspath(gtf_path)
    for _ in range(level):
        path = os.path.dirname(path)
    return os.path.basename(path)


def extract_transcript_contigs(gtf_path):
    """Returns {transcript_id: contig} from a GTF's 'transcript' lines."""
    mapping = {}
    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "transcript":
                continue
            contig = fields[0]
            transcript_id = fields[8].strip()
            mapping[transcript_id] = contig
    return mapping


def main():
    args = parse_args()

    with open(args.gtf_list) as f:
        gtf_paths = [l.strip() for l in f if l.strip()]

    with open(args.output, "w") as out:
        out.write("sample\tgene_id\tcontig\n")
        total = 0
        for gtf_path in gtf_paths:
            if not os.path.exists(gtf_path):
                print(f"WARNING: {gtf_path} not found — skipping", file=sys.stderr)
                continue
            sample = sample_from_path(gtf_path, args.sample_dir_level)
            mapping = extract_transcript_contigs(gtf_path)
            print(f"  {sample}: {len(mapping)} transcript(s) from {gtf_path}",
                  file=sys.stderr)
            for gene_id, contig in mapping.items():
                out.write(f"{sample}\t{gene_id}\t{contig}\n")
                total += 1

    print(f"\nDone. {total} gene->contig mapping(s) -> {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()

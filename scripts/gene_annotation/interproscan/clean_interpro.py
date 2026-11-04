#!/usr/bin/env python
"""
Parses raw InterProScan TSV and produces a per-gene summary table including GO terms.

Usage: 
python clean_interpro.py \\ 
    path/to/fasta \\
    path/to/interproscan/tsv
"""

from Bio import SeqIO
import sys


def main(argv):
    # Column index mapping for per-database hit rows (col 3 = database name)
    ipr_pos = {
        "CDD":         0,
        "Coils":       1,
        "FunFam":      2,
        "Gene3D":      3,
        "MobiDBLite":  4,
        "NCBIfam":     5,
        "PANTHER":     6,
        "Pfam":        7,
        "PIRSF":       8,
        "PRINTS":      9,
        "SFLD":        10,
        "SMART":       11,
        "SUPERFAMILY": 12,
        # slot 13 = InterPro summary (IPR accession + name)
        # slot 14 = GO terms
    }

    # Build seqs dict from FASTA — use only the first word of the header
    # to match InterProScan's gene ID format (handles ">gene1 description")
    seqs = {}
    for seq in SeqIO.parse(sys.argv[1], 'fasta'):
        gene_id = seq.id.split()[0]
        seqs[gene_id] = [[] for _ in range(15)]  # slots 0-14

    skipped = 0
    for line in open(sys.argv[2]):
        line = line.strip('\n')
        if not line or line.startswith('#'):
            continue

        cols = line.split('\t')
        # Use first word of gene ID to match FASTA keys        
        gene = cols[0].split()[0]

        if gene not in seqs:
            skipped += 1
            continue

        # InterPro summary entry (cols 11=IPR accession, col 12=IPR name)
        if len(cols) >= 13 and cols[11] != '-' and cols[11].startswith('IPR'):
            seqs[gene][13].append('-'.join(cols[11:13]))

        # GO terms (col 14)
        # Format: "GO:0005887(PANTHER)|GO:0016020(InterPro)|..."
        if len(cols) >= 14 and cols[13] != '-' and 'GO:' in cols[13]:
            for term in cols[13].split('|'):
                term = term.strip()
                if term and term != '-':
                    seqs[gene][14].append(term)

        # Per-database hits (col 3 = database name)
        if len(cols) >= 6:
            ident = cols[3]
            if ident in ipr_pos:
                pos = ipr_pos[ident]
                entry = '-'.join(cols[4:6]).strip('-')
                if entry:
                    seqs[gene][pos].append(entry)

    if skipped > 0:
        print(f"# WARNING: {skipped} rows skipped (gene ID not in FASTA)", file=sys.stderr)

    # Print header 
    print("Gene\tCDD-IPR\tCoils\tFunFam\tGene3D\tMobiDBLite\tNCBIfam\tPANTHER\t"
          "Pfam\tPIRSF\tPRINTS\tSFLD\tSMART\tSUPERFAMILY\tInterpro\tGO")

    # Print per-gene rows 
    for seq_id, entries in seqs.items():
        cols_out = [seq_id]
        for e in entries:
            if len(e) > 0 and len(e[0]) > 1:
                uniq = list(set(e))
                if len(uniq) < len(e):
                    add = ','.join(f"{j} ({e.count(j)})" for j in uniq)
                else:
                    add = ','.join(e)
                cols_out.append(add)
            else:
                cols_out.append('-')
        print('\t'.join(cols_out))


if __name__ == "__main__":
    main(sys.argv)

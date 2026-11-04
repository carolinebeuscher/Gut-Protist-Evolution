#! /usr/bin/env python
from ete3 import NCBITaxa
import sys

"""
Grabs lineage information of given species from NCBITaxa
Usage: 
    ncbi_taxonomy.py species_list.txt > output.tsv
"""

def main(argv):

    ncbi = NCBITaxa()

    # add header
    desired_ranks = ["genus", "family", "order", "class", "phylum", "domain", "clade_1", "clade_2", "clade_3"]
    print("species\t" + "\t".join(desired_ranks))

    for line in open(sys.argv[1]):
        line = line.strip('\n')
        species = line.split("\t")[0]

        if len(list(ncbi.get_name_translator([species]).values())) == 0:
            taxid = 1
        else:
            taxid = list(ncbi.get_name_translator([species]).values())[0][0]

        # get lineage and ranks
        lineage = ncbi.get_lineage(taxid)
        ranks = ncbi.get_rank(lineage)          # {taxid: rank_string}
        names = ncbi.get_taxid_translator(lineage)  # {taxid: sci_name}

        # Walk lineage in order (from root to species) to preserve all clade occurrences
        clades = []
        genus = family = order = class_name = phylum = domain = "NA"

        for tid in lineage:
            rank = ranks.get(tid, "no rank")
            name = names.get(tid, "NA")

            if rank == "genus":
                genus = name
            elif rank == "family":
                family = name
            elif rank == "order":
                order = name
            elif rank == "class":
                class_name = name
            elif rank == "phylum":
                phylum = name
            elif rank == "domain":
                domain = name
            elif rank == "clade":
                clades.append(name)

        # Pad clades list to always have 3 entries
        while len(clades) < 3:
            clades.append("NA")

        clade_1, clade_2, clade_3 = clades[0], clades[1], clades[2]

        print(f"{species}\t{genus}\t{family}\t{order}\t{class_name}\t{phylum}\t{domain}\t{clade_1}\t{clade_2}\t{clade_3}")


if __name__ == "__main__":
    main(sys.argv)

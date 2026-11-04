#!/usr/bin/env python

"""
Summarizes individual kofam outputs into one multi-species file

Usage: 
python kofam_summary.py \
    kegg_pathway_ids
    kegg_orth_and_pathway_corr
    list_of_filtered_kofam_files
"""

import sys, os


def main(argv):

    """Load KEGG pathways into a dict: {pathway_id: [name, [orthologs]]}"""
    pathways = {}
    with open(sys.argv[1]) as f:
          pathways = {line.split('\t')[0]: [line.split('\t')[1].strip(), []] 
                for line in f if line.strip()}



    """Add orthologs to pathways and remove empty pathways"""
    for line in open(sys.argv[2]):
        parts = line.strip().split('\t')
        path = parts[0]
        orth = parts[1]
        pathways[path][1].append(orth)
    remove = []
    for path in pathways:
        if len(pathways[path][1]) == 0:
            remove.append(path)

    for r in remove:
        del pathways[r]



    """Build Table"""
    headers = [os.path.basename(line.strip()).split('_')[0] for line in open(sys.argv[3])]
    toprint = "Pathway\tAnnot\tOrthcounts\tAll_orths\t"

    for h in headers:
        toprint += h + "_counts\t" + h + "_percent\t" + h + "_ids\t"
    print(toprint.strip('\t'))



    """ Create kofam dictionary """
    kofam_data = {}
    with open(sys.argv[3]) as f:
        for line in f:
            filename = line.strip()
            prefix = os.path.basename(line.strip()).split('_')[0]
            orthologs = set()
            with open(filename) as ann:
                for l in ann:
                    parts = l.split('\t')
                    if len(parts) > 2:
                        orthologs.add(parts[2])
            kofam_data[prefix] = orthologs



    """Reads KofamScan filtered files into table"""
    for path in pathways:
        pathway_members = pathways[path][1]
        toprint = path + '\t' + pathways[path][0].replace('\t','-') + '\t' + str(len(pathway_members)) + '\t' + ",".join(pathway_members) + '\t'
        for prefix, orthologs in kofam_data.items():
            found = [k for k in pathway_members if k in orthologs]
            toprint += str(len(found)) + '\t' + str(float(len(found)/len(pathway_members)) * 100) + '\t' + ",".join(found) + '\t'
        print(toprint.strip('\t'))

if __name__ == "__main__":
    main(sys.argv)


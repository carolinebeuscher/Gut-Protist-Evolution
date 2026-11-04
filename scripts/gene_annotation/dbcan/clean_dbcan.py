#! /usr/bin/env python
"""
processes dbCAN overview files (which contain predictions 
from three different CAZyme annotation tools).

takes the following input columns from the dbCAN overview file:
        geneid — the gene/protein identifier
        ec — EC number (read in but not used in the output)
        hmmer — CAZyme family prediction from HMMER
        dbcan — CAZyme family prediction from dbCAN-sub
        diamond — CAZyme family prediction from DIAMOND/BLAST
And it outputs: 
        geneid — the gene ID
        gh — the consensus CAZyme family ID
        Error disagreement:
                If predictions differ

Usage: 
clean_dbcan.py dbcan_out/$prefix"_overview.txt" > $prefix"_dbcan.txt"
"""

import sys
def main():

        for line in open(sys.argv[1]):
                line = line.strip('\n')
                if "#" not in line:

                        #gene id, ec, hmmer, dbcan, diamond

                        geneid = line.split('\t')[0]
                        ec = line.split('\t')[1]
                        hmmer = line.split('\t')[2]
                        dbcan = line.split('\t')[3]
                        diamond = line.split('\t')[4]

                        if "+" in line:
                                hmmer_parsed = [e.split('(')[0] for e in hmmer.split('+')]
                                dbcan_parsed = [e.split('_')[0] for e in dbcan.split('+')]
                                diamond_parsed = [e for e in diamond.split('+') if "." not in e]
                                all_ids = list(set(hmmer_parsed + dbcan_parsed + diamond_parsed))
                                if "-" in all_ids:
                                        all_ids.remove('-')

                                gh = "+".join(all_ids)
                        else:
                                hmmer_parsed = hmmer.split('(')[0]
                                dbcan_parsed = dbcan.split('_')[0]

                                all_ids = list(set([hmmer_parsed, dbcan_parsed, diamond]))
                                if "-" in all_ids:
                                        all_ids.remove('-')
                                if len(all_ids) > 1:
                                        #check if subfam
                                        major_fam = []
                                        for l in all_ids:
                                                if "_" in l:
                                                        major_fam.append(l.split('_')[0])
                                                else:
                                                        major_fam.append(l)
                                        if len(list(set(major_fam))) == 1:
                                                gh = list(set(major_fam))[0]
                                        else:
                                                print("ERROR: Disagreement " + line)
                                                continue
                                else:
                                        if not all_ids:
                                                continue
                                        gh = all_ids[0]
                        
                        print(geneid + '\t' + gh)
                        
if __name__ == "__main__":
  main()

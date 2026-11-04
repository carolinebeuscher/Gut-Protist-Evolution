#!/usr/bin/env python


from Bio import SeqIO
import sys
def main(argv):
    for seq in SeqIO.parse(sys.argv[1], 'fasta'):
        newname = "|".join(seq.description.split('|')[0:4])
        sequence = str(seq.seq)
        #print(newname)
        print(">" + newname + '\n' + sequence)

if __name__ == "__main__":
    main(sys.argv)
#! /usr/bin/env python

"""
Runs an alien index to hypothesize if HGT has occurred. 

requires a simple diamond self blast (bn):
    JUST USE OUTFMT 6 no extras. It's important the bitscore is the last column of the self blast or this doesn't work
and a more structured diamond blast with the following format: 
    diamond blastp \\
        --query [qfile] \\
        --db [db] \\
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen slen qlen nident positive staxids qcovhsp \\
        --threads [threads] \\
        -b12 \\
        -c1 \\
        --out [output] \\
        --max-target-seqs 500
    
    note: 
        max target seqs can be larger if desired but I do not recommend smaller.
        Default of 25 is much too low for valid inference.

Usage: 
    python alien_index_multiskip.py \\
        -bs "path/to/self/blast" \\
        -bn "path/to/nr or nr_cluster_seq blast" \\
        -r "Recipient taxonomy ID(s)" \\
        -s "Self taxonomy ID(s) to skip" \\
        -t "optional/path/to/ete3/ncbi/taxonomy/sqlite/db" \\
        -o "path/to/output/file" \\
        -e "path/to/file/with/taxids/to/exclude"

"""


from ete3 import NCBITaxa
import sys, os, re, argparse


_lineage_cache = {}

def get_lineage_cached(taxid, ncbi):
    if taxid not in _lineage_cache:
        _lineage_cache[taxid] = ncbi.get_lineage(int(taxid))
    return _lineage_cache[taxid]


def resolve_lineage_taxids(staxids_str, ncbi):
    lineage_taxids = []
    for taxid in staxids_str.split(';'):
        if taxid == "N/A" or len(taxid) == 0:
            return ["N/A"]
        for l in get_lineage_cached(taxid, ncbi):
            if l not in lineage_taxids:
                lineage_taxids.append(l)
    return lineage_taxids



def get_phylum_for_taxid(taxid, ncbi):
    try:
        lineage = get_lineage_cached(str(taxid), ncbi)
        ranks = ncbi.get_rank(lineage)
        rank_ids = {ranks[v]: v for v in ranks}
        if "phylum" in rank_ids:
            return list(ncbi.get_taxid_translator([rank_ids["phylum"]]).values())[0]
        elif "order" in rank_ids:
            return list(ncbi.get_taxid_translator([rank_ids["order"]]).values())[0]
        else:
            return "No phylum or order in NCBI taxid"
    except Exception:
        return "Taxid lookup error"


def format_skipped_hits(skipped_hits, ncbi):
    if not skipped_hits:
        return "-"
    entries = []
    for sseqid, bitscore, taxid in skipped_hits:
        try:
            species = list(ncbi.get_taxid_translator([int(taxid)]).values())[0]
        except Exception:
            species = "Unknown"
        phylum = get_phylum_for_taxid(taxid, ncbi)
        entries.append(f"{sseqid}|{species}|{phylum}|{bitscore}")
    return ",".join(entries)


def parse_blast(bls, skips, recipients, ncbi, exclude):
    max_bitscore = 0
    bitscores = []
    in_group_hits = []
    out_group_hits = []
    skipped_hits = []  #hits belonging to the skips lineage(s)

    for line in bls:
        line = line.strip("\n")
        qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, \
        sstart, send, evalue, bitscore, saccver, slen, qlen, nident, positive, \
        staxids, qcovs = line.split('\t')

        lineage_taxids = resolve_lineage_taxids(staxids, ncbi)

        if float(bitscore) > max_bitscore:
            max_bitscore = float(bitscore)

        excluded = len([e for e in lineage_taxids if e in exclude]) > 0

        if "N/A" in lineage_taxids or excluded:
            pass  #silently drop N/A and explicitly excluded hits as before
        elif any(s in lineage_taxids for s in skips):
            #collect skipped hits for the new output column
            skipped_hits.append([sseqid, float(bitscore), staxids.split(';')[0]])
        else:
            bitscores.append([float(bitscore), staxids.split(';')[0]])
            info = [sseqid, float(bitscore), float(evalue), staxids.split(';')[0]]
            if any(r in lineage_taxids for r in recipients):
                in_group_hits.append(info)
            else:
                out_group_hits.append(info)

    sorted_bitscores = sorted(bitscores, key=lambda x: float(x[0]), reverse=True)
    if len(sorted_bitscores) < 5:
        top_bitscores = sorted_bitscores
    else:
        top_bitscores = sorted_bitscores[0:4]

    classes = []
    for hit in top_bitscores:
        lineage = get_lineage_cached(hit[1], ncbi)
        ranks = ncbi.get_rank(lineage)
        rank_ids = {ranks[v]: v for v in ranks}
        if "phylum" in rank_ids:
            classes.append(list(ncbi.get_taxid_translator([rank_ids["phylum"]]).values())[0])
        elif "order" in rank_ids:
            classes.append(list(ncbi.get_taxid_translator([rank_ids["order"]]).values())[0])
        else:
            classes.append("No phylum or order in NCBI taxid")

    classes = set(classes)
    return float(max_bitscore), classes, in_group_hits, out_group_hits, skipped_hits


def calculate_alien_info(max_bitscore, in_lineage_hits, out_lineage_hits, ncbi):
    best_inlineage = ["", 0, float("inf")]
    best_outlineage = ["", 0, float("inf")]

    for hit in in_lineage_hits:
        if hit[1] > best_inlineage[1]:
            best_inlineage = hit
    for hit in out_lineage_hits:
        if hit[1] > best_outlineage[1]:
            best_outlineage = hit

    alien = float(float(best_outlineage[1] / max_bitscore) - float(best_inlineage[1] / max_bitscore))

    if len(in_lineage_hits) == 0:
        best_inlineage = "No in group hits"
        best_inlineage_species = "-"
        in_class_name = '-'
    else:
        best_inlineage_species = list(ncbi.get_taxid_translator([best_inlineage[3]]).values())[0]
        in_lineage = get_lineage_cached(best_inlineage[3], ncbi)
        in_ranks = ncbi.get_rank(in_lineage)
        in_rank_ids = {in_ranks[v]: v for v in in_ranks}
        if "phylum" in in_rank_ids:
            in_class_name = list(ncbi.get_taxid_translator([in_rank_ids["phylum"]]).values())[0]
        elif "order" in in_rank_ids:
            in_class_name = list(ncbi.get_taxid_translator([in_rank_ids["order"]]).values())[0]
        else:
            in_class_name = "No phylum in NCBI taxid"

    if len(out_lineage_hits) == 0:
        best_outlineage = "No non-group hits"
        best_outlineage_species = "-"
        out_class_name = "-"
    else:
        best_outlineage_species = list(ncbi.get_taxid_translator([best_outlineage[3]]).values())[0]
        out_lineage = get_lineage_cached(best_outlineage[3], ncbi)
        out_ranks = ncbi.get_rank(out_lineage)
        out_rank_ids = {out_ranks[v]: v for v in out_ranks}
        if "phylum" in out_rank_ids:
            out_class_name = list(ncbi.get_taxid_translator([out_rank_ids["phylum"]]).values())[0]
        elif "order" in out_rank_ids:
            out_class_name = list(ncbi.get_taxid_translator([out_rank_ids["order"]]).values())[0]
        else:
            out_class_name = "No phylum in NCBI taxid"

    return alien, best_inlineage, best_inlineage_species, in_class_name, \
           best_outlineage, best_outlineage_species, out_class_name


def process_and_write(gene, parse, self_bitscores, skips, recipients, recipient_name,
                      ncbi, exclude, dest):
    if gene not in self_bitscores:
        return

    max_bitscore = self_bitscores[gene]
    _max_blast_bitscore, best_lineages, in_lineage_hits, out_hits, skipped_hits = \
        parse_blast(parse, skips, recipients, ncbi, exclude)

    alien, best_inlineage, best_inlineage_species, in_class_name, \
    best_outlineage, best_outlineage_species, out_class_name = \
        calculate_alien_info(max_bitscore, in_lineage_hits, out_hits, ncbi)

    skipped_hits_str = format_skipped_hits(skipped_hits, ncbi)

    dest.write(gene + '\t' + str(alien) + '\t' + recipient_name + '\t' +
               str(best_inlineage[0]) + '\t' + best_inlineage_species + '\t' +
               in_class_name + '\t' + str(best_inlineage[1]) + '\t' +
               best_outlineage[0] + '\t' + best_outlineage_species + '\t' +
               out_class_name + '\t' + str(best_outlineage[1]) + '\t' +
               ",".join(best_lineages) + '\t' + str(max_bitscore) + '\t' +
               skipped_hits_str + '\n')


def main(argv):

    parser = argparse.ArgumentParser(description="Process command line arguments")
    parser.add_argument("-bs", "--blsself",   required=True,  help="Self blast")
    parser.add_argument("-bn", "--blsnr",     required=True,
                        help="Diamond (or blast) nr or nr_cluster_seq blast. "
                             "Command is: diamond blastp --query [qfile] --db [db] "
                             "--outfmt 6 qseqid sseqid pident length mismatch gapopen "
                             "qstart qend sstart send evalue bitscore slen slen qlen "
                             "nident positive staxids qcovhsp -b12 -c1 --max-target-seqs 500")
    parser.add_argument("-r", "--recipient",  required=True,
                        help="Recipient taxonomy ID(s). A single taxid (e.g. 2759 for Eukaryotes) "
                             "or a comma-separated list (e.g. 2759,4751). A hit is classified as "
                             "in-group if it belongs to ANY of the supplied taxids.")
    parser.add_argument("-s", "--self",       required=True,
                        help="Self taxonomy ID(s) to skip — hits belonging to these taxa are excluded "
                             "from both in-group and out-group scoring. A single taxid (e.g. 1234567) "
                             "or a comma-separated list (e.g. 1234567,9606). Usually species- or "
                             "genus-level.")
    parser.add_argument("-t", "--taxdb",
                        help="Optional argument - location for ete3 NCBI taxonomy sqlite database "
                             "(taxa.sqlite). If you are using a non current version of nr or similar "
                             "this is necessary, must provide taxonomy database that corresponds to nr version")
    parser.add_argument("-o", "--output",     required=True,  help="Output file")
    parser.add_argument("-e", "--exclude",    required=False,
                        help="file of taxIDs to exclude - useful if there is a spurious taxon "
                             "skewing results (e.g., contaminated MAG or other contaminated assembly)")

    args = parser.parse_args()

    blsself     = args.blsself
    blsnr       = args.blsnr
    taxdb       = args.taxdb
    destfile    = args.output
    excludefile = args.exclude

    #Parse self and recipient both as comma-separated lists; single values still work unchanged
    skips      = [int(s.strip()) for s in args.self.split(',')]
    recipients = [int(r.strip()) for r in args.recipient.split(',')]

    if taxdb != None:
        ncbi = NCBITaxa(taxdb)
    else:
        ncbi = NCBITaxa()

    if excludefile == None:
        exclude = []
    else:
        exclude = [int(line.strip('\n')) for line in open(excludefile)]

    #Build display names: single taxid resolves to its name; multiple taxids
    #are joined so the output column is human-readable
    recipient_names = [list(ncbi.get_taxid_translator([r]).values())[0] for r in recipients]
    recipient_name = ";".join(recipient_names)  #e.g. "Eukaryota" or "Eukaryota;Fungi"

    #Parse self-blast to get per-gene max bitscore
    self_bitscores = {}
    for line in open(blsself):
        line = line.strip('\n')
        q = line.split('\t')[0]
        s = line.split('\t')[1]
        if q == s:
            self_bitscores[q] = float(line.split('\t')[-1])

    parse = []
    currgene = ""

    dest = open(destfile, 'w')
    dest.write("query\talien_score\trecipient_name\t"
               "best_ingroup_hit\tbest_ingroup_species\tbest_ingroup_class\t"
               "best_ingroup_bitscore\tbest_outgroup_hit\tbest_outgroup_species\t"
               "best_outgroup_class\tbest_outgroup_bitscore\t"
               "classes_of_top_4_hits\tself_bitscore\tskipped_hits\n")

    for line in open(blsnr):
        line = line.strip('\n')
        gene = line.split('\t')[0]

        if gene != currgene:
            if currgene != "":
                print("parsing:", currgene)
                process_and_write(currgene, parse, self_bitscores, skips, recipients,
                                  recipient_name, ncbi, exclude, dest)
            currgene = gene
            parse = [line]
        else:
            parse.append(line)

    #once more for the last guy
    if currgene != "":
        print("parsing:", currgene)
        process_and_write(currgene, parse, self_bitscores, skips, recipients,
                          recipient_name, ncbi, exclude, dest)

    dest.close()


if __name__ == "__main__":
    main(sys.argv)

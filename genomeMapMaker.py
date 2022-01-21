"""
Author: Nathan Hall
File: genomeMapMaker.py

Program takes in .coords file and the orginal fasta files and produces a genome map file for use by ExpanGe.py
Currently, ExpanGe.py only makes 1 to 1 comparisons.
"""
import os
import sourmash
import argparse


class Accession :
    def __init__(self):
        self.length = None
        self.homolog = None
        self.kmerhash = None
        self.containment_scores = []


def parse_fasta(fasta):
    with open(fasta) as f:
        genome = f.read().lstrip('>').split('>')
    return genome

def parse_accesion(acc,accession_lens):
    accession = acc.splitlines()
    header = accession[0].split(' ')[0]
    assert header not in accession_lens, f'{header} is repeated twice in input file. Accession names must be unique'
    seq = ''.join(accession[1:])
    return (header,seq)



def update_accession_lens(acc, accession_lens,cutoff=1e5):
    header,seq = parse_accesion(acc,accession_lens)
    seq_len = len(seq)
    if seq_len >= cutoff:

        seqhash = sourmash.MinHash(n=0, ksize=51,scaled=1000)
        seqhash.add_sequence(seq.upper().replace('N',''))
        a = Accession()
        a.length = seq_len
        a.kmerhash = seqhash
        accession_lens[header] = a
    return accession_lens

def parse_query(query,accession_lens,query_cutoff):
    total = [['refChr','refLen','queryChr','queryLen','JaccardSimilarity']]
    genome_map = [['#refChr','refLen','queryChr','queryLen']]
    radded = []
    qadded = []
    for acc in query:
        max_key = []
        max_jaccard_value = 0
        header,seq = parse_accesion(acc,accession_lens)
        seq_len = len(seq)
        if seq_len > query_cutoff :
            qadded.append(header)
            query_hash = sourmash.MinHash(n=0,ksize=51,scaled=1000)
            query_hash.add_sequence(seq.upper().replace('N',''))
            for ref in list(accession_lens.keys()):
                jaccard_value = accession_lens[ref].kmerhash.jaccard(query_hash)
                total.append([ref,str(accession_lens[ref].length),header,str(seq_len),str(jaccard_value)])
                if jaccard_value > max_jaccard_value :
                    max_key = [ref]
                    max_jaccard_value = jaccard_value
                elif jaccard_value == max_key:
                    max_key.append(ref)
                radded += max_key
            for key in max_key:
                genome_map.append([
                  key,
                  str(accession_lens[key].length),
                  header,
                  str(seq_len)
                ])

    return (genome_map,total)


def main():
    parser = argparse.ArgumentParser(description='Take in a reference and query fasta'
                                                 'and return best hit for each query')

    parser.add_argument('-r,--reference',
                        help='fasta formatted genome file',
                        required=True,
                        dest='ref')
    parser.add_argument('-q,--query',
                        help='fasta formatted genome file',
                        required=True,
                        dest='query')
    parser.add_argument('-o,--out',
                        help='base name for output files. Two files with be written'
                             '{out_base}_total.tsv and {out_base}_genome_map.tsv',
                        required=True,
                        dest='out')
    parser.add_argument('--query_cutoff',
                        help='minimum length for query sequence to be included in map',
                        required=False,
                        default=1e6,
                        dest='query_cutoff')
    parser.add_argument('--reference_cutoff',
                        help='minimum length for reference sequence to be included in map',
                        required=False,
                        default=1e6,
                        dest='ref_cutoff')

    arg = parser.parse_args()
    query_cutoff = arg.query_cutoff
    out_base = arg.out.rstrip('_')
    rfasta = arg.ref
    qfasta = arg.query
    ref = parse_fasta(rfasta)
    query = parse_fasta(qfasta)
    """
    Add accessions to 
    """
    accession_lens = {}
    for acc in ref:
        accession_lens = update_accession_lens(acc, accession_lens)

    genome_map, total = parse_query(query, accession_lens,query_cutoff)
    with open(f'{out_base}_total.tsv', 'w') as f:
        for line in total:
            f.write('\t'.join(line) + '\n')

    with open(f'{out_base}_genome_map.tsv', 'w') as f:
        for line in genome_map:
            f.write('\t'.join(line) + '\n')

if __name__ == '__main__':
    main()





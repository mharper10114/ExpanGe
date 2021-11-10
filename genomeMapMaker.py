"""
Author: Nathan Hall
File: genomeMapMaker.py

Program takes in .coords file and the orginal fasta files and produces a genome map file for use by ExpanGe.py
Currently, ExpanGe.py only makes 1 to 1 comparisons.
"""

class Accession :
    def __init__(self):
        self.length = None
        self.homolog = None





def parse_fasta(fasta):
    with open(fasta) as f:
        genome = f.read().rstrip('>').split('>')
    return genome

def get_accession_lengths(genome):
    accession_lens = {}
    for accession in genome :
        header = accession[0]
        seq_len = len(''.join(accession[1:]))
        assert header not in accession_lens, f'{header} is repeated twice in input file. Accession names must be unique'
        accession = Accession()
        accession.length = seq_len
        accession_lens[header] = accession



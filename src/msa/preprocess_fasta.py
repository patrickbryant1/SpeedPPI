import argparse
import sys
import os
import numpy as np
import pandas as pd
import pdb

parser = argparse.ArgumentParser(description = '''Read a fasta with >1 seq.
                                                Write individual fasta files.
                                                Save all sequence information to a csv.
                                                ''')

parser.add_argument('--fasta_file', nargs=1, type= str, default=sys.stdin, help = 'Fasta.')
parser.add_argument('--outdir', nargs=1, type= str, default=sys.stdin, help = 'Path to output csv')

##############FUNCTIONS##############
def read_fasta(fasta_file):
    """Read a fasta file
    """

    ids = []
    seqs = []

    with open(fasta_file, 'r') as file:
        for line in file:
            line = line.rstrip()
            if '>' in line:
                ids.append(line[1:])
            else:
                seqs.append(line)

    return ids, seqs

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

def write_fasta(fasta_df, outdir):
    """Write individual fasta files
    """

    for ind, row in fasta_df.iterrows():
        with open(outdir+row.ID+'.fasta', 'w') as file:
            file.write('>'+row.ID+'\n')
            file.write(row.sequence)


##################MAIN#######################

#Parse args
args = parser.parse_args()
ids, seqs = read_fasta(args.fasta_file[0])
outdir = args.outdir[0]

#Create a df
fasta_df = pd.DataFrame()
fasta_df['ID'] = ids
fasta_df['sequence'] = seqs
#Save df
fasta_df.to_csv(outdir+'id_seqs.csv', index=None)
#Write individual fastas
write_fasta(fasta_df, outdir)

import argparse
import sys
import os
import numpy as np
import pandas as pd
import numpy as np
import time
import pdb

parser = argparse.ArgumentParser(description = '''Pick seqs from the human uniprot Fasta.
                                                ''')
parser.add_argument('--protein_fasta', nargs=1, type= str, default=sys.stdin, help = 'Path to fasta file with all proteins.')
parser.add_argument('--n_seqs', nargs=1, type= int, default=sys.stdin, help = 'How many seqs to pick.')
parser.add_argument('--outdir', nargs=1, type= str, default=sys.stdin, help = 'Path to output csv')

##############FUNCTIONS##############
def read_fasta(fasta_file):
    """Read a fasta file
    """

    ids, seqs = [], []
    with open(fasta_file, 'r') as file:
        for line in file:
            line = line.rstrip()
            if line[0]=='>':
                ids.append(line[1:])
                if len(ids)>1:
                    seqs.append(seq)
                seq=''
            else:
                seq+=line

    seqs.append(seq)
    return ids, seqs

##################MAIN#######################

#Parse args
args = parser.parse_args()
protein_fasta = args.protein_fasta[0]
n_seqs = args.n_seqs[0]
outdir = args.outdir[0]
#Read fasta
ids, seqs = read_fasta(protein_fasta)
#Sample
sel_inds = np.random.choice(len(ids), n_seqs, replace=False)
seq_df = pd.DataFrame()
seq_df['ID']=np.array(ids)[sel_inds]
seq_df['Sequence']=np.array(seqs)[sel_inds]
seq_df.to_csv(outdir+'selected_human_seqs.csv',index=None)

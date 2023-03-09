import argparse
import sys
import numpy as np
import pandas as pd
sys.path.insert('../msa/', 0)
import pdb
#Custom
from pair_msas import read_a3m, pair_msas, analyse_paired_msa


parser = argparse.ArgumentParser(description = '''Get the Neff and gap fractions for all MSA pairs in a csv.''')

parser.add_argument('--selected_ints', nargs=1, type= str, default=sys.stdin, help = 'Path to csv file with selected interacting sequences.')
parser.add_argument('--msa_dir', nargs=1, type= str, default=sys.stdin, help = 'Path to dir with MSAs.')



##################MAIN#######################

#Parse args
args = parser.parse_args()
#Data
selected_ints = pd.read_csv(args.selected_ints[0])
msa_dir = args.msa_dir[0]

for ind, row in selected_ints.iterrows():
    pdb.set_trace()
#Read MSAS
msa1='../../data/dev/4G4S_O.a3m'
msa2='../../data/dev/4G4S_P.a3m'
msa1, ox1 = read_a3m(msa1)
msa2, ox2 = read_a3m(msa2)
#Get the unique ox seqs from the MSAs
u_ox1, inds1 = np.unique(ox1, return_index=True)
u_msa1 = msa1[inds1]
u_ox2, inds2 = np.unique(ox2, return_index=True)
u_msa2 = msa2[inds2]
#Pair MSAs
paired_msa = pair_msas(u_ox1, u_ox2, u_msa1, u_msa2)
#Analyse
Neff, gap_fraction = analyse_paired_msa(paired_msa)
pdb.set_trace()

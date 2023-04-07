import argparse
import sys
import os
import numpy as np
import pandas as pd
import glob
import pdb

parser = argparse.ArgumentParser(description = '''Builds a PPI network from AlphaFold predictions.
                                                ''')

parser.add_argument('--pred_dir', nargs=1, type= str, default=sys.stdin, help = 'Path to directory with predictions.')
parser.add_argument('--pdockq_t', nargs=1, type= float, default=sys.stdin, help = 'pDockQ threshold for considering TP interactions.')
parser.add_argument('--outdir', nargs=1, type= str, default=sys.stdin, help = 'Path to output csv')

##############FUNCTIONS##############


##################MAIN#######################

#Parse args
args = parser.parse_args()
pred_dir = args.pred_dir[0]
pdockq_t = args.pdockq_t[0]
outdir = args.outdir[0]

all_preds = glob.glob(pred_dir+'pred*/*_metrics.csv')
all_ppis = []
for pred in all_preds:
    all_ppis.append(pd.read_csv(pred))
#Concat dfs
ppi_net = pd.concat(all_ppis)
#Save
ppi_net.to_csv(outdir+'all_ppis_unfiltered.csv', index=None)
print('Saved all PPIs before filtering on pDockQ to', outdir+'all_ppis_unfiltered.csv')
#Filter
ppi_net = ppi_net[ppi_net.pdockq>pdockq_t]
print('Filtered PPI network on pDockQ>',pdockq_t,'resulting in',len(ppi_net),'interactions.')
ppi_net.to_csv(outdir+'ppis_filtered.csv', index=None)
print('Saved all PPIs after filtering on pDockQ to', outdir+'ppis_filtered.csv')

# SpeedPPI

This repository contains a code for predicting the pairwise interaction network from a set of protein sequences. /
The procedure is based on MSA creation and evaluation with AlphaFold2 adapted for pairwise interactions aka [FoldDock](https://www.nature.com/articles/s41467-022-28865-w) \
\
The number of possible pairs grows exponentially with the number of input sequences according to: \
n*(n-1)/2, where n is the number of input sequences.
\
\
As this number rapidly becomes too large to handle, we apply a number of techniques to speed up the
network generation. We also reduce the memory footprint of the dependencies for the protein structure prediction
and generated MSAs.
\
\
Run the pipeline: \
Input:
1. A fasta file with sequences for all proteins you want to analyse \
2. Path to HHblits \
3. Output directory \
bash create_ppi.sh ./data/dev/test.fasta HHblits ./data/dev/
\
\
#Note
If you have a computational cluster available, it will be much faster to run your predictions in parallel. This requires some knowledge of computational infrastructure, however. Steps 2-4 in create_ppi.sh are written in individual scripts assuming a SLURM infrastructure. You can copy these, modify the paths and variables and queue them at your cluster to make the predictions even more efficient!

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




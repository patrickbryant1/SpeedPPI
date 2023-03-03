#Plot MSA pair
import argparse
import sys
import jax
import jax.numpy as jnp
import numpy as np
import pandas as pd

import pdb




parser = argparse.ArgumentParser(description = '''Plot an MSA pair.''')

parser.add_argument('--selected_ints', nargs=1, type= str, default=sys.stdin, help = 'Path to csv file with selected interacting sequences.')
parser.add_argument('--fasta_dir', nargs=1, type= str, default=sys.stdin, help = 'Path to dir with fasta seqs.')


##########FUNCIONS##########
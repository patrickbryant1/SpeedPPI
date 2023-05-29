# Copyright 2021 DeepMind Technologies Limited
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Full AlphaFold protein structure prediction script.
This has been modified to do all-vs-all predictions for PPI network creation.
"""
import json
import os
import warnings
import pathlib
import pickle
import random
import sys
import time
from typing import Dict, Optional

import argparse
from alphafold.common import protein
from alphafold.common import residue_constants
from alphafold.data import msaonly
from alphafold.data import foldonly
from alphafold.data import pipeline
from alphafold.data import templates
from alphafold.data import pair_msas
from alphafold.model import data
from alphafold.model import config
from alphafold.model import model
import pandas as pd
import numpy as np
#Data loading
from tinyloader import DataLoader
# Internal import (7716).
from collections import defaultdict
import pdb




parser = argparse.ArgumentParser(description = '''Predict a set of putative PPIs given two lists of single chain protein sequences.''')
parser.add_argument('--complex_id', nargs=1, type= str, default=sys.stdin, help = 'ID of complex.')
parser.add_argument('--msa1', nargs=1, type= str, default=sys.stdin, help = 'Path to dir with single chain MSAs.')
parser.add_argument('--msa2', nargs=1, type= str, default=sys.stdin, help = 'Path to dir with single chain MSAs.')
parser.add_argument('--data_dir', nargs=1, type= str, default=sys.stdin, help = 'Path to directory of supporting data (params).')
parser.add_argument('--max_recycles', nargs=1, type= int, default=10, help = 'Number of recyles through the model.')
parser.add_argument('--output_dir', nargs=1, type= str, default=sys.stdin, help = 'Path to a directory that will store the results.')


#######################FUNCTIONS#######################

##########INPUT DATA#########
def load_input_example(msa1, msa2, complex_id):
    """Load for a single example
    """


    # Get features. The features are prefetched on CPU.
    #Pair and block MSAs
    msa1, ox1 = pair_msas.read_a3m(msa1)
    msa2, ox2 = pair_msas.read_a3m(msa2)
    #Get the unique ox seqs from the MSAs
    u_ox1, inds1 = np.unique(ox1, return_index=True)
    u_msa1 = msa1[inds1]
    u_ox2, inds2 = np.unique(ox2, return_index=True)
    u_msa2 = msa2[inds2]
    #This is a list with seqs
    paired_msa = pair_msas.pair_msas(u_ox1, u_ox2, u_msa1, u_msa2)
    #Block the MSA
    blocked_msa = []
    gaps1 = '-'*len(msa2[0])
    for seq in msa1:
        blocked_msa.append(seq+gaps1)
    gaps2 = '-'*len(msa1[0])
    for seq in msa2:
        blocked_msa.append(gaps2+seq)

    # The msas must be str representations of the blocked+paired MSAs here
    #Define the data pipeline
    data_pipeline = foldonly.FoldDataPipeline()

    #Get the sequences from the msa
    cat_sequence = msa1[0]+msa2[0]

    #Get features
    feature_dict = data_pipeline.process(
          input_sequence=cat_sequence,
          input_description=complex_id,
          input_msas=[paired_msa,blocked_msa],
          template_search=None)

    # Introduce chain breaks for oligomers
    idx_res = feature_dict['residue_index']
    idx_res[len(msa1[0]):] += 200
    feature_dict['residue_index'] = idx_res #This assignment is unnecessary (already made?)
    # Add the id
    feature_dict['ID'] = complex_id
    return feature_dict, len(msa1[0])


#############Run PPI evaluation#############
def score_PPI(CB_coords, plddt, l1):
    """Score the PPI
    """

    #Cβs within 8 Å from each other from different chains are used to define the interface.
    CB_dists = np.sqrt(np.sum((CB_coords[:,None]-CB_coords[None,:])**2,axis=-1))

    #Get contacts
    contact_dists = CB_dists[:l1,l1:] #upper triangular --> first dim = chain 1
    contacts = np.argwhere(contact_dists<=8)

    #Get plddt per chain
    plddt1 = plddt[:l1]
    plddt2 = plddt[l1:]

    if contacts.shape[0]<1:
        pdockq=0
        avg_if_plddt=0
        n_if_contacts=0
    else:
        #Get the average interface plDDT
        avg_if_plddt = np.average(np.concatenate([plddt1[np.unique(contacts[:,0])], plddt2[np.unique(contacts[:,1])]]))
        #Get the number of interface contacts
        n_if_contacts = contacts.shape[0]
        x = avg_if_plddt*np.log10(n_if_contacts)
        pdockq = 0.724 / (1 + np.exp(-0.052*(x-152.611)))+0.018


    return pdockq, avg_if_plddt, n_if_contacts

def parse_atm_record(line):
    '''Get the atm record
    '''
    record = defaultdict()
    record['name'] = line[0:6].strip()
    record['atm_no'] = int(line[6:11])
    record['atm_name'] = line[12:16].strip()
    record['atm_alt'] = line[17]
    record['res_name'] = line[17:20].strip()
    record['chain'] = line[21]
    record['res_no'] = int(line[22:26])
    record['insert'] = line[26].strip()
    record['resid'] = line[22:29]
    record['x'] = float(line[30:38])
    record['y'] = float(line[38:46])
    record['z'] = float(line[46:54])
    record['occ'] = float(line[54:60])
    record['B'] = float(line[60:66])

    return record

def save_design(pdb_info, output_name, l1):
    '''Save the resulting protein-peptide design to a pdb file
    '''

    chain_name = 'A'
    with open(output_name, 'w') as f:
        pdb_contents = pdb_info.split('\n')
        for line in pdb_contents:
            try:
                record = parse_atm_record(line)
                if record['res_no']>l1:
                    chain_name='B'
                outline = line[:21]+chain_name+line[22:]
                f.write(outline+'\n')
            except:
                f.write(line+'\n')

def main(num_ensemble,
        complex_id,
        max_recycles,
        data_dir,
        msa1,
        msa2,
        output_dir):

  """Predict the structure of all possible interacting pairs to the protein in the target row.
  """

  #Define the data pipeline
  data_pipeline = foldonly.FoldDataPipeline()

  #Define the model runner - only once
  model_runners = {}
  for model_name in ['model_1']:
    model_config = config.model_config(model_name)
    model_config.data.eval.num_ensemble = num_ensemble
    model_config.data.common.num_recycle = max_recycles
    model_config.model.num_recycle = max_recycles
    model_params = data.get_model_haiku_params(model_name=model_name, data_dir=data_dir)
    model_runner = model.RunModel(model_config, model_params)
    model_runners[model_name] = model_runner

  #Get a seed
  random_seed = random.randrange(sys.maxsize)

  #Metrics
  metrics = {'ID':[], 'num_contacts':[], 'avg_if_plddt':[], 'pdockq':[]}


  #Merge fasta and predict the structure for each of the sequences.

  # Load an input example - on CPU
  feature_dict, l1 = load_input_example(msa1, msa2, complex_id)

  print('Evaluating pair', feature_dict['ID'])
  # Run the model - on GPU
  t0 = time.time()
  for model_name, model_runner in model_runners.items():
    processed_feature_dict = model_runner.process_features(
      feature_dict, random_seed=random_seed)
    prediction_result = model_runner.predict(processed_feature_dict)
  t1 = time.time()
  print('It took',t1-t0,'s to predict the interaction.')
  # Get pLDDT confidence metric.
  plddt = prediction_result['plddt']

  # Add the predicted LDDT in the b-factor column.
  # Note that higher predicted LDDT value means higher model confidence.
  plddt_b_factors = np.repeat(
    plddt[:, None], residue_constants.atom_type_num, axis=-1)
  unrelaxed_protein = protein.from_prediction(
    features=processed_feature_dict,
    result=prediction_result,
    b_factors=plddt_b_factors)

  #Get the pdb and CB coords
  pdb_info, CB_coords = protein.to_pdb(unrelaxed_protein)
  #Score - calculate the pDockQ (number of interface residues and average interface plDDT)
  pdockq, avg_if_plddt, n_if_contacts = score_PPI(CB_coords, plddt, l1)
  output_name = output_dir+feature_dict['ID']+'.pdb'
  save_design(pdb_info, output_name, l1)
  #Add to df
  metrics['ID'].append(feature_dict['ID'])
  metrics['num_contacts'].append(n_if_contacts)
  metrics['avg_if_plddt'].append(avg_if_plddt)
  metrics['pdockq'].append(pdockq)
  #Save df
  metric_df = pd.DataFrame.from_dict(metrics)
  metric_df.to_csv(output_dir+target_id+'_metrics.csv', index=None)
  print(feature_dict['ID'], pdockq)


################MAIN################
#Parse args
args = parser.parse_args()

main(num_ensemble=1,
    complex_id=args.complex_id[0],
    max_recycles=args.max_recycles[0],
    data_dir=args.data_dir[0],
    msa1=args.msa1[0],
    msa2=args.msa2[0],
    output_dir=args.output_dir[0])

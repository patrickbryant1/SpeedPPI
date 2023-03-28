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

from absl import app
from absl import flags
from absl import logging
from alphafold.common import protein
from alphafold.common import residue_constants
from alphafold.data import msaonly
from alphafold.data import foldonly
from alphafold.data import pipeline
from alphafold.data import templates
from alphafold.model import data
from alphafold.model import config
from alphafold.model import model
import pandas as pd
import pdb
import numpy as np
# Internal import (7716).

##### essential flags #####
flags.DEFINE_string('protein_csv', None,
    'Path to a csv with all proteins to be evaluated in an all-vs-all fashion.')

flags.DEFINE_int('target_row', None,
    'What row index to use to compare to all others in the protein csv.')

flags.DEFINE_list('fasta_dir', None,
    'Paths to FASTA files, each containing '
    'one sequence. Paths should be separated by commas. '
    'All FASTA paths must have a unique basename as the '
    'basename is used to name the output directories for '
    'each prediction.')

flags.DEFINE_string('output_dir', None,
    'Path to a directory that will store the results.')

flags.DEFINE_list('msa_dir', None,
    'Comma separated list of msa paths')

##### databases flags #####
flags.DEFINE_string('data_dir', None,
    'Path to directory of supporting data (params).')

flags.DEFINE_integer('max_recycles', 10,
    'Number of recyles through the model')

FLAGS = flags.FLAGS

#############Run PPI evaluation#############
def main(num_ensemble,
        max_recycles,
        data_dir,
        fasta_dir,
        msa_dir,
        output_dir,
        protein_csv,
        target_row):

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

  #Get length of the first chain
  target_row = protein_csv.loc[target_row]
  target_id = target_row.id
  target_seq = target_row.sequence
  chain_break=len(target_seq)
  #Get the remaining rows - only use the subsequent rows (upper-triangular)
  remaining_rows = np.arange(len(protein_csv))[target_row:]
  #Check the previous preds
  if os.path.exists(output_dir+target_id+'.csv'):
      metric_df = pd.read_csv(output_dir+target_id+'.csv')
      metrics = {'ID1': metric_df.ID1.values, 'ID2':metric_df.ID2.values,
                'num_contacts':metric_df.num_contacts.values, 'avg_if_plddt':metric_df.avg_if_plddt.values}
  else:
      metrics = {'ID1':[], 'ID2':[], 'num_contacts':[], 'avg_if_plddt':[]}

  #Merge fasta and predict the structure for each of the sequences.
  for i in remaining_rows[len(metrics):]:
    pdb.set_trace()
    row_i = protein_csv.loc[i]

    # Get features. The features are prefetched on CPU.
    # The msas must be str representations of the blocked+paired MSAs here
    feature_dict = data_pipeline.process(
          input_fasta_path=fasta_path,
          input_msas=msas,
          template_search=None)

    # Introduce chain breaks for oligomers
    idx_res = feature_dict['residue_index']
    idx_res[chain_break:] += 200
    feature_dict['residue_index'] = idx_res #This assignment is unnecessary (already made?)

    # Run the model - on GPU
    for model_name, model_runner in model_runners.items():
      processed_feature_dict = model_runner.process_features(
          feature_dict, random_seed=random_seed)
      prediction_result = model_runner.predict(processed_feature_dict)

    #Calculate the if contacts and if plDDT

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
    #Cβs within 8 Å from each other from different chains are used to define the interface.

    pdb.set_trace()




################MAIN################
main(num_ensemble=1,
    max_recycles=FLAGS.max_recycles,
    data_dir=FLAGS.data_dir,
    fasta_dir=FLAGS.fasta_dir,
    msa_dir=FLAGS.msa_dir,
    output_dir=FLAGS.output_dir,
    protein_csv=pd.read_csv(FLAGS.protein_csv),
    target_row=FLAGS.target_row)

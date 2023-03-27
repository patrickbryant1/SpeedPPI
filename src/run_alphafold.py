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

"""Full AlphaFold protein structure prediction script."""
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


def predict_structure(
    fasta_path: str,
    fasta_name: str,
    output_dir_base: str,
    data_pipeline: pipeline.DataPipeline,
    random_seed: int,
    model_runners: Optional[Dict[str, model.RunModel]]
    ):

  """Predicts structure using AlphaFold for the given sequence."""
  timings = {}
  output_dir = os.path.join(output_dir_base, fasta_name)
  if not os.path.exists(output_dir):
    os.makedirs(output_dir)
  msa_output_dir = os.path.join(output_dir, 'msas')
  if not os.path.exists(msa_output_dir):
    os.makedirs(msa_output_dir)

  # Get features.
  feature_dict = data_pipeline.process(
        input_fasta_path=fasta_path,
        input_msas=FLAGS.msas,
        template_search=FLAGS.template_search,
        msa_output_dir=msa_output_dir)
  timings['features'] = time.time() - t_0

  # Introduce chain breaks for oligomers ########## NEW!
  idx_res = feature_dict['residue_index']
  prev_overlay = 0
  for chain_break in FLAGS.chain_break_list:
    idx_res[chain_break:] += 200

  feature_dict['residue_index'] = idx_res

  relaxed_pdbs = {}
  plddts = {}

  # Run the models.
  for model_name, model_runner in model_runners.items():
    logging.info('Running model %s', model_name)
    t_0 = time.time()
    processed_feature_dict = model_runner.process_features(
        feature_dict, random_seed=random_seed)
    prediction_result = model_runner.predict(processed_feature_dict)

    #Calculate the if contacts and if plDDT

    # Get mean pLDDT confidence metric.
    plddt = prediction_result['plddt']


    # Add the predicted LDDT in the b-factor column.
    # Note that higher predicted LDDT value means higher model confidence.
    plddt_b_factors = np.repeat(
        plddt[:, None], residue_constants.atom_type_num, axis=-1)
    unrelaxed_protein = protein.from_prediction(
        features=processed_feature_dict,
        result=prediction_result,
        b_factors=plddt_b_factors)

    unrelaxed_pdb_path = os.path.join(output_dir, f'unrelaxed_{model_name}.pdb')
    with open(unrelaxed_pdb_path, 'w') as f:
      f.write(protein.to_pdb(unrelaxed_protein))

    #Score
    


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

  #Define the model runner
  model_runners = {}
  for model_name in ['model_1']:
    model_config = config.model_config(model_name)
    model_config.data.eval.num_ensemble = num_ensemble
    model_config.data.common.num_recycle = max_recycles
    model_config.model.num_recycle = max_recycles
    model_params = data.get_model_haiku_params(
    model_name=model_name, data_dir=data_dir)
    model_runner = model.RunModel(model_config, model_params)
    model_runners[model_name] = model_runner

  #Get a seed
  random_seed = random.randrange(sys.maxsize)

  #Go through all
  #Merge fasta and get length of the first chain

  # Predict structure for each of the sequences.
  for fasta_path, fasta_name in zip(FLAGS.fasta_paths, fasta_names):
    predict_structure(
        fasta_path=fasta_path,
        fasta_name=fasta_name,
        output_dir_base=FLAGS.output_dir,
        data_pipeline=data_pipeline,
        model_runners=model_runners,
        amber_relaxer=amber_relaxer,
        benchmark=FLAGS.benchmark,
        random_seed=random_seed)



################MAIN################

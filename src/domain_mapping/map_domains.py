import argparse
import sys
import os
import numpy as np
import pandas as pd
import json
import numpy as np
import tensorflow.compat.v1 as tf
from collections import Counter
import matplotlib.pyplot as plt
import tqdm

# Suppress noisy log messages.
from tensorflow.python.util import deprecation
deprecation._PRINT_DEPRECATION_WARNINGS = False
import pdb

parser = argparse.ArgumentParser(description = '''Map all Pfam domains in a protein sequence.
                                                This script assumes you have the model weights and mapping in
                                                the places specified in the setup.
                                                ''')
parser.add_argument('--protein_csv', nargs=1, type= str, default=sys.stdin, help = 'Path to csv file with all proteins to be evaluated in an all-vs-all fashion.')
parser.add_argument('--target_row', nargs=1, type= int, default=sys.stdin, help = 'What row index to use to compare to all others in the protein csv.')
parser.add_argument('--model_dir', nargs=1, type= str, default=sys.stdin, help = 'Path to model dir.')
parser.add_argument('--pfam_vocab', nargs=1, type= str, default=sys.stdin, help = 'Path to Pfam vocabulary.')
parser.add_argument('--outdir', nargs=1, type= str, default=sys.stdin, help = 'Path to output csv')


# Copyright 2021 Google Inc.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#From https://colab.research.google.com/github/google-research/google-research/blob/master/using_dl_to_annotate_protein_universe/neural_network/Neural_network_accuracy_on_random_seed_split.ipynb#scrollTo=Lp-3Ccx09IJ7
##############FUNCTIONS##############

def sequence_to_window(sequence, min_len=50, max_len=200):
    """Split a sequence into windows
    """
    seq_len = len(sequence)
    min_len = min(seq_len, min_len)
    max_len = min(seq_len, max_len)

    #If too short - use entire seq
    if min_len == max_len:
        return [sequence]

    #Go through and slice
    seq_slices, start_points, end_points = [], [], []
    for i in range(min_len, max_len+1):
        for j in range(seq_len-i):
            seq_slices.append(sequence[j:i+j])
            start_points.append(j)
            end_points.append(i+j)
    return seq_slices, start_points, end_points


AMINO_ACID_VOCABULARY = [
    'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R',
    'S', 'T', 'V', 'W', 'Y'
]
def residues_to_one_hot(amino_acid_residues):
  """Given a sequence of amino acids, return one hot array.

  Supports ambiguous amino acid characters B, Z, and X by distributing evenly
  over possible values, e.g. an 'X' gets mapped to [.05, .05, ... , .05].

  Supports rare amino acids by appropriately substituting. See
  normalize_sequence_to_blosum_characters for more information.

  Supports gaps and pads with the '.' and '-' characters; which are mapped to
  the zero vector.

  Args:
    amino_acid_residues: string. consisting of characters from
      AMINO_ACID_VOCABULARY

  Returns:
    A numpy array of shape (len(amino_acid_residues),
     len(AMINO_ACID_VOCABULARY)).

  Raises:
    ValueError: if sparse_amino_acid has a character not in the vocabulary + X.
  """
  to_return = []
  normalized_residues = amino_acid_residues.replace('U', 'C').replace('O', 'X')
  for char in normalized_residues:
    if char in AMINO_ACID_VOCABULARY:
      to_append = np.zeros(len(AMINO_ACID_VOCABULARY))
      to_append[AMINO_ACID_VOCABULARY.index(char)] = 1.
      to_return.append(to_append)
    elif char == 'B':  # Asparagine or aspartic acid.
      to_append = np.zeros(len(AMINO_ACID_VOCABULARY))
      to_append[AMINO_ACID_VOCABULARY.index('D')] = .5
      to_append[AMINO_ACID_VOCABULARY.index('N')] = .5
      to_return.append(to_append)
    elif char == 'Z':  # Glutamine or glutamic acid.
      to_append = np.zeros(len(AMINO_ACID_VOCABULARY))
      to_append[AMINO_ACID_VOCABULARY.index('E')] = .5
      to_append[AMINO_ACID_VOCABULARY.index('Q')] = .5
      to_return.append(to_append)
    elif char == 'X':
      to_return.append(
          np.full(len(AMINO_ACID_VOCABULARY), 1. / len(AMINO_ACID_VOCABULARY)))
    elif char == _PFAM_GAP_CHARACTER:
      to_return.append(np.zeros(len(AMINO_ACID_VOCABULARY)))
    else:
      raise ValueError('Could not one-hot code character {}'.format(char))
  return np.array(to_return)

def pad_one_hot_sequence(sequence: np.ndarray,
                         target_length: int) -> np.ndarray:
  """Pads one hot sequence [seq_len, num_aas] in the seq_len dimension."""
  sequence_length = sequence.shape[0]
  pad_length = target_length - sequence_length
  if pad_length < 0:
    raise ValueError(
        'Cannot set a negative amount of padding. Sequence length was {}, target_length was {}.'
        .format(sequence_length, target_length))
  pad_values = [[0, pad_length], [0, 0]]
  return np.pad(sequence, pad_values, mode='constant')

def batch_iterable(iterable, batch_size):
  """Yields batches from an iterable.

  If the number of elements in the iterator is not a multiple of batch size,
  the last batch will have fewer elements.

  Args:
    iterable: a potentially infinite iterable.
    batch_size: the size of batches to return.

  Yields:
    array of length batch_size, containing elements, in order, from iterable.

  Raises:
    ValueError: if batch_size < 1.
  """
  if batch_size < 1:
    raise ValueError(
        'Cannot have a batch size of less than 1. Received: {}'.format(
            batch_size))

  current = []
  for item in iterable:
    if len(current) == batch_size:
      yield current
      current = []
    current.append(item)

  # Prevent yielding an empty batch. Instead, prefer to end the generation.
  if current:
    yield current

def predict_domains(sequences,
                    start_points,
                    end_points,
                    model_dir,
                    pfam_vocab,
                    batch_size=32):
    """Predict domains
    """

    sess = tf.Session()
    graph = tf.Graph()

    with graph.as_default():
      saved_model = tf.saved_model.load(sess, ['serve'], model_dir)

    #top_pick_signature = saved_model.signature_def['serving_default']
    #top_pick_signature_tensor_name = top_pick_signature.outputs['output'].name

    class_confidence_signature = saved_model.signature_def['confidences']
    class_confidence_signature_tensor_name = class_confidence_signature.outputs['output'].name

    sequence_input_tensor_name = saved_model.signature_def['confidences'].inputs['sequence'].name
    sequence_lengths_input_tensor_name = saved_model.signature_def['confidences'].inputs['sequence_length'].name

    with open(pfam_vocab) as f:
      vocab = json.loads(f.read())

    inference_results = []
    inference_confidences = []
    batches = list(batch_iterable(sequences, batch_size))
    for batch in tqdm.tqdm(batches[:10], position=0):
      seq_lens = [len(seq) for seq in batch]
      one_hots = [residues_to_one_hot(seq) for seq in batch]
      padded_sequence_inputs = [pad_one_hot_sequence(seq, max(seq_lens)) for seq in one_hots]
      with graph.as_default():
        ir = sess.run(
            class_confidence_signature_tensor_name,
            {
                sequence_input_tensor_name: padded_sequence_inputs,
                sequence_lengths_input_tensor_name: seq_lens,
            })

      inference_results.extend(np.argmax(ir,axis=-1))
      inference_confidences.extend(np.max(ir,axis=-1))


    #Find domains
    #Go through all non-overlapping regions and find the highest confidences
    #Report these as domains following: https://www.nature.com/articles/s41587-021-01179-w#Sec4
    pdb.set_trace()

    #Map to Pfam vocab
    #Load vocab
    with open(pfam_vocab, 'r') as f:
        vocab = json.loads(f.read())
    pfam_mapping = [vocab[i] for i in inference_results]
    pdb.set_trace()


def find_max_confidence_regions(confidences, domain_pos, start_points, end_points):
    """Go through all the predicted confidences and annotate the non-overlapping regions
    """

    annotated_domains = []

    #Get the first domain
    fdp = np.argmax(confidences)
    annotated_domains.append(domain_pos[fdp])
    ads = start_points[fdp]
    ade = end_points[fdp]
    #Go through all domains and see if they are in another region
    for i in range(len(confidences)):
        ds_i = start_points[i]
        de_i = end_points[i]
        #Check overap (before start or after end)
        if de_i < ads or ds_i > ade):
            pdb.set_trace()


##################MAIN#######################

#Parse args
args = parser.parse_args()
protein_csv = pd.read_csv(args.protein_csv[0])
target_row = args.target_row[0]
model_dir = args.model_dir[0]
pfam_vocab = args.pfam_vocab[0]
outdir = args.outdir[0]

target_row = protein_csv.loc[target_row]
target_id, target_seq = target_row.ID, target_row.sequence
#Slice the sequence to search for domains
seq_slices, start_points, end_points = sequence_to_window(target_seq)
#Predict
mapped_domains = predict_domains(seq_slices, start_points, end_points, model_dir, pfam_vocab)
pdb.set_trace()

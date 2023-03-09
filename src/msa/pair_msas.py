import sys
import os
import numpy as np
import time
import pdb

def read_a3m(infile,max_gap_fraction=0.9):
    '''Read a3m MSA'''
    mapping = {'-': 21, 'A': 1, 'B': 21, 'C': 2, 'D': 3, 'E': 4, 'F': 5,
             'G': 6,'H': 7, 'I': 8, 'K': 9, 'L': 10, 'M': 11,'N': 12,
             'O': 21, 'P': 13,'Q': 14, 'R': 15, 'S': 16, 'T': 17,
             'V': 18, 'W': 19, 'Y': 20,'U': 21, 'Z': 21, 'X': 21, 'J': 21}

    parsed = []#Save extracted msa
    species = []
    seqlen = 0
    lc = 0
    with open(infile, 'r') as file:
        for line in file:
            line = line.rstrip()

            if line.startswith('>'): #OX=OrganismIdentifier
                if 'OX=' in line:
                    OX= line.split('OX=')[1]
                    if len(OX)>0:
                        try:
                            species.append(int(OX.split(' ')[0]))
                        except:
                            species.append(0)
                    else:
                        species.append(0)
                elif 'TaxID=' in line:
                    OX= line.split('TaxID=')[1]
                    if len(OX)>0:
                        try:
                            species.append(int(OX.split(' ')[0]))
                        except:
                            species.append(0)
                    else:
                        species.append(0)
                else:
                    species.append(0)
                continue
            line = line.rstrip()
            gap_fraction = line.count('-') / float(len(line))
            if gap_fraction <= max_gap_fraction:#Only use the lines with less than 90 % gaps
                parsed.append([mapping.get(ch, 22) for ch in line if not ch.islower()])
            else:
                if len(species)>1:
                    species = species[:-1] #Remove the previously stored species
                    continue
            #Check that the lengths match
            if len(parsed[-1])!=seqlen and lc>=1:
                parsed = parsed[:-1]
                species = species[:-1]
                continue
            seqlen = len(parsed[-1])
            lc+=1


    return np.array(parsed, dtype=np.int8, order='F'), np.array(species)

def pair_msas(ox1, ox2, msa1, msa2):
    '''Select the top ox match (first match) in each MSA and
    merge the sequences to a final MSA file.
    '''

    #Get the matches
    matching_ox = np.intersect1d(ox1,ox2)
    #Get the matching inds - only unique ox so can take directly
    #MSA1
    match_inds1 = np.argwhere(np.isin(ox1, matching_ox)).ravel()
    ox1, msa1 = ox1[match_inds1], msa1[match_inds1]
    #MSA2
    match_inds2 = np.argwhere(np.isin(ox2, matching_ox)).ravel()
    ox2, msa2 = msa2[match_inds2], msa2[match_inds2]

    #Return a concatenated MSA
    return np.concatenate([msa1, msa2], axis=1)

def analyse_paired_msa():
    """Analyse the paired MSA to see if it is likely
    to obtain an accurate prediction.
    1. Neff
    2. Gaps
    """

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

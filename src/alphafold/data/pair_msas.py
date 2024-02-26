import sys
import os
import numpy as np
import time
import copy
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
                parsed.append(''.join([ch for ch in line if not ch.islower()]))
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

    return np.array(parsed), np.array(species)

def pair_msas(ox1, ox2, msa1, msa2):
    '''Select the top ox match (first match) in each MSA and
    merge the sequences to a final MSA file.
    '''

    #Get the matches
    matching_ox = np.intersect1d(ox1,ox2)
    #Go through all matching and select the first (top) hit
    ind1 = [] #Index to select from the individual MSAs
    ind2 = []
    for ox in matching_ox:
        ind1.append(min(np.argwhere(ox1==ox)[:,0]))
        ind2.append(min(np.argwhere(ox2==ox)[:,0]))

    #Select from MSAs and merge
    cat_msa = np.array([x1+x2 for x1,x2 in zip(msa1[ind1],msa2[ind2])])

    return cat_msa

def analyse_paired_msa(paired_msa, seqid=0.62):
    """Analyse the paired MSA to see if it is likely
    to obtain an accurate prediction.
    1. Neff - 62% seqid
    2. Gaps
    """

    #Neff. Cluster values on 62% seqid
    Neff=0
    remaining_msa = copy.deepcopy(paired_msa)
    t = paired_msa.shape[1]*seqid #Threshold
    while remaining_msa.shape[0]>0:
        msa_diff = np.count_nonzero(remaining_msa-remaining_msa[0],axis=1)
        #Select
        remaining_msa = remaining_msa[np.argwhere(msa_diff>t)[:,0]]
        Neff+=1
        #print(remaining_msa.shape, Neff)

    #Gaps
    gap_fraction = np.argwhere(paired_msa==21).shape[0]/(paired_msa.shape[0]*paired_msa.shape[1])

    return Neff, gap_fraction


def write_a3m(merged_msa, outfile):
    '''Write a3m MSA'''
    backmap = { 1:'A', 2:'C', 3:'D', 4:'E', 5:'F',6:'G' ,7:'H',
               8:'I', 9:'K', 10:'L', 11:'M', 12:'N', 13:'P',14:'Q',
               15:'R', 16:'S', 17:'T', 18:'V', 19:'W', 20:'Y', 21:'-'} #Here all unusual AAs and gaps are set to the same char (same in the GaussDCA script)

    with open(outfile,'w') as file:
        for i in range(len(merged_msa)):
            file.write('>'+str(i)+'\n')
            file.write(''.join([backmap[ch] for ch in merged_msa[i]])+'\n')

    return None

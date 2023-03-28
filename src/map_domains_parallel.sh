#!/bin/bash -x
#SBATCH --output="path to out"
#SBATCH --error="path to error"
#SBATCH -A "your allocation id"
#SBATCH -t 00:10:00 #Set at 10 minutes
#SBATCH --array=1-N #N is the number of total proteins
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2 #How many CPUs to use
#SBATCH --mem-per-cpu=2500 #Memory in Mb per CPU

#Load the necessary modules (e.g. python)

#This checks if an offset is provided - useful if more than
#the max amount of allowed jobs is used
if [ -z $1 ]
then
        offset=0
else
        offset=$1
fi
LN=$(($SLURM_ARRAY_TASK_ID+$offset))
IDS="path to ids"
#Get ID
ID=$(sed -n $LN'p' $IDS)
echo $ID
#Fasta with sequences to use for the PPI network
FASTADIR="path to fasta seqs"
FASTA=$FASTADIR/$ID'.fasta'

#Annotate the domains

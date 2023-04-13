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
#the max amount of allowed jobs is used (0 if not provided)
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
OUTDIR="path to output"
MSADIR=$OUTDIR/msas/
mkdir $MSADIR
FASTADIR="path to individual fasta seqs created in step 1"
FASTA=$FASTADIR/$ID'.fasta'
HHBLITS="path to HHblits"
UNICLUST="path to Uniclust30, ../../data/uniclust30/uniclust30_2018_08 according to setup"


# Run HHblits to create MSA
echo $ID
if [ -f "$MSADIR/$ID.a3m" ]; then
  echo $MSADIR/$ID.a3m exists
else
  echo Creating MSA for $ID
  $HHBLITS -i $FASTA -d $UNICLUST -E 0.001 -all -oa3m $MSADIR/$ID'.a3m'
fi

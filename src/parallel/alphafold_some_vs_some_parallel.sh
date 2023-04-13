#!/bin/bash -x
#SBATCH --output="path to out"
#SBATCH --error="path to error"
#SBATCH -A "your allocation id"
#SBATCH -t 08:00:00 #Set at 8 hours
#SBATCH --array=1-N #N is the number of total proteins
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2 #How many CPUs to use
#SBATCH --mem-per-cpu=2500 #Memory in Mb per CPU
#SBATCH --export=ALL,CUDA_VISIBLE_DEVICES
#SBATCH --gpus=1


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
FASTADIR1="path to individual fasta seqs created in step 1 for list 1"
FASTADIR2="path to individual fasta seqs created in step 1 for list 2"

# Predict the structure using a modified version of AlphaFold2 (FoldDock)
# The predictions continue where they left off if the run is timed out
PR_CSV1=$FASTADIR1/id_seqs.csv
PR_CSV2=$FASTADIR2/id_seqs.csv
DATADIR="Path to where the AlphaFold2 parameteres are stored, default ../../data/"
RECYCLES=10
PDOCKQ_T=0.5
NUM_CPUS=2

mkdir $OUTDIR'/pred'$LN'/'
echo Running pred $c out of $NUM_PREDS
python3 "path to SpeedPPI"/src/run_alphafold_all_vs_all.py --protein_csv1 $PR_CSV1 \
--protein_csv2 $PR_CSV2 \
--target_row $LN \
--msa_dir $MSADIR \
--data_dir $DATADIR \
--max_recycles $RECYCLES \
--pdockq_t $PDOCKQ_T \
--num_cpus $NUM_CPUS \
--output_dir $OUTDIR'/pred'$LN'/'

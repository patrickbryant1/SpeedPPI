#Create a PPI network

#ARGS
#INPUT
FASTA_SEQS=$1 #All fasta seqs
HHBLITS=$2 #Path to HHblits
OUTDIR=$3
#DEFAULT
UNICLUST=./data/uniclust30/uniclust30_2018_08 #Assume path according to setup


#The pipeline starts here
#1. Create individual fastas
FASTADIR=$OUTDIR/fasta/
if [ -f "$FASTADIR/id_seqs.csv" ]; then
  echo Fastas exist...
  echo "Remove the directory $FASTADIR if you want to write new fastas."
else
mkdir $FASTADIR
python3 ./src/preprocess_fasta.py --fasta_file $FASTA_SEQS \
--outdir $FASTADIR
echo "Writing fastas of each sequence to $FASTADIR"
fi
wait

#2. Run HHblits for all fastas to create MSAs
MSADIR=$OUTDIR/msas/
if [ -d "$MSADIR" ]; then
  echo MSAs exists...
  echo Checking if all are present
else
  mkdir $MSADIR
fi
for FASTA in $FASTADIR/*.fasta
do
  ID=$(basename $FASTA)
  ID=$(echo $ID|cut -d '.' -f 1)
  echo $ID
  if [ -f "$MSADIR/$ID.a3m" ]; then
    echo $MSADIR/$ID.a3m exists
  else
    echo Creating MSA for $ID
    $HHBLITS -i $FASTA -d $UNICLUST -E 0.001 -all -oa3m $MSADIR/$ID'.a3m'
  fi
done

#3. Predict the structure using a modified version of AlphaFold2 (FoldDock)
PR_CSV=$FASTADIR/id_seqs.csv
NUM_PREDS=$(wc -l $PR_CSV|cut -d ' ' -f 1)
NUM_PREDS=$(($NUM_PREDS-1))
DATADIR=./data/
RECYCLES=10
PDOCKQ_T=0.5
NUM_CPUS=1
for (( c=1; c<=$NUM_PREDS; c++ ))
do
  mkdir $OUTDIR'/pred'$c'/'
  echo Running pred $c out of $NUM_PREDS
  python3 ./src/run_alphafold.py --protein_csv $PR_CSV \
    --target_row $c \
    --msa_dir $MSADIR \
    --data_dir $DATADIR \
    --max_recycles $RECYCLES \
    --pdockq_t $PDOCKQ_T \
    --num_cpus $NUM_CPUS \
    --output_dir $OUTDIR'/pred'$c'/'
done

#4. Merge all predictions to construct a PPI network.
#When the pDockQ > 0.5, the PPV is >0.9 (https://www.nature.com/articles/s41467-022-28865-w, https://www.nature.com/articles/s41594-022-00910-8)
#The default threshold to construct edges (links) is 0.5

#Create a PPI network

#ARGS
#INPUT
FASTA_SEQS1=$1 #All fasta seqs from list 1
FASTA_SEQS2=$2 #All fasta seqs from list 2
HHBLITS=$3 #Path to HHblits
PDOCKQ_T=$4
OUTDIR=$5
#DEFAULT
UNICLUST=./data/uniclust30/uniclust30_2018_08 #Assume path according to setup


#The pipeline starts here
#1. Create individual fastas

#List 1
FASTADIR=$OUTDIR/fasta1/
if [ -f "$FASTADIR/id_seqs.csv" ]; then
  echo Fastas exist...
  echo "Remove the directory $FASTADIR if you want to write new fastas."
else
mkdir -p $FASTADIR
python3 ./src/preprocess_fasta.py --fasta_file $FASTA_SEQS1 \
--outdir $FASTADIR
echo "Writing fastas of each sequence to $FASTADIR"
fi

#List 2
FASTADIR=$OUTDIR/fasta2/
if [ -f "$FASTADIR/id_seqs.csv" ]; then
  echo Fastas exist...
  echo "Remove the directory $FASTADIR if you want to write new fastas."
else
mkdir $FASTADIR
python3 ./src/preprocess_fasta.py --fasta_file $FASTA_SEQS2 \
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
#List 1
for FASTA in $OUTDIR/fasta1/*.fasta
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

#List 2
for FASTA in $OUTDIR/fasta2/*.fasta
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
PR_CSV1=$OUTDIR/fasta1/id_seqs.csv
PR_CSV2=$OUTDIR/fasta2/id_seqs.csv
NUM_PREDS=$(wc -l $PR_CSV1|cut -d ' ' -f 1)
NUM_PREDS=$(($NUM_PREDS-1))
DATADIR=./data/
RECYCLES=10
NUM_CPUS=1
for (( c=1; c<=$NUM_PREDS; c++ ))
do
  mkdir $OUTDIR'/pred'$c'/'
  echo Running pred $c out of $NUM_PREDS
  python3 ./src/run_alphafold_some_vs_some.py --protein_csv1 $PR_CSV1 \
  --protein_csv2 $PR_CSV2 \
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
python3 ./src/build_ppi.py --pred_dir $OUTDIR/ \
--pdockq_t $PDOCKQ_T --outdir $OUTDIR/

#5. Move all high-confidence predictions to a dir
mkdir $OUTDIR'/high_confidence_preds/'
mv $OUTDIR/pred*/*.pdb $OUTDIR'/high_confidence_preds/'
echo Moved all high confidence predictions to $OUTDIR'/high_confidence_preds/'

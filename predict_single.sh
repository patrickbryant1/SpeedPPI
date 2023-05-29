#Run a single PPI prediction

#ARGS
#INPUT
FASTA1=$1 #Fasta sequence 1
FASTA2=$2 #Fasta sequence 2
HHBLITS=$3 #Path to HHblits
PDOCKQ_T=$4
OUTDIR=$5
#DEFAULT
UNICLUST=./data/uniclust30_2018_08/uniclust30_2018_08 #Assume path according to setup

#2. Run HHblits for all fastas to create MSAs
MSADIR=$OUTDIR/msas/
if [ -d "$MSADIR" ]; then
  echo MSAs exists...
  echo Checking if all are present
else
  mkdir $MSADIR
fi
#Create MSAs
for FASTA in $FASTA1 $FASTA2
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

#Predict a single example
ID1=$(basename $FASTA1)
ID2=$(basename $FASTA2)
MSA1=$MSADIR/$ID1'.a3m'
MSA2=$MSADIR/$ID2'.a3m'
DATADIR=./data/
RECYCLES=10
NUM_CPUS=1

echo Predicting...
python3 ./src/run_alphafold_some_vs_some.py --protein_csv1 $PR_CSV1 \
--protein_csv2 $PR_CSV2 \
--target_row $c \
--msa_dir $MSADIR \
--data_dir $DATADIR \
--max_recycles $RECYCLES \
--pdockq_t $PDOCKQ_T \
--num_cpus $NUM_CPUS \
--output_dir $OUTDIR'/pred'$c'/'

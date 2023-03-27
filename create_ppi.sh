#Create a PPI network

#Path to settings
FASTA_SEQS=$1 #All fasta seqs
HHBLITS=$2
OUTDIR=$3 #Path to HHblits

UNICLUST=./data/uniclust30/uniclust30_2018_08 #Assume path according to setup


#The pipeline starts here
#Create individual fastas
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

#Run HHblits for all MSAs
MSADIR=$OUTDIR/msas/
if [ -d "$MSADIR" ]; then
  echo MSAS exists...
  echo Remove $MSADIR if you want to make new MSAs.
else
  mkdir $MSADIR
for FASTA in $FASTADIR/*.fasta
do
  ID=$(basename $FASTA)
  ID=$(echo $ID|cut -d '.' -f 1)
  echo $ID
  if [ -f "$MSADIR/$ID.a3m" ]; then
    echo $MSADIR/$ID.a3m exists
  else
    $HHBLITS -i $FASTA -d $UNICLUST -E 0.001 -all -oa3m $MSADIR/$ID'.a3m'
  fi
done
fi

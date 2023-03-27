#Create a PPI network


#Run HHblits for all MSAs
hhblits -i $FASTA1 -d $BASEDIR/uniclust30_2018_08/uniclust30_2018_08 -E 0.001 -all -oa3m ${FASTA1%.fasta}.a3m

#Get the Neff and gap fractions of the MSAs in the Marks set
SEL_INTS=../../data/dev/AF_dockQ.csv
MSA_DIR=/group/ag_cmb/pbryant/data/PPI/hhblits_msas/ #Local
python3 ./get_msa_metrics.py --selected_ints $SEL_INTS \
--msa_dir $MSA_DIR

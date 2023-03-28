#Get requirements
#conda install -c conda-forge -c bioconda hhsuite
#git clone https://github.com/soedinglab/hh-suite.git
#mkdir -p hh-suite/build && cd hh-suite/build
#cmake -DCMAKE_INSTALL_PREFIX=. ..
#make -j 4 && make install

#Download  Uniclust30
echo "Getting unlclust30..."
wget http://wwwuser.gwdg.de/~compbiol/uniclust/2018_08/uniclust30_2018_08_hhsuite.tar.gz
mkdir data/uniclust30
mv uniclust30_2018_08_hhsuite.tar.gz data/uniclust30
tar -zxvf data/uniclust30/uniclust30_2018_08_hhsuite.tar.gz

#Download AF2 parameters
echo "Getting AlphaFold parameters (v. 2021.07.14)..."
mkdir data/params
wget https://storage.googleapis.com/alphafold/alphafold_params_2021-07-14.tar
mv alphafold_params_2021-07-14.tar data/params
tar -xf data/params/alphafold_params_2021-07-14.tar

#Cleanup
echo "Cleaning up unnecessary files..."
rm data/uniclust30/uniclust30_2018_08_hhsuite.tar.gz
rm data/params/alphafold_params_2021-07-14.tar

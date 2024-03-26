#!/bin/bash

printf "\n################ Build dependency files ################\n" 


rootdir=`realpath "$0"`
rootdir=`dirname $rootdir`
rootdir="${rootdir/"/LigExtract/bin"/}"

# download BeEM from https://github.com/kad-ecoli/BeEM/
wget https://github.com/kad-ecoli/BeEM/releases/download/v1.0.1/BeEM.linux -O "$rootdir"/LigExtract/data/bin/BeEM.linux
chmod a+x "$rootdir"/LigExtract/data/bin/BeEM.linux

# Dependency files - one-off

mkdir -p "$rootdir"/LigExtract/data


# https://www.ebi.ac.uk/pdbe/docs/sifts/quick.html
echo "Get pdb_chain_uniprot.csv.gz"
wget ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/pdb_chain_uniprot.csv.gz -O $rootdir/LigExtract/data/pdb_chain_uniprot.csv.gz --quiet
echo "Get entries.idx"
wget https://files.wwpdb.org/pub/pdb/derived_data/index/entries.idx -O $rootdir/LigExtract/data/allpdbs.txt --quiet
echo "Get prd-all.cif.gz"
wget ftp://ftp.wwpdb.org/pub/pdb/data/bird/prd/prd-all.cif.gz -O $rootdir/LigExtract/data/prd-all.cif.gz --quiet
gunzip $rootdir/LigExtract/data/prd-all.cif.gz

python $rootdir/LigExtract/bin/get_prd2pdb.py $rootdir/LigExtract/data/prd-all.cif
echo "Get cc-counts.td from ligand Expo"
wget http://ligand-expo.rcsb.org/dictionaries/cc-counts.tdd -O $rootdir/LigExtract/data/all_pdbligs.txt --quiet




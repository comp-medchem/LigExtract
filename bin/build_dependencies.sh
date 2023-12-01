#!/bin/bash

printf "\n################ Build dependency files ################\n" 

# Dependency files - one-off
mkdir -p ~/LigExtract/data
# https://www.ebi.ac.uk/pdbe/docs/sifts/quick.html
wget ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/pdb_chain_uniprot.csv.gz -O ~/LigExtract/data/pdb_chain_uniprot.csv.gz --quiet
wget https://files.wwpdb.org/pub/pdb/derived_data/index/entries.idx -O ~/LigExtract/data/allpdbs.txt --quiet
wget http://ftp.wwpdb.org/pub/pdb/data/bird/prd/prd-all.cif.gz -O ~/LigExtract/data/prd-all.cif.gz --quiet
gunzip ~/LigExtract/data/prd-all.cif.gz

python ~/LigExtract/bin/get_prd2pdb.py ~/LigExtract/data/prd-all.cif
wget http://ligand-expo.rcsb.org/dictionaries/cc-counts.tdd -O ~/LigExtract/data/all_pdbligs.txt --quiet

#!/bin/bash

printf "\n################ Build dependency files ################\n" 

# Dependency files - one-off
mkdir -p ../data
# https://www.ebi.ac.uk/pdbe/docs/sifts/quick.html
wget ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/pdb_chain_uniprot.csv.gz -O ../data/pdb_chain_uniprot.csv.gz --quiet
wget https://ftp.wwpdb.org/pub/pdb/derived_data/index/entries.idx -O ~/data/allpdbs.txt --quiet
# Download the full PRD cif file at http://ftp.wwpdb.org/pub/pdb/data/bird/prd/prd-all.cif.gz; gunzip prd-all.cif.gz
python ~/scripts/get_prd2pdb.py ~/data/prd-all.cif
wget http://ligand-expo.rcsb.org/dictionaries/cc-counts.tdd -O ~/data/all_pdbligs.txt --quiet
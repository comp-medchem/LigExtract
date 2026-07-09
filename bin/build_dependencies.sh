#!/bin/bash
set -e

length=90; padding=$(printf '%*s' "$length" '' | tr ' ' '#')

title=" Build dependency files "
printf "%.*s %s %.*s\n" "$(((length - 1 - ${#title}) / 2))" "$padding" "$title" "$(((length - ${#title}) / 2))" "$padding"


rootdir=`realpath "$0"`
rootdir=`dirname $rootdir`
rootdir="${rootdir/"/LigExtract/bin"/}"


# download BeEM from https://github.com/kad-ecoli/BeEM/
echo "Get BeEM"
wget https://github.com/kad-ecoli/BeEM/releases/download/v1.0.1/BeEM.linux -O "$rootdir"/LigExtract/bin/BeEM.linux --quiet
chmod a+x "$rootdir"/LigExtract/bin/BeEM.linux

# Dependency files - one-off

mkdir -p "$rootdir"/LigExtract/data


# https://www.ebi.ac.uk/pdbe/docs/sifts/quick.html
echo "Get pdb_chain_uniprot.csv.gz"
wget ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/pdb_chain_uniprot.csv.gz -O $rootdir/LigExtract/data/pdb_chain_uniprot.csv.gz --quiet
echo "Get entries.idx"
wget https://files.wwpdb.org/pub/pdb/derived_data/index/entries.idx -O $rootdir/LigExtract/data/allpdbs.txt --quiet
echo "Get prd-all.cif.gz"
wget https://files.wwpdb.org/pub/pdb/data/bird/prd/prd-all.cif.gz -O $rootdir/LigExtract/data/prd-all.cif.gz --quiet
gunzip $rootdir/LigExtract/data/prd-all.cif.gz

python $rootdir/LigExtract/bin/get_prd2pdb.py $rootdir/LigExtract/data/prd-all.cif
echo "Get ligand counts across all PDB" # Ligand Expo is no longer maitained
# from https://github.com/rcsb/rcsb-training-resources/blob/master/example-use-cases/pdb-ligand-composition
python $rootdir/LigExtract/bin/generate_pdb_ligand_mappings.py --chem_comp_types nonpolymer # --generate_cc_extra_file

python -c "import pandas as pd; d = [x.strip().split('\t') for x in open('cc-to-pdb.tsv').readlines()]; d=pd.DataFrame([[x[0], len(x[1].split(' '))] for x in d], columns=['id', 'count']); d.to_csv('all_pdbligs.txt', index=False, sep='\t')"
mv all_pdbligs.txt $rootdir/LigExtract/data/.
rm pdb-to-cc.tsv cc-to-pdb.tsv #cc-counts-extra.tsv 

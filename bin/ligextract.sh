#!/bin/bash
set -e

# /path/to/ligextract.sh -d myproteins -r 3.5 -o cluster -c nan
# /path/to/ligextract.sh -d myproteins -r 3.5 -o filter -c clean

while getopts ':h?d:r:o:c:' opts; do
  case ${opts} in

  d) d="$OPTARG";;
	r) res="$OPTARG";;
	o) filter_option="$OPTARG";;
	c) cleanup="$OPTARG";;
    \?|h|*) 
		echo "Please use:"
        echo "   bash ligextract.sh -d [directory with pdbs] -r [maximum PDB resolution accepted] -o ['filter' or 'cluster'] -c ['clean']"
        echo "   e.g. bash ligextract.sh -d path/to/pdbsdir -r 3 -o filter"
        exit 0
      ;;
  esac
done

rootdir=`realpath "$0"`
rootdir=`dirname $rootdir`
rootdir="${rootdir/"/LigExtract/bin"/}"

# clear out all log files before a new run
rm -f *.log
#rm -fr $d
rm -fR "$d"_LIGS

length=100; padding=$(printf '%*s' "$length" '' | tr ' ' '#')

title=" Processing [${d}] Query "
printf "%.*s %s %.*s\n\n\n" "$(((length - 1 - ${#title}) / 2))" "$padding" "$title" "$(((length - ${#title}) / 2))" "$padding"


######################################################################################################################
###################################################### MODULE 1 ######################################################
title=" MODULE 1: Structure retrieval and processing "
printf "%.*s %s %.*s\n" "$(((length - 1 - ${#title}) / 2))" "$padding" "$title" "$(((length - ${#title}) / 2))" "$padding"

##download all PDBs
python $rootdir/LigExtract/bin/getPdbsFromUniprot.py --outputDir $d --uniprots "$d"_uniprot_list.txt --allPdbs $rootdir/LigExtract/data/allpdbs.txt --maxResol $res

rm -f pdbs2download.txt
mkdir -p cifs
python $rootdir/LigExtract/bin/downloadPdbs.py --outputDir $d 

touch pdb_download.log
if [[ -f pdbs2download.txt ]]; then 
	  bash $rootdir/LigExtract/bin/pdb_batch_download.sh -f pdbs2download.txt -c
fi

# catch failed download and update *_pdb_uniprot_filteredlist.txt
python $rootdir/LigExtract/bin/bypass_missing_pdbs.py

### RENAMING PDBs to lowercase: no longer needed

if compgen -G "*.gz" > /dev/null; then
    gunzip *.gz
    #for file in *.pdb; do mv -- $file $(echo $file | tr 'A-Z' 'a-z'); done
    mv *.cif cifs/.
fi

# convert all cifs to PDB to access ATOM and HETATM
printf "\nConvert mmCIF to pdb\n"
for cif in cifs/*cif; do
	pdbname=$(basename "$cif" .cif).pdb
	if [ ! -f "$d/$pdbname" ]; then
		$rootdir/LigExtract/bin/BeEM.linux $cif >> cifpdbconvert.log
		mv $pdbname $d/.	
	fi
done

if compgen -G "*ligand-id-mapping.tsv" > /dev/null; then
    mv *ligand-id-mapping.tsv cifs/.
fi

#mv *chain-id-mapping.txt $d/.

# handle all 5-letter cases
printf "\nProcess 5-letter codes\n"
python $rootdir/LigExtract/bin/processLongCodes.py


> "$d"_process_uniprot_chains.err
> "$d"_process_uniprot_chains.txt
python $rootdir/LigExtract/bin/process_chains.py --pdbspath $d


## check if some pdbs failed
redo_pdbs=`grep "Request to Uniprot failed" "$d"_process_uniprot_chains.err | cut -f1 -d":"`
echo `grep "Request to Uniprot failed" "$d"_process_uniprot_chains.err | wc -l` "failed PDBs."
sed -i '/Request to Uniprot failed after/d' "$d"_process_uniprot_chains.err
if [[ $redo_pdbs != "" ]]; then 
	echo "LigExtract will repeat chain processing for failed PDBs"
	python $rootdir/LigExtract/bin/process_chains.py --pdbspath $d --pdbsToRedo $redo_pdbs; 
fi 


######################################################################################################################
###################################################### MODULE 2 ######################################################

padding=$(printf '%*s' "$length" '' | tr ' ' '#')

title=" MODULE 2: Ligand Extraction "
printf "\n\n%.*s %s %.*s\n\n" "$(((length - 1 - ${#title}) / 2))" "$padding" "$title" "$(((length - ${#title}) / 2))" "$padding"


python $rootdir/LigExtract/bin/extract_ligands.py --pdbPath $d --outputPath "$d"_LIGS --uniprot2pdbFile "$d"_pdb_uniprot_filteredlist.txt > ligand_extraction.log


ls "$d"_LIGS > rawlist_extraction.txt

echo "$filter_option option selected"


if [ $filter_option == "filter" ]; then
	#############################################   MODULE 3  ##################################################
	title=" MODULE 3: Ligand Selection and Curation "
    printf "%.*s %s %.*s\n" "$(((length - 1 - ${#title}) / 2))" "$padding" "$title" "$(((length - ${#title}) / 2))" "$padding"
    printf "\n First Ligand Clean-up (crystallography additives, solvents, etc) & Pocket detection \n"

	python $rootdir/LigExtract/bin/find_ligands.py --pdbPath $d --ligandsPath "$d"_LIGS --dist 6 --uniprot2pdbFile "$d"_pdb_uniprot_filteredlist.txt --keeprepeats n > find_ligands.log
	

    #############################################   MODULE 4  ##################################################
	title=" MODULE 4: Final Ligand Selection "
    printf "%.*s %s %.*s\n" "$(((length - 1 - ${#title}) / 2))" "$padding" "$title" "$(((length - ${#title}) / 2))" "$padding"

	python $rootdir/LigExtract/bin/filter_ligands.py --pdbPath $d --ligandsPath "$d"_LIGS --prdCif $rootdir/data/prd-all.cif > filter_ligands.log
	
	python $rootdir/LigExtract/bin/assemble_finalreport.py --pdbPath $d --ligandsPath "$d"_LIGS
fi



if [ $filter_option == "cluster" ]; then
	#############################################   MODULE 3  ##################################################
	title=" MODULE 3: Ligand Selection and Curation "
    printf "%.*s %s %.*s\n" "$(((length - 1 - ${#title}) / 2))" "$padding" "$title" "$(((length - ${#title}) / 2))" "$padding"
    printf "\n First Ligand Clean-up (crystallography additives, solvents, etc) & Pocket detection \n"

	python $rootdir/LigExtract/bin/find_ligands.py --pdbPath $d --ligandsPath "$d"_LIGS --dist 6 --uniprot2pdbFile "$d"_pdb_uniprot_filteredlist.txt --keeprepeats y > find_ligands.log
	

    #############################################   MODULE 4  ##################################################
	title=" MODULE 4: Ligands clustering "
    printf "%.*s %s %.*s\n" "$(((length - 1 - ${#title}) / 2))" "$padding" "$title" "$(((length - ${#title}) / 2))" "$padding"

	rm -f *_pockets_hierarch-clusters.txt
	python $rootdir/LigExtract/bin/cluster_ligands_hierarchical.py --pdbPath $d --ligandsPath "$d"_LIGS --prdCif $rootdir/data/prd-all.cif --uniprot2pdbFile "$d"_pdb_uniprot_filteredlist.txt > cluster_ligands.log
	
fi


# cleaup
if [ $cleanup == "clean" ]; then
    rm *.cif
fi


#!/bin/bash

# ligextract.sh -d myproteins -r 3.5 -o cluster -c nan
# ligextract.sh -d myproteins -r 3.5 -o filter -c clean

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

# clear out all log files before a new run
rm -f *.log

printf "\n################ Processing [${d}] ################\n" 

######################################################################################################################
###################################################### MODULE 1 ######################################################

##download all PDBs
python ~/LigExtract/bin/getPdbsFromUniprot.py --outputDir $d --uniprots "$d"_uniprot_list.txt --allPdbs ~/LigExtract/data/allpdbs.txt --maxResol $res
if [[ $? = 123 ]]; then printf "\nAbort workflow! inspect the logs to fix any issues."; exit 1 ; fi
rm -f pdbs2download.txt
python ~/LigExtract/bin/downloadPdbs.py --outputDir $d 
#if [[ $? = 123 ]]; then echo "Abort workflow! inspect the logs to fix any issues."; exit 1 ; fi
if [[ -f pdbs2download.txt ]]; then 
	  bash ~/LigExtract/bin/pdb_batch_download.sh -f pdbs2download.txt -p > pdb_download.log
fi

# catch failed download and update *_pdb_uniprot_filteredlist.txt
python ~/LigExtract/bin/bypass_missing_pdbs.py

### RENAMING PDBs to lowercase: no longer needed
#if compgen -G "*.gz" > /dev/null; then
#    gunzip *.gz
#    for file in *.pdb; do mv -- $file $(echo $file | tr 'A-Z' 'a-z'); done
#    mv *.pdb $d/.
#fi


> "$d"_process_uniprot_chains.err
> "$d"_process_uniprot_chains.txt
python ~/LigExtract/bin/process_chains.py --pdbspath $d
if [[ $? = 123 ]]; then printf "\nAbort workflow! inspect the logs to fix any issues."; exit 1 ; fi

## check if some pdbs failed
redo_pdbs=`grep "Request to Uniprot failed" "$d"_process_uniprot_chains.err | cut -f1 -d":"`
echo `grep "Request to Uniprot failed" "$d"_process_uniprot_chains.err | wc -l` "failed PDBs."
sed -i '/Request to Uniprot failed after/d' "$d"_process_uniprot_chains.err
if [[ $redo_pdbs != "" ]]; then 
	echo "LigExtract will repeat chain processing for failed PDBs"
	python ~/LigExtract/bin/process_chains.py --pdbspath $d --pdbsToRedo $redo_pdbs; 
fi 


######################################################################################################################
###################################################### MODULE 2 ######################################################
printf "\n---------------------------  Extracting all possible ligands from PDBs  ---------------------------\n\n"
rm -fR "$d"_LIGS
python ~/LigExtract/bin/extract_ligands.py --pdbPath $d --outputPath "$d"_LIGS --uniprot2pdbFile "$d"_pdb_uniprot_filteredlist.txt > ligand_extraction.log
if [[ $? = 123 ]]; then printf "\nAbort workflow! inspect the logs to fix any issues."; exit 1 ; fi
ls "$d"_LIGS > rawlist_extraction.txt



if [ $filter_option == "filter" ]; then
	#############################################   MODULE 3  ##################################################
	printf "\n------  First Ligand Clean-up (crystallography additives, solvents, etc) & Pocket detection  ------\n\n"
	python ~/LigExtract/bin/find_ligands.py --pdbPath $d --ligandsPath "$d"_LIGS --dist 6 --uniprot2pdbFile "$d"_pdb_uniprot_filteredlist.txt --keeprepeats n > find_ligands.log
	if [[ $? = 123 ]]; then printf "\nAbort workflow! inspect the logs to fix any issues."; exit 1 ; fi

  #############################################   MODULE 4  ##################################################
	printf "\n------------------------------------  Final Ligand Selection  ------------------------------------\n\n"
	python ~/LigExtract/bin/filter_ligands.py --pdbPath $d --ligandsPath "$d"_LIGS --prdCif ~/data/prd-all.cif > filter_ligands.log
	if [[ $? = 123 ]]; then printf "\nAbort workflow! inspect the logs to fix any issues."; exit 1 ; fi
	python ~/LigExtract/bin/assemble_finalreport.py --pdbPath $d --ligandsPath "$d"_LIGS
fi



if [ $filter_option == "cluster" ]; then
	#############################################   MODULE 3  ##################################################
	printf "\n------  First Ligand Clean-up (crystallography additives, solvents, etc) & Pocket detection  ------\n\n"
	python ~/LigExtract/bin/find_ligands.py --pdbPath $d --ligandsPath "$d"_LIGS --dist 6 --uniprot2pdbFile "$d"_pdb_uniprot_filteredlist.txt --keeprepeats y > find_ligands.log
	if [[ $? = 123 ]]; then printf "\nAbort workflow! inspect the logs to fix any issues."; exit 1 ; fi


  #############################################   MODULE 4  ##################################################
	printf "\n---------------------------------------  Pockets clustering  ---------------------------------------\n\n"
	rm -f *_pockets_hierarch-clusters.txt
	python ~/LigExtract/bin/cluster_ligands_hierarchical.py --pdbPath $d --ligandsPath "$d"_LIGS --prdCif ~/data/prd-all.cif --uniprot2pdbFile "$d"_pdb_uniprot_filteredlist.txt > cluster_ligands.log
	if [[ $? = 123 ]]; then printf "\nAbort workflow! inspect the logs to fix any issues."; exit 1 ; fi
fi


# cleaup
if [ $cleanup == "clean" ]; then
    rm *.cif
fi


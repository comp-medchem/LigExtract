import pandas as pd
import numpy as np
from glob import glob

pdbslist = glob("*_pdb_uniprot_filteredlist.txt")[0]

failed = open("pdb_download.log").readlines()
failed = [ln.split("/")[-1].split(".")[0].upper() for ln in failed if ln.startswith("Failed to download")]

allpdbs = pd.read_csv(pdbslist, sep="\t")

lines2rm = np.in1d(allpdbs.pdb, failed)
if lines2rm.sum()>0:
	print(f"Removing ** {lines2rm.sum()} ** lines from {pdbslist} due to missing PDF-format file")
	allpdbs = allpdbs[~lines2rm]
	allpdbs.to_csv(pdbslist, sep="\t", index=False)
	outfile = open("failed_pdbs_missingFile.log", "w")
	outfile.write("\n".join(failed))
	outfile.close()
	print("\nthe full list of failed PDBs is available in failed_pdbs_missingFile.log\n\n")
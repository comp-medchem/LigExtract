# python3

#
# --- DESCRIPTION ----
# This script converts all 5-letter codes in cif files into the short letter code attributed by BeEM during conversion
# back into PDB
# This is a temporary solution, because ideally the 5-letter code should be handled as is during the whole workflow

from glob import glob
from biopandas.pdb import PandasPdb
import pandas as pd
import os

for pdb_convert in glob("cifs/*ligand-id-mapping.tsv"):
    pdb = pdb_convert.split("-")[0].split("/")[1]
    #check that there are no idmapping-to-pdb mismatch
    if len(glob(f"cifs/{pdb}.cif"))==0: 
        os.remove(f"cifs/{pdb}-ligand-id-mapping.tsv")
        continue
    
    # get name conversions
    
    pdb_lig_convert = [ln.strip().split("\t") for ln in open(pdb_convert).readlines()]
    pdb_lig_convert = pd.DataFrame(pdb_lig_convert[1:], columns=pdb_lig_convert[0])
    
    # replace all long codes with the respective short code in the cif file
    # pdb will have the extra space for 2-letter codes, the mmcif can handle 2 and 3-letter codes simultaneously
    # all "long" codes in cif will be translated into the "short" codes
    
    cif = open(f'cifs/{pdb}.cif').read()
    for short,long in pdb_lig_convert.values:
        #print(short,log)
        cif = cif.replace(long, short)
    ciffile = open(f'cifs/{pdb}.cif', "w")
    ciffile.write(cif)
    ciffile.close()


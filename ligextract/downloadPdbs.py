import sys, os
import pandas as pd
import numpy as np
from progress.bar import ChargingBar, Bar
from time import sleep
import requests
import subprocess
import json
import argparse

parser = argparse.ArgumentParser(description='Download all PDBs corresponding to a list of Uniprot IDs')
parser.add_argument('--outputDir', type=str, required=True, dest="targetdir", help='name of the directory that will be created/used to receive the PDB files')
args = parser.parse_args()

targetdir = args.targetdir

uniprot_pdb_dict = pd.read_csv(f'{targetdir}_pdb_uniprot_filteredlist.txt', sep="\t")

try: 
    os.mkdir(targetdir)
except FileExistsError:
    None

pdbs2download = uniprot_pdb_dict["pdb"].unique()
bar = Bar('Fetching PDB structures...', max=len(pdbs2download))
for pdb in pdbs2download:
    bar.next()
    if f'{pdb.lower()}.pdb' in os.listdir(targetdir):
        continue
    sleep(5)
    try: o=subprocess.check_output(f"wget https://files.rcsb.org/download/{pdb.lower()}.pdb --quiet -O {targetdir}/{pdb.lower()}.pdb", shell=True)
    except subprocess.CalledProcessError:
        print(f"failed {pdb}")

sys.exit(f"\n\n**** Finished collecting {len(pdbs2download)} PDBs for the user-provided Uniprot IDs!")

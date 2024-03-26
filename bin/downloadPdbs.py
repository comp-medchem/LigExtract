import sys, os
import pandas as pd
import numpy as np
from progress.bar import ChargingBar, Bar
from time import sleep
import requests
import subprocess
import json
import argparse
home = os.path.realpath(__file__)
home = home.split("/LigExtract")[0]

parser = argparse.ArgumentParser(description='Download all PDBs corresponding to a list of Uniprot IDs')
parser.add_argument('--outputDir', type=str, required=True, dest="targetdir", help='name of the directory that will be created/used to receive the PDB files')
args = parser.parse_args()

targetdir = args.targetdir

uniprot_pdb_dict = pd.read_csv(f'{targetdir}_pdb_uniprot_filteredlist.txt', sep="\t")

try: 
    os.mkdir(targetdir)
except FileExistsError:
    None

pdbs2download = [x.lower() for x in uniprot_pdb_dict["pdb"].unique()]


# check is pdb is already downloaded 
already_downloaded = [x.split(".")[0] for x in os.listdir(home)]

pdbs2download = np.setdiff1d(pdbs2download, already_downloaded)

if len(pdbs2download)>0:
    print(len(pdbs2download), "PDBs queued to download")
    outfile = open("pdbs2download.txt", "w")
    outfile.write(",".join(pdbs2download))
    outfile.close()
else:
    print("All PDBs already downloaded")


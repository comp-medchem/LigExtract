import sys, os
import pandas as pd
import numpy as np
from progress.bar import ChargingBar, Bar
from time import sleep
import requests
from requests.adapters import HTTPAdapter, Retry
from xml.etree import ElementTree
from urllib.parse import urlparse, parse_qs, urlencode
import subprocess
import json
import argparse
import urllib
from pathlib import Path
import gzip
#home = str(Path.home())
home = os.path.realpath(__file__)
home = home.split("/LigExtract")[0]
uniprot_fun = f'{home}/LigExtract/bin/'
sys.path.insert(1, uniprot_fun)
from uniprot_map_api import *

parser = argparse.ArgumentParser(description='Find all PDB IDs corresponding to a list of Uniprot IDs using the RCSB PDB Search API')
parser.add_argument('--outputDir', type=str, required=True, dest="targetdir",
    help='name of the directory that will be created/used to receive the PDB files')
parser.add_argument('--uniprots', type=str, required=True, dest="uniprot_lst", help='File containing just a list of Uniprot IDs of interest')
parser.add_argument('--allPdbs', type=str, required=True, dest="all_pdbs", help='File downloaded from ftp.wwpdb.org/pub/pdb/derived_data/index listing all PDBs.')
parser.add_argument('--maxResol', type=float, required=True, dest="resolution_lim", help='maximum resolution value allowed to download a PDB')
parser.add_argument('--pdbFilter', type=str, required=False, dest="pdbFilter", help='list of PDBs to consider, otherwise all PDBs mapped to the uniprotIDs will be considered.')

args = parser.parse_args()

targetdir = args.targetdir
uniprot_lst = args.uniprot_lst
all_pdbs_file = args.all_pdbs
resolution_lim = args.resolution_lim
pdbFilter = args.pdbFilter



length = 90
pad_char = '-'

title = "Collecting PDBs from Uniprot IDs provided"
padding_total = length - len(title) - 2
print(f"\n{pad_char * (padding_total // 2)} {title} {pad_char * (padding_total - (padding_total // 2))}\n")



# Uniprot mapping no longer retrieves anything (excludes 5-letter codes in ligands)
#def fastuniprot2pdb(uniprotLstFile, allpdbsTabl, resolution_limit):

# This is the function to use now
def getuniprot2pdb(uniprotLstFile, allpdbsTabl, resolution_limit):
    response = []
    line_cnt= 0
    #['PDB', 'CHAIN', 'SP_PRIMARY', 'RES_BEG', 'RES_END', 'PDB_BEG', 'PDB_END', 'SP_BEG', 'SP_END']
    with gzip.open(f"{home}/LigExtract/data/pdb_chain_uniprot.csv.gz") as f:
        for ln in f: # PDB  CHAIN   SP_PRIMARY
            line_cnt+=1
            ln=ln.decode("utf-8").strip().split(",")
            if line_cnt==2 and ln[2]!= 'SP_PRIMARY':
                sys.exit("file not structured as originally expected")
            if line_cnt<2: continue
            if ln[2] in uniprot_lst:
                response.append([ln[2], ln[0].upper()])
    response = pd.DataFrame(response, columns=["From", "To"])
    print(f"{len(response)} retrieved PDBs")
    response = response[np.isin(response.To, allpdbsTabl.query(f"RESOLUTION < {resolution_limit}").pdb)]
    print(f"{len(response)} PDBs under the set resolution")
    response = [list(x) for x in response.values]
    return(response)



############ Map Uniprot to PDB ###########

all_pdbs = [ln.strip() for ln in open(all_pdbs_file).readlines()]
all_pdbs_head = all_pdbs[0].split(", ")
all_pdbs = pd.DataFrame([ln.split("\t") for ln in all_pdbs[2:]], columns = all_pdbs_head)
all_pdbs = all_pdbs.rename(columns={"IDCODE":"pdb"})[["pdb", "HEADER","RESOLUTION","EXPERIMENT TYPE (IF NOT X-RAY)"]]
all_pdbs = all_pdbs.query("RESOLUTION != 'NOT'")
resol_lst = [np.nanmax(np.array(x.replace("NOT","nan").strip().split(",")).astype(float)) if "," in x else x for x in all_pdbs.RESOLUTION]

filter_lines = []
for val in all_pdbs.RESOLUTION:
    try:
        float(val)
        filter_lines.append(True)
    except ValueError:
        filter_lines.append(False)

all_pdbs = all_pdbs[filter_lines]
resol_lst = np.array(resol_lst)[filter_lines]

all_pdbs.loc[:,"RESOLUTION"] = np.array(resol_lst).astype(float)

uniprot_lst = np.unique([x.strip() for x in open(uniprot_lst).readlines()])
print(f'There are {len(uniprot_lst)} unique Uniprot IDs that will be processed.\n')


bar = Bar('Fetching PDB codes from Uniprot IDs... ', max=len(uniprot_lst))

uniprot_pdb_dict = getuniprot2pdb(uniprot_lst, all_pdbs, resolution_lim)

uniprot_pdb_dict = pd.DataFrame(uniprot_pdb_dict, columns = ["uniprot","pdb"])
missing = np.setdiff1d(uniprot_pdb_dict.pdb.values, all_pdbs.pdb.values)
if len(missing)>0:
    sys.stderr.write("\nThere are pdbs associated with some of the input Uniprot IDs which are missing in allpdbs.txt. Try to re-download a more recent version of this file.\n")
    sys.stderr.write(f"\t{','.join(missing)}\n\n")
    sys.exit(123)
uniprot_pdb_dict = uniprot_pdb_dict.merge(all_pdbs, on="pdb")
uniprot_pdb_dict = uniprot_pdb_dict.drop_duplicates()
print("\nFinished gathering PDB codes.\n")

############ Validate Experiments ###########

unique_experiment_counts = uniprot_pdb_dict.groupby("EXPERIMENT TYPE (IF NOT X-RAY)")["pdb"].nunique().reset_index().values
accepted_experiment_types = []
for param,n in unique_experiment_counts:
    txt = input(f"   Do you accept '{param.upper()}' entries (with {n} PDBs)? (type y/n) ")
    txt = txt.strip()
    if txt=="y": accepted_experiment_types.append(param)
    if txt not in ["n","y"]:
        print(f"Please type either 'y' or 'n' !")
        txt = input(f"   Do you accept '{param.upper()}' entries (used in {n} PDBs)? (type y/n) ")
        txt = txt.strip()
        if txt=="y": accepted_experiment_types.append(param)
        if txt not in ["n","y"]:
            print(f"Aborted due to invalid choice '{txt}'")
            sys.exit(123)


uniprot_pdb_dict = uniprot_pdb_dict[np.isin(uniprot_pdb_dict["EXPERIMENT TYPE (IF NOT X-RAY)"].values, accepted_experiment_types)]

# optional filter with the manual list
if args.pdbFilter is not None:          # user supplied a file
    pdb_ids_manualinput = [x.strip().upper() for x in open(args.pdbFilter).readlines()]
    uniprot_pdb_dict = uniprot_pdb_dict[uniprot_pdb_dict.pdb.isin(pdb_ids_manualinput)]


uniprot_pdb_dict.to_csv(f'{targetdir}_pdb_uniprot_filteredlist.txt', sep="\t", index=False)

accepted_experiment_types = [f"'{x}'" for x in accepted_experiment_types]
print("\nAccepted experiment types:",", ".join(accepted_experiment_types), "\n")


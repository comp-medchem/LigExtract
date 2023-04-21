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
home = str(Path.home())
uniprot_fun = f'{home}/LigExtract/bin/'
sys.path.insert(1, uniprot_fun)
from uniprot_map_api import *

parser = argparse.ArgumentParser(description='Find all PDB IDs corresponding to a list of Uniprot IDs using the RCSB PDB Search API')
parser.add_argument('--outputDir', type=str, required=True, dest="targetdir",
    help='name of the directory that will be created/used to receive the PDB files')
parser.add_argument('--uniprots', type=str, required=True, dest="uniprot_lst", help='File containing just a list of Uniprot IDs of interest')
parser.add_argument('--allPdbs', type=str, required=True, dest="all_pdbs", help='File downloaded from ftp.wwpdb.org/pub/pdb/derived_data/index listing all PDBs.')
parser.add_argument('--maxResol', type=float, required=True, dest="resolution_lim", help='maximum resolution value allowed to download a PDB')

args = parser.parse_args()

targetdir = args.targetdir
uniprot_lst = args.uniprot_lst
all_pdbs_file = args.all_pdbs
resolution_lim = args.resolution_lim


print("\n------------------  Collecting PDBs from Uniprot IDs provided  ------------------\n")



def fastuniprot2pdb(uniprotLstFile, allpdbsTabl, resolution_limit):
    for ntry in range(5):
        try:
            job_id = submit_id_mapping(from_db="UniProtKB_AC-ID", to_db="PDB", ids=uniprotLstFile)
            break
        except requests.HTTPError:
            sleep(3)
            print("retry uniprot mapping...")
            continue
    if check_id_mapping_results_ready(job_id):
        link = get_id_mapping_results_link(job_id)
        response = get_id_mapping_results_search(link)
    else:
        print("There was an issue with the server. Failed after 5 requests.")
    response = [[x["from"], x["to"]] for x in response["results"]]
    if len(response)==0:
       sys.exit("No retrieved PDBs")
    response = pd.DataFrame(response, columns=["From", "To"])
    print(f"{len(response)} retrieved PDBs")
    response = response[np.in1d(response.To, allpdbsTabl.query(f"RESOLUTION < {resolution_limit}").pdb)]
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

uniprot_pdb_dict = fastuniprot2pdb(uniprot_lst, all_pdbs, resolution_lim)

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


uniprot_pdb_dict = uniprot_pdb_dict[np.in1d(uniprot_pdb_dict["EXPERIMENT TYPE (IF NOT X-RAY)"].values, accepted_experiment_types)]
uniprot_pdb_dict.to_csv(f'{targetdir}_pdb_uniprot_filteredlist.txt', sep="\t", index=False)

accepted_experiment_types = [f"'{x}'" for x in accepted_experiment_types]
print("\nAccepted experiment types:",", ".join(accepted_experiment_types), "\n")


# python ~/scripts/getPdbsFromUniprot.py --outputDir targetsofinterest --uniprots uniprot_list.txt --allPdbs allpdbs.txt --maxResol 2.5

import sys, os
import pandas as pd
import numpy as np
from progress.bar import ChargingBar, Bar
from time import sleep
import requests
import subprocess
import json
import argparse
import urllib


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

def uniprot2pdb(uniprotID, resolution_limit):
    sleep(5)
    rowlim = 100
    params = {"query": 
    {"type": "group","logical_operator": "and","nodes": 
    [
    {"type": "terminal","service": "text", "parameters": {"operator": "exact_match", "value": f"{uniprotID.upper()}" , "attribute":"rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession"}},
    {"type": "terminal", "service": "text","parameters": { "operator": "exact_match", "value": "UniProt","attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_name"}},
    {"type": "terminal","service": "text","parameters": {"operator": "less_or_equal","value": resolution_limit , "attribute":"rcsb_entry_info.resolution_combined"}}
    ]},
    "request_options": {"pager": {"start": 0,"rows": rowlim}},
    "return_type": "polymer_entity"}
    result = requests.get("https://search.rcsb.org/rcsbsearch/v1/query", {"json": json.dumps(params, separators=(",", ":"))})
    res_stat = result.status_code
    result = result.content.decode("utf-8")
    if result == "":
        if res_stat not in [200,204]: pdbmaplog.write(f"{uniprotID}: unssuccessful retrieval of uniprot ID\n")  
        return []
    result = eval(result)
    cnt = result["total_count"]
    if cnt > rowlim: 
        params["request_options"] = {'pager': {'start': 0, 'rows': cnt}}
        result = requests.get("https://search.rcsb.org/rcsbsearch/v1/query", {"json": json.dumps(params, separators=(",", ":"))})
        result = eval(result.content.decode("utf-8"))
    result=[x['identifier'] for x in result['result_set']]
    result = [x.split("_")[0] for x in result]
    return result


def fastuniprot2pdb(uniprotLstFile, allpdbsTabl, resolution_limit):
    url = 'https://www.uniprot.org/uploadlists/'
    unip_q = " ".join(uniprotLstFile)
    params = {
    'from': 'ACC+ID',
    'to': 'PDB_ID',
    'format': 'tab',
    'query': f'{unip_q}'
    }
    
    data = urllib.parse.urlencode(params)
    data = data.encode('utf-8')
    req = urllib.request.Request(url, data)
    with urllib.request.urlopen(req) as f:
       response = f.read()
    
    response = [x.split("\t") for x in response.decode('utf-8').split("\n")]
    response = response[:-1] #trim empty last line
    
    if len(response)==1:
       sys.exit("No retrieved PDBs")
    
    response = pd.DataFrame(response[1:], columns=response[0])
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
all_pdbs.loc[:,"RESOLUTION"] = np.array(resol_lst).astype(float)

uniprot_lst = np.unique([x.strip() for x in open(uniprot_lst).readlines()])
print(f'There are {len(uniprot_lst)} unique Uniprot IDs that will be processed.\n')

pdbmaplog = open("pdb_mapping.log","w")

bar = Bar('Fetching PDB codes from Uniprot IDs... ', max=len(uniprot_lst))

uniprot_pdb_dict = fastuniprot2pdb(uniprot_lst, all_pdbs, resolution_lim)

pdbmaplog.close()

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


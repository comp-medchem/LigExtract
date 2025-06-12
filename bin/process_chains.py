import sys,os
import numpy as np
import urllib
from urllib import parse,request
from time import sleep
from progress.bar import ChargingBar, Bar
import argparse
import pandas as pd
import gzip
from pathlib import Path
from glob import glob
from pdbecif.mmcif_io import CifFileReader
#home = str(Path.home())
home = os.path.realpath(__file__)
home = home.split("/LigExtract")[0]
sys.path.append(f'{home}/LigExtract/bin')
from uniprot_map_api import *

parser = argparse.ArgumentParser(description='Process all PDBs to annotate all chains with a possible protein identifier.')
parser.add_argument('--pdbspath', type=str, required=True, dest="pdbpath", help='Path to directory containing all PDBs to process')
parser.add_argument('--pdbsToRedo', type=str, required=False, nargs='+', default=["None"], dest="pdbs_redo", help='a list of space-separated pdb files to re-run (only used in re-runs of this script)')
args = parser.parse_args() #   nargs='+'

pdbpath = args.pdbpath
pdbs_redo = args.pdbs_redo
if pdbpath[-1]=="/": pdbpath = pdbpath[:-1]

# 1ai8 has a chain with P28501 in the CIF file, which is not the primary uniprot P01050
# pdb_chain_uniprot.csv.gz provides the updated primary uniprot

length = 90; pad_char = '-'

title="Processing PDBs to gather chain information"
padding_total = length - len(title) - 2
print(f"\n{pad_char * (padding_total // 2)} {title} {pad_char * (padding_total - (padding_total // 2))}\n")


pdbs = [x for x in os.listdir(pdbpath) if x.endswith(".pdb")]

# old2new
testprotein_lst = []

for pdbname in pdbs:
    pdbcode = pdbname.split(".")[0]
    ciffile = f"cifs/{pdbcode}.cif"
    if glob(ciffile) == 0:
        print("file does not exit:", ciffile)
    data = CifFileReader().read(ciffile)
    data = data[pdbcode.upper()]
    if "_struct_ref_seq" in data:
        struct_ref = pd.DataFrame.from_dict(data["_struct_ref"], orient="index").T
        struct_ref = struct_ref[["id", "db_name"]]
        struct_ref = struct_ref.query("db_name == 'UNP'")
        struct_ref = struct_ref.rename(columns = {"id":"ref_id"})
        testprotein = pd.DataFrame.from_dict(data["_struct_ref_seq"], orient="index").T
        # ref_id is a pointer to _struct_ref.id
        testprotein = testprotein[["ref_id", "pdbx_PDB_id_code", "pdbx_db_accession", "pdbx_strand_id"]]
        testprotein = testprotein.merge(struct_ref, on="ref_id")
        testprotein_lst.append(testprotein)

testprotein_lst = pd.concat(testprotein_lst)
testprotein_lst = testprotein_lst[testprotein_lst.pdbx_PDB_id_code != testprotein_lst.pdbx_db_accession]
l = testprotein_lst.pdbx_PDB_id_code.str.lower()
testprotein_lst.loc[:,"pdbx_PDB_id_code"] = l

new_uniprots = []
list_pdbs = [pdbname.split(".")[0].lower() for pdbname in pdbs]
with gzip.open(f"{home}/LigExtract/data/pdb_chain_uniprot.csv.gz") as f:
    for ln in f: # PDB  CHAIN   SP_PRIMARY
        ln=ln.decode("utf-8").strip().split(",")
        if ln[0] in list_pdbs:
            new_uniprots.append(ln[0:3])

new_uniprots = pd.DataFrame(new_uniprots, columns = ["pdbs", "chainName", "uniprot"])
new_uniprots = new_uniprots[[x.startswith("NOR")==False for x in new_uniprots.uniprot.values]]
new_uniprots = new_uniprots.drop_duplicates()

old2new_lst = []
chains2solve = []
for pdb, ch in testprotein_lst[["pdbx_PDB_id_code", "pdbx_strand_id"]].drop_duplicates().values:
    res = new_uniprots.query(f"pdbs == '{pdb}' and chainName == '{ch}'").drop_duplicates()[['chainName', 'uniprot']]
    #comes from the raw cif file
    q = testprotein_lst.query(f"pdbx_PDB_id_code == '{pdb}' and pdbx_strand_id == '{ch}'").drop_duplicates()[['pdbx_db_accession', 'pdbx_strand_id']]
    if len(res) == 0:
        # there is nothing to do as there are no uniprot codes for this chain in the SIFTS file
        q.loc[:,"pdb"] = pdb
        q.loc[:,"uniprot"] = q.pdbx_db_accession.values
        q = q[['uniprot', "pdbx_db_accession", 'pdbx_strand_id', "pdb"]]
        q=q.rename(columns={"pdbx_strand_id":"chainName"})
        old2new_lst.append(q)
        continue
    if len(q) > 1:
        # More than one uniprot is mapped to a chain, currently in the file. We cannot assert old-2-new dict
        # probably a chimera
        for q_i in q.pdbx_db_accession.unique(): chains2solve.append([pdb,ch, q_i])
        continue

    uniprot_pairs = res.merge(q, right_on="pdbx_strand_id", left_on="chainName").drop_duplicates()[['uniprot', "pdbx_db_accession", 'chainName']]
    uniprot_pairs.loc[:,"pdb"] = pdb
    uniprot_pairs = uniprot_pairs[['uniprot', "pdbx_db_accession", 'chainName', "pdb"]]
    old2new_lst.append(uniprot_pairs)

old2new_lst = pd.concat(old2new_lst)
old2new_lst.columns = ["updated_uniprot", "uniprot_in_pdbfile", "chain", "pdb"]
old2new_lst = old2new_lst.drop_duplicates()
old2new_lst = old2new_lst[[ "pdb", "chain", "uniprot_in_pdbfile", "updated_uniprot"]]


if len(chains2solve)>0:
    uniprots2solve = np.unique([x[-1] for x in chains2solve])
    job_id = submit_id_mapping(
        from_db="UniProtKB_AC-ID", to_db="UniProtKB", ids=list(uniprots2solve)
    )
    if check_id_mapping_results_ready(job_id):
        link = get_id_mapping_results_link(job_id)
        results = get_id_mapping_results_search(link)

    results = {x["from"]:x["to"]["primaryAccession"] for x in results["results"]}

    #chains2solve.append([pdb,ch, q_i])
    old2new_lst2 = []
    for pdb, ch, uniprot in chains2solve:
        new_unp = results[uniprot]
        old2new_lst2.append([pdb, ch, uniprot, new_unp])
    old2new_lst2 = pd.DataFrame(old2new_lst2, columns = [ "pdb", "chain", "uniprot_in_pdbfile", "updated_uniprot"])
    old2new_lst = pd.concat([old2new_lst, old2new_lst2])


old2new_lst.to_csv(f"{pdbpath.split('/')[-1]}_process_uniprot_chains.txt", sep="\t", index=False)


print("\n\nFinished.\n\n")

# python3
#
# --- DESCRIPTION ----
# This script runs through a directory containing protein PDB files and searches for ligands.
# In order to avoid removing small ligands by accident this script only filters out molecules that
# have 3 atoms, maximum. ligands that are larger than this but exist in hundreds of other PDBs will
# be extracted into a file but will raise a warning shown in the LOG file.
# These cases should be submitted to manual inspection as molecules such as ATP exist in thousands
# of PDBs but may be the main ligand in some particular cases.
#
# This script creates a results folder if it does not exist already. Information of the run is stored in the .log file
#

import pandas as pd
import numpy as np
from biopandas.pdb import PandasPdb
import string
from copy import copy
import sys
import os
from pathlib import Path
import networkx as nx
import urllib
from time import sleep
from progress.bar import Bar
import argparse
import gzip
import subprocess
import shlex
import json
import requests
import re
from concurrent.futures import ProcessPoolExecutor, as_completed
from pdbecif.mmcif_io import CifFileReader
from glob import glob
from scipy.spatial.distance import euclidean
from itertools import product,combinations

#home = str(Path.home())
home = os.path.realpath(__file__)
home = home.split("/LigExtract")[0]

sys.path.append(f'{home}/LigExtract/bin')
from ligextract_utils import *

parser = argparse.ArgumentParser(description='Search and extract all *possible* ligands in each PDB. These will be submitted to further filtration.')
parser.add_argument('--pdbPath', type=str, required=True, dest="pdbpath", help='Path to directory containing all PDBs to process')
parser.add_argument('--outputPath', type=str, required=True, dest="outpath", help='Path to a directory that will created/re-used to put the extracted ligand files.')
parser.add_argument('--uniprot2pdbFile', type=str, required=True, dest="uniprot2pdb_file", help='File with Uniprot-to-Pdb mappings, ending with *pdb_uniprot_filteredlist.txt')
args = parser.parse_args() #                                 

pdbpath = args.pdbpath
outpath = args.outpath
uniprot2pdb_file = args.uniprot2pdb_file

if pdbpath[-1]=="/": pdbpath = pdbpath[:-1]
if outpath[-1]=="/": outpath = outpath[:-1]


length = 90; pad_char = '-'

title="Extracting all possible ligands from PDBs"
padding_total = length - len(title) - 2
print(f"\n{pad_char * (padding_total // 2)} {title} {pad_char * (padding_total - (padding_total // 2))}\n")


# First handle CIFs that are to large (remove now; to include later)
largecifs = [x.split("-")[0] for x in glob("*-chain-id-mapping.txt")]

for pdb in largecifs:
    for pdbsplit in glob(f'{pdbpath}/{pdb}*.pdb'):
        os.remove(pdbsplit)
    if len(glob(f"cifs/{pdb}.cif"))>0:
        os.remove(f"cifs/{pdb}.cif")

originpath = f'{os.getcwd()}/{pdbpath}'.replace("/./","/")
targetpath = f'{os.getcwd()}/{outpath}'.replace("/./","/")

if ("pdb_chain_uniprot.csv.gz" in os.listdir(f"{home}/LigExtract/data")) == False: 
    print("Expecting pdb_chain_uniprot.csv.gz to be in the 'data' directory but file cannot be found. Abort!")
    sys.exit(123)

whitelisted_ligs = [ln.strip() for ln in open(f"{home}/LigExtract/docs/whitelist_ligands.txt").readlines()]

uniprot2pdb = pd.read_csv(uniprot2pdb_file, sep="\t") 
missing_downloads = np.setdiff1d(uniprot2pdb.pdb.str.lower(), [x.split(".")[0] for x in os.listdir(pdbpath) if x.endswith(".pdb")])
extra_downloads = np.setdiff1d([x.split(".")[0] for x in os.listdir(pdbpath) if x.endswith(".pdb")], uniprot2pdb.pdb.str.lower())
pdbs = [x for x in os.listdir(pdbpath) if x.endswith(".pdb")]

if len(missing_downloads)>0:
    if len(missing_downloads)>0: sys.stderr.write(f"\nWARNING!!!  There are {len(missing_downloads)} pdbs missing ({' '.join(missing_downloads)}). Please re-run the mmCIF-to-Pdb conversion block\n")
    sys.exit(123)


if len(extra_downloads)>0:
    for pdbx in extra_downloads:
        for pdb_f in glob(f"{pdbpath}/{pdbx}*"): os.remove(pdb_f)
        for cif_f in glob(f"cifs/{pdbx}*"): os.remove(cif_f)


uniprot_pdb = pd.read_csv(f"{pdbpath.split('/')[-1]}_process_uniprot_chains.txt", sep="\t")
uniprot2pdb_secondary = copy(uniprot_pdb)
# pdb  uniprot_in_pdbfile  updated_uniprot

uniprot_pdb = uniprot_pdb[["uniprot_in_pdbfile","updated_uniprot"]].drop_duplicates()
uniprot_pdb["updated_uniprot"] = uniprot_pdb["updated_uniprot"]+","
uniprot_pdb = uniprot_pdb.groupby("uniprot_in_pdbfile")["updated_uniprot"].sum().reset_index()
uniprot_pdb["updated_uniprot"]=[x.split(",")[:-1] for x in uniprot_pdb.updated_uniprot]



try : os.listdir(outpath)
except FileNotFoundError: os.makedirs(outpath)


keep_others = ["HEADER","TITLE","LINK","CONECT", "MASTER","END"]

if len(pdbs)==0:
    print(f" There are no PDB files in the target folder {targetpath}.\n Make sure your PDB files have a '.pdb' extension.")
    sys.exit(123)


def detectBoundLigands(file1, file2):
    """
    detect minimum distance between two molecules (HETATM group) to indicate possible covalent bond; confirmed with PyMOL dist function
    """
    lig1 = PandasPdb().read_pdb(file1)
    lig2 = PandasPdb().read_pdb(file2)
    lig1_coor = lig1.df['HETATM'][["x_coord","y_coord","z_coord"]].values
    lig2_coor = lig2.df['HETATM'][["x_coord","y_coord","z_coord"]].values
    ed = []
    for coor1, coor2 in product(lig1_coor, lig2_coor):
        ed.append(euclidean(coor1, coor2))
    #resnum1 = lig1.df["HETATM"].residue_number.unique()[0]
    #resnum2 = lig2.df["HETATM"].residue_number.unique()[0]
    return(min(ed)) #1.7: print("might be covalently bound")



# After all cifs have been downloaded, screen them to grab PRD IDs
bar = Bar('Extract all chains with PRD code... ', max=len(pdbs))
all_prds_chains = []
for pdb in pdbs:
    bar.next()
    pdbcode = pdb.split(".")[0]
    ciffile = f"cifs/{pdbcode}.cif"
    cifdata = CifFileReader().read(ciffile)
    data = cifdata[pdbcode.upper()]
    if "_pdbx_molecule" in data:
        data = pd.DataFrame.from_dict(data["_pdbx_molecule"], orient="index").T
        data.loc[:,"pdb"] = pdbcode
        all_prds_chains.append(data)

bar.finish()

if len(all_prds_chains)>0:
    all_prds_chains = pd.concat(all_prds_chains)
else:
    all_prds_chains = pd.DataFrame([],columns = ["instance_id","prd_id","asym_id","pdb"])


# gather a new-2-old dictionary of chains
bar = Bar('Building chains dictionary... ', max=len(pdbs))
pdb_auth_chains_dict = {}

for pdb in pdbs:
    bar.next()
    pdbcode = pdb.split(".")[0]
    ciffile = f"cifs/{pdbcode.lower()}.cif"
    cifdata = CifFileReader().read(ciffile)
    data = cifdata[pdbcode.upper()]
     # translate new chain ID (asym_id) to original chain ID (author_asym_id)
    if "_atom_site" in data:
        data = pd.DataFrame.from_dict(data["_atom_site"], orient="index").T
        pdb_auth_chains = data[["auth_asym_id", "label_asym_id"]].drop_duplicates()
        pdb_auth_chains_dict[pdbcode] = {newlabel:authlabel for authlabel,newlabel in pdb_auth_chains[["auth_asym_id", "label_asym_id"]].values}
bar.finish()

# translate PRD's chains (reported as label_asym_id) to auth_asym_id
all_prds_chains_NEW = []
for instance_id,prd_id,asym_id,pdb in all_prds_chains.values:
    asym_id = pdb_auth_chains_dict[pdb][asym_id]
    all_prds_chains_NEW.append([instance_id,prd_id,asym_id,pdb ])

all_prds_chains_NEW = pd.DataFrame(all_prds_chains_NEW, columns = all_prds_chains.columns)
all_prds_chains_NEW.to_csv(f"{home}/LigExtract/data/prdID_chain_pdb.txt", sep="\t", index=False)


# update prd_to_pdb_IDs.txt
prd2pdb = pd.read_csv(f"{home}/LigExtract/data/prd_to_pdb_IDs.txt", sep="\t")
save_headers = prd2pdb.columns
prd2pdb = np.vstack([prd2pdb.values, all_prds_chains_NEW[["prd_id", "pdb"]].values])
prd2pdb = pd.DataFrame(prd2pdb, columns = save_headers)
prd2pdb.to_csv(f"{home}/LigExtract/data/prd_to_pdb_IDs.txt", sep="\t", index=False)


nonpolymer_cnt = {}




#for pdbname in pdbs:
def extractorSplit(pdbname):
    storeLog=[]
    print(f"#######  Finding ligands in  {pdbname.split('.')[0]}  #######")
    pdbcode = pdbname.upper().split(".")[0]
    pdb = open(f"{pdbpath}/{pdbname}").readlines()
    if len(pdb) == 0:
        storeLog.append(f"! {pdbname} has zero lines. Please verify the file is valid.")
        storeLog = "\n".join(storeLog)
        print(storeLog)
        return(pdbcode)
    
    try: uniprotid = uniprot2pdb.query(f"pdb == '{pdbcode.upper()}'").uniprot.values[0]
    except IndexError:
        storeLog.append(f"Unable to find {pdbcode.upper()} in {uniprot2pdb_file}. Make sure this pdb code is in the file, in all capitals.")
        sys.exit(123)
    # remove PDBs where the protein actually is not found within DBREF
    uniprotid_alternates = uniprot2pdb_secondary.query(f"pdb == '{pdbcode.lower()}'").uniprot_in_pdbfile.values

    chains_w_Uniprot = uniprot2pdb_secondary.query(f"pdb == '{pdbcode}'").chain.unique()
    
    if len(chains_w_Uniprot)==0 and uniprotid not in uniprot2pdb_secondary.updated_uniprot.values:
        # Finally only delete if this pdb is not being used by other Uniprot IDs
        other_uniprots = uniprot2pdb_secondary.query(f"pdb == '{pdbcode.lower()}'").updated_uniprot.values
        other_uniprots = np.setdiff1d(other_uniprots, uniprotid)
        if np.isin(uniprot2pdb.uniprot.unique(), other_uniprots).any() == False:
            storeLog.append(f"{pdbcode} does not have any chain for ANY of the desired uniprotIDs in the full collection of targets! REMOVED protein")
        os.remove(f'{pdbpath}/{pdbname}')
        storeLog = "\n".join(storeLog)
        print(storeLog)
        return(pdbcode)
    
    # initialise file to store as a log
    logfile = open(f'{outpath}/{pdbcode.lower()}_ligand_extraction.log',"w")
    logfile.write(f" ###### LOG OF LIGAND EXTRACTION FOR {pdbname} ######\n\n")
    logfile.write(f"results obtained from PDBs in {originpath}\n")
    logfile.write(f"results stored in {targetpath}\n\n")
    
    ciffile = f"cifs/{pdbcode.lower()}.cif"
    cifdata = CifFileReader().read(ciffile)
    data = cifdata[pdbcode.upper()]
    if "_pdbx_entity_nonpoly" in data:
        data = pd.DataFrame.from_dict(data["_pdbx_entity_nonpoly"], orient="index").T
        nonpoly = list(data.comp_id.values)
    else:
        nonpoly = []


    if "_struct_conn" in cifdata[pdbcode.upper()]:
        links = pd.DataFrame.from_dict(cifdata[pdbcode.upper()]["_struct_conn"], orient="index").T
        #dist = links.pdbx_dist_value.astype(float) < 2.15 # max covalent distance
        #links = links[dist]
        links = links.query("ptnr1_label_comp_id != 'HOH'")
        links = links.query("conn_type_id == 'covale'")
    else:
        links = pd.DataFrame([])
    
    # what chains map to uniprot
    isprotein = uniprot2pdb_secondary.query(f"pdb == '{pdbcode.lower()}'").chain.unique()
    if len(isprotein) == 0:
        # try to populate with chains directly from the SIFTS mapping
        isprotein = []
        with gzip.open(f"{home}/LigExtract/data/pdb_chain_uniprot.csv.gz") as f:
            for ln in f: # PDB  CHAIN   SP_PRIMARY
                ln=ln.decode("utf-8").strip().split(",")
                if ln[0] in [pdbcode.lower()]:
                    isprotein.append(ln[1])
        isprotein = np.unique(isprotein)

    # isprotein cannot be empty
    if len(isprotein) == 0:
        sys.stderr.write("Interrupt process: no proteins mapped to PDB [",pdbcode, "] which is unexpected")

    # get all chains
    cifdata = CifFileReader().read(ciffile)
    data = cifdata[pdbcode.upper()]
    data = pd.DataFrame.from_dict(data["_struct_ref_seq"], orient="index").T
    allchains = data.pdbx_strand_id.values
    notaprotein = np.setdiff1d(allchains, isprotein)

    # add found peptides to the chain lists
    #peptidechains = findpeptides(f"{pdbpath}/{pdbname}")
    peptidechains = findpeptide(pdbcode).split(",")
    if peptidechains == [""]: peptidechains = []
    norine = [x.startswith("NOR") and len(x) == 8 for x in data.pdbx_db_accession.values]
    peptidechains_norine = data[norine].pdbx_strand_id.values
    peptidechains = np.union1d(peptidechains, peptidechains_norine)
    storeLog.append(f"peptide chains for {pdbcode}: {','.join(peptidechains)}")
    isprotein = np.setdiff1d(isprotein, peptidechains)
    notaprotein = np.union1d(notaprotein, peptidechains)
    
    #seq = [x for x in pdb if x.startswith('SEQRES')]
    data = cifdata[pdbcode.upper()]
    chainEntityID_AuthChainName_dict = pd.DataFrame.from_dict(data["_atom_site"], orient="index").T
    chainEntityID_AuthChainName_dict = chainEntityID_AuthChainName_dict[["label_entity_id","auth_asym_id"]].drop_duplicates()
    #AuthChainName_chainEntityID_dict = {y:x for x,y in chainEntityID_AuthChainName_dict.values}
    chainEntityID_AuthChainName_dict = {x:chainEntityID_AuthChainName_dict.query(f"label_entity_id == '{x}'").auth_asym_id.values for x in chainEntityID_AuthChainName_dict.label_entity_id.unique()}

    data = cifdata[pdbcode.upper()]
    pdbseqs = []
    if "_entity_poly_seq" in data:
        pdbseqs = pd.DataFrame.from_dict(data["_entity_poly_seq"], orient="index").T
        pdbseqs = pdbseqs[["entity_id","mon_id"]]
        #pdbseqs = pdbseqs.replace(chainEntityID_AuthChainName_dict)
    #pd.DataFrame.from_dict(data["_pdbx_nonpoly_scheme"], orient="index").T
    
    # gather ligands; # ARG in 4kx9 is a ligand 
    # ligand among the protein sequence (e.g. 4H9M) cannot be considered; warning: a standard residue code might appear in HETNAM
    
    #ciffile = f"cifs/{pdbcode.lower()}.cif"
    #cifdata = CifFileReader().read(ciffile)
    if "_pdbx_entity_nonpoly" in data:
        ligands = pd.DataFrame.from_dict(data["_pdbx_entity_nonpoly"], orient="index").T #_chem_comp
        ligands = ligands.comp_id.values
        ligands = np.setdiff1d(ligands, "HOH")
    else:
        ligands = np.array([])

    ###########  EXTRACT LIGAND - METHOD 1  ###########
    # This method extracts ligands that may exist (partly or entirely) among the residues section.
    # This happens when the ligand is a peptide or has aminoacids in its structure.
    #print(notaprotein)
    if len(notaprotein)==0:
        ligres=[]
    else:
        ligres = []
        for entityID in pdbseqs.entity_id.unique():
            chains = chainEntityID_AuthChainName_dict[entityID]
            for chainName in chains:
                if chainName in notaprotein:
                    ligres.append(pdbseqs.query(f"entity_id == '{entityID}'").mon_id.unique())
        ligres = np.unique(np.hstack(ligres))
        
        for chain in notaprotein:
            lig = PandasPdb().read_pdb(f'{originpath}/{pdbname}')
            lig.df["HETATM"] = lig.df["HETATM"].query(f'chain_id=="{chain}"')
            lig.df["ATOM"] = lig.df["ATOM"].query(f'chain_id=="{chain}"')
            # clean from solvents; not exhaustive 
            lig.df["HETATM"] = lig.df["HETATM"].query(f'residue_name!="HOH"')
            res_atoms = lig.df["HETATM"].groupby(["residue_name","residue_number"])["atom_number"].nunique()
            res2keep = res_atoms[res_atoms >= 1].reset_index().residue_number # some residues are just the carbon
            
            lig.df["HETATM"] = lig.df["HETATM"][np.isin(lig.df["HETATM"].residue_number,res2keep)]

            res_atom = lig.df['ATOM'][["residue_name","chain_id","residue_number"]].drop_duplicates().values
            res_hetatm = lig.df['HETATM'][["residue_name","chain_id","residue_number"]].drop_duplicates().values
            
            if len(res_atom) > 0: # part/all of the ligand is in ATOM
                res_atom =[f"{a}.{b}.{c}" for a,b,c in res_atom]
                res_hetatm =[f"{a}.{b}.{c}" for a,b,c in res_hetatm]
                # see which res connect from HETATM to ATOM
                
                # At the moment this does not cover all exceptions as some cases fail to mention ATOM-HETATM links that exist (e.g. 1de7)
                all_linked_atms = np.unique(findAllLinks(res_atom, "ATOM", links) + findAllLinks(res_hetatm,"HETATM", links))
                if len(all_linked_atms)==0: continue
                storeLog.append(f"----------######### {all_linked_atms}")
                filtered_lig_het = [lig.df["HETATM"].query(f'residue_name=="{ln.split(".")[0]}" and chain_id=="{ln.split(".")[1]}" and residue_number=={ln.split(".")[2]}') for ln in all_linked_atms]
                filtered_lig_atm = [lig.df["ATOM"].query(f'residue_name=="{ln.split(".")[0]}" and chain_id=="{ln.split(".")[1]}" and residue_number=={ln.split(".")[2]}') for ln in all_linked_atms]
                
                lig.df['ATOM'] = pd.concat(filtered_lig_atm)
                lig.df['HETATM'] = pd.concat(filtered_lig_het)

                
                if (len(lig.df["ATOM"][["residue_name", "chain_id", "residue_number"]].drop_duplicates())+len(lig.df["HETATM"][["residue_name", "chain_id", "residue_number"]].drop_duplicates())) != len(all_linked_atms): sys.exit("all_linked_atoms is not mapping correctly to ATOM/HETATM.")
                lig.to_pdb(path=f'{outpath}/{pdbname.split(".")[0]}_lig_chain-{chain}.pdb', records=['ATOM', 'HETATM'], gz=False, append_newline=True)
                storeLog.append(f'lig_chain-{chain}')
                logfile.write(f"rebuilt ligand detected in chain {chain}\n")

            if len(res_atom) == 0: # ligand is entirely in HETATM but has multiple residues
                res_hetatm =[f"{a}.{b}.{c}" for a,b,c in res_hetatm]
                linked_hetatm = findAllLinks(res_hetatm, "HETATM", links)
                lig.df['ATOM'] = lig.df['ATOM'][lig.df['ATOM'].chain_id==chain]
                lig.df['HETATM'] = lig.df['HETATM'][lig.df['HETATM'].chain_id==chain]
                if len(linked_hetatm)>0:
                    lig.df['HETATM'] = lig.df['HETATM'][np.isin(lig.df['HETATM'].residue_number,linked_hetatm)]
                if len(linked_hetatm)==0:
                    storeLog.append("[LOG WARN] no linked hetatm residues - bypass")
                    # untie with HET
                    #het = [x.split() for x in pdb if x.startswith('HET ')]
                    #het = [x[1] for x in het if x[2]==chain]
                    #lig.df['HETATM'] = lig.df['HETATM'][np.isin(lig.df['HETATM'].residue_name,het)]
                lig.to_pdb(path=f'{outpath}/{pdbname.split(".")[0]}_lig_chain-{chain}.pdb', records=['ATOM', 'HETATM', "OTHERS"], gz=False, append_newline=True)
                storeLog.append(f'lig_chain-{chain}')
                logfile.write(f"rebuilt ligand detected in chain {chain}; exclusively found in HETATM records\n")
    
    
    ##########  EXTRACT LIGAND - METHOD 2  ###########
    # This method extracts ligands from the typical section where ligands are found
    #cmpds_all = [cmpd[12:15].strip() for cmpd in formulas]
    cifdata = CifFileReader().read(ciffile)
    seqs = pd.DataFrame.from_dict(cifdata[pdbcode.upper()]["_entity_poly"], orient="index").T 
    res_not_cmpds = [x for x in seqs.pdbx_seq_one_letter_code.values if "(" in x]
    if len(res_not_cmpds)>0:
        res_not_cmpds = np.hstack([re.findall(r"\(.*?\)", x) for x in res_not_cmpds])
        res_not_cmpds = [x[1:-1] for x in res_not_cmpds]
    else:
        res_not_cmpds = []

    lig_ids_in_bird = pd.read_csv(f"{home}/LigExtract/data/prd_to_pdb_ligIDs.txt", sep="\t").pdblig_ID.values
    cmpds_ok, cmpds_removed = countMolsAtoms(pdbcode)
    # save cmpds_removes according to BIRD
    cmpds_removed = np.setdiff1d(cmpds_removed, lig_ids_in_bird)

    if len(cmpds_removed)>0: storeLog.append(f"The following compounds do not comply wih the inclusion criteria of minimum 3 heavy atoms and maximum of 10 repeated molecules: {cmpds_removed}")
    cmpds = np.union1d(cmpds_ok, ligands)

    cmpds = np.setdiff1d(cmpds, cmpds_removed)
    cmpds = np.setdiff1d(cmpds, res_not_cmpds)


    if len(cmpds)>0:
        for cpd in cmpds:
            message = ""
            with open(f"{home}/LigExtract/data/all_pdbligs.txt") as f: 
                for line in f:
                    if line == "id\tcount\n": continue
                    pdbcnt = int(line.strip().split("\t")[1])
                    if pdbcnt>250 and line.split("\t")[0]==cpd and line.split("\t")[0] not in whitelisted_ligs:
                        # see if compound has been already tested for non-polymer count
                        if cpd in nonpolymer_cnt.keys():
                            pdbcnt = nonpolymer_cnt[cpd]
                            if pdbcnt>250 : message = f"(!) Found in {pdbcnt} PDBs! Might not be a ligand."
                        else:
                            # calculate real non-polymer count
                            params = json.loads(open(f"{home}/LigExtract/bin/pdb_api_counts.q").read().replace("LIGAND", cpd))
                            result = requests.get("https://search.rcsb.org/rcsbsearch/v2/query", {"json": json.dumps(params, separators=(",", ":"))})
                            res_stat = result.status_code
                            result = result.content.decode("utf-8")
                            if result == "" and res_stat not in [200,204]: #,204
                                storeLog.append(f"{cpd}: unssuccessful request for PDB count\n") 
                                message = f"(!) Found in {pdbcnt} PDBs! Might not be a ligand."
                            elif res_stat==204: # has no hits
                                pdbcnt = 0
                                message = f"(!) Found in {pdbcnt} PDBs! Might not be a ligand."
                                nonpolymer_cnt[cpd]=pdbcnt
                            else:
                                result = eval(result)
                                pdbcnt = result["total_count"]
                                #FINAL result
                                if pdbcnt>250 :
                                    message = f"(!) Found in {pdbcnt} PDBs! Might not be a ligand."
                                nonpolymer_cnt[cpd]=pdbcnt
            storeLog.append(cpd)
            links_warning = []
            for link in links:
                if f' {cpd} ' in link:
                    links_warning.append(link)
            
            if len(links_warning)>0:
                logfile.write(f"{cpd}   ligand connected to other molecules. {message}\n")
                storeLog.append(f"    WARNING!  Ligand {cpd} appears to be connected to other molecules (see .log file).\n")
                storeLog.append(f"    This means it (1) is just a portion of the ligand or (2) is covalently bound (which requires further struture correction to return both ligand and protein to original unbound structure)\n")
                logfile.write("\t"+"            ---- res 1 ----               ---- res 2 ----                 bond (A)\n")
                for w in links_warning:
                    logfile.write("\t"+w+"\n")
            else:
                logfile.write(f"{cpd}   ligand detected. {message}\n")
            
            lig = PandasPdb().read_pdb(f'{originpath}/{pdbname}')
            lig_instances = lig.df["HETATM"].query(f"residue_name=='{cpd}'")[["residue_name", "residue_number", "chain_id"]].drop_duplicates()
            #since this compound is a isolated entity, there should be no case where it is a multi-residue compound.
            # separate multiple residues
            for resname, resx, chain in lig_instances.values:
                lig = PandasPdb().read_pdb(f'{originpath}/{pdbname}')
                lig.df["HETATM"] = lig.df["HETATM"].query(f'residue_name=="{cpd}" and chain_id=="{chain}" and residue_number=={resx}')
                lig.df["ATOM"] = lig.df["ATOM"].query(f'residue_name=="NOTHING"') # should not exist here
                lig.to_pdb(path=f'{outpath}/{pdbname.split(".")[0]}_chain-{chain}_lig-{cpd}-{resx}.pdb', records=['ATOM', 'HETATM', "OTHERS"], 
                                gz=False, append_newline=True)
                storeLog.append(f'{pdbname.split(".")[0]}_chain-{chain}_lig-{cpd}-{resx}') #cases like 8uv1 ligands will start replacing each other.


    # compare all ligands in cmpds, if links is empty. Rebuilding from links will happen in next module
    lig_files = glob(f'{outpath}/{pdbcode.lower()}*_lig-*.pdb')
    if len(links)> 0: 
        mismatched_chains_links = (links.ptnr1_label_asym_id != links.ptnr2_label_asym_id).any()
    else: mismatched_chains_links = False
    if (len(links)==0 or mismatched_chains_links==True) and len(lig_files)>1:# empty links or residues are of different chains
        G = nx.Graph()
        G.add_nodes_from(lig_files)
        for f1, f2 in combinations(lig_files, 2):
            mindist = detectBoundLigands(f1, f2)
            if mindist <= 1.7 and mindist>0.1: # avoid merging alternate states of the ligand (i.e. isomers in 1RGT for example)
                G.add_edge(f1, f2)

        # Find connected components
        connected_ligand_groups = list(nx.connected_components(G))

        for group in connected_ligand_groups:
            if len(group) == 1:
                continue  # single-res ligands, skip
            lig_merged = []
            atom_numbers = []
            for file in group:
                lig = PandasPdb().read_pdb(file)
                lig_merged.append(lig.df["HETATM"])
                atom_numbers.extend(lig.df["HETATM"].atom_number.values)
            pair_bound = pd.concat(lig_merged)

            # Save merged ligand
            chainName = pair_bound.chain_id.unique()[0]
            new_filename = f'{outpath}/{pdbcode.lower()}_lig_chain-h{chainName}.pdb'
            
            lig = PandasPdb().read_pdb(f'{originpath}/{pdbcode.lower()}.pdb')
            lig.df["HETATM"] = lig.df["HETATM"][lig.df["HETATM"].atom_number.isin(atom_numbers)]
            lig.to_pdb(path=new_filename, records=['HETATM'], gz=False, append_newline=True)

            # Remove original ligand files
            for file in group: os.remove(file)

            # log
            names = ', '.join([f.split('_lig-')[-1].split('.pdb')[0] for f in group])
            logfile.write(f"Connected ligands merged ({len(group)}): {names} (individual ligand files removed)\n")


    ########## METHOD 3: collecting oligosacharide chain ligands ##########

    # is glycosylation group - might be oligo or not
    glycosylation_chains = np.array([])
    cifdata = CifFileReader().read(ciffile)
    data = cifdata[pdbcode.upper()]
    if "_struct_conn" in data:
        data = pd.DataFrame.from_dict(data["_struct_conn"], orient="index").T
        if "pdbx_role" in data.columns:
            has_sugar = ["osylation" in x for x in data.pdbx_role]
            #partner chains in the connection
            glycosylation_chains = data[has_sugar][["ptnr1_auth_asym_id","ptnr2_auth_asym_id"]]
            glycosylation_chains = np.unique(glycosylation_chains.values.flatten())
            # remove chains that are the receiving proteins
            glycosylation_chains = np.setdiff1d(glycosylation_chains, isprotein)

    # is oligo ligand
    oligo_entities = np.array([])
    cifdata = CifFileReader().read(ciffile)
    data = cifdata[pdbcode.upper()]
    if "_pdbx_entity_branch" in data:
        oligo_entities = pd.DataFrame.from_dict(data["_pdbx_entity_branch"], orient="index").T
        oligo_entities = oligo_entities.query("type == 'oligosaccharide'").entity_id.values

    new2old_chainID = pdb_auth_chains_dict[pdbcode.lower()]

    oligo_entities_chains = np.array([])
    if "_pdbx_branch_scheme" in data:
        oligo_entities_chains = pd.DataFrame.from_dict(data["_pdbx_branch_scheme"], orient="index").T
        oligo_entities_chains = [new2old_chainID[x] for x in oligo_entities_chains.asym_id.values]
        oligo_entities_chains = np.unique(oligo_entities_chains)
        
    if len(oligo_entities_chains)>0: 
        non_glycosylation_oligos_chains = np.setdiff1d(oligo_entities_chains, glycosylation_chains)
    else: 
        non_glycosylation_oligos_chains = np.array([])

    for chain in non_glycosylation_oligos_chains:
        storeLog.append(f'lig_chain-{chain} (oligosaccharide)')
        lig = PandasPdb().read_pdb(f'{originpath}/{pdbname}')
        lig.df["HETATM"] = lig.df["HETATM"].query(f'chain_id=="{chain}"')
        lig.df["ATOM"] = lig.df["ATOM"].query(f'chain_id=="{chain}"')
        lig.to_pdb(path=f'{outpath}/{pdbname.split(".")[0]}_lig_chain-{chain}.pdb', records=['ATOM', 'HETATM'], 
                        gz=False, append_newline=True)
        storeLog.append(f'{pdbname.split(".")[0]}_lig_chain-{chain}')
        # ADD TO INDIVIDUAL LOG


    #############  EXTRACT LIGAND - METHOD 4  #############
    # Directly extract PRD from CIF
    data = cifdata[pdbcode.upper()]
    prd_chains = []
    if "_pdbx_molecule" in data:
        prd_chains = pd.DataFrame.from_dict(data["_pdbx_molecule"], orient="index").T
        prd_chains = prd_chains.asym_id.values
        new2old_chainID = pdb_auth_chains_dict[pdbcode.lower()]
        prd_chains = [new2old_chainID[x] for x in prd_chains]

    # Check if there are PRDs within a protein chain and not on its own chain - this is not correct
    # A PRD ID should automatically be associated with a ligand
    #prd_chains = np.setdiff1d(prd_chains, isprotein)
    
    for chain in prd_chains:
        storeLog.append(f'lig_chain-{chain} (oligosaccharide/peptide found from BIRD)')
        lig = PandasPdb().read_pdb(f'{originpath}/{pdbname}')
        lig.df["HETATM"] = lig.df["HETATM"].query(f'chain_id=="{chain}"')
        lig.df["ATOM"] = lig.df["ATOM"].query(f'chain_id=="{chain}"')
        lig.to_pdb(path=f'{outpath}/{pdbname.split(".")[0]}_lig_chain-{chain}.pdb', records=['ATOM', 'HETATM'], 
                        gz=False, append_newline=True)

    # extract PRD from sequence matching

    prd_seq_list = [ln.strip().split("\t") for ln in open(f"{home}/LigExtract/data/prd_cmpds_seq_full_experim.txt").readlines()]

    chainEntityID_AuthChainName_dict = pd.DataFrame.from_dict(data["_atom_site"], orient="index").T
    chainEntityID_AuthChainName_dict = chainEntityID_AuthChainName_dict[["label_entity_id","auth_asym_id"]].drop_duplicates()
    chainEntityID_AuthChainName_dict = {x:y for x,y in chainEntityID_AuthChainName_dict.values}

    data = cifdata[pdbcode.upper()]
    pdbseqs = []
    if "_entity_poly_seq" in data:
        pdbseqs = pd.DataFrame.from_dict(data["_entity_poly_seq"], orient="index").T
        pdbseqs = pdbseqs[["entity_id","mon_id"]]
        pdbseqs = pdbseqs.replace(chainEntityID_AuthChainName_dict)

    lig_seq_all = []
    # only current proteins are eligible for this. All other chains are already classed as non-proteins
    for protchainName in isprotein:
        ligseq = pdbseqs.query(f"entity_id == '{protchainName}'").mon_id.values
        if len(ligseq)> 0:
            ligseq = " ".join(ligseq) 
            lig_seq_all.append([protchainName, ligseq])
            #lig_seq_all = np.array(lig_seq_all)
    if len(lig_seq_all) == 0:
        lig_seq_all = np.array([[np.nan,np.nan]])

    found_match = np.intersect1d([a[1] for a in lig_seq_all],[x[1] for x in prd_seq_list])
    if len(found_match) == 0:
        storeLog = "\n".join(storeLog)
        print(storeLog)
        return(pdbcode)
    
    match_ids = []
    for m in found_match:
        chains = [a[0] for a in lig_seq_all if a[1]==m]
        prdID = [x[0] for x in prd_seq_list if x[1]==m][0]
        for chain in chains:
            lig = PandasPdb().read_pdb(f'{originpath}/{pdbname}')
            lig.df["HETATM"] = lig.df["HETATM"].query(f'residue_name!="HOH"')
            lig.df["HETATM"] = lig.df["HETATM"].query(f'chain_id=="{chain}"')
            lig.df["ATOM"] = lig.df["ATOM"].query(f'chain_id=="{chain}"')
            # remove 1-atom residues
            res_atoms = lig.df["HETATM"].groupby(["residue_name","residue_number"])["atom_number"].nunique()
            res2keep = res_atoms[res_atoms > 1].reset_index().residue_number
            lig.df["HETATM"] = lig.df["HETATM"][np.isin(lig.df["HETATM"].residue_number,res2keep)]

            if f'{pdbname.split(".")[0]}_lig_chain-{chain}.pdb' not in os.listdir(outpath):
                lig.to_pdb(path=f'{outpath}/{pdbname.split(".")[0]}_lig_chain-{chain}.pdb', records=['ATOM', 'HETATM'], gz=False, append_newline=True)
                storeLog.append(f'lig_chain-{chain}')
            else:
                storeLog.append(f'lig_chain-{chain} found again')
            logfile.write(f"rebuilt ligand detected in chain {chain}; associated with PRD ID {prdID} which means it is likely active\n")
    storeLog = "\n".join(storeLog)
    print(storeLog)
    return(pdbcode)

# multithreading
num_workers = os.cpu_count()-1


with ProcessPoolExecutor(max_workers=num_workers) as executor:
    #executor.map(extractorSplit, pdbs)
    futures = [executor.submit(extractorSplit, pdb_q) for pdb_q in pdbs]
    bar = Bar('Extracting Ligands... ', max=len(futures))
    results = []
    for f in as_completed(futures):
        results.append(f.result())
        bar.next()
    bar.finish()


print("\n\nFinished Ligand extraction.")
sys.stderr.write("\n\nFinished Ligand extraction.\n\n")

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
import itertools
import urllib
from urllib import parse,request
from time import sleep
from progress.bar import ChargingBar, Bar
import argparse
import gzip
import subprocess
import shlex
import json
import requests
home = str(Path.home())

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

print("\n------------------  Extracting all possible ligands from PDBs  ------------------\n")

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

if len(missing_downloads)>0 or len(extra_downloads)>0:
    if len(missing_downloads)>0: sys.stderr.write(f"\nWARNING!!!  There are {len(missing_downloads)} pdbs missing ({' '.join(missing_downloads)}). Please re-run the downloadPdbs.py command\n")
    if len(extra_downloads)>0: sys.stderr.write(f"\nWARNING!!!  There are {len(extra_downloads)} extra pdbs ({' '.join(extra_downloads)}). Please remove them from the PDBs directory you are passing to -d\n")
    sys.exit(123)

uniprot_pdb = pd.read_csv(f"{pdbpath.split('/')[-1]}_process_uniprot_chains.txt", sep="\t")
uniprot2pdb_secondary = copy(uniprot_pdb)
# pdb  uniprot_in_pdbfile  updated_uniprot

uniprot_pdb = uniprot_pdb[["uniprot_in_pdbfile","updated_uniprot"]].drop_duplicates()
uniprot_pdb["updated_uniprot"] = uniprot_pdb["updated_uniprot"]+","
uniprot_pdb = uniprot_pdb.groupby("uniprot_in_pdbfile")["updated_uniprot"].sum().reset_index()
uniprot_pdb["updated_uniprot"]=[x.split(",")[:-1] for x in uniprot_pdb.updated_uniprot]


def findpeptides(pdbfile):
    pdb_raw_text = PandasPdb().read_pdb(pdbfile)
    source_dat = " ".join(pdb_raw_text.df['OTHERS'].query("record_name == 'SOURCE'").entry.values).split("MOL_ID:")
    source_synth = [x for x in source_dat if "SYNTHESI" in x and "CHEMICAL" in x and "ORGANISM_SCIENTIFIC:" not in x]
    cmpd_dat = " ".join(pdb_raw_text.df['OTHERS'].query("record_name == 'COMPND'").entry.values).split("MOL_ID:")
    cmpd_inhib = [x for x in cmpd_dat if "PEPTIDE" in x]
    covered_id_1 = [x.split(";")[0].strip() for x in source_synth]
    covered_id_2 = [x.split(";")[0].strip() for x in cmpd_inhib]
    covered_id = np.intersect1d(covered_id_1, covered_id_2)
    if len(covered_id)==0: return([])
    peptide_chains = [x.split(" CHAIN: ")[1].split(";")[0].split(",") for x in cmpd_dat if x.split(";")[0].strip() in covered_id]
    peptide_chains = [x.strip() for x in np.hstack(peptide_chains)]
    return(peptide_chains)
    

def findlink(residues_atom, residues_hetatm):
    linked_hetatm = []
    for link in links:
        for res1 in residues_hetatm:
            for res2 in residues_atom:
                if res1 in link and res2 in link:
                    l = [x for x in link.split(" ") if x!=""]
                    linked_hetatm.append(int(res1[-4:].strip()))
                    break
    return(linked_hetatm)


def findAllLinks(residues_lst, restypes):
    res_chain_q = [x[4] for x in residues_lst]
    allLinks = []
    for link in links:
        if link[21] in res_chain_q and link[51] in res_chain_q:
            l1 = link[17:26]
            l2 = link[47:56]
            foundlink = [l1, l2]
            allLinks.append(foundlink)
    allLinks_ids=[(x,y) for x,y in allLinks]
    if len(allLinks_ids)==0:
        return allLinks_ids
    G = nx.Graph()
    G.add_edges_from(allLinks_ids)

    if restypes == "HETATM" and len(residues_lst)==1: 
        # means there is no connections to be made as there is only one HETATM residue in this chain
        return([])

    if restypes == "HETATM" and len(residues_lst)>1:
        all_connect = [int(x[-4:].strip()) for x in list(G.nodes)]
        return(all_connect)
    
    if restypes == "ATOM" and len(residues_lst)>1: 
        prot_link = zip(residues_lst[:-1],residues_lst[1:])
        prot_link_ids=[(x,y) for x,y in prot_link]
        G.add_edges_from(prot_link_ids)
    else: prot_link_ids=["0x0"]

    catch_neigh = residues_lst[:]
    expand=True
    while expand == True:
        new_neigh = np.unique(np.hstack([list(G.neighbors(x)) for x in catch_neigh]))
        new_neigh = np.setdiff1d(new_neigh,catch_neigh)
        if len(new_neigh) > 0:
            catch_neigh = np.hstack([catch_neigh,new_neigh])
        else:
            expand=False

    catch_neigh = [int(x[-4:].strip()) for x in catch_neigh]
    return(catch_neigh)


# some of these might be residues and others HETATM
# go back to the link list and get more info so that this list can be passed to both ATOM and HETATM


def countMolsAtoms(formula):
    # takes chemical formula from "FORMUL" strings
    formula=formula.replace("*","")
    wheremolcount = formula.find("(")
    if wheremolcount == -1: molcount = 1
    elif wheremolcount == 0: molcount = 1
    else: molcount = int(formula[:wheremolcount])
    # clean formula
    formula = formula.split("(")[-1].split(")")[0]
    # some molecules have empty molformula; avoid excluding them at this stage
    if formula.strip()=="": return(np.nan,molcount)
    formula = formula.split(" ")
    if "+" in formula[-1] or "-" in formula[-1]: formula = formula[:-1]
    formula = (" ").join(formula)
    formula = "".join([x for x in formula if x in string.digits or x==" "]).split(" ")
    formula = [x if len(x)>0 else 1 for x in formula]
    formula = np.array(formula).astype(float).sum()
    return(formula, molcount)


try : os.listdir(outpath)
except FileNotFoundError: os.makedirs(outpath)
pdbs = [x for x in os.listdir(pdbpath) if x.endswith(".pdb")]


keep_others = ["HEADER","TITLE","LINK","CONECT", "MASTER","END"]

if len(pdbs)==0:
    print(f" There are no PDB files in the target folder {targetpath}.\n Make sure your PDB files have a '.pdb' extension.")
    sys.exit(123)


subprocess.run("mkdir -p cifs", shell=True)
pdbs2extr = np.setdiff1d([x.split(".")[0] for x in pdbs], [x.split(".")[0] for x in os.listdir("cifs")])
if len(pdbs2extr)>0: 
    bar = Bar('Extracting CIFs... ', max=len(pdbs2extr))
    for pdbname in pdbs2extr:
        bar.next()
        pdbcode = pdbname.split(".")[0]
        if f'{pdbcode}.cif' not in os.listdir("cifs"):
            print(f'Download {pdbcode}.cif')
            sleep(2)
            subprocess.run(f"wget https://files.rcsb.org/view/{pdbcode.lower()}.cif --quiet -O cifs/{pdbcode}.cif", shell=True)
    
    sys.stderr.write("\n")


# After all cifs have been downloaded, screen them to grab PRD IDs
bar = Bar('Extract all chains with PRD code... ', max=len(pdbs))
all_prds_chains = []
for pdb in pdbs:
    bar.next()
    pdbcode = pdb.split(".")[0]
    ciffile = open(f"cifs/{pdbcode}.cif").readlines()
    ciffile = "".join(ciffile).split("# \n")
    ciffile = [x for x in ciffile if "_pdbx_molecule." in x]
    if len(ciffile)>0:
        ciffile = ciffile[0].strip().split("\n")
        if ciffile[0].startswith("loop_"):
            headers = [x.split(".")[-1].strip() for x in ciffile[1:] if x.startswith("_pdbx_molecule")]
            body = ciffile[len(headers)+1:]
            body = pd.DataFrame([shlex.split(x.strip(), posix=False) for x in body], columns=headers)
            body.loc[:,"pdb"] = pdbcode
            # translate asym_id to author ID
            all_prds_chains.append(body) # asym_id = label_asym_id (the new label)

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
    ciffile = open(f"cifs/{pdbcode}.cif").readlines()
    ciffile = "".join(ciffile).split("# \n")
    ciffile = [x for x in ciffile if "_atom_site." in x and "'_atom_site." not in x]
    if len(ciffile)>0:
        ciffile = ciffile[0].strip().split("\n")
        if ciffile[0].startswith("loop_"):
            headers = [x.split(".")[-1].strip() for x in ciffile[1:] if x.startswith("_atom_site")]
            body = ciffile[len(headers)+1:]
            line_length = max([len(ln.split()) for ln in body])
            lines2mergeup = []
            for i,ln in enumerate(body[:-1]):
                if len(ln.split())<line_length/2:
                    lines2mergeup.append(i)
            for i in lines2mergeup:
                body[i-1] = body[i-1]+body[i]
                body[i] = ''
            # protect spaces inside quotation marks
            body = pd.DataFrame([shlex.split(ln.strip(), posix=False) for ln in body], columns=headers)
            pdb_auth_chains = body[["auth_asym_id", "label_asym_id"]].drop_duplicates()
            pdb_auth_chains_dict[pdbcode] = {newlabel:authlabel for authlabel,newlabel in pdb_auth_chains[["auth_asym_id", "label_asym_id"]].values}


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


bar = Bar('Extracting Ligands... ', max=len(pdbs))

for pdbname in pdbs:
    bar.next()
    print(f"#######  Finding ligands in  {pdbname}  #######")
    pdbcode = pdbname.upper().split(".")[0]
    pdb = open(f"{pdbpath}/{pdbname}").readlines()
    if len(pdb) == 0:
        print(f"! {pdbname} has zero lines. Please verify the file is valid.")
        continue
    
    try: uniprotid = uniprot2pdb.query(f"pdb == '{pdbcode.upper()}'").uniprot.values[0]
    except IndexError:
        print(f"Unable to find {pdbcode.upper()} in {uniprot2pdb_file}. Make sure this pdb code is in the file, in all capitals.")
        sys.exit(123)
    # remove PDBs where the protein actually is not found within DBREF
    uniprot_refs = [x for x in pdb if x.startswith('DBREF')]
    uniprotid_alternates = uniprot2pdb_secondary.query(f"pdb == '{pdbcode.lower()}' and (updated_uniprot == '{uniprotid}' or uniprot_in_pdbfile == '{uniprotid}')").uniprot_in_pdbfile.values

    chains_w_Uniprot = np.unique([x[12] for x in uniprot_refs if x[33:41].strip() in uniprotid_alternates])
    extra_prot_chains = []
    # extend the mapping beyond what is in the file using the mapping from SIFTS
    with gzip.open(f"{home}/LigExtract/data/pdb_chain_uniprot.csv.gz") as f:
        for ln in f:
            ln=ln.decode("utf-8").strip().split(",")
            if ln[0]==pdbcode.lower():
                extra_prot_chains.append(ln[1:3])
    extra_prot_chains = np.unique([c for c,unp in extra_prot_chains]) 
    
    if len(chains_w_Uniprot)==0 and uniprotid not in uniprot2pdb_secondary.updated_uniprot.values:
        # Finally only delete if this pdb is not being used by other Uniprot IDs
        other_uniprots = uniprot2pdb_secondary.query(f"pdb == '{pdbcode.lower()}'").updated_uniprot.values
        other_uniprots = np.setdiff1d(other_uniprots, uniprotid)
        if np.in1d(uniprot2pdb.uniprot.unique(), other_uniprots).any() == False:
            print(f"{pdbcode} does not have any chain for ANY of the desired uniprotIDs in the full collection of targets! REMOVED protein")
        os.remove(f'{pdbpath}/{pdbname}')
        continue
    
    # initialise file to store as a log
    outfile = open(f'{outpath}/{pdbname.split(".")[0]}_ligand_extraction.log',"w")
    outfile.write(f" ###### LOG OF LIGAND EXTRACTION FOR {pdbname} ######\n\n")
    outfile.write(f"results obtained from PDBs in {originpath}\n")
    outfile.write(f"results stored in {targetpath}\n\n")
    
    cif = open(f"cifs/{pdbcode.lower()}.cif").read().split("# \n")
    if "_pdbx_entity_nonpoly." in "".join(cif):
        block = [x for x in cif if "_pdbx_entity_nonpoly" in x if x.split("\n")[1].startswith("_pdbx_entity_nonpoly")][0].replace(";","").split("\n")[:-1]
        if block[0]=="loop_":
            headers = [x.replace("_pdbx_entity_nonpoly.","").strip() for x in block if x.startswith("_pdbx")]
            body = []
            for x in block[1:]:
                if x.startswith("_pdbx") == False:
                    body.append(shlex.split(x.strip(), posix=False))
            nonpoly = [x[-1] for x in body if len(x)>0] 
        else:
            nonpoly = [x.strip().split()[-1] for x in block if x.startswith("_pdbx_entity_nonpoly.comp_id")]
    else:
        nonpoly = []

    
    links = [x.strip() for x in pdb if x.startswith("LINK")]
    links = [x.strip() for x in links if float(x.split(" ")[-1])<2.1 and " HOH " not in x] # max covalent distance
    
    
    chaintype = [x for x in pdb if x.startswith('DBREF ')]
    notaprotein = [x[12] for x in chaintype if x[33:41].strip() not in uniprot_pdb.uniprot_in_pdbfile.values] # OK
    notaprotein = np.setdiff1d(notaprotein,extra_prot_chains)
    isprotein = [x[12] for x in chaintype if x[33:41].strip() in uniprot_pdb.uniprot_in_pdbfile.values] # chain
    isprotein = np.union1d(isprotein,extra_prot_chains)
    
    # add chains that may not appear to be a protein from DBREF but have E.C. code in the COMPND section
    moltype = [x for x in pdb if x.startswith('COMPND')]
    chunks = np.argwhere(["MOL_ID:" in x for x in moltype]).flatten()
    chunks = list(zip(chunks,list(chunks[1:])+[len(moltype)]))
    moltype = [moltype[a:b] for a,b in chunks]
    moltype_prot = []
    for molblock in moltype:
        #if "CHAIN:" in "".join(molblock) and "EC:" in "".join(molblock):
        #    moltype_prot.append([molline[10:].split("CHAIN:")[-1].split(";")[0].strip() for molline in molblock if "CHAIN:" in molline and molline[10:].strip().startswith("CHAIN")])
        for line_n, molline in enumerate(molblock): 
            if "CHAIN:" in molline and molline[10:].strip().startswith("CHAIN"):
                molline = molline.strip()
                currline = copy(line_n)
                while molline.endswith(","):
                    molline_add = molblock[currline+1][10:].strip()
                    molline = molline+molline_add
                    currline+=1
                moltype_prot.append(molline[10:].split("CHAIN:")[-1].split(";")[0].strip())

    if len(moltype_prot)>0:
        moltype_prot = [x.split(",") for x in np.hstack(moltype_prot)]
        moltype_prot = np.hstack(moltype_prot)
        moltype_prot = [x.strip() for x in moltype_prot]
    isprotein = np.union1d(isprotein,moltype_prot)

    notaprotein = np.setdiff1d(notaprotein,isprotein)
    
    # add found peptides to the chain lists
    peptidechains = findpeptides(f"{pdbpath}/{pdbname}")
    print(f"peptide chains for {pdbcode}:", ",".join(peptidechains))
    isprotein = np.setdiff1d(isprotein, peptidechains)
    notaprotein = np.union1d(notaprotein, peptidechains)
    
    seq = [x for x in pdb if x.startswith('SEQRES')]

    # gather protein residues
    protres = []

    for chain in isprotein:
        r = [x.strip()[19:].split() for x in seq if chain == x[11]]
        r = np.hstack(r)
        r = [x for x in r if x!=""]
        protres.append(r)

    if len(protres)>0: protres = np.hstack(protres)

    # gather ligands 
    ligands = [x.strip() for x in pdb if x.startswith('HETNAM ')]
    ligands = list(np.unique([x[11:14].strip() for x in ligands]))
    # ligand among the protein sequence (e.g. 4H9M) cannot be considered; warning: a standard residue code might appear in HETNAM
    formulas = [x for x in pdb if x.startswith('FORMUL ')]
    modres = [x[12:15] for x in pdb if x.startswith('MODRES ')]
    modres2 = [x[12:15] for x in pdb if x.startswith('SEQADV ') and "MODIFIED RESIDUE" in x]
    modres = list(np.unique(np.hstack([modres,modres2])))
    ligands = np.setdiff1d(ligands,modres)
    filtered_ligs = []
    for x in ligands:
        if x in protres and x not in nonpoly: continue
        else: filtered_ligs.append(x) # e.g. ARG in 4kx9 which is a ligand 
            
    ligands = copy(filtered_ligs)

    ###########  EXTRACT LIGAND - METHOD 1  ###########
    # This method extracts ligands that may exist (partly or entirely) among the residues section.
    # This happens when the ligand is a peptide or has aminoacids in its structure.
    #print(notaprotein)
    if len(notaprotein)==0:
        ligres=[]
    else:
        ligres = []
        for chain in notaprotein:
            r = [x.strip() for x in seq if f' {chain}  ' in x]
            r = [x[x.find(f' {chain}  ')+8:].split(" ") for x in r]
            r = np.hstack(r)
            r = [x for x in r if x!=""]
            ligres.append(r)
            
        ligres = np.unique(np.hstack(ligres))
        
        for chain in notaprotein:
            lig = PandasPdb().read_pdb(f'{originpath}/{pdbname}')
            lig.df["HETATM"] = lig.df["HETATM"].query(f'chain_id=="{chain}"')
            lig.df["ATOM"] = lig.df["ATOM"].query(f'chain_id=="{chain}"')
            # clean from solvents; not exhaustive 
            lig.df["HETATM"] = lig.df["HETATM"].query(f'residue_name!="HOH"')
            res_atoms = lig.df["HETATM"].groupby(["residue_name","residue_number"])["atom_number"].nunique()
            res2keep = res_atoms[res_atoms >= 1].reset_index().residue_number # some residues are just the carbon
            
            lig.df["HETATM"] = lig.df["HETATM"][np.in1d(lig.df["HETATM"].residue_number,res2keep)]

            res_atom = lig.df['ATOM'][["residue_name","chain_id","residue_number"]].drop_duplicates().values
            res_hetatm = lig.df['HETATM'][["residue_name","chain_id","residue_number"]].drop_duplicates().values
            
            if len(res_atom) > 0: # part/all of the ligand is in ATOM
                res_atom =["{:>3}".format(a) + "{:>2}".format(b) + "{:>4}".format(c) for a,b,c in res_atom.astype(str)]
                res_hetatm =["{:>3}".format(a) + "{:>2}".format(b) + "{:>4}".format(c) for a,b,c in res_hetatm]
                # see which res connect from HETATM to ATOM
                
                # At the moment this does not cover all exceptions as some cases fail to mention ATOM-HETATM links that exist (e.g. 1de7)
                all_linked_atms = np.unique(findAllLinks(res_atom, "ATOM") + findAllLinks(res_hetatm,"HETATM"))
                # add additional residues that are in the sequence but not in all_linked_atms
                chainres = np.hstack([x.strip()[19:].split() for x in seq if chain == x[11]])
                chainres = [int(x[-4:].strip()) for x in res_hetatm if x[:3] in chainres]
                all_linked_atms = np.hstack([all_linked_atms, chainres])

                
                lig.df['ATOM'] = lig.df['ATOM'][lig.df['ATOM'].chain_id==chain]
                lig.df['HETATM'] = lig.df['HETATM'][lig.df['HETATM'].chain_id==chain]
                lig.df['HETATM'] = lig.df['HETATM'][np.in1d(lig.df['HETATM'].residue_number,all_linked_atms)] # linked_hetatm
                lig.df['ANISOU'] = lig.df['ANISOU'][lig.df['ANISOU'].chain_id==chain]
                lig.df["OTHERS"] = lig.df["OTHERS"][np.in1d(lig.df["OTHERS"].record_name, keep_others)]
                lig.to_pdb(path=f'{outpath}/{pdbname.split(".")[0]}_lig_chain-{chain}.pdb', records=['ATOM', 'HETATM', "OTHERS"], gz=False, append_newline=True)
                print(f'lig_chain-{chain}')
                outfile.write(f"rebuilt ligand detected in chain {chain}\n")

            if len(res_atom) == 0: # ligand is entirely in HETATM but has multiple residues
                res_hetatm =["{:>3}".format(a) + "{:>2}".format(b) + "{:>4}".format(c) for a,b,c in res_hetatm]
                linked_hetatm = findAllLinks(res_hetatm, "HETATM")
                lig.df['ATOM'] = lig.df['ATOM'][lig.df['ATOM'].chain_id==chain]
                lig.df['HETATM'] = lig.df['HETATM'][lig.df['HETATM'].chain_id==chain]
                if len(linked_hetatm)>0:
                    lig.df['HETATM'] = lig.df['HETATM'][np.in1d(lig.df['HETATM'].residue_number,linked_hetatm)]
                if len(linked_hetatm)==0:
                    # untie with HET
                    het = [x.split() for x in pdb if x.startswith('HET ')]
                    het = [x[1] for x in het if x[2]==chain]
                    lig.df['HETATM'] = lig.df['HETATM'][np.in1d(lig.df['HETATM'].residue_name,het)]
                lig.df['ANISOU'] = lig.df['ANISOU'][lig.df['ANISOU'].chain_id==chain]
                lig.df["OTHERS"] = lig.df["OTHERS"][np.in1d(lig.df["OTHERS"].record_name, keep_others)]
                lig.to_pdb(path=f'{outpath}/{pdbname.split(".")[0]}_lig_chain-{chain}.pdb', records=['ATOM', 'HETATM', "OTHERS"], gz=False, append_newline=True)
                print(f'lig_chain-{chain}')
                outfile.write(f"rebuilt ligand detected in chain {chain}; exclusively found in HETATM records\n")

    
    ##########  EXTRACT LIGAND - METHOD 2  ###########
    # This method extracts ligands from the typical section where ligands are found
    cmpds_all = [cmpd[12:15].strip() for cmpd in formulas]
    
    cmpds = [cmpd[12:15].strip() for cmpd in formulas if ( countMolsAtoms(cmpd[19:70].strip())[0]>3 or np.isfinite(countMolsAtoms(cmpd[19:70].strip())[0])==False) and countMolsAtoms(cmpd[19:70].strip())[1]<10]
    # restore cmpds which are in BIRD
    lig_ids_in_bird = pd.read_csv(f"{home}/LigExtract/data/prd_to_pdb_ligIDs.txt", sep="\t").pdblig_ID.values
    cmpds = np.union1d(cmpds, np.intersect1d(cmpds_all,lig_ids_in_bird))
    cmpds_removed = np.setdiff1d(cmpds_all,cmpds)
    if len(cmpds_removed)>0: print("The following compounds do not comply wih the inclusion criteria of minimum 3 heavy atoms and maximum of 10 repeated molecules:",cmpds_removed)
    cmpds = np.setdiff1d(cmpds, modres) # ligres
    cmpds = np.intersect1d(cmpds, ligands)

    if len(cmpds)>0:
        for cpd in cmpds:
            message = ""
            if cpd in modres:
                continue
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
                            params = {"query": {
                                      "type": "terminal","service": "text","parameters": {
                                      "attribute": "rcsb_nonpolymer_instance_feature_summary.comp_id",
                                      "operator": "exact_match",
                                      "value": f"{cpd}"}
                                      },
                                     "request_options": {"return_counts": True},"return_type": "entry"}
                            result = requests.get("https://search.rcsb.org/rcsbsearch/v1/query", {"json": json.dumps(params, separators=(",", ":"))})
                            res_stat = result.status_code
                            result = result.content.decode("utf-8")
                            if result == "" and res_stat not in [200,204]: #,204
                                print(f"{cpd}: unssuccessful request for PDB count\n") 
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
            print(cpd)
            links_warning = []
            for link in links:
                if f' {cpd} ' in link:
                    links_warning.append(link)
            
            if len(links_warning)>0:
                outfile.write(f"{cpd}   ligand connected to other molecules. {message}\n")
                print(f"    WARNING!  Ligand {cpd} appears to be connected to other molecules (see .log file).\n")
                print(f"    This means it (1) is just a portion of the ligand or (2) is covalently bound (which requires further struture correction to return both ligand and protein to original unbound structure)\n")
                outfile.write("\t"+"            ---- res 1 ----               ---- res 2 ----                 bond (A)\n")
                for w in links_warning:
                    outfile.write("\t"+w+"\n")
            else:
                outfile.write(f"{cpd}   ligand detected. {message}\n")
            
            lig = PandasPdb().read_pdb(f'{originpath}/{pdbname}')
            chains_lig = lig.df["HETATM"][lig.df["HETATM"].residue_name==cpd].chain_id.unique()
            for chain in chains_lig:
                lig = PandasPdb().read_pdb(f'{originpath}/{pdbname}')
                lig.df["HETATM"] = lig.df["HETATM"].query(f'residue_name=="{cpd}" and chain_id=="{chain}"')
                lig.df["ATOM"] = lig.df["ATOM"].query(f'residue_name=="NOTHING" and chain_id=="{chain}"')
                lig.df["ANISOU"] = lig.df["ANISOU"].query(f'residue_name=="NOTHING" and chain_id=="{chain}"')
                lig.df["OTHERS"] = lig.df["OTHERS"][np.in1d(lig.df["OTHERS"].record_name, keep_others)]
                lig.to_pdb(path=f'{outpath}/{pdbname.split(".")[0]}_chain-{chain}_lig-{cpd}.pdb', records=['ATOM', 'HETATM', "OTHERS"], 
                                gz=False, append_newline=True)
                print(f'{pdbname.split(".")[0]}_chain-{chain}_lig-{cpd}')

    ########## METHOD 3: collecting oligosacharide chain ligands ##########
    cif = open(f"cifs/{pdbcode.lower()}.cif").read().split("# \n")

    if "_pdbx_entity_branch.entity_id" in "".join(cif):
        block = [x for x in cif if "_pdbx_entity_branch." in x and x.split("\n")[1].startswith("_pdbx_entity_branch.")][0].split("\n")[:-1]
        if block[0]=="loop_":
            headers = [x.replace("_pdbx_entity_branch.","").strip() for x in block if x.startswith("_pdbx")]
            body = [x.strip().split() for x in block[1:] if x.startswith("_pdbx") == False]
            oligo_entities = pd.DataFrame(body, columns=headers)
            oligo_entities = oligo_entities.query("type == 'oligosaccharide'").entity_id.values
        else:
            is_oligo = [x.strip().split()[-1] for x in block if "_pdbx_entity_branch.type" in x][0]
            if is_oligo == "oligosaccharide":
                oligo_entities = [x.strip().split()[-1] for x in block if "_pdbx_entity_branch.entity_id" in x]
            else:
                oligo_entities = []
    else:
        oligo_entities = []


    if "_struct_conn.pdbx_role" in "".join(cif):
        block = [x for x in cif if "_struct_conn." in x and x.split("\n")[1].startswith("_struct_conn.")][0].split("\n")[:-1]
        if block[0]=="loop_":
            headers = [x.replace("_struct_conn.","").strip() for x in block if x.startswith("_struct_conn")]
            body = [x for x in block[1:] if x.startswith("_struct_conn") == False]
            join_ln = np.argwhere([len(shlex.split(x))<len(headers) for x in body]).flatten()
            for x in join_ln:
                # if current line is empty, this is a continuation line: bypass.
                if len(body[x].split())==0: continue
                body[x] = body[x]+body[x+1]
                body[x+1]=""
            body = [shlex.split(x.strip(), posix=False) for x in body if "Glycosylation" in x]
            glycosylation_chains = pd.DataFrame(body, columns=headers)
            glycosylation_chains = glycosylation_chains[["ptnr1_auth_asym_id","ptnr2_auth_asym_id"]] 
            glycosylation_chains = glycosylation_chains.values.flatten()
            glycosylation_chains = np.setdiff1d(glycosylation_chains, isprotein)
        else:
            glycosylation_chains = np.array([x.strip().split() for x in block])
            glycosylation_chains = pd.DataFrame(glycosylation_chains[:,1].reshape(1,-1), columns = [x.replace("_struct_conn.","") for x in glycosylation_chains[:,0]])
            if "Glycosylation" not in glycosylation_chains.pdbx_role.values[0]:
                glycosylation_chains = []
            else:
                glycosylation_chains = glycosylation_chains[["ptnr1_auth_asym_id","ptnr2_auth_asym_id"]].values.flatten()
    else:
        glycosylation_chains = []

    block = [x for x in cif if "_atom_site." in x][0].split("\n")[:-1]
    headers = [x.replace("_atom_site.","").strip() for x in block if x.startswith("_atom_site.")]
    body = [x.strip().split() for x in block[1:] if x.startswith("_atom_site.") == False]
    old_new_chainID = pd.DataFrame(body, columns=headers)
    new2old_chainID = old_new_chainID[["label_asym_id","auth_asym_id"]].drop_duplicates()
    new2old_chainID = {x:y for x,y in new2old_chainID.values}
    
    if "_pdbx_branch_scheme." in "".join(cif):
        block = [x for x in cif if "_pdbx_branch_scheme." in x and x.split("\n")[1].startswith("_pdbx_branch_scheme.")][0].split("\n")[:-1]
        if block[0]=="loop_":
            headers = [x.replace("_pdbx_branch_scheme.","").strip() for x in block if x.startswith("_pdbx")]
            body = [shlex.split(x.strip(), posix=False) for x in block[1:] if x.startswith("_pdbx") == False]
            oligo_entities_chains = pd.DataFrame(body, columns=headers)
            oligo_entities_chains = oligo_entities_chains[np.in1d(oligo_entities_chains.entity_id, oligo_entities)]
        else:
            oligo_entities_chains = np.array([x.strip().split() for x in block])
            oligo_entities_chains = pd.DataFrame(oligo_entities_chains[:,1].reshape(1,-1), columns = [x.replace("_pdbx_branch_scheme.","") for x in oligo_entities_chains[:,0]])
        oligo_entities_chains = [new2old_chainID[x] for x in oligo_entities_chains.asym_id.values]
    else:
        oligo_entities_chains = []
    
    
    if len(oligo_entities_chains)>0: non_glycosylation_oligos_chains = np.setdiff1d(oligo_entities_chains, glycosylation_chains)
    else: non_glycosylation_oligos_chains = []

    for chain in non_glycosylation_oligos_chains:
        print(f'lig_chain-{chain} (oligosaccharide)')
        lig = PandasPdb().read_pdb(f'{originpath}/{pdbname}')
        lig.df["HETATM"] = lig.df["HETATM"].query(f'chain_id=="{chain}"')
        lig.df["ATOM"] = lig.df["ATOM"].query(f'chain_id=="{chain}"')
        lig.df["ANISOU"] = lig.df["ANISOU"].query(f'chain_id=="{chain}"')
        lig.df["OTHERS"] = lig.df["OTHERS"][np.in1d(lig.df["OTHERS"].record_name, keep_others)]
        lig.to_pdb(path=f'{outpath}/{pdbname.split(".")[0]}_lig_chain-{chain}.pdb', records=['ATOM', 'HETATM', "OTHERS"], 
                        gz=False, append_newline=True)
        print(f'{pdbname.split(".")[0]}_lig_chain-{chain}')
        # ADD TO INDIVIDUAL LOG


    #############  EXTRACT LIGAND - METHOD 4  #############
    # Directly extract PRD from CIF
    if "_pdbx_molecule.prd_id" in "".join(cif):
        block = [x for x in cif if "_pdbx_molecule." in x and x.split("\n")[1].startswith("_pdbx_molecule.")][0].split("\n")[:-1]
        if block[0]=="loop_":
            headers = [x.replace("_pdbx_molecule.","").strip() for x in block if x.startswith("_pdbx")]
            body = [x.strip().split() for x in block[1:] if x.startswith("_pdbx") == False]
            prd_chains = pd.DataFrame(body, columns=headers)
        else:
            prd_chains = np.array([x.strip().split() for x in block])
            
            prd_chains = pd.DataFrame(prd_chains[:,1].reshape(1,-1), columns = [x.replace("_pdbx_molecule.","") for x in prd_chains[:,0]])
        prd_chains = prd_chains.asym_id.values

        # translate new chain ID (asym_id) to original chain ID (author_asym_id)
        block = [x for x in cif if "_atom_site." in x][0].split("\n")[:-1]
        headers = [x.replace("_atom_site.","").strip() for x in block if x.startswith("_atom_site.")]
        body = [x.strip().split() for x in block[1:] if x.startswith("_atom_site.") == False]
        old_new_chainID = pd.DataFrame(body, columns=headers)
        new2old_chainID = old_new_chainID[["label_asym_id","auth_asym_id"]].drop_duplicates()
        new2old_chainID = {x:y for x,y in new2old_chainID.values}

        # apply translation of new-2-old chain IDs
        prd_chains = [new2old_chainID[x] for x in prd_chains]

    else:
        prd_chains = []

    # Check if there are PRDs within a protein chain and not no its own chain
    prd_chains = np.setdiff1d(prd_chains, isprotein)
    
    for chain in prd_chains:
        print(f'lig_chain-{chain} (oligosaccharide found from BIRD)')
        lig = PandasPdb().read_pdb(f'{originpath}/{pdbname}')
        lig.df["HETATM"] = lig.df["HETATM"].query(f'chain_id=="{chain}"')
        lig.df["ATOM"] = lig.df["ATOM"].query(f'chain_id=="{chain}"')
        lig.df["ANISOU"] = lig.df["ANISOU"].query(f'chain_id=="{chain}"')
        lig.df["OTHERS"] = lig.df["OTHERS"][np.in1d(lig.df["OTHERS"].record_name, keep_others)]
        lig.to_pdb(path=f'{outpath}/{pdbname.split(".")[0]}_lig_chain-{chain}.pdb', records=['ATOM', 'HETATM', "OTHERS"], 
                        gz=False, append_newline=True)

    # extract PRD from sequence matching

    prd_seq_list = [ln.strip().split("\t") for ln in open(f"{home}/LigExtract/data/prd_cmpds_seq_full_experim.txt").readlines()]
    pdbseqs = [ln.strip() for ln in pdb if ln.startswith("SEQRES ")]
    if len(pdbseqs) > 0:
        lig_seq_all = [[seq[0:12][-1]," ".join(seq[12:].split()[1:])+" "] for seq in pdbseqs]
        lig_seq_all = pd.DataFrame(lig_seq_all)
        lig_seq_all = lig_seq_all.groupby(0)[1].sum()
        lig_seq_all = [[x,y.strip()] for x,y in lig_seq_all.reset_index().values]
    if len(pdbseqs) == 0:
        lig_seq_all = [[np.nan,np.nan]]

    found_match = np.intersect1d([a[1] for a in lig_seq_all],[x[1] for x in prd_seq_list])
    if len(found_match) == 0:
        continue
    
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
            lig.df["HETATM"] = lig.df["HETATM"][np.in1d(lig.df["HETATM"].residue_number,res2keep)]

            lig.df["ANISOU"] = lig.df["ANISOU"].query(f'chain_id=="{chain}"')
            lig.df["OTHERS"] = lig.df["OTHERS"][np.in1d(lig.df["OTHERS"].record_name, keep_others)]
            if f'{pdbname.split(".")[0]}_lig_chain-{chain}.pdb' not in os.listdir(outpath):
                lig.to_pdb(path=f'{outpath}/{pdbname.split(".")[0]}_lig_chain-{chain}.pdb', records=['ATOM', 'HETATM', "OTHERS"], gz=False, append_newline=True)
                print(f'lig_chain-{chain}')
            else:
                print(f'lig_chain-{chain} found again')
            outfile.write(f"rebuilt ligand detected in chain {chain}; associated with PRD ID {prdID} which means it is likely active\n")


print("\n\nFinished Ligand extraction.")
sys.stderr.write("\n\nFinished Ligand extraction.")

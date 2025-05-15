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
import re
from pdbecif.mmcif_io import CifFileReader
from glob import glob
home = os.path.realpath(__file__)
home = home.split("/LigExtract")[0]


def findpeptide(pdbcode):
    # first get chain2uniprot correspondence
    chain2unip = {}
    # extend the mapping beyond what is in the file using the mapping from SIFTS
    with gzip.open(f"{home}/LigExtract/data/pdb_chain_uniprot.csv.gz") as f:
        for ln in f:
            ln=ln.decode("utf-8").strip().split(",")
            if ln[0]==pdbcode.lower():
                chain2unip[ln[1]]=ln[2]
    # Next text according to sequence lengths
    ciffile = f"cifs/{pdbcode.lower()}.cif"
    cifdata = CifFileReader().read(ciffile)
    is_polym_syn = pd.DataFrame.from_dict(cifdata[pdbcode.upper()]["_entity"], orient="index").T
    # has to be syn, to reduce the risk of falsely finding a protein like thrombin light chain
    is_polym_syn = is_polym_syn.query("type == 'polymer' and src_method == 'syn'")
    if len(is_polym_syn) == 0:
        return("") #"no peptides"
    session = requests.Session()
    peptides = []
    for line in range(len(is_polym_syn)):
        resolve_chain = is_polym_syn.iloc[line].id
        resolve_chain_seq = pd.DataFrame.from_dict(cifdata[pdbcode.upper()]["_entity_poly"], orient="index").T
        chain_theor_seq = resolve_chain_seq.query(f"entity_id == '{resolve_chain}'")
        chain_name = chain_theor_seq.pdbx_strand_id.values[0].split(",")
        if "_pdbx_entity_src_syn" in cifdata[pdbcode.upper()]:
            synconst = pd.DataFrame.from_dict(cifdata[pdbcode.upper()]["_pdbx_entity_src_syn"], orient="index").T
            synconst = synconst.query(f"entity_id == '{resolve_chain}' and organism_scientific == 'synthetic construct'")
            if len(synconst) > 0: 
                peptides.append([resolve_chain, ",".join(chain_name)])
                continue
        if np.isin(chain_name,list(chain2unip.keys())).all() == False:
            peptides.append([resolve_chain, ",".join(chain_name)])
            continue
        uniprot = chain2unip[chain_name[0]]
        response = session.get(f"https://rest.uniprot.org/uniprotkb/search?query=accession_id:{uniprot}")
        prot_length = json.loads(response.text)["results"][0]['sequence']["length"]
        theor_len = len(chain_theor_seq.pdbx_seq_one_letter_code_can.values[0])
        #real_len = np.abs(eval(is_polym_syn.pdbx_fragment.values[0].split()[-1]))+1
        allres = pd.DataFrame.from_dict(cifdata[pdbcode.upper()]["_atom_site"], orient="index").T
        allres = allres.query("group_PDB == 'ATOM'")
        if "_pdbx_struct_mod_residue" in cifdata[pdbcode.upper()]:
            mod_res = pd.DataFrame.from_dict(cifdata[pdbcode.upper()]["_pdbx_struct_mod_residue"], orient="index").T
            mod_res = mod_res.query(f"auth_asym_id == '{chain_name[0]}'")
        else:
            mod_res = []
        real_len = allres.query(f"label_entity_id == '{resolve_chain}'").label_seq_id.nunique()
        real_len += len(mod_res)
        # for very short peptides
        if theor_len <= 20 and prot_length>=50:
            if real_len >= 3: 
                peptides.append([resolve_chain, ",".join(chain_name)])
        # larger peptides; more that 70% resolved
        elif real_len/theor_len >= 0.7 and theor_len/prot_length < 0.1:
            peptides.append([resolve_chain, ",".join(chain_name)])
        # larger peptides; fairly incomplete
        elif theor_len<40 and real_len>7 and theor_len/prot_length < 0.1:
            peptides.append([resolve_chain, ",".join(chain_name)])
            print(f"Warning: chain {chain_name} in {pdbcode} looks like a peptide but the majority of expected sequence is unresolved")
        else:
            peptides.append([None,""])
    #peptides stores chain number and chain name. output only the second
    return( ",".join([x[-1] for x in peptides if len(x[-1])>0]) )
    



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




def findAllLinks(residues_lst, restypes, links_table):
    if len(residues_lst) == 0:
        return([])
    if len(links_table) == 0:
        return([x for x in residues_lst])
    res_chain_q = np.unique([x.split("-")[1] for x in residues_lst]) # order: "residue_name","chain_id","residue_number"
    allLinks = []
    for ch in res_chain_q:
        foundlink = links_table.query(f"ptnr1_auth_asym_id == '{ch}' and ptnr2_auth_asym_id == '{ch}'")
        if len(foundlink) == 0:
            return([])
        foundlink1 = ["-".join(x) for x in foundlink[['ptnr1_auth_comp_id', 'ptnr1_auth_asym_id', 'ptnr1_auth_seq_id']].values]
        foundlink2 = ["-".join(x) for x in foundlink[['ptnr2_auth_comp_id', 'ptnr2_auth_asym_id', 'ptnr2_auth_seq_id']].values]
        foundlink = np.vstack(list(zip(foundlink1, foundlink2)))
        allLinks.append(list(foundlink))
    allLinks_ids = np.vstack(allLinks)
    if len(allLinks_ids)==0:
        return allLinks_ids
    G = nx.Graph()
    G.add_edges_from(allLinks_ids)
    if restypes == "HETATM" and len(residues_lst)==1: 
        # means there is no connections to be made as there is only one HETATM residue in this chain
        return([])
    if restypes == "HETATM" and len(residues_lst)>1:
        all_connect = [x for x in list(G.nodes)]
        return(all_connect)
    if restypes == "ATOM" and len(residues_lst)>1: 
        prot_link = copy(residues_lst)
        prot_link = zip(prot_link[:-1],prot_link[1:]) # pair all res to res+1
        prot_link_ids=[(x,y) for x,y in prot_link]
        G.add_edges_from(prot_link_ids)
    else: prot_link_ids=["0x0"]
    groups = nx.connected_components(G)
    all_groups_found =[]
    for group in groups:
        all_groups_found.append(sorted(group, key=lambda x: int(x.split('-')[-1])))
    largest_group=max(all_groups_found, key=lambda df: len(df))
    return([str(x) for x in largest_group]) #ABC-A-23





# some of these might be residues and others HETATM
# go back to the link list and get more info so that this list can be passed to both ATOM and HETATM

def countMolsAtoms(pdbcode):
    ciffile = f"cifs/{pdbcode.lower()}.cif"
    cifdata = CifFileReader().read(ciffile)
    formula = pd.DataFrame.from_dict(cifdata[pdbcode.upper()]["_chem_comp"], orient="index").T
    formula = formula.query("mon_nstd_flag != 'y'")[["id","name","formula"]]
    n_mols = pd.DataFrame.from_dict(cifdata[pdbcode.upper()]["_entity"], orient="index").T
    n_mols = n_mols.query(f"type != 'polymer' and type != 'water'")
    n_mols = n_mols[["pdbx_description", "pdbx_number_of_molecules"]]
    n_mols = n_mols.merge(formula, right_on="name", left_on="pdbx_description") #= int(n_mols)
    n_mols.loc[:,"pdbx_number_of_molecules"] = n_mols.pdbx_number_of_molecules.values.astype(int)
    n_mols.loc[:,"atom_types"] = [len(np.intersect1d(list(string.ascii_uppercase),list(x))) for x in n_mols.formula.values]
    clean_formul =[]
    for ln in n_mols.formula:
        isformul = [len(np.intersect1d(list(string.ascii_uppercase + string.ascii_lowercase), list(x) )) for x in ln.split()]
        isformul = np.array(ln.split())[np.array(isformul)>0]
        clean_formul.append(" ".join(isformul))  
    n_mols.loc[:,"clean_formula"] = clean_formul
    #formula = "".join([x for x in clean_formul if x in string.digits or x==" "]).split(" ")
    atom_counts_lst = []
    for ln in clean_formul:
        atom_counts = "".join([x for x in ln if x in string.digits or x==" "]).split(" ")
        atom_counts = [x if len(x)>0 else 1 for x in atom_counts]
        atom_counts_lst.append(np.array(atom_counts).astype(float).sum())
    n_mols.loc[:,"atom_counts"] = atom_counts_lst
    n_mols_ok = n_mols.query("atom_counts > 3 and pdbx_number_of_molecules < 10")
    n_mols_refuse = n_mols.loc[np.setdiff1d(n_mols.index, n_mols_ok.index)]
    return(n_mols_ok.id.values, n_mols_refuse.id.values)





def groupresidues(residues_lst, linkstable):
    res_chain_q = np.unique([x[4] for x in residues_lst])
    allLinks = []
    foundlinks_idx = [linkstable.query(f"ptnr1_auth_asym_id == '{ch}' and ptnr2_auth_asym_id == '{ch}'").index for ch in res_chain_q]
    foundlinks = linkstable.loc[np.unique(np.hstack(foundlinks_idx))]
    for linksLine in foundlinks.index:
        link = foundlinks.loc[linksLine]
        #if link[21] in res_chain_q and link[51] in res_chain_q:
        # emulate the spacing in PDB files, LINKS sections
        l1 = link[['ptnr1_auth_comp_id', 'ptnr1_auth_asym_id', 'ptnr1_auth_seq_id']].values[0] #link[17:26]
        l2 = link[['ptnr2_auth_comp_id', 'ptnr2_auth_asym_id', 'ptnr2_auth_seq_id']].values[0] #link[47:56]
        l1 = "{:>3}".format(l1[0])+"{:>2}".format(l1[1])+"{:>4}".format(l1[2])
        l2 = "{:>3}".format(l2[0])+"{:>2}".format(l2[1])+"{:>4}".format(l2[2])
        foundlink = [l1, l2]
        allLinks.append(foundlink)
    allLinks_ids=[(x,y) for x,y in allLinks]
    if len(allLinks_ids)==0:
        g = [x[-4:].strip() for x in residues_lst]
        g = pd.DataFrame(zip(g,range(1,len(g)+1)), columns=["resn","group"])
        g.loc[:,"resn"] = g.resn.values.astype(int)
        return g
    G = nx.Graph()
    G.add_edges_from(allLinks_ids)
    groups = []
    man_group = 0
    for r in residues_lst:
        if r not in G.nodes:
            groups.append(r[-4:].strip())
            continue
        expand=True
        catch_neigh = [r]
        while expand == True:
            new_neigh = np.unique(np.hstack([list(G.neighbors(x)) for x in catch_neigh]))
            new_neigh = np.setdiff1d(new_neigh,catch_neigh)
            new_neigh = [x for x in new_neigh if r[:3] in x]
            if len(new_neigh) > 0:
                catch_neigh = np.unique(np.hstack([catch_neigh,new_neigh]))
            else:
                expand=False
        catch_neigh = ",".join([x[-4:].strip() for x in catch_neigh])
        groups.append(catch_neigh)
    g_n = 0
    unique_groups = []
    for g in np.unique(groups):
        g_n+=1
        g = g.split(",")
        unique_groups.append(pd.DataFrame(zip(g, [g_n]*len(g)), columns=["resn","group"]))
    unique_groups = pd.concat(unique_groups)
    unique_groups.loc[:,"resn"] = unique_groups.resn.values.astype(int)
    return unique_groups
import pandas as pd
import numpy as np
from biopandas.pdb import PandasPdb
import sys
import os
from rdkit import Chem
from rdkit.Chem.Descriptors import ExactMolWt
from rdkit.Chem.rdmolfiles import MolFromPDBFile
from glob import glob
import subprocess
from time import sleep
from progress.bar import ChargingBar, Bar
import argparse
from string import ascii_uppercase
from rdkit import RDLogger
import gzip
from pathlib import Path
import networkx as nx
from pdbecif.mmcif_io import CifFileReader
from concurrent.futures import ProcessPoolExecutor, as_completed
#HOME=Path.home()
HOME = os.path.realpath(__file__)
HOME = HOME.split("/LigExtract")[0]
RDLogger.DisableLog('rdApp.*')

sys.path.append(f'{HOME}/LigExtract/bin')
from ligextract_utils import *


parser = argparse.ArgumentParser(description='This does a preliminary ligand clean-up removing the most obvious non ligands. A file with the pockets around each ligand will be produced.')
parser.add_argument('--pdbPath', type=str, required=True, dest="protein_dir", help='Path to directory containing all PDBs to process.')
parser.add_argument('--ligandsPath', type=str, required=True, dest="ligands_dir", help='Path to directory with extracted ligand files.')
parser.add_argument('--dist', type=float, required=True, dest="distmax", help='maximum distance from every ligand atom, used to detect near residues. A reasonable distance is 5-7 Angstrom.')
parser.add_argument('--uniprot2pdbFile', type=str, required=True, dest="uniprot2pdb_file", help='File with Uniprot-to-Pdb mappings, ending with *pdb_uniprot_filteredlist.txt')
parser.add_argument('--keeprepeats', type=str, required=True, dest="keeprepeats", help="Use '--keeprepeats y' to keep all ligands when multiple instances of the same ligand exist (e.g. 2C0-res357 and 2C0-res358 split into their individual files and both kept). Use '--keeprepeats n' to keep the instance of the same ligand with the largest pocket (e.g. From 2C0-res357 and 2C0-res358 keep just 2C0-res357).")

args = parser.parse_args()

protein_dir = args.protein_dir
ligands_dir = args.ligands_dir
distmax = args.distmax
uniprot2pdb_file = args.uniprot2pdb_file
keeprepeats = args.keeprepeats

if protein_dir[-1]=="/": protein_dir = protein_dir[:-1]
if ligands_dir[-1]=="/": ligands_dir = ligands_dir[:-1]


length = 90; pad_char = '-'

title="First Ligand Clean-up (crystallography additives, solvent, etc)"
padding_total = length - len(title) - 2
print(f"\n{pad_char * (padding_total // 2)} {title} {pad_char * (padding_total - (padding_total // 2))}\n")


uniprot2pdb = pd.read_csv(uniprot2pdb_file, sep="\t")
# this mapping has the most updated uniprot

if "pdb" not in uniprot2pdb.columns:
    print("ERROR!  the uniprot2pdb file must have a column called 'pdb'")
    sys.exit(123)

if "uniprot" not in uniprot2pdb.columns:
    print("ERROR!  the uniprot2pdb file must have a column called 'uniprot'")
    sys.exit(123)

uniprot2pdb_secondary = pd.read_csv(f"{protein_dir.split('/')[-1]}_process_uniprot_chains.txt", sep="\t")

# SIFTS pdb chain mapping
sifts_pdb2uniprot = []
with gzip.open(f"{home}/LigExtract/data/pdb_chain_uniprot.csv.gz") as f:
    for ln in f: # PDB  CHAIN   SP_PRIMARY
        ln=ln.decode("utf-8").strip().split(",")
        if ln[0].upper() in uniprot2pdb.pdb.values:
            sifts_pdb2uniprot.append(ln[0:3])

sifts_pdb2uniprot = pd.DataFrame(sifts_pdb2uniprot, columns = ["pdbs", "chainName", "uniprot"])
sifts_pdb2uniprot = sifts_pdb2uniprot[[x.startswith("NOR")==False for x in sifts_pdb2uniprot.uniprot.values]]


## Removing x-ray additives and other solutes
print("removing x-ray additives and other solutes listed in generic_solute_ligands.txt...")
exclusionligs = [ln.strip() for ln in open(f"{HOME}/LigExtract/docs/generic_solute_ligands.txt").readlines()]

rm_trash = []
for exclud in exclusionligs:
    for f in glob(f'{ligands_dir}/*_lig-{exclud}-*'): os.remove(f); rm_trash.append(f)


print("\n".join([", ".join(rm_trash[i:i+5]) for i in range(0, len(rm_trash), 5)]))

# Clean solvent from chain ligands; small ligs already only have their only residue

removed_solvent_res = []
bar1 = Bar('Cleaning all chain ligands... ', max=len(glob(f'{ligands_dir}/*_lig_chain-*.pdb')))
for f in glob(f'{ligands_dir}/*_lig_chain-*.pdb'):
    bar1.next()
    if "_lig_chain-h" in f: continue #bypass ligands that have been assembled based on distance in the previous module
    original_lig = PandasPdb().read_pdb(f)
    #cleaned_lig = subprocess.run(f"obabel -i pdb {f} -r -o pdb", shell=True, capture_output=True)
    resname_counts = original_lig.df["HETATM"].drop_duplicates(["residue_name", "residue_number"]).value_counts("residue_name").to_frame()
    cleaned_lig_het = resname_counts.query("count <= 3").index
    removed_res_solvent_het = original_lig.df['HETATM'][~np.isin(original_lig.df['HETATM'].residue_name, cleaned_lig_het)][["residue_name", "residue_number"]]
    removed_solvent_res.append(removed_res_solvent_het.residue_name.unique())
    original_lig.df['HETATM'] = original_lig.df['HETATM'][np.isin(original_lig.df['HETATM'].residue_name, cleaned_lig_het)]
    original_lig.to_pdb(path=f, records=None, gz=False, append_newline=True) 

bar1.finish()
sys.stderr.write("\ndone cleaning ligands.\n")

try: 
    removed_solvent_res = np.hstack(removed_solvent_res)
except ValueError:
    removed_solvent_res = np.array(removed_solvent_res)

if len(removed_solvent_res)>0: removed_solvent_res = np.unique(removed_solvent_res)

print("\n\nThe following salt/solvent residues were striped from chain ligands:", removed_solvent_res)



list_proteins = [x for x in os.listdir(protein_dir) if x.endswith(".pdb")]
processbar = Bar('Processing all potential ligands... ', max=len(list_proteins))
chains_w_Uniprot_lst = {}
pockets = []


for proteinfile in list_proteins:
    print(proteinfile)
    processbar.next()
    if proteinfile.endswith(".pdb") == False:
        continue
    
    pdbcode = proteinfile.split(".pdb")[0]
    
    print(f"\n############### {pdbcode} ###############\n")

    if len(glob(f'{ligands_dir}/{pdbcode}_*.pdb'))==0:
        print("There are no ligands to consider.")
        continue
    
    try: uniprotid = uniprot2pdb.query(f"pdb == '{pdbcode.upper()}'").uniprot.values[0]
    except IndexError:
        print(f"Unable to find {pdbcode.upper()} in {uniprot2pdb_file}. Make sure this pdb code is in the file, in all capitals.")
        sys.exit(123)
    
    prot = PandasPdb().read_pdb(f'{protein_dir}/{proteinfile}')
    
    ciffile = f"cifs/{pdbcode.lower()}.cif"
    cifdata = CifFileReader().read(ciffile)
    if "_struct_conn" in cifdata[pdbcode.upper()]:
        links = pd.DataFrame.from_dict(cifdata[pdbcode.upper()]["_struct_conn"], orient="index").T
        links = links.query("ptnr1_label_comp_id != 'HOH'")
        links = links.query("conn_type_id == 'covale'")
    else:
        links = pd.DataFrame([])
    
    prd_ligs = [x.split(";")[0].split("chain ")[1] for x in open(f"{ligands_dir}/{pdbcode}_ligand_extraction.log").readlines() if "associated with PRD ID" in x]
    # Get all uniprot that map to the most recent uniprot ID
    # some PDBs have no uniprot associated with any of the chains (e.g. 3in9), so they will not be in *process_uniprot_chains.txt
    # See if they are in SIFTS and only if not, raise error
    fail1 = False
    fail2 = False
    if len(uniprot2pdb_secondary.query(f"pdb == '{pdbcode}'"))==0:
        fail1 = True
        print(f"{protein_dir.split('/')[-1]}_process_uniprot_chains.txt does not include {pdbcode}\n") #. Please re-run process_chains.py
    if len(sifts_pdb2uniprot.query(f"pdbs == '{pdbcode.lower()}'"))==0:#SIFTS
        fail2=True
    if fail1 and fail2:
        print(f"SIFTS mapping also does not include any uniprot mapping to {pdbcode} chains. Abort\n\n")
        sys.exit(f"{pdbcode} has no uniprot mapping, which is not expected (this PDB was retrieved from a uniprot query)")
    
    # only allow chains that map to the query uniprotid
    chains_w_Uniprot1 = uniprot2pdb_secondary.query(f"pdb == '{pdbcode}' and updated_uniprot == '{uniprotid}'").chain.unique()
    chains_w_Uniprot2 = sifts_pdb2uniprot.query(f"pdbs == '{pdbcode.lower()}' and uniprot == '{uniprotid}'").chainName.unique()
    chains_w_Uniprot = np.union1d(chains_w_Uniprot1, chains_w_Uniprot2)
    

    # update dict
    if len(chains_w_Uniprot)> 0: chains_w_Uniprot_lst[pdbcode]=chains_w_Uniprot
    else: chains_w_Uniprot_lst[pdbcode]=[] # to avoid raising the "elementwise comparison failed" warning
    
    allligands = [x for x in os.listdir(ligands_dir) if pdbcode in x and x.endswith(".pdb")]
    
    # locate and remove ions
    cifdata = CifFileReader().read(ciffile)
    ions = []
    if "_chem_comp" in cifdata[pdbcode.upper()]:
        ions = pd.DataFrame.from_dict(cifdata[pdbcode.upper()]["_chem_comp"], orient="index").T
        ions = ions.query("type == 'non-polymer'")
        if len(ions)>0:
            ions = ions[[" ION" in x for x in ions.name.values]].id.values
        else:
            ions = np.array([])
    
    for ion in ions:
        for ionfile in glob(f'{ligands_dir}/{pdbcode}_*_lig-{ion}-*'): 
            os.remove(ionfile) 
            print(f"remove {ion} ions by detecting the word ION in the _chem_comp section")
    
    
    # remove any chains that are other proteins
    # inspect all entities/individual molecules in this PDB using PDBe
    chain_ligs = [x for x in allligands if "lig_chain" in x]

    # get all peptide chains
    peptide_lig_chains = [x.split(":")[-1].strip().split(",") for x in open("ligand_extraction.log").readlines() if x.startswith(f"peptide chains for {pdbcode.upper()}:")]
    if len(peptide_lig_chains)>0: peptide_lig_chains = np.hstack(peptide_lig_chains)
    
    if len(chain_ligs)>0:
        # cross-check with the chain-to-uniprot mapping provided by uniprot
        prot_chains = []
        with gzip.open(f"{HOME}/LigExtract/data/pdb_chain_uniprot.csv.gz") as f:
            for ln in f:
                ln=ln.decode("utf-8").strip().split(",")
                if ln[0]==pdbcode.lower():
                    prot_chains.append(ln[1:3])
        if len(prot_chains)>0: prot_chains = np.array([[c,unp] for c,unp in prot_chains])
        else: prot_chains = np.array([["np.nan", "np.nan"]])
        
        for chainlig in chain_ligs:
            chainlig_chain = chainlig.split("_chain-")[1].split(".pdb")[0]
            if (chainlig_chain in prot_chains[:,0] or chainlig_chain in chains_w_Uniprot_lst[pdbcode]) and chainlig_chain not in prd_ligs:
                chainlig_uniprot = np.unique(prot_chains[prot_chains[:,0]==chainlig_chain,1])
                if chainlig_chain not in peptide_lig_chains:
                    print(f"removing {chainlig} as it is a protein ({','.join(np.unique(chainlig_uniprot))})")
                    os.remove(f'{ligands_dir}/{chainlig}')
                
    # update
    allligands = [x for x in os.listdir(ligands_dir) if pdbcode in x and x.endswith(".pdb")]
    for ligandfile in allligands:
        print(ligandfile)
        if ligandfile not in chain_ligs:
            # inspect whether the ligand is part of the solvent; if so: exclude
            try:
                m = Chem.MolFromPDBFile(f'{ligands_dir}/{ligandfile}')
                if m is None: 
                    print(f"{ligandfile} appears to have some issues in its original structure (as deposited by the author). Please make sure you inspect this ligand if you relly on it having a chemically valid struture.")
                    
                else:
                    smi = Chem.MolToSmiles(m)
                    mols = smi.split(".")
                    mol_count = len(mols)
                    if len(np.unique(mols))==1 and mol_count>1:
                        solv = np.unique(mols)[0]
                        nheavy_solv = sum([x in ascii_uppercase for x in list(solv)])
                        if nheavy_solv<3:
                            print(f"removed {ligandfile} as it has been detected as part of the solvent (multiple molecules with the same structure in the same chain, each molecule with < 3 heavy atoms)")
                            os.remove(f'{ligands_dir}/{ligandfile}')
                            continue
            except OSError: None # Pdb File cannot be interpreted by RDKit
        ligand = PandasPdb().read_pdb(f'{ligands_dir}/{ligandfile}')
        res_hetatm = ligand.df['HETATM'][["residue_name","chain_id","residue_number"]].drop_duplicates().values
        res_hetatm =["{:>3}".format(a) + "{:>2}".format(b) + "{:>4}".format(c) for a,b,c in res_hetatm]
        heav_cnt = (ligand.df["ATOM"].element_symbol!="H").sum() + (ligand.df["HETATM"].element_symbol!="H").sum()
        res_atom = ligand.df['ATOM'][["residue_name","chain_id","residue_number"]].drop_duplicates().values
        res_atom =["{:>3}".format(a) + "{:>2}".format(b) + "{:>4}".format(c) for a,b,c in res_atom]

        ligand = pd.concat([ligand.df["HETATM"], ligand.df["ATOM"]])
        
        # Detect groups of multi-residue ligands // this only works for same-chain residues
        if len(res_hetatm)>1 and len(links)>0:
            resgroups = groupresidues(res_hetatm, links)
        else:
            #len hetatm is 0,1
            # can either be a 1-res small lig alone or a chain lig
            resgroups = [x[-4:].strip() for x in res_hetatm + res_atom]
            resgroups = pd.DataFrame(zip(resgroups,range(1,len(resgroups)+1)), columns=["resn","group"])
            resgroups.loc[:,"resn"] = resgroups.resn.values.astype(int)

        # accomodate for chain ligands
        
        # go through all ligands
        contacts_chainlig = []
        contacts_chainlig_reschain = []
        for resgroup in resgroups.group.unique():
            close_res = []
            close_res_chain = []
            for res in resgroups.query(f'group=={resgroup}').resn.values:
                ligres = ligand[ligand.residue_number == res]
                for ref_coord in ligres[['x_coord', 'y_coord', 'z_coord']].values:
                    keepclose = prot.distance(ref_coord, records=('ATOM',))<distmax
                    closeatomschain = prot.df["ATOM"][keepclose].residue_name + prot.df["ATOM"][keepclose].residue_number.astype(str) + "-" + prot.df["ATOM"][keepclose].chain_id
                    closeatoms = prot.df["ATOM"][keepclose].residue_name + prot.df["ATOM"][keepclose].residue_number.astype(str)
                    # Remove self contact
                    if "_lig_chain-" in ligandfile:
                        chain_lig = ligandfile.split("_lig_chain-")[1].split(".")[0]
                        lines2keep = prot.df["ATOM"][keepclose].chain_id != chain_lig
                        closeatomschain = closeatomschain[lines2keep]
                        closeatoms = closeatoms[lines2keep]
                    close_res.append(closeatoms.drop_duplicates())
                    close_res_chain.append(closeatomschain.drop_duplicates())
            # store contacts per residue
            close_res = np.unique(np.hstack(close_res))
            close_res_chain = np.unique(np.hstack(close_res_chain))
            # if ligand is a chain the contacts in the current residue must be appended to a larger list
            contacts_chainlig.append(close_res)
            contacts_chainlig_reschain.append(close_res_chain)
        
        if len(contacts_chainlig)>0: contacts_chainlig = np.unique(np.hstack(contacts_chainlig))
        else: contacts_chainlig = np.array([])
        if len(contacts_chainlig_reschain)>0: contacts_chainlig_reschain = np.unique(np.hstack(contacts_chainlig_reschain))
        else: contacts_chainlig_reschain = np.array([])
        if len(contacts_chainlig)==0: 
            pockets.append([ligandfile,np.nan,"","",0,0,pdbcode,""])
        else:
            chains = ";".join(np.unique([x.split("-")[-1] for x in contacts_chainlig_reschain]))
            pockets.append([ligandfile,np.nan,";".join(contacts_chainlig_reschain), ";".join(contacts_chainlig), 
                            len(contacts_chainlig_reschain),len(contacts_chainlig), pdbcode, chains])




print ("\n\n--------- Processing the pockets ---------\n")

pockets = pd.DataFrame(pockets)

allres = []
for pocket in pockets.values:
    if pocket[4]==0:
        print(f"There are no closeby residues for {pocket[0]} @ {pocket[-2]}")
        continue
    allres+=[x.split("-")[0] for x in pocket[2].split(";") if x.split("-")[1] in chains_w_Uniprot_lst[pocket[6]]]

if len(allres)==0:
    print("there are no contact residues between the identified ligands and the respective target protein in each PDB. Please check your Uniprot code(s)")
    sys.exit(123)

allres = np.hstack(allres)
pockets.columns = ["ligandfile","resnumber","pocketres_chain","pocketres","pocketres_chain_size","pocketres_size","pdbcode","chain_name"]

# At least one chain of desired uniprot ID must be in contact with the ligand
acceptpockets = []
for pdbi,ci in pockets[["pdbcode","chain_name"]].values:
    okchains = np.intersect1d(ci.split(";"),chains_w_Uniprot_lst[pdbi])
    if len(okchains)>0:
        acceptpockets.append(True)
    else:
        acceptpockets.append(False)
        


rejected = pockets[~np.array(acceptpockets)]
for ligi,pdbi,ci in rejected[["ligandfile","pdbcode","chain_name"]].values:
    print(f"{ligi} excluded: bound to chain(s) {ci} which does not coincide with chains of the desired protein {';'.join(chains_w_Uniprot_lst[pdbi])}")

pockets = pockets[acceptpockets]
if len(pockets) == 0:
    print(f"There are no pockets where the ligand is interacting with the protein {uniprotid}")
    sys.exit(123)

print ("\n\n--------- Solving cases of multiple similar ligands in the same chain (which appear as one ligand record) ---------\n")

# Untie cases of same ligand, same chain according to the largest pocket. Clean the ligand file accordingly
resnum_count = pockets.groupby(["ligandfile"])["resnumber"].nunique()
ligs2untie = resnum_count[resnum_count>1].reset_index()
print (f"\n {len(ligs2untie)} ligands to untie \n")

# reset indices to allow safely removing instances by index ID
pockets.index = range(len(pockets))

ligs2drop = []
for lig,n_repeat in ligs2untie.values:
    lig_lines = pockets.query(f"ligandfile == '{lig}'")
    lig_lines = lig_lines.sort_values("pocketres_chain_size", ascending=False)
    if keeprepeats=="n":
        print("(Keeping only the best ligand from the same ligand ID)")
        selected_lig = lig_lines[0:1]
        rejected_ligs = lig_lines[1:]
        ligcode = lig.split('lig-')[1].split('.')[0]
        rejected_abbrev = ", ".join([f'{ligcode}-res{int(x)}' for x in rejected_ligs.resnumber])
        print(f"{lig} contains {n_repeat} distinct ligands: {ligcode}-res{int(selected_lig.resnumber.values[0])}  will be SELECTED for having the maximum number of pocket residues {selected_lig.pocketres_chain_size.values[0]}.")
        print(f"                                                      {rejected_abbrev}  will be REJECTED.")
        pockets = pockets.drop(rejected_ligs.index, axis="index")
        # Correct the ligand to contain just the desired residue
        ligand = PandasPdb().read_pdb(f'{ligands_dir}/{lig}')
        ligand.df["HETATM"] = ligand.df["HETATM"].query(f'residue_number == {int(selected_lig.resnumber.values[0])}')
        ligand.to_pdb(path=f'{ligands_dir}/{lig}', records=None, gz=False, append_newline=True)
    if keeprepeats=="y":
        print("(Keeping all ligands from the same ligand ID)")
        for i in range(len(lig_lines)): # one line per residue
            selected_lig = lig_lines.iloc[i]
            #print(selected_lig)
            ligand = PandasPdb().read_pdb(f'{ligands_dir}/{lig}')
            r=int(selected_lig.resnumber)
            ligand.df["HETATM"] = ligand.df["HETATM"].query(f'residue_number == {r}')
            print(f'created ligand file:',f'{lig.replace("_lig",f"_res{r}_lig")}')
            ligand.to_pdb(path=f'{ligands_dir}/{lig.replace("_lig",f"_res{r}_lig")}', records=None, gz=False, append_newline=True)
            selected_lig = lig_lines.iloc[[i]].replace({lig:lig.replace("_lig",f"_res{r}_lig")})
            ligs2drop.append(lig)
            pockets = pd.concat([pockets,selected_lig])

pockets = pockets[np.isin(pockets.ligandfile, ligs2drop, invert=True)]
pockets.to_csv(f'pockets_{protein_dir.split("/")[-1]}.txt', sep="\t", index=False)

# final clean-up
curr_ligands = [x for x in os.listdir(ligands_dir) if x.endswith(".pdb")]
for l in curr_ligands:
    if l not in pockets.ligandfile.values: os.remove(f'{ligands_dir}/{l}')


pockets_file = f'pockets_{protein_dir.split("/")[-1]}.txt'
sys.stderr.write(f"\n\nFinished. All ligands and their pockets are listed in {pockets_file}\n\n")
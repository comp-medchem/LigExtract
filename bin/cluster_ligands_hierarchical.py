import pandas as pd
import numpy as np
from biopandas.pdb import PandasPdb
import sys
import os
import linecache
from time import sleep
import subprocess
from progress.bar import ChargingBar, Bar
import argparse
from glob import glob
import shlex
from itertools import combinations
from sklearn.cluster import AgglomerativeClustering
from scipy.spatial.distance import euclidean
from copy import copy
from pathlib import Path
HOME = str(Path.home())

parser = argparse.ArgumentParser(description='Final procedure do decide which ligands are the most likely in a given protein, using different criteria that looks into ligand structure and annotation, ligand placement in the pocket, etc.')
parser.add_argument('--pdbPath', type=str, required=True, dest="prot_dir", help='Path to directory containing all PDBs to process.')
parser.add_argument('--ligandsPath', type=str, required=True, dest="lig_dir", help='Path to directory with extracted ligand files.')
parser.add_argument('--prdCif', type=str, required=True, dest="cif_file", help='CIF File with all ligands in BIRD')
parser.add_argument('--uniprot2pdbFile', type=str, required=True, dest="uniprot2pdbFile", help='File containing the PDB-to-UniprotID mapping')
args = parser.parse_args()

prot_dir = args.prot_dir
lig_dir = args.lig_dir
cif_file = args.cif_file
uniprot2pdbFile = args.uniprot2pdbFile

uniprot2pdbFile = pd.read_csv(uniprot2pdbFile, sep="\t")

if prot_dir[-1]=="/": prot_dir = prot_dir[:-1]
if lig_dir[-1]=="/": lig_dir = lig_dir[:-1]

print("\n------------------  Filtering Ligands  ------------------\n")


pockets = pd.read_csv(f"pockets_{prot_dir.split('/')[-1]}.txt", sep="\t")
oligos_file = open(f"{lig_dir}_oligosaccharides.txt","w")
oligos_file.write("WARNING! Some protein residues might be added here due to the nature of the N-glycosylation group detection procedure.\n")
oligos_file.write("         However this has no impact on the whole procedure as this list is used to exclude ligands for being N-glycosylation\n")
oligos_file.write("         groups. Even if a protein residue is within the list of ligands (i.e. a modified residue in HETATM),\n")
oligos_file.write("         it would simply remove that supposed ligand (correctly, as it is actually not a ligand)\n\n\n")

pdbs_in_pockets = pockets.pdbcode.drop_duplicates()

save_clean_pockets_list = []

bar = Bar('Filtering ligands in each PDB... ', max=len(pdbs_in_pockets))


for pdb in pdbs_in_pockets:
    #bar.next()
    print(f"\n############# {pdb} ###########\n")
    
    if len(glob(f'{lig_dir}/{pdb}_*.pdb'))==0:
        print("There are no ligands to consider.")
        continue
    # Remove ligands that are oligosacharides
    cif_f = open(f"cifs/{pdb}.cif").read()
    if "_pdbx_entity_branch_list" in cif_f: 
        cmpd_where = np.argwhere([".comp_id" in x for x in cif_f.split("\n") if x.startswith("_pdbx_entity_branch_list.")]).flatten()[0]
        oligos = cif_f.split("_pdbx_entity_branch_list.hetero \n")[1].split("#")[0]
        oligos = [x for x in oligos.split("\n")[:-1]]
        oligos = np.unique([shlex.split(x)[cmpd_where] for x in oligos])
        if len(oligos)>0: 
            o='\t'.join(oligos)
            oligos_file.write(f"{pdb}\t{o}\n")
    else: oligos = []
    
    
    if "_struct_conn.pdbx_role" in cif_f:
        block = [x for x in cif_f.split("# \n") if "_struct_conn." in x and x.split("\n")[1].startswith("_struct_conn.")][0].split("\n")[:-1]
        if block[0]=="loop_":
            headers = [x.replace("_struct_conn.","").strip() for x in block if x.startswith("_struct_conn")]
            body = [x for x in block[1:] if x.startswith("_struct_conn") == False]
            join_ln = np.argwhere([len(shlex.split(x))<len(headers) for x in body]).flatten()
            for x in join_ln:
                # if current line is empty, this is a continuation line: bypass.
                if len(body[x].split())==0: continue
                body[x] = body[x]+body[x+1]
                body[x+1]=""
            body = [shlex.split(x) for x in body if "Glycosylation" in x]
            glycosylation_het = pd.DataFrame(body, columns=headers)
            add_oligos = glycosylation_het[["ptnr1_label_comp_id","ptnr2_label_comp_id"]].values.flatten() # "pdbx_role"
        else:
            glycosylation_het = np.array([x.strip().split() for x in block])
            glycosylation_het = pd.DataFrame(glycosylation_het[:,1].reshape(1,-1), columns = [x.replace("_struct_conn.","") for x in glycosylation_het[:,0]])
            if "Glycosylation" not in glycosylation_het.pdbx_role.values[0]:
                add_oligos = []
            else:
                add_oligos = glycosylation_het[["ptnr1_label_comp_id","ptnr2_label_comp_id"]].values.flatten()
        if len(add_oligos)>0: 
            # ! some of these will be residues but this will have no impact on the process
            # since residues will never be a possible individual ligand to begin with.
            o='\t'.join(add_oligos)
            oligos_file.write(f"{pdb}\t{o}\n")
    
    else: add_oligos = []
    
    oligos = list(oligos) + list(add_oligos)
    
    del cif_f 
    # get all identified ligs
    logfile = [ln.strip() for ln in open(f"{lig_dir}/{pdb}_ligand_extraction.log").readlines()][5:]
    chain2prdIDs = np.array([[x.split("rebuilt ligand detected in chain ")[1][0],x.split("PRD ID ")[1].split(" ")[0]] for x in logfile if "associated with PRD ID" in x])
    has_rebuilt_lig = any(["lig_chain-" in x for x in os.listdir(lig_dir) if pdb in x])
    rare_ligand = [x.split()[0] for x in logfile if "ligand " in x and "Might not be a ligand." not in x and "rebuilt" not in x]
    
    # clean all non-rare ligands
    currentfiles = [x for x in os.listdir(lig_dir) if pdb in x and x.endswith("pdb") and "_lig-" in x]
    trash_ligands = [x.split()[0] for x in logfile if "ligand " in x and "Might not be a ligand." in x]
    if len(currentfiles)>0 and len(trash_ligands)>0:
        print("Clean all ligands that are likely to be crystallographic additives/reagents:")
        for ligkey in trash_ligands: 
            for filex in currentfiles:
                if ligkey in filex:
                    print("\t",filex)
                    os.remove(f'{lig_dir}/{filex}')
                    pockets = pockets[~np.in1d(pockets.ligandfile, filex)]
        print("-----------------------------\n")
    
    # update rare ligands from the log file
    currentfiles = [x for x in os.listdir(lig_dir) if pdb in x and x.endswith("pdb")]
    currentfiles = [x.split("_lig-")[-1].split(".")[0] for x in currentfiles]
    if np.intersect1d(currentfiles,oligos).shape[0]>0: print("The current ligands are glycosylation groups, and will be excluded:", " ".join(np.intersect1d(currentfiles,oligos)))
    for x in np.intersect1d(currentfiles,oligos): 
    	for i in glob(f'{lig_dir}/{pdb}*-{x}.pdb'):
    		os.remove(i)
    currentfiles = np.setdiff1d(currentfiles, oligos)
    rare_ligand = np.intersect1d(rare_ligand,currentfiles)
    print("small ligand(s):", rare_ligand)
    
    # if has_rebuilt_lig is True and len(rare_ligand)>0, check if rebuilt lig contains rare_lig
    if has_rebuilt_lig == True and len(rare_ligand)>0:
        rebuilts_lig_files = [x for x in os.listdir(lig_dir) if "_lig_chain-" in x and pdb in x]
        files2save = []
        for f in rebuilts_lig_files:
            lig = PandasPdb().read_pdb(f"{lig_dir}/{f}")
            all_res_cnt = len(lig.df["ATOM"].residue_number.drop_duplicates()) + len(lig.df["HETATM"].residue_number.drop_duplicates())
            all_res = np.hstack([lig.df["ATOM"].residue_name.unique(), lig.df["HETATM"].residue_name.unique()])
            rareLigs_in_rebuiltLigs = np.intersect1d(rare_ligand,all_res)
            if all_res_cnt>10: ligtype = "large"
            if all_res_cnt<=10: ligtype = "small"
            chain = "chain-"+f.split("chain-")[1][0]
            if len(rareLigs_in_rebuiltLigs) > 0:
                decision="rebuilt"
                files2remove = np.hstack([[x for x in os.listdir(lig_dir) if (l in x and chain in x and x.startswith(pdb))] for l in rareLigs_in_rebuiltLigs])
                print("remove", "; ".join(files2remove), f"--> contained in {f}")
                files2save.append(f)
                for x in files2remove: os.remove(f'{lig_dir}/{x}')
                # update pockets file
                pockets = pockets[~np.in1d(pockets.ligandfile, files2remove)]
                remaining_ligs = [x for x in os.listdir(lig_dir) if pdb in x and chain in x and x!=f and x.endswith(".pdb")]
                # update pockets file
                pockets = pockets[~np.in1d(pockets.ligandfile, remaining_ligs)]
                
    
    if has_rebuilt_lig == True:
        # cross ref with PRDs
        prd2pdb = pd.read_csv(f"{HOME}/LigExtract/data/prd_to_pdb_IDs.txt", sep="\t")
        prd2pdb = prd2pdb.query(f"pdb == '{pdb}'").prd_ID.drop_duplicates().values
        rebuilt_ligs = [x for x in os.listdir(lig_dir) if pdb in x and "lig_chain-" in x]
        if len(prd2pdb)==0:
            print(f"There are {len(rebuilt_ligs)} outstanding ligands with no further annotation in BIRD.")
            for a in rebuilt_ligs: print("\t"+a)
            found_match = None
        else:
            print("Searching for matches in BIRD...")
            # try to use what was gathered first
            if len(chain2prdIDs)>0: lig2search_bird = [x for x in rebuilt_ligs if x.split("chain-")[1].split(".")[0] not in chain2prdIDs[:,0]]
            else: lig2search_bird = []
            # save the already found macthes from the log
            if len(chain2prdIDs)>0: ligfound_bird = [[x,x.split("chain-")[1].split(".")[0]] for x in rebuilt_ligs if x.split("chain-")[1].split(".")[0] in chain2prdIDs[:,0]]
            else: ligfound_bird = []
            
            for l,ch in ligfound_bird:
                p = chain2prdIDs[chain2prdIDs[:,0]==ch][:,1][0]
                print(f"ligand  {l}  is annotated in BIRD under {p}  (full sequence match). This is likely the ligand of focus in this PDB entry.")
                
            
            if len(lig2search_bird)>0:
                print("searching ligands that weren't found in ligand extraction...")
                prd_block = []
                titles = []
                c = 0
                with open(cif_file) as prd:
                    for line in prd:
                        if line.startswith("data_PRD_"):
                            prd_block.append(c)
                            titles.append(line.strip().replace("data_",""))
                        c +=1
                
                prd_block = np.array(list(zip(prd_block,prd_block[1:]+[c])))
                get_blocks = prd_block[np.in1d(titles, prd2pdb)]
                
                seq_list = []
                for blockStart,blockStop in get_blocks:
                    block = [linecache.getline(cif_file,x) for x in range(blockStart, blockStop)]
                    obs_seq = "".join(block).split('_pdbx_reference_entity_poly_seq.observed \n')[1].split("# \n")[0]
                    obs_seq = [x.split() for x in obs_seq.split("\n")]
                    obs_seq = pd.DataFrame(obs_seq[0:-1])
                    headers = "".join(block).split("loop_\n_pdbx_reference_entity_poly_seq.")[1].split(".observed \n")[0] + ".observed"
                    headers = headers.replace("_pdbx_reference_entity_poly_seq.","").split(" \n")
                    obs_seq.columns = headers
                    full_seq = obs_seq.mon_id.values
                    exp_seq = obs_seq.query("observed == 'Y'").mon_id.values
                    name = obs_seq.prd_id.drop_duplicates().values[0]
                    seq_list.append([name," ".join(full_seq), " ".join(exp_seq)])
                
                pdbseqs = [ln.strip() for ln in open(f'{prot_dir}/{pdb}.pdb').readlines() if ln.startswith("SEQRES ")]
                found_match = []
                for lig in rebuilt_ligs:
                    lig_chain = lig.split("lig_chain-")[1].split(".")[0]
                    lig_seq_match = []
                    for seq in pdbseqs:
                        if lig_chain in seq[0:13]:
                            lig_seq_match.append(" ".join(seq[12:].split()[1:]))
                    lig_seq_match = " ".join(lig_seq_match)
                    for prd_seq in seq_list:
                        if lig_seq_match == prd_seq[1]:
                            found_match.append([lig,prd_seq[0]])
                
                if len(found_match) == 0: print("None of the rebuilt-chain ligands are in BIRD.")
                else: 
                    for ln in found_match: 
                        print(f"ligand  {ln[0]}  is annotated in BIRD under {ln[1]}  (full sequence match). This is likely the ligand of focus in this PDB entry.")
                        
    
    
    # First remove cross-contacts between ligands; assuming all surviving ligands at this stage are real ligands
    pdb_pockets = pockets.query(f"pdbcode == '{pdb}'")
    allchainligs = pdb_pockets.ligandfile.unique()
    save_clean_pockets = []
    for lig_line,pocketres_chain in pdb_pockets[["ligandfile","pocketres_chain"]].values:
        chain2rm = [x.split("lig_chain-")[1].split(".")[0] for x in np.setdiff1d(allchainligs,lig_line) if "lig_chain-" in x]
        pocketres_chain = [x for x in pocketres_chain.split(";") if x.split("-")[1] not in chain2rm]
        pocketres_chain_size = len(pocketres_chain)
        chain_name = ";".join(np.unique([x.split("-")[1] for x in pocketres_chain]))
        pocketres_chain = ";".join(pocketres_chain)
        save_clean_pockets.append([lig_line, pocketres_chain, pocketres_chain_size, chain_name])
    
    
    save_clean_pockets = pd.DataFrame(save_clean_pockets, columns = ["ligandfile","pocketres_chain", "pocketres_chain_size", "chain_name"])
    save_clean_pockets["ligtype"] = np.where(["lig_chain" in x for x in save_clean_pockets.ligandfile], "chain ligand", "small-molecule ligand")
    save_clean_pockets["lig_ID"] = [x.split("lig-")[-1].split(".")[0] for x in save_clean_pockets.ligandfile]
    save_clean_pockets = save_clean_pockets[~np.in1d(save_clean_pockets.lig_ID, oligos)]
    
    if len(oligos)>0: print(f"The following oligosacharide residues were found and will be excluded from the ligands: {', '.join(oligos)}")
    
    
    # inspect the size of the pocket
    sparse_pockets = save_clean_pockets.query("pocketres_chain_size<6").ligandfile
    if len(sparse_pockets)>0:
        print("WARNING! There are ligands with very sparse pockets [<6 close res] (which indicates this might not be a ligand). These will be excluded")
        for x in sparse_pockets: 
            print("\t"+x)
            os.remove(f"{lig_dir}/{x}")
        save_clean_pockets = save_clean_pockets[np.in1d(save_clean_pockets.ligandfile, sparse_pockets, invert=True)]
    
    save_clean_pockets_list.append(save_clean_pockets)


oligos_file.close()

# Clean ligands to keep only the final selected list
all_ligands = [x for x in os.listdir(lig_dir) if x.endswith(".pdb")]


# save clean pockets, updated with the selected ligands
save_clean_pockets = pd.concat(save_clean_pockets_list, sort=False)
save_clean_pockets.loc[:,"pdbcode"] = [x.split("_")[0] for x in save_clean_pockets.ligandfile]

save_clean_pockets.to_csv(f"cleanpockets_{prot_dir.split('/')[-1]}.txt", sep="\t", index=False)


print("\n------------------  Clustering Ligands  ------------------\n")
subprocess.run(f"rm -Rf {prot_dir}/pdbs_filtered_chains; mkdir {prot_dir}/pdbs_filtered_chains", shell=True)

prot_lst = uniprot2pdbFile.uniprot.unique()

sys.stderr.write("\n")

for p_i, prot in enumerate(prot_lst):
    print(f'**** Protein {prot}')
    pdbs = uniprot2pdbFile.query(f"uniprot == '{prot}'").pdb.str.lower().values
    pockets_prot = save_clean_pockets[np.in1d(save_clean_pockets.pdbcode, pdbs)]
    if len(pockets_prot) == 0:
        sys.stderr.write(f"\n\n*** Protein {prot} ({p_i+1}/{len(prot_lst)}) : bypass as it has no found ligands in any PDBs.\n")
        print(f"bypass {prot} as it has no found ligands in any PDBs.")
        continue
    # reset folder in case this was already run before for this protein
    subprocess.run(f"rm -Rf {prot_dir}/pdbs_filtered_chains/{prot}; mkdir {prot_dir}/pdbs_filtered_chains/{prot}", shell=True)
    sys.stderr.write(f"\n\n*** protein {prot} ({p_i+1}/{len(prot_lst)}) ***\n")
    bar = Bar("Aligning PDBs ... ", max=pockets_prot.pdbcode.nunique()) 
    
    # prior to alignment chains need to be filtered and split when appropriate
    for pdb in pockets_prot.pdbcode.unique():
        bar.next()
        pdb_ligs = pockets_prot.query(f"pdbcode == '{pdb}'")
        ligs_inmultichain = [";" in x for x in pdb_ligs.chain_name.values]
        
        # Split chains to alow self alignment; only keep multiple chains if at least one ligand has then simultaneouly
        allcontactchains = pdb_ligs.chain_name.unique()
        for chainset in allcontactchains:
            chainset = chainset.split(";")
            #ligchains = [x.split("chain-")[1].split(".")[0] for x in pdb_ligs.ligandfile.values if "lig_chain" in x]
            # get ligands that have all its chains in the current chainset (i.e. a ligand with chain K will be included in chainset K;L)
            ligchains = [np.in1d(x.split(";"),chainset).all() for x in pdb_ligs.chain_name]
            ligchains = [x.split("chain-")[1].split(".")[0] for x in pdb_ligs[ligchains].ligandfile.values]
            allcontactchains = np.hstack([chainset,ligchains])
            protein_pdb = PandasPdb().read_pdb(f"{prot_dir}/{pdb}.pdb")
            protein_pdb.df["ATOM"] = protein_pdb.df["ATOM"][np.in1d(protein_pdb.df["ATOM"].chain_id, allcontactchains)]
            protein_pdb.df["HETATM"] = protein_pdb.df["HETATM"][np.in1d(protein_pdb.df["HETATM"].chain_id, allcontactchains)]
            protein_pdb.to_pdb(path=f"{prot_dir}/pdbs_filtered_chains/{prot}/{pdb}_keychain{'-'.join(chainset)}.pdb", records=None, gz=False, append_newline=True)

    rmsd=subprocess.run(f'pymol -cq ~/LigExtract/bin/align_pdbs_pockets.py -- {prot_dir}/pdbs_filtered_chains/{prot}', shell=True, capture_output=True)
    refpdb_align = [x for x in os.listdir(f'{prot_dir}/pdbs_filtered_chains/{prot}') if x.endswith(".pdb")][0]
    if len(pockets_prot.pdbcode.unique())>1:
        rmsd = rmsd.stdout.decode("utf=8")
        rmsd = eval("["+rmsd[rmsd.find("["):].strip().replace("\n",",")+"]")
        print(f"RMSDs from alignment against {refpdb_align}:")
        for i in rmsd:
        	print("\t",i)
        bad_align = [x for x in rmsd if x[1]>3.5]
    else:
        rmsd = []
        bad_align = []
    if len(bad_align)>0:
        print(f"After aligning {len(rmsd)+1} structures using {refpdb_align} as reference, {len(bad_align)} showed poor alignment to the rest. Please consider manual alignment and inspection of the following structures:")
        for x in bad_align:
            print("\t",x)
    else:
        print(f"Successful alignment of {len(rmsd)+1} structures")
    
    ligand_centroid = []
    sys.stderr.write("\nBuilding pockets clusters...\n")
    for pdbalign in os.listdir(f"{prot_dir}/pdbs_filtered_chains/{prot}/aligned_pdbs"):#pockets_prot.pdbcode.unique():
        print("PROTEIN: ", pdbalign)
        #bypass = False
        pdb = pdbalign.split("_")[0]
        if pdbalign.replace(".pdb","") in bad_align: continue #bypass=True; ; break
        #if bypass == True:
        #    continue
        protein_pdb = PandasPdb().read_pdb(f"{prot_dir}/pdbs_filtered_chains/{prot}/aligned_pdbs/{pdbalign}")
        pdb_ligs = pockets_prot.query(f"pdbcode == '{pdb}'")
        print(pdb_ligs.ligandfile.values)
        for lig_pdb in pdb_ligs.ligandfile.values:
            print(lig_pdb, "-----------------")
            lig_pdb_file = PandasPdb().read_pdb(f"{lig_dir}/{lig_pdb}")
            lig_resn = np.hstack([lig_pdb_file.df["ATOM"].residue_number.drop_duplicates().values, lig_pdb_file.df["HETATM"].residue_number.drop_duplicates().values])
            save_aligned_lig = PandasPdb().read_pdb(f"{prot_dir}/pdbs_filtered_chains/{prot}/aligned_pdbs/{pdbalign}")
            if "_lig-" in lig_pdb:
                resname = lig_pdb.split("lig-")[1].split(".")[0]
                chainid = lig_pdb.split("chain-")[1].split("_")[0]
                atom_coord = protein_pdb.df["HETATM"].query(f"residue_name == '{resname}' and chain_id == '{chainid}' and residue_number == {lig_resn[0]}")
                atom_coord = atom_coord[['x_coord', 'y_coord', 'z_coord']]
                save_aligned_lig.df["HETATM"] = save_aligned_lig.df["HETATM"].query(f"residue_name == '{resname}' and chain_id == '{chainid}' and residue_number == {lig_resn[0]}")
                save_aligned_lig.df["ATOM"] = save_aligned_lig.df["ATOM"].query(f"residue_name == 'NORESIDUES'")
            
            if "lig_chain" in lig_pdb:
                chainid = lig_pdb.split("_chain-")[1].split(".")[0]
                atom_coord = protein_pdb.df["ATOM"].query(f"chain_id == '{chainid}'")
                atom_coord = atom_coord[np.in1d(atom_coord.residue_number, lig_resn)]
                atom_coord = atom_coord[['x_coord', 'y_coord', 'z_coord']]
                save_aligned_lig.df["ATOM"] = save_aligned_lig.df["ATOM"].query(f"chain_id == '{chainid}'")
                save_aligned_lig.df["ATOM"] = save_aligned_lig.df["ATOM"][np.in1d(save_aligned_lig.df["ATOM"].residue_number,lig_resn)]
                
                hetatom_coord = protein_pdb.df["HETATM"].query(f"chain_id == '{chainid}'")
                hetatom_coord = hetatom_coord[np.in1d(hetatom_coord.residue_number, lig_resn)]
                hetatom_coord = hetatom_coord[['x_coord', 'y_coord', 'z_coord']]
                atom_coord = pd.concat([hetatom_coord,atom_coord])
                save_aligned_lig.df["HETATM"] = save_aligned_lig.df["HETATM"].query(f"chain_id == '{chainid}'")
                save_aligned_lig.df["HETATM"] = save_aligned_lig.df["HETATM"][np.in1d(save_aligned_lig.df["HETATM"].residue_number,lig_resn)]
            
            liglines = save_aligned_lig.df["HETATM"].shape[0]+save_aligned_lig.df["ATOM"].shape[0]
            if liglines == 0:
                continue
            newligname = lig_pdb.replace('.pdb','_aligned_LIG.pdb')
            # if we are considering different sets of chains, we must avoid overwriting a previous ligand file with the same name
            if newligname in os.listdir(f'{prot_dir}/pdbs_filtered_chains/{prot}/aligned_pdbs/'):
                newligname = newligname.replace("_aligned_", "_1_aligned_")
            save_aligned_lig.to_pdb(path=f"{prot_dir}/pdbs_filtered_chains/{prot}/aligned_pdbs/{newligname}", records=None, gz=False, append_newline=True)
            centroid = atom_coord.mean().values.round(3)
            print(newligname, centroid)
            ligand_centroid.append(np.hstack([newligname, centroid]))
    
    # Cluster controids
    ligand_centroid = pd.DataFrame(ligand_centroid)
    if ligand_centroid[0].nunique() != len(ligand_centroid):
        print("WARNING!!! Some ambiguities detected in the ligands that will be used in the clustering. Please inspect inspect_ligand_centroids.txt")
        ligand_centroid.to_csv("inspect_ligand_centroids.txt", sep="\t", index=False)
    centroids_coords = ligand_centroid[[1,2,3]].astype(float)#.round(3)
    if len(centroids_coords) > 1:
        clustering = AgglomerativeClustering(n_clusters=None, distance_threshold=10).fit(centroids_coords)
        cluster_classes = clustering.labels_
        clust_ec = {cl:max([round(euclidean(i,j),3) for i,j in combinations(centroids_coords[cluster_classes==cl].values,2)]) if sum(cluster_classes==cl)>1 else 0 for cl in np.unique(cluster_classes)}
    else:
        cluster_classes = np.array([0])
        clust_ec = {0:0}
    sys.stderr.write(f"--- {len(np.unique(cluster_classes))} different clusters were found.\n\n")
    print(f"--- {len(np.unique(cluster_classes))} different clusters were found")
    
    ligand_centroid.columns = ["ligandfile","centroidX","centroidY", "centroidZ"]
    ligand_centroid.loc[:,"cluster"] = cluster_classes
    sorted_clusters = ligand_centroid["cluster"].value_counts().sort_values(ascending=False).index
    ligand_centroid_sorted = []
    for s in sorted_clusters:
        l = copy(ligand_centroid.query(f"cluster == {s}"))
        l.loc[:,"maxdist"]=clust_ec[s]
        ligand_centroid_sorted.append(l)
    
    ligand_centroid_sorted = pd.concat(ligand_centroid_sorted)
    ligand_centroid_sorted.to_csv(f"{prot}_pockets_hierarch-clusters.txt", sep="\t", index=False)
    
    # Produce images for clustered pockets
    for c in ligand_centroid_sorted.cluster.unique():
        pocket_ligs = ligand_centroid_sorted.query(f"cluster == '{c}'").ligandfile.values
        pocket_ligs = " ".join(pocket_ligs)
        subprocess.run(f"pymol -cq ~/LigExtract/bin/align_ligs_figures.py -- pocketcluster{c} {prot_dir}/pdbs_filtered_chains/{prot}/aligned_pdbs {pocket_ligs}", shell=True, capture_output=True)
    
    sys.stderr.write(f"Clustered pockets have been stored in {prot}_pockets_hierarch-clusters.txt. Pictures of the different pocket clusters created in {prot_dir}/pdbs_filtered_chains/{prot}/aligned_pdbs")
    
    # Produce one global image where each pocket has a color
    subprocess.run(f"pymol -cq ~/LigExtract/bin/allpockets_figure.py -- {prot_dir}/pdbs_filtered_chains/{prot}/aligned_pdbs {prot}_pockets_hierarch-clusters.txt", shell=True, capture_output=True)

sys.stderr.write(f"\n\nFinished Clustering pockets")
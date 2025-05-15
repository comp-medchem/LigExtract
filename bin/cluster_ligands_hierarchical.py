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
from concurrent.futures import ProcessPoolExecutor, as_completed
from pdbecif.mmcif_io import CifFileReader
#HOME = str(Path.home())
HOME = os.path.realpath(__file__)
HOME = HOME.split("/LigExtract")[0]

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


length = 90; pad_char = '-'

title="Filtering Ligands"
padding_total = length - len(title) - 2
print(f"\n{pad_char * (padding_total // 2)} {title} {pad_char * (padding_total - (padding_total // 2))}\n")


pockets = pd.read_csv(f"pockets_{prot_dir.split('/')[-1]}.txt", sep="\t")
oligos_file = open(f"{lig_dir}_oligosaccharides.txt","w")
oligos_file.write("WARNING! Some protein residues might be added here due to the nature of the N-glycosylation group detection procedure.\n")
oligos_file.write("         However this has no impact on the whole procedure as this list is used to exclude ligands for being N-glycosylation\n")
oligos_file.write("         groups. Even if a protein residue is within the list of ligands (i.e. a modified residue in HETATM),\n")
oligos_file.write("         it would simply remove that supposed ligand (correctly, as it is actually not a ligand)\n\n\n")

pdbs_in_pockets = pockets.pdbcode.drop_duplicates()

prd2pdb_data = pd.read_csv(f"{HOME}/LigExtract/data/prd_to_pdb_IDs.txt", sep="\t")

save_clean_pockets_list = []

bar = Bar('Filtering ligands in each PDB... ', max=len(pdbs_in_pockets))


for pdb in pdbs_in_pockets:
    bar.next()
    print(f"\n############# {pdb} ###########\n")
    
    if len(glob(f'{lig_dir}/{pdb}_*.pdb'))==0:
        print("There are no ligands to consider.")
        continue

    # Remove ligands that are oligosacharides
    cifdata = CifFileReader().read(f"cifs/{pdb}.cif")
    data = cifdata[pdb.upper()]
    if "_pdbx_entity_branch_list" in data:
        data = pd.DataFrame.from_dict(data["_pdbx_entity_branch_list"], orient="index").T
        oligos = data.comp_id.unique()
        if len(oligos)>0: 
            o='\t'.join(oligos)
            oligos_file.write(f"{pdb}\t{o}\n")
    else: oligos = []
    
    
    data = cifdata[pdb.upper()]
    if "_struct_conn" in data:
        data = pd.DataFrame.from_dict(data["_struct_conn"], orient="index").T
        data = data[["Glycosylation" in x for x in data.pdbx_role]]
        add_oligos = np.unique(data[["ptnr1_label_comp_id","ptnr2_label_comp_id"]].values.flatten())
        if len(add_oligos)>0: 
            o='\t'.join(add_oligos)
            oligos_file.write(f"{pdb}\t{o}\n")
    else: add_oligos = []
    
    oligos = list(oligos) + list(add_oligos)
    
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
                    pockets = pockets[~np.isin(pockets.ligandfile, filex)]
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
                pockets = pockets[~np.isin(pockets.ligandfile, files2remove)]
                remaining_ligs = [x for x in os.listdir(lig_dir) if pdb in x and chain in x and x!=f and x.endswith(".pdb")]
                # update pockets file
                pockets = pockets[~np.isin(pockets.ligandfile, remaining_ligs)]
                
    
    if has_rebuilt_lig == True:
        # cross ref with PRDs
        prd2pdb = prd2pdb_data.query(f"pdb == '{pdb}'").prd_ID.drop_duplicates().values
        rebuilt_ligs = [x for x in os.listdir(lig_dir) if pdb in x and "lig_chain-" in x]
        if len(prd2pdb)==0:
            print(f"There are {len(rebuilt_ligs)} outstanding ligands with no further annotation in BIRD.")
            for a in rebuilt_ligs: print("\t"+a)
            found_match = None
        else:
            print("Searching for matches in BIRD...")
            # METHOD 1: DIRECT MAPPING
            data = cifdata[pdb.upper()]
            prd_map = copy(rebuilt_ligs)
            if "_pdbx_molecule" in data:
                data = pd.DataFrame.from_dict(data["_pdbx_molecule"], orient="index").T
                for l in rebuilt_ligs:
                    l_chain = l.split("lig_chain-")[-1].split(".")[0]
                    rebuilt_ligs_PRD = data.query(f"asym_id == '{l_chain}'")
                    if len(rebuilt_ligs_PRD)> 0:
                        rebuilt_ligs_PRD = rebuilt_ligs_PRD.prd_id.values[0]
                        prd_map.remove(l)
                        print(f"ligand  {l}  is annotated in BIRD under {rebuilt_ligs_PRD}  ; This is likely a ligand of focus in this PDB entry.")
            
            
            # METHOD 2: INDIRECT MAPPING through sequence
            if len(prd_map)>0: # still chain ligands to be mapped
                print("searching ligands that weren't found in ligand extraction...")
                prd =  CifFileReader().read(cif_file)
                
                found_match = []
                for lig in rebuilt_ligs:
                    lig_structure = PandasPdb().read_pdb(f"{lig_dir}/{lig}").df
                    lig_seq = " ".join(lig_structure["ATOM"].residue_name.drop_duplicates())
                    lig_seq = pd.concat([lig_structure["ATOM"][["residue_number", "residue_name"]], lig_structure["HETATM"][["residue_number", "residue_name"]]]).drop_duplicates().sort_values("residue_number")
                    lig_seq = " ".join(lig_seq.residue_name)
                    
                    for prd_id in prd:
                        if "_pdbx_reference_entity_poly_seq" in prd[prd_id]:
                            obs_seq = pd.DataFrame.from_dict(prd[prd_id]["_pdbx_reference_entity_poly_seq"], orient="index").T
                            full_seq = " ".join(obs_seq.mon_id.values)
                            exp_seq = " ".join(obs_seq.query("observed == 'Y'").mon_id.values)
                            if lig_seq == full_seq or lig_seq == exp_seq:
                                # found a seq-based match
                                found_match.append([lig,exp_seq, prd_id])
                                break 
                
                if len(found_match) == 0: print("None of the rebuilt-chain ligands are in BIRD.")
                else: 
                    for ln in found_match: 
                        print(f"ligand  {ln[0]}  is annotated in BIRD under {ln[1]}  (full sequence match to {ln[2]}). This is likely the ligand of focus in this PDB entry.")
                        
    
    
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
    save_clean_pockets = save_clean_pockets[~np.isin(save_clean_pockets.lig_ID, oligos)]
    
    if len(oligos)>0: print(f"The following oligosacharide residues were found and will be excluded from the ligands: {', '.join(oligos)}")
    
    
    # inspect the size of the pocket
    sparse_pockets = save_clean_pockets.query("pocketres_chain_size<6").ligandfile
    if len(sparse_pockets)>0:
        print("WARNING! There are ligands with very sparse pockets [<6 close res] (which indicates this might not be a ligand). These will be excluded")
        for x in sparse_pockets: 
            print("\t"+x)
            os.remove(f"{lig_dir}/{x}")
        save_clean_pockets = save_clean_pockets[np.isin(save_clean_pockets.ligandfile, sparse_pockets, invert=True)]
    
    save_clean_pockets_list.append(save_clean_pockets)

bar.finish()

oligos_file.close()

# Clean ligands to keep only the final selected list
all_ligands = [x for x in os.listdir(lig_dir) if x.endswith(".pdb")]


# save clean pockets, updated with the selected ligands
save_clean_pockets = pd.concat(save_clean_pockets_list, sort=False)
save_clean_pockets.loc[:,"pdbcode"] = [x.split("_")[0] for x in save_clean_pockets.ligandfile]

save_clean_pockets.to_csv(f"cleanpockets_{prot_dir.split('/')[-1]}.txt", sep="\t", index=False)


title="Clustering Ligands"
padding_total = length - len(title) - 2
print(f"\n{pad_char * (padding_total // 2)} {title} {pad_char * (padding_total - (padding_total // 2))}\n")

subprocess.run(f"rm -Rf {prot_dir}/pdbs_filtered_chains; mkdir {prot_dir}/pdbs_filtered_chains", shell=True)

prot_lst = uniprot2pdbFile.uniprot.unique()

sys.stderr.write("\n")


#for p_i, prot in enumerate(prot_lst):
def clusteringSplit(prot, p_i):
    print(f'**** Protein {prot}')
    pdbs = uniprot2pdbFile.query(f"uniprot == '{prot}'").pdb.str.lower().values
    pockets_prot = save_clean_pockets[np.isin(save_clean_pockets.pdbcode, pdbs)]
    if len(pockets_prot) == 0:
        sys.stderr.write(f"\n\n*** Clustering Protein {prot} ({p_i+1}/{len(prot_lst)} proteins) : bypass as it has no found ligands in any PDBs.\n")
        print(f"bypass {prot} as it has no found ligands in any PDBs.")
        return(prot) #continue
    # reset folder in case this was already run before for this protein
    subprocess.run(f"rm -Rf {prot_dir}/pdbs_filtered_chains/{prot}; mkdir {prot_dir}/pdbs_filtered_chains/{prot}", shell=True)
    #sys.stderr.write(f"\n\n*** Clustering protein {prot} ({p_i+1}/{len(prot_lst)} proteins) ***\n")
    #bar = Bar("Aligning PDBs ... ", max=pockets_prot.pdbcode.nunique()) 
    
    # prior to alignment chains need to be filtered and split when appropriate
    print("CLEANING STRUCTURES BEFORE ALIGNMENT:")
    for pdb in pockets_prot.pdbcode.unique():
        #bar.next()
        pdb_ligs = pockets_prot.query(f"pdbcode == '{pdb}'")
        #ligs_inmultichain = [";" in x for x in pdb_ligs.chain_name.values]
        
        # modres
        ciffile = f"cifs/{pdb}.cif"
        cifdata = CifFileReader().read(ciffile)
        data = cifdata[pdb.upper()]
        if "_pdbx_struct_mod_residue" in data:
            data = pd.DataFrame.from_dict(data["_pdbx_struct_mod_residue"], orient="index").T
            modres = data.label_comp_id.unique()
        else:
            modres = []

        # Split chains to alow self alignment; only keep multiple chains if at least one ligand has then simultaneouly
        allcontactchains = pdb_ligs.chain_name.unique()
        for chainset in allcontactchains:
            chainset = chainset.split(";")
            #ligchains = [x.split("chain-")[1].split(".")[0] for x in pdb_ligs.ligandfile.values if "lig_chain" in x]
            # get ligands that have all its chains in the current chainset (i.e. a ligand with chain K will be included in chainset K;L)
            complexingchains = [np.isin(x.split(";"),chainset).all() for x in pdb_ligs.chain_name]
            ligchains = pdb_ligs[complexingchains].ligandfile.values
            chainfromlig_lst = []
            ligs2keep = []
            for ligFile in ligchains:
                protein_pdb = PandasPdb().read_pdb(f"{lig_dir}/{ligFile}").df
                chainfromlig = np.unique(np.hstack([protein_pdb["ATOM"].chain_id.unique(), protein_pdb["HETATM"].chain_id.unique()]))
                chainfromlig_lst.append(chainfromlig)
                ligs2keep.append(protein_pdb["HETATM"].residue_name.unique())
            chainfromlig_lst = np.unique(np.hstack(chainfromlig_lst))
            allcontactchains = np.unique(np.hstack([chainset,chainfromlig_lst]))
            # save PDB with all participating chains (ligand or protein(s) )
            protein_pdb = PandasPdb().read_pdb(f"{prot_dir}/{pdb}.pdb")
            protein_pdb.df["ATOM"] = protein_pdb.df["ATOM"][protein_pdb.df["ATOM"].chain_id.isin(allcontactchains)]
            protein_pdb.df["HETATM"] = protein_pdb.df["HETATM"][protein_pdb.df["HETATM"].chain_id.isin(allcontactchains)]
            # remove all lines that are HETATM and not in pdb_ligs and NOT is modres
            all_ligs = protein_pdb.df["HETATM"].residue_name.unique()
            ligs2keep = np.unique(np.hstack(ligs2keep))
            ligs2keep = np.setdiff1d(ligs2keep, modres)
            lig2discard = protein_pdb.df["HETATM"][~protein_pdb.df["HETATM"].residue_name.isin(ligs2keep)].residue_name.unique()
            protein_pdb.df["HETATM"] = protein_pdb.df["HETATM"][protein_pdb.df["HETATM"].residue_name.isin(ligs2keep)]
            print(f"removing resname from {pdb}: {';'.join(lig2discard)}")
            # save cleaned file
            protein_pdb.to_pdb(path=f"{prot_dir}/pdbs_filtered_chains/{prot}/{pdb}_keychain{'-'.join(chainset)}.pdb", records=None, gz=False, append_newline=True)
    
    print("\nALIGNMENT --------------------------")
    rmsd=subprocess.run(f'pymol -cq {HOME}/LigExtract/bin/align_pdbs_pockets.py -- {prot_dir}/pdbs_filtered_chains/{prot}', shell=True, capture_output=True)
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
    #sys.stderr.write("\nBuilding ligand clusters...\n")
    unique_pdbs_2alignligs = np.unique([x.split("/")[-1].split("_")[0] for x in glob(f"{prot_dir}/pdbs_filtered_chains/{prot}/aligned_pdbs/*_align.pdb")])
    ligs2cluster = pockets_prot[pockets_prot.pdbcode.isin(unique_pdbs_2alignligs)].ligandfile.values
    for ligfile in ligs2cluster:
        print(ligfile, "-----------------")
        pdb = ligfile.split("_")[0]
        # ref protein with ligands (aligned) - get the corresponding protein for the currently selected ligand
        pdbsref = glob(f"{prot_dir}/pdbs_filtered_chains/{prot}/aligned_pdbs/{pdb}*align.pdb")
        #if pdb == "8uzz": break
        for pdbalign in pdbsref:
            fullPdb = PandasPdb().read_pdb(pdbalign)
            
            # original ligand (not aligned)
            lig_pdb_file = PandasPdb().read_pdb(f"{lig_dir}/{ligfile}")
            
            # get atoms to extract from the original ligand
            lig_lines2keep_het = lig_pdb_file.df["HETATM"].sort_values(by='occupancy', ascending=False).drop_duplicates(subset=["atom_name", "residue_name", "chain_id", "residue_number"]).sort_values("atom_number") # keep first alt_loc
            lig_lines2keep_atm = lig_pdb_file.df["ATOM"].sort_values(by='occupancy', ascending=False).drop_duplicates(subset=["atom_name", "residue_name", "chain_id", "residue_number"]).sort_values("atom_number")
            lig_lines2keep_het = lig_lines2keep_het[["atom_name", "alt_loc", "residue_name", "chain_id", "residue_number"]].values
            lig_lines2keep_atm = lig_lines2keep_atm[["atom_name", "alt_loc", "residue_name", "chain_id", "residue_number"]].values

            #if pdb == "8uv1": print(asfdhasfdh)
            #lig_atomnum2keep = np.hstack([lig_atomnum2keep_het, lig_atomnum2keep_atm])
            # save atoms from the aligned full protein file
            filtered_lig_het = [fullPdb.df["HETATM"].query(f'atom_name=="{ln[0]}" and alt_loc=="{ln[1]}" and residue_name=="{ln[2]}" and chain_id=="{ln[3]}" and residue_number=={ln[4]}') for ln in lig_lines2keep_het]
            filtered_lig_atm = [fullPdb.df["ATOM"].query(f'atom_name=="{ln[0]}" and alt_loc=="{ln[1]}" and residue_name=="{ln[2]}" and chain_id=="{ln[3]}" and residue_number=={ln[4]}') for ln in lig_lines2keep_atm]
    
            if len(filtered_lig_het)==0:
                fullPdb.df["HETATM"] = fullPdb.df["HETATM"].query("atom_name=='NOTHING'")
            else:
                filtered_lig_het = pd.concat(filtered_lig_het)
                if len(lig_lines2keep_het) != len(filtered_lig_het): continue # it is not getting the lines from correct reference pdb
                fullPdb.df["HETATM"] = filtered_lig_het
            
            if len(filtered_lig_atm)==0:
                fullPdb.df["ATOM"] = fullPdb.df["ATOM"].query("atom_name=='NOTHING'")
            else:
                filtered_lig_atm = pd.concat(filtered_lig_atm)
                if len(lig_lines2keep_atm) != len(filtered_lig_atm): continue # it is not getting the lines from correct reference pdb
                fullPdb.df["ATOM"] = filtered_lig_atm
            #liglines = save_aligned_lig.df["HETATM"].shape[0]+save_aligned_lig.df["ATOM"].shape[0]
            newligname = ligfile.replace('.pdb','_aligned_LIG.pdb')
            fullPdb.to_pdb(path=f"{prot_dir}/pdbs_filtered_chains/{prot}/aligned_pdbs/{newligname}", records=['ATOM', 'HETATM'], gz=False, append_newline=True)
            centroid = pd.concat([fullPdb.df["HETATM"], fullPdb.df["ATOM"]])[["x_coord", "y_coord","z_coord"]].mean().values.round(3)
            print(newligname, centroid)
            ligand_centroid.append(np.hstack([newligname, centroid]))
    
    # Cluster controids
    ligand_centroid = pd.DataFrame(ligand_centroid)
    if ligand_centroid[0].nunique() != len(ligand_centroid):
        print("WARNING!!! Some ambiguities detected in the ligands that will be used in the clustering. Please inspect inspect_ligand_centroids.txt")
        ligand_centroid.to_csv("inspect_ligand_centroids.txt", sep="\t", index=False)
    centroids_coords = ligand_centroid[[1,2,3]].astype(float)#.round(3)
    if len(centroids_coords) > 1:
        clustering = AgglomerativeClustering(n_clusters=None, distance_threshold=10, linkage="single").fit(centroids_coords)
        cluster_classes = clustering.labels_
        clust_ec = {cl:max([round(euclidean(i,j),3) for i,j in combinations(centroids_coords[cluster_classes==cl].values,2)]) if sum(cluster_classes==cl)>1 else 0 for cl in np.unique(cluster_classes)}
    else:
        cluster_classes = np.array([0])
        clust_ec = {0:0}
    #sys.stderr.write(f"--- {len(np.unique(cluster_classes))} different clusters were found.\n\n")
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
    #for c in ligand_centroid_sorted.cluster.unique():
    #    pocket_ligs = ligand_centroid_sorted.query(f"cluster == '{c}'").ligandfile.values
    #    pocket_ligs = " ".join(pocket_ligs)
    #    subprocess.run(f"pymol -cq {HOME}/LigExtract/bin/align_ligs_figures.py -- pocketcluster{c} {prot_dir}/pdbs_filtered_chains/{prot}/aligned_pdbs {pocket_ligs}", shell=True, capture_output=True)
    
    #sys.stderr.write(f"Clustered pockets have been stored in {prot}_pockets_hierarch-clusters.txt. Pictures of the different pocket clusters created in {prot_dir}/pdbs_filtered_chains/{prot}/aligned_pdbs")
    
    # Produce one global PyMOL session where each cluster has a color
    subprocess.run(f"pymol -cq {HOME}/LigExtract/bin/allpockets_figure.py -- {prot_dir}/pdbs_filtered_chains/{prot}/aligned_pdbs {prot}_pockets_hierarch-clusters.txt", shell=True, capture_output=True)
    return(prot)


# multithreading
num_workers = os.cpu_count()-1

with ProcessPoolExecutor(max_workers=num_workers) as executor: #num_workers
    #executor.map(extractorSplit, pdbs)
    futures = [executor.submit(clusteringSplit, protein, pi) for pi,protein in enumerate(prot_lst)]
    barthread = Bar('Aligning chains and Clustering ligands... ', max=len(futures))
    results = []
    for f in as_completed(futures):
        results.append(f.result())
        barthread.next()
    barthread.finish()

sys.stderr.write(f"\n\nFinished Clustering pockets.\n\n")
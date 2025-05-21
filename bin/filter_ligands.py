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
import gc


parser = argparse.ArgumentParser(description='Final procedure do decide which ligands are the most likely in a given protein, using different criteria that looks into ligand structure and annotation, ligand placement in the pocket, etc.')
parser.add_argument('--pdbPath', type=str, required=True, dest="prot_dir", help='Path to directory containing all PDBs to process.')
parser.add_argument('--ligandsPath', type=str, required=True, dest="lig_dir", help='Path to directory with extracted ligand files.')
parser.add_argument('--prdCif', type=str, required=True, dest="cif_file", help='CIF File with all ligands in BIRD')
args = parser.parse_args()

prot_dir = args.prot_dir
lig_dir = args.lig_dir
cif_file = args.cif_file

if prot_dir[-1]=="/": prot_dir = prot_dir[:-1]
if lig_dir[-1]=="/": lig_dir = lig_dir[:-1]

length = 90; pad_char = '-'

title="Final Ligand Selection"
padding_total = length - len(title) - 2
print(f"\n{pad_char * (padding_total // 2)} {title} {pad_char * (padding_total - (padding_total // 2))}\n")


pockets = pd.read_csv(f"pockets_{prot_dir.split('/')[-1]}.txt", sep="\t")
keep_ligands_file = open(f"{lig_dir}_final_ligands2keep.txt","w")
oligos_file = open(f"{lig_dir}_oligosaccharides.txt","w")
oligos_file.write("WARNING! Some protein residues might be added here due to the nature of the N-glycosylation group detection procedure.\n")
oligos_file.write("         However this has no impact on the whole procedure as this list is used to exclude ligands for being N-glycosylation\n")
oligos_file.write("         groups. Even if a protein residue is within the list of ligands (i.e. a modified residue in HETATM),\n")
oligos_file.write("         it would simply remove that supposed ligand (correctly, as it is actually not a ligand)\n\n\n")

pdbs_in_pockets = pockets.pdbcode.drop_duplicates()
#pdbs_in_pockets = pdbs_in_pockets[np.argwhere(pdbs_in_pockets=="1bb0").flatten()[0]:]

save_clean_pockets_list = []

bar = Bar('Filtering ligands... ', max=len(pdbs_in_pockets))

for pdb in pdbs_in_pockets:#["3ilr"]:#
    bar.next()
    print(f"\n############# {pdb} ###########\n")

    if len(glob(f'{lig_dir}/{pdb}_*.pdb'))==0:
        print("There are no ligands to consider.")
        continue
    # Remove ligands that are glycosylation groups
    #sleep(2)
    #subprocess.run(f"wget https://files.rcsb.org/view/{pdb}.cif --quiet -O {pdb}.cif", shell=True)
    cif_f = open(f"cifs/{pdb}.cif").read()
    prd_chain_pdb = pd.read_csv("prdID_chain_pdb.txt", sep="\t")
    if "_pdbx_entity_branch_list" in cif_f: 
        cmpd_where = np.argwhere([".comp_id" in x for x in cif_f.split("\n") if x.startswith("_pdbx_entity_branch_list.")]).flatten()[0]
        oligos = cif_f.split("_pdbx_entity_branch_list.hetero \n")[1].split("#")[0]
        oligos = [x for x in oligos.split("\n")[:-1]]
        oligos = np.unique([shlex.split(x)[cmpd_where] for x in oligos])
        o='\t'.join(oligos)
        oligos_file.write(f"{pdb}\t{o}\n")
    else: oligos = []

    if "_struct_conn.pdbx_role" in cif_f:
        headers_cif = [x for x in cif_f.split("\n") if x.startswith("_struct_conn.")]
        where_resnames = np.argwhere(["ptnr1_label_comp_id" in x or "ptnr2_label_comp_id" in x for x in headers_cif]).flatten()
        add_oligos = cif_f.split("_struct_conn.pdbx_role \n")[1].split("#")[0]
        add_oligos = [x for x in add_oligos.split("\n")[:-1]]
        join_ln = np.argwhere([len(shlex.split(x))<len(headers_cif) for x in add_oligos]).flatten()
        for x in join_ln:
            # if current line is empty, this is a continuation line: bypass.
            if len(add_oligos[x].split())==0: continue
            add_oligos[x] = add_oligos[x]+add_oligos[x+1]
            add_oligos[x+1]=""
        add_oligos = [shlex.split(x) for x in add_oligos if "Glycosylation" in x]
        if len(add_oligos)>0: add_oligos = np.hstack(np.array(add_oligos)[:,where_resnames]) 
        # ! some of these will be residues but this will have no impact on the process
        # since residues will never be a possible individual ligand to begin with.
        o='\t'.join(add_oligos)
        oligos_file.write(f"{pdb}\t{o}\n")
    else: add_oligos = []
    
    oligos = list(oligos) + list(add_oligos)

    del cif_f # no need to have this huge file in memory
    gc.collect()
    #os.remove(f"{pdb}.cif")
    # get all identified ligs
    logfile = [ln.strip() for ln in open(f"{lig_dir}/{pdb}_ligand_extraction.log").readlines()][5:]
    chain2prdIDs = np.array([[x.split("rebuilt ligand detected in chain ")[1][0],x.split("PRD ID ")[1].split(" ")[0]] for x in logfile if "associated with PRD ID" in x])
    has_rebuilt_lig = any(["lig_chain-" in x for x in os.listdir(lig_dir) if pdb in x])
    rare_ligand = [x.split()[0] for x in logfile if "ligand " in x and "Might not be a ligand." not in x and "rebuilt" not in x]

    # update rare ligands from the log file
    currentfiles = [x for x in os.listdir(lig_dir) if pdb in x and x.endswith("pdb")]
    currentfiles = [x.split("_lig-")[-1].split(".")[0] for x in currentfiles]
    if np.intersect1d(currentfiles,oligos).shape[0]>0: print("The current ligands are glycosylation groups, and will be excluded:", " ".join(np.intersect1d(currentfiles,oligos)))
    currentfiles = np.setdiff1d(currentfiles, oligos)
    rare_ligand = np.intersect1d(rare_ligand,currentfiles)

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
                files2remove = np.hstack([[x for x in os.listdir(lig_dir) if l in x and chain in x] for l in rareLigs_in_rebuiltLigs])
                print("remove", "; ".join(files2remove), f"--> contained in {f}")
                files2save.append(f)
                for x in files2remove: os.remove(f'{lig_dir}/{x}')
                # update pockets file
                pockets = pockets[~np.isin(pockets.ligandfile, files2remove)]
                remaining_ligs = [x for x in os.listdir(lig_dir) if pdb in x and chain in x and x!=f and x.endswith(".pdb")]
                if len(remaining_ligs)>0: print("remove remaining ligands that are likely to be crystallography elements and N-glycolysation branches:")
                for a in remaining_ligs: print("\t",a)
                for x in remaining_ligs: os.remove(f'{lig_dir}/{x}')
                # update pockets file
                pockets = pockets[~np.isin(pockets.ligandfile, remaining_ligs)]

            if len(rareLigs_in_rebuiltLigs) == 0:
                # both the small and the larger molecule can be a ligand
                decision="small_mol"
                remaining_ligs = np.array([x for x in os.listdir(lig_dir) if pdb in x and chain in x and x!=f and x.endswith(".pdb")])
                rare_lig_file = np.hstack([[x for x in os.listdir(lig_dir) if pdb in x and chain in x and x!=f and r in x] for r in rare_ligand])
                remaining_ligs = np.setdiff1d(remaining_ligs, rare_lig_file)
                print(f"Consider as likely ligands in chain {chain[-1]}: {', '.join(np.hstack([f,rare_lig_file]))}")
                files2save+=[f]+list(rare_lig_file)
                if ligtype == "large": print(f"However {f} is larger than 10 aminoacids")
                if len(remaining_ligs)>0: print("remove ligands that are likely to be crystallography elements and N-glycolysation branches:")
                for a in remaining_ligs: print("\t",a)
                for x in remaining_ligs: os.remove(f'{lig_dir}/{x}')
                # update pockets file
                pockets = pockets[~np.isin(pockets.ligandfile, remaining_ligs)]
            
        # catch ligands left out, after loop finishes; only remove ligands from current chain; otherwise all other ligands in other chains will be removed
        clean_leftout = [x for x in os.listdir(lig_dir) if pdb in x and x.endswith(".pdb") and x not in files2save and f'_chain-{chain[-1]}_' in x]
        for x in clean_leftout: os.remove(f'{lig_dir}/{x}')
        # update pockets file
        pockets = pockets[~np.isin(pockets.ligandfile, clean_leftout)]
        if len(clean_leftout)>0: print("remove ligands that are likely to be crystallography additives and N-glycolysation branches:")
        for a in clean_leftout: print("\t",a)
    
    if has_rebuilt_lig == False and len(rare_ligand)>0: # consider the rare ligand
        remaining_ligs = np.unique([x for x in os.listdir(lig_dir) if pdb in x and x.endswith(".pdb")])
        rare_lig_file = np.hstack([[x for x in os.listdir(lig_dir) if pdb in x and r in x] for r in rare_ligand])
        remaining_ligs = np.setdiff1d(remaining_ligs, rare_lig_file)
        print(f"Consider as the most likely ligand(s): {' '.join(rare_lig_file)}")
        if len(rare_lig_file) == 1: keep_ligands_file.write(rare_lig_file[0]+"\n")
        if len(remaining_ligs)>0: print("remove ligands that are likely to be crystallography additives and N-glycolysation branches:")
        for a in remaining_ligs: print("\t",a)
        for x in remaining_ligs: os.remove(f'{lig_dir}/{x}')
        # update pockets file
        pockets = pockets[~np.isin(pockets.ligandfile, remaining_ligs)]

    if has_rebuilt_lig == False and len(rare_ligand)==0: # no ligands of any kind
        print("There are no chain ligands and all small molecules are likely crystallographic additives/reagents as they occur in hundreds other proteins")
        print("All current ligands will be therefore removed:")
        remaining_ligs = np.unique([x for x in os.listdir(lig_dir) if pdb in x and x.endswith(".pdb")])
        for a in remaining_ligs: 
            print("\t",a)
            os.remove(f'{lig_dir}/{a}')
        continue
        pockets = pockets[~np.isin(pockets.ligandfile, remaining_ligs)]

    if has_rebuilt_lig == True:
        # clean non-rare ligs
        if len(rare_ligand)>0: 
            rare_lig_file = np.hstack([[x for x in os.listdir(lig_dir) if pdb in x and r in x] for r in rare_ligand])
        else:
            rare_lig_file = []
        remaining_ligs = np.unique([x for x in os.listdir(lig_dir) if pdb in x and x.endswith(".pdb") and "_lig-" in x])
        remaining_ligs = np.setdiff1d(remaining_ligs, rare_lig_file)
        if len(remaining_ligs)>0:
            print("Clean all ligands that are likely to be crystallographic additives/reagents:")
            pockets = pockets[~np.isin(pockets.ligandfile, remaining_ligs)]
        for a in remaining_ligs: 
            print("\t",a)
            os.remove(f'{lig_dir}/{a}')
            
        # cross ref with PRDs
        #sleep(3)
        #subprocess.run(f"wget https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/{pdb} --quiet -O {pdb}.dict", shell=True)
        prd2pdb = pd.read_csv("prd_to_pdb_IDs.txt", sep="\t")
        prd2pdb = prd2pdb.query(f"pdb == '{pdb}'").prd_ID.drop_duplicates().values
        rebuilt_ligs = [x for x in os.listdir(lig_dir) if pdb in x and "lig_chain-" in x]
        if len(prd2pdb)==0:
            print(f"There are {len(rebuilt_ligs)} outstanding ligands with no basic annotation in BIRD.")
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
                keep_ligands_file.write(l+"\n")

            # add chain ligands from prdID_chain_pdb.txt
            for x in rebuilt_ligs: 
                ch = x.split("chain-")[1].split(".")[0] 
                match_prd2chain = prd_chain_pdb.query(f"pdb=='{pdb}' and asym_id=='{ch}'")
                if len(match_prd2chain)>0:
                    prd_match = match_prd2chain.prd_id.values[0]
                    print(f"ligand  {x}  is annotated in BIRD under {prd_match}  (chain-to-prdID match). This is likely the ligand of focus in this PDB entry.")
                    keep_ligands_file.write(x+"\n")


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
                get_blocks = prd_block[np.isin(titles, prd2pdb)]

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

                #deprecated SEQRES
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
                        keep_ligands_file.write(ln[0]+"\n")
    
    
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
    save_clean_pockets.loc[:,"ligtype"] = np.where(["lig_chain" in x for x in save_clean_pockets.ligandfile], "chain ligand", "small-molecule ligand")
    save_clean_pockets.loc[:,"lig_ID"] = [x.split("lig-")[-1].split(".")[0] for x in save_clean_pockets.ligandfile]
    save_clean_pockets = save_clean_pockets[~np.isin(save_clean_pockets.lig_ID, oligos)]
    save_clean_pockets.loc[:,"pdbcode"] = [x[0:4] for x in save_clean_pockets.ligandfile]
    #save_clean_pockets_list.append(save_clean_pockets)
    if len(oligos)>0: print(f"The following oligosacharide residues were found and will be excluded from the ligands: {', '.join(oligos)}")

    # Untie chain among the surviving ligands
    save_clean_pockets.index=range(len(save_clean_pockets))
    # infer duplicates from same PRD_ID
    for prd in prd_chain_pdb.query(f"pdb=='{pdb}'").prd_id.unique():
        prd_duplicates = prd_chain_pdb.query(f"pdb=='{pdb}' and prd_id=='{prd}'")
        if len(prd_duplicates)<2: continue
        prd_duplicates = prd_duplicates.asym_id.values
        rep_lig = [x for x in save_clean_pockets.query(f"pdbcode=='{pdb}'").ligandfile.values if "lig_chain-" in x and x.split("_lig_chain-")[-1][0] in prd_duplicates]
        if len(rep_lig)==0: continue
        sorted_reps = save_clean_pockets[np.isin(save_clean_pockets.ligandfile, rep_lig)].sort_values("pocketres_chain_size", ascending=False).ligandfile.values
        # first item in sorted_reps is the one to keep; the others are repeated
        save_clean_pockets = save_clean_pockets[~np.isin(save_clean_pockets.ligandfile, sorted_reps[1:])]
        print(f"deduplicate from {rep_lig} ---> keep {sorted_reps[0]}")
    unique_ligs_idx = save_clean_pockets.sort_values("pocketres_chain_size", ascending=False).drop_duplicates(subset=["lig_ID"],keep="first").index
    removed_idx = np.setdiff1d(save_clean_pockets.index,unique_ligs_idx)
    print("deduplication by lig_ID removed ligands:")
    for x in save_clean_pockets.loc[removed_idx].ligandfile.values:
        print(f"\t{x}")

    save_clean_pockets = save_clean_pockets.loc[unique_ligs_idx]
    print("deduplication by lig_ID kept ligands:")
    for x in save_clean_pockets.ligandfile.values:
        print(f"\t{x}")

    # inspect the size of the pocket
    sparse_pockets = save_clean_pockets.query("pocketres_chain_size<6").ligandfile
    if len(sparse_pockets)>0:
        print("WARNING! There are ligands with very sparse pockets [<6 close res] (which indicates this might not be a ligand). These will be excluded")
        for x in sparse_pockets: 
            print("\t"+x)
            os.remove(f"{lig_dir}/{x}")

    # Scenario 1 : all ligands are small molecules
    if len(save_clean_pockets)>1 and (save_clean_pockets.ligtype == "small-molecule ligand").all() and save_clean_pockets.lig_ID.nunique()==1:
        # this can only be done if the two ligands are equivalent (i.e. typical ligands)
        print(f"There are {len(save_clean_pockets)} identified ligands, and all are the same molecule ({save_clean_pockets.lig_ID.unique()[0]}) ")
        print(f"The ligand with the most nearby residues (according to set distmax) will be selected...")
        max_pocket = save_clean_pockets.pocketres_chain_size.max()
        save_clean_pockets = save_clean_pockets.query(f"pocketres_chain_size == '{max_pocket}'")
        if len(save_clean_pockets)>1: print(f"there are {len(save_clean_pockets)} binding sites with max residue count ({max_pocket}). The first one will be kept.")
        print("The ligand kept is:")
        save_clean_pockets = save_clean_pockets[0:1]
        print(save_clean_pockets.ligandfile.values[0])
        keep_ligands_file.write(save_clean_pockets.ligandfile.values[0]+"\n")

    # Scenario 2 : there are ligands of both types - small molecules and chain ligands
    elif len(save_clean_pockets)>1 and save_clean_pockets.ligtype.nunique()==2:
        print(f"There are {len(save_clean_pockets)} identified ligands, and there are molecules of different types ({save_clean_pockets.ligtype.unique()}) ")
        # keep unique small mols 
        unique_smallcpds = [x for x,y in save_clean_pockets[["ligandfile","lig_ID"]].values if sum(save_clean_pockets.lig_ID.values==y)==1 and "lig_chain-" not in y]
        if len(unique_smallcpds)>0: 
            print("keep unique small ligands:")
            for uniq_lig in unique_smallcpds:
                print(f"\t{uniq_lig}")
                keep_ligands_file.write(uniq_lig+"\n")
        dupe_cadidates = np.unique([x for x in save_clean_pockets.lig_ID.values if sum(save_clean_pockets.lig_ID.values==x)>1])
        # deduplicate small-ligs
        for dupe_lig in dupe_cadidates:
            unique_lig = save_clean_pockets.query(f"lig_ID == '{dupe_lig}'").sort_values("pocketres_chain_size", ascending=False)[0:1]
            unique_lig = unique_lig.ligandfile.values[0]
            rep_lig = unique_lig.ligandfile.values[1:]
            save_clean_pockets = save_clean_pockets[~np.isin(save_clean_pockets.ligandfile, rep_lig)]
            print(f"keep deduped ligand {unique_lig}")
            keep_ligands_file.write(unique_lig+"\n")

        # dedupe chain ligs
        chain_pockets = save_clean_pockets[["lig_chain-" in x for x in save_clean_pockets["lig_ID"].values]].copy()
        print(chain_pockets.shape,"################")
        if len(chain_pockets)>1:
            print("De-duplication of chain ligands by inspecting MOL_ID and keeping the ligand of each MOL_ID with the most residues closeby.")
            pdb_mols = [ln.strip() for ln in open(f'{prot_dir}/{pdb}.pdb').readlines() if ln.startswith("COMPND ")]
            chunks = [x for x,_ in enumerate(pdb_mols) if "MOL_ID: " in pdb_mols[x]]
            mol_samechain = []
            for a,b in zip(chunks, chunks[1:]+[len(pdb_mols)]): mol_samechain.append("".join(pdb_mols[a:b]).split("CHAIN: ")[1].split(";")[0])
            mol_samechain = np.vstack([[[k,i+1] for k in ln.split(", ")] for i,ln in enumerate(mol_samechain)] )
            mol_samechain = {k:i for k,i in mol_samechain}
            chain_pockets.loc[:,"same_molid"] = [mol_samechain[ch[-1]] if ch[-1] in mol_samechain else ch for ch in chain_pockets["lig_ID"].values ]
            chain_pockets = chain_pockets.sort_values("pocketres_chain_size", ascending=False).drop_duplicates("same_molid")
            rep_lig = [x for x in save_clean_pockets.ligandfile.values if "lig_chain-" in x and x not in chain_pockets.ligandfile.values]
            save_clean_pockets = save_clean_pockets[~np.isin(save_clean_pockets.ligandfile, rep_lig)]
            # infer duplicates from same PRD_ID
            for prd in prd_chain_pdb.query(f"pdb=='{pdb}'").prd_id.unique():
                prd_duplicates = prd_chain_pdb.query(f"pdb=='{pdb}' and prd_id=='{prd}'")
                if len(prd_duplicates)<2: continue
                prd_duplicates = prd_duplicates.asym_id.values
                rep_lig = [x for x in save_clean_pockets.query(f"pdbcode=='{pdb}'").ligandfile.values if "lig_chain-" in x and x.split("_lig_chain-")[-1][0] in prd_duplicates]
                sorted_reps = save_clean_pockets[np.isin(save_clean_pockets.ligandfile, rep_lig)].sort_values("pocketres_chain_size", ascending=False).ligandfile.values
                # first item in sorted_reps is the one to keep; the others are repeated
                save_clean_pockets = save_clean_pockets[~np.isin(save_clean_pockets.ligandfile, sorted_reps[1:])]
                print(f"deduplicate from {rep_lig} ---> keep {sorted_reps[0]}")


            print("\tRemoved duplicate ligands",rep_lig)

        print("\n     ------ All ligands kept:")
        for x in save_clean_pockets.ligandfile: 
            print("\t"+x)
            keep_ligands_file.write(x+"\n")

    # Scenario 3 : all ligands are chain ligands
    elif len(save_clean_pockets)>1 and (save_clean_pockets.ligtype == "chain ligand").all():
        print(f"There are {len(save_clean_pockets)} identified ligands, and all molecules are of type ({save_clean_pockets.ligtype.unique()[0]}) ")
        # check if they are the same molecule
        pdb_mols = [ln.strip() for ln in open(f'{prot_dir}/{pdb}.pdb').readlines() if ln.startswith("COMPND ")]
        chunks = [x for x,_ in enumerate(pdb_mols) if "MOL_ID: " in pdb_mols[x]]
        mol_samechain = []
        for a,b in zip(chunks, chunks[1:]+[len(pdb_mols)]): mol_samechain.append("".join(pdb_mols[a:b]).split("CHAIN: ")[1].split(";")[0])
        mol_samechain = np.vstack([[[k,i+1] for k in ln.split(", ")] for i,ln in enumerate(mol_samechain)] )
        mol_samechain = {k:i for k,i in mol_samechain}
        save_clean_pockets.loc[:,"same_molid"] = [mol_samechain[ch[-1]] if ch[-1] in mol_samechain else ch for ch in save_clean_pockets["lig_ID"].values ]
        
        # dedupe molecules of the same molecule ID, selecting the largest pocket
        save_clean_pockets = save_clean_pockets.sort_values("pocketres_chain_size",ascending=False).drop_duplicates("same_molid")

        print("keep ligands, after deduplication according to MOL_ID and selecting the largest close residue count (for tied counts, the first in the list is selected):")
        for x in save_clean_pockets.ligandfile: 
            print("\t"+x)
            keep_ligands_file.write(x+"\n")

    else:
        print("keep all remaining unique ligands that survived the filtration process")
        for x in save_clean_pockets.ligandfile: 
            print("\t"+x)
            keep_ligands_file.write(x+"\n")

    save_clean_pockets_list.append(save_clean_pockets)

keep_ligands_file.close()
oligos_file.close()

# Clean ligands to keep only the final selected list
keep_ligands = [ln.strip() for ln in open(f"{lig_dir}_final_ligands2keep.txt").readlines()]
all_ligands = [x for x in os.listdir(lig_dir) if x.endswith(".pdb")]
for f in np.setdiff1d(all_ligands,keep_ligands): os.remove(f"{lig_dir}/{f}")

# save clean pockets, updated with the selected ligands
save_clean_pockets_list = pd.concat(save_clean_pockets_list, sort=False)

k = np.setdiff1d(keep_ligands, save_clean_pockets_list.ligandfile)
if k.shape[0]>0: print(f"\nINSPECT LIGANDS! Some ligands in {lig_dir}_final_ligands2keep.txt are meant to be kept but are no longer in the pockets list: {','.join(k)}")
save_clean_pockets_list = save_clean_pockets_list[np.isin(save_clean_pockets_list.ligandfile, keep_ligands)]
save_clean_pockets_list.to_csv(f"cleanpockets_{prot_dir.split('/')[-1]}.txt", sep="\t", index=False)


sys.stderr.write("\n\nFinished.")
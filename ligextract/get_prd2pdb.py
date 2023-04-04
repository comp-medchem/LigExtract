# python ~/scripts/get_prd2pdb.py prd_cif_file
# e.g.
# python ~/scripts/get_prd2pdb.py prd-all.cif


import pandas as pd
import numpy as np
import sys
import linecache

cif_file = sys.argv[1] # prd-all.cif

prd2pdb = []
prd_block =[]
prd2lig=[]
with open(cif_file) as prd:
    for line in prd:
        if line.startswith("data_PRD_"):
            #print(prd_block)
            prd_id = [x.strip().split()[-1] for x in prd_block if x.startswith("_pdbx_reference_molecule.prd_id")]
            pdb_code = [x.strip().split()[-1].lower() for x in prd_block if "representative_PDB_id_code" in x]
            lig_code = [x.strip().split()[-1] for x in prd_block if "_pdbx_reference_molecule.chem_comp_id" in x]
            prd2pdb.append(list(np.hstack([prd_id,pdb_code])))
            prd2lig.append(list(np.hstack([prd_id,lig_code])))
            prd_block = [] # start new block
        prd_block.append(line)

prd2lig = prd2lig[1:]
prd2pdb = prd2pdb[1:]
prd2pdb = pd.DataFrame(prd2pdb, columns=["prd_ID","pdb"])
prd2lig = pd.DataFrame(prd2lig, columns=["prd_ID","pdblig_ID"])
prd2lig = prd2lig.query("pdblig_ID != '?'")
prd2pdb.to_csv("prd_to_pdb_IDs.txt", sep="\t", index=False)
prd2lig.to_csv("prd_to_pdb_ligIDs.txt", sep="\t", index=False)

print("extracted PRD_IDs to PDB_IDS --> stored into prd_to_pdb_IDs.txt")
print("extracted PRD_IDs to ligand IDS --> stored into prd_to_pdb_ligIDs.txt")


# Get all PRD instances and their sequences
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
    if '_pdbx_reference_entity_poly_seq' not in "".join(block):
    	continue
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


seq_list = pd.DataFrame(seq_list)
seq_list.to_csv("../data/prd_cmpds_seq_full_experim.txt", sep="\t", index=False, header=None)
# python ~/scripts/process_chains.py --pdbspath targetsofinterest --pdbsToRedo None

import sys,os
import numpy as np
import urllib
from urllib import parse,request
from time import sleep
from progress.bar import ChargingBar, Bar
import argparse
import pandas as pd
import gzip

parser = argparse.ArgumentParser(description='Process all PDBs to annotate all chains with a possible protein identifier.')
parser.add_argument('--pdbspath', type=str, required=True, dest="pdbpath", help='Path to directory containing all PDBs to process')
parser.add_argument('--pdbsToRedo', type=str, required=False, nargs='+', default=["None"], dest="pdbs_redo", help='a list of space-separated pdb files to re-run (only used in re-runs of this script)')
args = parser.parse_args() #                                 nargs='+'

pdbpath = args.pdbpath
pdbs_redo = args.pdbs_redo
if pdbpath[-1]=="/": pdbpath = pdbpath[:-1]


print("\n------------------  Processing PDBs to gather chain information  ------------------\n")

pdbs_lst = [x for x in os.listdir(pdbpath) if x.endswith(".pdb")]

catch_html_err = open(f"{pdbpath.split('/')[-1]}_process_uniprot_chains.err","a")

save_uniprot = open(f"{pdbpath.split('/')[-1]}_process_uniprot_chains.txt","a")

if "None" in pdbs_redo:
    save_uniprot.write("pdb\tuniprot_in_pdbfile\tupdated_uniprot\n")
elif len(pdbs_redo)>0 and pdbs_redo[0]!="None":
    pdbs_lst = pdbs_redo[:] #prevent aliasing
else:
    sys.exit("uncaught scenario. This appears to be an empty list provided to --pdbsToRedo")


# if mapping from old to new uniprot IDs
bar = Bar('Collecting all uniprot IDs in PDB chains...', max=len(pdbs_lst))

testprotein_lst = []
for pdbname in pdbs_lst: 
    bar.next()
    pdb = open(f"{pdbpath}/{pdbname}").readlines()
    if len(pdb) == 0:
        catch_html_err.write(f"{pdbname}: PDB has zero lines. Please verify the file is valid.\n")
        continue
    chaintype = [x for x in pdb if x.startswith('DBREF ')]
    testprotein = np.unique([x[33:41].strip() for x in chaintype if x[26:32].strip()=="UNP"])
    testprotein_lst.append(testprotein)

sys.stderr.write("\n")

testprotein_lst = np.unique(np.hstack(testprotein_lst))
chunks_idx = list(range(0,len(testprotein_lst),5000))
chunks_idx = list(zip(chunks_idx,chunks_idx[1:]+[len(testprotein_lst)]))
print(f'{len(testprotein_lst)} proteins found in PDB chains. Uniprot mapping will be done in {len(chunks_idx)} chunk(s)...')

old2new=[]
for chunk_i,chunk_f in chunks_idx:
    retries = 0
    success = False
    chunk = testprotein_lst[chunk_i:chunk_f]
    params = {'from': 'ACC+ID', 'to': 'ACC', 'format': 'tab', 'query': f'{" ".join(chunk)}'}
    data = parse.urlencode(params).encode('utf-8')
    while success == False and retries < 3:
        try:
            sleep(6)
            req = request.Request('https://www.uniprot.org/uploadlists/', data)
            r = request.urlopen(req) # Doing this in case old Uniprot IDs are being used
            success = True
        except urllib.error.HTTPError as err:
            # retry
            sleep(10)
            retries+=1
    if success == False:
        catch_html_err.write(f"{pdbname}: Request to Uniprot failed after {retries} retries! PDB {pdbname} will be bypassed.\n")
        continue
    try:
        r = r.read().decode("utf-8").split("\n")[1:-1]
    except ConnectionResetError:
        catch_html_err.write(f"{pdbname}: Request to Uniprot failed after {retries} retries! PDB {pdbname} will be bypassed.\n")
        continue
    if len(r)>0:
        r = [x.strip().split("\t") for x in r]
        for q,retrieved_unip in r:
            old2new.append([q,retrieved_unip])


unip_dict = {} # old:new
old2new = pd.DataFrame(old2new, columns=["old","new"]).drop_duplicates()
for old in old2new["old"].values:
    new = ",".join(old2new.query(f"old=='{old}'").new.values)
    unip_dict[old]=new


bar = Bar('Processing...', max=len(pdbs_lst))

for pdbname in pdbs_lst: 
    bar.next()
    pdb = open(f"{pdbpath}/{pdbname}").readlines()
    if len(pdb) == 0:
        catch_html_err.write(f"{pdbname}: PDB has zero lines. Please verify the file is valid.\n")
        continue

    pdbcode=pdbname.split(".")[0]
    chaintype = [x for x in pdb if x.startswith('DBREF ')]
    testprotein_chain = pd.DataFrame([[x[12].strip(),x[33:41].strip()] for x in chaintype if x[26:32].strip()=="UNP"])
    if len(testprotein_chain)>0:
        testprotein_chain = testprotein_chain.groupby([0])[1].agg(lambda x: ','.join(x))

    chain2uniprot = {}
    with gzip.open(f"../data/pdb_chain_uniprot.csv.gz") as f:
        for ln in f: # PDB  CHAIN   SP_PRIMARY
            ln=ln.decode("utf-8").strip().split(",")
            if ln[0]==pdbcode.lower():
                if ln[1] in chain2uniprot.keys():
                    chain2uniprot = {}
                    # return an empty dict because of N-to-N mapping of chains to uniprot (it is already found)
                    break
                chain2uniprot[ln[1]]=ln[2]

    premapping_old_new = []
    if len(testprotein_chain) > 0:
        for c in testprotein_chain.reset_index()[0].unique():
            uniprots = testprotein_chain.loc[c].split(",")
            if len(uniprots)>1:
                # ambiguous N-to-N mapping
                continue
            if c in chain2uniprot.keys():
                premapping_old_new.append([uniprots[0], chain2uniprot[c]])

    allchains = pd.DataFrame([x[12].strip() for x in chaintype])
    for c in np.unique(allchains):
        if c in chain2uniprot.keys():
            premapping_old_new.append([chain2uniprot[c], chain2uniprot[c]])

    if len(premapping_old_new)>0:
        premapping_old_new = list(pd.DataFrame(premapping_old_new).drop_duplicates().values)
        premapping_old_new = ["\t".join(x) for x in premapping_old_new]
    
    if len(premapping_old_new)>0:
        for a in premapping_old_new: save_uniprot.write(pdbcode+"\t"+a+"\n")
    
    testprotein = np.unique([x[33:41].strip() for x in chaintype if x[26:32].strip()=="UNP"])

    for q in testprotein:
        if q in unip_dict.keys():
            new_unip=unip_dict[q]
        else:
            continue
        save_uniprot.write(f"{pdbname.split('.')[0]}\t{q}\t{retrieved_unip}\n")

save_uniprot.close()

sys.stderr.write("\n")
# check any missing proteins

save_uniprot = pd.read_csv(f"{pdbpath.split('/')[-1]}_process_uniprot_chains.txt", sep="\t")
bar = Bar('Checking mappings...', max=len(pdbs_lst))
for pdbname in pdbs_lst: 
    bar.next()
    pdbcode = pdbname.split(".")[0]
    pdb = open(f"{pdbpath}/{pdbname}").readlines()
    if len(pdb) == 0:
        catch_html_err.write(f"{pdbname}: PDB has zero lines. Please verify the file is valid.\n")
        continue
    if pdbcode not in save_uniprot.pdb.values:
        outfile = open(f"{pdbpath.split('/')[-1]}_process_uniprot_chains.txt", "a")
        outfile.write(f"{pdbcode}\tno_uniprot\tno_uniprot\n")
        outfile.close()


print("\n\nFinished.")
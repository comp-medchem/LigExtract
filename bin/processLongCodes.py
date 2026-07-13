# python3

#
# --- DESCRIPTION ----
# This script converts all 5-letter codes in cif files into the short letter code attributed by BeEM during conversion
# back into PDB
# This is a temporary solution, because ideally the 5-letter code should be handled as is during the whole workflow

from glob import glob
#from biopandas.pdb import PandasPdb
import pandas as pd
import os, sys
import string
from pdbecif.mmcif_io import CifFileReader, CifFileWriter
import numpy as np
from progress.bar import Bar

# handle ligand codes (e.g. 9f77)
# convert large Ligand names into BeEM names (e.g. A1AZR --> 01)
print("Screening for ligand-id-mapping files to convert IDs into short form.")

bar = Bar('Processing IDs... ', max=len(glob("cifs/*ligand-id-mapping.tsv")))
for pdb_convert in glob("cifs/*ligand-id-mapping.tsv"):
    bar.next()
    pdb = pdb_convert.split("-")[0].split("/")[-1]
    print(pdb)

    pdb_lig_convert = [ln.strip().split("\t") for ln in open(pdb_convert).readlines()]
    pdb_lig_convert = pd.DataFrame(pdb_lig_convert[1:], columns=pdb_lig_convert[0])
    conversionDict = {lng:short for short,lng in pdb_lig_convert.values}

    cifdata = CifFileReader().read(f'cifs/{pdb}.cif')
    data = cifdata[pdb.upper()]
    found4replacement = []
    for tbl in data:
        data_tbl = pd.DataFrame.from_dict(data[tbl], orient="index").T
        #if tbl == "_chem_comp_atom": break
        #sanity check: are the terms to be replaced actually there
        if np.all([(data_tbl == x).any().any() for x in conversionDict.keys()]):
            found4replacement.append(tbl)
            data_tbl = data_tbl.replace(conversionDict)
            # replace in the original data subsection
            data[tbl] = data_tbl.to_dict(orient='list')
    # replace the full data, but first check if replacement has already been done
    if len(glob(f'cifs/{pdb}.log'))>0:
        if pdb in open(f'cifs/{pdb}.log').read().strip():
            continue # no replacement is needed - already done in a previous pass
    if len(found4replacement)==0:
        # no points of replacement were found. HALT
        sys.stderr.write(f"interrupt process. {pdb} calls for code replacement in ligands but no place for replacement was found")
    cifdata[pdb.upper()] = data
    # overwrite the cif file
    newciffile = CifFileWriter(f'cifs/{pdb}.cif')
    newciffile.write(cifdata)
    #print(f"ligand ID(s) converted for {pdb}")
    logfile = open(f'cifs/{pdb}.log',"w")
    logfile.write(f"ligand ID(s) converted for {pdb}")
    logfile.close()


# handle chain codes (Currently does nothing - in process...)

# from direct screening from 335 files from BeEM
place2replaceChain = ['_atom_site-auth_asym_id','_pdbx_modification_feature-auth_asym_id', '_pdbx_nonpoly_scheme-pdb_strand_id',
       '_pdbx_poly_seq_scheme-pdb_strand_id', '_pdbx_struct_mod_residue-auth_asym_id',
       '_pdbx_struct_sheet_hbond-range_1_auth_asym_id', '_pdbx_struct_sheet_hbond-range_2_auth_asym_id',
       '_pdbx_unobs_or_zero_occ_residues-auth_asym_id', '_pdbx_validate_close_contact-auth_asym_id_1',
       '_pdbx_validate_close_contact-auth_asym_id_2','_pdbx_validate_torsion-auth_asym_id',
       '_refine_ls_restr_ncs-pdbx_auth_asym_id','_struct_conf-beg_auth_asym_id', '_struct_conf-end_auth_asym_id',
       '_struct_conn-ptnr1_auth_asym_id', '_struct_conn-ptnr2_auth_asym_id', '_struct_ncs_dom-details',
       '_struct_ncs_dom_lim-beg_auth_asym_id', '_struct_ncs_dom_lim-end_auth_asym_id', '_struct_ref_seq-pdbx_strand_id',
       '_struct_ref_seq_dif-pdbx_pdb_strand_id', '_struct_sheet_range-beg_auth_asym_id',
       '_struct_sheet_range-end_auth_asym_id', "_entity_poly-pdbx_strand_id", "_pdbx_poly_seq_scheme-auth_mon_id",
       '_atom_site-auth_atom_id', '_atom_site-id','_atom_site-label_alt_id', '_atom_site-label_asym_id',
       '_atom_site-label_atom_id', '_atom_site-type_symbol','_atom_site_anisotrop-id',
       '_atom_site_anisotrop-pdbx_auth_asym_id','_atom_site_anisotrop-pdbx_auth_atom_id',
       '_atom_site_anisotrop-pdbx_label_alt_id','_atom_site_anisotrop-pdbx_label_asym_id',
       '_atom_site_anisotrop-pdbx_label_atom_id','_atom_site_anisotrop-type_symbol', '_atom_type-symbol',
       '_chem_comp_atom-atom_id', '_chem_comp_atom-pdbx_stereo_config','_chem_comp_atom-type_symbol', '_chem_comp_bond-atom_id_1',
       '_chem_comp_bond-atom_id_2','_ndb_struct_na_base_pair-i_auth_asym_id',
       '_ndb_struct_na_base_pair-j_auth_asym_id','_ndb_struct_na_base_pair_step-i_auth_asym_id_1',
       '_ndb_struct_na_base_pair_step-i_auth_asym_id_2','_ndb_struct_na_base_pair_step-j_auth_asym_id_1',
       '_ndb_struct_na_base_pair_step-j_auth_asym_id_2','_pdbx_branch_scheme-auth_asym_id',
       '_pdbx_branch_scheme-pdb_asym_id','_pdbx_distant_solvent_atoms-auth_asym_id','_pdbx_modification_feature-label_alt_id',
       '_pdbx_modification_feature-label_asym_id','_pdbx_modification_feature-modified_residue_auth_asym_id',
       '_pdbx_modification_feature-modified_residue_label_alt_id','_pdbx_modification_feature-modified_residue_label_asym_id',
       '_pdbx_nonpoly_scheme-asym_id', '_pdbx_nonpoly_scheme-auth_mon_id','_pdbx_poly_seq_scheme-asym_id',
       '_pdbx_refine_tls_group-beg_auth_asym_id','_pdbx_refine_tls_group-end_auth_asym_id',
       '_pdbx_struct_conn_angle-ptnr1_auth_asym_id','_pdbx_struct_conn_angle-ptnr1_label_alt_id',
       '_pdbx_struct_conn_angle-ptnr1_label_asym_id', '_pdbx_struct_conn_angle-ptnr2_auth_asym_id',
       '_pdbx_struct_conn_angle-ptnr3_auth_asym_id',  '_pdbx_struct_conn_angle-ptnr3_label_alt_id',
       '_pdbx_struct_conn_angle-ptnr3_label_asym_id','_pdbx_struct_sheet_hbond-range_1_label_asym_id',
       '_pdbx_struct_sheet_hbond-range_2_label_asym_id', '_pdbx_struct_special_symmetry-auth_asym_id',
       '_pdbx_struct_special_symmetry-label_asym_id', '_pdbx_unobs_or_zero_occ_atoms-auth_asym_id',
       '_pdbx_unobs_or_zero_occ_atoms-label_alt_id', '_pdbx_unobs_or_zero_occ_atoms-label_asym_id',
       '_pdbx_unobs_or_zero_occ_residues-label_asym_id', '_pdbx_validate_chiral-auth_asym_id',
       '_pdbx_validate_close_contact-label_alt_id_1','_pdbx_validate_main_chain_plane-auth_asym_id',
       '_pdbx_validate_peptide_omega-auth_asym_id_1','_pdbx_validate_peptide_omega-auth_asym_id_2',
       '_pdbx_validate_planes-auth_asym_id','_pdbx_validate_polymer_linkage-auth_asym_id_1',
       '_pdbx_validate_polymer_linkage-auth_asym_id_2','_pdbx_validate_rmsd_angle-auth_asym_id_1',
       '_pdbx_validate_rmsd_angle-auth_asym_id_2','_pdbx_validate_rmsd_angle-auth_asym_id_3',
       '_pdbx_validate_rmsd_bond-auth_asym_id_1','_pdbx_validate_rmsd_bond-auth_asym_id_2',
       '_pdbx_validate_symm_contact-auth_asym_id_1','_pdbx_validate_symm_contact-auth_asym_id_2',
       '_pdbx_validate_torsion-label_alt_id', '_struct_asym-id',
       '_struct_conf-beg_label_asym_id', '_struct_conf-end_label_asym_id',
       '_struct_conn-pdbx_ptnr1_label_alt_id','_struct_conn-pdbx_ptnr2_label_alt_id',
       '_struct_conn-ptnr1_label_asym_id','_struct_conn-ptnr2_label_asym_id', '_struct_mon_prot_cis-auth_asym_id',
       '_struct_mon_prot_cis-label_asym_id','_struct_mon_prot_cis-pdbx_auth_asym_id_2',
       '_struct_mon_prot_cis-pdbx_label_asym_id_2', '_struct_ncs_dom_lim-beg_label_asym_id',
       '_struct_ncs_dom_lim-end_label_asym_id','_struct_sheet_range-beg_label_asym_id',
       '_struct_sheet_range-end_label_asym_id','_struct_site-pdbx_auth_asym_id', '_struct_site_gen-auth_asym_id']
place2replaceChain = [x.split("-") for x in place2replaceChain]

tbl_col_replace = []
'''
for chains2convert in glob("cifs/*chain-id-mapping.txt"):
    # WARNING - need to handle files that get split into two pdbs!
    pdb = chains2convert.split("-")[0].split("/")[-1]
    chains2convert = [ln.strip().split() for ln in open(chains2convert).readlines()]
    chains2convert = pd.DataFrame(chains2convert[3:], columns=["NewchainID","OriginalchainID"])
    chains2convert = chains2convert.query("NewchainID != OriginalchainID")
    chainconversionDict = {orig:new for new,orig in chains2convert.values}
    cifdata = CifFileReader().read(f'cifs/{pdb}.cif')
    data = cifdata[pdb.upper()]
    for tbl,col in place2replaceChain:
        if tbl in data:
            data_tbl = pd.DataFrame.from_dict(data[tbl], orient="index").T
            newcol = data_tbl[[col]].replace(chainconversionDict)
'''
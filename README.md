# LigExtract: Automated Ligand Identification and Extraction from PDB Structures

![](docs/sources/images/ligextract_logo.png)


Software that allows large-scale ligand extraction from UniProt ID queries. 

Developed by **Natália Aniceto** (ORCID 0000-0001-7039-0022), **Nuno Martinho** (ORCID 0000-0001-5102-4756) and **Rita Guedes** (ORCID 0000-0002-5790-9181).

Below is the overall workflow of LigExtract:


![](docs/sources/images/scheme_app_nologo.png)


## Dependencies

LigExtract was developed on Linux and written in Python 3. It requires PyMol to be installed and, to avoid issues with properly interfacing with a system-wide PyMOL, the ligextract env has its own pymol installation.

Once you have installed the conda env provided and have activated it, you can run the command below to check if PyMOL can be used properly within the LigExtract environment.

    pymol -cq

If this command returns an error indicating that PyMOL cannot be properly used, go to the troubleshooting section.


## Running LigExtract

**1.** git clone LigExtract into your home directory. Make ligextract.sh and build_dependencies.sh inside of bin executable (i.e. chmod 755 *.sh). Add path to LigExtract/bin to your $PATH.

**2.** create a conda environment from ligextract.yml.

**3.** cd into your working directory. This is where all PDBs will be downloaded and processed. This directory must contain a file with a name following the format 

        <projectname>_uniprot_list.txt

For example, my project is called "myproteins" so the file will be named
        
        myproteins_uniprot_list.txt

This file will contain a list of UniProt IDs (see example in docs)

**4.** Build dependencies (one-off downloads):
        
        build_dependencies.sh

**5.** Run LigExtract for your query proteins in your *_uniprot_list.txt file (without cleaning-up at the end).

        ligextract.sh -d myproteins -r 3.5 -o cluster -c nan

with clean-up at the end (i.e. removing all *.cif files):

        ligextract.sh -d myproteins -r 3.5 -o filter -c clean

  In this example, ligextract will only consider PDBs up to 3.5 Angstrom resolution and will employ the "cluster" mode (i.e., all ligands that survive filtration are kept, even if duplicated)
  
  Notice how "myproteins" is the name provided to the -d argument, as this must correspond to the prefix of the *_uniprot_list.txt file.
  

#### Arguments of ligextract.sh:
     -h    usage information
     -d    name of the directory that will be created with all PDBs. This will be the prefix for multiple files.
     -r    maximum PDB resolution accepted
     -o    selected ligand selection mode: can be 'filter' or 'cluster'
     -c    cleaning outcome: 'clean' will delete all *.cif files; 'nan' will keep all *.cif files.

     
#### Usage Notes:

- LigExtract will prompt the user to manually reject or accept each of all experiment methods used to obtain the PDBs associated with the UniProt queries.
- LigExtract will bypass entries with no PDB-format file available. Read more on what causes this situation on https://www.rcsb.org/docs/general-help/structures-without-legacy-pdb-format-files. Eventually we plan to accommodate these exceptions but at the moment LigExtract is not capable on handling them.


#### Troubleshooting

If you cannot run PyMOL from the command line once LigExtract's environment is activate, you can try these alternative conditions:

1) You should replace `cluster_ligands_hierarchical.py` inside `bin` by `cluster_ligands_hierarchical_fixPymol.py` by simply renaming the first script or moving it out of the bin directory, and removing the "_fixPymol" portion of the second script. In order for `cluster_ligands_hierarchical_fixPymol.py` to work, you must make sure that `source deactivate` deactivates LigExtract's environment without reverting into any base environment.

2) If you still have an error when running `cluster_ligands_hierarchical.py` ("FileNotFoundError: No such file or directory: 'myproteins/pdbs_filtered_chains/PXXXX/aligned'"), try to replace `deactivate` by `{HOME}/anaconda3/bin/deactivate` inside that script (3 locations in the script). Note that this command should be used exactly as is, because `HOME` is already a variable in the python script that grabs your own home path. This is the typical path to the `deactivate` executable, and if your path is different you must replace it.

3) If you still get that error when running `cluster_ligands_hierarchical.py`, replace `cluster_ligands_hierarchical.py` by `cluster_ligands_hierarchical_fixPymol_hideconda.py`, making sure you rename the second script with the name of the first.

#### _________
If you encounter any errors or issues, or if you have any suggestions, please email Natalia Aniceto at: nataliaaniceto[at]ff.ul.pt

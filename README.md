# LigExtract

### Automated Ligand Identification and Extraction from PDB Structures

- **Version: 1.1 (22 May 2025)**

![](docs/sources/images/ligextract_logo.png)


Software that allows large-scale ligand extraction from UniProt ID queries. 



The original publication describing this tool in depth, including performance benchmarks, is [here](https://academic.oup.com/gpb/advance-article/doi/10.1093/gpbjnl/qzaf018/8046017)


Below is the overall workflow of LigExtract:


![](docs/sources/images/scheme_app_nologo.png)


## Dependencies

LigExtract was developed on Linux and written in Python 3. It requires PyMol to be installed and, to avoid issues with properly interfacing with a system-wide PyMOL, the ligextract env has its own pymol installation.

Once you have installed the conda env provided and have activated it, you can run the command below to check if PyMOL can be used properly within the LigExtract environment.

    pymol -cq


## Running LigExtract

**1.** git clone LigExtract into your home directory. Make ligextract.sh and build_dependencies.sh inside of bin executable (i.e. chmod 755 *.sh). Add path to LigExtract/bin to your $PATH environment.

**2.** create a conda environment from ligextract.yml, and activate it.

**3.** cd into your working directory. This is where all PDBs will be downloaded and processed. This directory must contain a file with a name following the format 

        <projectname>_uniprot_list.txt

For example, my project is called "myproteins" so the file will be named
        
        myproteins_uniprot_list.txt

This file will contain a list of UniProt IDs (see example in docs)

**4.** Build dependencies (one-off downloads):
        
        build_dependencies.sh

**5.** Run LigExtract for your query proteins in your *_uniprot_list.txt file (without cleaning-up at the end). You can use the *cluster* mode if you want to keep all ligands (recommended for molecular docking, binding sites study, etc):

        ligextract.sh -d myproteins -r 2.5 -o cluster -c no


Cluster mode with clean-up at the end (i.e. removing all raw *.pdb files converted from *.cif):

        ligextract.sh -d myproteins -r 2.5 -o cluster -c yes

  In this second example, ligextract will only consider PDBs up to 2.5 Angstrom resolution and will employ the "cluster" mode (i.e., all ligands that survive filtration are kept, even if duplicated). All raw PDB files will be removed after the process is finished with **-c yes**.
  
  Notice how "myproteins" is the name provided to the -d argument, as this must correspond to the prefix of the *_uniprot_list.txt file.
  

Additionally you can provide your own list of PDBs. This is meant to make the query more efficient in cases where you do not want to consider/process all PDBs that map to your UniProt ID(s). This should be a simple *.txt file with one PDB code per line (not case sensitive).

        ligextract.sh -d myproteins -r 3 -f myPdbQueries.txt


You can inspect the arguments available with:

        ligextract.sh -h


*In the example query provided in docs/myproteins_uniprot_list.txt, LigExtract processes 267 structures in 6m32s using 7 CPUs.*

## Outputs

#### Cluster mode:

Ligextract produces a table called **projectname_ligandsList.txt** with all ligands and some data characterising them, looking like this:

ligandfile | pocketres_chain | pocketres_chain_size | chain_name | ligtype | lig_ID | pdbcode 
--- | --- | --- | --- | --- | --- | --- 
1sb1_lig_chain-I.pdb | ARG67-H;(...);TYR76-H | 18 | H | chain ligand | 1sb1_lig_chain-I | 1sb1 
1sb1_chain-H_lig-165-1001.pdb | ALA190-H;(...);VAL213-H | 30 | H | small-molecule ligand | 165-1001 | 1sb1


In the example shown above, structure [1SB1](https://www.rcsb.org/structure/1SB1) has two ligands: one peptide ligand annotated under chain I, and a small-molecule ligand with ID 165 (and residue number 1001). Both ligands are bound to chain H.

The of identified ligands in the output file, and their corresponding cleaned proteins, are stored under **projectname/pdbs_filtered_chains/uniprotQuery/aligned_pdbs**. Here, *projectname* corresponds to the name you provided earlier, and *uniprotQuery* will correspond to each query in your input file. In this directory, a **pymol session file** (*.pse) is also saved with all ligands clustered (each cluster has a color and a code registered in the **clusters** directory).

All structures inside a given uniprot query are aligned and saved separately (ligands and proteins in their own individual files). The **unaligned complexes** are also saved one level up, in **projectname/pdbs_filtered_chains/uniprotQuery**.

The original, raw list of ligands after the first pass (module 1) of ligand identification is saved in **rawlist_extraction.txt**.



## Additional Information:
  

#### Arguments of ligextract.sh:
     -h    usage information
     -d    name of the directory that will be created with all PDBs. This will be the prefix for multiple files. (required)
     -r    maximum PDB resolution accepted (default=2.5)
     -o    selected ligand selection mode: can be 'filter' or 'cluster' (default='cluster')
     -c    cleaning outcome: 'yes' will delete cifs directory; 'no' will keep all *.cif files. (default='yes')
     -f    file with a list of PDB IDs to use (optional)
     -v    See installed version

     
#### Usage Notes:

- LigExtract will prompt the user to manually reject or accept each of all experiment method used to obtain the PDBs associated with the UniProt queries.
- LigExtract will replace 5-character ligand codes with a short numerical code (e.g.,A1AIM -> 01) everywhere, through the use of BeEM (mmCIF-to-pdb converter). The user can easily access the original ID in the *-ligand-id-mapping.tsv files for the corresponding PDB, inside cifs.



#### _________

Developed by **Natália Aniceto** (ORCID 0000-0001-7039-0022), **Nuno Martinho** (ORCID 0000-0001-5102-4756) and **Rita Guedes** (ORCID 0000-0002-5790-9181).

If you encounter any errors or issues, or if you have any suggestions, please email Natalia Aniceto at: nataliaaniceto [at]ff.ul.pt

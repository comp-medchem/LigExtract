# LigExtract: Automated Ligand Indentification and Extraction from PDB Structures

![](docs/sources/images/ligextract_logo.png)


Software that allows large-scale ligand extraction from UniProt ID queries. Below is the overall workflow of LigExtract:


![](docs/sources/images/scheme_app_nologo.png)


## Running LigExtract

1. Add path to LigExtract/ligextract directory to $PATH

2. cd into your working directory. This is where all PDBs will be downloaded and processed. This directory must contain a file with a name following the format 

        <projectname>_uniprot_list.txt

For example, my project is called "myproteins" so the file will be named
        
        myproteins_uniprot_list.txt

3. Build dependencies (one-off downloads):
        
        bash build_dependencies.sh

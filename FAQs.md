# Frequently Asked Questions

### What is considered to be a ligand in LigExtract?

In theory any molecule bound to a protein is a ligand. However, LigExtract is aimed at facilitating drug design efforts, so it only considers as ligands molecules that fall the profile a molecule of pharmaceutical interest, namely, small molecules that are non-ubiquitous accross the Protein Data Bank (i.e. do not occur in more than 250 PDBs), that have more than 3 heavy atoms and which do not occur more than 10 times within a given PDB structure under analysis. Some exceptions to this are molecules such as ADP that are ubiquitous but still could fall under the drug-like type. A few of these exceptions are gathered in a list whose ligands will bypass removal during the filtration stage.

### How does LigExtract handle 5-character ligand IDs?

LigExtract still relies on the PDB file format, so all cif files (originally downloaded in that format) are converted to the pdb using BeEM. BeEM converts these long codes to short codes such as "01", and stores this conversion in a file ending with *-ligand-id-mapping.tsv*, that is stored inside the cifs directory. To keep everything consistent, the long codes in original cif file are also converted to the short code assigned by BeEM. However the final output file (ending with \*_ligandsList.txt) has the original code.


### my protein of interest contains an ion that plays a modulatory role (therefore it should be considered a ligand of pharmaceutical interest). Why is it not considered as a ligand? 
Several cases exist of ions that play a modulatory role. In some cases an ion ligand even acts as an inhibitor. For example, Fluoride inhibits enolase. However, the main focus of LigExtract if to aid drug design (often aimed at gathering ligands for molecular docking, for instance) so, it considers only drug-like molecules, excluding all others. One can easily access the ions listed for each PDB using RSBC's own lists (built in the dependencies of LigExtract, under data).


### How does LigExtract handle structures with no ligands of pharmaceutical interest?

As this tool is focused on ligands if, at any stage, a given pdb is found to have no ligands it is removed.

### Why does LigExtract sometimes keep more than one small molecule?

In some cases it is not possible to know, without reading the original publication and/or having knowledge on crystallography procedures, whether a given molecule in a ligand (of pharmaceutical interest). If we consider the case of 5F1A, LigExtract will classify SAL, AKR and COH as possible ligands.
- In this case SAL (salicylic acid) is the key ligand, according to the title of the PDB entry ("The Crystal Structure of Salicylate Bound to Human Cyclooxygenase-2"). 
- Regarding AKR (acrylic acid), this is a molecule that occurs in fewer PDBs than SAL and is complexed in just 3 locations in the PDB structure (which can legitimatelly happen for inhibitors of certain kinases, for instance), as is not overly small (similar size to SAL). As a result, there is no reason to exclude AKR as, by all accounts, it looks like it could be a ligand.
- COH (protoporphyrin IX) can be considered a second legitimate ligand, as it is a cofactor (and could be used to drive drug design, since cofactor-competitive inhibitors are a viable route of drug development).
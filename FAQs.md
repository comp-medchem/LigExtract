# Frequently Asked Questions


### How does LigExtract handle 5-character ligand IDs?

LigExtract still relies on the PDB file format, so all cif files (originally downloaded in that format) are converted to the pdb using BeEM. BeEM converts these long codes to short codes such as "01", and stores this conversion in a file ending with *-ligand-id-mapping.tsv*, that is stored inside the cifs directory. To keep everything consistent, the long codes in original cif file are also converted to the short code assigned by BeEM. However the final output file (ending with \*_ligandsList.txt) has the original code.


### An ion is a ligand that plays a modulatory role in my protein of interest. Why is it not considered as a ligand? 
Several cases exist of ions that play a modulatory role. In some cases an ion ligand even acts as an inhibitor. For example, Fluoride inhibits enolase. However, the main focus of LigExtract if to aid drug design (often aimed at gathering ligands for molecular docking, for instance) so, it considers only drug-like molecules, excluding all others.


### How does LigExtract handle structures with no ligands of pharmaceutical interest?

As this tool is focused on ligands if, at any stage, a given pdb is found to have no ligands it is removed.


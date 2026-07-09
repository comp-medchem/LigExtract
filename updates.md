## Updates


### 23 April 2025
	- Fix bugs related to absent fields in the CIF file
	- Add multithreading to Module 2 and 4 to increase efficiency
	- Accomodate cases where multiple HETATM residues are covalently bound but there is no covalent instances in _struct_conn (e.g 6FJ2) - bind if minimum distance <= 1.7A



### 9 July 2026
    - replaced Ligand Expo (no longer maintained) with ligand information directly from PDB (https://github.com/rcsb/rcsb-training-resources/blob/master/example-use-cases/pdb-ligand-composition/generate_pdb_ligand_mappings.py) as per this PDB post: https://www.rcsb.org/docs/general-help/ligandexpo-shutdown

    - Add functionality to allow for PDB queries (DONE)
    - multiple fixes to cater to exceptions found


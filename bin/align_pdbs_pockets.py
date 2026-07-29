import sys, copy, os
import numpy
filesloc = sys.argv[1] 
targetdir = sys.argv[1]
refpdb = sys.argv[2]
if targetdir[-1]=="/": targetdir = targetdir[:-1]

if targetdir=='':
	targetdir = "aligned_pdbs/"
else:
	targetdir = targetdir+"/aligned_pdbs/"

try: os.mkdir(targetdir)
except OSError: None


listpdbs = [x for x in os.listdir(filesloc) if (x.endswith(".pdb") and x!=refpdb)]

# aligment
#refpdb = listpdbs[0]
cmd.load(filesloc+"/"+refpdb)
# ref pdb is already one of the files that will be aligned.
#cmd.save(targetdir+refpdb.split(".pdb")[0]+"_align.pdb", refpdb.split(".pdb")[0])

rms_align = []
for f in listpdbs:
	cmd.load(filesloc+"/"+f)
	al=cmd.align(f.split(".pdb")[0],refpdb.split(".pdb")[0])
	rmsd = al[0]
	if rmsd>3.5:
		al=cmd.cealign(refpdb.split(".pdb")[0],f.split(".pdb")[0])
		rmsd = al["RMSD"]
	rms_align.append([f.split(".pdb")[0],rmsd])
	# save aligned PDBs
	cmd.save(targetdir+f.split(".pdb")[0]+"_align.pdb", f.split(".pdb")[0])

# Save picture
for i in rms_align: print(i)




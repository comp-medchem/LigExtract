import sys, copy, os
import numpy as np

workdir = sys.argv[1]
pockets_file = sys.argv[2]

outfile = open("{}_pocket_clusters_colorcodes.txt".format(pockets_file.split("_")[0]),"w")

if workdir[-1]=="/": workdir = workdir[:-1]
color_lst = [x[0] for x in cmd.get_color_indices()[2:26] if x[0]!="dash"]

pockets_file = [x.strip().split("\t") for x in open(pockets_file).readlines()]
pockets_file_head = np.array(pockets_file[0])
column_idx = {x:i for i,x in enumerate(pockets_file_head)}
pockets_file = pockets_file[1:]
unique_pockets = np.unique([x[column_idx["cluster"]] for x in pockets_file])

outfile.write("color\tcluster_id\n")
for color_i,cluster in enumerate(unique_pockets):
	cluster_ligs = [x[column_idx["ligandfile"]] for x in pockets_file if x[column_idx["cluster"]]==cluster]
	for l in cluster_ligs:
		cmd.load(workdir+"/"+l, l.split(".")[0])
		cmd.do("color {}, {}".format(color_lst[color_i], l.split(".")[0]))
	# group clusters in PyMOL
	ligs2group=" ".join([l.split(".")[0] for l in cluster_ligs])
	cmd.do("group cluster{}, {}".format(cluster,ligs2group))
	outfile.write("{}\t{}\n".format(color_lst[color_i],cluster))

outfile.close()

# PyMOL session
cmd.show("sticks")
cmd.center()
#cmd.zoom(complete=1)
cmd.do("zoom all, 1")
cmd.save("{}/all_clusters.pse".format(workdir))

# rotate 45* right
#cmd.do("turn x,45")
#cmd.png("{}/all_pockets_2.png".format(workdir), ray=1, dpi=600)
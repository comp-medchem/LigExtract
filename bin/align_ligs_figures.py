import sys, copy, os
import numpy

pocketname = sys.argv[1]
workdir = sys.argv[2]
files = sys.argv[3:]

if workdir[-1]=="/": workdir = workdir[:-1]

for f in files:
	cmd.load(workdir+"/"+f)

cmd.show("sticks")
cmd.center()
#cmd.zoom(complete=1)
cmd.do("zoom all, 2")

cmd.png("{}/{}.png".format(workdir,pocketname), ray=1, dpi=100)

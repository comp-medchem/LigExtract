#!/bin/bash

# THIS SCRIPT WAS TAKEN FROM https://www.rcsb.org/docs/programmatic-access/batch-downloads-with-shell-script
# Script to download files from RCSB http file download services.
# Use the -h switch to get help on usage.

if ! command -v curl &> /dev/null
then
    echo "'curl' could not be found. You need to install 'curl' for this script to work."
    exit 1
fi

PROGNAME=$0
BASE_URL="https://files.rcsb.org/download"

usage() {
  cat << EOF >&2
Usage: $PROGNAME -f <file> [-o <dir>] [-c] [-p]

 -f <file>: the input file containing a comma-separated list of PDB ids
 -o  <dir>: the output dir, default: current dir
 -c       : download a cif.gz file for each PDB id
 -p       : download a pdb.gz file for each PDB id (not available for large structures)
 -a       : download a pdb1.gz file (1st bioassembly) for each PDB id (not available for large structures)
 -x       : download a xml.gz file for each PDB id
 -s       : download a sf.cif.gz file for each PDB id (diffraction only)
 -m       : download a mr.gz file for each PDB id (NMR only)
 -r       : download a mr.str.gz for each PDB id (NMR only)
EOF
  exit 1
}

function ProgressBar {
    let _progress=(${1}*100/${2}*100)/100
    let _done=(${_progress}*4)/10
    let _left=40-$_done
    _fill=$(printf "%${_done}s")
    _empty=$(printf "%${_left}s")
    printf "\rProgress : [${_fill// /#}${_empty// /-}] ${_progress}%%"

}



download() {
  url="$BASE_URL/$1"
  out=$2/$1
  #echo "Downloading $1 to $out"
  curl -s -f $url -o $out || echo "Failed to download $url" >> pdb_download.log
}

listfile=""
outdir="."
cif=false
pdb=false
pdb1=false
xml=false
sf=false
mr=false
mrstr=false
while getopts f:o:cpaxsmr o
do
  case $o in
    (f) listfile=$OPTARG;;
    (o) outdir=$OPTARG;;
    (c) cif=true;;
    (p) pdb=true;;
    (a) pdb1=true;;
    (x) xml=true;;
    (s) sf=true;;
    (m) mr=true;;
    (r) mrstr=true;;
    (*) usage
  esac
done
shift "$((OPTIND - 1))"

if [ "$listfile" == "" ]
then
  echo "Parameter -f must be provided"
  exit 1
fi
contents=$(cat $listfile)

# see https://stackoverflow.com/questions/918886/how-do-i-split-a-string-on-a-delimiter-in-bash#tab-top
IFS=',' read -ra tokens <<< "$contents"


#for number in $(seq ${_start} ${_end})
#do ProgressBar ${number} ${_end}
#done
numpdbs=`echo $contents | tr -cd , | wc -c`
_start=1
_end=$numpdbs

> pdb_download.log
number=0

for token in "${tokens[@]}"
do
  if [ "$cif" == true ]
  then
    number=$(($number + 1))
    ProgressBar ${number} ${_end}
    download ${token}.cif.gz $outdir
  fi
done









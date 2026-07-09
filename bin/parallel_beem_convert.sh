printf "\nConvert mmCIF to pdb\n"

rootdir=$1
d=$2

# how many conversions to run at once
JOBS=$(nproc --all)
JOBS=$((JOBS * 3 / 4))
((JOBS<1)) && JOBS=1

running=0

for cif in cifs/*.cif; do
(
    pdbcode=$(basename "$cif" .cif)
    pdbname="${pdbcode}.pdb"

    if [ ! -f "$d/$pdbname" ]; then
        "$rootdir/LigExtract/bin/BeEM.linux" "$cif"

        # if BeEM produces *chain-id-mapping.txt file, change output pdb file's name back to <pdbname>.pdb
        if [ -f "${pdbcode}-chain-id-mapping.txt" ]; then
            # if more than 1 pdb bundle, bypass this for now
            if [ "$(ls "${pdbcode}"-pdb-bundle*.pdb 2>/dev/null | wc -l)" -gt 1 ]; then
                rm -f "${pdbcode}"-*
                exit 0
            fi
            mv "${pdbcode}-pdb-bundle1.pdb" "$pdbname"
        fi

        mv "$pdbname" "$d/."
    fi
) &

    running=$((running+1))
    if [ "$running" -ge "$JOBS" ]; then
        # wait for any one job to finish, then continue launching more
        wait -n
        running=$((running-1))
    fi
done

# wait for remaining jobs
wait

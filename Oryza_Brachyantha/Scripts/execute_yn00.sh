for d in cluster_3526/*.ali.phy; do
    a=$(basename "$d")
    pair="${a%.ali.phy}"
    newDir="cluster_3526/${pair}"

    mkdir -p "$newDir"

    awk -v file="$d" -v out="$newDir/yn" '
        {
            gsub("XXXXX", file)
            gsub("outfile *= *yn", "outfile = " out)
            print
        }
    ' yn00.ctl_master > "$newDir/yn00.ctl"
    yn00 $newDir/yn00.ctl
    mv -i 2YN.dN 2YN.dS 2YN.t $newDir
    mv -i rst rst1 rub $newDir
done


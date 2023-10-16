#!/bin/bash

#####################
#
# This script takes the flag files meticulously created by Omar
# and applies them to the scaled, SPC'ed, and SFXC archives with 'paz'.
#
# This script assumes that r67_generate_archives.sh was run previously
#
# Once all the flags are applied, we can start to compute fluences, S/N, widths, TOAs...
#
####################


pwait() {
    while [ $(jobs -p | wc -l) -ge $1 ]; do
	sleep 0.33
    done
}
njobs=16

flagrundir=./flags/
# 
outbasedir=<root-dir-where-archives-live>

runs=( \
 exp dish scan flag
"r67008 o8 17" \
"r67011 o8 17" \
"r67014 o8 46" \
"r67014 o8 13" \
"r67022 o8 27" \
"r67022 o8 28" \
"r67028 o8 20" \
"r67l01 tr 29" \
"r67030 o8 16" \
"r67079 o8 27" \
"r67079 o8 31" \
"r67081 o8 05" \
"r67081 o8 07" \
"p67112 wb 20" \
"pre027 wb 18" \
"pre027 wb 22" \
"pre028 wb 11" \
"p67113 wb 24" \
"p67113 wb 29" \
"p67114 wb 17" \
"p67116 wb 32" \
"p67117 wb 39" \
"p67118 wb 34" \
"p67118 o8 34" \
"p67118 wb 46" \
"p67118 o8 46" \
"p67118 wb 47" \
"p67118 o8 47" \
"p67119 o8 03" \
"p67120 wb 18" \
"p67120 wb 27" \
"p67121 wb 39" \
"p67124 o8 35" \
"p67124 o8 39" \
"p67125 wb 18" \
"p67125 wb 22" \
"pre033 wb 32" \
"pre033 wb 38" \
"p67128 wb 21" \
"p67128 wb 24" \
"p67128 wb 40" \
"p67138 wb 25" \
"p67138 wb 40" \
"r67l29 tr 19" \
"stocke st 01" \
"stocke st 02" \
"stocke st 03" \
"stocke st 04" \
"stocke st 05" \
"stocke st 06" \
"stocke st 07" \
"stocke st 08" \
"stocke st 09" \
"stocke st 10" \
"stocke st 11" \
"stocke st 12" \
"stocke st 13" \
"stocke st 14" \
"stocke st 15" \
"stocke st 16" \
)
for run in "${runs[@]}";do
    IFS=" " read -r -a info <<< "${run}"
    exp="${info[0]}"
    dish="${info[1]}"
    scan="${info[2]}"
    arbase=${exp}_${dish}_no00${scan}
    flagfile=${arbase}.sh
    expbasedir=${outbasedir}/${exp}
    sfxc_file=${expbasedir}/sfxc/${arbase}.ds.ar
    scale_file=${expbasedir}/scale_spc/${arbase}_allIFs.vdif_pol2.fil.scale.ds.ar
    spc_file=${expbasedir}/scale_spc/${arbase}_allIFs.vdif_pol2.fil.scale.spc.ds.ar
    fs="${sfxc_file} ${scale_file} ${spc_file}"
    for f in ${fs}; do
        if ! [ -e ${f} ];then
            echo "${f} doesn't exist. Moving on..."
            continue
        fi
        cmd="${flagrundir}/${flagfile} ${f}"
        echo "running ${cmd}"
        eval ${cmd}
    done
done

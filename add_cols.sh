#!/bin/bash
set -e

sefd_o8=310
sefd_wb=420
sefd_tr=250
sefd_st=1100
sefd_st2=385

runs=( \
 exp dish scan sefd ID
"r67008 o8 17 ${sefd_o8} B02-o8" \
"r67011 o8 17 ${sefd_o8} B03-o8" \
"r67014 o8 13 ${sefd_o8} B04-o8" \
"r67014 o8 46 ${sefd_o8} B05-o8" \
"r67022 o8 27 ${sefd_o8} B06-o8" \
"r67022 o8 28 ${sefd_o8} B07-o8" \
"r67028 o8 20 ${sefd_o8} B08-o8" \
"r67l01 tr 29 ${sefd_tr} B08-tr" \
"r67030 o8 16 ${sefd_o8} B09-o8" \
"r67079 o8 27 ${sefd_o8} B11-o8" \
"r67079 o8 31 ${sefd_o8} B12-o8" \
"r67081 o8 05 ${sefd_o8} B13-o8" \
"r67081 o8 07 ${sefd_o8} B14-o8" \
"p67112 wb 20 ${sefd_wb} B15-wb" \
"pre027 wb 18 ${sefd_wb} B16-wb" \
"pre027 wb 22 ${sefd_wb} B17-wb" \
"pre028 wb 11 ${sefd_wb} B18-wb" \
"p67113 wb 24 ${sefd_wb} B19-wb" \
"p67113 wb 29 ${sefd_wb} B20-wb" \
"p67114 wb 17 ${sefd_wb} B21-wb" \
"p67119 o8 03 ${sefd_o8} B22-o8" \
"p67116 wb 32 ${sefd_wb} B23-wb" \
"p67117 wb 39 ${sefd_wb} B24-wb" \
"p67118 wb 34 ${sefd_wb} B25-wb" \
"p67118 o8 34 ${sefd_o8} B25-o8" \
"p67118 wb 46 ${sefd_wb} B26-wb" \
"p67118 o8 46 ${sefd_o8} B26-o8" \
"p67118 wb 47 ${sefd_wb} B27-wb" \
"p67118 o8 47 ${sefd_o8} B27-o8" \
"p67120 wb 18 ${sefd_wb} B28-wb" \
"p67120 wb 27 ${sefd_wb} B29-wb" \
"p67121 wb 39 ${sefd_wb} B30-wb" \
"p67124 o8 35 ${sefd_o8} B31-o8" \
"p67124 o8 39 ${sefd_o8} B32-o8" \
"p67125 wb 18 ${sefd_wb} B33-wb" \
"p67125 wb 22 ${sefd_wb} B34-wb" \
"pre033 wb 32 ${sefd_wb} B35-wb" \
"pre033 wb 38 ${sefd_wb} B36-wb" \
"p67128 wb 21 ${sefd_wb} B37-wb" \
"p67128 wb 24 ${sefd_wb} B38-wb" \
"p67128 wb 40 ${sefd_wb} B39-wb" \
"p67138 wb 25 ${sefd_wb} B41-wb" \
"p67138 wb 40 ${sefd_wb} B42-wb" \
"r67l29 tr 19 ${sefd_tr} B42-tr" \
"stocke st 01 ${sefd_st} B01-st" \
"stocke st 02 ${sefd_st} B06-st" \
"stocke st 03 ${sefd_st} B10-st" \
"stocke st 04 ${sefd_st2} B31-st" \
"stocke st 05 ${sefd_st2} B32-st" \
"stocke st 06 ${sefd_st2} B33-st" \
"stocke st 07 ${sefd_st2} B34-st" \
"stocke st 08 ${sefd_st2} B35-st" \
"stocke st 09 ${sefd_st2} B37-st" \
"stocke st 10 ${sefd_st2} B38-st" \
"stocke st 11 ${sefd_st2} B39-st" \
"stocke st 12 ${sefd_st2} B40-st" \
"stocke st 13 ${sefd_st2} B43-st" \
"stocke st 14 ${sefd_st2} B44-st" \
"stocke st 15 ${sefd_st2} B45-st" \
"stocke st 16 ${sefd_st2} B46-st" \
)

dbdir=./dbs/
dm=410.775
distance=453  # Mpc

db="${dbdir}/burst_info.pickle"  # everything goes here
#db="${dbdir}/burst_info_test.pickle"  # everything goes here
add_id='./python/add_column_2_db.py'
for run in "${runs[@]}";do
    IFS=" " read -r -a info <<< "${run}"
    exp="${info[0]}"
    dish="${info[1]}"
    scan=no00"${info[2]}"
    sefd="${info[3]}"
    id="${info[4]}"
    scan_sfxc=${scan}.ds.pazi
    if [[ ${dish} == 'st' ]];then
        scan=${scan}.ds.pazi
    fi
    colnames=( 'id' 'sefd_jy' 'distance_Mpc' 'dm' )
    values=( ${id} ${sefd} ${distance} ${dm} )
    for cnt in `seq 0 3`;do
        colname="${colnames[${cnt}]}"
        value="${values[${cnt}]}"
        unique="${uniques[${cnt}]}"
        cmd="python3 ${add_id} ${db} --exp ${exp} --dish ${dish} --scan ${scan} --colname ${colname} --val ${value}"
        echo "running ${cmd}"
        eval ${cmd}
        cmd="python3 ${add_id} ${db} --exp ${exp} --dish ${dish} --scan ${scan_sfxc} --colname ${colname} --val ${value}"
        echo "running ${cmd}"
        eval ${cmd}
    done
done

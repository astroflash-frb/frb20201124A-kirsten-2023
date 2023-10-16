#!/bin/bash
set -e

##############
#
# This is a wrapper around ./python/get_MJD.py which
# combines the info from the ar-files and the pandas data
# frames to fit the bursts in the SFXC/Stockert-filterbanks with
# Guassians to determine the TOAs. For SFXC-filterbanks these
# are geocentric TOAs, while for Stockert they are topocentric
# at Stockert.
#
#############


ref_freq_o8_4g=1698.0
ref_freq_o8_2g=1608.0
ref_freq_o8_1g=1480.0
ref_freq_o8_1g_off=1578.0
ref_freq_wb_1g=1379.49
ref_freq_tr_2g=1580.0  # r67l01
ref_freq_tr_1g=1470.0  # r67l29
ref_freq_st_1g='none'

# exp dish scan sefd ID ref_freq DM
runs=( \
"r67008 o8 17 coe B02-o8 ${ref_freq_o8_1g} 0.0" \
"r67011 o8 17 coe B03-o8 ${ref_freq_o8_1g} 0.0" \
"r67014 o8 13 o8 B04-o8 ${ref_freq_o8_1g} 410.775" \
"r67014 o8 46 coe B05-o8 ${ref_freq_o8_1g} 0.0" \
"r67022 o8 27 coe B06-o8 ${ref_freq_o8_1g_off} 0.0" \
"r67022 o8 28 coe B07-o8 ${ref_freq_o8_1g_off} 0.0" \
"r67028 o8 20 coe B08-o8 ${ref_freq_o8_1g} 0.0" \
"r67l01 tr 29 coe B08-tr ${ref_freq_tr_2g} 0.0" \
"r67030 o8 16 coe B09-o8 ${ref_freq_o8_1g} 0.0" \
"r67079 o8 27 coe B11-o8 ${ref_freq_o8_1g} 0.0" \
"r67079 o8 31 coe B12-o8 ${ref_freq_o8_1g} 0.0" \
"r67081 o8 05 coe B13-o8 ${ref_freq_o8_1g} 0.0" \
"r67081 o8 07 coe B14-o8 ${ref_freq_o8_1g} 0.0" \
"p67112 wb 20 coe B15-wb ${ref_freq_wb_1g} 0.0" \
"pre027 wb 18 coe B16-wb ${ref_freq_wb_1g} 0.0" \
"pre027 wb 22 coe B17-wb ${ref_freq_wb_1g} 0.0" \
"pre028 wb 11 coe B18-wb ${ref_freq_wb_1g} 0.0" \
"p67113 wb 24 coe B19-wb ${ref_freq_wb_1g} 0.0" \
"p67113 wb 29 coe B20-wb ${ref_freq_wb_1g} 0.0" \
"p67114 wb 17 coe B21-wb ${ref_freq_wb_1g} 0.0" \
"p67119 o8 03 coe B22-o8 ${ref_freq_o8_2g} 0.0" \
"p67116 wb 32 coe B23-wb ${ref_freq_wb_1g} 0.0" \
"p67117 wb 39 coe B24-wb ${ref_freq_wb_1g} 0.0" \
"p67118 wb 34 coe B25-wb ${ref_freq_wb_1g} 0.0" \
"p67118 o8 34 coe B25-o8 ${ref_freq_o8_4g} 0.0" \
"p67118 wb 46 coe B26-wb ${ref_freq_wb_1g} 0.0" \
"p67118 o8 46 coe B26-o8 ${ref_freq_o8_4g} 0.0" \
"p67118 wb 47 coe B27-wb ${ref_freq_wb_1g} 0.0" \
"p67118 o8 47 coe B27-o8 ${ref_freq_o8_4g} 0.0" \
"p67120 wb 18 coe B28-wb ${ref_freq_wb_1g} 0.0" \
"p67120 wb 27 coe B29-wb ${ref_freq_wb_1g} 0.0" \
"p67121 wb 39 coe B30-wb ${ref_freq_wb_1g} 0.0" \
"p67124 o8 35 coe B31-o8 ${ref_freq_o8_4g} 0.0" \
"p67124 o8 39 coe B32-o8 ${ref_freq_o8_4g} 0.0" \
"p67125 wb 18 coe B33-wb ${ref_freq_wb_1g} 0.0" \
"p67125 wb 22 coe B34-wb ${ref_freq_wb_1g} 0.0" \
"pre033 wb 32 coe B35-wb ${ref_freq_wb_1g} 0.0" \
"pre033 wb 38 coe B36-wb ${ref_freq_wb_1g} 0.0" \
"p67128 wb 21 coe B37-wb ${ref_freq_wb_1g} 0.0" \
"p67128 wb 24 coe B38-wb ${ref_freq_wb_1g} 0.0" \
"p67128 wb 40 coe B39-wb ${ref_freq_wb_1g} 0.0" \
"p67138 wb 25 coe B41-wb ${ref_freq_wb_1g} 0.0" \
"p67138 wb 40 coe B42-wb ${ref_freq_wb_1g} 0.0" \
"r67l29 tr 19 coe B42-tr ${ref_freq_tr_1g} 0.0" \
"stocke st 01 st B01-st ${ref_freq_st_1g} 410.775" \
"stocke st 02 st B06-st ${ref_freq_st_1g} 410.775" \
"stocke st 03 st B10-st ${ref_freq_st_1g} 410.775" \
"stocke st 04 st B31-st ${ref_freq_st_1g} 410.775" \
"stocke st 05 st B32-st ${ref_freq_st_1g} 410.775" \
"stocke st 06 st B33-st ${ref_freq_st_1g} 410.775" \
"stocke st 07 st B34-st ${ref_freq_st_1g} 410.775" \
"stocke st 08 st B35-st ${ref_freq_st_1g} 410.775" \
"stocke st 09 st B37-st ${ref_freq_st_1g} 410.775" \
"stocke st 10 st B38-st ${ref_freq_st_1g} 410.775" \
"stocke st 11 st B39-st ${ref_freq_st_1g} 410.775" \
"stocke st 12 st B40-st ${ref_freq_st_1g} 410.775" \
"stocke st 13 st B43-st ${ref_freq_st_1g} 410.775" \
"stocke st 14 st B44-st ${ref_freq_st_1g} 410.775" \
"stocke st 15 st B45-st ${ref_freq_st_1g} 410.775" \
"stocke st 16 st B46-st ${ref_freq_st_1g} 410.775" \
)

basedir=<root-dir-where-the-data-live>
dbdir=./dbs/
#Fscrunch=4
#Tscrunch=2

db="${dbdir}/burst_info.pickle"  # everything goes here

get_MJD=./python/get_MJD.py

for run in "${runs[@]}";do
    IFS=" " read -r -a info <<< "${run}"
    exp="${info[0]}"
    dish="${info[1]}"
    scan="${info[2]}"
    loc="${info[3]}"
    id="${info[4]}"
    refFreq="${info[5]}"
    dm="${info[6]}"
    src='sfxc'
    src_dir='sfxc'
    freqref="--ref_freq ${refFreq}"
    if [[ ${refFreq} == 'none' ]];then
        freqref=''
    fi

    arbase=${exp}_${dish}_no00${scan}
    if [[ "$arbase" == 'r67014_o8_no0013' ]]; then
        # for that scan there are not baseband data, only old filterbanks.
        # Thus we take the 'scale' version of all
        src='scale'
        src_dir='scale_spc'
    fi
    expbasedir=${basedir}/${exp}
    arfile=${expbasedir}/${src_dir}/${arbase}.ds.pazi
    filfile=${expbasedir}/${src_dir}/${arbase}.ds.fil

    cmd="python3 ${get_MJD} -d ${dm} --db ${db} --src ${src} --arfile ${arfile} \
                            --ID ${id} ${freqref} --location ${loc} \
                            ${filfile}"
    echo "running ${cmd}"
    eval ${cmd}
done

#!/bin/bash
#
#########
#
# This script assumes that the required filterbank files live in ${basedir}/${experiment}
# - using 'dspsr', it will create archives files, attempting to have the pulse in the middle
# - each ar-file (one per subband) will first be scaled
# - the scaled ar-files with then be corrected using psrchive's SPC alogorithm
# - the scaled ar-files will be co-added to have the full band in one archive
# - the first scaled and then SPC'ed ar-files will be co-added to have the full band in one archive
# - a diff-archive of the scaled and the SPC'ed file will be generated
# - diagnostic plots of the scaled, SPC'ed, diff'ed arcives are generated
#
# - note that burst B04-o8 from r67014 scan 0013 is not taken care of here (no baseband data)
#
############
# the below works with dspsr: on the vdif-devel at commit 97ad7c68a08cac0b172e407a61c0e22f64267eb9
# and with psrchive on master branch at commit 18312bb996cdef8ea3e1f10f0160dba41c63784e
###########
pwait() {
    while [ $(jobs -p | wc -l) -ge $1 ]; do
	sleep 0.33
    done
}
njobs=16
expcounters=( 3 ) # this is zero based!

bw128=128000000
bw256=256000000
bw512=512000000
subbw16=16000000
subbw32=32000000

runs=( \
# exp dish scan bw     subbandbw  cepoch    scaling factor for spc
"r67008 o8 17 ${bw128} ${subbw16} 0.0000135 0.50" \
"r67011 o8 17 ${bw128} ${subbw16} 0.0000230 0.50" \
"r67014 o8 46 ${bw128} ${subbw16} 0.0000225 0.50" \
"r67014 o8 13 ${bw128} ${subbw16} 0.0000225 0.50" \
"r67022 o8 27 ${bw128} ${subbw16} 0.0000155 0.50" \
"r67022 o8 28 ${bw128} ${subbw16} 0.0000210 0.50" \
"r67028 o8 20 ${bw128} ${subbw16} 0.0000235 0.50" \
"r67l01 tr 29 ${bw256} ${subbw32} 0.0000225 0.50" \
"r67030 o8 16 ${bw128} ${subbw16} 0.0000150 0.50" \
"r67079 o8 27 ${bw128} ${subbw16} 0.0000230 0.50" \
"r67079 o8 31 ${bw128} ${subbw16} 0.0000125 0.50" \
"r67081 o8 05 ${bw128} ${subbw16} 0.0000215 0.50" \
"r67081 o8 07 ${bw128} ${subbw16} 0.0000165 0.50" \
"p67112 wb 20 ${bw128} ${subbw16} 0.0000260 0.50" \
"pre027 wb 18 ${bw128} ${subbw16} 0.0000080 0.50" \
"pre027 wb 22 ${bw128} ${subbw16} 0.0000215 0.50" \
"pre028 wb 11 ${bw128} ${subbw16} 0.0000090 0.50" \
"p67113 wb 24 ${bw128} ${subbw16} 0.0000095 0.50" \
"p67113 wb 29 ${bw128} ${subbw16} 0.0000080 0.50" \
"p67114 wb 17 ${bw128} ${subbw16} 0.0000225 0.50" \
"p67116 wb 32 ${bw128} ${subbw16} 0.0000180 0.50" \
"p67117 wb 39 ${bw128} ${subbw16} 0.0000200 0.50" \
"p67118 wb 34 ${bw128} ${subbw16} 0.0000195 0.50" \
"p67118 o8 34 ${bw512} ${subbw32} 0.0000170 0.65" \
"p67118 wb 46 ${bw128} ${subbw16} 0.0000075 0.50" \
"p67118 o8 46 ${bw512} ${subbw32} 0.0000060 0.60" \
"p67118 wb 47 ${bw128} ${subbw16} 0.0000205 0.50" \
"p67118 o8 47 ${bw512} ${subbw32} 0.0000185 0.60" \
"p67119 o8 03 ${bw256} ${subbw16} 0.0000200 0.50" \
"p67120 wb 18 ${bw128} ${subbw16} 0.0000155 0.50" \
"p67120 wb 27 ${bw128} ${subbw16} 0.0000225 0.50" \
"p67121 wb 39 ${bw128} ${subbw16} 0.0000100 0.50" \
"p67124 o8 35 ${bw512} ${subbw32} 0.0000175 0.60" \
"p67124 o8 39 ${bw512} ${subbw32} 0.0000125 0.61" \
"p67125 wb 18 ${bw128} ${subbw16} 0.0000175 0.50" \
"p67125 wb 22 ${bw128} ${subbw16} 0.0000215 0.50" \
"pre033 wb 32 ${bw128} ${subbw16} 0.0000195 0.50" \
"pre033 wb 38 ${bw128} ${subbw16} 0.0000105 0.50" \
"p67128 wb 21 ${bw128} ${subbw16} 0.0000085 0.50" \
"p67128 wb 24 ${bw128} ${subbw16} 0.0000125 0.50" \
"p67128 wb 40 ${bw128} ${subbw16} 0.0000075 0.50" \
"p67138 wb 25 ${bw128} ${subbw16} 0.0000050 0.50" \
"p67138 wb 40 ${bw128} ${subbw16} 0.0000220 0.50" \
"r67l29 tr 19 ${bw128} ${subbw16} 0.0000210 0.50" \
)



dm=410.775

#freq_resolution=31250.0  #Hz; equiv of 16e6/512; allows for 32e-6s time bins for 16MHz subbands
freq_resolution=62500.0   # however, we have to go to 64e-6s time bins for spc to work; the DM-smearing
                          # in the lowest channel is more than 32e-6s.
#freq_resolution=31250.0  # spc doesn't work on all subbands with this freq setting
#time_res=0.000064
time_res=0.000128
#time_res=0.000256 # For B04-o8
length=2.097152

nbins=`echo "${length}/${time_res}" | bc -l | cut -d '.' -f1`

basedir=/data1/franz/
dir_in_outdir=scale_spc/

dspsr=`which dspsr`
pam=`which pam`

# we create archives on a per subband level, then scale them, spc them, psradd them.
# need to have NPOL 2 in hdr files for vdif-devel dspsr to work.

for expcounter in "${expcounters[@]}";do
    run="${runs[${expcounter}]}"
#for run in "${runs[@]}";do
    IFS=" " read -r -a info <<< "${run}"
    exp="${info[0]}"
    dish="${info[1]}"
    scan="${info[2]}"
    bw="${info[3]}"
    subbandbw="${info[4]}"
    cepoch="${info[5]}"
    scale="${info[6]}"
    nif=`echo "${bw}/${subbandbw}" | bc | cut -d '.' -f1`
    echo "${exp} ${dish} ${scan} ${bw} ${subbandbw} ${cepoch} ${nif}"
    outdir=${basedir}/${exp}/${dir_in_outdir}
    if ! [ -d ${outdir} ];then
	mkdir -p ${outdir}
    fi
    nchans=`echo "${subbandbw}/${freq_resolution}" | bc | cut -d '.' -f1`
    leakage_factor=`echo "${nchans}*4" | bc | cut -d '.' -f1`
    for IF in `seq 1 ${nif}`;do
	hdr=${exp}_${dish}_no00${scan}_IF${IF}.vdif_pol2.hdr
	ar=${outdir}/${hdr%.hdr}
        cmd=${cmd}"${dspsr} -b ${nbins} -D ${dm} -c ${length} -T ${length} -cepoch ${cepoch} ${ar}.fil -O ${ar}.fil && "
        cmd=${cmd}"./spc/scale_${scale}.sh -e scale.ar ${ar}.fil.ar && "
        cmd=${cmd}"./spc/spc.sh -e spc.ar ${ar}.fil.scale.ar && "
        cmd=${cmd}"${pam} ${fscrunch} -e ds.ar ${ar}.fil.scale.ar && "
        cmd=${cmd}"${pam} ${fscrunch} -e ds.ar ${ar}.fil.scale.spc.ar "
        echo "$cmd"
        eval ${cmd} &
        pwait ${njobs}
    done
    wait
    # add the subbands into one
    ar_scale=${exp}_${dish}_no00${scan}_allIFs.vdif_pol2.fil.scale.ar
    ar_spc=${exp}_${dish}_no00${scan}_allIFs.vdif_pol2.fil.scale.spc.ar
    ar_scale_ds=${exp}_${dish}_no00${scan}_allIFs.vdif_pol2.fil.scale.ds.ar
    ar_spc_ds=${exp}_${dish}_no00${scan}_allIFs.vdif_pol2.fil.scale.spc.ds.ar
    psradd -R -o ${outdir}/${ar_scale} ${outdir}/${exp}_${dish}_no00${scan}_IF*.vdif_pol2.fil.scale.ar
    psradd -R -o ${outdir}/${ar_scale_ds} ${outdir}/${exp}_${dish}_no00${scan}_IF*.vdif_pol2.fil.scale.ds.ar
    psradd -R -o ${outdir}/${ar_spc} ${outdir}/${exp}_${dish}_no00${scan}_IF*.vdif_pol2.fil.scale.spc.ar
    psradd -R -o ${outdir}/${ar_spc_ds} ${outdir}/${exp}_${dish}_no00${scan}_IF*.vdif_pol2.fil.scale.spc.ds.ar
    psrplot -pfreq+ -c x:unit=s -D /CPS -jpDT -j 'F x4' -j 'B x8' ${outdir}/${ar_scale}
    mv pgplot.ps ${outdir}/${ar_scale}.ps
    psrplot -pfreq+ -c x:unit=s -D /CPS -jpDT -j 'F x4' -j 'B x8' ${outdir}/${ar_spc}
    mv pgplot.ps ${outdir}/${ar_spc}.ps
    psrdiff ${outdir}/${ar_scale} ${outdir}/${ar_spc}
    mv psrdiff.out ${outdir}/${exp}_${dish}_no00${scan}_scale_spc_diff.ar
    psrplot -pD -c x:unit=s -D /CPS -jpFD ${outdir}/${exp}_${dish}_no00${scan}_scale_spc_diff.ar
    mv pgplot.ps ${outdir}/${exp}_${dish}_no00${scan}_scale_spc_diff.ar.ps
done

# Intro

This repo contains the scripts and databases that were used to generate the plots and results found in [Kirsten+2023](https://arxiv.org/abs/2306.15505). The data needed to run the scripts are available at [this Zenodo link](https://zenodo.org/uploads/10006350). The tarball contains the readily processed and flagged data, as well as the unprocessed filterbanks. Thus, it is possible to use the plotting scripts directly or reprocess the data first.

The database that contains all of the numbers generated via the scripts in this repo is in `./dbs/burst_info.[pickle,cvs]`. For convenience, the data from [Xu et al. 2022](https://www.nature.com/articles/s41586-022-05071-8) are also provided in `./dbs/FRB20201124A_2021AprMay_FAST_BurstInfo.txt` 

If you are using this code, please cite [Kirsten+2023](https://arxiv.org/abs/2306.15505)

# Run the code
Most scripts contain somewhat of an explanation right at the top of the code and further comments throught. In order to reprocess the data, the following scripts need to be run in order:

- r67_generate_archives.sh
  - uses the per-subband filterbanks to generate psrchive-archives,
  - creates readily scaled and also SPC'ed archives, diff archives and diagnostic plots
- flag_archives.sh
  - applies manually found RFI-flags to all archive-files (generated both from digifil-created filterbanks and SFXC-generated filterbanks)
- get_fluences.sh
  - this is a wrapper around `./python/ACF_archive.py` which helps to interactively determine the fluences per burst
- get_MJDs.sh
  - is a wrapper around `./python/get_MJD.py` which uses some of the output from the previous step to interactively determine the TOAs per burst
  


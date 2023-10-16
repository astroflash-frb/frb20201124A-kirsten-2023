#!/usr/bin/env python
# coding: utf-8
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from matplotlib import rcParams

rcParams['ytick.direction'] = 'in'
rcParams['xtick.direction'] = 'in'
rcParams['xtick.labelsize'] = 10
rcParams['ytick.labelsize'] = 10
rcParams['axes.labelsize'] = 10
rcParams['legend.fontsize'] = 8
rcParams['legend.handlelength'] = 1.0

def load_df(file):
    "Load in the dataframe"
    fluence_df = pd.read_pickle(file)
    print(fluence_df['src'].unique())
    ##Split the dataframe into one without Stockert for the comparission plot
    fluence_df_diff = fluence_df[fluence_df['experiment'] != 'stocke']
    return fluence_df_diff


def fluence(fluence_df, src):
    "return the array with fluence values"
    fluence_list = []
    exps = fluence_df['id'].unique()
    fluence_df_config_temp = fluence_df[(fluence_df['src'] == src)]

    for exp in exps:
        fluence_df_exp_temp = fluence_df_config_temp[fluence_df_config_temp['id'] == exp]
        fluence = sum(fluence_df_exp_temp['fluence_jyms'].values)
        fluence_list.append(fluence)

    fluence_array = np.asarray(fluence_list)
    return fluence_array


def sort_fluence(fluence_df_diff):
    #list of all experiments over which we can loop
    ids = fluence_df_diff['id'].unique()
    configs = fluence_df_diff['src'].unique()

    fluence_scale = fluence(fluence_df_diff, src='scale')
    fluence_spc = fluence(fluence_df_diff, src='spc')
    fluence_sfxc = fluence(fluence_df_diff, src='sfxc')

    #Order burst from large to small Fluence
    zip_ordered = sorted(zip(fluence_scale, fluence_spc, fluence_sfxc, ids), reverse=True)
    fluences_sorted = list(zip(*zip_ordered))

    fluences_scale = np.array(fluences_sorted[0])
    fluences_spc = np.array(fluences_sorted[1])
    fluences_sfxc = np.array(fluences_sorted[2])
    burst_id = np.array(fluences_sorted[3])
    print("check len", len(fluences_scale), len(fluences_spc), len(fluences_sfxc), len(burst_id))
    return fluences_scale, fluences_spc, fluences_sfxc, burst_id


def fluence_plot(fluences_scale, fluences_spc, fluences_sfxc, burst_id):
    x_grid = np.arange(1, len(fluences_scale)+1)
    plt.figure(figsize=(5, 4))
    plt.scatter(fluences_spc, fluences_scale/fluences_spc, label="digifil/SPC", marker="^", s=60)
    plt.scatter(fluences_spc, fluences_sfxc/fluences_spc, label="SFXC/SPC", marker="o", s=60)
    plt.axhline(1, color='red', linestyle='-.')
    plt.ylabel("Fluence ratio")
    plt.xlabel(r"SPC-corrected Fluence [Jy$\,$ms]")
    plt.legend()
    plt.savefig("../plots/fluence_ratio.png", bbox_inches='tight',
                facecolor='white', transparent=False, dpi=600)
    try:
        plt.savefig(f"{os.environ['HOME']}/git/overleaf/r67-single-dish/figures/fluence_ratio.png",
                    bbox_inches='tight',
                    facecolor='white', transparent=False, dpi=600)
    except:
        print('Could not save to overleaf. Got only a local copy.')
    #plt.show()

fluence_df_diff = load_df(file="../dbs/burst_info.pickle")
fluences_scale, fluences_spc, fluences_sfxc, burst_id = sort_fluence(fluence_df_diff)
fluence_plot(fluences_scale, fluences_spc, fluences_sfxc, burst_id)

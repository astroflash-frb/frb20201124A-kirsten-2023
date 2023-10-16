#!/usr/bin/env python
# coding: utf-8
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import pandas as pd
import os
from astropy.time import Time
from matplotlib.patches import Rectangle
from matplotlib.patches import Patch
from station_colors import get_station_colors


def burst_mjds_csv():
    "Function to retrieve the topocentric arrival times in UTC from the bursts, from the main csv file"

    #file
    csv_file = "../dbs/burst_info.csv"
    #Read in the main csv file
    mjd_df = pd.read_csv(csv_file, sep=",", header=0, index_col=0, na_values='NA', engine='python')

    #Remove Stockert
    #mjd_df = mjd_df[mjd_df['experiment'] != 'stocke']

    #Remove the nan entries
    mjd_df = mjd_df[~mjd_df['toa_mjd_utc'].isna()]

    mjd_list = []
    exps = mjd_df['id'].unique()

    for exp in exps:
        mjd_df_exp = mjd_df[mjd_df['id'] == exp]

        #Just get the average of the MJD, the difference between average and actual arrival time is not visible in the plot
        mjd = np.average(mjd_df_exp['toa_mjd_utc'].values)
        mjd_list.append(mjd)

    data_mjd = {'mjd': mjd_list, 'exp': exps}
    df_mjd_small = pd.DataFrame(data=data_mjd)
    df_mjd_small = df_mjd_small.sort_values(by=['mjd'])
    #print(df_mjd_small)
    return df_mjd_small


def burst_mjds(breakpoint):

    #Retrieve all the average MJD in an array
    df_mjd_small = burst_mjds_csv()
    #mjd_arr = df_mjd_small['mjd'].values
    #mjd_arr.sort()

    #Since we have burst detected at the same MJD, this will filter out those dupplicates.
    #mjd_arr is of lenght 47
    #syntax from: https://stackoverflow.com/questions/58202837/drop-element-in-numpy-array-or-pandas-series-if-difference-to-previous-element
    df_mjd_small= df_mjd_small[~(df_mjd_small['mjd'].diff()<=0.002)]
    mjd_arr = df_mjd_small['mjd'].values

    #Will create two arrays with MJD, based on breakpoint
    mjd_1 = mjd_arr[np.where(mjd_arr < breakpoint)]
    mjd_2 = mjd_arr[np.where(mjd_arr > breakpoint)]
    print(len(mjd_1), len(mjd_2))
    return mjd_1, mjd_2


def dummy_maker(telescope, obs_code, df_temp):
    #"Function to make an overview of begin, end and frequency setup per observations"
    #"duplicates are allowed in this df"

    dummy_data = {}
    dummy_data["Telescope"] = telescope
    dummy_data["obs_code"] = obs_code

    dummy_data["mjd_start"] = np.min(df_temp["MJD_start_time"])
    dummy_data["mjd_end"] = np.max(df_temp["MJD_end_time"])

    bandwidth = df_temp["Total_Bandwidth_(MHz)"].values[0]
    nchan = df_temp["Number_of_channels"].values[0]

    #Calculate half the channelwidth
    half_channel_width = (bandwidth/nchan) / 2

    #frequency setup of observation
    dummy_data["bandwidth"] = bandwidth
    #dummy_data["nchan"] = nchan

    dummy_data["cent_f"] = df_temp["Central_freq_(MHz)"].values[0]
    dummy_data["low_c_l"] = df_temp["Low_channel_(MHz)"].values[0] - half_channel_width
    dummy_data["high_c_h"] = df_temp["High_channel_(MHz)"].values[0] + half_channel_width

    rec_time_s = np.sum(df_temp["total_obs_time_s"])
    dummy_data["rec_time_s"] = rec_time_s
    dummy_data["rec_time_m"] = rec_time_s / 60.
    dummy_data["rec_time_h"] = rec_time_s / 3600.
    dummy_data["rec_time_d"] = rec_time_s / (3600. * 24.)
    return dummy_data

def small_r67(big_df):
    #"Function to load in the big r67 dataframe and reduce it to a smaller one"
    #"This smaller dataframe is based per observation instead of per scan"

    dict_list = []

    #list of all telescopes that were used
    telescopes = list(set(big_df["Telescope"]))
    for single_telescope in telescopes:
        big_df_telescope = big_df[(big_df["Telescope"] == single_telescope)]
        #list of all the observations code that were used for that specific telescope
        observations = list(set(big_df_telescope["scan_name"]))
        #observations.sort()
        for s_observation in observations:
            #Create a temp df for a single telescope for a single observation
            #Even though this df only had source r67, I filter on it as a failsave
            df_temp = big_df_telescope[(big_df_telescope["scan_name"] == s_observation)
                                       & (big_df_telescope["Source_Name"] == "R67")]

            # some runs have more than one frequency setup, need to account for this
            freqs = df_temp["Central_freq_(MHz)"].unique()

            #If there are more than 1 frequency setup in a single observations -- switching frequencies
            if len(df_temp["Central_freq_(MHz)"].unique()) > 1 :

                print(f"!!!Splitting the dataframe {s_observation} because of different frequency setup!!!")
                #When the frequency changes it will get a new numbers in the dataframe, this way we can utilze the groupby method of pandas to loop over the differt dfs
                #example: https://stackoverflow.com/questions/64163984/how-to-split-a-dataframe-each-time-a-string-value-changes-in-a-column
                df_temp['group'] = df_temp['Central_freq_(MHz)'].ne(df_temp['Central_freq_(MHz)'].shift()).cumsum()
                df_temp = df_temp.groupby('group')

                for name, df_temp_freqs in df_temp:

                    dummy_data = dummy_maker(single_telescope, s_observation, df_temp_freqs)
                    dict_list.append(dummy_data)

            elif s_observation == "p67116":
                #Only in the case of p17116, and also dirty hardcoded, because I don't want to break my head about it
                #In the case of this observations, we stopped it, to later resume it causing a gap in the schedule.
                #This happens at: 59611.648009 and before: 59611.896111 , so breakpoint at 59611.65
                #This is extremly ugly and hardcoded and I am sorry for those who have to look at it. it is what it is.

                dummy_data = dummy_maker(single_telescope, s_observation, df_temp[df_temp["MJD_end_time"] < 59611.65])
                dict_list.append(dummy_data)

                dummy_data = dummy_maker(single_telescope, s_observation, df_temp[df_temp["MJD_end_time"] > 59611.65])
                dict_list.append(dummy_data)
            else:

                dummy_data = dummy_maker(single_telescope, s_observation, df_temp)
                dict_list.append(dummy_data)

    df_r67_small = pd.DataFrame(dict_list)
    df_r67_small = df_r67_small.sort_values(by=['mjd_start'])
    df_r67_small.to_pickle('../dbs/df_r67_reduced.pkl')
    print(df_r67_small)

    return df_r67_small

def observing_plot(df_reduced, colours):
    # set width and height of figure in inches
    width, height = 18, 12
    fig, (ax1, ax2, ax3, ax4, ax5, ax6, ax7) = \
        plt.subplots(7, 1, figsize=(width, height), sharex=False,
                     gridspec_kw={'height_ratios': (1/1.86, 1, 1/5.3, 0.5, 1/1.86, 1, 1/5.3)})
    #ratio = maxbandwith + space / (min_bandwidth +space)
    #'height_ratios': (1/1.27, 1, 1/3.611, 0.5, 1/1.27, 1, 1/3.611)}
    fig.subplots_adjust(hspace=0)
    ax1.set_yticks([4600, 4800])
    ax5.set_yticks([4600, 4800])
    ax2.set_yticks([1300, 1450, 1600])
    ax6.set_yticks([1300, 1450, 1600])
    #ax3.set_yticks([300,350])
    #ax7.set_yticks([300,350])
    ax3.set_yticks([330])
    ax7.set_yticks([330])

    space = 20
    ax1.set_ylim(4550.0-space, 4806.0+space)
    ax2.set_ylim(1202-space,1714+space)
    ax3.set_ylim(290-space,364+space)
    ax5.set_ylim(4550.0-space, 4806.0+space)
    ax6.set_ylim(1202-space,1714+space)
    ax7.set_ylim(300-space,364+space)

    ax1.spines['bottom'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax2.spines['bottom'].set_visible(False)
    ax3.spines['top'].set_visible(False)
    ax5.spines['bottom'].set_visible(False)
    ax6.spines['top'].set_visible(False)
    ax6.spines['bottom'].set_visible(False)
    ax7.spines['top'].set_visible(False)
    ax4.set_visible(False)
    ax2.tick_params(axis='x', length=0)

    ax1.tick_params(labeltop=False)  # don't put tick labels at the top
    ax3.tick_params(labeltop=False)
    ax2.tick_params(labelbottom=False)
    ax5.tick_params(labeltop=False)  # don't put tick labels at the top
    ax6.tick_params(labelbottom=False)
    ax7.tick_params(labeltop=False)
    ax6.tick_params(axis='x', length=0)
    #ax7.tick_params(axis='x', length=0)
    ax3.xaxis.tick_bottom()
    ax1.xaxis.tick_top()
    ax5.xaxis.tick_top()
    ax7.xaxis.tick_bottom()

    fig.text(0.06, 0.5, 'Frequency (MHz)', va='center', rotation='vertical', fontsize=16)
    ax7.set_xlabel("Time (MJD)", fontsize=16)
    for ax in [ax1, ax2, ax3, ax4, ax5, ax6, ax7]:
        ax.tick_params(labelsize=14)

    # Draw the diagonal lines to show broken axes
    d = 0.  # proportion of vertical to horizontal extent of the slanted line
    kwargs = dict(marker=[(-1, -d), (1, d)], markersize=16,
              linestyle="none", color='k', mec='k', mew=1, clip_on=False)
    ax1.plot([0, 1], [0, 0], transform=ax1.transAxes, **kwargs)
    ax2.plot([1, 0], [1, 1], transform=ax2.transAxes, **kwargs)
    ax2.plot([0, 1], [0, 0], transform=ax2.transAxes, **kwargs)
    ax3.plot([1, 0], [1, 1], transform=ax3.transAxes, **kwargs)

    ax5.plot([0, 1], [0, 0], transform=ax5.transAxes, **kwargs)
    ax6.plot([1, 0], [1, 1], transform=ax6.transAxes, **kwargs)
    ax6.plot([0, 1], [0, 0], transform=ax6.transAxes, **kwargs)
    ax7.plot([1, 0], [1, 1], transform=ax7.transAxes, **kwargs)
    fig.subplots_adjust(hspace=0.1)
    legend_elements = [Patch(facecolor=colours['tr'], label=r'Toru$\mathrm{\'n}$'),
                       Patch(facecolor=colours['on'], label='Onsala'),
                       Patch(facecolor=colours['st'], label='Stockert'),
                       Patch(facecolor=colours['wb'], label='Westerbork')]
    ax1.legend(handles=legend_elements, loc='upper right', fontsize=13)
    for i in tqdm(range(len(df_reduced["mjd_start"]))):
        tele = df_reduced["Telescope"][i]
        obs_code = df_reduced["obs_code"][i]
        mjd_s = df_reduced["mjd_start"][i]
        mjd_e = df_reduced["mjd_end"][i]
        low_c = df_reduced["low_c_l"][i]
        high_c = df_reduced["high_c_h"][i]
        bandwidth = df_reduced["bandwidth"][i]
        obs_time_d = df_reduced["rec_time_d"][i]
        if tele == "tr":
            lab, cc = "Torun", colours['tr']
        if tele == "stk":
            lab, cc = "Stockert", colours['st']
        if tele == "o8":
            lab, cc = "Onsala", colours['on']
        if tele == "wb":
            lab, cc = "Westerbork", colours['wb']
        #Draw the figure in all three axes in one panel
        for ax in (ax1, ax2, ax3, ax5, ax6, ax7):
            rect = Rectangle((mjd_s, low_c), obs_time_d, height=bandwidth,
                             facecolor=cc, lw=0.0, edgecolor='k')
            ax.add_patch(rect)
    ax3.get_shared_x_axes().join(ax1, ax2, ax3)
    ax7.get_shared_x_axes().join(ax5, ax6, ax7)

    date_right_bound = 59600
    mjd_1, mjd_2 = burst_mjds(date_right_bound)
    #list_mjd_panel1 = [59326.6459, 59329.5213, 59337.80098, 59349.57226, 59349.57695, 59356.47362, \
    #            59358.4335]
    #list_mjd_panel2 = [59485.2013, 59485.2363, 59485.9682]

    for burst_mjd in mjd_1:
        ax1.axvline(x=burst_mjd, ymax=1, ymin=-1.5, clip_on=False, ls='--', color='k')
        ax2.axvline(x=burst_mjd, ymax=1, ymin=0, clip_on=False, ls='--', color='k')
        ax3.axvline(x=burst_mjd, ymax=1.5, ymin=0, clip_on=False, ls='--', color='k')

    for burst_mjd in mjd_2:
        ax5.axvline(x=burst_mjd, ymax=1, ymin=-1.5, clip_on=False, ls='--', color='k')
        ax6.axvline(x=burst_mjd, ymax=1, ymin=0, clip_on=False, ls='--', color='k')
        ax7.axvline(x=burst_mjd, ymax=1.5, ymin=0, clip_on=False, ls='--', color='k')

    #Set the limits on the plots, -10 because the plot will get some overla
    ax1.set_xlim(59300,date_right_bound)
    ax5.set_xlim((date_right_bound),59643)
    #Set dates to be first of month
    ax1.set_xticks([59300, 59400, 59500, 59600])
    ax3.set_xticks([59300, 59400, 59500, 59600])

    times_top = ax1.get_xticks()
    t_t = Time(times_top, format='mjd', scale='utc')
    mjd_isot_t = t_t.to_value('isot', 'date')
    ax1.set_xticklabels(mjd_isot_t)
    ax1.tick_params(labelrotation=0)

    #ax5 is for the date isot ticks on top
    #ax7 is for mjd at the bottom
    ax5.set_xticks([59600,59620, 59640])
    ax7.set_xticks([59600,59620, 59640])

    times_bottom = ax5.get_xticks()
    t_b = Time(times_bottom, format='mjd', scale='utc')
    mjd_isot_b = t_b.to_value('isot', 'date')
    ax5.set_xticklabels(mjd_isot_b)
    ax5.tick_params(labelrotation=0)

    # highlight time range that FAST observed
    yminvspan=0
    ymaxvspan=0.2
    #ax1.axvspan(xmin=59305, xmax=59363, ymin=yminvspan, ymax=ymaxvspan, alpha=0.15, color='gray')
    #ax2.axvspan(xmin=59305, xmax=59363, ymin=yminvspan, ymax=ymaxvspan, alpha=0.15, color='gray')
    ax3.axvspan(xmin=59305, xmax=59363, ymin=yminvspan, ymax=ymaxvspan, alpha=0.5, color='black')

    lbl0=ax1.text(0.03,0.75,(r'${\rm {\bf %s}}$')%'a', ha='center', va='center', fontsize=14, transform = ax1.transAxes,
                  bbox=dict(facecolor='grey', alpha=0.2, edgecolor='black'))
    lbl1=ax5.text(0.03,0.75,(r'${\rm {\bf %s}}$')%'b', ha='center', va='center', fontsize=14, transform = ax5.transAxes,
                  bbox=dict(facecolor='grey', alpha=0.2, edgecolor='black'))


    plt.savefig("../plots/Observingplot_V4_2panels.pdf", bbox_inches='tight')
    #plt.savefig("../Observingplot_V4_2panels" + ".png", bbox_inches='tight',
    #            facecolor='white', transparent=False, dpi=600)
    try:
        figpath = f"{os.environ['HOME']}/git/overleaf/r67-single-dish/figures/"
        plt.savefig(f"{figpath}/Observingplot_V4_2panels.pdf", bbox_inches='tight')
    except:
        print('Could not save to overleaf git. Got a local copy though.')
    #plt.show()


df_r67 = pd.read_pickle('../dbs/df_r67.pkl')
df_r67_small = small_r67(df_r67)
#df_r67_small[df_r67_small["obs_code"] == "r67008"]

mjds = df_r67_small["mjd_start"].values
print(min(mjds), max(mjds))
observing_plot(df_r67_small, colours=get_station_colors())

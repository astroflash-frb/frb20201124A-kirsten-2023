#!/usr/bin/env python
# coding: utf-8

import argparse
import numpy as np
import pandas as pd


def pickle_open():

    r67_df = pd.read_pickle('../dbs/df_r67_reduced.pkl')
    r67_df = r67_df.sort_values(by=['mjd_start'])

    r67_df = r67_df[r67_df['mjd_start'].between(59611, 59618)]
    print(r67_df)

    return

def dummy_maker(telescope, obs_code, df_temp):
    "Function to make an overview of begin, end and frequency setup per observations"
    "duplicates are allowed in this df"

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

    
#Based on code from Mark Snelders from the SGR paper
def small_r67():

    big_df = pd.read_pickle('../dbs/df_r67.pkl')

    "Function to load in the big r67 dataframe and reduce it to a smaller one"
    "This smaller dataframe is based per observation instead of per scan"
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
    
                dummy_data = dummy_maker(single_telescope, s_observation, df_temp[df_temp["MJD_end_time"] < 59611.65])
                dict_list.append(dummy_data)

                dummy_data = dummy_maker(single_telescope, s_observation, df_temp[df_temp["MJD_end_time"] > 59611.65])
                dict_list.append(dummy_data)
   
            else:

                dummy_data = dummy_maker(single_telescope, s_observation, df_temp)
                dict_list.append(dummy_data)
  
    df_r67_small = pd.DataFrame(dict_list)
    df_r67_small = df_r67_small.sort_values(by=['mjd_start'])
    df_r67_smaller = df_r67_small[df_r67_small['mjd_start'].between(59605, 59618)]
    #print(r67_df)
    print(df_r67_smaller)

    return

small_r67()


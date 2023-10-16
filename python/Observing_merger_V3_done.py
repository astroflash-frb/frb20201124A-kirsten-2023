#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import glob
import sys
#Need this command for the _merge function for recursion
sys.setrecursionlimit(10**5)


def pandas_grab(path):

    #extension to find all pickle files in directory
    glob_path = path + "/**/*.pkl"
    pandas_files = glob.glob(glob_path, recursive=True)

    df_all = pd.concat([pd.read_pickle(fp) for fp in pandas_files], ignore_index=True)

    return df_all

df_all = pandas_grab(path = "../scripts_omar/observing_hours/r67/")

#dictionary with frequency setups
def freq_dict_maker():

    #acces with: freq_dict[332]
    freq_dict = {
        332 : "P",
        1424 : "L1",
        1524 : "L1_offset",
        1446 : "L8",
        1323.49 : "L10",
        1468 : "L11",
        1418 : "L12",
        1381.46484375 : "L13",
        1414 : "L15",
        1488 : "L17_2g",
        1458 : "L16_4g",
        4678 : "C4",
        4664 : "C5"
    }

    return freq_dict


def _merge(head, tail):
    if tail == []:
        return head

    a, b = head[-1]
    x, y = tail[0]

    do_merge = b > x
    if do_merge:
        head_ = head[:-1] + [(a, max(b, y))]
        tail_ = tail[1:]
        return _merge(head_, tail_)
    else:
        head_ = head + tail[:1]
        tail_ = tail[1:]
        return _merge(head_, tail_)


def merge_intervals(lst):
    if len(lst) <= 1:
        return lst

    #sort on the first element of the tuples for the whole list
    lst = sorted(lst, key=lambda x: x[0])

    return _merge(lst[:1], lst[1:])


def get_overlap_times(df, coln1="MJD_start", coln2="MJD_end"):

    # get the right columns
    start = np.array(df["MJD_start_time"])
    end = np.array(df["MJD_end_time"])

    # check if input is correct:
    for s, e in zip(start, end):
        if s >= e:
            print("WARNING: start time should not exceed end time")

    # create tuple pairs:
    times = [(x, y) for x, y in zip(start, end)]
    times_reduced = merge_intervals(times)

    return times_reduced


def calc_overlap_times(times, unit="days"):
    """ Times is a list of tuple pairs, e.g.:
    [(1, 3), (5.5, 6.0)]
    Will return the sum of the differences of the pairs:
    [(1, 3), (5.5, 6.0)] --> (3-1) + (6.0-5.5) = 2 + 0.5 = 2.5 """
    ttotal = 0 # total time
    for timepair in times:
        s, e = timepair[0], timepair[1]
        ttotal += (e - s)
    print(f"The total non overlapping observing time is {ttotal:.2f} {unit} or {(ttotal*24):.2f} hours")


def pickle_opener(df_all, total_time=False,
                  activity_phase=False, mjdstart=59305, mjdstop=59363):

    #call the freq dictionary
    freq_dict = freq_dict_maker()

    #Conmbine all the pickles files into one big dataframe
    df_all["total_obs_time_s"] = (df_all["MJD_end_time"] - df_all["MJD_start_time"]) * 24 * 3600

    #Create mask for only R67, apply mask.
    mask_r67 = df_all['Source_Name'] == 'R67'
    df_r67 = df_all[mask_r67]

    ##If you want to limit in time:
    if activity_phase is True:
        #59300 start of observing and 59363 1 day after Stk
        df_r67 = df_r67[(df_r67['MJD_start_time'] >= mjdstart) &
                        (df_r67['MJD_end_time'] < mjdstop)]

        #59602 start of wb; 59644 last Stk obs
        #df_r67 = df_r67[(df_r67['MJD_start_time'] > 59602) & (df_r67['MJD_end_time'] < 59644)]

    #print setup of observing frequenies
    freqs = df_r67['Central_freq_(MHz)'].values
    print(set(freqs))

    #Westerbork
    df_wb = df_r67[df_r67["Telescope"] == 'wb']
    wb_obs_h = sum(df_wb["total_obs_time_s"].values)/3600

    #Onsala
    df_o8 = df_r67[df_r67["Telescope"] == 'o8']
    o8_obs_h = sum(df_o8["total_obs_time_s"].values)/3600

    #Torun
    df_tr = df_r67[df_r67["Telescope"] == 'tr']
    tr_obs_h = sum(df_tr["total_obs_time_s"].values)/3600

    #Stockert
    df_stk = df_r67[df_r67["Telescope"] == 'stk']
    stk_obs_h = sum(df_stk["total_obs_time_s"].values)/3600

    #Print observing hours
    print(f"total obs Wb {wb_obs_h:.2f} h")
    print(f"total obs o8 {o8_obs_h:.2f} h")
    print(f"total obs Tr {tr_obs_h:.2f} h")
    print(f"total obs Stk {stk_obs_h:.2f} h")

    total_obs = sum(df_r67["total_obs_time_s"].values) / 3600
    print(f"total observing time: {total_obs:.2f} h")
    print("--------------------")

    print("-P band-")
    #Pband
    df_r67_pband = df_r67[df_r67["Central_freq_(MHz)"] == 332.0]
    times_p = get_overlap_times(df_r67_pband)
    calc_overlap_times(times_p)
    print(f'Total P-band: {sum(df_r67_pband["total_obs_time_s"].values)/3600:.2f} hrs.')

    print("-C band-")
    #Cband
    df_r67_cband = df_r67[(df_r67["Central_freq_(MHz)"] == 4664.0) | (df_r67["Central_freq_(MHz)"] == 4678.0)]
    times_c = get_overlap_times(df_r67_cband)
    calc_overlap_times(times_c)
    print(f'total C-band: {sum(df_r67_cband["total_obs_time_s"].values)/3600:.2f} hrs.')

    print("-L band-")
    #Lband
    df_r67_lband = df_r67[(df_r67["Central_freq_(MHz)"] != 332.0) &
                          (df_r67["Central_freq_(MHz)"] != 4664.0) &
                          (df_r67["Central_freq_(MHz)"] != 4678.0)]
    times_l = get_overlap_times(df_r67_lband)
    calc_overlap_times(times_l)
    print(f'Total L-band: {sum(df_r67_lband["total_obs_time_s"].values)/3600:.2f} hrs.')

    #Calculate the unique hours in the observations
    if total_time is True:
        print("--------------------")
        times = get_overlap_times(df_r67)
        calc_overlap_times(times)
        print("--------------------")

    #Counters for observing time
    total_h_obs = 0
    total_h_obs_lband = 0

    #Loop over all entries in the freq_dict
    for obser_freq in freq_dict:

        #Create mask per frequency setup
        df_freq_setup_mask = (df_r67['Central_freq_(MHz)'] == obser_freq)
        df_freq_setup = df_r67[(df_freq_setup_mask)]

        #As a double check; get the set of all telescopes in this setup (should be unique)
        part_telescope = set(df_freq_setup["Telescope"])

        #Amount of hours per frequency setup in hours
        df_freq_setup_h = sum(df_freq_setup["total_obs_time_s"].values)/3600
        print(f"{freq_dict[obser_freq]}: time = {df_freq_setup_h:.2f} h -- central_freq {obser_freq} -- {part_telescope}")

        #Add to total
        total_h_obs += df_freq_setup_h

        #Add if observation is on Lband
        if freq_dict[obser_freq].startswith("L"):
            total_h_obs_lband += df_freq_setup_h

    print(f"Total observing on L-band: {total_h_obs_lband:.2f}")
    print(f"as a double check, the total is over the freq setup is: {total_h_obs:.2f}")
    print("--------------------")

    df_r67 = df_r67.sort_values(by=['MJD_start_time'])

    return df_r67

# FAST window
mjdstart = 59305
mjdstop = 59363

## 2nd activity window
mjdstart = 59602
mjdstop =59642

df_r67 = pickle_opener(df_all, total_time=True, activity_phase=True,
                       mjdstart=mjdstart, mjdstop=mjdstop)
#df_r67.to_pickle("/data1/omar/sfxc_r67/Analysis/obsering_hours/merged_files/" + "df_r67" + '.pkl')
#df_r67.to_pickle('../dbs/df_r67.pkl')

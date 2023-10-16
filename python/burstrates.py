#!/usr/bin/env python3
# coding: utf-8

# ## Cumulative burst distributions
#
# Version8 of the burst distributions <br>
#
# In this version I will try to implement the burst rate the activity is between: <br>
# full first activity window: 59264 and 59363 <br>
# FAST window : 59307 and 59363 <br>
# third activity window: 59600 and 59645 <br>
#
# Determine the alpha
# - (done) threshold: MLE code, power-law code -> determine the alpha.
# - look by eye
# - (done) Threshold by SNR in Fluence equation
# - Use the cutoff as input parameter to determine the actual alpha with curve-fit
#
# Error on the alpha
# - bootstrapping to calculate a better error on the alpha. take a subsample and using curve fit; add the alpha to a list and then std on the list of alpha's.
# - 'right plot' exclude from low-high; determine the alpha -> see Kelly paper

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib.patheffects as PathEffects
import scipy.optimize as spopt
import numpy.ma as ma
import pandas as pd
import powerlaw

from tqdm import tqdm
from matplotlib import transforms
import random
from scipy.optimize import curve_fit
from station_colors import get_station_colors

def fluence_limit(sefd=420, width=1, bw=100, snr=20):
    '''
    sefd in Jy, width in ms, bw in MHz
    '''
    fluence_limit = snr * sefd * np.sqrt( (width/1000.) / (2 * bw * 1e6))
    # return things in Jy~ms
    return fluence_limit * 1000


def fluence_looper(df_fluence, sumcomponents=True):
    "function to sum over the components returns the burst id and fluences in a pd df"

    #the burst names to loop over
    exps = df_fluence['id'].unique()

    #Create a list with the telescope names for easy filtering later on
    exps_tel = [i.split('-')[1] for i in exps]

    #loop over every burst to sum the fluence
    fluence_list = []
    for exp in exps:
        fluence_df_exp_temp = df_fluence[df_fluence['id'] == exp]
        fluence = fluence_df_exp_temp['fluence_jyms'].values
        if sumcomponents:
            fluence = sum(fluence)
        fluence_list.append(fluence)

    #Convert the fluence list to array
    if sumcomponents:
        fluence_array = np.asarray(fluence_list)
    else:
        fluence_array = [np.asarray(arr) for arr in fluence_list]

    #Convert the Fluences to dataframe
    new_data = {'id': exps, 'fluence': fluence_array, 'telescope': exps_tel}
    df_new_fluence = pd.DataFrame(data=new_data)

    return df_new_fluence


def load_df(file, sumcomponents=True):
    "Load in the dataframe and return the fluences"

    #load in the file
    fluence_df = pd.read_pickle(file)
    #display(fluence_df)

    ##Wb, On and Tr Fluences
    #Only take the spc fluences
    df_spc = fluence_df[fluence_df['src'] == 'spc']
    df_spc_short = fluence_looper(df_spc, sumcomponents=sumcomponents)

    ##Stk Fluences
    df_stk = fluence_df[fluence_df['dish']=='st']
    df_stk_short = fluence_looper(df_stk, sumcomponents=sumcomponents)

    #Combine the two
    df_combined = pd.concat([df_spc_short, df_stk_short], ignore_index=True)
    #display(df_combined)
    return df_combined


pd.set_option('display.max_columns', 500)
pd.set_option('display.max_rows', 500)
fluence_df = load_df("../dbs/burst_info.pickle", sumcomponents=True)

def flatten_array(arr):
    try:
        arr = np.asarray([component for burst in arr for component in burst])
    except:
        pass
    return arr

def fluences_arrays(df_fluences):
    "return Fluences arrays from different activity windows"

    #### Onsala
    #Onsala window 1
    onsala_burst_window1 = ['B02-o8', 'B03-o8', 'B04-o8', 'B05-o8', 'B06-o8', 'B07-o8', 'B08-o8', 'B09-o8']
    onsala_burst_window2 = ['B11-o8', 'B12-o8', 'B13-o8', 'B14-o8']
    onsala_burst_window12 = np.concatenate([onsala_burst_window1, onsala_burst_window2])

    #Fluences window 1 for Onsala
    fluences_on_w1 = df_fluences[df_fluences['id'].isin(onsala_burst_window1)]
    #display(fluences_on_w1)
    fluences_on_w1_arr = flatten_array(fluences_on_w1['fluence'].values)
    #sort from high to low
    fluences_on_w1_arr[::-1].sort()


    #Fluences window 3 for Onsala
    fluences_on_w3 = df_fluences[~df_fluences['id'].isin(onsala_burst_window12)]
    fluences_on_w3 = fluences_on_w3[fluences_on_w3['telescope'] == 'o8']
    fluences_on_w3_arr = flatten_array(fluences_on_w3['fluence'].values)
    fluences_on_w3_arr[::-1].sort()
    #display(fluences_on_w3)

    #### Stockert
    stockert_burst_window1 = ['B01-st', 'B06-st', 'B10-st']

    #Fluences window 1 for Stk
    fluences_stk_w1 = df_fluences[df_fluences['id'].isin(stockert_burst_window1)]
    fluences_stk_w1_arr = flatten_array(fluences_stk_w1['fluence'].values)
    fluences_stk_w1_arr[::-1].sort()
    #display(fluences_stk_w1)

    #Fluences window 3 for Stk
    fluences_stk_w3 = df_fluences[~df_fluences['id'].isin(stockert_burst_window1)]
    fluences_stk_w3 = fluences_stk_w3[fluences_stk_w3['telescope'] == 'st']
    fluences_stk_w3_arr = flatten_array(fluences_stk_w3['fluence'].values)
    fluences_stk_w3_arr[::-1].sort()

    #display(fluences_stk_w3)

    #print(fluences_on_w1_arr, fluences_on_w3_arr)

    #### Westerbork
    fluences_wb_w3 = df_fluences[df_fluences['telescope'] == 'wb']
    #display(fluences_wb_w3)
    fluences_wb_w3_arr = flatten_array(fluences_wb_w3['fluence'].values)
    fluences_wb_w3_arr[::-1].sort()

    return fluences_on_w1_arr, fluences_on_w3_arr, fluences_stk_w1_arr, fluences_stk_w3_arr,\
        fluences_wb_w3_arr

fluences_on_w1_arr, fluences_on_w3_arr, fluences_stk_w1_arr, fluences_stk_w3_arr, fluences_wb_w3_arr = fluences_arrays(fluence_df)

def FAST_read_data():
    "function to make a pd from the FAST paper"
    #source: https://psr.pku.edu.cn/index.php/publications/frb20201124a/

    col = ['Barycentric arrival times','S/N','PeakFlux','PeakFlux_error','Fluence', 'Fluence_error', 'eq_width',
           'eq_width_err', 'bandwidth', 'DM', 'DM_rr', 'RM', 'RM_err+', 'RM_err-', 'RM_corr_io', 'RM_corr_io_err',
           'dol', 'dol_err', 'doc', 'doc_err', 'scin_bw', 'scin_bw_err', 'pl_index', 'pl_index_err', 'rm_index',
           'rm_index_err-', 'rm_index_err+']
    data = pd.read_csv('../dbs/FRB20201124A_2021AprMay_FAST_BurstInfo.txt', header=None, skiprows=13,
                       sep='\s+', engine='python', names=col)
    data.sort_values(by='Barycentric arrival times', axis=0, inplace=True)
    fluence=data['Fluence'].values
    times = data['Barycentric arrival times']
    return fluence, times


fast_fluences, fast_times = FAST_read_data()

def FAST_sum_components(fluences, times, threshold=0.1):
    '''
    Function to sum up fluences of bursts that occur within threshold seconds as we
    consider them as components of one burst.
    '''
    # wait time between bursts in seconds
    diff = np.diff(times) * 86400.
    components = diff < threshold
    print(f'Bursts affected by summing: {components.sum()}')
    bursts = []
    added = False
    for i, addnext in enumerate(components):
        if not added:
            bursts.append(fluences[i])
        if addnext:
            bursts[-1] += fluences[i+1]
            added = True
        else:
            added = False
    return np.array(bursts)

fast_fluences_summed = FAST_sum_components(fast_fluences, fast_times, threshold=0.1)
#fast_fluences_summed = fast_fluences

def isotropic_energies(fluences, distance=453):
    "Convert list of fluences to spectral densities"
    "Method based on fluence.py by K.Nimmo"

    #convert from jyms to jys
    fluence_jys = np.array(fluences) * 1e-3

    #mpc to cm
    distance_lum_cm = 3.086e24*distance

    #redshift correction, see evernote for details
    z = 0.098
    red_cor = 1 / (1 + z)**2

    energy_iso = fluence_jys * 4*np.pi*(distance_lum_cm**2) * 1e-23 * red_cor
    #print(energy_iso)

    return energy_iso


def func_powerlaw(x, scale=1, gamma=1):
    return x**gamma * scale


def slope_est(energies_array):
    "Function to do a first order estimate of the slope and error"
    "This method comes from Dante and an older paper, see his slack and paper he mentions"

    energies_array.sort()
    M=len(energies_array)
    si=np.array(energies_array)/energies_array[0]

    corrfact=(M-1)/M
    a_star=corrfact*1/( (-1/M) * np.sum(np.log(si)) )
    err= M*a_star / ( (M-1)*(M-2)**0.5 )

    return a_star,err


def power_law_fit(energies, observed_time, scale=1, limit=(1, 1),
                  limit_b=False, fit=True, ignore_turnover=False,
                  time_is_weights=False):
    "Function to return the fitted powerlaw to the datapoints"

    energies.sort()
    if fit:
        #Determine the turning point in the powerlaw with the Powerlaw module
        if not limit_b:
            pw_results = powerlaw.Fit(energies)
            print("Turning point:", pw_results.xmin)
        else:
            pw_results = powerlaw.Fit(energies,  xmin=(limit[0], limit[1]))
            print("No limit_b, Turning point:", pw_results.xmin)
        if ignore_turnover:
            pw_results.xmin = limit[0]

        #Values that will be included in the fit
        energies_fit = energies[energies>=pw_results.xmin]

        #Values that are excluded in the fit
        energies_ex = energies[energies<pw_results.xmin]

        if not time_is_weights:
            #The burst rate, devide by the observed_time to get a rate per hours
            #ypoints for the fitted data
            y_points_fit = np.arange(1, len(energies_fit)+1, 1) / observed_time

            #ypoints for the not fitted data; scaled
            y_points_ex = np.arange(len(energies_fit)+1,
                                    len(energies_fit)+len(energies_ex)+1, 1) / observed_time
        else:
            y_points_fit = observed_time
            y_points_ex = []
        ###
        #Curve_fit for the slope
        a_star, a_err = slope_est(energies)
        #print(f"estimate: {a_star:.4f} {a_err:.4f}")

        #apply curve_fit to do the fit
        popt, pcov = curve_fit(func_powerlaw, energies_fit, np.flip(y_points_fit),
                               p0 = np.asarray([scale, a_star]), maxfev=50000)

        #flip x and y to use the x-error
        popt_flip, pcov_flip = curve_fit(func_powerlaw, np.flip(y_points_fit),energies_fit,
                                         p0=[popt[0]**(1./-popt[1]), 1./popt[1]],
                                         sigma=np.asarray(energies_fit*0.2),\
                                         absolute_sigma=True)

        #Calculate the popt values in the xy frame, use the errors on the fit
        popt_back=[popt_flip[0]**(1./-popt_flip[1]), 1./popt_flip[1]]
        sol2 = [popt_back, pcov_flip]

        sol2_err = np.sqrt(np.diag(sol2[1]))
        xx = np.linspace(np.min(energies_fit),np.max(energies_fit),500)
        p_y = func_powerlaw(xx, *sol2[0])

        ##Bootstrapping
        subset_int = int(len(energies_fit) * 0.9)
        print(f"Bootstrapping, only using {subset_int} points out of {len(energies_fit)} points")

        N=1000
        alpha_list = []
        for trial in tqdm(range(N)):

            #subset_list = random.sample(list(energies_fit), k = subset_int)

            #Get a random sample of 90% from the fitted energies
            ss_energie, ss_ypoint = zip(*random.sample(list(zip(energies_fit, y_points_fit)),
                                                       k=subset_int))

            #Sort them again from big to small
            ss_energie_1, ss_ypoint_1 = zip(*sorted(zip(ss_energie, ss_ypoint)))
            popt_bs_yx, pcov_bs_yx=curve_fit(func_powerlaw,np.flip(ss_ypoint_1),ss_energie_1,
                                             p0=[popt[0]**(1./-popt[1]),1./popt[1]],
                                             sigma=np.asarray(ss_energie_1)*0.2,absolute_sigma=True)

            #Covert back to values, bootstramp reversed
            popt_bs_r=[popt_bs_yx[0]**(1./-popt_bs_yx[1]),1./popt_bs_yx[1]]

            #Append the value for the slope the list
            alpha_list.append(popt_bs_r[1])

        #Calculate the boot_alpha
        boot_alpha = np.std(alpha_list)
        print(f"boot_alpha: {boot_alpha}")

        return energies_fit, y_points_fit, energies_ex, y_points_ex, \
            sol2, sol2_err, pw_results.xmin, boot_alpha
    else:
        if not time_is_weights:
            y_points = np.arange(1, len(energies)+1, 1) / observed_time
        else:
            y_points = observed_time
        return energies, y_points


def joint_powerlaw_fit(energies, y_points_fit, scale):
    # sorting things by energy
    tmp = np.sort(np.vstack([energies, y_points_fit]), axis=1)
    energies, y_points_fit = tmp
    ###
    #Curve_fit for the slope
    a_star, a_err = slope_est(energies)
    #print(f"estimate: {a_star:.4f} {a_err:.4f}")

    #apply curve_fit to do the fit
    popt, pcov = curve_fit(func_powerlaw, energies, np.flip(y_points_fit),
                           p0 = np.asarray([scale, a_star]), maxfev=50000)

    #flip x and y to use the x-error
    popt_flip, pcov_flip = curve_fit(func_powerlaw, np.flip(y_points_fit), energies,
                                     p0=[popt[0]**(1./-popt[1]), 1./popt[1]],
                                     sigma=np.asarray(energies*0.2),\
                                     absolute_sigma=True)

    #Calculate the popt values in the xy frame, use the errors on the fit
    popt_back=[popt_flip[0]**(1./-popt_flip[1]), 1./popt_flip[1]]
    sol2 = [popt_back, pcov_flip]

    sol2_err = np.sqrt(np.diag(sol2[1]))
    xx = np.linspace(np.min(energies), np.max(energies), 500)
    p_y = func_powerlaw(xx, *sol2[0])

    ##Bootstrapping
    subset_int = int(len(energies) * 0.9)
    print(f"Bootstrapping, only using {subset_int} points out of {len(energies)} points")

    N=1000
    alpha_list = []
    for trial in tqdm(range(N)):

        #subset_list = random.sample(list(energies_fit), k = subset_int)

        #Get a random sample of 90% from the fitted energies
        ss_energie, ss_ypoint = zip(*random.sample(list(zip(energies, y_points_fit)),                                                                k=subset_int))

        #Sort them again from big to small
        ss_energie_1, ss_ypoint_1 = zip(*sorted(zip(ss_energie, ss_ypoint)))

        popt_bs_yx, pcov_bs_yx = curve_fit(func_powerlaw,np.flip(ss_ypoint_1),ss_energie_1,                            p0=[popt[0]**(1./-popt[1]),1./popt[1]],                                sigma=np.asarray(ss_energie_1)*0.2,absolute_sigma=True)

        #Covert back to values, bootstramp reversed
        popt_bs_r=[popt_bs_yx[0]**(1./-popt_bs_yx[1]),1./popt_bs_yx[1]]

        #Append the value for the slope the list
        alpha_list.append(popt_bs_r[1])

    #Calculate the boot_alpha
    boot_alpha = np.std(alpha_list)
    print(f"boot_alpha: {boot_alpha}")

    return energies, y_points_fit, sol2, sol2_err, boot_alpha



def plot_window1_fast():
    "Cum disstribution of burst densities of window 1 for Onsala and Stk"

    SNR_limit = 15
    clrs = get_station_colors()

    energies_array_on = isotropic_energies(fluences_on_w1_arr)
    energies_array_stk = isotropic_energies(fluences_stk_w1_arr)
    energies_array_fast = isotropic_energies(fast_fluences)
    energies_array_fast_summed = isotropic_energies(fast_fluences_summed)

    #total amount of hours observed between 59300 and 59363
    #last mjd fast is 59360 and stk 59362.4

    f_limit_on = fluence_limit(sefd=310, width=3, bw=100, snr=SNR_limit)
    f_limit_st = fluence_limit(sefd=1100, width=3, bw=80, snr=SNR_limit)
    f_limit_fast = 0.053 # in Jy~ms
    energie_threshold_on = isotropic_energies([f_limit_on])
    energie_threshold_st = isotropic_energies([f_limit_st])
    energie_threshold_fast = isotropic_energies([f_limit_fast])

    #Onsala
    ob_time_on_l = 257.23  # 269.71
    #Onsala
    print('Onsala')
    energies_fit_on, y_points_fit_on, energies_ex_on, y_points_ex_on, \
        sol2_on, sol2_err_on, turning_point_on, b_alpha_on = \
            power_law_fit(energies_array_on, observed_time=ob_time_on_l,
                          scale=10**18, limit=(energie_threshold_on, 1e33),
                          limit_b=True, ignore_turnover=True)

    #Stockert
    ob_time_stk_l = 390.4 # 457.82 # 368.6
    print('Stockert')
    energies_stk, y_points_stk = \
        power_law_fit(energies_array_stk, observed_time=ob_time_stk_l,
                      fit=False)

    #Onsala + Stockert fit
    energies_array_on_st = np.concatenate([energies_fit_on, energies_array_stk])
    y_points_on_stk = np.concatenate([y_points_fit_on, y_points_stk])
    print('Joint O8+St')
    energies_fit_on_st, y_points_fit_on_st, \
        sol2_on_st, sol2_err_on_st, b_alpha_on_st = \
            joint_powerlaw_fit(energies_array_on_st, y_points_on_stk,
                               scale=10**18)


    #Fast
    ob_time_fast = 88.0 # 96.9 # 88
    #energies_fast, y_points_fast = power_law_fit(energies_array_fast, observed_time=ob_time_fast, fit=False)
    print('Doing FAST with Xu et al. fit limit')
    energies_fit_fast, y_points_fit_fast, energies_ex_fast, y_points_ex_fast,\
        sol2_fast, sol2_err_fast, turning_point_fast, b_alpha_f = \
            power_law_fit(energies_array_fast, observed_time=ob_time_fast, scale=10**18,
                          limit=(5.9e29, 1e33), limit_b=True, ignore_turnover=True)
    print('Doing FAST fitting for a fit limit')
    energies_fit_fast_0, y_points_fit_fast_0, energies_ex_fast_0, y_points_ex_fast_0,\
        sol2_fast_0, sol2_err_fast_0, turning_point_fast_0, b_alpha_f_0 = \
            power_law_fit(energies_array_fast, observed_time=ob_time_fast, scale=10**18)

    print('Doing FAST with Xu et al. fit limit --summed fluences')
    energies_fit_fast_summed, y_points_fit_fast_summed, energies_ex_fast_summed, y_points_ex_fast_summed,\
        sol2_fast_summed, sol2_err_fast_summed, turning_point_fast_summed, b_alpha_f_summed = \
            power_law_fit(energies_array_fast_summed, observed_time=ob_time_fast, scale=10**18,
                          limit=(5.9e29, 1e33), limit_b=True, ignore_turnover=True)
    print('Doing FAST fitting for a fit limit -- summed fluences')
    energies_fit_fast_0_summed, y_points_fit_fast_0_summed, energies_ex_fast_0_summed, y_points_ex_fast_0_summed,\
        sol2_fast_0_summed, sol2_err_fast_0_summed, turning_point_fast_0_summed, b_alpha_f_0_summed = \
            power_law_fit(energies_array_fast_summed, observed_time=ob_time_fast, scale=10**18)

    #Figure
    #plt.figure(figsize=(12,8))
    rcParams['ytick.direction'] = 'in'
    rcParams['xtick.direction'] = 'in'
    rcParams['xtick.labelsize'] = 12
    rcParams['ytick.labelsize'] = 12
    rcParams['axes.labelsize'] = 12
    rcParams['legend.fontsize'] = 10
    rcParams['legend.handlelength'] = 1.0

    fig = plt.figure(figsize=(18, 6))
                     #label_f)
    spec = fig.add_gridspec(ncols=2, nrows=1, wspace=0.1, hspace=4.0)

    ax0 = fig.add_subplot(spec[0, 0])

    #Onsala data
    #plt.plot(energies_fit_on,func_powerlaw(energies_fit_on,*sol2_on[0]),color='blue',lw=0.8)
    ax0.scatter(np.sort(energies_ex_on), np.flip(y_points_ex_on),
                color=clrs['on'], marker="s", alpha=0.5)
    ax0.scatter(np.sort(energies_fit_on), np.flip(y_points_fit_on),
                color=clrs['on'], marker = "s",
                label = fr"O8 2021")
                #label = fr"Onsala $\gamma$ = {sol2_on[0][1]:.2f} $\pm$ {sol2_err_on[1]:.2f} $\pm$ {b_alpha_on:.2f}")


    #Stockert data
    ax0.scatter(np.sort(energies_array_stk), np.flip(y_points_stk), color=clrs['st'], marker = "s",
                label = fr"St 2021")

    # FAST data -- counting each component as burst
    ##ax0.axvline(complete_fast_energy)
    #ax0.scatter(np.sort(energies_ex_fast), np.flip(y_points_ex_fast),
    #            color='red', marker="2", s=16, alpha=0.5)
    ##ax0.scatter(np.sort(energies_fit_fast), np.flip(y_points_fit_fast), color="#e66101", marker = "2",s=16,
    ##            label= fr"Fast")
    #ax0.scatter(np.sort(energies_fit_fast), np.flip(y_points_fit_fast), color='red',
    #            marker="2", s=16,
    #            label = 'FAST 2021 -- per component')

    #ax0.axvline(complete_fast_energy)
    ## FAST data -- summing bursts within 100ms and count them as one burst
    ax0.scatter(np.sort(energies_ex_fast_summed), np.flip(y_points_ex_fast_summed),
                color=clrs['fast'], marker="2", s=16, alpha=0.5)
    #ax0.scatter(np.sort(energies_fit_fast), np.flip(y_points_fit_fast), color="#e66101", marker = "2",s=16,
    #            label= fr"Fast")
    ax0.scatter(np.sort(energies_fit_fast_summed), np.flip(y_points_fit_fast_summed), color=clrs['fast'],
                marker="2", s=16,
                label = 'FAST 2021')

    #Onsala + Stocker fit
    ax0.plot(energies_fit_on_st, func_powerlaw(energies_fit_on_st, *sol2_on_st[0]),
             color='gray', lw=1.8,
             label = fr"O8+St $\gamma$ = {sol2_on_st[0][1]:.2f} $\pm$ {sol2_err_on_st[1]:.2f} $\pm$ {b_alpha_on_st:.2f}")


    ##Fast fits -- per component
    #fast_fit_plot_range = np.concatenate([energies_ex_fast[-200:], energies_fit_fast])
    ## FAST fit with fit limit determinded using powerlaw.Fit
    #ax0.plot(fast_fit_plot_range, func_powerlaw(fast_fit_plot_range, *sol2_fast_0[0]),
    #         color='red', lw=1.8, linestyle='-',
    #         label = r"$\mathrm{\gamma_{FAST}}$"+f" = {sol2_fast_0[0][1]:.3f} $\pm$ {sol2_err_fast_0[1]:.3f} $\pm$ {b_alpha_f_0:.3f}")
    ## FAST fit with break point from Xu et al. 2022 used as fitting limit
    #ax0.plot(fast_fit_plot_range, func_powerlaw(fast_fit_plot_range, *sol2_fast[0]),
    #         color='red', lw=1.8, linestyle='--', dashes=(2, 2),
    #         label= r"$\mathrm{\gamma_{FAST}}$"+f" = {sol2_fast[0][1]:.3f} $\pm$ {sol2_err_fast[1]:.3f} $\pm$ {b_alpha_f:.3f}")

    #Fast fits -- summing bursts within 100ms and count them as one
    fast_fit_plot_range = np.concatenate([energies_ex_fast_summed[-200:], energies_fit_fast_summed])
    # FAST fit with fit limit determinded using powerlaw.Fit
    ax0.plot(fast_fit_plot_range, func_powerlaw(fast_fit_plot_range, *sol2_fast_0_summed[0]),
             color=clrs['fast'], lw=1.8, linestyle='-',
             label = r"$\mathrm{\gamma_{FAST}}$"+f" = {sol2_fast_0_summed[0][1]:.3f} $\pm$ {sol2_err_fast_0_summed[1]:.3f} $\pm$ {b_alpha_f_0_summed:.3f}")
    # FAST fit with break point from Xu et al. 2022 used as fitting limit
    ax0.plot(fast_fit_plot_range, func_powerlaw(fast_fit_plot_range, *sol2_fast_summed[0]),
             color=clrs['fast'], lw=1.8, linestyle='--', dashes=(2, 2),
             label= r"$\mathrm{\gamma_{FAST}}$"+f" = {sol2_fast_summed[0][1]:.3f} $\pm$ {sol2_err_fast_summed[1]:.3f} $\pm$ {b_alpha_f_summed:.3f}")

    # completeness thresholds for all
    ax0.axvline(energie_threshold_on, color=clrs['on'], linestyle='dotted')
    ax0.axvline(energie_threshold_st, color=clrs['st'], linestyle='dotted')
    #ax0.axvline(energie_threshold_fast, color=clrs['fast'], linestyle='--')
    ax0.axvline(turning_point_fast_0, color=clrs['fast'], linestyle='dotted')
    ax0.axvline(turning_point_fast, color=clrs['fast'], linestyle='dashdot')
    print(f'Turning point fast: {turning_point_fast}')

    #ax0.axvline(turning_point_fast, color='#e66101', linestyle='dashdot')
    #without fit
    #ax0.scatter(np.sort(energies_array_fast), np.flip(y_points_fast), color="#e66101", marker = "2", s=16,\
     #           label = "FAST")

    ax0.set_ylabel(r"R (>E$_{\nu}$) [h$^{-1}$]")
    ax0.set_xlabel(r"E$_{\nu}$ [erg$\,$Hz$^{-1}$]")
    ax0.set_xscale('log')
    ax0.set_yscale('log')
    ax0.legend(loc=1)
    ax0.set_ylim(0.001, 40)
    ax0.set_xlim(3*1e27, 1e33)
    #ax0.savefig("../plots/burstrate_o8_fast_1.png", bbox_inches='tight', facecolor='white', transparent=False, dpi=600)
    #ax0.show()
    #ax0.clf()

    #########################
    #
    # PLOT PART two
    #
    #######################
#def plot_window3():
    "Cum disstribution of burst densities of window 3 for Onsala,Stk and Westerbork"

    #total amount of hours observed between 59602 and 59643
    energies_array_on = isotropic_energies(fluences_on_w3_arr)
    energies_array_stk2 = isotropic_energies(fluences_stk_w3_arr)
    energies_array_wb = isotropic_energies(fluences_wb_w3_arr)

    f_limit_wb = fluence_limit(sefd=420, width=3, bw=100, snr=SNR_limit)
    f_limit_stk2 = fluence_limit(sefd=385, width=3, bw=80, snr=SNR_limit)

    energie_threshold_wb = isotropic_energies([f_limit_wb])
    energie_threshold_stk2 = isotropic_energies([f_limit_stk2])

    #Onsala
    print('Fitting Onsala Epoch 2')
    ob_time_on_l2 = 9.84+40.89+7.22 # 57.95
    energies_on2, y_points_on2 = \
        power_law_fit(energies_array_on, observed_time=ob_time_on_l2, fit=False)

    #fit
#     energies_fit_on, y_points_fit_on, energies_ex_on, y_points_ex_on, xx_on, p_y_on, \
#     sol2_on, sol2_err_on, turning_point_on = power_law_fit(energies_array_on, \
#                                                            observed_time=ob_time_on_l, scale=10**18)

    #Stockert
    print('Fitting Stockert Epoch 2')
    ob_time_stk_l2 = 177.97 # 177.97
    energies_fit_stk2, y_points_fit_stk2, energies_ex_stk2, y_points_ex_stk2,\
        sol2_stk2, sol2_err_stk2, turning_point_stk2, b_alpha_stk2 = \
            power_law_fit(energies_array_stk2, observed_time=ob_time_stk_l2,
                          scale=10**18, limit=(energie_threshold_stk2, 1e33),
                          limit_b=True, ignore_turnover=True)
    #power_law_fit(energies_array_stk, observed_time=ob_time_stk_l, scale=10**18)

    #Westerbork
    print('Fitting Westerbork Epoch 2')
    ob_time_wb_l = 199.01 # 214?
    energies_fit_wb, y_points_fit_wb, energies_ex_wb, y_points_ex_wb, sol2_wb, \
        sol2_err_wb, turning_point_wb, b_alpha_wb = \
            power_law_fit(energies_array_wb, observed_time=ob_time_wb_l,
                          scale=10**18, limit=(energie_threshold_wb,9e32),
                          limit_b=True, ignore_turnover=True)

    # global fit to data from all dishes and all epochs plus obstimes for weights
    # from 1st epoch
    #energies_fit_on, y_points_fit_on, ob_time_on_l
    #energies_stk, y_points_stk, ob_time_stk_l
    # from 2nd epoch
    #energies_fit_wb, y_points_fit_wb, ob_time_wb_l
    #energies_fit_stk2, y_points_fit_stk2, ob_time_stk_l2
    #energies_on2, y_points_on2, ob_time_on_l2

    # collecting all energies and observing times for a global histogram
    # the inverse of the observing time will be the weights for the cumulative histogram
    energies_on2_limit = energies_on2[energies_on2 > energie_threshold_on]
    energies_for_hist = np.concatenate([energies_fit_on,
                                        energies_on2_limit,
                                        #energies_stk,
                                        energies_fit_stk2,
                                        energies_fit_wb])
    o8_weights_1 = np.array([1./ob_time_on_l] * len(energies_fit_on))
    o8_weights_2 = np.array([1./ob_time_on_l2] * len(energies_on2_limit))
    st_weights_1 = np.array([1./ob_time_stk_l] * len(energies_stk))
    st_weights_2 = np.array([1./ob_time_stk_l2] * len(energies_fit_stk2))
    wb_weights = np.array([1./ob_time_wb_l] * len(energies_fit_wb))

    weights_for_hist = np.concatenate([o8_weights_1, o8_weights_2,
                                       #st_weights_1,
                                       st_weights_2,
                                       wb_weights])

    # collect energies for global fit on data points.
    energies_array_all = np.concatenate([energies_fit_on_st,
                                         energies_on2[energies_on2 > energie_threshold_on],
                                         energies_fit_stk2, energies_fit_wb])
    y_points_all = np.concatenate([y_points_fit_on_st, y_points_on2[energies_on2 > energie_threshold_on],
                                   y_points_fit_stk2, y_points_fit_wb])
    print('Global fit.')
    print(f'{np.sum(energies_array_all >= 3e30)}')
    print(f'{np.sum(energies_array_all >= 1.3e32)}')
    energies_fit_all, y_points_fit_all, \
        sol2_all, sol2_err_all, b_alpha_all = \
            joint_powerlaw_fit(energies_array_all, y_points_all,
                               scale=10**18)

    #Figure
    #fig = plt.figure(figsize=(12, 10))
    ax1 = fig.add_subplot(spec[0, 1])
    #spec = fig.add_gridspec(ncols=1, nrows=2, wspace=0.0, hspace=0.0)

    #Stockert data points
    #ax0 = fig.add_subplot() #spec[0, 0])
    ax1.scatter(np.sort(energies_ex_stk2), np.flip(y_points_ex_stk2), color=clrs['st'],
                marker="o", alpha=0.5)
    ax1.scatter(np.sort(energies_fit_stk2), np.flip(y_points_fit_stk2),
                color=clrs['st'], marker = "o",
                label = "St 2022")
    #ax1.axvline(turning_point_stk, color='red', linestyle='--')
    ax1.set_ylim(0.001, 0.3)
    ax1.set_xlim(1e30, 8e32)
    #ax1.tick_params(labelbottom=False)
    ax1.axvline(energie_threshold_stk2, color=clrs['st'], linestyle='--')
    #ax1.axvline(energie_threshold_st, color=clrs['st'], linestyle='-.')

    ########

    #Westerbork data points
    #ax1 = fig.add_subplot(spec[1, 0], sharex=ax0, sharey=ax0)
    ax1.scatter(np.sort(energies_ex_wb), np.flip(y_points_ex_wb),
                color=clrs['wb'], marker="o", alpha=0.5)
    ax1.scatter(np.sort(energies_fit_wb), np.flip(y_points_fit_wb),
                color=clrs['wb'], marker = "o",
                label = "Wb 2022")

    # STocket and Wb fits
    ax1.plot(energies_fit_stk2,func_powerlaw(energies_fit_stk2,*sol2_stk2[0]),
             color=clrs['st'], lw=1.8,
             label = r"$\mathrm{\gamma_{St}}$"+f" = {sol2_stk2[0][1]:.2f} $\pm$ {sol2_err_stk2[1]:.2f} $\pm$ {b_alpha_stk2:.2f}")
    ax1.plot(energies_fit_wb, func_powerlaw(energies_fit_wb,*sol2_wb[0]),
             color=clrs['wb'], lw=1.8,
             label=r"$\mathrm{\gamma_{Wb}}$"+f" = {sol2_wb[0][1]:.2f} $\pm$ {sol2_err_wb[1]:.2f} $\pm$ {b_alpha_wb:.2f}")
    #ax1.axvline(turning_point_wb, color='green', linestyle='--')
    ax1.axvline(energie_threshold_wb, color=clrs['wb'], linestyle='--')

    # Onsala, no fit
    ax1.scatter(np.sort(energies_on2), np.flip(y_points_on2), color=clrs['on'], marker = "o",
                label = fr"O8 2022")

    #Onsala - first window (FAST window)
    #ax1.plot(energies_fit_on,func_powerlaw(energies_fit_on,*sol2_on[0]),color='blue',lw=0.8)
    ax1.scatter(np.sort(energies_ex_on), np.flip(y_points_ex_on), color=clrs['on'], marker="s", alpha=0.2)
    ax1.scatter(np.sort(energies_fit_on), np.flip(y_points_fit_on), color=clrs['on'], marker = "s",
                label = fr"O8 2021", alpha=0.2)
    ax1.axvline(energie_threshold_on, color=clrs['on'], linestyle='--')

    #Stockert - first window (FAST window)
    ax1.scatter(np.sort(energies_array_stk), np.flip(y_points_stk), color=clrs['st'], marker = "s",
                label = fr"St 2021", alpha=0.2)

    # Global fit to all data from all epochs
    #ax1.plot(energies_fit_all, func_powerlaw(energies_fit_all,*sol2_all[0]),
    #         color='gray', lw=1.8,
    #         label=r"$\mathrm{\gamma_{all}}$"+f" = {sol2_all[0][1]:.2f} $\pm$ {sol2_err_all[1]:.2f} $\pm$ {b_alpha_all:.2f}")
    ax1.set_xlabel(r"E$_{\nu}$ [erg$\,$Hz$^{-1}$]")
    ax1.legend(loc='best')
    ax1.set_xscale('log')
    ax1.set_yscale('log')

    lbl0=ax0.text(0.05,0.85,(r'${\rm {\bf %s}}$')%'a', ha='center', va='center', fontsize=12, transform = ax0.transAxes,
                  bbox=dict(facecolor='grey', alpha=0.2, edgecolor='black'))
    lbl1=ax1.text(0.05,0.85,(r'${\rm {\bf %s}}$')%'b', ha='center', va='center', fontsize=12, transform = ax1.transAxes,
                  bbox=dict(facecolor='grey', alpha=0.2, edgecolor='black'))
    #lbl0.set_path_effects([PathEffects.withStroke(linewidth=5, foreground='w')])

    #plotlabel=ax0.legend((0,0), ('a',), loc='upper left',handlelength=0,bbox_to_anchor=(0.02, 0.99), handletextpad=0,frameon=True,markerscale=0)
    #plt.gca().add_artist(plotlabel)

    #for ax in [ax0, ax1]:
    #    ax.tick_params(labelsize=20)
    #    ax.legend(loc='best', fontsize=15)
    #    ax.set_xscale('log')
    #    ax.set_yscale('log')

    #fig.text(0.03, 0.5, "R (>E; per hour)", rotation="vertical", va="center", fontsize=25)
    #plt.savefig("../plots/burstrate_wb_stk_3.png",
    #            bbox_inches='tight', facecolor='white', transparent=False, dpi=600)
    plt.savefig("../plots/burstrates.pdf",
                bbox_inches='tight', facecolor='white', transparent=False, dpi=600)
    try:
        plt.savefig("/home/pharao/git/overleaf/r67-single-dish/figures/burstrates.pdf",
                    bbox_inches='tight', facecolor='white', transparent=False, dpi=600)
    except:
        print('Could not save to overleaf git. Got a local copy though.')

    #plt.show()

    # separate plot for global average of rates; we collect all energies and their
    # respective weights into one array, sort by energies and compute the cumulative sum
    # of the weights. This we then fit with a power law.
    global_rates = np.vstack([energies_for_hist, weights_for_hist]).T # shape = (nburst, 2)
    global_rates = global_rates[global_rates[:, 0].argsort()] # sort by energies, lowest first
    energies = global_rates[:, 0]
    weights = global_rates[:, 1] # these are sorted by lowest energy first, we want the cumsum > E,
                                 # i.e. we sum the rates starting with largest E
    rates = weights[::-1].cumsum() # cumulative sum with highest E first
    # in the power_law_fit function, energies are sorted and rates are computed by
    # summing equal weights per burst, i.e. with rate for highest E first too. Flipping done in
    # the function.
    print('New global fit. (Excludes St bursts before their upgrade)')
    energies_fit_global, y_points_fit_global, energies_ex_global, \
        y_points_ex_global, sol2_global, sol2_err_global, \
        turning_point_global, b_alpha_global = \
            power_law_fit(energies=energies, observed_time=rates,
                          scale=10**18, limit_b=False, fit=True,
                          ignore_turnover=True,
                          time_is_weights=True)

    print(f'{np.sum(energies_fit_global >= 1.3e32)}')
    fig = plt.figure(figsize=(9, 6))
    #spec = fig.add_gridspec(ncols=2, nrows=1, wspace=0.1, hspace=4.0)
    ax0 = fig.add_subplot()
    ax0.scatter(np.sort(energies_fit_global), np.flip(y_points_fit_global),
                color='black', marker = "s",
                label = fr"Combined data")
    ax0.plot(energies_fit_global,
             func_powerlaw(energies_fit_global, *sol2_global[0]),
             color='black', lw=1.8,
             label = r"$\mathrm{\gamma_{global}}$"+ f" = {sol2_global[0][1]:.2f} $\pm$ {sol2_err_global[1]:.2f} $\pm$ {b_alpha_global:.2f}")
    ax0.set_xlabel(r"E$_{\nu}$ [erg$\,$Hz$^{-1}$]")
    ax0.legend(loc='best')
    ax0.set_xscale('log')
    ax0.set_yscale('log')
    ax0.set_ylabel(r"R (>E$_{\nu}$) [h$^{-1}$]")
    ax0.set_xlabel(r"E$_{\nu}$ [erg$\,$Hz$^{-1}$]")
    #ax0.legend(loc=1)
    ax1.set_ylim(0.001, 0.3)
    ax1.set_xlim(1e30, 8e32)
    plt.savefig("../plots/burstrates_global.pdf",
                bbox_inches='tight', facecolor='white', transparent=False, dpi=600)
    try:
        plt.savefig("/home/pharao/git/overleaf/r67-single-dish/figures/burstrates_global.pdf",
                    bbox_inches='tight', facecolor='white', transparent=False, dpi=600)
    except:
        print('Could not save to overleaf git. Got a local copy though.')

    return energies, weights

energies, weights = plot_window1_fast()
#plot_window3()

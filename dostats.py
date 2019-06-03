#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on %(date)s

@author: %(username)s
"""

#%%
#%pylab
#%load_ext autoreload
#%autoreload 2
#%reload_ext autoreload

#%%

#import sys
#year = np.int(sys.argv[1])
#assert 2018<=year<=2020

#%%




import os
import numpy as np

#%%

""" read sun-moon data """

import re
from datetime import datetime
#import time
from astropy.time import Time
#from scipy.stats import binned_statistic
from scipy import signal
from scipy.stats import binned_statistic_2d
from matplotlib import rcParams
from matplotlib import pyplot as plt
rcParams.update({"font.size":15})

#%%

def read_sunmoon(fp="./sqm/lhsunmoon.dat"):
        
    with open(fp, "r+", encoding="gb2312") as f:
        lines = f.readlines()
    
    data_tstr = []
    data_mask = []
    for this_line in lines[3:]:
        this_data = re.split(" +", this_line.strip())
        this_row0 = []
        this_mask = []
        ymd = this_data[0][0:4],this_data[0][4:6],this_data[0][6:8]
        for data in this_data[1:]:
            hh, mm = data.split(":")
            if hh=="55":
                # bad data
                this_mask.append(True)
                this_row0.append("{}-{}-{}T{}:{}:00.000".format(*ymd, "00", mm))
            else:
                this_mask.append(False)
                this_row0.append("{}-{}-{}T{}:{}:00.000".format(*ymd, hh, mm))
        this_row0.insert(0, "{}-{}-{}T00:00:00.000".format(*ymd))
        this_mask.insert(0, False)
        data_tstr.append(this_row0)
        data_mask.append(this_mask)
    data_tstr = np.array(data_tstr)
    data_mask = np.array(data_mask)
    data = np.ma.MaskedArray(data_tstr, data_mask)
    
    t0 = Time(data[:,0])
    t1 = Time(data[:,[1,6]])
    t2 = Time(data[:,[2,5]])
    t3 = Time(data[:,[3,4]])
    tmoon = Time(data[:,[7,8]])
    return t0, t1, t2, t3, tmoon


#t0, t1, t2, t3, tmoon = read_sunmoon()
"""
t0: local date
t1: astronomy
t2: sailing
t3: civil
tmoon: moon rise and set
"""

#%%
""" day / night """
def isdaytime(x, t1):
    """ given t1, return True if elements in x is in daytime """
    medata = t1.mjd
    xmjd = x.mjd
    ind_day = np.array([np.any((medata[:,0]<xmjd_)&(medata[:,1]>xmjd_)) for xmjd_ in xmjd])
    return ind_day

#isdaytime(Time([57755.5], format="mjd"), t1)

#%%
""" low pass filter """
def lfilter(x, y, N=8, Wn=0.05, btype='low', analog=False, output='ba'):
    b, a = signal.butter(N, Wn, btype=btype, analog=analog, output=output)
    padlen=np.int(0.3*len(x))
    yf = signal.filtfilt(b, a, y, padlen=padlen)
    return yf

#%%
""" moving std """
def moving_std(x, n_pix=4):
    xstd = np.zeros_like(x, float)
    xlen = len(x)
    for i in range(xlen):
        xstd[i] = np.std(x[np.max((0, i-n_pix)):np.min((i+n_pix, xlen))])
    return xstd

#%%

""" read sky data """
from astropy.table import Table
#sky = Table.read("/home/cham/lh/sqm/SQMReadings_DLH.txt", format="ascii")
def read_sky(fp_sky="/home/cham/lh/sqm/SQMReadings_20181205.txt"):
    """ read SQM data """
    with open(fp_sky, "r+") as f:
        lines_sky = f.readlines()
        
    lines_sky[1] = lines_sky[1].replace("Year/Month/Day", "YMD").replace("Hour/Minute/Second", "HMS").replace("(C)", "")
    sky = Table.read(lines_sky[1:], format="ascii.basic")
    return sky


#%%
""" sky brightness """
def plot_sky_brightness(tsky, sky, figfp_sky_brightness):
    fig = plt.figure(figsize=(10,7))
    ax = fig.add_subplot(111)
    # all data
    l1 = ax.plot(np.mod(tsky.jd,1), sky["MPSAS"],'k.', alpha=0.8, ms=0.2, label="all data")
    
#    # yesterday
#    t_lastnoon = Time(np.round(Time(datetime.now()).jd)-1, format="jd")
#    ind_lastday = tsky.jd>t_lastnoon.jd
#    ax.plot(np.mod(tsky[ind_lastday].jd,1), sky["MPSAS"][ind_lastday],'-m',alpha=0.8, ms=0.2)
    
    # last day
    fjd = np.floor(tsky.jd)
    if np.mod(tsky.jd[-1],1) > 0.5:
        jd_last = np.unique(fjd)[-1]
    else:
        jd_last = np.unique(fjd)[-2]
    ind_lastday = fjd==jd_last
    date_last = Time(jd_last,format="jd").isot[:10]
    l2 = ax.plot(np.mod(tsky[ind_lastday].jd,1), sky["MPSAS"][ind_lastday],'-', color="r", alpha=0.8, ms=0.2, label=date_last)
    ax.legend(loc="lower left", framealpha=0.1)
        
    _xticks = np.linspace(0,1,13)
    _xticklabels = ["{}".format(_+4) for _ in np.arange(0, 25, 2)]
    ax.set_xticks(_xticks)
    ax.set_xticklabels(_xticklabels)
    
    ax.set_ylabel("SQM sky brightness (mag)")
    ax.set_xlabel("Time (UT hours)")
    
    ax.set_ylim(24,6)
    ax.set_xlim(.2, 0.93)
    ax.set_title("SQM")
    afontsize = 20
    ax.annotate(Time.now().isot[:10]+"  @SST", xy=(0.5,0.7), xycoords="axes fraction", fontsize=afontsize, horizontalalignment="center")
    
    tstr_now = Time.now().isot[:10]
    tstr_min = tsky.min().isot[:10]
    tstr_max = tsky.max().isot[:10]
    ax.annotate("{} - {}".format(tstr_min, tstr_max), xy=(0.5,0.8), xycoords="axes fraction", fontsize=afontsize, horizontalalignment="center")
    
    fig.tight_layout()
    fig.savefig(figfp_sky_brightness)
    
    return

#%%
""" sky goodness """
def plot_sky_goodness(tsky, sky, figfp_sky_brightness):
    
    # day / night
    ind_day = isdaytime(tsky, t1)
    ind_night = np.logical_not(ind_day)
    
    # only need night sky data
    tsky = tsky[ind_night]
    sky = sky[ind_night]
    
    
    fjd0 = np.floor(Time("{:04d}-01-01T01:01:00.000".format(year), format="isot").jd)
    fjd1 = np.floor(Time("{:04d}-01-01T01:01:00.000".format(year+1), format="isot").jd)
    
    
    fig = plt.figure(figsize=(10,8))
    ax = fig.add_subplot(111)
    
    lw = 3
    threshold = 0.03
    
    time_total = 0
    time_tbd = 0
    time_down = 0
    time_work = 0
    time_good = 0
    
    jd_now = Time.now().jd
    
    for this_jd in np.arange(fjd0, fjd1):
        
        # look for evening & morning time
        this_ev, this_mn = t1[(t1.jd>this_jd)&(t1.jd<this_jd+1)]
        
        # count sky data in this night
        ind_this_night = (tsky>this_ev)&(tsky<this_mn)
        #this_night_bin = (np.linspace(0,1,10)*(this_mn-this_ev)+this_ev).jd
        #this_night_center = np.diff(this_night_bin)+this_night_bin[:-1]
        
        if 0<np.sum(ind_this_night)<10:
            # worked, but down
            this_time_total = this_mn.jd-this_ev.jd
            time_total+=this_time_total
            time_down+=this_time_total
            
            # plot background
            plt.plot([this_jd,this_jd], [this_ev.jd-this_jd,this_mn.jd-this_jd], 'gray',lw=lw)
        
        elif 10<=np.sum(ind_this_night):
            # plot background
            plt.plot([this_jd,this_jd], [this_ev.jd-this_jd,this_mn.jd-this_jd], 'gray',lw=lw)
            # low pass filter & moving std
            x = tsky[ind_this_night].jd
            y = sky["MPSAS"][ind_this_night].data
            yf = lfilter(x, y)
            ystd = moving_std(yf-y, n_pix=10)
            x_plot = this_jd*np.ones_like(x)
            ind_clear = ystd<threshold
            y_plot = np.where(ind_clear, x-this_jd, np.nan)
            
            # plot good time
            plt.plot(x_plot, y_plot, 'cyan', lw=lw)
            
            # worked
            this_time_total = this_mn.jd-this_ev.jd
            time_total+=this_time_total
            time_work+=this_time_total
            
            dx = np.hstack((np.diff(x[:2]), (np.diff(x[:-1])+np.diff(x[1:]))/2, np.diff(x[-2:])))
            time_good+=np.sum(dx[ind_clear])
            
        else:
            # no observation
            this_time_total = this_mn.jd-this_ev.jd
            time_total+=this_time_total
            if this_ev.jd > jd_now:
                time_tbd += this_time_total
            else:
                time_down += this_time_total
                
            # plot background
            if this_ev.jd > jd_now:
                plt.plot([this_jd,this_jd], [this_ev.jd-this_jd,this_mn.jd-this_jd], 'b',lw=lw)
            else:
                plt.plot([this_jd,this_jd], [this_ev.jd-this_jd,this_mn.jd-this_jd], 'r',lw=lw)
    
    ax.set_xlabel("Month")
    ax.set_ylabel("Hour(UT)")
    
    ytick_hours = np.arange(10,28,2)
    
    ax.set_yticks((ytick_hours/24)-12/24+8/24)
    ax.set_yticklabels(["{}".format(_) for _ in ytick_hours])
    ax.set_ylim(0.2,.9)
    
    #ax.vlines(jd_now, 0.2,.9, linestyle="dashed", colors="k", zorder=4, alpha=0.5)
    
    xtick_times = Time(["{:04d}-{:02d}-01T01:01:00.000".format(year, _) for _ in np.arange(1,13)], format="isot")
    xtick_fjd = np.floor(xtick_times.jd)
    ax.set_xticks(xtick_fjd)
    xtick_labels = [_[:7] for _ in xtick_times.isot]
    ax.set_xticklabels(xtick_labels, rotation=45)
    ax.set_xlim(xtick_fjd[0]-2,xtick_fjd[-1]+32)
    
    ax.set_title("Observing time (Astronomical twilights) stat @ SST")
    afontsize = 20
    ax.annotate("Done  : {:02.2f}%".format(100*(1-time_tbd/time_total)), xy=(0.1, 0.1), xycoords="axes fraction", fontsize=afontsize)
    ax.annotate("Down  : {:02.2f}%".format(100*(time_down/(time_work+time_down))), xy=(0.65, 0.1), xycoords="axes fraction", fontsize=afontsize)
    ax.annotate("Clear  : {:02.2f}%".format(100*time_good/time_work),  xy=(0.1, 0.9), xycoords="axes fraction", fontsize=afontsize)
    ax.annotate("Cloudy: {:02.2f}%".format(100*(1-time_good/time_work)), xy=(0.65, 0.9), xycoords="axes fraction", fontsize=afontsize)
    fig.tight_layout()
    
    fig.savefig(figfp_sky_goodness)
    
    return 

#%%
from matplotlib import transforms
def plot_wind(ws, wd, ttws, figfp_wind=None):
    
    ind_day = isdaytime(ttws, t1)
    # last day
    fjd = np.floor(ttws.jd)
    if np.mod(ttws.jd[-1],1) > 0.5:
        jd_last = np.unique(fjd)[-1]
    else:
        jd_last = np.unique(fjd)[-2]
    ind_lastday = fjd==jd_last
    date_last = Time(jd_last,format="jd").isot[:10]
    jd_lastmidnight = jd_last+0.5
    
    wsbins = np.arange(0,26)
    wdbins = np.linspace(0,2*np.pi,18)
    
    fig = plt.figure(figsize=(15,10))
    
    ax = fig.add_subplot(2,3,1)
    ax.hist(ws[ind_lastday], bins=wsbins, density=False, histtype="stepfilled", lw=5, color="gray", label="all data")
    ax.hist(ws[ind_lastday&ind_day], bins=wsbins, density=False, histtype="step", lw=5, color="red", label="daytime")
    ax.hist(ws[ind_lastday&~ind_day], bins=wsbins, density=False, histtype="step",  lw=5, color="blue", label="nighttime")
    ax.set_xlabel("Wind speed [m/s]")
    ax.set_ylabel("Counts")
    ax.set_title(date_last)
    ax.set_xlim(wsbins[[0,-1]])
    ax.legend()
    
    ax = fig.add_subplot(2,3,4)
    ax.hist(ws, bins=wsbins, density=False, histtype="stepfilled", lw=5, color="gray", label="all data")
    ax.hist(ws[ind_day], bins=wsbins, density=False, histtype="step", lw=5, color="red", label="daytime")
    ax.hist(ws[~ind_day], bins=wsbins, density=False, histtype="step",  lw=5, color="blue", label="nighttime")
    ax.set_xlabel("Wind speed [m/s]")
    ax.set_ylabel("Counts")
    ax.set_title("All")
    ax.set_xlim(wsbins[[0,-1]])
    ax.legend()
    
    ax = fig.add_subplot(2,3,2, projection="polar")
    plt.scatter(wd[ind_lastday&ind_day],ws[ind_lastday&ind_day],s=10, c=np.abs(ttws.jd[ind_lastday&ind_day]-jd_lastmidnight)*24, cmap=plt.cm.jet, alpha=0.8, vmin=0, vmax=12)
    ca = plt.colorbar()
    ca.set_ticks([0,12])
    ca.set_ticklabels(["nighttime","daytime"])#, rotation=90)
    ax.set_ylim(0, wsbins[-1])
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    
    ax = fig.add_subplot(2,3,5, projection="polar")
    grid_ws, grid_wd = wsbins,wdbins
    mesh_wd,mesh_ws = np.meshgrid(grid_wd,grid_ws)
    mesh_wc_day = binned_statistic_2d(
            wd[ind_day], 
            ws[ind_day], 
            ws[ind_day], 
            bins=(grid_wd,grid_ws),
            statistic="count")
    l = ax.pcolormesh(mesh_wd, mesh_ws, np.log10(mesh_wc_day[0].T), cmap=plt.cm.hot_r, vmin=0, alpha=1)#, extent=(*grid_wd[[0,-1]], *grid_ws[[0,-1]]))
    plt.colorbar(l)
    ax.set_title("log10(counts) [daytime]")
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    
    ax = fig.add_subplot(2,3,6, projection="polar")
    grid_ws, grid_wd = wsbins,wdbins
    mesh_wd,mesh_ws = np.meshgrid(grid_wd,grid_ws)
    mesh_wc_night = binned_statistic_2d(
            wd[~ind_day], 
            ws[~ind_day], 
            ws[~ind_day], 
            bins=(grid_wd,grid_ws),
            statistic="count")
    l = ax.pcolormesh(mesh_wd, mesh_ws, np.log10(mesh_wc_night[0].T), cmap=plt.cm.hot_r, vmin=0, alpha=1)#, extent=(*grid_wd[[0,-1]], *grid_ws[[0,-1]]))
    plt.colorbar(l)
    ax.set_title("log10(counts) [nighttime]")
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    
    fig.tight_layout()
    if figfp_wind is not None:
        fig.savefig(figfp_wind)
    return
    
#%%

if __name__ == "__main__":
    
    year = 2019
    
    if os.uname()[1] == "T7610":
        
        # working dir
        dir_work = "/home/cham/lh" 
        
        # sunmoon data
        datafp_sunmoon = "/home/cham/lh/sqm/lhsunmoon.dat"
        # sky brightness data
        datafp_sky = "/home/cham/lh/2019/SQMReadings.txt"
        # wind data
        datafp_wind = "/home/cham/lh/2019/weather2019.csv"
        
        # figure paths
        figfp_sky_brightness = "/home/cham/PycharmProjects/lhstat/figs/latest_sky_brightness.png"
        figfp_sky_goodness = "/home/cham/PycharmProjects/lhstat/figs/latest_sky_goodness_{}.png".format(year)
        # wind figure
        figfp_wind = "/home/cham/PycharmProjects/lhstat/figs/latest_wind_stat.png"
    
    else: # on ali server
        # working dir
        dir_work = "/root/lhstat"
        
        # sunmoon data
        datafp_sunmoon = "./data/lhsunmoon.dat"
        # sky brightness data
        datafp_sky = "./latest_data/SQMReadings.txt"
        # wind data
        datafp_wind = "./latest_data/weather2019.csv"
        
        # figure paths
        figfp_sky_brightness = "./figs/latest_sky_brightness.png"
        figfp_sky_goodness = "./figs/latest_sky_goodness_{}.png".format(year)
        # wind figure
        figfp_wind = "./figs/latest_wind_stat.png"
        
    os.chdir(dir_work)   
    
    # read sunmoon data
    t0, t1, t2, t3, tmoon = read_sunmoon(datafp_sunmoon)
    
    # read sky data
    sky = read_sky(datafp_sky)
        
    # construct t_str for sky data
    sky_tstr = [(sky["YMD"][i]+"T"+sky["HMS"][i]).replace("/","-") for i in range(len(sky))]
    tsky = Time(sky_tstr)
    
    plot_sky_brightness(tsky, sky, figfp_sky_brightness)
    plot_sky_goodness(tsky, sky, figfp_sky_goodness)
        
    # wind data
    tws = Table.read(datafp_wind, format="ascii.commented_header")
    
    ttws = Time(["{}T{}:00".format(tws["date"][i],tws["time"][i]) for i in range(len(tws))])
    
    ws = tws["wind_speed_2mins"]
    wd = tws["wind_direction"]/180*np.pi
    plot_wind(ws, wd, ttws, figfp_wind)


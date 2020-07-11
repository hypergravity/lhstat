# -*- coding: utf-8 -*-

import datetime
import glob
import os
import re

import numpy as np
from astroML.stats import binned_statistic
from astropy import table
from astropy.table import Table, Column
from astropy.time import Time
from matplotlib import colors
from matplotlib import pyplot as plt
from matplotlib import rcParams
from scipy import signal
from scipy.stats import binned_statistic_2d

rcParams.update({"font.size": 15})


def read_sunmoon(fp="./sqm/lhsunmoon.dat"):
    with open(fp, "r+", encoding="gb2312") as f:
        lines = f.readlines()

    data_tstr = []
    data_mask = []
    for this_line in lines[3:]:
        this_data = re.split(" +", this_line.strip())
        this_row0 = []
        this_mask = []
        ymd = this_data[0][0:4], this_data[0][4:6], this_data[0][6:8]
        for data in this_data[1:]:
            hh, mm = data.split(":")
            if hh == "55":
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

    t0 = Time(data[:, 0])
    t1 = Time(data[:, [1, 6]])
    t2 = Time(data[:, [2, 5]])
    t3 = Time(data[:, [3, 4]])
    tmoon = Time(data[:, [7, 8]])
    return t0, t1, t2, t3, tmoon


def read_whitelist(fp_whitelist):
    # the whitelist contains jd for good days
    with open(fp_whitelist, "r+") as f:
        s_wdate = [_.strip() for _ in f.readlines()]
    wjd0 = Time(["{}T12:00:00".format(_) for _ in s_wdate], format="isot").jd
    return wjd0


# sky = Table.read("/home/cham/lh/sqm/SQMReadings_DLH.txt", format="ascii")
def read_sky(fp_sky="/home/cham/lh/sqm/SQMReadings_20181205.txt"):
    """ read SQM data """
    with open(fp_sky, "r+") as f:
        lines_sky = f.readlines()

    lines_sky[1] = lines_sky[1].replace("Year/Month/Day", "YMD").replace("Hour/Minute/Second", "HMS").replace("(C)", "")

    # pop null lines because of adding old data
    for i in np.flipud(np.arange(len(lines_sky))):
        if len(lines_sky[i]) > 100:
            lines_sky.pop(i)
    sky = Table.read(lines_sky[1:], format="ascii.basic")
    return sky


def count_delta(t, flag, teps=1e-10):
    """
    Parameters
    ----------
    t:
        time
    flag:
        1 for good, 0 for bad
    teps:
        for pending head and tail

    Return:
    -------
    list of continued time
    """
    if np.sum(flag) > 0:
        t_pended = np.hstack((t[0] - teps, np.array(t).flatten(), t[-1] + teps))
        flag = np.array(flag, int)
        flag_pended = np.hstack((0, flag, 0))
        flag_diff = np.hstack((np.diff(flag_pended), 0))
        t_start = t_pended[flag_diff > 0]
        t_stop = t_pended[flag_diff < 0]
        t_delta = t_stop - t_start
        t_deltamax = np.max(t_delta)
        return t_start, t_stop, t_delta, t_deltamax
    else:
        return np.nan, np.nan, np.nan, 0


def isdaytime(x, t1):
    """ day / night
    given t1, return True if elements in x is in daytime
    isdaytime(Time([57755.5], format="mjd"), t1)
    """
    medata = t1.mjd
    xmjd = x.mjd
    ind_day = np.array([np.any((medata[:, 0] < xmjd_) & (medata[:, 1] > xmjd_)) for xmjd_ in xmjd])
    return ind_day


def lfilter(x, y, N=8, Wn=0.05, btype='low', analog=False, output='ba'):
    """ low pass filter """
    b, a = signal.butter(N, Wn, btype=btype, analog=analog, output=output)
    padlen = np.int(0.3 * len(x))
    yf = signal.filtfilt(b, a, y, padlen=padlen)
    return yf


def moving_std(x, n_pix=4):
    """ moving std """
    xstd = np.zeros_like(x, float)
    xlen = len(x)
    for i in range(xlen):
        xstd[i] = np.std(x[np.max((0, i - n_pix)):np.min((i + n_pix, xlen))])
    return xstd


def plot_sky_brightness(tsky, sky, figfp_sky_brightness):
    fig = plt.figure(figsize=(8, 7))
    ax = fig.add_subplot(111)
    # all data
    l1 = ax.plot(np.mod(tsky.jd, 1), sky["MPSAS"], 'k.', alpha=0.8, ms=0.2, label="all data")

    #    # yesterday
    #    t_lastnoon = Time(np.round(Time(datetime.now()).jd)-1, format="jd")
    #    ind_lastday = tsky.jd>t_lastnoon.jd
    #    ax.plot(np.mod(tsky[ind_lastday].jd,1), sky["MPSAS"][ind_lastday],'-m',alpha=0.8, ms=0.2)

    # last day
    fjd = np.floor(tsky.jd)
    if np.mod(tsky.jd[-1], 1) > 0.5:
        jd_last = np.unique(fjd)[-1]
    else:
        jd_last = np.unique(fjd)[-2]
    ind_lastday = fjd == jd_last
    date_last = Time(jd_last, format="jd").isot[:10]

    l2 = ax.plot(np.mod(tsky[ind_lastday].jd, 1), sky["MPSAS"][ind_lastday], '-', color="r", alpha=0.8, ms=0.2,
                 label=date_last)
    ax.legend(loc="lower left", framealpha=0.1)

    _xticks = np.linspace(0, 1, 13)
    _xticklabels = ["{}".format(_ + 4) for _ in np.arange(0, 25, 2)]
    ax.set_xticks(_xticks)
    ax.set_xticklabels(_xticklabels)

    ax.set_ylabel("SQM sky brightness (mag/arcsec$^2$)")
    ax.set_xlabel("Time (UT hours)")

    ax.set_ylim(24, 6)
    ax.set_xlim(.2, 0.93)
    ax.set_title("SQM")
    afontsize = 20
    # ax.annotate(Time.now().isot[:10]+"  @SST", xy=(0.5,0.7), xycoords="axes fraction", fontsize=afontsize, horizontalalignment="center")

    tstr_now = date_last  # Time.now().isot[:10]
    tstr_min = tsky.min().isot[:10]
    tstr_max = date_last  # tsky.max().isot[:10]
    ax.annotate("{} - {}".format(tstr_min, tstr_max), xy=(0.5, 0.8), xycoords="axes fraction", fontsize=afontsize,
                horizontalalignment="center")

    fig.tight_layout()
    fig.savefig(figfp_sky_brightness)

    return


def plot_sky_goodness(tsky_, sky_, year=2020, figfp_sky_goodness_fmt="./figs/latest_sky_goodness_{}.png", wjd0=[],
                      dt_filter=0, fjd_dates=None):
    figfp_sky_goodness = figfp_sky_goodness_fmt.format(year)

    # start & end of the year, say 2019
    # fjd_dates = "2018-04-01","2020-05-05"
    if fjd_dates is None:
        if year == 2018:
            fjd0 = np.floor(Time("{:04d}-03-30T12:00:00.000".format(year), format="isot").jd)
            fjd1 = np.floor(Time("{:04d}-01-01T12:00:00.000".format(year + 1), format="isot").jd)
        else:
            fjd0 = np.floor(Time("{:04d}-01-01T12:00:00.000".format(year), format="isot").jd)
            fjd1 = np.floor(Time("{:04d}-01-01T12:00:00.000".format(year + 1), format="isot").jd)
    else:
        fjd0 = np.floor(Time("{}T12:00:00.000".format(fjd_dates[0]), format="isot").jd)
        fjd1 = np.floor(Time("{}T12:00:00.000".format(fjd_dates[1]), format="isot").jd)

    # day / night in tsky
    ind_day = isdaytime(tsky_, t1)
    ind_night = np.logical_not(ind_day)

    # only need night sky data
    tsky = tsky_[ind_night]
    sky = sky_[ind_night]

    # only this year
    ind_this_year = (tsky.jd > fjd0) & (tsky.jd < fjd1)
    tsky = tsky[ind_this_year]
    sky = sky[ind_this_year]

    # flag for sky
    # down : -1
    # bad  :  0
    # good :  1
    flag_good = np.zeros(len(tsky), int)

    # init figure
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111)
    # plotting parameters
    lw = 2

    # smoothing parameters
    threshold = 0.03

    # count good/bad/down/tbd time & days
    time_total = 0
    time_tbd = 0
    time_down = 0
    time_work = 0
    time_good = 0
    count_total = 0
    count_tbd = 0
    count_down = 0
    count_good = 0
    count_bad = 0

    # now time
    jd_now = Time.now().jd

    print("------- date-time ----- src counts selected")
    for this_jd in np.arange(fjd0, fjd1):
        # look for evening & morning time
        this_ev, this_mn = t1[(t1.jd > this_jd) & (t1.jd < this_jd + 1)]

        # the whitelist case [added on 2019-09-29]
        if this_jd in wjd0:
            # plot good time
            x_plot = [this_jd, this_jd]
            y_plot = [this_ev.jd - this_jd, this_mn.jd - this_jd]
            plt.plot(x_plot, y_plot, 'cyan', lw=lw)

            this_time_total = this_mn.jd - this_ev.jd
            time_total += this_time_total
            time_work += this_time_total
            time_good += this_time_total

            # count time & days
            count_good += 1
            count_total += 1
            continue

        # count sky data in this night
        ind_this_night = (tsky > this_ev) & (tsky < this_mn)
        sub_this_night = np.where(ind_this_night)[0]
        npts = np.sum(ind_this_night)

        if npts > 0:
            # note: use sqmsrc to selecte SQM source
            u_sqmsrc, c_sqmsrc = np.unique(sky[ind_this_night]["sqmsrc"].data, return_counts=True)
            o_sqmsrc = u_sqmsrc[np.argmax(c_sqmsrc)]
            ind_this_night = (tsky > this_ev) & (tsky < this_mn) & (sky["sqmsrc"] == o_sqmsrc)
            sub_this_night = np.where(ind_this_night)[0]
            print(Time(this_jd, format="jd").isot, u_sqmsrc, c_sqmsrc, o_sqmsrc)

        if 0 < npts < 10:
            # worked, but down
            this_time_total = this_mn.jd - this_ev.jd
            time_total += this_time_total
            time_down += this_time_total

            # plot background
            plt.plot([this_jd, this_jd], [this_ev.jd - this_jd, this_mn.jd - this_jd], 'gray', lw=lw)

            # set flag
            flag_good[sub_this_night] = -1

            # count time & days
            count_down += 1
            count_total += 1

        elif 10 <= npts:
            # worked, worked

            # plot background
            plt.plot([this_jd, this_jd], [this_ev.jd - this_jd, this_mn.jd - this_jd], 'gray', lw=lw)
            # low pass filter & moving std
            x = tsky[ind_this_night].jd
            y = sky["MPSAS"][ind_this_night].data.data
            yf = lfilter(x, y)
            ystd = moving_std(yf - y, n_pix=10)
            x_plot = this_jd * np.ones_like(x)
            ind_clear = ystd < threshold
            y_plot = np.where(ind_clear, x - this_jd, np.nan)

            # count time delta
            t_start, t_stop, t_delta, t_deltamax = count_delta(x, ind_clear * 1)
            if t_deltamax > dt_filter / 24:
                # the longest duration is longer than 3 hours
                # set flag
                flag_good[sub_this_night[ind_clear]] = 1

                # plot good time
                plt.plot(x_plot, y_plot, 'cyan', lw=lw)

                # worked
                this_time_total = this_mn.jd - this_ev.jd
                time_total += this_time_total
                time_work += this_time_total

                dx = np.hstack((np.diff(x[:2]),
                                (np.diff(x[:-1]) + np.diff(x[1:])) / 2,
                                np.diff(x[-2:])))
                time_good += np.sum(dx[ind_clear])

                # count time & days
                count_good += 1
                count_total += 1

            else:
                # longest duration is shorter than 3 hours, set to bad
                # set flag
                # flag_good[sub_this_night[ind_clear]] = 0

                # plot good time
                # plt.plot(x_plot, y_plot, 'cyan', lw=lw)

                # worked
                this_time_total = this_mn.jd - this_ev.jd
                time_total += this_time_total
                time_work += this_time_total

                # dx = np.hstack((np.diff(x[:2]),
                #                 (np.diff(x[:-1]) + np.diff(x[1:])) / 2,
                #                 np.diff(x[-2:])))
                # time_good += np.sum(dx[ind_clear])

                # count time & days
                count_bad += 1
                count_total += 1

        else:
            # no observation
            this_time_total = this_mn.jd - this_ev.jd
            time_total += this_time_total
            if this_ev.jd > jd_now:
                time_tbd += this_time_total
                # count time & days
                count_tbd += 1
                count_total += 1
            else:
                time_down += this_time_total
                # count time & days
                count_down += 1
                count_total += 1

            # plot background
            if this_ev.jd > jd_now:
                plt.plot([this_jd, this_jd], [this_ev.jd - this_jd, this_mn.jd - this_jd], 'b', lw=lw)
            else:
                plt.plot([this_jd, this_jd], [this_ev.jd - this_jd, this_mn.jd - this_jd], 'r', lw=lw)

    ax.set_xlabel("Month")
    ax.set_ylabel("Hour(UT)")

    ytick_hours = np.arange(10, 28, 2)

    ax.set_yticks((ytick_hours / 24) - 12 / 24 + 8 / 24)
    ax.set_yticklabels(["{}".format(_) for _ in ytick_hours])
    ax.set_ylim(0.2, .9)

    # ax.vlines(jd_now, 0.2,.9, linestyle="dashed", colors="k", zorder=4, alpha=0.5)
    # xtick_times = Time(["{:04d}-{:02d}-01T01:01:00.000".format(year, _) for _ in np.arange(1, 13)], format="isot")

    monthticks = Time(np.arange(fjd0, fjd1 + 1), format="jd")
    ind_01 = np.array([_[8:10]=="01" for _ in monthticks.isot])

    # xtick_fjd = np.floor(xtick_times.jd)
    xtick_fjd = monthticks[ind_01]
    ax.set_xticks(xtick_fjd.jd)
    xtick_labels = [_[:7] for _ in xtick_fjd.isot]
    ax.set_xticklabels(xtick_labels, rotation=45)
    ax.set_xlim(xtick_fjd[0].jd - 2, xtick_fjd[-1].jd + 2)

    ax.set_title("Observing time (Astronomical twilights) stat @ SST")
    afontsize = 20
    ax.annotate("Done  : {:02.2f}%".format(100 * (1 - time_tbd / time_total)), xy=(0.1, 0.1), xycoords="axes fraction",
                fontsize=afontsize)
    ax.annotate("Down  : {:02.2f}%".format(100 * (time_down / (time_work + time_down))), xy=(0.65, 0.1),
                xycoords="axes fraction", fontsize=afontsize)
    ax.annotate("Clear  : {:02.2f}%".format(100 * time_good / time_work), xy=(0.1, 0.9), xycoords="axes fraction",
                fontsize=afontsize)
    ax.annotate("Cloudy: {:02.2f}%".format(100 * (1 - time_good / time_work)), xy=(0.65, 0.9), xycoords="axes fraction",
                fontsize=afontsize)
    ax.annotate("N(photometric/bad/down/tbd): {}/{}/{}/{}".format(count_good, count_bad, count_down, count_tbd),
                xy=(0.2, 0.02), xycoords="axes fraction", fontsize=afontsize)
    # ax.annotate("N(photometric/bad/down): {}/{}/{}".format(count_good, count_bad, count_down, count_tbd),
    #             xy=(0.2, 0.02), xycoords="axes fraction", fontsize=afontsize)

    # add legend
    lgood = ax.plot([0, 0], [1, 1], "-", lw=lw, color="cyan", label="good")
    lbad = ax.plot([0, 0], [1, 1], "-", lw=lw, color="gray", label="bad")
    ldown = ax.plot([0, 0], [1, 1], "-", lw=lw, color="red", label="down")
    ltbd = ax.plot([0, 0], [1, 1], "-", lw=lw, color="blue", label="tbd")
    ax.legend([lgood[0], lbad[0], ldown[0], ltbd[0]],
              ["photometric", "bad", "down", "tbd"],
              loc="upper center", framealpha=0, fontsize=afontsize * 0.6)
    # ax.legend([lgood[0], lbad[0], ldown[0]],
    #           ["photometric", "bad", "down"],
    #           loc="upper center", framealpha=0, fontsize=afontsize * 0.6)

    fig.tight_layout()

    fig.savefig(figfp_sky_goodness)

    tsky_flagged = Table([Column(tsky.jd, "jd"),
                          # Column(ind_night, "isnight"),
                          Column(flag_good, "isgood")])
    return tsky_flagged


def plot_wind(ws, wd, ttws, figfp_wind=None):
    ind_day = isdaytime(ttws, t3)
    # last day
    fjd = np.floor(ttws.jd)
    if np.mod(ttws.jd[-1], 1) > 0.5:
        jd_last = np.unique(fjd)[-1]
    else:
        jd_last = np.unique(fjd)[-2]
    ind_lastday = fjd == jd_last
    date_last = Time(jd_last, format="jd").isot[:10]
    jd_lastmidnight = jd_last + 0.5 + 2 / 24.  # @ Lenghu, midnight is 2:00 am

    wmax = np.max((np.max(ws) + 1, 26))
    wsbins = np.arange(0, wmax)
    wdbins = np.linspace(0, 2 * np.pi, 18)

    fig = plt.figure(figsize=(15, 10))

    ax = fig.add_subplot(2, 3, 1)
    ax.hist(ws[ind_lastday], bins=wsbins, density=False, histtype="stepfilled", lw=5, color="gray",
            label="all data")
    ax.hist(ws[ind_lastday & ind_day], bins=wsbins, density=False, histtype="step", lw=5, color="red",
            label="daytime")
    ax.hist(ws[ind_lastday & ~ind_day], bins=wsbins, density=False, histtype="step", lw=5, color="blue",
            label="nighttime")
    ax.set_xlabel("Wind speed [m/s]")
    ax.set_ylabel("Counts")
    ax.set_title(date_last)
    ax.set_xlim(wsbins[[0, -1]])
    ax.legend()

    ax = fig.add_subplot(2, 3, 4)
    ax.hist(ws, bins=wsbins, density=False, histtype="stepfilled", lw=5, color="gray", label="all data")
    ax.hist(ws[ind_day], bins=wsbins, density=False, histtype="step", lw=5, color="red", label="daytime")
    ax.hist(ws[~ind_day], bins=wsbins, density=False, histtype="step", lw=5, color="blue", label="nighttime")
    ax.set_xlabel("Wind speed [m/s]")
    ax.set_ylabel("Counts")
    ax.set_title("All")
    ax.set_xlim(wsbins[[0, -1]])
    ax.legend()

    ax = fig.add_subplot(2, 3, 2, projection="polar")
    plt.scatter(wd[ind_lastday & ind_day], ws[ind_lastday & ind_day], s=10,
                c=np.abs(ttws.jd[ind_lastday & ind_day] - jd_lastmidnight) * 24, cmap=plt.cm.jet, alpha=0.8, vmin=0,
                vmax=12)
    ca = plt.colorbar()
    ca.set_ticks([0, 12])
    ca.set_ticklabels(["nighttime", "daytime"])  # , rotation=90)
    ax.set_ylim(0, wsbins[-1])
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)

    ax = fig.add_subplot(2, 3, 5, projection="polar")
    grid_ws, grid_wd = wsbins, wdbins
    mesh_wd, mesh_ws = np.meshgrid(grid_wd, grid_ws)
    mesh_wc_day = binned_statistic_2d(
        wd[ind_day],
        ws[ind_day],
        ws[ind_day],
        bins=(grid_wd, grid_ws),
        statistic="count")
    l = ax.pcolormesh(mesh_wd, mesh_ws, np.log10(mesh_wc_day[0].T), cmap=plt.cm.hot_r, vmin=0,
                      alpha=1)  # , extent=(*grid_wd[[0,-1]], *grid_ws[[0,-1]]))
    plt.colorbar(l)
    ax.set_title("log10(counts) [daytime]")
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)

    ax = fig.add_subplot(2, 3, 6, projection="polar")
    grid_ws, grid_wd = wsbins, wdbins
    mesh_wd, mesh_ws = np.meshgrid(grid_wd, grid_ws)
    mesh_wc_night = binned_statistic_2d(
        wd[~ind_day],
        ws[~ind_day],
        ws[~ind_day],
        bins=(grid_wd, grid_ws),
        statistic="count")
    l = ax.pcolormesh(mesh_wd, mesh_ws, np.log10(mesh_wc_night[0].T), cmap=plt.cm.hot_r, vmin=0,
                      alpha=1)  # , extent=(*grid_wd[[0,-1]], *grid_ws[[0,-1]]))
    plt.colorbar(l)
    ax.set_title("log10(counts) [nighttime]")
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)

    fig.tight_layout()
    if figfp_wind is not None:
        fig.savefig(figfp_wind)
    return


def plot_wind_sub(ws, wd, ttws, nwdbins=12, figfp_wind=None):
    _figsize = (8, 7)

    ind_day = isdaytime(ttws, t3)
    ind_night = ~ind_day
    # last day
    fjd = np.floor(ttws.jd)
    if np.mod(ttws.jd[-1], 1) > 0.5:
        jd_last = np.unique(fjd)[-1]
    else:
        jd_last = np.unique(fjd)[-2]
    ind_lastday = fjd == jd_last
    date_last = Time(jd_last, format="jd").isot[:10]
    jd_lastmidnight = jd_last + 0.5

    wmax = np.max((np.max(ws) + 1, 26))
    wsbins = np.arange(0, wmax)
    wdp = np.copy(wd)
    wdshift = 2 * np.pi / nwdbins / 2
    wdp[wdp > (2 * np.pi - wdshift)] -= 2 * np.pi
    wdpbins = np.linspace(0, 2 * np.pi, nwdbins + 1) - wdshift

    # -----------------------------
    fig = plt.figure(figsize=_figsize)
    ax = fig.add_subplot(111)
    ax.hist(ws[ind_lastday], bins=wsbins, density=False, histtype="stepfilled", lw=5, color="gray",
            label="all data")
    # ax.hist(ws[ind_lastday&ind_day], bins=wsbins, density=False, histtype="step", lw=5, color="red", label="daytime")
    ax.hist(ws[ind_lastday & ~ind_day], bins=wsbins, density=False, histtype="step", lw=5, color="blue",
            label="nighttime")
    ax.set_xlabel("Wind speed [m/s]")
    ax.set_ylabel("Counts")
    ax.set_title(date_last)
    ax.set_xlim(wsbins[[0, -1]])
    ax.legend()
    fig.tight_layout()
    fig.savefig(figfp_wind.replace(".png", "_1.png"))

    # -----------------------------
    fig = plt.figure(figsize=_figsize)
    ax = fig.add_subplot(111)
    ax.hist(ws, bins=wsbins, density=False, histtype="stepfilled", lw=5, color="gray", label="all data")
    # ax.hist(ws[ind_day], bins=wsbins, density=False, histtype="step", lw=5, color="red", label="daytime")
    ax.hist(ws[~ind_day], bins=wsbins, density=False, histtype="step", lw=5, color="blue", label="nighttime")
    ax.set_xlabel("Wind speed [m/s]")
    ax.set_ylabel("Counts")
    ax.set_title("All")
    ax.set_xlim(wsbins[[0, -1]])
    ax.legend()
    fig.tight_layout()
    fig.savefig(figfp_wind.replace(".png", "_2.png"))

    # -----------------------------
    fig = plt.figure(figsize=_figsize)
    ax = fig.add_subplot(111, projection="polar")
    plt.scatter(wd[ind_lastday], ws[ind_lastday], s=10, c=np.abs(ttws.jd[ind_lastday] - jd_lastmidnight) * 24,
                cmap=plt.cm.jet, alpha=0.8, vmin=0, vmax=12)
    ca = plt.colorbar()
    ca.set_ticks([0, 12])
    ca.set_ticklabels(["nighttime", "daytime"])  # , rotation=90)
    ax.set_ylim(0, wsbins[-1])
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_title("Wind speed & direction [{}]".format(date_last))
    fig.tight_layout()
    fig.savefig(figfp_wind.replace(".png", "_3.png"))

    # -----------------------------
    fig = plt.figure(figsize=_figsize)
    ax = fig.add_subplot(111, projection="polar")
    grid_ws, grid_wd = wsbins, wdpbins
    mesh_wd, mesh_ws = np.meshgrid(grid_wd, grid_ws)
    mesh_wc_day = binned_statistic_2d(
        wdp[ind_day],
        ws[ind_day],
        ws[ind_day],
        bins=(grid_wd, grid_ws),
        statistic="count")
    l = ax.pcolormesh(mesh_wd, mesh_ws, np.log10(mesh_wc_day[0].T), cmap=plt.cm.hot_r, vmin=0,
                      alpha=1)  # , extent=(*grid_wd[[0,-1]], *grid_ws[[0,-1]]))
    plt.colorbar(l)
    ax.set_title("log10(counts) [daytime]")
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    fig.tight_layout()
    fig.savefig(figfp_wind.replace(".png", "_4.png"))

    # -----------------------------
    fig = plt.figure(figsize=_figsize)
    ax = fig.add_subplot(111, projection="polar")
    grid_ws, grid_wd = wsbins, wdpbins
    mesh_wd, mesh_ws = np.meshgrid(grid_wd, grid_ws)
    mesh_wc_night = binned_statistic_2d(
        wdp[~ind_day],
        ws[~ind_day],
        ws[~ind_day],
        bins=(grid_wd, grid_ws),
        statistic="count")
    l = ax.pcolormesh(mesh_wd, mesh_ws, np.log10(mesh_wc_night[0].T), cmap=plt.cm.hot_r, vmin=0,
                      alpha=1)  # , extent=(*grid_wd[[0,-1]], *grid_ws[[0,-1]]))
    plt.colorbar(l)
    ax.set_title("log10(counts) [nighttime]")
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    fig.tight_layout()
    fig.savefig(figfp_wind.replace(".png", "_5.png"))


if __name__ == "__main__":

    if os.uname()[1] in ["T7610", "MBP16.local", "localhost"]:
        # working dir
        dir_work = os.getenv("HOME") + "/PycharmProjects/lhstat"
    else:  # on ali server
        # working dir
        dir_work = "/root/lhstat"
    os.chdir(dir_work)

    # sunmoon data
    datafp_sunmoon = "./data/lhsunmoon.dat"
    # sky brightness data
    datafp_skys = ["./latest_data/SQMReadings_20180923.txt",
                   "./latest_data/SQMReadings.txt",
                   "./latest_data/sqm_ext.txt",
                   ]
    # seeing data
    datafp_seeing = "./latest_data/Seeing_Data.txt"
    # whitelist added on 2019-09-29
    datafp_whitelist = "./latest_data/whitelist"

    # flagged tsky data
    datafp_tsky_flagged_fmt = "./figs/tsky_flagged_{}.csv"

    # figure paths
    figfp_sky_brightness = "./figs/latest_sky_brightness.png"
    figfp_sky_goodness_fmt = "./figs/latest_sky_goodness_{}.png"

    """ read white list """
    wjd0 = read_whitelist(datafp_whitelist)

    """ read sunmoon data
    t0: local date
    t1: astronomical
    t2: nautical
    t3: civil
    tmoon: moon rise and set
    """

    t0, t1, t2, t3, tmoon = read_sunmoon(datafp_sunmoon)

    """ sky stats"""
    sky_list = []
    for i, _ in enumerate(datafp_skys):
        print("reading sqm {}".format(_))
        this_sky = read_sky(_)
        this_sky.add_column(Column(np.ones(len(this_sky), dtype=int) * i, "sqmsrc"))
        sky_list.append(this_sky)
    sky = table.vstack(sky_list)

    sky_tstr = [(sky["YMD"][i] + "T" + sky["HMS"][i]).replace("/", "-") for i in range(len(sky))]

    tsky = Time(sky_tstr)
    # sort data
    indsort = np.argsort(tsky.jd)
    tsky = tsky[indsort]
    sky = sky[indsort]

    # log info
    # with open("./{}.log".format(datetime.datetime.now().isoformat()), "w+") as f:
    #     f.write(" now: " + datetime.datetime.now().isoformat() + "\n")
    #     f.write(" last entry: " + tsky[-1].isot + "\n")

    plot_sky_brightness(tsky, sky, figfp_sky_brightness)
    for year in [2018, 2019, 2020]:
        print("processing sky goodness of year ", year)
        tsky_flagged = plot_sky_goodness(tsky, sky, year, figfp_sky_goodness_fmt, wjd0=wjd0, dt_filter=0)
        tsky_flagged.write(datafp_tsky_flagged_fmt.format(year), overwrite=True)

    """ wind stats """
    # wind data
    datafp_wind = "./latest_data/weather2019.csv"
    # wind figure
    figfp_wind = "./figs/latest_wind_stat.png"
    tws = Table.read(datafp_wind, format="ascii.commented_header")
    ttws = Time(["{}T{}:00".format(tws["date"][i], tws["time"][i]) for i in range(len(tws))])
    ws = tws["wind_speed"]
    wd = tws["wind_direction"] / 180 * np.pi
    plot_wind(ws, wd, ttws, figfp_wind)
    plot_wind_sub(ws, wd, ttws, nwdbins=24, figfp_wind=figfp_wind)

    """ close all figures """
    plt.close("all")
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


def plot_seeing(sws, tsws, figfp_seeing):
    _figsize = (8, 7)

    fig = plt.figure(figsize=_figsize)
    ax = fig.add_subplot(111)
    ax.hist(sws["seeing"], bins=np.linspace(0, 3, 31), histtype="bar")
    ax.set_title("Seeing stat of SST [all data]")
    ax.set_xlabel("Seeing (arcsec)")
    ax.set_ylabel("Counts")
    ax.set_xlim(0, 3)
    _ylim = ax.get_ylim()
    ax.set_ylim(*_ylim)
    ax.annotate("Median = {:.2f}\"".format(np.nanmedian(sws["seeing"])), xy=(0.5, 0.7), xycoords="axes fraction",
                fontsize=20, horizontalalignment="left")
    ax.annotate("Ntotal = {:d}".format(np.sum(np.isfinite(sws["seeing"]))), xy=(0.5, 0.6), xycoords="axes fraction",
                fontsize=20, horizontalalignment="left")
    ax.vlines(np.nanmedian(sws["seeing"]), *_ylim, colors="r", linestyles="solid")
    fig.tight_layout()
    fig.savefig(figfp_seeing.replace(".png", "_hist.png"))

    fig = plt.figure(figsize=_figsize)
    ax = fig.add_subplot(111)
    # last day
    fjd = np.floor(tsws.jd)
    # fixed bug: figure not updated if observation is not till midnight
    if np.mod(tsws.jd[-1], 1) > 0.1:
        jd_last = np.unique(fjd)[-1]
    else:
        jd_last = np.unique(fjd)[-2]
    ind_lastday = fjd == jd_last
    date_last = Time(jd_last, format="jd").isot[:10]
    ax.plot(tsws.jd[ind_lastday] - jd_last, sws["seeing"][ind_lastday], "s-", alpha=0.5, label="seeing data")
    ax.set_xticks(np.linspace(0, 1, 25))
    int_hours_ut = np.arange(25) + 4
    int_hours_ut[int_hours_ut > 24] -= 24
    ax.set_xticklabels(["{}".format(_) for _ in int_hours_ut])
    ax.set_xlim(0.3, 0.85)
    ax.set_ylim(0, 5)
    ax.set_title("Seeing stat of SST [{}]".format(date_last))
    ax.set_xlabel("Hour (UT)")
    ax.set_ylabel("Seeing (arcsec)")
    seeing_med_last = np.nanmedian(sws["seeing"][ind_lastday])
    ax.hlines(seeing_med_last, *ax.get_xlim(), colors="r", linestyle="dashed",
              label="Median={:.2f}\"".format(seeing_med_last))
    ax.legend()
    fig.tight_layout()
    fig.savefig(figfp_seeing.replace(".png", "_last_seeing.png"))

    fig = plt.figure(figsize=_figsize)
    ax = fig.add_subplot(111)
    # last day
    fjd = np.floor(tsws.jd)
    if np.mod(tsws.jd[-1], 1) > 0.1:
        jd_last = np.unique(fjd)[-1]
    else:
        jd_last = np.unique(fjd)[-2]
    ind_lastday = fjd == jd_last
    date_last = Time(jd_last, format="jd").isot[:10]
    ax.plot(tsws.jd[ind_lastday] - jd_last, sws["col4"][ind_lastday], "s-", alpha=0.5)
    ax.set_xticks(np.linspace(0, 1, 25))
    int_hours_ut = np.arange(25) + 4
    int_hours_ut[int_hours_ut > 24] -= 24
    ax.set_xticklabels(["{}".format(_) for _ in int_hours_ut])
    ax.set_xlim(0.3, 0.85)
    # ax.set_ylim(0,3)
    ax.set_title("DIMM target flux [{}]".format(date_last))
    ax.set_xlabel("Hour (UT)")
    ax.set_ylabel("Flux")
    fig.tight_layout()
    fig.savefig(figfp_seeing.replace(".png", "_last_flux.png"))

    return


def plot_dust():
    """ plot dust """
    _ystday = Time(datetime.datetime.now()) - 1
    datafp_dust = "./latest_data/dust/measurement_{}-dM.dat".format(
        _ystday.isot[:10])
    figfp_dust = "./figs/latest_dust.png"

    # read data (yesterday)
    _dust_size = np.array([
        0.25, 0.28, 0.3, 0.35, 0.4, 0.45, 0.5, 0.58, 0.65,
        0.7, 0.8, 1., 1.3, 1.6, 2., 2.5, 3., 3.5,
        4., 5., 6.5, 7.5, 8.5, 10., 12.5, 15., 17.5,
        20., 25., 30., 32.])
    with open(datafp_dust, "r+") as _f:
        s = _f.readlines()
    t_dust = Table.read(s, format="ascii.no_header", delimiter="\t")
    t_dust.remove_column("col1")
    data_dust = np.array(t_dust.to_pandas())
    dust_mean = np.mean(data_dust, axis=0)
    dust_err = np.abs(np.percentile(data_dust, q=[25, 75], axis=0) - dust_mean)

    # plot figure
    fig = plt.figure(figsize=(8, 7))
    ax = fig.add_subplot(111)
    ax.errorbar(np.log10(_dust_size), dust_mean, yerr=dust_err,
                marker='s', ms=10, mfc='#1f77b4', mec='#1f77b4',
                ecolor="gray", elinewidth=5, label="daily mean$\\pm1$quantile")
    ax.legend(loc="upper right")
    ax.set_ylim(-0.01, ax.get_ylim()[1])
    _xticks = [0.2, 0.5, 1, 2, 5, 10, 20]
    ax.set_xticks(np.log10(_xticks))
    ax.set_xticklabels(["{}".format(_) for _ in _xticks])
    ax.set_xlabel("Particle Size [$\\mu$m]")
    ax.set_ylabel("Particle Counts [$\\mu$g m$^{-3}$]")
    ax.set_title("Particle size spectrum @SST [{}]".format(_ystday.isot[:10]))
    fig.tight_layout()
    # savefig
    if os.path.exists(figfp_dust):
        os.remove(figfp_dust)
    fig.savefig(figfp_dust)


def plot_aqi_daily():
    """ plot dust PM10 """
    # glob PM10 data
    _today = Time(datetime.datetime.now())
    datafp_pm10 = glob.glob("./latest_data/dust/measurement_*-M.dat")
    datafp_pm10.sort()
    figfp_pm10 = "./figs/latest_aqi_daily.png"

    data_pm10 = Table.read(datafp_pm10[-2], format="ascii.no_header", delimiter="\t",
                           names=["t", "pm10", "pm2p5", "pm1p0"])
    t_isot = []
    for i in range(len(data_pm10)):
        s = data_pm10["t"][i]
        s_splitted = [np.int(str(_)) for _ in re.split(r"[\s\b-\/:]", s)]
        t_isot.append("{:04d}-{:02d}-{:02d}T{:02d}:{:02d}:{:02d}".format(*s_splitted))
    tmjd = Time(t_isot, format="isot").mjd
    tmjd1 = tmjd - np.floor(tmjd)
    thour = tmjd1 * 24

    fig = plt.figure(figsize=(16 * .7, 9 * .7))
    ax = fig.add_subplot(111)
    l100, = ax.semilogy(thour, data_pm10["pm10"], ".-", alpha=0.5, label="PM10")
    l025, = ax.semilogy(thour, data_pm10["pm2p5"], ".-", alpha=0.5, label="PM2.5")
    l010, = ax.semilogy(thour, data_pm10["pm1p0"], ".-", alpha=0.5, label="PM1.0")
    ax.legend()
    ax.set_xlim(0, 24)
    ax.set_xticks(np.linspace(0, 24, 13))
    ax.set_ylim(0.1, 1000)
    ax.set_xlabel("Hour")
    ax.set_ylabel("Dust grain concentration [$\\mu$g m$^{-3}$]")
    ax.set_title("AQI of SST: {}".format(os.path.basename(datafp_pm10[-2])))
    fig.tight_layout()

    # savefig
    if os.path.exists(figfp_pm10):
        os.remove(figfp_pm10)
    fig.savefig(figfp_pm10)
    return


def plot_aqi_stats():
    """ plot dust PM10 """
    # glob PM10 data
    _today = Time(datetime.datetime.now())
    datafp_pm10 = glob.glob("./latest_data/dust/measurement_*-M.dat")
    datafp_pm10.sort()
    figfp_pm10 = "./figs/latest_aqi_stats.png"

    # read PM10 data
    data_pm10 = []
    for _ in datafp_pm10:
        try:
            data_pm10.append(Table.read(_, format="ascii.no_header", delimiter="\t"))
        except Exception:
            # 2019-12-27 file has a header
            if os.path.basename(_) == "measurement_2019-12-27-M.dat":
                with open(_, "r+") as f:
                    s = f.readlines()
                    data_pm10.append(Table.read(s[14:], format="ascii.no_header", delimiter="\t"))
            else:
                raise RuntimeError("error occurs when parsing", _)
            
    data_pm10 = table.vstack(data_pm10)
    data_pm10.rename_columns(["col1", "col2", "col3", "col4", ], ["t", "pm10", "pm2p5", "pm1p0"])
    # fix time string
    t_isot = []
    for i in range(len(data_pm10)):
        s = data_pm10["t"][i]
        s_splitted = [np.int(str(_)) for _ in re.split(r"[\s\b-\/:]", s)]
        t_isot.append("{:04d}-{:02d}-{:02d}T{:02d}:{:02d}:{:02d}".format(*s_splitted))
    data_pm10.add_column(Time(t_isot), name="t_isot")

    pm10 = data_pm10["pm10"].data
    tjd = data_pm10["t_isot"].jd
    jd_min = np.ceil(np.min(tjd))
    jd_max = np.floor(np.max(tjd))
    jd_x = np.arange(jd_min, jd_max + 1, 1, dtype=float)
    jd_edges = np.arange(jd_min - .5, jd_max + 1.5, 1, dtype=float)

    # quantiles
    bs_median = binned_statistic(tjd, pm10, statistic="median", bins=jd_edges)[0]
    bs_mean = binned_statistic(tjd, pm10, statistic="mean", bins=jd_edges)[0]
    bs_max = binned_statistic(tjd, pm10, statistic=np.max, bins=jd_edges)[0]
    bs_min = binned_statistic(tjd, pm10, statistic=np.min, bins=jd_edges)[0]
    bs_16 = binned_statistic(tjd, pm10, statistic=lambda x: np.percentile(x, 16), bins=jd_edges)[0]
    bs_84 = binned_statistic(tjd, pm10, statistic=lambda x: np.percentile(x, 84), bins=jd_edges)[0]

    # read lenghu PM10 data
    data_pm10_lh = Table.read("./latest_data/lhdust02.dat", format="ascii.no_header")
    tjd_lh = Time(["{}-{}-{}T12:00:00".format(str(_)[:4], str(_)[4:6], str(_)[6:8])
                   for _ in data_pm10_lh["col1"]], format="isot").jd
    pm10_lh = data_pm10_lh["col2"].data

    fig = plt.figure(figsize=(16 * .7, 9 * .7))
    ax = fig.add_subplot(111)
    # plot sst aqi
    llh, = ax.plot(tjd_lh, np.log10(pm10_lh), "rs-", alpha=.5)
    # fit a line
    p = np.polyfit(tjd_lh, np.log10(pm10_lh), deg=1)
    lfit_lh, = ax.plot(tjd_lh, np.polyval(p, tjd_lh), "r--", lw=3)

    # plot lenghu aqi
    lmedian, = ax.plot(jd_x, np.log10(bs_median), 's-', color="k", alpha=.9)
    lmean, = ax.plot(jd_x, np.log10(bs_mean), 's-', color="b", alpha=.5)
    l16 = ax.fill_between(jd_x, np.log10(bs_16), np.log10(bs_84), color="k", alpha=.4)
    l00 = ax.fill_between(jd_x, np.log10(bs_min), np.log10(bs_16), color="k", alpha=.1)
    l00 = ax.fill_between(jd_x, np.log10(bs_84), np.log10(bs_max), color="k", alpha=.1)
    # fit a line
    p = np.polyfit(jd_x, np.log10(bs_median), deg=1)
    lfit_median, = ax.plot(jd_x, np.polyval(p, jd_x), "k--", lw=3)
    p = np.polyfit(jd_x, np.log10(bs_mean), deg=1)
    lfit_mean, = ax.plot(jd_x, np.polyval(p, jd_x), "b--", lw=3, alpha=.5)

    ax.set_xlabel("Date")
    ax.set_ylabel("AQI(PM10 $\\mu$g m$^{-3}$)")
    ax.legend([llh, lmedian, lmean, l16, l00],
              ["LH daily mean", "SST daily median", "SST daily mean",
               "SST daily 16/84th pct", "SST daily min/max"],
              framealpha=0, fontsize=10, loc="upper left")
    ax.set_xlim(jd_min - .2, jd_max + .2)

    jd2dates = np.array([_[:10] for _ in Time(jd_x, format="jd").isot])
    jd_ticks_ind = [True if _[-2:] in ["01", "11", "21"] else False for _ in jd2dates]
    ax.set_xticks(jd_x[jd_ticks_ind])
    ax.set_xticklabels(jd2dates[jd_ticks_ind], rotation=45, fontsize=10)
    ax.set_yticks(np.linspace(-1, 3, 5))
    ax.set_yticklabels(["0.1", "1.0", "10", "100", "1000"])
    ax.set_title("Dust daily stats: Lenghu vs SST-site [{}]".format(_today.isot[:10]))
    fig.tight_layout()

    # savefig
    if os.path.exists(figfp_pm10):
        os.remove(figfp_pm10)
    fig.savefig(figfp_pm10)


def plot_dust_ts():
    """ plot dust DEPRECATED """
    _ystday = Time(datetime.datetime.now()) - 1
    datafp_dust = "./latest_data/dust/measurement_{}-dM.dat".format(
        _ystday.isot[:10])
    figfp_dust = "./figs/latest_dust_ts.png"

    # read data (yesterday)
    _dust_size = np.array([
        0.25, 0.28, 0.3, 0.35, 0.4, 0.45, 0.5, 0.58, 0.65,
        0.7, 0.8, 1., 1.3, 1.6, 2., 2.5, 3., 3.5,
        4., 5., 6.5, 7.5, 8.5, 10., 12.5, 15., 17.5,
        20., 25., 30., 32.])
    with open(datafp_dust, "r+") as _f:
        s = _f.readlines()
    t_dust = Table.read(s, format="ascii.no_header", delimiter="\t")
    t_time = Time([_.replace("/", "-") for _ in t_dust["col1"].data], format="iso").mjd - np.floor(_ystday.mjd)
    t_dust.remove_column("col1")
    data_dust = np.array(t_dust.to_pandas())
    dust_mean = np.mean(data_dust, axis=0)
    dust_err = np.abs(np.percentile(data_dust, q=[25, 75], axis=0) - dust_mean)
    scalar_map = plt.cm.ScalarMappable(
        norm=colors.Normalize(vmin=np.log10(np.min(_dust_size)), vmax=np.log10(np.max(_dust_size))),
        cmap=plt.cm.RdYlBu_r)
    fig = plt.figure(figsize=(8, 7))
    ax = fig.add_subplot(111)
    for i in range(len(_dust_size)):
        this_dust = data_dust[:, i]
        plt.plot(t_time, this_dust, "-", c=this_dust, color=scalar_map.to_rgba(np.log10(_dust_size[i])), alpha=0.3)

    ax.set_xlim(-0.001, 1.001)
    ax.set_ylim(-0.01, 0.5)
    _xticks = [0.2, 0.5, 1, 2, 5, 10, 20]
    # ax.set_xticks(np.log10(_xticks))
    # ax.set_xticklabels(["{}".format(_) for _ in _xticks])
    ax.set_xlabel("Time/Day")
    ax.set_ylabel("Particle Counts [$\\mu$g m$^{-3}$]")
    ax.set_title("Particle time series @SST [{}]".format(_ystday.isot[:10]))
    fig.tight_layout()
    # savefig
    if os.path.exists(figfp_dust):
        os.remove(figfp_dust)
    fig.savefig(figfp_dust)


if __name__ == "__main__":

    year = 2019

    if os.uname()[1] in ["T7610", "MBP16.local"]:
        # working dir
        dir_work = os.getenv("HOME") + "/PycharmProjects/lhstat"
    else:  # on ali server
        # working dir
        dir_work = "/root/lhstat"
    os.chdir(dir_work)

    # #################### common fps #####################
    # seeing data
    datafp_seeing = "./latest_data/Seeing_Data.txt"
    # seeing figure
    figfp_seeing = "./figs/latest_seeing_stat.png"
    # ######################################################

    """ seeing stats """
    sws = Table.read(datafp_seeing, format="ascii.no_header", delimiter="|",
                     names=['datetime1', 'datetime2', 'col3', 'col4', 'col5', 'seeing', 'col7', 'col8', 'col9'])
    tsws = Time(
        [sws["datetime2"][i].replace(" ", "T").replace("/", "-")
         for i in range(len(sws))])
    plot_seeing(sws, tsws, figfp_seeing)

    """ dust stats """
    plot_dust()
    plot_aqi_daily()
    plot_aqi_stats()

    """ close all figures """
    plt.close("all")

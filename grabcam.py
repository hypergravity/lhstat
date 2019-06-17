import os
import numpy as np

# %%

""" read sun-moon data """

import re
from datetime import datetime
# import time
from astropy.time import Time, TimeDelta
# from scipy.stats import binned_statistic
from scipy import signal
from scipy.stats import binned_statistic_2d
from matplotlib import rcParams
from matplotlib import pyplot as plt
from astropy.table import Table
import shutil

rcParams.update({"font.size": 15})


# %%

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


def picname2jd(fp):
    return Time(os.path.basename(fp).replace(".jpg",""), format="iso").jd


if __name__ == "__main__":

    if os.uname()[1] == "T7610":
        # working dir
        dir_work = "/home/cham/PycharmProjects/lhstat"

        # sunmoon data
        datafp_sunmoon = "./data/lhsunmoon.dat"

    elif os.uname()[1] == "50bindataprocess":
        # working dir
        dir_work = "/home/cham/lhstat"

        # sunmoon data
        datafp_sunmoon = "./data/lhsunmoon.dat"

    else:  # on ali server
        # working dir
        dir_work = "/root/lhstat"

        # sunmoon data
        datafp_sunmoon = "./data/lhsunmoon.dat"

    os.chdir(dir_work)
    t0, t1, t2, t3, tmoon = read_sunmoon(datafp_sunmoon)

    # make directory
    dir_temp = "./tempcam"
    if os.path.exists(dir_temp):
        os.system("rm -rf {}".format(dir_temp))
    os.mkdir(dir_temp)

    # today, noon
    dt_today = datetime.now()
    t_today = Time("{:04d}-{:02d}-{:02d}T12:00:00".format(
        dt_today.year, dt_today.month, dt_today.day), format="isot")
    t_yesterday = t_today- TimeDelta(1, format="jd")
    dt_yesterday = t_yesterday.to_datetime()

    # jd of ev & momonts of last night
    from bisect import bisect
    i_ = bisect(t1.jd[:,0], t_today.jd)
    jd_lastnight = (t1.jd[i_-1, 1], t1.jd[i_, 0])

    # glob pics
    from glob import glob
    fps = glob("/data/lh/allskycam/{:04d}/{:02d}/{:02d}/*.jpg".format(dt_today.year,dt_today.month,dt_today.day))
    fps.extend(glob("/data/lh/allskycam/{:04d}/{:02d}/{:02d}/*.jpg".format(dt_yesterday.year,dt_yesterday.month,dt_yesterday.day)))

    # time of cam
    jd_pics = np.array([picname2jd(fp) for fp in fps])

    # validate time
    jd_pics_flag = (jd_pics>=jd_lastnight[0])&(jd_pics<=jd_lastnight[1])

    for i in range(len(jd_pics_flag)):
        if jd_pics_flag[i]:
            shutil.copyfile(fps[i],"{}/{}".format(dir_temp, os.path.basename(fps[i])))
            print("@lhstat: copying file [{}] ...".format(fps[i]))

    print("@lhstat: Done, pics stored in {}!".format(dir_temp))



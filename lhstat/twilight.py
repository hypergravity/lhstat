# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 19:13:09 2020

@author: cham

adapted from:
    https://stackoverflow.com/questions/2637293/calculating-dawn-and-sunset-times-using-pyephem
    
PyEphem doc:
    https://rhodesmill.org/pyephem/rise-set.html
"""

# %%
import ephem
import numpy as np
from astropy.time import Time, TimeDelta
from astropy import units
from astropy.table import Table
from collections import OrderedDict


def eval_twilight(date="2020-03-10 04:00:00", lon=120, lat=30, elev=0., gmt=8, kind="astronomical"):
    # Make an observer
    fred = ephem.Observer()

    # PyEphem takes and returns only UTC times. 15:00 is noon in Fredericton
    # 4:00 is noon in Beijing
    fred.date = date

    # Location of Fredericton, Canada
    fred.lon = str(lon)  # Note that lon should be in string format
    fred.lat = str(lat)  # Note that lat should be in string format
    # Elevation of Fredericton, Canada, in metres
    fred.elev = elev

    # To get U.S. Naval Astronomical Almanac values, use these settings
    fred.pressure = 0
    fred.horizon = '-0:34'

    # sunrise = fred.previous_rising(ephem.Sun())  # Sunrise
    # noon = fred.next_transit(ephem.Sun(), start=sunrise)  # Solar noon
    # sunset = fred.next_setting(ephem.Sun())  # Sunset

    # We relocate the horizon to get twilight times
    # -6=civil twilight, -12=nautical, -18=astronomical
    assert kind in ["civil", "nautical", "astronomical"]
    if kind == "civil":
        fred.horizon = '-6'
    if kind == "nautical":
        fred.horizon = '-12'
    if kind == "astronomical":
        fred.horizon = '-18'

    beg_twilight = fred.previous_rising(ephem.Sun(), use_center=True)  # Begin civil twilight
    end_twilight = fred.next_setting(ephem.Sun(), use_center=True)  # End civil twilight

    # print(beg_twilight.datetime(), end_twilight.datetime())
    return (Time([beg_twilight.datetime(), end_twilight.datetime()]) + TimeDelta(gmt*units.h)).isot


def generate_sunmoon(year0=2017, year1=2022, lon=93.8961, lat=38.6068, elev=4200, gmt=8):
    assert year0 <= year1
    td = TimeDelta(gmt*units.h)
    jd0 = (Time("{}-01-01T12:00:00".format(year0), ) - td).jd
    jd1 = (Time("{}-12-31T12:00:00".format(year1), ) - td).jd
    jd_array = np.arange(jd0, jd1+1, 1, dtype=float)
    t_array = Time(jd_array, format="jd")
    t_local_array = t_array + td

    results = []
    for i in range(len(t_array)):
        this_result = OrderedDict()
        this_result["noon"] = t_local_array[i].isot
        this_result["sunrise_civil"], this_result["sunset_civil"] = \
            eval_twilight(date=t_array[i].datetime, lon=lon, lat=lat, elev=elev, gmt=gmt, kind="civil")
        this_result["sunrise_nauti"], this_result["sunset_nauti"] = \
            eval_twilight(date=t_array[i].datetime, lon=lon, lat=lat, elev=elev, gmt=gmt, kind="nautical")
        this_result["sunrise_astro"], this_result["sunset_astro"] = \
            eval_twilight(date=t_array[i].datetime, lon=lon, lat=lat, elev=elev, gmt=gmt, kind="astronomical")
        results.append(this_result)
    return Table(results)


if __name__ == "__main__":
    # print(eval_twilight())
    sunmoon = generate_sunmoon(2017, 2022,)
    sunmoon.write()
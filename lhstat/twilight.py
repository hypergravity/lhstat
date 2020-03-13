# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 19:13:09 2020

@author: cham

adapted from:
    https://stackoverflow.com/questions/2637293/calculating-dawn-and-sunset-times-using-pyephem
    
PyEphem doc:
    https://rhodesmill.org/pyephem/rise-set.html
"""

#%% 
import ephem


def twilight_gmtp8(date="2020-03-10 04:00:00", lon=120, lat=30, elev=0.,kind="astronomical"):

    #Make an observer
    fred = ephem.Observer()
    
    #PyEphem takes and returns only UTC times. 15:00 is noon in Fredericton
    fred.date = date
    
    #Location of Fredericton, Canada
    fred.lon  = str(lon) #Note that lon should be in string format
    fred.lat  = str(lat)      #Note that lat should be in string format
    #Elevation of Fredericton, Canada, in metres
    fred.elev = elev
    
    #To get U.S. Naval Astronomical Almanac values, use these settings
    fred.pressure= 0
    fred.horizon = '-0:34'
    
    sunrise=fred.previous_rising(ephem.Sun()) #Sunrise
    noon =fred.next_transit(ephem.Sun(), start=sunrise) #Solar noon
    sunset =fred.next_setting(ephem.Sun()) #Sunset
    
    # We relocate the horizon to get twilight times
    # -6=civil twilight, -12=nautical, -18=astronomical
    assert kind in ["civial", "nautical", "astronomical"]
    if kind == "civil":
        fred.horizon = '-6'
    if kind == "nautical":
        fred.horizon = '-12'
    if kind == "astronomical":
        fred.horizon = '-18'
        
    beg_twilight=fred.previous_rising(ephem.Sun(), use_center=True) #Begin civil twilight
    end_twilight=fred.next_setting(ephem.Sun(), use_center=True) #End civil twilight
    
    # print(beg_twilight.datetime(), end_twilight.datetime())
    return beg_twilight.datetime(), end_twilight.datetime()

if __name__ == "__main__":
    print(twilight_gmtp8())

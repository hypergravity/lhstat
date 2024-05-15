# -*- codeing = utf-8 -*-
# @Time : 2021/12/8 7:31
# @Author : Ava
# @File : pwv_new.py
# @Software : PyCharm
import pandas as pd
from datetime import datetime,timedelta
import numpy as np
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import pyecharts.options as opts
from pyecharts.charts import Line
from astral.sun import sun
from astral import LocationInfo
import time

filename = r'/root/lhstat/latest_data/weather2019.csv'
df = pd.read_csv(filename, sep='\\s+', skiprows=[0], header=None)
df.columns = ['date', 'time', 'temperature', 'humidity', 'dew_temperature', 'pressure', 'wind_speed',
              'wind_speed_2mins', 'wind_speed_10mins', 'wind_direction']
df['date&time'] = df['date'] + " " + df['time']
df['date&time'] = df['date&time'].apply(lambda x: datetime.strptime(x, '%Y-%m-%d %H:%M'))

latest2=datetime.strptime(df.iloc[-1]['date'], '%Y-%m-%d')
latest1=latest2-timedelta(days=1)
city= LocationInfo('LengHu', 'China', 'Asia/ShangHai',  38.6068,93.8961)
s1 = sun(city.observer, date=datetime.date(latest1), dawn_dusk_depression=18.0,tzinfo=city.timezone)
s2 = sun(city.observer, date=datetime.date(latest2), dawn_dusk_depression=18.0,tzinfo=city.timezone)
dusk=s1["dusk"].strftime('%Y-%m-%d %H:%M')
dawn=s2["dawn"].strftime('%Y-%m-%d %H:%M')
print('Dusk: '+str(s1["dusk"])+'\n'+'Dawn: '+str(s2["dawn"]))
df=df[(df['date&time'] >= dusk) & (df['date&time'] <= dawn)]
df = df.reset_index(drop=True)

df['date&time']=df['date&time']-timedelta(hours=8)
df['date'] = pd.to_datetime(df['date&time']).dt.date
df['time'] = pd.to_datetime(df['date&time']).dt.time
print(df)

x=df['date&time']
temperature=df['temperature']
pressure=df['pressure']
humidity=df['humidity']
p=np.array(pressure)
t=np.array(temperature)
RH=np.array(humidity)
PWV=np.zeros(len(p))
g=9.80665
rou=1000.0
cin1=[i for i in range(1,681)]
cin=np.array(cin1)
tcp=680
for n in range (0,len(p)):
    PWV[n] = 0
    tloop=t[n]+273.15-6.5*(cin-1)/100
    ploop=p[n]*(tloop/(t[n]+273.15))**(9.8/6.5/0.28705287)
    cdloop = 10.79574* (1 - 273.16/ tloop) - 5.02800* np.log10(tloop / 273.16) + 1.50475 * 0.0001* (1 - 10**(-8.2969 * (tloop / 273.16 - 1))) + 0.42873 * 0.001 * (10** (4.76955 * (1 - 273.16 / tloop)) - 1) + 0.78614
    esdloop = 10**cdloop * RH[n]
    pwvloop=0.622*esdloop/(ploop-0.378*esdloop)/rou/g*1000.0
    PWV[n]=-0.5*sum((pwvloop[1:tcp]+pwvloop[0:tcp-1])*(ploop[1:tcp]-ploop[0:tcp-1]))
    # print(n)

PWV=pd.DataFrame(PWV)
a=df['date&time'][0].strftime('%Y')
b=df['date&time'][0].strftime('%m')
c=df['date&time'][0].strftime('%d')
# print(PWV)
pd.set_option("display.max_columns",300)
df['PWV']=PWV
y=df['PWV'].apply(lambda x:round(x,5))
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.xaxis.set_major_formatter(mdates.DateFormatter('%H'))
ax.xaxis.set_major_locator(mdates.HourLocator())
plt.plot(x,y,marker='s',alpha = 0.5,markersize=5)
plt.xlabel('hour (UT)')
plt.ylabel('PWV (mm)')
plt.title('Precipitable water vapour [%s-%s-%s]'%(a,b,c))
plt.savefig("/root/lhstat/figs/pwv")
plt.show()


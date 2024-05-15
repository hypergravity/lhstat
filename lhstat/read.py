#!/root/anaconda3/bin/python3
# -*- codeing = utf-8 -*-
# @Time : 2021/12/20 6:12
# @Author : Ava
# @File : read.py
# @Software : PyCharm
import pandas as pd
from datetime import datetime,timedelta
from astral.sun import sun
from astral import LocationInfo


def read(filename):
    df = pd.read_csv(filename, sep='\\s+', skiprows=[0], header=None)
    # df.drop(df.head(800000).index,inplace=True)
    # df=df.reset_index(drop=True)
    df.columns = ['date', 'time', 'temperature', 'humidity', 'dew_temperature', 'pressure', 'wind_speed',
                  'wind_speed_2mins', 'wind_speed_10mins', 'wind_direction']
    df['date&time'] = df['date'] + " " + df['time']
    df['date&time'] = df['date&time'].apply(lambda x: datetime.strptime(x, '%Y-%m-%d %H:%M'))
    latest2 = datetime.strptime(df.iloc[-1]['date'], '%Y-%m-%d')
    latest1 = latest2 - timedelta(days=1)
    df = df[(df['date&time'] >= latest1 + timedelta(hours=12)) & (df['date&time'] <= latest2 + timedelta(hours=12))]
    df = df.reset_index(drop=True)
    df_utc = df.copy(deep=True)
    # df_utc['date&time'] = df_utc['date&time'] - timedelta(hours=8)
    df_utc['date'] = pd.to_datetime(df_utc['date&time']).dt.date
    df_utc['time'] = pd.to_datetime(df_utc['date&time']).dt.time
    df_utc.to_csv('latest_data_all.csv', sep=' ', index=None)

    city = LocationInfo('LengHu', 'China', 'Asia/ShangHai', 38.6068, 93.8961)
    s1 = sun(city.observer, date=datetime.date(latest1), dawn_dusk_depression=18.0, tzinfo=city.timezone)
    s2 = sun(city.observer, date=datetime.date(latest2), dawn_dusk_depression=18.0, tzinfo=city.timezone)
    dusk = s1["dusk"].strftime('%Y-%m-%d %H:%M')
    dawn = s2["dawn"].strftime('%Y-%m-%d %H:%M')
    df = df[(df['date&time'] >= dusk) & (df['date&time'] <= dawn)]
    df = df.reset_index(drop=True)
    # df['date&time'] = df['date&time'] - timedelta(hours=8)
    df['date'] = pd.to_datetime(df['date&time']).dt.date
    df['time'] = pd.to_datetime(df['date&time']).dt.time
    df.to_csv('latest_data.csv', sep=' ', index=None)

filename = r'/root/lhstat/latest_data/weather2019.csv'
read(filename)
#!/usr/bin/env bash
cd /root/lhstat/
rm ./latest_data/*
scp root@159.226.170.49:/data/lh/sqm/weather2019.csv ./latest_data/
scp root@159.226.170.49:/data/lh/sqm/SQMReadings.txt ./latest_data/
scp root@159.226.170.49:/data/lh/seeing/Seeing_Data.txt ./latest_data/

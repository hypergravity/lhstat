#!/usr/bin/env bash
cd /root/lhstat/
# rm ./latest_data/*
rsync -avz root@159.226.170.49:/data/lh/sqm/weather2019.csv ./latest_data/
# SQM
rsync -avz root@159.226.170.49:/data/lh/sqm/SQMReadings.txt ./latest_data/
rsync -avz root@159.226.170.49:/data/lh/sqm/SQMReadings_20180923.txt ./latest_data/
rsync -avz root@159.226.170.49:/data/lh/sqm/sqm_ext.txt ./latest_data/
# SQM white list
rsync -avz root@159.226.170.49:/data/lh/sqm/whitelist ./latest_data/
# seeing
rsync -avz root@159.226.170.49:/data/lh/seeing/Seeing_Data.txt ./latest_data/
rsync -avz root@159.226.170.49:/data/lh/seeing/Seeing_Data.txt ./latest_data/
# dust
rsync -avrz root@159.226.170.49:/data/lh/dust ./latest_data/
rsync -avrz root@159.226.170.49:/data/lh/lhdust02.dat ./latest_data/

#!/bin/bash

# cd working directory
cd /root/lhstat

# update weather data
bash update_data_ecs.sh

# do stats
ipython dostats.py > ./figs/`date +%Y%m%d`/`date +%Y%m%d`.log

# save a copy for today
mkdir -p ./figs/`date +%Y%m%d`
cp ./figs/*.png ./figs/`date +%Y%m%d`/
cp ./latest_data/* ./figs/`date +%Y%m%d`/

#!/bin/bash

# cd working directory
cd /root/lhstat

# update weather data
bash update_data_ecs.sh

# do stats ---

# AQI, dust, seeing
ipython ./lhstat/dostats.py
# SQM & wind
ipython ./lhstat/proc_sqm_wind.py

# new scripts
ipython ./lhstat/pwv_plt.py
ipython ./lhstat/read.py

echo "DONE!!!"
# save a copy for today
#mkdir -p ./figs/`date +%Y%m%d`
#cp ./figs/*.png ./figs/`date +%Y%m%d`/
#cp ./latest_data/* ./figs/`date +%Y%m%d`/
#tail ./latest_data/* > ./log/`date +%Y%m%d%H%M%S`.log

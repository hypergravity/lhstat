#!/bin/bash

# cd working directory
cd /root/lhstat

# update weather data
bash update_data_ecs.sh

# do stats
ipython dostats.py

# save a copy for today
mkdir -p ./figs/`date +%Y%m%d`
cp ./figs/*.png ./figs/`date +%Y%m%d`/
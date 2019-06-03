cd /root/lhstat/
rm ./latest_data/*
scp -p 10415 root@159.226.170.49:/data/lh/sqm/weather2019.csv ./latest_data/
scp -p 10415 root@159.226.170.49:/data/lh/sqm/SQMReadings.txt ./latest_data/

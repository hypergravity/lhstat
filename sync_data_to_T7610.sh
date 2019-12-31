#rsync -avzr --port=10415 root@101.201.57.57:/root/lhstat/latest_data /home/cham/PycharmProjects/lhstat/
#rsync -avzr --port=10284 root@ucas.china-vo.org:/root/lhstat/latest_data /home/cham/PycharmProjects/lhstat/
scp -P 10284 -r root@ucas.china-vo.org:/root/lhstat/latest_data /home/cham/PycharmProjects/lhstat/

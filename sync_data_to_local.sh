#rsync -avzr --port=10415 root@101.201.57.57:/root/lhstat/latest_data /home/cham/PycharmProjects/lhstat/
#rsync -avzr --port=10284 root@ucas.china-vo.org:/root/lhstat/latest_data /home/cham/PycharmProjects/lhstat/
#scp -P 10284 -r root@ucas.china-vo.org:/root/lhstat/latest_data /home/cham/PycharmProjects/lhstat/
#rsync -e 'ssh -p 10284' -avzr root@ucas.china-vo.org:/root/lhstat/latest_data ~/PycharmProjects/lhstat/
rsync -e 'ssh -p 12898' -avzr root@cloudvm.china-vo.org:/root/lhstat/latest_data ~/PycharmProjects/lhstat/

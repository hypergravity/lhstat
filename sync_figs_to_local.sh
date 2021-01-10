#scp -P 10415 -r root@101.201.57.57:/root/lhstat/figs/* /home/cham/PycharmProjects/lhstat/figs/
#rsync -e 'ssh -p 10284' -avzr root@ucas.china-vo.org:~/lhstat/figs/* ~/PycharmProjects/lhstat/figs/
rsync -e 'ssh -p 12898' -avzr root@cloudvm.china-vo.org:~/lhstat/figs/* ~/PycharmProjects/lhstat/figs/

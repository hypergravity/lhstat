## lhstat

冷湖监测数据网址: [http://lenghu.china-vo.org/](http://lenghu.china-vo.org/)

下载运行方式：
```commandline
git clone https://github.com/hypergravity/lhstat.git
cd lhstat
sh allin1.sh
```

## 模块说明

- do_stats.py:  处理AQI，dust，seeing数据
- grabcam.py:  读取原晨昏蒙影数据文件，现已不用
- proc_sqm_wind.py: 处理SQM(`plot_sky_brightness`, `plot_sky_goodness`)，wind(`plot_wind_wind`)数据 
- twilight.py: 计算晨昏蒙影时间


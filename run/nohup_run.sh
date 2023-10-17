#/bin/bash

# option="fidx"
# option="qidxON"
option="q1"
# option="f"
datasets="DBLP_2"
input_dir="/mnt/data/wangshu/scan/$datasets"
output_path="/mnt/data/wangshu/scan/$datasets"
mu_values=(2 5 10 15 20)
e_values=(0.2 0.4 0.6 0.8)
mu=2
e=0.2
query_file_name="query_1.txt"
query_file="./$query_file_name"
# log file path
log_file="../debug_result/$datasets-${query_file_name%.txt}-$option.log"

nohup ./cs_hin_scan -$option $input_dir $query_file \
    > $log_file 2>&1 &

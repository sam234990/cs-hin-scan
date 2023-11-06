#/bin/bash

option="fidx"
# option="qidx"
# option="q"
# option="f"
# /mnt/data/wangshu/scan/amazon/
# datasets="DBLP_2"
datasets="amazon"
input_dir="/mnt/data/wangshu/scan/$datasets"
output_path="/mnt/data/wangshu/scan/$datasets"
mu_values=(2 5 10 15 20)
e_values=(0.2 0.4 0.6 0.8)
mu=2
e=0.2
query_file_name="query_1.txt"
query_file="./$query_file_name"
# log file path
# log_file="../debug_result/$datasets-${query_file_name%.txt}-$option.log"
log_file="../debug_result/fidx/$datasets-${query_file_name%.txt}-3-$option.log"


nohup ./cs_hin_scan -$option $input_dir $query_file 3 \
    > $log_file 2>&1 &

# k_values=(3 4 5 6 7 8 9 10)
# for k in "${k_values[@]}"; do
#     log_file="../debug_result/fidx/$datasets-${query_file_name%.txt}-$k-$option.log"

#     nohup ./cs_hin_scan -$option $input_dir $query_file $k \
#         > $log_file 2>&1 &
# done



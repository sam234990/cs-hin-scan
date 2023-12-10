#/bin/bash

option="fidx1"
# option="qidx"
# option="fidxmutree"
# option="qidxSCAN"
# option="qSCAN"
# option="f"
# /mnt/data/wangshu/scan/amazon/
# datasets="amazon"
# datasets="DBLP_2"
datasets="DBPedia"
# datasets="IMDB_2"
# datasets="foursquare"
input_dir="/mnt/data/wangshu/scan/$datasets"
output_path="/mnt/data/wangshu/scan/$datasets"
# query_file_name="foursquare_query_1.txt"
query_file_name="query_DBP.txt"
# query_file_name="query_1.txt"
# query_file_name="query_type0.txt"
query_file="./$query_file_name"

# 1. query run
# log_file="../debug_result/$datasets-${query_file_name%.txt}-$option.log"

# nohup ./cs_hin_scan -$option $input_dir $query_file \
#     > $log_file 2>&1 &

# 2. index run
start_k="3"
log_file="../debug_result/fidx/$datasets-${query_file_name%.txt}-$start_k-$option.log"


nohup ./cs_hin_scan -$option $input_dir $query_file $start_k \
    > $log_file 2>&1 &


# k_values=(7 9 11 13 15)
# for k in "${k_values[@]}"; do
#     log_file="../debug_result/fidx/$datasets-${query_file_name%.txt}-$k-$option.log"

#     nohup ./cs_hin_scan -$option $input_dir $query_file $k \
#         > $log_file 2>&1 &
# done



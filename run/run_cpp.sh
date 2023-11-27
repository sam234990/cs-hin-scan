#!/bin/bash

# option="fidx"
option="f"
# option="q3"
# option="qidx"
datasets="IMDB_2"
input_dir="/mnt/data/wangshu/scan/$datasets"
index_output_path="/mnt/data/wangshu/scan/$datasets"
query_file="/home/wangshu/graphs_a/hin_scan/run/query_1.txt"
result_output_path="/home/wangshu/graphs_a/hin_scan/scan_result"

./cs_hin_scan -$option $input_dir $input_dir \
    # -o $result_output_path
    # $index_output_path

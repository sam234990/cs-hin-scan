#!/bin/bash

# option="fidx"
# option="f"
option="q"
# option="qidx"
datasets="DBPedia"
input_dir="/mnt/data/wangshu/scan/$datasets"
index_output_path="/mnt/data/wangshu/scan/$datasets"
query_file_name="query_DBP.txt"
query_file="./$query_file_name"
result_output_path="/home/wangshu/graphs_a/hin_scan/scan_result"

./cs_hin_scan -$option $input_dir $input_dir \
    # -o $result_output_path
    # $index_output_path

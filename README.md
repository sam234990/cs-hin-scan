# SSCS & SSCD

This repository contains code and data processing for the paper "Searching and Detecting Structurally Similar Communities in Large Heterogeneous Information Networks"

## Compile the code
```sh
$ cd run
$ make clean
$ make
```

## process the graph
```sh
$ cd run
$ ./hin_sscs -f {input_directory} {outpout_directory}
```


## run the code
```sh
$ cd run
$ ./hin_sscs {option} {input_directory} {query_file}
```

For example,
```sh
$ cd run
$ ./hin_sscs -q ../dataset/amazon ./query_1
```

## Usage Instructions
```sh
$ ./hin_sscs --help
# SSCS & SSCD

This repository contains code and data processing for the paper "Searching and Detecting Structurally Similar Communities in Large Heterogeneous Information Networks"

## Compile the code
```sh
$ cd run
$ make clean
$ make
```

## Process the graph
```sh
$ cd run
$ ./hin_sscs -f {input_directory} {outpout_directory}
```


## Run the code
```sh
$ cd run
$ ./hin_sscs {option} {input_directory} {query_file}
```

For example,
```sh
$ cd run
$ ./hin_sscs -q ../dataset/amazon ./query_1
```

## Usage instructions
```sh
$ ./hin_sscs --help
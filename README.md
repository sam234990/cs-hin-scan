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

For example, SSCS on ramdom 100 verices on amazon:

```sh
$ cd run                                    
$ ./hin_sscs -f ../dataset/amazon ../dataset/amazon  # first, process the graph
$ ./hin_sscs -q ../dataset/amazon ./query_1.txt  # then, SSCS on amazon
```

## Usage instructions
To run the code, please refer to the help: 
```sh
$ ./hin_sscs --help
```

For query file format, please refer to the help_query: 
```sh
$ ./hin_sscs -help_query
```

## Data format
Please refer to README.txt in /dataset/amazon
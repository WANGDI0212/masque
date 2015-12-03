#!/bin/bash


#Check arguments
if [ $# -ne 3  ]
then
    echo "$0 <reads_dir> <raw_reads_dir> <reads_dir> <output_file>"
    exit
fi



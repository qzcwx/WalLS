#!/bin/bash

for i in `cat machines.lst`
do
    echo ${i} 
    ssh ${i} mpstat 1 2 | tail -1 | awk '{ print $NF }'
done
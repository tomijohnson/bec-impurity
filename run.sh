#!/bin/bash

for i in `seq $1 $2`;
do
./load.sh isingpaired_evol_init.mat ./data/TomiCalcs $i 1
sleep 0.3
done
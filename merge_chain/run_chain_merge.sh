#!/bin/bash

python chain_merge.py -s 360 -n 10 &> merge_0.log &
pross_id_0=$!

python chain_merge.py -s 370 -n 10 &> merge_1.log &
pross_id_1=$!

python chain_merge.py -s 380 -n 4 &> merge_2.log &
pross_id_2=$!


wait $pross_id_0
wait $pross_id_1
wait $pross_id_2


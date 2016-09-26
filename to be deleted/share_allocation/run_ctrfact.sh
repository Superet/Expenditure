#!/bin/bash

cd /home/brgordon/ccv103/Exercise/run/estrun_8
NOW=$(date +"%m-%d-%Y")
nohup Rscript 3_ctrfact_sim_seq.r > ctrfact_seq1_$NOW.out 2>&1 "seg_id = 1" &
nohup Rscript 3_ctrfact_sim_seq.r > ctrfact_seq2_$NOW.out 2>&1 "seg_id = 2" &
nohup Rscript 3_ctrfact_sim_seq.r > ctrfact_seq3_$NOW.out 2>&1 "seg_id = 3" &

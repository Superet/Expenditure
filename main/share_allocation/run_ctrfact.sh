#!/bin/bash

cd /home/brgordon/ccv103/Exercise/run/estrun_4
NOW=$(date +"%m-%d-%Y")
nohup Rscript 3_ctrfact_sim_newformat.r > ctrfact_newf1_$NOW.out 2>&1 "seg_id = 1" &
nohup Rscript 3_ctrfact_sim_newformat.r > ctrfact_newf2_$NOW.out 2>&1 "seg_id = 2" &
nohup Rscript 3_ctrfact_sim_newformat.r > ctrfact_newf3_$NOW.out 2>&1 "seg_id = 3" &

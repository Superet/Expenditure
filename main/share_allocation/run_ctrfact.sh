#!/bin/bash

cd /home/brgordon/ccv103/Exercise/run/estrun_2
NOW=$(date +"%m-%d-%Y")
nohup Rscript 3_ctrfact_sim_all.r > ctrfact_sim1_$NOW.out 2>&1 "seg_id = 1" &
nohup Rscript 3_ctrfact_sim_all.r > ctrfact_sim2_$NOW.out 2>&1 "seg_id = 2" &
nohup Rscript 3_ctrfact_sim_all.r > ctrfact_sim3_$NOW.out 2>&1 "seg_id = 3" &

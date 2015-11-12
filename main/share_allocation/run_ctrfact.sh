#!/bin/bash

cd /home/brgordon/ccv103/Exercise/run/estrun_1
nohup Rscript 3_ctrfact_sim_all.r > ctrfact_sim1.out 2>&1 "seg_id = 1" &
nohup Rscript 3_ctrfact_sim_all.r > ctrfact_sim2.out 2>&1 "seg_id = 2" &
nohup Rscript 3_ctrfact_sim_all.r > ctrfact_sim3.out 2>&1 "seg_id = 3" &

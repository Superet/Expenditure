#!/bin/bash

cd /home/brgordon/ccv103/Exercise/run/estrun_11
nohup Rscript 2_mdcev_share_estimation.r > mdcev_est1.out 2>&1 "seg_id = 1" &
nohup Rscript 2_mdcev_share_estimation.r > mdcev_est2.out 2>&1 "seg_id = 2" &
nohup Rscript 2_mdcev_share_estimation.r > mdcev_est3.out 2>&1 "seg_id = 3" &

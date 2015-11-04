#!/bin/bash

cd /home/brgordon/ccv103/Exercise/run/run_5
nohup Rscript 1_nielsen_7_marginal.r > 7_effect_seg1.out 2>&1 "run_id=5" "seg_id=1" "vid=1" &
nohup Rscript 1_nielsen_7_marginal.r > 7_effect_seg2.out 2>&1 "run_id=5" "seg_id=2" "vid=1" &
nohup Rscript 1_nielsen_7_marginal.r > 7_effect_seg3.out 2>&1 "run_id=5" "seg_id=3" "vid=1" &
nohup Rscript 1_nielsen_7_marginal.r > 7_effect_seg4.out 2>&1 "run_id=5" "seg_id=4" "vid=1" &
nohup Rscript 1_nielsen_7_marginal.r > 7_effect_seg5.out 2>&1 "run_id=5" "seg_id=5" "vid=1" &
nohup Rscript 1_nielsen_7_marginal.r > 7_effect_seg6.out 2>&1 "run_id=5" "seg_id=6" "vid=1" &

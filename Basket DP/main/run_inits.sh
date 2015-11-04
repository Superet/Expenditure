#!/bin/bash

cd /home/brgordon/ccv103/Exercise/run/run_5
nohup Rscript 1_nielsen_6_DP_findinits_v3.r > inits_seg1.out 2>&1 "run_id=5" "seg_id=1" "vid=2" &
nohup Rscript 1_nielsen_6_DP_findinits_v3.r > inits_seg2.out 2>&1 "run_id=5" "seg_id=2" "vid=2" &
nohup Rscript 1_nielsen_6_DP_findinits_v3.r > inits_seg3.out 2>&1 "run_id=5" "seg_id=3" "vid=2" &
nohup Rscript 1_nielsen_6_DP_findinits_v3.r > inits_seg4.out 2>&1 "run_id=5" "seg_id=4" "vid=2" &
nohup Rscript 1_nielsen_6_DP_findinits_v3.r > inits_seg5.out 2>&1 "run_id=5" "seg_id=5" "vid=2" &
nohup Rscript 1_nielsen_6_DP_findinits_v3.r > inits_seg6.out 2>&1 "run_id=5" "seg_id=6" "vid=2" &


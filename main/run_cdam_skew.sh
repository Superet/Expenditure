#!/bin/bash

cd /kellogg/users/marketing/2661703/Exercise/run/run_2
nohup Rscript 1_nielsen_6_DP.r > 6_est_seg2_v1.out 2>&1 "run_id=2" "seg_id=2" "vid=1" &
nohup Rscript 1_nielsen_6_DP.r > 6_est_seg6_v1.out 2>&1 "run_id=2" "seg_id=6" "vid=1" &
nohup Rscript 1_nielsen_6_DP.r > 6_est_seg10_v1.out 2>&1 "run_id=2" "seg_id=10" "vid=1" &

nohup Rscript 1_nielsen_6_DP.r > 6_est_seg3_v1.out 2>&1 "run_id=2" "seg_id=3" "vid=1" &
nohup Rscript 1_nielsen_6_DP.r > 6_est_seg7_v1.out 2>&1 "run_id=2" "seg_id=7" "vid=1" &
nohup Rscript 1_nielsen_6_DP.r > 6_est_seg11_v1.out 2>&1 "run_id=2" "seg_id=11" "vid=1" &

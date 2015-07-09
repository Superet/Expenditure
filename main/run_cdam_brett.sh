#!/bin/bash

cd /home/brgordon/ccv103/Exercise/run/run_1
nohup Rscript 1_nielsen_6_DP_v2.r > 6_est_seg1_v2_chain1.out 2>&1 "arr_idx = 1" &
nohup Rscript 1_nielsen_6_DP_v2.r > 6_est_seg1_v2_chain2.out 2>&1 "arr_idx = 2" &
nohup Rscript 1_nielsen_6_DP_v2.r > 6_est_seg1_v2_chain3.out 2>&1 "arr_idx = 3" &
nohup Rscript 1_nielsen_6_DP_v2.r > 6_est_seg1_v2_chain4.out 2>&1 "arr_idx = 4" &

nohup Rscript 1_nielsen_6_DP_v2.r > 6_est_seg2_v2_chain1.out 2>&1 "arr_idx = 5" &
nohup Rscript 1_nielsen_6_DP_v2.r > 6_est_seg2_v2_chain2.out 2>&1 "arr_idx = 6" &
nohup Rscript 1_nielsen_6_DP_v2.r > 6_est_seg2_v2_chain3.out 2>&1 "arr_idx = 7" &
nohup Rscript 1_nielsen_6_DP_v2.r > 6_est_seg2_v2_chain4.out 2>&1 "arr_idx = 8" &

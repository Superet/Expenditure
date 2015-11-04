#!/bin/bash

cd /kellogg/users/marketing/2661703/Exercise/run/run_5
nohup Rscript 1_nielsen_7_marginal.r > 7_effect_seg7.out 2>&1 "run_id=5" "seg_id=7" "vid=1" &
nohup Rscript 1_nielsen_7_marginal.r > 7_effect_seg8.out 2>&1 "run_id=5" "seg_id=8" "vid=1" &
nohup Rscript 1_nielsen_7_marginal.r > 7_effect_seg9.out 2>&1 "run_id=5" "seg_id=9" "vid=1" &
nohup Rscript 1_nielsen_7_marginal.r > 7_effect_seg10.out 2>&1 "run_id=5" "seg_id=10" "vid=1" &
nohup Rscript 1_nielsen_7_marginal.r > 7_effect_seg11.out 2>&1 "run_id=5" "seg_id=11" "vid=1" &
nohup Rscript 1_nielsen_7_marginal.r > 7_effect_seg12.out 2>&1 "run_id=5" "seg_id=12" "vid=1" &

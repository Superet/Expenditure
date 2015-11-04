#!/bin/bash

cd /kellogg/users/marketing/2661703/Exercise/run/estrun_1
nohup Rscript 1_nielsen_6_counterfactural_v2.r > seg1.out 2>&1 "seg_id=1" &
nohup Rscript 1_nielsen_6_counterfactural_v2.r > seg2.out 2>&1 "seg_id=2" &
nohup Rscript 1_nielsen_6_counterfactural_v2.r > seg3.out 2>&1 "seg_id=3" &

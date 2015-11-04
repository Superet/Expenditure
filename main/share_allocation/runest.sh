#!/bin/bash

cd /home/brgordon/ccv103/Exercise/run/estrun_1
nohup Rscript 2_mdcev_share_inclusive_value.r > inclusive1.out 2>&1 "seg_id = 1" &
nohup Rscript 2_mdcev_share_inclusive_value.r > inclusive2.out 2>&1 "seg_id = 2" &
nohup Rscript 2_mdcev_share_inclusive_value.r > inclusive3.out 2>&1 "seg_id = 3" &

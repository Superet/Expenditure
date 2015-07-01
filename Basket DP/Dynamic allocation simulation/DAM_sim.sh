#!/bin/bash

cd /home/brgordon/ccv103/processed\ data/run_sim
nohup Rscript DAM_3.R > DAM_w0t0.out 2>&1 'w0t0'&
nohup Rscript DAM_3.R > DAM_w0t1.out 2>&1 'w0t1'&
nohup Rscript DAM_3.R > DAM_w0t2.out 2>&1 'w0t2'&

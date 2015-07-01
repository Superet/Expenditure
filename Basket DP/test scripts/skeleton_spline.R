rm(list=ls())
library(Rcpp)
library(RcppGSL)
library(devtools)

mycpp	<- "~/Documents/Research/Store switching/Exercise/Basket DP/test scripts/skeleton_spline.cpp"
sourceCpp(mycpp)


Rcpp.package.skeleton(name="myspl", attributes=TRUE, list=c("splineC"), 
cpp_files=mycpp, example_code=FALSE)
compileAttributes("myspl")

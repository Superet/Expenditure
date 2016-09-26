# This script estimates model. The moduels follow as: 
# 1. Segment households based on initial income level and subset data; 
# 2. Prepare estimation data; 
# 3. Estimate the conditional allocation model at different initial values; 
# 4. Simulate inclusive values; 
# 5. Estimate upper level expenditure model. 

cat("This program begins to run at", as.character(Sys.time()), ".\n")

library(ggplot2)
library(reshape2)
library(Rcpp)
library(RcppArmadillo)
library(maxLik)
library(evd)
library(data.table)
library(doParallel)
library(foreach)
library(nloptr)
library(mgcv)
options(error = quote({dump.frames(to.file = TRUE)}))

# setwd("~/Documents/Research/Store switching/Processed_data")
# plot.wd	<- '~/Desktop'
# sourceCpp(paste("../Exercise/Multiple_discrete_continuous_model/", model_name, ".cpp", sep=""))
# source("../Exercise/Multiple_discrete_continuous_model/0_Efficient_Allocation_function.R")
# source("../Exercise/main/ctrfact_sim_functions.r")

# setwd("/home/brgordon/ccv103/Exercise/run")
# setwd("/kellogg/users/marketing/2661703/Exercise/run")
setwd("/sscc/home/c/ccv103/Exercise/run")
# setwd("U:/Users/ccv103/Documents/Research/Store switching/run")

########################
# Read estimation data # 
########################
run_id		<- 1
ver.date	<- "2016-09-11"
nseg 		<- 9 
ehat1.ls	<- vector("list", length = nseg)
ehat.lm.ls	<- vector("list", length = nseg)
fit.ls		<- vector("list", length = nseg)
sim.unq.ful	<- data.frame()
tmp.ls		<- ls()

for(ii in 1:nseg){
	fname	<- paste("estrun_",run_id,"/comp_MDCEV_est_seg", ii, "_", ver.date,".rdata",sep="")
	cat(fname,"\n")
	load(fname)
	
	ehat1.ls[[ii]]		<- ehat1
	ehat.lm.ls[[ii]]	<- ehat.lm
	fit.ls[[ii]]		<- c(OneStage.ll = sol.n$maximum, lmAIC = sum(sapply(lm.fit, AIC)))
	sim.unq.ful			<- rbind(sim.unq.ful, cbind(segment = ii, sim.unq))
	
	rm(list = setdiff(ls(), c("ii", "tmp.ls",tmp.ls)))
}

load("estrun_1/marginal_effect_sim500_2016-09-12.rdata")

identical(sim.unq$household_code, sim.unq.ful$household_code)
selyr	<- 2007
sel		<- mydata.full$year == selyr & !duplicated(mydata.full[,c("household_code", "year")]) # & mydata.full$household_code < 2250000
identical(mydata.full[sel,"household_code"], sim.unq.ful$household_code)

###########################
# Summarize model fitness #
###########################
e.true	<- as.matrix(mydata.full[sel,paste("DOL_", gsub("\\s", "_", fmt_name),sep="")])
esign.true	<- 1*(e.true > 0)
ehat.ls	<- list(Propose = do.call(rbind, lapply(sim.base07, function(x) x$Average[,fmt_name])), 
				Simultaenous = do.call(rbind, lapply(ehat1.ls, function(x) apply(x, c(2,3), mean, na.rm=T)))[,fmt_name], 
				Regression = do.call(rbind, ehat.lm.ls) )
esignhat.ls	<- lapply(ehat.ls, function(x) 1*(x!=0))
sel		<- esign.true == 0
tmp.tab <- cbind(mse = sapply(ehat.ls, function(x) mean((x-e.true)^2)), 
			 	hit.positive = sapply(esignhat.ls, function(x) sum((x==esign.true)[!sel])/sum(!sel) ), 
			 	hit.zero = sapply(esignhat.ls, function(x) sum((x==esign.true)[sel])/sum(sel) ) 	)
cat("Model fitness:\n"); print(tmp.tab); cat("\n")





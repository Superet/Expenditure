library(ggplot2)
library(reshape2)
library(Rcpp)
library(RcppArmadillo)
library(maxLik)
library(evd)
library(data.table)
library(doParallel)
library(foreach)
# library(chebpol)
library(nloptr)
library(mgcv)
options(error = quote({dump.frames(to.file = TRUE)}))

args <- commandArgs(trailingOnly = TRUE)
print(args)
if(length(args)>0){
    for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
    }
}

# setwd("~/Documents/Research/Store switching/processed data")
# plot.wd	<- '~/Desktop'
# source("../Exercise/Multiple_discrete_continuous_model/0_Allocation_function.R")
# source("../Exercise/main/share_allocation/ctrfact_sim_functions.r")

# setwd("/home/brgordon/ccv103/Exercise/run")
# setwd("/kellogg/users/marketing/2661703/Exercise/run")
setwd("/sscc/home/c/ccv103/Exercise/run")
run_id		<- 4
make_plot	<- TRUE
ww			<- 10
ar			<- .6

###########################
# Read counterfactal data # 
###########################
# Read estimation from all segments #
nseg	<- 3
sim.ls	<- setNames(vector("list", length = nseg), c("Low", "Med", "High"))
my.simdata	<- data.frame()
numsim		<- 1000
ver.date1	<- "2016-03-28"
cpi.adj		<- TRUE
annual		<- TRUE
annual.week	<- 26
tmpls	<- ls()

for(ii in 1:nseg){
	fname			<- paste("ctrfact_ref_seg",ii,sep="")
	loadf			<- paste("estrun_",run_id, "/",fname, "_sim", numsim, "_", ver.date1, ".rdata", sep="")
	load(loadf)
	
	# Convert biweekly simulations to annual simulations
	if(annual){
		sim.base07$Allocation	<- sim.base07$Allocation * annual.week
		sim.base08$Allocation	<- sim.base08$Allocation * annual.week
		for(j in 1:length(ref.sim)){
			ref.sim[[j]]$Allocation 		<- ref.sim[[j]]$Allocation * annual.week
		}
	}
	
	tmp	<- lapply(ref.sim, function(x) x$Allocation)
	sim.ls[[ii]]<- c(list(Base07 = sim.base07$Allocation, Base08 = sim.base08$Allocation), tmp)
	my.simdata	<- rbind(my.simdata, data.frame(sim.data, IncGrp = names(sim.ls)[ii]))
	
	# Delete unused objects
	rm(list = setdiff(ls(), c(tmpls, "tmpls", "ii", "fmt_name", "R", "ret.a", "ret.b")))
	print(ii)
}

# Set small simulation value to zero 
trim.alpha		<- 0
zero.max		<- .01	# 1 cent
for(i in 1:nseg){
	for(j in 1:length(sim.ls[[1]])){
		sim.ls[[i]][[j]]	<- 1*(sim.ls[[i]][[j]] >= zero.max) * sim.ls[[i]][[j]]
	}
}

#############
# Functions #
#############
weight_summary	<- function(x, freq, fun, ...){
	x.full	<- unlist(lapply(1:length(x), function(i) rep(x[i], freq[i])))
	out		<- fun(x.full, ...)
	return(out)
}

mean_qt <- function(x, alpha = .25){
	out	<- c(mean(x), quantile(x, c(alpha, 1-alpha)))
	names(out)	<- c("y", "ymin", "ymax")
	return(out)
}

weight_se	<- function(se, freq){
# This function computes the standard error of x.bar = sum_i (freq_i*x_i)/sum_i freq_i
	wt	<- freq/sum(freq)
	v	<- sum(wt^2 * se^2, na.rm =T)
	return(sqrt(v))
}

mytrimfun	<- function (x, alpha = 0.05){
    if (alpha == 0) {
        return(x)
    }else {
        n <- length(x)
        lo <- floor(n * alpha) + 1
        hi <- n + 1 - lo
		if(sum(is.na(x)) >= (lo-1)){
			out	<- x
		}else{
			out <- sort.int(x, partial = unique(c(lo, hi)), na.last = FALSE)[lo:hi]
		}      
        return(out)
    }
}

#################
# Income effect #
#################
# Income frequency
tab.inc	<- data.table(my.simdata)
tab.inc	<- tab.inc[, list(nhh = length(unique(household_code))), by = list(IncGrp, income)]
setnames(tab.inc, "income", "income2007")
setkeyv(tab.inc, c("IncGrp", "income2007"))
tab.inc$income2007	<- round(tab.inc$income2007, 2)

# Summarize the the allocation at each income level within each segment
sel.base	<- "Base08"
tmp	<- setNames(vector("list", length(sim.ls)), names(sim.ls))
for(i in 1:nseg){
	tmp[[i]]	<- setNames(vector("list", length(ret.b)), ret.b)
	for(j in ret.b){
		tmp[[i]][[j]]	<- apply(sim.ls[[i]][[j]] - sim.ls[[i]][[sel.base]], c(2, 3), mean, trim = trim.alpha)
	}
}
dol.dt	<- melt(tmp)
names(dol.dt)	<- c("lnInc", "retailer", "Difference", "retail.b", "IncGrp")
dol.dt$retailer	<- factor(fmt_name[dol.dt$retailer], levels = fmt_name)
dol.dt$income2007	<- round(exp(dol.dt$lnInc)/.9, 2)
dim(dol.dt)
dol.dt	<- merge(dol.dt, tab.inc, by = c("IncGrp", "income2007"), all.x = TRUE)
dim(dol.dt)
dol.dt	<- data.table(dol.dt)

# Summarize the overall expenditure change 
tmp.tab	<- dol.dt[, list(Difference = weight_summary(Difference, nhh, mean), sd = weight_summary(Difference, nhh, sd)), 
					by = list(retail.b, retailer)]
tmp.tab1<- dol.dt[, list(Difference = weight_summary(Difference, nhh, mean), sd = weight_summary(Difference, nhh, sd)), 
					by = list(retail.b, retailer, IncGrp)]
setkeyv(tmp.tab, c("retail.b", "retailer"))
setkeyv(tmp.tab1, c("retail.b", "IncGrp", "retailer"))					
for(i in ret.b){
	cat("---------------------------------------------------------------------------\n")
	cat("Overall annual expenditure change of ", ret.a,"reformat into", i, ":\n"); print(tmp.tab[i,]); cat("\n")
	cat("Annual expenditure change of ", ret.a,"reformat into", i, " by income group:\n"); print(tmp.tab1[i,]); cat("\n")		
}

# -------------------------#
# Compute share difference #
# Convert dollar expenditure to share  
sim.shr	<- sim.ls
for(i in 1:nseg){
	for(j in 1:length(sim.ls[[1]])){
		sim.shr[[i]][[j]]	<- apply(sim.shr[[i]][[j]], c(1,2), function(x) x/sum(x)*100 )
	}
}

tmp	<- setNames(vector("list", length(sim.shr)), names(sim.shr))
for(i in 1:nseg){
	tmp[[i]]	<- setNames(vector("list", length(ret.b)), ret.b)
	for(j in ret.b){
		tmp[[i]][[j]]	<- apply(sim.shr[[i]][[j]] - sim.shr[[i]][[sel.base]], c(1, 3), mean, trim = trim.alpha)
	}
}
shr.dt	<- melt(tmp)
names(shr.dt)	<- c("retailer", "lnInc", "Difference", "retail.b", "IncGrp")
shr.dt$retailer	<- factor(fmt_name[shr.dt$retailer], levels = fmt_name)
shr.dt$income2007	<- round(exp(shr.dt$lnInc)/.9, 2)
dim(shr.dt)
shr.dt	<- merge(shr.dt, tab.inc, by = c("IncGrp", "income2007"), all.x = TRUE)
dim(shr.dt)
shr.dt	<- data.table(shr.dt)

# Summarize the overall expenditure change 
tmp.tab	<- shr.dt[, list(Difference = weight_summary(Difference, nhh, mean), sd = weight_summary(Difference, nhh, sd)), 
					by = list(retail.b, retailer)]
tmp.tab1<- shr.dt[, list(Difference = weight_summary(Difference, nhh, mean), sd = weight_summary(Difference, nhh, sd)), 
					by = list(retail.b, retailer, IncGrp)]
setkeyv(tmp.tab, c("retail.b", "retailer"))
setkeyv(tmp.tab1, c("retail.b", "IncGrp", "retailer"))					
for(i in ret.b){
	cat("---------------------------------------------------------------------------\n")
	cat("Overall annual expenditure change of ", ret.a,"reformat into", i, ":\n"); print(tmp.tab[i,]); cat("\n")
	cat("Annual expenditure change of ", ret.a,"reformat into", i, " by income group:\n"); print(tmp.tab1[i,]); cat("\n")		
}

library(ggplot2)
library(reshape2)
library(maxLik)
library(evd)
library(data.table)
library(doParallel)
library(foreach)
# library(chebpol)
library(nloptr)
library(mgcv)
options(error = quote({dump.frames(to.file = TRUE)}))

# seg_id	<- as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))

args <- commandArgs(trailingOnly = TRUE)
print(args)
if(length(args)>0){
    for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
    }
}


model_name 	<- "MDCEV_share"
# setwd("~/Documents/Research/Store switching/Processed_data/processed data_20160809")
# plot.wd	<- '~/Desktop'
# sourceCpp(paste("../../Exercise/Multiple_discrete_continuous_model/", model_name, ".cpp", sep=""))
# source("../../Exercise/Multiple_discrete_continuous_model/0_Efficient_Allocation_function.R")
# source("../../Exercise/main/ctrfact_sim_functions.r")

# setwd("/home/brgordon/ccv103/Exercise/run")
# setwd("/kellogg/users/marketing/2661703/Exercise/run")
# setwd("/sscc/home/c/ccv103/Exercise/run")
setwd("U:/Users/ccv103/Documents/Research/Store switching/run")

run_id		<- 1
plot.wd		<- getwd()
make_plot	<- TRUE
ww			<- 10
ar			<- .6

source("0_Efficient_Allocation_function.R")
source("ctrfact_sim_functions.r")

# Load estimation data 
ver.date	<- "2016-08-22"
cpi.adj		<- TRUE

(loadf	<-  paste("estrun_",run_id,"/MDCEV_annual_est_", ver.date,".rdata",sep=""))
load(loadf)
rm(list = intersect(ls(), c("gamfit", "shr","model_name", "tmpdat")))

# Set simulation parameters
week.price		<- FALSE
interp.method	<- "spline"					# Spline interpolation or Chebyshev interpolation
exp.method		<- "Utility"				# Utility maximization or finding roots for first order condition
trim.alpha		<- 0
numsim			<- 500
draw.par		<- FALSE
sim.omega		<- FALSE
# fname			<- paste("marginal_annual_sim_", numsim, "_", as.character(Sys.Date()), sep="")
fname			<- paste("marginal_", seldpt, "_sim_", numsim, "_", as.character(Sys.Date()), sep="")
cat("Output file name is", fname, ".\n")

###############################
# Prepare simulation elements # 
###############################
# Data required: parameters, income level, price, retail attributes, and random draws
# For each simulation scenario, if income, price or retail attributes change, we need to re-simulate inclusive values.

# Set simulation parameters
lambda 	<- coef(lsol)	#coef(sol.top2)
shr.par	<- coef(sol)
# ret.idx	<- 26:30		# Column index flagging the retail attribugtes in X1

#-----------------------#
# Construct income data # 
# Take the households' income in 2007 as basis 
selyr	<- 2007
sel		<- mydata$year == selyr & !duplicated(mydata[,c("household_code", "year")])# & mydata$household_code < 2250000
sim.unq	<- cbind(household_code = mydata[sel,c("household_code")], Z[sel,])
sim.unq$price	<- as.matrix(mydata[sel,paste("PRC_", gsub("\\s", "_", fmt_name), sep="")])
tmp		<- do.call(cbind, lapply(X1, function(x) x[sel,ret.idx]))
sim.unq$X<- as.matrix(cbind(Z1[sel,setdiff(colnames(Z1), c("1","ln_inc"))], tmp))
cat("dim(sim.unq) =", dim(sim.unq), "\n")

# Retail attributes
X1.07	<- lapply(X1, function(x) x[sel,])
X2.07	<- lapply(X2, function(x) x[sel,])
X2.08	<- X1.08 	<- vector("list", length=length(fmt_name))
for(i in 1:length(fmt_name)){
	tmp0<- rep(0, R-1)
	if(i<beta0_base){
		tmp0[i]	<- 1
	}else if(i > beta0_base){
		tmp0[(i-1)]	<- 1
	}
	tmp1<- as.matrix(cbind(1, rec=1, Z[sel,]))	
	tmp2<- kronecker(tmp1, t(tmp0))
	colnames(tmp2)	<- paste(rep(c("Intercept", "Rec", demo.col), each = R-1), fmt_name[-beta0_base], sep=":")
	tmp3<- rep(0, R)
	tmp3[i]	<- 1
	tmp3<- kronecker(tmp1, t(tmp3))
	colnames(tmp3)	<- paste(rep(c("Intercept", "Rec", demo.col), each = R), fmt_name[-beta0_base], sep=":")
	sel1	<- attr_mat$channel_type == fmt_name[i]
	tmp	<- as.matrix(attr_mat[sel1,selcol][sel,])
	X1.08[[i]]	<- cbind(tmp2, tmp, tmp*1)
	X2.08[[i]]	<- cbind(tmp3, tmp, tmp*1)
	colnames(X1.08[[i]])	<- c(colnames(tmp2), colnames(tmp), paste(colnames(tmp), ":Rec", sep=""))
	colnames(X2.08[[i]])	<- c(colnames(tmp3), colnames(tmp), paste(colnames(tmp), ":Rec", sep=""))
}

#-------------------#
# Take random draws #
set.seed(666)
eps_draw	<- matrix(rgev(numsim*R, scale = exp(shr.par["ln_sigma"])), numsim, R)
par_draw <- NULL

##############
# Simulation #
##############
# Register parallel computing
mycore 	<- 3
# cl		<- makeCluster(mycore, type = "FORK")
cl		<- makeCluster(mycore, type = "PSOCK")			# Register cores in Windows
registerDoParallel(cl)
cat("Register", mycore, "core parallel computing. \n")

# ------------------ #
# Basline simulation #
# Simulate expenditure and expenditure share. 
pct				<- proc.time()
sim.base07		<- SimWrapper_fn(omega_deriv, lambda, shr.par, base = beta0_base, 
								X1 = X1.07, X2=X2.07, price = sim.unq$price, 
								demo_mat = as.matrix(cbind(Intercept = 1, Rec = 0, sim.unq[,demo.col])), demo.x = TRUE, 
								eps_draw = eps_draw, ret.idx = ret.idx, sim.y = NULL, ret.sim = TRUE)
use.time		<- proc.time() - pct						
cat("2007 Baseline simulation finishes with", use.time[3]/60, "min.\n")

pct				<- proc.time()
sim.base08		<- SimWrapper_fn(omega_deriv, lambda, shr.par, base = beta0_base, 
								X1 = X1.08, X2=X2.08, price = sim.unq$price, 
								demo_mat = as.matrix(cbind(Intercept = 1, Rec = 1, sim.unq[,demo.col])), demo.x = TRUE, 
								eps_draw = eps_draw, ret.idx = ret.idx, sim.y = NULL, ret.sim = TRUE)
use.time		<- proc.time() - pct						
cat("2008 Baseline simulation finishes with", use.time[3]/60, "min.\n")

iv.change <- .1
(sel.retail	<- which(fmt_name == "Discount Store"))

#-------------------------------- #
# Simulation of fixed expenditure #
pct			<- proc.time()
sim.ls07	<- foreach(v = c(selcol, "price"), .errorhandling = "remove", .packages=c("mgcv", "nloptr")) %dopar% {
	X1_new 		<- X1.07
	X2_new		<- X2.07
	price_new	<- sim.unq$price
	if(v == "price"){
		price_new[,sel.retail]		<- price_new[,sel.retail]*(1 + iv.change)
	}else if(v %in% c("ln_upc_per_mod", "ln_num_module")){
		X1_new[[sel.retail]][,v]	<- X1_new[[sel.retail]][,v] + log(1 + iv.change)
		X2_new[[sel.retail]][,v]	<- X2_new[[sel.retail]][,v] + log(1 + iv.change)
	}else{
		X1_new[[sel.retail]][,v]	<- X1_new[[sel.retail]][,v] * (1 + iv.change)
		X2_new[[sel.retail]][,v]	<- X2_new[[sel.retail]][,v] * (1 + iv.change)
	}
	
	out1	<- SimWrapper_fn(omega_deriv, lambda, shr.par, base = beta0_base, 
							X1 = X1_new, X2=X2_new, price = price_new, 
							demo_mat = as.matrix(cbind(Intercept = 1, Rec = 0, sim.unq[,demo.col])), demo.x = TRUE, 
							eps_draw = eps_draw, ret.idx = ret.idx, sim.y = sim.base07$y, ret.sim = TRUE)
	out1.avg	<- data.frame(out1$Average, retailer = fmt_name[sel.retail], Var = v, change = iv.change )
	out1$Average<- out1.avg
	return(out1)
}
use.time	<- proc.time() - pct
cat("Parallel simulation of retail price finishes with", use.time[3]/60, "min.\n")

pct			<- proc.time()
sim.ls08	<- foreach(v = c(selcol, "price"), .errorhandling = "remove", .packages=c("mgcv", "nloptr")) %dopar% {
	X1_new 		<- X1.08
	X2_new		<- X2.08
	price_new	<- sim.unq$price
	if(v == "price"){
		price_new[,sel.retail]	<- price_new[,sel.retail] * (1 + iv.change)
	}else if(v %in% c("ln_upc_per_mod", "ln_num_module")){
		X1_new[[sel.retail]][,v]	<- X1_new[[sel.retail]][,v] + log(1 + iv.change)
		X1_new[[sel.retail]][,paste(v,":Rec",sep="")]	<- X1_new[[sel.retail]][,paste(v,":Rec",sep="")] + log(1 + iv.change)
		X2_new[[sel.retail]][,v]	<- X2_new[[sel.retail]][,v] + log(1 + iv.change)
		X2_new[[sel.retail]][,paste(v,":Rec",sep="")]	<- X2_new[[sel.retail]][,paste(v,":Rec",sep="")] + log(1 + iv.change)
	}else{
		X1_new[[sel.retail]][,v]	<- X1_new[[sel.retail]][,v] * (1 + iv.change)
		X1_new[[sel.retail]][,paste(v,":Rec",sep="")]	<- X1_new[[sel.retail]][,paste(v,":Rec",sep="")] * (1 + iv.change)
		X2_new[[sel.retail]][,v]	<- X2_new[[sel.retail]][,v] * (1 + iv.change)
		X2_new[[sel.retail]][,paste(v,":Rec",sep="")]	<- X2_new[[sel.retail]][,paste(v,":Rec",sep="")] * (1 + iv.change)
	}
	
	out1	<- SimWrapper_fn(omega_deriv, lambda, shr.par, base = beta0_base, 
							X1 = X1_new, X2=X2_new, price = price_new, 
							demo_mat = as.matrix(cbind(Intercept = 1, Rec = 1, sim.unq[,demo.col])), demo.x = TRUE, 
							eps_draw = eps_draw, ret.idx = ret.idx, sim.y = sim.base08$y, ret.sim = TRUE)
	out1.avg	<- data.frame(out1$Average, retailer = fmt_name[sel.retail], Var = v, change = iv.change )
	out1$Average<- out1.avg
	return(out1)
}
use.time	<- proc.time() - pct
cat("Parallel simulation of retail price finishes with", use.time[3]/60, "min.\n")

#--------------------------------------------------------- #
# Simulation of effects on both expenditure and allocation #
pct			<- proc.time()
sim.2ls07	<- foreach(v = c(selcol, "price"), .errorhandling = "remove", .packages=c("mgcv", "nloptr")) %dopar% {
	X1_new 		<- X1.07
	X2_new		<- X2.07
	price_new	<- sim.unq$price
	if(v == "price"){
		price_new[,sel.retail]		<- price_new[,sel.retail]*(1 + iv.change)
	}else if(v %in% c("ln_upc_per_mod", "ln_num_module")){
		X1_new[[sel.retail]][,v]	<- X1_new[[sel.retail]][,v] + log(1 + iv.change)
		X2_new[[sel.retail]][,v]	<- X2_new[[sel.retail]][,v] + log(1 + iv.change)
	}else{
		X1_new[[sel.retail]][,v]	<- X1_new[[sel.retail]][,v] * (1 + iv.change)
		X2_new[[sel.retail]][,v]	<- X2_new[[sel.retail]][,v] * (1 + iv.change)
	}
	
	out1	<- SimWrapper_fn(omega_deriv, lambda, shr.par, base = beta0_base, 
							X1 = X1_new, X2=X2_new, price = price_new, 
							demo_mat = as.matrix(cbind(Intercept = 1, Rec = 0, sim.unq[,demo.col])), demo.x = TRUE, 
							eps_draw = eps_draw, ret.idx = ret.idx, sim.y = NULL, ret.sim = TRUE)
	out1.avg	<- data.frame(out1$Average, retailer = fmt_name[sel.retail], Var = v, change = iv.change )
	out1$Average<- out1.avg
	return(out1)
}
use.time	<- proc.time() - pct
cat("Parallel simulation of retail price finishes with", use.time[3]/60, "min.\n")

pct			<- proc.time()
sim.2ls08	<- foreach(v = c(selcol, "price"), .errorhandling = "remove", .packages=c("mgcv", "nloptr")) %dopar% {
	X1_new 		<- X1.08
	X2_new		<- X2.08
	price_new	<- sim.unq$price
	if(v == "price"){
		price_new[,sel.retail]	<- price_new[,sel.retail] * (1 + iv.change)
	}else if(v %in% c("ln_upc_per_mod", "ln_num_module")){
		X1_new[[sel.retail]][,v]	<- X1_new[[sel.retail]][,v] + log(1 + iv.change)
		X1_new[[sel.retail]][,paste(v,":Rec",sep="")]	<- X1_new[[sel.retail]][,paste(v,":Rec",sep="")] + log(1 + iv.change)
		X2_new[[sel.retail]][,v]	<- X2_new[[sel.retail]][,v] + log(1 + iv.change)
		X2_new[[sel.retail]][,paste(v,":Rec",sep="")]	<- X2_new[[sel.retail]][,paste(v,":Rec",sep="")] + log(1 + iv.change)
	}else{
		X1_new[[sel.retail]][,v]	<- X1_new[[sel.retail]][,v] * (1 + iv.change)
		X1_new[[sel.retail]][,paste(v,":Rec",sep="")]	<- X1_new[[sel.retail]][,paste(v,":Rec",sep="")] * (1 + iv.change)
		X2_new[[sel.retail]][,v]	<- X2_new[[sel.retail]][,v] * (1 + iv.change)
		X2_new[[sel.retail]][,paste(v,":Rec",sep="")]	<- X2_new[[sel.retail]][,paste(v,":Rec",sep="")] * (1 + iv.change)
	}
	
	out1	<- SimWrapper_fn(omega_deriv, lambda, shr.par, base = beta0_base, 
							X1 = X1_new, X2=X2_new, price = price_new, 
							demo_mat = as.matrix(cbind(Intercept = 1, Rec = 1, sim.unq[,demo.col])), demo.x = TRUE, 
							eps_draw = eps_draw, ret.idx = ret.idx, sim.y = NULL, ret.sim = TRUE)
	out1.avg	<- data.frame(out1$Average, retailer = fmt_name[sel.retail], Var = v, change = iv.change )
	out1$Average<- out1.avg
	return(out1)
}
use.time	<- proc.time() - pct
cat("Parallel simulation of retail price finishes with", use.time[3]/60, "min.\n")

stopCluster(cl)
cat("Stop clustering. \n")

# Check the marginal effects. 
selcol2	<- c(selcol, "price")
tmp0	<- colSums(sim.base07$Average[, fmt_name], na.rm=T)
tmp1	<- sapply(sim.ls07, function(x) colSums(x$Average[,gsub(" ", ".", fmt_name)], na.rm=T))
tmp		<- setNames(tmp1[sel.retail,] - tmp0[sel.retail], selcol2)
cat("Impact on exp. share in 2007:\n"); print(tmp); cat("\n")
cat("Elasticity in 2007 when expenditure is fixed:\n"); print(tmp/tmp0[sel.retail]/iv.change); cat("\n")
tmp1	<- sapply(sim.2ls07, function(x) colSums(x$Average[,gsub(" ", ".", fmt_name)], na.rm=T))
tmp		<- setNames(tmp1[sel.retail,] - tmp0[sel.retail], selcol2)
cat("Impact on exp. share in 2007:\n"); print(tmp); cat("\n")
cat("Elasticity in 2007 when expenditure is free:\n"); print(tmp/tmp0[sel.retail]/iv.change); cat("\n")

tmp0	<- colSums(sim.base08$Average[, fmt_name], na.rm=T)
tmp1	<- sapply(sim.ls08, function(x) colSums(x$Average[,gsub(" ", ".", fmt_name)], na.rm=T))
tmp		<- setNames(tmp1[sel.retail,] - tmp0[sel.retail], selcol2)
cat("Impact on exp. share in 2008:\n"); print(tmp); cat("\n")
cat("Elasticity in 2008 when expenditure is fixed:\n"); print(tmp/tmp0[sel.retail]/iv.change); cat("\n")
tmp1	<- sapply(sim.2ls08, function(x) colSums(x$Average[,gsub(" ", ".", fmt_name)], na.rm=T))
tmp		<- setNames(tmp1[sel.retail,] - tmp0[sel.retail], selcol2)
cat("Impact on exp. share in 2008:\n"); print(tmp); cat("\n")
cat("Elasticity in 2008 when expenditure is free:\n"); print(tmp/tmp0[sel.retail]/iv.change); cat("\n")

################
# Save results #
################
rm(list  = intersect(ls(), c("ar", "args", "cl", "ggtmp", "ggtmp1", "ggtmp2", "i", "lastFuncGrad", "lastFuncParam", "make_plot", "mycore", 
			"myfix", "plot.wd", "s1_index", "sel", "selyr", "price", "sol", "sol.top", "sol.top2", "tmp", "tmp_coef", 
			"tmp_price", "tmp1", "tmp2", "use.time", "ver.date", "var1", "ww", "f", 
			"numnodes", "out", "out1", "pct", "tmpd1", "tmpd2", "tmpdat", "u", "W", "y", "y.nodes",
			"Allocation_constr_fn", "Allocation_fn", "Allocation_nlop_fn", "cheb.1d.basis", "cheb.basis", "chebfun", 
			"exp_fn", "expFOC_fn", "incl_value_fn", "mysplfun", "mytrimfun", "param_assignR", "simExp_fn", "SimOmega_fn", 
			"SimWrapper_fn", "solveExp_fn", "spl2dfun", "uP_fn", "uGrad_fn")))

save.image(paste("estrun_",run_id, "/", fname, ".rdata",sep=""))

cat("This program is done.\n")

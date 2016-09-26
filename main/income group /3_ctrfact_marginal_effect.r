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
seg_id <- 1
cat("seg_id =", seg_id, ".\n")

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

run_id		<- 2
plot.wd		<- getwd()
make_plot	<- TRUE
ww			<- 10
ar			<- .6

source("0_Efficient_Allocation_function.R")
source("ctrfact_sim_functions.r")

# Load estimation data 
ver.date	<- "2016-08-28"
cpi.adj		<- TRUE

loadf	<- paste("estrun_",run_id,"/MDCEV_est_seg",seg_id,"_", ver.date,".rdata",sep="")
loadf
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
fname			<- paste("marginal_seg",seg_id,sep="")
if(draw.par){ fname <- paste(fname, "_pardraw", sep="") }
if(sim.omega){fname	<- paste(fname, "_simomega", sep="") }
fname			<- paste(fname, "_sim", numsim, "_", as.character(Sys.Date()), sep="")
cat("Output file name is", fname, ".\n")

###############################
# Prepare simulation elements # 
###############################
# Data required: parameters, income level, price, retail attributes, and random draws
# For each simulation scenario, if income, price or retail attributes change, we need to re-simulate inclusive values.

# Set simulation parameters
lambda 	<- coef(lsol)	#coef(sol.top2)
shr.par	<- coef(sol)
ret.idx	<- 11:15		# Column index flagging the retail attribugtes in X1

#-----------------------#
# Construct income data # 
selyr	<- 2007
sel		<- mydata$year == selyr & !duplicated(mydata[,c("household_code", "year")]) # & mydata$household_code < 2250000
sim.unq	<- mydata[sel,c("household_code", "ln_inc")]
sim.unq$price	<- as.matrix(mydata[sel,paste("PRC_", gsub("\\s", "_", fmt_name), sep="")])
tmp		<- do.call(cbind, lapply(X1, function(x) x[sel,ret.idx]))
sim.unq$X<- as.matrix( tmp)
cat("dim(sim.unq) =", dim(sim.unq), "\n")

# attr_mat	<- merge(mydata[,c("household_code", "year", "month", "scantrack_market_descr")], 
# 					hh_dist[,c("household_code", "year", "channel_type", "distance_haver")], 
# 				by = c("household_code", "year"), all.x=T)
# Normalize distance by zipcode
tmp		<- data.table(hh_dist[,c("household_code", "panelist_zip_code","year", "channel_type", "distance_haver")])
tmp		<- tmp[, maxd:=max(distance_haver, na.rm=T), by = list(panelist_zip_code)]
tmp		<- tmp[, distance_haver := distance_haver/maxd]
tmp		<- data.frame(tmp)
attr_mat	<- merge(mydata[,c("household_code", "year", "month", "scantrack_market_descr")], 
					tmp[,c("household_code", "year", "channel_type", "distance_haver")], 
					by = c("household_code", "year"), all.x=T)
attr_mat	<- merge(attr_mat, fmt_attr[,c("scantrack_market_descr", "year", "channel_type", selcol)], 
				by = c("scantrack_market_descr", "year", "channel_type"))				
attr_mat	<- attr_mat[order(attr_mat$household_code, attr_mat$year, attr_mat$month),]	
cat("dim(mydata) =", dim(mydata), ".\n")
cat("dim(attr_mat) =", dim(attr_mat), ", dim(mydata)*R =",nrow(mydata)*R, ".\n")
selcol1		<- c("distance_haver", selcol)

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
	tmp1<- as.matrix(cbind(1, rec=rep(1, sum(sel))))	
	tmp2<- kronecker(tmp1, t(tmp0))
	colnames(tmp2)	<- paste(rep(c("Intercept", "Rec"), each = R-1), fmt_name[-beta0_base], sep=":")
	tmp3<- rep(0, R)
	tmp3[i]	<- 1
	tmp3<- kronecker(tmp1, t(tmp3))
	colnames(tmp3)	<- paste(rep(c("Intercept", "Rec"), each = R), fmt_name, sep=":")
	sel1	<- attr_mat$channel_type == fmt_name[i]
	tmp	<- as.matrix(attr_mat[sel1,selcol1][sel,])
	X1.08[[i]]	<- cbind(tmp2, tmp, tmp*1)
	X2.08[[i]]	<- cbind(tmp3, tmp, tmp*1)
	colnames(X1.08[[i]])	<- c(colnames(tmp2), colnames(tmp), paste(colnames(tmp), ":Rec", sep=""))
	colnames(X2.08[[i]])	<- c(colnames(tmp3), colnames(tmp), paste(colnames(tmp), ":Rec", sep=""))
	
	colnames(X1.07[[i]])	<- c(colnames(tmp2), colnames(tmp), paste(colnames(tmp), ":Rec", sep=""))
	colnames(X2.07[[i]])	<- c(colnames(tmp3), colnames(tmp), paste(colnames(tmp), ":Rec", sep=""))
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
								demo_mat = as.matrix(cbind(Intercept = 1, rec = 0, ln_inc = sim.unq$ln_inc)), 
								eps_draw = eps_draw, ret.idx = ret.idx, sim.y = NULL, ret.sim = TRUE, method = exp.method)
use.time		<- proc.time() - pct						
cat("2007 Baseline simulation finishes with", use.time[3]/60, "min.\n")

pct				<- proc.time()
sim.base08		<- SimWrapper_fn(omega_deriv, lambda, shr.par, base = beta0_base, 
								X1 = X1.08, X2=X2.08, price = sim.unq$price, 
								demo_mat = as.matrix(cbind(Intercept = 1, rec = 1, ln_inc = sim.unq$ln_inc)), 
								eps_draw = eps_draw, ret.idx = ret.idx, sim.y = NULL, ret.sim = TRUE, method = exp.method)
use.time		<- proc.time() - pct						
cat("2008 Baseline simulation finishes with", use.time[3]/60, "min.\n")

iv.change <- .1
(sel.retail	<- which(fmt_name == "Discount Store"))

#-------------------------------- #
# Simulation of fixed expenditure #
pct			<- proc.time()
sim.ls07	<- foreach(v = c(selcol1, "price"), .errorhandling = "remove", .packages=c("mgcv", "nloptr")) %dopar% {
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
							demo_mat = as.matrix(cbind(Intercept = 1, rec = 0, ln_inc = sim.unq$ln_inc)), 
							eps_draw = eps_draw, ret.idx = ret.idx, sim.y = sim.base07$y, ret.sim = TRUE)
	out1.avg	<- data.frame(out1$Average, retailer = fmt_name[sel.retail], Var = v, change = iv.change )
	out1$Average<- out1.avg
	return(out1)
}
use.time	<- proc.time() - pct
cat("Parallel simulation of retail price finishes with", use.time[3]/60, "min.\n")

pct			<- proc.time()
sim.ls08	<- foreach(v = c(selcol1, "price"), .errorhandling = "remove", .packages=c("mgcv", "nloptr")) %dopar% {
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
							demo_mat = as.matrix(cbind(Intercept = 1, rec = 1, ln_inc = sim.unq$ln_inc)), 
							eps_draw = eps_draw, ret.idx = ret.idx, sim.y = sim.base08$y, ret.sim = TRUE, method = exp.method)
	out1.avg	<- data.frame(out1$Average, retailer = fmt_name[sel.retail], Var = v, change = iv.change )
	out1$Average<- out1.avg
	return(out1)
}
use.time	<- proc.time() - pct
cat("Parallel simulation of retail price finishes with", use.time[3]/60, "min.\n")

# #--------------------------------------------------------- #
# # Simulation of effects on both expenditure and allocation #
# pct			<- proc.time()
# sim.2ls07	<- foreach(v = c(selcol1, "price"), .errorhandling = "remove", .packages=c("mgcv", "nloptr")) %dopar% {
# 	X1_new 		<- X1.07
# 	X2_new		<- X2.07
# 	price_new	<- sim.unq$price
# 	if(v == "price"){
# 		price_new[,sel.retail]		<- price_new[,sel.retail]*(1 + iv.change)
# 	}else if(v %in% c("ln_upc_per_mod", "ln_num_module")){
# 		X1_new[[sel.retail]][,v]	<- X1_new[[sel.retail]][,v] + log(1 + iv.change)
# 		X2_new[[sel.retail]][,v]	<- X2_new[[sel.retail]][,v] + log(1 + iv.change)
# 	}else{
# 		X1_new[[sel.retail]][,v]	<- X1_new[[sel.retail]][,v] * (1 + iv.change)
# 		X2_new[[sel.retail]][,v]	<- X2_new[[sel.retail]][,v] * (1 + iv.change)
# 	}
# 	
# 	out1	<- SimWrapper_fn(omega_deriv, lambda, shr.par, base = beta0_base, 
# 							X1 = X1_new, X2=X2_new, price = price_new, 
# 							demo_mat = as.matrix(cbind(Intercept = 1, rec = 0, ln_inc = sim.unq$ln_inc)), 
# 							eps_draw = eps_draw, ret.idx = ret.idx, sim.y = NULL, ret.sim = TRUE, method = exp.method)
# 	out1.avg	<- data.frame(out1$Average, retailer = fmt_name[sel.retail], Var = v, change = iv.change )
# 	out1$Average<- out1.avg
# 	return(out1)
# }
# use.time	<- proc.time() - pct
# cat("Parallel simulation of retail price finishes with", use.time[3]/60, "min.\n")
# 
# pct			<- proc.time()
# sim.2ls08	<- foreach(v = c(selcol1, "price"), .errorhandling = "remove", .packages=c("mgcv", "nloptr")) %dopar% {
# 	X1_new 		<- X1.08
# 	X2_new		<- X2.08
# 	price_new	<- sim.unq$price
# 	if(v == "price"){
# 		price_new[,sel.retail]	<- price_new[,sel.retail] * (1 + iv.change)
# 	}else if(v %in% c("ln_upc_per_mod", "ln_num_module")){
# 		X1_new[[sel.retail]][,v]	<- X1_new[[sel.retail]][,v] + log(1 + iv.change)
# 		X1_new[[sel.retail]][,paste(v,":Rec",sep="")]	<- X1_new[[sel.retail]][,paste(v,":Rec",sep="")] + log(1 + iv.change)
# 		X2_new[[sel.retail]][,v]	<- X2_new[[sel.retail]][,v] + log(1 + iv.change)
# 		X2_new[[sel.retail]][,paste(v,":Rec",sep="")]	<- X2_new[[sel.retail]][,paste(v,":Rec",sep="")] + log(1 + iv.change)
# 	}else{
# 		X1_new[[sel.retail]][,v]	<- X1_new[[sel.retail]][,v] * (1 + iv.change)
# 		X1_new[[sel.retail]][,paste(v,":Rec",sep="")]	<- X1_new[[sel.retail]][,paste(v,":Rec",sep="")] * (1 + iv.change)
# 		X2_new[[sel.retail]][,v]	<- X2_new[[sel.retail]][,v] * (1 + iv.change)
# 		X2_new[[sel.retail]][,paste(v,":Rec",sep="")]	<- X2_new[[sel.retail]][,paste(v,":Rec",sep="")] * (1 + iv.change)
# 	}
# 	
# 	out1	<- SimWrapper_fn(omega_deriv, lambda, shr.par, base = beta0_base, 
# 							X1 = X1_new, X2=X2_new, price = price_new, 
# 							demo_mat = as.matrix(cbind(Intercept = 1, rec = 1, ln_inc = sim.unq$ln_inc)), 
# 							eps_draw = eps_draw, ret.idx = ret.idx, sim.y = NULL, ret.sim = TRUE, method = exp.method)
# 	out1.avg	<- data.frame(out1$Average, retailer = fmt_name[sel.retail], Var = v, change = iv.change )
# 	out1$Average<- out1.avg
# 	return(out1)
# }
# use.time	<- proc.time() - pct
# cat("Parallel simulation of retail price finishes with", use.time[3]/60, "min.\n")

stopCluster(cl)
cat("Stop clustering. \n")

# Check the marginal effects. 
selcol2	<- c(selcol1, "price")
tmp0	<- colSums(sim.base07$Average[, fmt_name], na.rm=T)
tmp1	<- sapply(sim.ls07, function(x) colSums(x$Average[,gsub(" ", ".", fmt_name)]))
tmp		<- setNames(tmp1[sel.retail,] - tmp0[sel.retail], selcol2)
cat("Impact on exp. share in 2007:\n"); print(tmp); cat("\n")
cat("Elasticity in 2007 when expenditure is fixed:\n"); print(tmp/tmp0[sel.retail]/iv.change); cat("\n")
# tmp1	<- sapply(sim.2ls07, function(x) colSums(x$Average[,gsub(" ", ".", fmt_name)]))
# tmp		<- setNames(tmp1[sel.retail,] - tmp0[sel.retail], selcol2)
# cat("Impact on exp. share in 2007:\n"); print(tmp); cat("\n")
# cat("Elasticity in 2007 when expenditure is free:\n"); print(tmp/tmp0[sel.retail]/iv.change); cat("\n")

tmp0	<- colSums(sim.base08$Average[, fmt_name], na.rm=T)
tmp1	<- sapply(sim.ls08, function(x) colSums(x$Average[,gsub(" ", ".", fmt_name)]))
tmp		<- setNames(tmp1[sel.retail,] - tmp0[sel.retail], selcol2)
cat("Impact on exp. share in 2008:\n"); print(tmp); cat("\n")
cat("Elasticity in 2008 when expenditure is fixed:\n"); print(tmp/tmp0[sel.retail]/iv.change); cat("\n")
# tmp1	<- sapply(sim.2ls08, function(x) colSums(x$Average[,gsub(" ", ".", fmt_name)]))
# tmp		<- setNames(tmp1[sel.retail,] - tmp0[sel.retail], selcol2)
# cat("Impact on exp. share in 2008:\n"); print(tmp); cat("\n")
# cat("Elasticity in 2008 when expenditure is free:\n"); print(tmp/tmp0[sel.retail]/iv.change); cat("\n")

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

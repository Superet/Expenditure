library(ggplot2)
library(reshape2)
library(maxLik)
library(evd)
library(data.table)
library(doParallel)
library(foreach)
library(nloptr)
library(mgcv)
options(error = quote({dump.frames(to.file = TRUE)}))

# seg_id	<- as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))

# args <- commandArgs(trailingOnly = TRUE)
# print(args)
# if(length(args)>0){
#     for(i in 1:length(args)){
#       eval(parse(text=args[[i]]))
#     }
# }
# cat("seg_id =", seg_id, ".\n")

# setwd("~/Documents/Research/Store switching/Processed_data/processed data_20160809")
# plot.wd	<- '~/Desktop'
# sourceCpp(paste("../../Exercise/Multiple_discrete_continuous_model/", model_name, ".cpp", sep=""))
# source("../../Exercise/Multiple_discrete_continuous_model/0_Efficient_Allocation_function.R")
# source("../../Exercise/main/ctrfact_sim_functions.r")

# setwd("/home/brgordon/ccv103/Exercise/run")
# setwd("/kellogg/users/marketing/2661703/Exercise/run")
setwd("/sscc/home/c/ccv103/Exercise/run")
# setwd("U:/Users/ccv103/Documents/Research/Store switching/run")

run_id		<- 1
plot.wd		<- getwd()
make_plot	<- TRUE
ww			<- 10
ar			<- .6

source("0_Efficient_Allocation_function.R")
source("ctrfact_sim_functions.r")

# Load estimation data 
ver.date	<- "2016-09-11"
loadf	<- paste("estrun_",run_id,"/MDCEV_est_", ver.date,".rdata",sep="")
loadf
load(loadf)

# Set simulation parameters
week.price		<- FALSE
interp.method	<- "spline"					# Spline interpolation or Chebyshev interpolation
exp.method		<- "Utility"				# Utility maximization or finding roots for first order condition
trim.alpha		<- 0
numsim			<- 500
draw.par		<- FALSE
sim.omega		<- FALSE
fname			<- paste("price_elst",sep="")
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
lambda 	<- sapply(exp.est.ls, coef)
shr.par	<- sapply(shr.est.ls, coef)
ret.idx	<- 11:15		# Column index flagging the retail attribugtes in X1

#-----------------------#
# Construct income data # 
selyr	<- 2007
sel		<- mydata.full$year == selyr & !duplicated(mydata.full[,c("household_code", "year")]) # & mydata.full$household_code < 2250000
mydata.full$ln_inc	<- log(mydata.full$income_mid)
sim.unq	<- mydata.full[sel,c("household_code", "ln_inc","rec", "segment")]
sim.unq$price	<- as.matrix(mydata.full[sel,paste("PRC_", gsub("\\s", "_", fmt_name), sep="")])
# sim.unq	<- sim.unq[sample(1:nrow(sim.unq), 100),]

# Normalize distance by zipcode
tmp		<- data.table(hh_dist[,c("household_code", "panelist_zip_code","year", "channel_type", selcol)])
tmp		<- tmp[, maxd:=max(distance_haver, na.rm=T), by = list(panelist_zip_code)]
tmp		<- tmp[, distance_haver := distance_haver/maxd]
tmp		<- data.frame(subset(tmp, year == selyr))
attr_mat	<- merge(sim.unq[,c("household_code", "ln_inc")], tmp[,c("household_code", "channel_type", selcol)], 
				by = c("household_code"), all.x=T)
attr_mat	<- attr_mat[order(attr_mat$household_code),]	
cat("dim(sim.unq) =", dim(sim.unq), ".\n")
cat("dim(attr_mat) =", dim(attr_mat), ", dim(sim.unq)*R =",nrow(sim.unq)*R, ".\n")
beta0_base 	<- which(fmt_name == "Grocery")

X2	<- X1 	<- vector("list", length=length(fmt_name))
for(i in 1:length(fmt_name)){
	tmp0<- rep(0, R-1)
	if(i<beta0_base){
		tmp0[i]	<- 1
	}else if(i > beta0_base){
		tmp0[(i-1)]	<- 1
	}
	tmp1<- as.matrix(cbind(1, rec=sim.unq$rec))	
	tmp2<- kronecker(tmp1, t(tmp0))
	colnames(tmp2)	<- paste(rep(c("Intercept", "Rec"), each = R-1), fmt_name[-beta0_base], sep=":")
	tmp3<- rep(0, R)
	tmp3[i]	<- 1
	tmp3<- kronecker(tmp1, t(tmp3))
	colnames(tmp3)	<- paste(rep(c("Intercept", "Rec"), each = R), fmt_name, sep=":")
	sel	<- attr_mat$channel_type == fmt_name[i]	
	tmp	<- as.matrix(attr_mat[sel,selcol])
	X1[[i]]	<- cbind(tmp2, tmp, tmp*sim.unq$rec)
	X2[[i]]	<- cbind(tmp3, tmp, tmp*sim.unq$rec)
	colnames(X1[[i]])	<- c(colnames(tmp2), colnames(tmp), paste(colnames(tmp), ":Rec", sep=""))
	colnames(X2[[i]])	<- c(colnames(tmp3), colnames(tmp), paste(colnames(tmp), ":Rec", sep=""))
}
nx 		<- ncol(X1[[1]])
cat("nx =", nx,"\n")

sim.unq$X<- as.matrix( do.call(cbind, lapply(X1, function(x) x[,ret.idx])))
cat("dim(sim.unq) =", dim(sim.unq), "\n")

# Retail attributes for simulation
X1.07	<- X1
X2.07	<- X2
X2.08	<- X1.08 	<- vector("list", length=length(fmt_name))
for(i in 1:length(fmt_name)){
	tmp0<- rep(0, R-1)
	if(i<beta0_base){
		tmp0[i]	<- 1
	}else if(i > beta0_base){
		tmp0[(i-1)]	<- 1
	}
	tmp1<- as.matrix(cbind(1, rec=rep(1, nrow(sim.unq))))	
	tmp2<- kronecker(tmp1, t(tmp0))
	colnames(tmp2)	<- paste(rep(c("Intercept", "Rec"), each = R-1), fmt_name[-beta0_base], sep=":")
	tmp3<- rep(0, R)
	tmp3[i]	<- 1
	tmp3<- kronecker(tmp1, t(tmp3))
	colnames(tmp3)	<- paste(rep(c("Intercept", "Rec"), each = R), fmt_name, sep=":")
	sel1	<- attr_mat$channel_type == fmt_name[i]
	tmp	<- as.matrix(attr_mat[sel1,selcol])
	X1.08[[i]]	<- cbind(tmp2, tmp, tmp*1)
	X2.08[[i]]	<- cbind(tmp3, tmp, tmp*1)
	colnames(X1.08[[i]])	<- c(colnames(tmp2), colnames(tmp), paste(colnames(tmp), ":Rec", sep=""))
	colnames(X2.08[[i]])	<- c(colnames(tmp3), colnames(tmp), paste(colnames(tmp), ":Rec", sep=""))
	
	# colnames(X1.07[[i]])	<- c(colnames(tmp2), colnames(tmp), paste(colnames(tmp), ":Rec", sep=""))
	# colnames(X2.07[[i]])	<- c(colnames(tmp3), colnames(tmp), paste(colnames(tmp), ":Rec", sep=""))
}

#-------------------#
# Take random draws #
set.seed(666)
eps_draw	<- lapply(1:nseg, function(i) matrix(rgev(numsim*R, scale = exp(shr.par["ln_sigma",i])), numsim, R))
par_draw <- NULL

rm(list = intersect(ls(), c("hh_dist", "hh_dpt","interp.method", "loadf", "X1", "X2", "shr")))

##############
# Simulation #
##############
# Register parallel computing
mycore 	<- 2
cl		<- makeCluster(mycore, type = "FORK")
# cl		<- makeCluster(mycore, type = "PSOCK")			# Register cores in Windows
registerDoParallel(cl)
cat("Register", mycore, "core parallel computing. \n")

# ------------------ #
# Basline simulation #
# Simulate expenditure and expenditure share. 
pct				<- proc.time()
sim.base07		<- foreach(i = 1:nseg, .errorhandling = "remove", .packages=c("mgcv", "nloptr")) %dopar% {
	sel		<- sim.unq$segment == i
	out		<- SimWrapper_fn(omega.ls[[i]], lambda[,i], shr.par[,i], base = beta0_base, 
				X1 = lapply(X1.07, function(x) x[sel,]), X2=lapply(X2.07, function(x) x[sel,]), price = sim.unq$price[sel,], 
				demo_mat = as.matrix(cbind(Intercept = 1, rec = 0, ln_inc = sim.unq[sel,"ln_inc"])), 
				eps_draw = eps_draw[[i]], ret.idx = ret.idx, sim.y = NULL, ret.sim = TRUE, method = exp.method)
	out$Average	<- cbind(segment = i, household_code = sim.unq[sel,"household_code"], out$Average)
	return(out)
}
use.time		<- proc.time() - pct						
cat("2007 Baseline simulation finishes with", use.time[3]/60, "min.\n")

iv.change <- .1

#-------------------------------- #
# Simulation of fixed expenditure #
pct			<- proc.time()
sim.ls07	<- foreach(sel.retail = 1:R, .errorhandling = "remove", .packages=c("mgcv", "nloptr")) %dopar% {
	v 		<- "price"
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
	
	out1	<- data.frame()
	for(i in 1:nseg){
		sel	<- sim.unq$segment == i
		tmp	<- SimWrapper_fn(omega.ls[[i]], lambda[,i], shr.par[,i], base = beta0_base, 
					X1 = lapply(X1_new, function(x) x[sel,]), X2=lapply(X2_new, function(x) x[sel,]), price = price_new[sel,], 
					demo_mat = as.matrix(cbind(Intercept = 1, rec = 0, ln_inc = sim.unq[sel,"ln_inc"])), 
					eps_draw = eps_draw[[i]], ret.idx = ret.idx, sim.y = sim.base07[[i]]$y, ret.sim = FALSE, method = exp.method)
		tmp	<- data.frame(tmp, retailer = fmt_name[sel.retail], Var = v, change = iv.change, segment = i, household_code = sim.unq[sel,"household_code"])	
		out1<- rbind(out1, tmp)
	}
	return(out1)
}
use.time	<- proc.time() - pct
cat("Parallel simulation of retail price finishes with", use.time[3]/60, "min.\n")
names(sim.ls07)	<- sapply(sim.ls07, function(x) as.character(unique(x$retailer)))

stopCluster(cl)
cat("Stop clustering. \n")

##############################
# Summarize marginal effects # 
##############################
fmt_name1	<- gsub("\\s", ".", fmt_name)
base.07	<- data.frame(do.call(rbind, lapply(sim.base07, function(x) x$Average)))
max(abs(base.07[,"household_code"] - as.numeric(sim.ls07[[1]]$household_code)))

# Mean own-elasticity across population 
tmp1	<- sapply(1:R, function(i) (sim.ls07[[i]][,fmt_name1[i]] - base.07[,fmt_name1[i]])/base.07[,fmt_name1[i]]/iv.change )
cat("Summary of individual elasticity:\n");print(summary(tmp1)); cat("\n")

# Overall Revenue
tmp0	<- colSums(base.07[,gsub("\\s", ".",fmt_name)])
tmp1	<- lapply(sim.ls07, function(x) colSums(x[,gsub("\\s", ".",fmt_name)]))
tmp1	<- sapply(tmp1, function(x) (x- tmp0)/tmp0/iv.change)
cat("Cross Elasticity on overall revenue:\n"); print(tmp1); cat("\n")

# Market share
tmp0	<- colSums(base.07[,gsub("\\s", ".",fmt_name)])
tmp0	<- tmp0/sum(tmp0)
tmp1	<- lapply(sim.ls07, function(x) colSums(x[,gsub("\\s", ".",fmt_name)]))
tmp1	<- lapply(tmp1, function(x) x/sum(x))
tmp.tab1	<- sapply(tmp1, function(x) x - tmp0)
tmp.tab2	<- sapply(tmp1, function(x) (x - tmp0)/tmp0/iv.change)
cat("Change in market share:\n"); print(tmp.tab1); cat("\n")
cat("Elasticity on market share:\n"); print(tmp.tab2); cat("\n")

################
# Save results #
################
rm(list  = intersect(ls(), c("ar", "args", "cl", "ggtmp", "ii", "i", "make_plot", "mycore", 
			"myfix", "plot.wd", "sel", "selyr", "price", "sol", "sol.top", "sol.top2", "tmp", "tmp_coef", 
			"tmp_price", "tmp1", "tmp2", "use.time", "ver.date", "var1", "ww", "f", 
			"pct", "tmpd1", "tmpd2", "tmpdat", "out", "out1", "tmp1", "tmp2", "tmp3","tmp.ls",
			"Allocation_constr_fn", "Allocation_fn", "Allocation_nlop_fn", "cheb.1d.basis", "cheb.basis", "chebfun", 
			"exp_fn", "expFOC_fn", "incl_value_fn", "mysplfun", "mytrimfun", "param_assignR", "simExp_fn", "SimOmega_fn", 
			"SimWrapper_fn", "solveExp_fn", "spl2dfun", "uP_fn", "uGrad_fn")))

save.image(paste("estrun_",run_id, "/", fname, ".rdata",sep=""))

cat("This program is done.\n")

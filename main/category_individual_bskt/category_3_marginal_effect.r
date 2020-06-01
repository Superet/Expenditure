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

run_id		<- 3
plot.wd		<- getwd()
make_plot	<- TRUE
ww			<- 10
ar			<- .6

source("0_Efficient_Allocation_function.R")
source("ctrfact_sim_functions.r")

# Load estimation data 
seldpt		<- "DG"
ver.date	<- "2016-09-11"
loadf	<- paste("estrun_",run_id,"/est_",seldpt,"_", ver.date,".rdata",sep="")
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
fname			<- paste(seldpt,"_marginal_effect",sep="")
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
sel		<- sapply(exp.est.ls, function(x) !is.null(x))
lambda 	<- sapply(exp.est.ls[sel], coef)
shr.par	<- sapply(shr.est.ls[sel], coef)
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
sel		<- hh_dpt$department == dpt_name[seldpt]
tmp		<- merge(hh_dpt[sel,], hh_dist[,c("household_code","year", "channel_type","panelist_zip_code","distance_haver")], 
				 by = c("household_code","year", "channel_type"))
tmp		<- data.table(tmp[,c("household_code", "panelist_zip_code","year", "channel_type", selcol)])
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
(exp.bound	<- quantile(mydata.full$dol, c(.01, .99)) )

rm(list = intersect(ls(), c("hh_dpt", "hh_month", "X1", "X2", "attr_mat", "pan")))

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
				eps_draw = eps_draw[[i]], ret.idx = ret.idx, sim.y = NULL, ret.sim = TRUE, method = exp.method, 
				exp.lower = exp.bound[1], exp.upper = exp.bound[2])
	out$Average	<- cbind(segment = i, household_code = sim.unq[sel,"household_code"], out$Average)
	return(out)
}
use.time		<- proc.time() - pct						
cat("2007 Baseline simulation finishes with", use.time[3]/60, "min.\n")

pct				<- proc.time()
sim.base08		<- foreach(i = 1:nseg, .errorhandling = "remove", .packages=c("mgcv", "nloptr")) %dopar% {
	sel		<- sim.unq$segment == i
	out		<- SimWrapper_fn(omega.ls[[i]], lambda[,i], shr.par[,i], base = beta0_base, 
				X1 = lapply(X1.08, function(x) x[sel,]), X2=lapply(X2.08, function(x) x[sel,]), price = sim.unq$price[sel,], 
				demo_mat = as.matrix(cbind(Intercept = 1, rec = 1, ln_inc = sim.unq[sel,"ln_inc"])), 
				eps_draw = eps_draw[[i]], ret.idx = ret.idx, sim.y = NULL, ret.sim = TRUE, method = exp.method, 
				exp.lower = exp.bound[1], exp.upper = exp.bound[2])
	out$Average	<- cbind(segment = i, household_code = sim.unq[sel,"household_code"], out$Average)
	return(out)				
}
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
names(sim.ls07)	<- sapply(sim.ls07, function(x) as.character(unique(x$Var)))

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
	
	out1	<- data.frame()
	for(i in 1:nseg){
		sel	<- sim.unq$segment == i
		tmp	<- SimWrapper_fn(omega.ls[[i]], lambda[,i], shr.par[,i], base = beta0_base, 
					X1 = lapply(X1_new, function(x) x[sel,]), X2=lapply(X2_new, function(x) x[sel,]), price = price_new[sel,], 
					demo_mat = as.matrix(cbind(Intercept = 1, rec = 1, ln_inc = sim.unq[sel,"ln_inc"])), 
					eps_draw = eps_draw[[i]], ret.idx = ret.idx, sim.y = sim.base08[[i]]$y, ret.sim = FALSE, method = exp.method)
		tmp	<- data.frame(tmp, retailer = fmt_name[sel.retail], Var = v, change = iv.change, segment = i, household_code = sim.unq[sel,"household_code"])	
		out1<- rbind(out1, tmp)
	}
		
	# out1	<- SimWrapper_fn(omega_deriv, lambda, shr.par, base = beta0_base, 
	# 						X1 = X1_new, X2=X2_new, price = price_new, 
	# 						demo_mat = as.matrix(cbind(Intercept = 1, rec = 1, ln_inc = sim.unq$ln_inc)), 
	# 						eps_draw = eps_draw, ret.idx = ret.idx, sim.y = sim.base08$y, ret.sim = TRUE, method = exp.method)
	# out1.avg	<- data.frame(out1$Average, retailer = fmt_name[sel.retail], Var = v, change = iv.change )
	# out1$Average<- out1.avg
	return(out1)
}
use.time	<- proc.time() - pct
cat("Parallel simulation of retail price finishes with", use.time[3]/60, "min.\n")
names(sim.ls08)	<- sapply(sim.ls08, function(x) as.character(unique(x$Var)))

stopCluster(cl)
cat("Stop clustering. \n")

##############################
# Summarize marginal effects # 
##############################
sel.retail	<- "Discount.Store" 
selcol1	<- c(selcol, "price")
base.07	<- data.frame(do.call(rbind, lapply(sim.base07, function(x) x$Average)))
base.08	<- data.frame(do.call(rbind, lapply(sim.base08, function(x) x$Average)))
max(abs(base.07[,"household_code"] - as.numeric(sim.ls07[[1]]$household_code)))

# Mean elasticity across population 
tmp1	<- lapply(sim.ls07, function(x) (x[,sel.retail] - base.07[,sel.retail])/base.07[,sel.retail]/x[,"change"])
tmp2	<- lapply(sim.ls08, function(x) (x[,sel.retail] - base.08[,sel.retail])/base.08[,sel.retail]/x[,"change"])
tmp		<- vector("list", length = length(tmp1))
for(i in 1:length(tmp1)){
	sel	<- !(is.na(tmp1[[i]]) | is.na(tmp2[[i]]) | abs(tmp1[[i]]) == Inf | abs(tmp2[[i]]) == Inf)
	tmp[[i]]	<- t.test(tmp2[[i]][sel], tmp1[[i]][sel])
}
tmp		<- do.call(rbind, lapply(tmp, function(x) c(x$estimate, x$p.value)))
colnames(tmp)	<- c("2008", "2007", "pvalue")
rownames(tmp)	<- names(tmp1)
cat("T test of mean elasticity:\n"); print(tmp); cat("\n")

# Overall Revenue
tmp0	<- colSums(base.07[,gsub("\\s", ".",fmt_name)])
tmp1	<- lapply(sim.ls07, function(x) colSums(x[,gsub("\\s", ".",fmt_name)]))
tmp1	<- sapply(tmp1, function(x) (x[sel.retail] - tmp0[sel.retail])/tmp0[sel.retail]/iv.change)
tmp0	<- colSums(base.08[,gsub("\\s", ".",fmt_name)])
tmp2	<- lapply(sim.ls08, function(x) colSums(x[,gsub("\\s", ".",fmt_name)]))
tmp2	<- sapply(tmp2, function(x) (x[sel.retail] - tmp0[sel.retail])/tmp0[sel.retail]/iv.change)
tmp.tab	<- cbind(tmp1, tmp2)
colnames(tmp.tab)	<- c("2007", "2008")
cat("Elasticity on overall revenue:\n"); print(tmp.tab); cat("\n")

# Market share
tmp0	<- colSums(base.07[,gsub("\\s", ".",fmt_name)])
tmp0	<- tmp0/sum(tmp0)
tmp1	<- lapply(sim.ls07, function(x) colSums(x[,gsub("\\s", ".",fmt_name)]))
tmp1	<- lapply(tmp1, function(x) x/sum(x))
tmp02	<- colSums(base.08[,gsub("\\s", ".",fmt_name)])
tmp02	<- tmp02/sum(tmp02)
tmp2	<- lapply(sim.ls08, function(x) colSums(x[,gsub("\\s", ".",fmt_name)]))
tmp2	<- lapply(tmp2, function(x) x/sum(x))
tmp.tab1	<- cbind(sapply(tmp1, function(x) (x[sel.retail] - tmp0[sel.retail])), 
					 sapply(tmp2, function(x) (x[sel.retail] - tmp02[sel.retail])) )
tmp.tab2	<- cbind(sapply(tmp1, function(x) (x[sel.retail] - tmp0[sel.retail])/tmp0[sel.retail]/iv.change), 
					 sapply(tmp2, function(x) (x[sel.retail] - tmp02[sel.retail])/tmp02[sel.retail]/iv.change))
colnames(tmp.tab1)	<- colnames(tmp.tab2) <- c("2007", "2008")
cat("Change in market share:\n"); print(tmp.tab1); cat("\n")
cat("Elasticity on market share:\n"); print(tmp.tab2); cat("\n")

# Distance marginal effect
# Normalize distance by zipcode
tmp		<- data.table(hh_dist[,c("household_code", "panelist_zip_code","year", "channel_type", "distance_haver")])
tmp		<- tmp[, maxd:=max(distance_haver, na.rm=T), by = list(panelist_zip_code)]
tmp		<- tmp[, distance_haver := distance_haver/maxd]
tmp		<- data.frame(tmp)

ggtmp	<- merge(sim.unq[,c("household_code", "segment")], subset(tmp, year == 2007 & channel_type == "Discount Store"), 
				by = "household_code", all.x=T)
tmp		<- setNames(base.07[,sel.retail], base.07[,"household_code"])
ggtmp$old<- tmp[as.character(ggtmp$household_code)]
tmp		<- setNames(sim.ls07[["distance_haver"]][,sel.retail], sim.ls07[["distance_haver"]][,"household_code"])
ggtmp$new<- tmp[as.character(ggtmp$household_code)]
ggtmp$d_d	<- with(ggtmp, .1*distance_haver*maxd)
ggtmp$d_exp	<- with(ggtmp, new - old)
summary(ggtmp$d_d)
summary(ggtmp$d_exp/ggtmp$d_d)
summary(with(ggtmp, d_exp/old/(d_d/(maxd*distance_haver)) ))

################
# Save results #
################
rm(list  = intersect(ls(), c("ar", "args", "cl", "ggtmp", "ii", "i", "make_plot", "mycore", "cpi.adj","sel1", 
			"myfix", "plot.wd", "sel", "selyr", "price", "sol", "sol.top", "sol.top2", "tmp", "tmp_coef", 
			"tmp_price", "tmp1", "tmp2", "use.time", "ver.date", "var1", "ww", "f", 
			"pct", "tmpd1", "tmpd2", "tmpdat", "out", "out1", "tmp1", "tmp2", "tmp3","tmp.ls","tmp0",
			"Allocation_constr_fn", "Allocation_fn", "Allocation_nlop_fn", "cheb.1d.basis", "cheb.basis", "chebfun", 
			"exp_fn", "expFOC_fn", "incl_value_fn", "mysplfun", "mytrimfun", "param_assignR", "simExp_fn", "SimOmega_fn", 
			"SimWrapper_fn", "solveExp_fn", "spl2dfun", "uP_fn", "uGrad_fn")))

save.image(paste("estrun_",run_id, "/", fname, ".rdata",sep=""))

cat("This program is done.\n")

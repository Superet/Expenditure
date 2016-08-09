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
# setwd("~/Documents/Research/Store switching/processed data/Estimation")
# plot.wd	<- '~/Desktop'
# source("../../Exercise/Multiple_discrete_continuous_model/0_Allocation_function.R")
# source("../../Exercise/main/share_allocation_month/ctrfact_sim_functions.r")

# setwd("/home/brgordon/ccv103/Exercise/run")
# setwd("/kellogg/users/marketing/2661703/Exercise/run")
# setwd("/sscc/home/c/ccv103/Exercise/run")
setwd("U:/Users/ccv103/Documents/Research/Store switching/run")

run_id		<- 9
plot.wd		<- getwd()
make_plot	<- TRUE
ww			<- 10
ar			<- .6

source("0_Allocation_function.R")
source("ctrfact_sim_functions.r")

# Load estimation data 
ver.date	<- "2016-07-13"
cpi.adj		<- TRUE

if(cpi.adj){
	loadf	<- paste("estrun_",run_id,"/MDCEV_cpi_est_seg",seg_id,"_", ver.date,".rdata",sep="")
}else{
	loadf	<- paste("estrun_",run_id,"/MDCEV_est_seg",seg_id,"_", ver.date,".rdata",sep="")
}
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
fname			<- paste("nostate_marginal_seg",seg_id,sep="")
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

#-----------------------#
# Construct income data # 
# Take the households' income in 2007 as basis 
selyr	<- 2007
tmp		<- data.table(subset(mydata, year %in% selyr))
tmp		<- tmp[,list(income = unique(income_midvalue)), by = list(household_code, year)]
sim.data<- data.frame(tmp)[,c("household_code","income")]
sim.unq		<- data.frame(income2007 = unique(sim.data[,-1]))

# Counterfactual scenario: income is lower by 10%. 
my.change		<- 0	#.1
lnInc_08		<- lnInc + log(1 - my.change)
sim.unq$income2008	<- (1 - my.change) * sim.unq$income2007		
sim.unq$Inc07	<- log(sim.unq[,"income2007"])
sim.unq$Inc08	<- log(sim.unq[,"income2008"])
sim.unq			<- sim.unq[order(sim.unq$income2007),]
cat("dim(sim.unq) =", dim(sim.unq), "\n")
sim.unq			<- sim.unq[c(1,nrow(sim.unq)),]

#----------------------------#
# Average price in year 2007
if(week.price){
	tmp 		<- dcast(subset(price_dat, year %in% selyr),	 
						scantrack_market_descr + year + biweek ~ channel_type, value.var = "bsk_price_paid_2004") 
	price.07	<- setNames( colMeans(as.matrix(tmp[,4:(3+R)]), na.rm=T), fmt_name)
}else{
	tmp 		<- dcast(subset(price_dat, year %in% selyr),	 
						scantrack_market_descr + year ~ channel_type, value.var = "bsk_price_paid_2004") 
	price.07	<- setNames( colMeans(as.matrix(tmp[,3:(2+R)]), na.rm=T), fmt_name)
}
cat("The average price level in 2007:\n"); print(price.07);cat("\n")

# Average retail attributes in 2007
selcol	<- c("size_index", "ln_upc_per_mod", "ln_num_module","overall_prvt")
X_list0<- setNames( lapply(fmt_name, function(x) colMeans(as.matrix(subset(fmt_attr, channel_type == x & year%in% selyr)[,selcol]))), 
					 fmt_name)
cat("The average retail attributes in 2007:\n"); print(do.call(rbind, X_list0)); cat("\n")

# Expand X_list and price to match the nobs of income 
price.07	<- rep(1, nrow(sim.unq)) %*% matrix(price.07, nrow = 1)
colnames(price.07)	<- fmt_name
X_list08	<- X_list07	<- setNames(vector("list", R), fmt_name)
for(i in 1:R){
	tmp2	<- tmp1	<- matrix(0, nrow(sim.unq), R-1)
	if(i<beta0_base){
		tmp1[,i]	<- 0
		tmp2[,i]	<- 1
	}else if(i>beta0_base){
		tmp1[,(i-1)]	<- 0
		tmp2[,(i-1)]	<- 1
	}
	X_list07[[i]]	<- cbind(tmp1, rep(1, nrow(sim.unq)) %*% matrix(X_list0[[i]], nrow = 1), matrix(0, nrow(sim.unq), length(selcol)))
	X_list08[[i]]	<- cbind(tmp2, rep(1, nrow(sim.unq)) %*% matrix(X_list0[[i]], nrow = 1), rep(1, nrow(sim.unq)) %*% matrix(X_list0[[i]], nrow = 1))
	colnames(X_list07[[i]])	<- colnames(X_list08[[i]])	<- c(fmt_name[-beta0_base], selcol, paste("rec*", selcol,sep=""))
}

#-------------------#
# Take random draws #
set.seed(666)
eps_draw	<- matrix(rgev(numsim*R, scale = exp(shr.par["ln_sigma"])), numsim, R)
if(draw.par){
	par_se	<- c(sqrt(diag(vcov(sol.top2))), sqrt(diag(vcov(sol))) )
	par_se[is.na(par_se)]	<- 0
	par_draw <- sapply(par_se, function(x) rnorm(numsim, mean = 0, sd = x))
}else{
	par_draw <- NULL
}

# Compare the deterministic vs. random components
tmp		<- param_assignR(shr.par, nx, R, beta0_base)
tmp1	<- sapply(1:R, function(i) X_list[[i]] %*% tmp$beta + tmp$beta0[i])
tmp2	<- sapply(1:R, function(i) X_list07[[i]] %*% tmp$beta + tmp$beta0[i])
tmp3	<- sapply(1:R, function(i) X_list08[[i]] %*% tmp$beta + tmp$beta0[i])
cat("The expected xbeta from the actual data:\n"); print((colMeans(tmp1))); cat("\n")
cat("The expected xbeta of the simulation X in year 2007:\n"); print(tmp2); cat("\n")
cat("The expected xbeta of the simulation X in year 2008:\n"); print(tmp3); cat("\n")
cat("Summary stats of random draws:\n"); print(summary(c(eps_draw))); cat("\n")

##############
# Simulation #
##############
if(interp.method == "spline"){
	y.nodes		<- quantile(mydata$dol, c(0:50)/50)
	y.nodes		<- sort(unique(c(y.nodes , seq(600, 1000, 100)) ))
}else{
	# Set up Chebyshev interpolation
	GH_num_nodes<- 100
	y.interval	<- c(.1, 1000)
	y.nodes		<- chebknots(GH_num_nodes, interval = y.interval)[[1]]
}
numnodes<- length(y.nodes)

# Register parallel computing
mycore 	<- 3
# cl		<- makeCluster(mycore, type = "FORK")
cl		<- makeCluster(mycore, type = "PSOCK")			# Register cores in Windows
registerDoParallel(cl)
cat("Register", mycore, "core parallel computing. \n")

# Simulate expenditure and expenditure share. 
pct				<- proc.time()
sim.base07		<- SimWrapper_fn(omega_deriv, ln_inc = sim.unq$Inc07, rec = rep(0, nrow(sim.unq)), lambda = lambda, param_est = shr.par, base = beta0_base, 
					X_list = X_list07, price = price.07, eps_draw = eps_draw, method = exp.method, ret.sim = TRUE, par.draw = par_draw,
					share.state = NULL)
use.time		<- proc.time() - pct						
cat("2007 Baseline simulation finishes with", use.time[3]/60, "min.\n")

pct				<- proc.time()
sim.base08		<- SimWrapper_fn(omega_deriv1, ln_inc = sim.unq$Inc08, rec = rep(1, nrow(sim.unq)), lambda = lambda, param_est = shr.par, 
					base = beta0_base, X_list = X_list08, price = price.07, eps_draw = eps_draw, method = exp.method, ret.sim = TRUE, 
					par.draw = par_draw, share.state = NULL)
use.time		<- proc.time() - pct						
cat("2008 Baseline simulation finishes with", use.time[3]/60, "min.\n")

# Check the baseline simulation 
tmp1	<- sim.base07$Average[,fmt_name]/sim.base07$Average[,"Exp"]
tmp2	<- sim.base08$Average[,fmt_name]/sim.base08$Average[,"Exp"]
tmp		<- tmp2-tmp1
cat("Summary of baseline income elasticity:\n"); print(tmp); cat("\n")

# Plot omega over y
tmp		<- seq(50, 250, 10)
tmpd1	<- data.frame(lnInc = rep(sim.unq$Inc07, each = length(tmp)), y = rep(tmp, length(sim.unq$Inc07)))
tmp1	<- omega_deriv(tmpd1)
tmpd2	<- data.frame(lnInc = rep(sim.unq$Inc08, each = length(tmp)), y = rep(tmp, length(sim.unq$Inc08)))
tmp2	<- omega_deriv(tmpd2)
ggtmp	<- rbind(cbind(tmpd1, omega = tmp1, year = 2007), cbind(tmpd2, omega = tmp2, year = 2008))
ggtmp$Income	<- exp(ggtmp$lnInc)

# Check if utility function at the top level is concave
tmp		<- seq(300, 400, 5)
u		<- matrix(NA, length(tmp), length(lnInc), dimnames = list(y = tmp, lnInc = lnInc))
for(i in 1:length(lnInc)){
	f 	<- function(y, ...){
		newx	<- data.frame(lnInc = lnInc[i], y = y)
		omega_deriv(newx, ...)
	}
	u[,i]<- -exp_fn(tmp, lambda, f, ln_inc = lnInc[i], rec = 0)
}
ggtmp2	<- melt(u, value.name = "utility")

if(make_plot){
	pdf(paste(plot.wd, "/estrun_",run_id,"/graph_checkbase_", fname, ".pdf", sep=""), width = ww/ar, height = ww)
	print(ggplot(ggtmp, aes(y, omega, col = factor(Income))) + geom_line() + 
			facet_wrap(~year) + 
			guides(col = FALSE) + 
			labs(title = "Inclusive value omega(Inc, y)")
		)
	print(ggplot(ggtmp2, aes(y, utility)) + geom_point() + geom_line() + facet_wrap(~lnInc, scales = "free"))
	dev.off()
}

#------------------------#
# Simulation of strategy #
# Assume only grocery stores change retail attributes. 
iv.change <- .1
(sel.retail	<- which(fmt_name == "Discount Store"))

#-------------------------------#
# Simulation of price reduction # 
pct			<- proc.time()
sim.ls07	<- foreach(v = selcol, .errorhandling = "remove", .packages=c("mgcv", "nloptr"), .combine = 'rbind') %dopar% {
	out1.avg	<- data.frame()
	out1.alc	<- array(NA, c(numsim, nrow(sim.unq), R), 
						dimnames = list(iter = 1:numsim, lnInc = sim.unq[,"Inc07"], retail = fmt_name))
	X_new 		<- X_list07
	if(v %in% c("ln_upc_per_mod", "ln_num_module")){
		X_new[[sel.retail]][,v]	<- X_new[[sel.retail]][,v] + log(1 + iv.change)
	}else{
		X_new[[sel.retail]][,v]	<- X_new[[sel.retail]][,v] * (1 + iv.change)
	}
	
	out1	<- SimWrapper_fn(omega_deriv = omega_deriv, ln_inc = sim.unq[,"Inc07"], rec = rep(0, nrow(sim.unq)),
							lambda = lambda, param_est = shr.par, base = beta0_base, 
							X_list = X_new, price = price.07, sim.y = sim.base07$y, 
							eps_draw = eps_draw, method = exp.method, ret.sim = TRUE, par.draw = par_draw, share.state = shr.stat[,-beta0_base])
	out1.avg	<- data.frame(out1$Average, lnInc = sim.unq[,"Inc07"], retailer = fmt_name[sel.retail], 
							Var = v, change = iv.change )
	return(out1.avg[1,])
}
use.time	<- proc.time() - pct
cat("Parallel simulation of retail price finishes with", use.time[3]/60, "min.\n")

pct			<- proc.time()
sim.ls08	<- foreach(v = selcol, .errorhandling = "remove", .packages=c("mgcv", "nloptr"), .combine = 'rbind') %dopar% {
	out1.avg	<- data.frame()
	out1.alc	<- array(NA, c(numsim, nrow(sim.unq), R), 
						dimnames = list(iter = 1:numsim, lnInc = sim.unq[,"Inc07"], retail = fmt_name))
	X_new 		<- X_list08
	if(v %in% c("ln_upc_per_mod", "ln_num_module")){
		X_new[[sel.retail]][,v]	<- X_new[[sel.retail]][,v] + log(1 + iv.change)
		X_new[[sel.retail]][,paste("rec*",v,sep="")]	<- X_new[[sel.retail]][,paste("rec*",v,sep="")] + log(1+ iv.change)
	}else{
		X_new[[sel.retail]][,v]	<- X_new[[sel.retail]][,v] * (1 + iv.change)
		X_new[[sel.retail]][,paste("rec*",v,sep="")]	<- X_new[[sel.retail]][,paste("rec*",v,sep="")] * (1 + iv.change)
	}
	
	out1	<- SimWrapper_fn(omega_deriv = omega_deriv1, ln_inc = sim.unq[,"Inc08"], rec = rep(1, nrow(sim.unq)),
							lambda = lambda, param_est = shr.par, base = beta0_base, 
							X_list = X_new, price = price.07, sim.y = sim.base08$y, 
							eps_draw = eps_draw, method = exp.method, ret.sim = TRUE, par.draw = par_draw, share.state = NULL)
	out1.avg	<- data.frame(out1$Average, lnInc = sim.unq[,"Inc08"], retailer = fmt_name[sel.retail], 
							Var = v, change = iv.change )
	return(out1.avg[1,])
}
use.time	<- proc.time() - pct
cat("Parallel simulation of retail price finishes with", use.time[3]/60, "min.\n")

stopCluster(cl)
cat("Stop clustering. \n")

# Check the marginal effects. 
tmp0	<- sim.base08$Average[1,fmt_name]/sim.base08$Average[,"Exp"]
tmp1	<- sim.ls08[,gsub(" ", ".",fmt_name)]/sim.ls08[,"Exp"]
tmp		<- ((tmp1[,sel.retail] - tmp0[sel.retail])/tmp0[sel.retail]) / iv.change
names(tmp)	<- selcol
cat("Elasticity for assortment variables is:\n"); print(tmp); cat("\n")

tmp		<- tmp1[,sel.retail] - tmp0[sel.retail]
names(tmp)	<- selcol
cat("Changes in share:\n"); print(round(tmp*100, 1)); cat("\n")

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

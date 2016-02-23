library(ggplot2)
library(reshape2)
library(Rcpp)
library(RcppArmadillo)
library(maxLik)
library(evd)
library(data.table)
library(doParallel)
library(foreach)
library(chebpol)
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

model_name 	<- "MDCEV_share"
# setwd("~/Documents/Research/Store switching/processed data")
# plot.wd	<- '~/Desktop'
# source("../Exercise/Multiple_discrete_continuous_model/0_Allocation_function.R")
# source("../Exercise/main/share_allocation/ctrfact_sim_functions_v2.r")

setwd("/home/brgordon/ccv103/Exercise/run")
# setwd("/kellogg/users/marketing/2661703/Exercise/run")
# setwd("/sscc/home/c/ccv103/Exercise/run")
run_id		<- 3
plot.wd		<- getwd()
make_plot	<- TRUE
ww			<- 10
ar			<- .6

source("0_Allocation_function.R")
source("ctrfact_sim_functions_v2.r")

# Load estimation data 
ver.date	<- "2016-02-16"
load(paste("estrun_",run_id,"/MDCEV_est_seg",seg_id,"_", ver.date,".rdata",sep=""))
rm(list = c("gamfit", "shr","model_name", "tmpdat"))

# Set simulation parameters
interp.method	<- "spline"					# Spline interpolation or Chebyshev interpolation
exp.method		<- "Utility"				# Utility maximization or finding roots for first order condition
trim.alpha		<- 0.05
numsim			<- numsim1		#<- 1000
draw.par		<- FALSE
sim.omega		<- FALSE
fname			<- paste("ctrfact_seg",seg_id,sep="")
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
lambda 	<- coef(sol.top2)
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
my.change		<- .1
lnInc_08		<- lnInc + log(1 - my.change)
sim.unq$income2008	<- (1 - my.change) * sim.unq$income2007		
sim.unq$Inc07	<- log(sim.unq[,"income2007"])
sim.unq$Inc08	<- log(sim.unq[,"income2008"])
sim.unq			<- sim.unq[order(sim.unq$income2007),]
cat("dim(sim.unq) =", dim(sim.unq), "\n")

#----------------------------#
# Average price in year 2007
tmp 		<- dcast(subset(price_dat, year %in% selyr),	 
					scantrack_market_descr + year + biweek ~ channel_type, value.var = "bsk_price_paid_2004") 
price.07	<- setNames( colMeans(as.matrix(tmp[,4:(3+R)]), na.rm=T), fmt_name)
cat("The average price level in 2007:\n"); print(price.07);cat("\n")

# Average retail attributes in 2007
selcol	<- c("size_index", "ln_upc_per_mod", "ln_num_module","overall_prvt")
X_list07<- setNames( lapply(fmt_name, function(x) colMeans(as.matrix(subset(fmt_attr, channel_type == x & year%in% selyr)[,selcol]))), 
					 fmt_name)
cat("The average retail attributes in 2007:\n"); print(do.call(rbind, X_list07)); cat("\n")

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
cl		<- makeCluster(mycore, type = "FORK")
registerDoParallel(cl)
cat("Register", mycore, "core parallel computing. \n")

# Simulate expenditure and expenditure share. 
pct				<- proc.time()
sim.base07		<- SimWrapper_fn(omega_deriv, ln_inc = sim.unq$Inc07, lambda = lambda, param_est = shr.par, base = beta0_base, 
					X_list = X_list07, price = price.07, eps_draw = eps_draw, method = exp.method, ret.sim = TRUE, par.draw = par_draw)
use.time		<- proc.time() - pct						
cat("2007 Baseline simulation finishes with", use.time[3]/60, "min.\n")

pct				<- proc.time()
sim.base08		<- SimWrapper_fn(omega_deriv, ln_inc = sim.unq$Inc08, lambda = lambda, param_est = shr.par, 
					base = beta0_base, X_list = X_list07, price = price.07, eps_draw = eps_draw, method = exp.method, ret.sim = TRUE, par.draw = par_draw)
use.time		<- proc.time() - pct						
cat("2008 Baseline simulation finishes with", use.time[3]/60, "min.\n")

# Check the baseline simulation 
tmp1	<- sim.base07$Allocation 
tmp2	<- sim.base08$Allocation
tmp		<- tmp2-tmp1
cat("Summary of baseline income elasticity:\n"); print(summary(c(tmp))); cat("\n")

# Plot omega over y
tmp		<- seq(50, 250, 10)
tmpd1	<- data.frame(lnInc = rep(sim.unq$Inc07, each = length(tmp)), y = rep(tmp, length(sim.unq$Inc07)))
tmp1	<- omega_deriv(tmpd1)
tmpd2	<- data.frame(lnInc = rep(sim.unq$Inc08, each = length(tmp)), y = rep(tmp, length(sim.unq$Inc08)))
tmp2	<- omega_deriv(tmpd2)
ggtmp	<- rbind(cbind(tmpd1, omega = tmp1, year = 2007), cbind(tmpd2, omega = tmp2, year = 2008))
ggtmp$Income	<- exp(ggtmp$lnInc)

# Plot the simulation baseline -- expenditure allocation by income 
tmp1	<- melt(sim.base07$Average[,fmt_name])
tmp2	<- melt(sim.base08$Average[,fmt_name])
ggtmp1	<- rbind(data.frame(tmp1, year = 2007), data.frame(tmp2, year = 2008))
names(ggtmp1)	<- c("inc_idx", "Retail", "Exp", "year")
ggtmp1$Inc	<- round(exp(ggtmp1$inc_idx), 0)
ggtmp1$Inc1	<- factor(ggtmp1$Inc, levels = sort(unique(ggtmp1$Inc)))

# Check if utility function at the top level is concave
tmp		<- seq(100, 110, .2)
u		<- matrix(NA, length(tmp), length(lnInc), dimnames = list(y = tmp, lnInc = lnInc))
for(i in 1:length(lnInc)){
	f 	<- function(y, ...){
		newx	<- data.frame(lnInc = lnInc[i], y = y)
		omega_deriv(newx, ...)
	}
	u[,i]<- -exp_fn(tmp, lambda, f, ln_inc = lnInc[i])
}
ggtmp2	<- melt(u, value.name = "utility")

if(make_plot){
	pdf(paste(plot.wd, "/estrun_",run_id,"/graph_checkbase_", fname, ".pdf", sep=""), width = ww/ar, height = ww)
	print(ggplot(ggtmp, aes(y, omega, col = factor(Income))) + geom_line() + 
			facet_wrap(~year) + 
			guides(col = FALSE) + 
			labs(title = "Inclusive value omega(Inc, y)")
		)
	print(ggplot(ggtmp1, aes(Inc1, Exp, fill = Retail)) + geom_bar(stat = "identity", position = "stack") + 
			facet_wrap(~year) + 
			labs(x = "Income", title = "Baseline simulation expenditure allocation") + 
			theme(axis.text.x = element_text(angle = 60))
		)
	print(ggplot(ggtmp2, aes(y, utility)) + geom_point() + geom_line() + facet_wrap(~lnInc, scales = "free"))
	dev.off()
}

save.image(paste("estrun_",run_id, "/", fname, ".rdata",sep=""))

#------------------------#
# Simulation of strategy #
# Assume only grocery stores change retail attributes. 
iv.change		<- .1			# Percentage increase of independent variable

sim.X.ls	<- setNames(vector("list", length(selcol)), selcol)
sim.X.df	<- data.frame()
for(i in 1:length(selcol)){
	pct			<- proc.time()
	var1		<- selcol[i]
	out	<- foreach(sel.retail = fmt_name, .errorhandling = "remove", .packages=c("mgcv", "nloptr")) %do% {
		X_list_new	<- X_list07
		if(var1 %in% c("ln_upc_per_mod", "ln_num_module")){
			X_list_new[[sel.retail]][var1]	<- X_list_new[[sel.retail]][var1] + log(1 + iv.change)
		}else{
			X_list_new[[sel.retail]][var1]	<- X_list_new[[sel.retail]][var1] * (1 + iv.change)
		}
		if(sim.omega){
			out1	<- SimOmega_fn(ln_inc = sim.unq[,"Inc08"], lambda = lambda, param_est = shr.par, 
						base = beta0_base, X_list = X_list_new, price = price.07, 
						lnInc_lv = lnInc_08, y.nodes = y.nodes, eps_draw = eps_draw, method = exp.method, 
						interp.method = interp.method, ret.sim = TRUE, alpha = trim.alpha, par.draw = par_draw)
		}else{
			out1 	<- SimWrapper_fn(omega_deriv = omega_deriv, ln_inc = sim.unq[,"Inc08"], lambda = lambda, param_est = shr.par, 
						base = beta0_base, X_list = X_list_new, price = price.07, 
						eps_draw = eps_draw, method = exp.method, ret.sim = TRUE, par.draw = par_draw)
		}
		
		out1$Average 	<- data.frame(out1$Average, lnInc = sim.unq[,"Inc08"], retailer = as.character(sel.retail), Var = var1)
		return(out1)
	}
	sim.X.ls[[i]]	<- out
	sim.X.df	<- rbind(sim.X.df, do.call(rbind, lapply(out, function(x) x$Average)))
	use.time	<- proc.time() - pct
	cat("Parallel simulation of retail attributes", var1, "finishes with", use.time[3]/60, "min.\n")
}

#-------------------------------#
# Simulation of price reduction # 
pct			<- proc.time()
sim.prc.ls	<- foreach(sel.retail = fmt_name, .errorhandling = "remove", .packages=c("mgcv", "nloptr")) %do% {
	price.new	<- price.07
	price.new[sel.retail]	<- price.07[sel.retail] * (1 - iv.change)
	if(sim.omega){
		out	<- SimOmega_fn(ln_inc = sim.unq[,"Inc08"], lambda = lambda, param_est = shr.par, 
					base = beta0_base, X_list = X_list07, price = price.new, 
					lnInc_lv = lnInc_08, y.nodes = y.nodes, eps_draw = eps_draw, method = exp.method, 
					interp.method = interp.method, ret.sim = TRUE, alpha = trim.alpha, par.draw = par_draw)
	}else{
		out	<- SimWrapper_fn(omega_deriv = omega_deriv, ln_inc = sim.unq[,"Inc08"], lambda = lambda, param_est = shr.par, 
					base = beta0_base, X_list = X_list07, price = price.new, 
					eps_draw = eps_draw, method = exp.method, ret.sim = TRUE, par.draw = par_draw)
	}
	out$Average		<- data.frame(out$Average, lnInc = sim.unq$Inc08, retailer = as.character(sel.retail))
	return(out)
}
use.time	<- proc.time() - pct
cat("Parallel simulation of retail price finishes with", use.time[3]/60, "min.\n")

stopCluster(cl)
cat("Stop clustering. \n")

################
# Save results #
################
rm(list  = c("ar", "args", "cl", "ggtmp", "ggtmp1", "ggtmp2", "i", "lastFuncGrad", "lastFuncParam", "make_plot", "mycore", 
			"myfix", "plot.wd", "s1_index", "sel", "selyr", "price", "sol", "sol.top", "sol.top2", "tmp", "tmp_coef", 
			"tmp_price", "tmp1", "tmp2", "use.time", "ver.date", "var1", "ww", "f", 
			"numnodes", "out", "out1", "pct", "tmpd1", "tmpd2", "tmpdat", "u", "W", "y", "y.nodes",
			"Allocation_constr_fn", "Allocation_fn", "Allocation_nlop_fn", "cheb.1d.basis", "cheb.basis", "chebfun", 
			"exp_fn", "expFOC_fn", "incl_value_fn", "mysplfun", "mytrimfun", "param_assignR", "simExp_fn", "SimOmega_fn", 
			"SimWrapper_fn", "solveExp_fn", "spl2dfun", "uP_fn", "uGrad_fn"))

save.image(paste("estrun_",run_id, "/", fname, ".rdata",sep=""))
# save.image(paste("estrun_",run_id,"/ctrfact_seg",seg_id, "_", Sys.Date(), ".rdata",sep=""))

cat("This program is done.\n")


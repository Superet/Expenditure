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
# source("../Exercise/main/share_allocation/ctrfact_sim_functions.r")

setwd("/home/brgordon/ccv103/Exercise/run")
# setwd("/kellogg/users/marketing/2661703/Exercise/run")
# setwd("/sscc/home/c/ccv103/Exercise/run")
run_id		<- 1
seg_id		<- 1
plot.wd		<- getwd()
make_plot	<- TRUE
ww			<- 6.5
ar			<- .6

source("0_Allocation_function.R")
source("ctrfact_sim_functions.r")

# Load estimation data 
ver.date	<- "2015-11-25"
load(paste("estrun_",run_id,"/MDCEV_est_seg",seg_id,"_", ver.date,".rdata",sep=""))
interp.method	<- "spline"					# Spline interpolation or Chebyshev interpolation
exp.method		<- "Utility"				# Utility maximization or finding roots for first order condition

#############
# Functions #
#############
solveExp_fn	<- function(lambda, omega_fn, ln_inc, method = "FOC", lower = NULL, upper = NULL){
# Solve for the optimal expediture given parameter lambda and a given omega function, at a single value of log(income)	

# lambda	... A vector of utility paraemter of purchase utility
# omega_fn	... Spline function for a given income level
# ln_inc	...	A scalor of log(income)
#==============================================================================#
	if(is.null(lower))	{lower = 50}
	if(is.null(upper))	{upper = 250}
	init.interval	<- c(lower, upper)

	switch(method, 
		FOC = {
			sol	<- try(uniroot(expFOC_fn, c(20, 100), lambda = lambda, omega_fn = omega_fn, ln_inc = ln_inc), silent = TRUE)
			# The root solution is sensitive to the solution range. 
			# So we search for wider range gradually
			if(class(sol) == "try-error"){
				maxpt	<- seq(150, 600, by = 50)
				i 		<- 1
				while(i <= length(maxpt) & class(sol) == "try-error"){
					sol	<- try(uniroot(expFOC_fn, c(1, maxpt[i]), lambda = lambda, omega_fn = omega_fn, ln_inc = ln_inc), silent = TRUE)
					i	<- i + 1
				}
				if(class(sol) == "try-error"){
					out <- 0
					cat("Expenditure solution was not found at ln_inc =", ln_inc, ".\n")
				}else{
					out	<- sol
				}
			}else{
				out	<- sol
			}
			out$obj	<- exp_fn(sol$root, lambda, omega_fn, ln_inc)
		}, 
		
		Utility = {
			sol	<- try(optimize(exp_fn, interval = init.interval, lambda = lambda, omega_fn = omega_fn, ln_inc = ln_inc), silent = TRUE)
			check	<- class(sol) == "try-error"
			
			if(class(sol) == "try-error"){
				interval	<- c(1, 500)
				sol	<- try(optimize(exp_fn, interval = interval, lambda = lambda, omega_fn = omega_fn, ln_inc = ln_inc), silent = TRUE)
				if(class(sol) == "try-error"){
					out <- 0
					cat("Expenditure solution was not found at ln_inc =", ln_inc, ".\n")
				}else{
					out <- sol		
				}
			}
			# else if (sol$minimum <= init.interval[1] + .1  | sol$minimum >= init.interval[2] - .1){
			# 		interval	<- c(1, 500)
			# 		sol	<- try(optimize(exp_fn, interval = interval, lambda = lambda, omega_fn = omega_fn, ln_inc = ln_inc), silent = TRUE)
			# 		if(class(sol) == "try-error"){
			# 			out <- 0
			# 			cat("Expenditure solution was not found at ln_inc =", ln_inc, ".\n")
			# 		}else{
			# 			out <- sol		
			# 		}
			# 	}
			else{
				out	<- sol
			}
		}
	)
	
	return(out)
}

simExp_fn	<- function(lambda, omega_fn_ls, ln_inc, method = "FOC", lower = NULL, upper = NULL){
# A function to simulate optimal expenditure for a vector log(income)

# lambda		... A vector of utility parameters for purchase utility
# omega_fn_ls	...	A list of spline function, omega_fn, with each element corresponding to an income level 
# ln_inc		... A vector of log(income) as input data
#==============================================================================#
	out		<- vector("list", length = length(ln_inc))		#rep(NA, length(ln_inc))		
	lvinc	<- unique(ln_inc)				# Levels of discrete income bin
	for(i in 1:length(lvinc)){
		sel <- which(ln_inc == lvinc[i])
		sel1<- which(names(omega_fn_ls) == lvinc[i])
		out[[sel]]	<- solveExp_fn(lambda, omega_fn_ls[[sel1]], lvinc[i], method)
	}
	return(out)
}


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
numsim 		<- 400
eps_draw	<- matrix(rgev(numsim*R), numsim, R)

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
mycore 	<- 2
cl		<- makeCluster(mycore, type = "FORK")
registerDoParallel(cl)
cat("Register", mycore, "core parallel computing. \n")

# --------------------#
# Simulate inclusive values 
omega_parallel	<- function(eps_draw, lnInc, y.nodes, X_list, price){
	numnodes<- length(y.nodes)
	out		<- matrix(NA, length(lnInc), numnodes)
	for(j in 1:length(lnInc)){
		tmpX_list1	<- lapply(X_list, function(x) rep(1, numnodes) %*% t(c(x, x*lnInc[j])))
		tmpsol 		<- incl_value_fn(param_est=shr.par, base= beta0_base, X_list=tmpX_list1, y=y.nodes, Q=Inf, price=rep(1, numnodes) %*% t(price), 
							R=R, Ra=R, qz_cons = 0, exp_outside = FALSE, quant_outside = FALSE, eps_draw = rep(1, numnodes) %*% t(eps_draw) )
		out[j,]		<- tmpsol$omega
	}
	return(out)
}

pct		<- proc.time()
omega_draw_ls		<- foreach(i = 1:numsim) %dopar% {
	omega_parallel(eps_draw = eps_draw[i,], lnInc = lnInc, y.nodes = y.nodes, X_list = X_list07, price = price.07)
}
use.time	<- proc.time() - pct
cat("Simulation of omega with random draws finishes with", use.time[3]/60, "min.\n")

#---------------#
# Interpolation #
omega.draw	<- array(NA, c(numsim, length(lnInc), numnodes), dimnames = list(NULL, lnInc, y.nodes))
for(i in 1:numsim){
	omega.draw[i,,]	<- omega_draw_ls[[i]]
}
omega.07	<- apply(omega.draw, c(2, 3), mean, trim = .05)

# Interpolation functions
if(interp.method == "spline"){
	omega.ls.itpl <- lapply(1:length(lnInc), function(i) splinefun(y.nodes, omega.07[i,], method = "natural"))
}else{
	omega.ls.itpl <- lapply(1:length(lnInc), function(i) chebfun(x = y.nodes, y = omega.07[i,], interval = y.interval))
}
names(omega.ls.itpl)	<- lnInc

# ------------------#
# Smoothing splines #
mysplfun	<- function(x, y){
	fit	<- smooth.spline(x, y)
	f	<- function(newx){
		return(predict(fit, newx)$y)
	}
	return(f)
}

mytrimfun	<- function(x, alpha = .05){
	sort.x	<- sort(x)
	dropn	<- floor(alpha * length(x))
	drop.idx<- c(1:dropn, (length(x) - dropn + 1):length(x))
	return(sort.x[-drop.idx])
}

tmp	<- apply(omega.draw, c(2,3), mytrimfun)
# omega.ls.smspl	<- lapply(1:length(lnInc), function(i) mysplfun(x = rep(y.nodes, each = numsim), y = as.vector(omega.draw[,i,])))
omega.ls.smspl	<- lapply(1:length(lnInc), function(i) 
							mysplfun(x = rep(y.nodes, each = dim(tmp)[1]), y = as.vector(tmp[,i,])))
names(omega.ls.smspl)	<- lnInc

#--------------------------------------------#
# Check the baseline simulation 
tmp		<- seq(50, 250, 10)
tmp1	<- sapply(omega.ls.itpl, function(f) f(tmp))
tmp2	<- sapply(omega.ls.smspl, function(f) f(tmp))
ggtmp	<- rbind(data.frame(melt(tmp1), method = "interpolation"), data.frame(melt(tmp2), method = "smoothspl"))
names(ggtmp)	<- c("yidx", "lnInc", "omega", "method")
ggtmp$y			<- tmp[ggtmp$yidx]
ggtmp$Income	<- exp(ggtmp$lnInc)

# Simulate expenditure and expenditure share #
pct				<- proc.time()
sim.itpl		<- simExp_fn(lambda, omega_fn_ls = omega.ls.itpl, ln_inc = lnInc, method = "Utility")
sim.smspl		<- simExp_fn(lambda, omega_fn_ls = omega.ls.smspl, ln_inc = lnInc, method = "Utility")
use.time		<- proc.time() - pct						
cat("2007 Baseline simulation finishes with", use.time[3]/60, "min.\n")

tmp1	<- sapply(sim.itpl, function(x) x$minimum) 
tmp2	<- sapply(sim.smspl, function(x) x$minimum)
tmp	<- seq(-20, 20, 1)
ggtmp1	<- data.frame()
for(i in 1:length(lnInc)){
	sel	<- which.min(abs(as.numeric(names(omega.ls.itpl)) - lnInc[i]))
	# u	<- (lambda[1] + lambda[2]*tmp.inc[i]) * omega_deriv[[sel]](tmp.opt[i] + tmp)
	u	<- sapply(tmp1[i] + tmp, function(y) -exp_fn(y, lambda, omega_fn = omega.ls.itpl[[sel]], ln_inc = lnInc[i]))
	ggtmp1	<- rbind(ggtmp1, data.frame(lnInc = lnInc[i], y = tmp1[i]+tmp, u = u, sol = 1*(tmp==0), method="interpolation"))
	u	<- sapply(tmp2[i] + tmp, function(y) -exp_fn(y, lambda, omega_fn = omega.ls.smspl[[sel]], ln_inc = lnInc[i]))
	ggtmp1	<- rbind(ggtmp1, data.frame(lnInc = lnInc[i], y = tmp2[i]+tmp, u = u, sol = 1*(tmp==0), method="smoothingspl"))
}
names(ggtmp1)	<- c("lnInc", "y", "u", "sol", "method")

if(make_plot){
	pdf(paste(plot.wd, "/estrun_",run_id,"/graph_ctrfact_test_seg", seg_id, ".pdf", sep=""), width = ww, height = ww)
	print(ggplot(ggtmp, aes(y, omega, col = factor(Income))) + geom_line() + 
			facet_wrap(~method) + 
			guides(col = FALSE) + 
			labs(title = "Inclusive value omega(Inc, y)")
		)
	print(	ggplot(ggtmp1, aes(y, u, linetype = method, shape = method)) + geom_point(aes(col=factor(sol))) + geom_line() + 
						facet_wrap(~lnInc, scales = "free") + 
						scale_color_manual(values = c("black", "red"))
		)
	dev.off()
}

stopCluster(cl)
cat("Stop clustering. \n")

################
# Save results #
################
save.image(paste("estrun_",run_id,"/test_smoothingspl.rdata", sep=""))

cat("This program is done.\n")


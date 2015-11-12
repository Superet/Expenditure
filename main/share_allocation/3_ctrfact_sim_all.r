library(ggplot2)
library(reshape2)
library(Rcpp)
library(RcppArmadillo)
library(maxLik)
library(evd)
library(data.table)
library(doParallel)
library(foreach)
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
# sourceCpp(paste("../Exercise/Multiple_discrete_continuous_model/", model_name, ".cpp", sep=""))
# source("../Exercise/Multiple_discrete_continuous_model/0_Allocation_function.R")

setwd("/home/brgordon/ccv103/Exercise/run")
# setwd("/kellogg/users/marketing/2661703/Exercise/run")
# setwd("/sscc/home/c/ccv103/Exercise/run")
run_id		<- 1
# seg_id		<- 1
plot.wd		<- getwd()
make_plot	<- TRUE
ww			<- 6.5
ar			<- .6

sourceCpp(paste(model_name, ".cpp", sep=""))
source("0_Allocation_function.R")

# Load estimation data 
load(paste("estrun_",run_id,"/MDCEV_est_seg",seg_id,".rdata",sep=""))

#############
# Functions #
#############
exp_fn <- function(y, lambda, omega_fn, ln_inc){
# Utility function expenditure for one single observation

# y			... Expenditure scalor
# lambda	...	A vector of utility parameters of purchase utility
# omega_fn	...	Spline function for a given income level within an income group
# ln_inc	... A scalor of log(income)
#==============================================================================#
	out	<- (lambda[1] + lambda[2] * ln_inc) * omega_fn(y) + (exp(ln_inc) - y)
	return(-out)
}

solveExp_fn	<- function(lambda, omega_fn, ln_inc){
# Solve for the optimal expediture given parameter lambda and a given omega function, at a single value of log(income)	

# lambda	... A vector of utility paraemter of purchase utility
# omega_fn	... Spline function for a given income level
# ln_inc	...	A scalor of log(income)
#==============================================================================#
	init.interval	<- c(10, 200)
	sol	<- try(optimize(exp_fn, interval = init.interval, lambda = lambda, omega_fn = omega_fn, ln_inc = ln_inc), silent = TRUE)
	if(class(sol) == "try-error" | sol$minimum <= init.interval[1] + .1  | sol$minimum >= init.interval[2] - .1){
		interval	<- c(1, 500)
		sol	<- try(optimize(exp_fn, interval = interval, lambda = lambda, omega_fn = omega_fn, ln_inc = ln_inc), silent = TRUE)
		if(class(sol) == "try-error"){
			out <- 0
			cat("Expenditure solution was not found at ln_inc =", ln_inc, ".\n")
		}else{
			out <- sol$minimum			
		}
	}else{
		out	<- sol$minimum
	}
	return(out)
}

simExp_fn	<- function(lambda, omega_fn_ls, ln_inc){
# A function to simulate optimal expenditure for a vector log(income)

# lambda		... A vector of utility parameters for purchase utility
# omega_fn_ls	...	A list of spline function, omega_fn, with each element corresponding to an income level 
# ln_inc		... A vector of log(income) as input data
#==============================================================================#
	out		<- rep(NA, length(ln_inc))		
	lvinc	<- unique(ln_inc)				# Levels of discrete income bin
	for(i in 1:length(lvinc)){
		sel <- ln_inc == lvinc[i]
		sel1<- which(names(omega_fn_ls) == lvinc[i])
		out[sel]	<- solveExp_fn(lambda, omega_fn_ls[[sel1]], lvinc[i])
	}
	return(out)
}

SimWrapper_fn	<- function(omega_fn_ls, ln_inc, lambda, param_est, base, X_list, price, eps_draw, sim.y = NULL){
# This function simulates expenditure and share, having known omega_fn_ls. 

# omega_fn_ls	...	A list of spline function, omega_fn, with each element corresponding to an income level 
# ln_inc		...	A vector of log(income) as simulation input data
# lambda		... A vector of utility parameters for purchase utility
# param_est		... A vector of parameters in the conditional allocation model
# base			... Integer index of which retailer has intercept of 0 in the conditional allocation model
# X_list			...	A list of retail attributes
# price			...	A vector of price 
# esp_draw		... A matrix of random draw
# sim.y			... If not NULL, then expenditure is given and only expenditure share is simulated. 
#=========================================================================================================#
	
	if(is.null(sim.y)){
		y		<- rep(NA, length(ln_inc))
	}else{
		y 		<- sim.y
	}
	numsim	<- nrow(eps_draw)
	e		<- array(NA, c(numsim, length(ln_inc), R))
	omega	<- matrix(NA, numsim, length(ln_inc))
	omega1	<- matrix(NA, numsim, length(ln_inc))
	lvinc	<- sort(unique(ln_inc))
	for(i in 1:length(lvinc)){
		sel <- ln_inc == lvinc[i]
		if(sum(sel) > 0){
			sel1<- which(abs(as.numeric(names(omega_fn_ls)) - lvinc[i]) < 1e-6)
			if(is.null(sim.y)){
				y[sel]	<- solveExp_fn(lambda, omega_fn_ls[[sel1]], lvinc[i])
			}
			tmpX_list1	<- lapply(X_list, function(x) rep(1, sum(sel)) %*% t(c(x, x*lvinc[i])))
			for(j in 1:numsim){
				tmpsol 		<- incl_value_fn(param_est, base, X_list=tmpX_list1, y=y[sel], Q=Inf, price=rep(1, sum(sel)) %*% t(price), 
							R=R, Ra=R, qz_cons = 0, exp_outside = FALSE, quant_outside = FALSE, 
							eps_draw = rep(1, sum(sel)) %*% t(eps_draw[j,]) )
				omega[j,sel] 	<- tmpsol$omega			
				omega1[j,sel]	<- omega_fn_ls[[sel1]](y[sel])
				e[j,sel,]		<- tmpsol$e
			}			
		}
	}
	out	<- cbind(y, colMeans(omega, na.rm = T), colMeans(omega1, na.rm = T), apply(e, c(2,3), mean, na.rm= T)/y)
	colnames(out)	<- c("Exp", "omega", "spl.omega", fmt_name)
	return(out)
}

SimOmega_fn	<- function(ln_inc, lambda, param_est, base, X_list, price, lnInc_lv, y.nodes, eps_draw){
# This function simulate inclusive values first, and then simulate expenditure and expenditure share; 

# ln_inc		...	A vector of log(income) as simulation input data
# lambda		... A vector of utility parameters for purchase utility
# param_est		... A vector of parameters in the conditional allocation model
# base			... Integer index of which retailer has intercept of 0 in the conditional allocation model
# X_list			...	A list of retail attributes
# price			...	A vector of price 
# lnInc_lv		... A vector of discrete log(income) level. 
# y.nodes			... A vector of nodes of expenditure to simulate omega
# esp_draw		... A matrix of random draw
#===================================================================================================#
	numsim		<- nrow(eps_draw)
	numnodes	<- length(y.nodes)
	
	# Simulate inclusive values under the new retail attributes
	omega_new	<- array(NA, c(numsim, length(lnInc_lv), numnodes), dimnames = list(NULL, lnInc_lv, y.nodes))
	for(i in 1:numsim){
		for(j in 1:length(lnInc_lv)){
			tmpX_list	<- lapply(X_list, function(x) rep(1, numnodes) %*% t(c(x, x*lnInc_lv[j])))
			tmpsol 		<- incl_value_fn(param_est=shr.par, base= beta0_base, X_list=tmpX_list, y=y.nodes, Q=Inf, price=rep(1,numnodes) %*% t(price), 
						R=R, Ra=R, qz_cons = 0, exp_outside = FALSE, quant_outside = FALSE, eps_draw = rep(1, numnodes) %*% t(eps_draw[i,]) )
			omega_new[i,j,] <- tmpsol$omega
		}	
	}	

	# Spline function of omega
	omega.ls.new <- lapply(1:length(lnInc_lv), function(i) splinefun(y.nodes, colMeans(omega_new[,i,], na.rm = T), method = "natural"))
	names(omega.ls.new)	<- lnInc_lv

	# Simulate expenditure and share
	out	<- SimWrapper_fn(omega_fn_ls = omega.ls.new, ln_inc = ln_inc, lambda = lambda, param_est = param_est, 
						base = base, X_list = X_list, price = price, eps_draw = eps_draw)
	return(out)
}

###############################
# Prepare simulation elements # 
###############################
# Data required: parameters, income level, price, retail attributes, and random draws
# For each simulation scenario, if income, price or retail attributes change, we need to re-simulate inclusive values.

# Set simulation parameters
lambda 	<- coef(rd.fit2)
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
numsim 		<- 100
eps_draw	<- matrix(rgev(numsim*R), numsim, R)

##############
# Simulation #
##############
# Register parallel computing
mycore 	<- 2
cl		<- makeCluster(mycore, type = "FORK")
registerDoParallel(cl)
cat("Register", mycore, "core parallel computing. \n")

# --------------------#
# Baseline simulation #
omega.07	<- array(NA, c(numsim, length(lnInc), numnodes), dimnames = list(NULL, lnInc, y.nodes))
omega.08	<- array(NA, c(numsim, length(lnInc), numnodes), dimnames = list(NULL, lnInc_08, y.nodes))

# Simulate inclusive values 
pct			<- proc.time()
omega.07	<- foreach(i = 1:length(lnInc), .combine = rbind) %dopar%{
	out		<- matrix(NA, numsim, numnodes, dimnames = list(NULL, y = y.nodes))
	for(j in 1:numsim){
		tmpX_list1	<- lapply(X_list07, function(x) rep(1, numnodes) %*% t(c(x, x*lnInc[i])))
		tmpsol 		<- incl_value_fn(param_est=shr.par, base= beta0_base, X_list=tmpX_list1, y=y.nodes, Q=Inf, price=rep(1,numnodes) %*% t(price.07), 
					R=R, Ra=R, qz_cons = 0, exp_outside = FALSE, quant_outside = FALSE, eps_draw = rep(1, numnodes) %*% t(eps_draw[j,]) )
		out[j,]		<- tmpsol$omega
	}
	return(colMeans(out, na.rm = T))			
}
rownames(omega.07)	<- lnInc

omega.08	<- foreach(i = 1:length(lnInc_08), .combine = rbind) %dopar%{
	out		<- matrix(NA, numsim, numnodes, dimnames = list(NULL, y = y.nodes))
	for(j in 1:numsim){
		tmpX_list1	<- lapply(X_list07, function(x) rep(1, numnodes) %*% t(c(x, x*lnInc_08[i])))
		tmpsol 		<- incl_value_fn(param_est=shr.par, base= beta0_base, X_list=tmpX_list1, y=y.nodes, Q=Inf, price=rep(1,numnodes) %*% t(price.07), 
					R=R, Ra=R, qz_cons = 0, exp_outside = FALSE, quant_outside = FALSE, eps_draw = rep(1, numnodes) %*% t(eps_draw[j,]) )		
		out[j,]		<- tmpsol$omega
	}
	return(colMeans(out, na.rm = T))			
}
rownames(omega.08)	<- lnInc_08
use.time	<- proc.time() - pct
cat("Simulation of omega with random draws finishes with", use.time[3]/60, "min.\n")

# Spline functions
omega.ls.07 <- lapply(1:length(lnInc), function(i) splinefun(y.nodes, omega.07[i,], method = "natural"))
omega.ls.08 <- lapply(1:length(lnInc), function(i) splinefun(y.nodes, omega.08[i,], method = "natural"))
names(omega.ls.07)	<- lnInc
names(omega.ls.08)	<- lnInc_08

# Simulate expenditure and expenditure share. 
pct				<- proc.time()
sim.base07		<- SimWrapper_fn(omega.ls.07, ln_inc = sim.unq$Inc07, lambda = lambda, param_est = shr.par, 
						base = beta0_base, X_list = X_list07, price = price.07, eps_draw = eps_draw)
use.time		<- proc.time() - pct						
cat("2007 Baseline simulation finishes with", use.time[3]/60, "min.\n")

pct				<- proc.time()
sim.base08		<- SimWrapper_fn(omega_fn_ls = omega.ls.08, ln_inc = sim.unq$Inc08, lambda = lambda, param_est = shr.par, 
						base = beta0_base, X_list = X_list07, price = price.07, eps_draw = eps_draw)
use.time		<- proc.time() - pct						
cat("2008 Baseline simulation finishes with", use.time[3]/60, "min.\n")

# Check the baseline simulation 
tmp		<- seq(10, 300, 20)
tmp1	<- sapply(omega.ls.07, function(f) f(tmp))
tmp2	<- sapply(omega.ls.08, function(f) f(tmp))
ggtmp	<- rbind(data.frame(melt(tmp1), year = 2007), data.frame(melt(tmp2), year = 2008))
names(ggtmp)	<- c("yidx", "lnInc", "omega", "year")
ggtmp$y			<- tmp[ggtmp$yidx]
ggtmp$Income	<- exp(ggtmp$lnInc)

tmp1	<- melt(sim.base07[,"Exp"] * sim.base07[,fmt_name])
tmp2	<- melt(sim.base08[,"Exp"] * sim.base08[,fmt_name])
ggtmp1	<- rbind(data.frame(tmp1, year = 2007), data.frame(tmp2, year = 2008))
names(ggtmp1)	<- c("inc_idx", "Retail", "Exp", "year")
ggtmp1$Inc	<- sim.unq[ggtmp1$inc_idx,"income2007"]
sel		<- ggtmp1$year == 2008
ggtmp1[sel,"Inc"]	<- sim.unq[ggtmp1[sel,"inc_idx"], "income2008"]

if(make_plot){
	pdf(paste(plot.wd, "/graph_ctrfact_checkbase_seg", seg_id, ".pdf", sep=""), width = ww, height = ww)
	print(ggplot(ggtmp, aes(y, omega, col = factor(Income))) + geom_line() + 
			facet_wrap(~year) + 
			guides(col = FALSE) + 
			labs(title = "Inclusive value omega(Inc, y)")
		)
	print(ggplot(ggtmp1, aes(as.character(Inc), Exp, fill = Retail)) + geom_bar(stat = "identity", position = "stack") + 
			facet_wrap(~year) + 
			labs(x = "Income", title = "Baseline simulation expenditure allocation") + 
			theme(axis.text.x = element_text(angle = 60))
		)
	dev.off()
}

#------------------------#
# Simulation of strategy #
# Assume only grocery stores change retail attributes. 
iv.change		<- .1			# Percentage increase of independent variable

pct			<- proc.time()
sim.X.ls	<- foreach(var = selcol) %dopar% {
	# Loop over retailers 
	out 	<- foreach(i = 1:R) %do% {
		sel.retail	<- fmt_name[i]
		X_list_new	<- X_list07
		if(var %in% c("ln_upc_per_mod", "ln_num_module")){
			X_list_new[[sel.retail]][var]	<- X_list_new[[sel.retail]][var] + log(1 + iv.change)
		}else{
			X_list_new[[sel.retail]][var]	<- X_list_new[[sel.retail]][var] * (1 + iv.change)
		}
		return(SimOmega_fn(ln_inc = sim.unq$Inc08, lambda = lambda, param_est = shr.par, 
					base = beta0_base, X_list = X_list_new, price = price.07, 
					lnInc_lv = lnInc_08, y.nodes = y.nodes, eps_draw = eps_draw) )
	}
	names(out)	<- fmt_name
	return(out)
}
names(sim.X.ls)	<- selcol
use.time	<- proc.time() - pct
cat("Parallel simulation of retail attributes finishes with", use.time[3]/60, "min.\n")

#-------------------------------#
# Simulation of price reduction # 
pct			<- proc.time()
sim.prc.ls	<- foreach(sel.retail = fmt_name) %dopar% {
	price.new	<- price.07
	price.new[sel.retail]	<- price.07[sel.retail] * (1 - iv.change)
	out		<- SimOmega_fn(ln_inc = sim.unq$Inc08, lambda = lambda, param_est = shr.par, 
				base = beta0_base, X_list = X_list07, price = price.new, 
				lnInc_lv = lnInc_08, y.nodes = y.nodes, eps_draw = eps_draw)
	return(out)
}
use.time	<- proc.time() - pct
names(sim.prc.ls)	<- fmt_name
cat("Parallel simulation of retail price finishes with", use.time[3]/60, "min.\n")

stopCluster(cl)
cat("Stop clustering. \n")

################
# Save results #
################
save.image(paste("estrun_",run_id,"/ctrfact_seg",seg_id,".rdata",sep=""))

cat("This program is done.\n")


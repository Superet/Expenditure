library(ggplot2)
library(reshape2)
library(Rcpp)
library(RcppArmadillo)
library(maxLik)
library(evd)
library(data.table)
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

# setwd("/home/brgordon/ccv103/Exercise/run")
setwd("/kellogg/users/marketing/2661703/Exercise/run")
# setwd("/sscc/home/c/ccv103/Exercise/run")
run_id		<- 1
# seg_id		<- 1
make_plot	<- FALSE

sourceCpp(paste(model_name, ".cpp", sep=""))
source("0_Allocation_function.R")

# Load estimation data 
load(paste("estrun_",run_id,"/MDCEV_est_seg",seg_id,".rdata",sep=""))

#############
# Functions #
#############
expFOC_fn <- function(y, lambda, omega_fn, ln_inc){
# First order condition of optimal expenditure for one single observation
# y			... Expenditure 
# lambda	...	Utility parameter of purchase utility
# omega_fn	...	Spline function for a given income level within an income group
	foc	<- (lambda[1] + lambda[2] * ln_inc) * omega_fn(y, deriv = 1) - 1
	return(foc)
}

solveExp_fn	<- function(lambda, omega_fn, ln_inc){
# Solve for the optimal expediture given parameter lambda and a given omega function	
	sol	<- try(uniroot(expFOC_fn, c(20, 100), lambda = lambda, omega_fn = omega_fn, ln_inc = ln_inc), silent = TRUE)
	if(class(sol) == "try-error"){
		maxpt	<- seq(200, 700, by = 100)
		i 		<- 1
		while(i <= length(maxpt) & class(sol) == "try-error"){
			sol	<- try(uniroot(expFOC_fn, c(1, maxpt[i]), lambda = lambda, omega_fn = omega_fn, ln_inc = ln_inc), silent = TRUE)
			i	<- i + 1
		}
		if(class(sol) == "try-error"){
			out <- 0
			cat("Expenditure solution was not found at ln_inc =", ln_inc, ".\n")
		}else{
			out	<- sol$root
		}
	}else{
		out	<- sol$root
	}
	return(out)
}

simExp_fn	<- function(lambda, omega_deriv, ln_inc){
# A function to simulate optimal expenditure
# lambda		... Utility parameter for purchase utility
# omega_deriv	...	A list of spline function, with each element corresponding to an income level 
# ln_inc		... A vector of log(income) as input data
	out		<- rep(NA, length(ln_inc))
	lvinc	<- unique(ln_inc)
	for(i in 1:length(lvinc)){
		sel <- ln_inc == lvinc[i]
		sel1<- which(names(omega_deriv) == lvinc[i])
		out[sel]	<- solveExp_fn(lambda, omega_deriv[[sel1]], lvinc[i])
	}
	return(out)
}

SimWrapper_fn	<- function(omega_deriv, ln_inc, lambda, param_est, base, X_list, price, eps_draw, sim.y = NULL){
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
			sel1<- which(abs(as.numeric(names(omega_deriv)) - lvinc[i]) < 1e-6)
			if(is.null(sim.y)){
				y[sel]	<- solveExp_fn(lambda, omega_deriv[[sel1]], lvinc[i])
			}
			tmpX_list1	<- lapply(X_list, function(x) rep(1, sum(sel)) %*% t(c(x, x*lvinc[i])))
			for(j in 1:numsim){
				tmpsol 		<- incl_value_fn(param_est, base, X_list=tmpX_list1, y=y[sel], Q=Inf, price=rep(1, sum(sel)) %*% t(price), 
							R=R, Ra=R, qz_cons = 0, exp_outside = FALSE, quant_outside = FALSE, 
							eps_draw = rep(1, sum(sel)) %*% t(eps_draw[j,]) )
				omega[j,sel] 	<- tmpsol$omega			
				omega1[j,sel]	<- omega_deriv[[sel1]](y[sel])
				e[j,sel,]		<- tmpsol$e
			}			
		}
	}
	out	<- cbind(y, colMeans(omega, na.rm = T), colMeans(omega1, na.rm = T), apply(e, c(2,3), mean, na.rm= T)/y)
	colnames(out)	<- c("Exp", "omega", "spl.omega", fmt_name)
	return(out)
}

#################################
# Prepare simulation input data # 
#################################
# Focus on year 07-08
selyr	<- 2007
tmp		<- data.table(subset(mydata, year %in% selyr))
tmp		<- tmp[,list(income = unique(income_midvalue)), by = list(household_code, year)]
sim.data<- data.frame(tmp)[,c("household_code","income")]
sim.data.all	<- sim.data
sim.data		<- data.frame(income2007 = unique(sim.data[,-1]))

# Counterfactual scenario: income is lower by 10%. 
my.change		<- .1
lnInc_new		<- lnInc + log(1 - my.change)
sim.data$income2008	<- (1 - my.change) * sim.data$income2007		
sim.data$Inc07	<- log(sim.data[,"income2007"])
sim.data$Inc08	<- log(sim.data[,"income2008"])
cat("dim(sim.data) =", dim(sim.data), "\n")

# Average price in year 2007
tmp 		<- dcast(subset(price_dat, year %in% selyr),	 
					scantrack_market_descr + year + biweek ~ channel_type, value.var = "bsk_price_paid_2004") 
tmp_price	<- colMeans(as.matrix(tmp[,4:(3+R)]), na.rm=T) 

# Average retail attributes in 2007
selcol	<- c("size_index", "ln_upc_per_mod", "ln_num_module","overall_prvt")
tmpX_list	<- lapply(fmt_name, function(x) colMeans(as.matrix(subset(fmt_attr, channel_type == x & year%in% selyr)[,selcol])))		

# Set simulation parameters
lambda 	<- coef(sol1)
shr.par	<- coef(sol)

#---------------------------------------------------------------- # 
# Compute the inclusive value under income levels with random draws
set.seed(666)
numsim 		<- 100
eps_draw	<- matrix(rgev(numsim*R), numsim, R)
omega_draw	<- array(NA, c(numsim, length(lnInc), numnodes), dimnames = list(NULL, lnInc, tmpy))
omega_draw_new	<- array(NA, c(numsim, length(lnInc), numnodes), dimnames = list(NULL, lnInc_new, tmpy))

pct			<- proc.time()
for(i in 1:numsim){
	for(j in 1:length(lnInc)){
		tmpX_list1	<- lapply(tmpX_list, function(x) rep(1, numnodes) %*% t(c(x, x*lnInc[j])))
		tmpsol 		<- incl_value_fn(param_est=shr.par, base= beta0_base, X_list=tmpX_list1, y=tmpy, Q=Inf, price=rep(1,numnodes) %*% t(tmp_price), 
					R=R, Ra=R, qz_cons = 0, exp_outside = FALSE, quant_outside = FALSE, eps_draw = rep(1, numnodes) %*% t(eps_draw[i,]) )
		omega_draw[i,j,] <- tmpsol$omega
		
		tmpX_list1	<- lapply(tmpX_list, function(x) rep(1, numnodes) %*% t(c(x, x*lnInc_new[j])))
		tmpsol 		<- incl_value_fn(param_est=shr.par, base= beta0_base, X_list=tmpX_list1, y=tmpy, Q=Inf, price=rep(1,numnodes) %*% t(tmp_price), 
					R=R, Ra=R, qz_cons = 0, exp_outside = FALSE, quant_outside = FALSE, eps_draw = rep(1, numnodes) %*% t(eps_draw[i,]) )
		omega_draw_new[i,j,] <- tmpsol$omega
	}	
}
use.time	<- proc.time() - pct
cat("Simulation of omega with random draws finishes with", use.time[3]/60, "min.\n")

# Spline functions
omega_deriv <- lapply(1:length(lnInc), function(i) splinefun(tmpy, colMeans(omega_draw[,i,], na.rm = T), method = "natural"))
omega_deriv_new <- lapply(1:length(lnInc), function(i) splinefun(tmpy, colMeans(omega_draw_new[,i,], na.rm = T), method = "natural"))
names(omega_deriv)		<- lnInc
names(omega_deriv_new)	<- lnInc_new

# ------------------- #
# Baseline simulation #
pct				<- proc.time()
sim.base07		<- SimWrapper_fn(omega_deriv, ln_inc = sim.data$Inc07, lambda = lambda, param_est = shr.par, 
						base = beta0_base, X_list = tmpX_list, price = tmp_price, eps_draw = eps_draw)
use.time		<- proc.time() - pct						
cat("2007 Baseline simulation finishes with", use.time[3]/60, "min.\n")

pct				<- proc.time()
sim.base08		<- SimWrapper_fn(omega_deriv = omega_deriv_new, ln_inc = sim.data$Inc08, lambda = lambda, param_est = shr.par, 
						base = beta0_base, X_list = tmpX_list, price = tmp_price, eps_draw = eps_draw)
use.time		<- proc.time() - pct						
cat("2008 Baseline simulation finishes with", use.time[3]/60, "min.\n")

# Fix expenditure level in 2008 at the level of 2007
pct				<- proc.time()
sim.fixy08		<- SimWrapper_fn(omega_deriv_new, ln_inc = sim.data$Inc08, lambda = lambda, param_est = shr.par, 
						base = beta0_base, X_list = tmpX_list, price = tmp_price, sim.y = sim.base07[,"Exp"], eps_draw = eps_draw)
use.time		<- proc.time() - pct						
cat("Counterfactual simulation with fixed expenditure level finishes with", use.time[3]/60, "min.\n")

# Keep prefrenece the same
pct				<- proc.time()
sim.fixp08		<- SimWrapper_fn(omega_deriv, ln_inc = sim.data$Inc07, lambda = lambda, param_est = shr.par, 
						base = beta0_base, X_list = tmpX_list, price = tmp_price, sim.y = sim.base08[,"Exp"], eps_draw = eps_draw)
use.time		<- proc.time() - pct						
cat("Counterfactual simulation with fixed preference finishes with", use.time[3]/60, "min.\n")

# Compare the population mean
tmp		<- 1:nrow(sim.data)
names(tmp)	<- sim.data$income2007
sel		<- tmp[as.character(sim.data.all$income)]
# sel		<- sim.data[,"income2008"] < sim.data[,"income2007"]
tmp1 <- cbind(	Base07 	= colMeans(sim.base07[sel,], na.rm = T), 
				Base08 	= colMeans(sim.base08[sel,], na.rm = T), 
				FixExp	= colMeans(sim.fixy08[sel,], na.rm = T), 
				FixPref	= colMeans(sim.fixp08[sel,], na.rm = T))						
cat("Mean of simuulation:\n"); print(round(tmp1, 4)); cat("\n")

# Compare the mean difference
tmp2 <- cbind(Base 		= colMeans(sim.base08[sel,] - sim.base07[sel,], na.rm = T), 
			  FixExp 	= colMeans(sim.fixy08[sel,] - sim.base07[sel,], na.rm = T), 
			  FixPref	= colMeans(sim.fixp08[sel,] - sim.base07[sel,], na.rm = T))						
cat("Mean of simuulation difference:\n"); print(round(tmp2, 4)); cat("\n")

# Save results 
save.image(paste("estrun_",run_id,"/counter_v2_seg",seg_id,".rdata",sep=""))

cat("This program is done.\n")

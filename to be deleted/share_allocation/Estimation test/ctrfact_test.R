library(ggplot2)
library(reshape2)
library(Rcpp)
library(RcppArmadillo)
library(maxLik)
library(evd)

model_name 	<- "MDCEV_share"
setwd("~/Documents/Research/Store switching/processed data")
plot.wd	<- '~/Desktop'
source("../Exercise/Multiple_discrete_continuous_model/0_Allocation_function.R")
source("../Exercise/main/share_allocation/ctrfact_sim_functions.r")

run_id		<- 1
# seg_id		<- 1
make_plot	<- FALSE

#############
# Functions #
#############
solveExp_fn	<- function(lambda, omega_fn, ln_inc, method = "FOC", lower = NULL, upper = NULL){
# Solve for the optimal expediture given parameter lambda and a given omega function, at a single value of log(income)	

# lambda	... A vector of utility paraemter of purchase utility
# omega_fn	... Spline function for a given income level
# ln_inc	...	A scalor of log(income)
#==============================================================================#
	if(is.null(lower))	{ lower <- 50 }
	if(is.null(upper))	{ upper <- 250}
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
			if(class(sol) == "try-error" | sol$minimum <= init.interval[1] + .001  | sol$minimum >= init.interval[2] - .001){
				interval	<- c(1, 500)
				sol	<- try(optimize(exp_fn, interval = interval, lambda = lambda, omega_fn = omega_fn, ln_inc = ln_inc), silent = TRUE)
				if(class(sol) == "try-error"){
					out <- 0
					cat("Expenditure solution was not found at ln_inc =", ln_inc, ".\n")
				}else{
					out <- sol		
				}
			}else{
				out	<- sol
			}
		}
	)
	
	return(out)
}

simExp_fn	<- function(lambda, omega_fn_ls, ln_inc, method = "FOC", use.bound = FALSE){
# A function to simulate optimal expenditure for a vector log(income)

# lambda		... A vector of utility parameters for purchase utility
# omega_fn_ls	...	A list of spline function, omega_fn, with each element corresponding to an income level 
# ln_inc		... A vector of log(income) as input data
#==============================================================================#
	out		<- vector("list", length = length(ln_inc))		#rep(NA, length(ln_inc))		
	lvinc	<- unique(ln_inc)				# Levels of discrete income bin
	lower 	<- NULL
	for(i in 1:length(lvinc)){
		if(use.bound & i>1){
			lower	<- max(sapply(out[1:(i-1)], function(x) x$minimum))
		}
		sel <- which(ln_inc == lvinc[i])
		sel1<- which(names(omega_fn_ls) == lvinc[i]) 
		out[[sel]]	<- solveExp_fn(lambda, omega_fn_ls[[sel1]], lvinc[i], method, lower = lower)
	}
	return(out)
}

##############
# Simulation #
##############
# Load estimation data 
load(paste("estrun_",run_id,"/MDCEV_est_seg",seg_id,".rdata",sep=""))
lambda	<- coef(sol.top2)
omega_deriv_save	<- omega_deriv
omega_deriv	<- omega.ls.07

# Check the concavity of omega function 
ggtmp	<- melt(apply(omega_draw, c(2, 3), mean, na.rm = T))
names(ggtmp)	<- c("lnInc", "y", "omega")
ggplot(ggtmp, aes(y, omega)) + geom_point() + geom_line() + facet_wrap(~lnInc, scales = "free_y") + 
	xlim(c(0, 500))

# Check the concavity of utility funciton 
# U = (lambda1 + lambda2 * lnInc)*omega(y, Inc) + (Inc - y)
# As long as omega(y, Inc) is concave in y and (lambda1 + lambda2*lnInc) >0, U will be concave
lambda[1] + lambda[2]*lnInc

# Plot the utility function to see if it is concave
tmp		<- seq(10, 200, 5)
ggtmp	<- data.frame()
for(i in 1:length(lnInc)){
	sel	<- which.min(abs(as.numeric(names(omega_deriv)) - lnInc[i]))
	for(j in 1:length(tmp)){
		u	<- -exp_fn(tmp[j], lambda, omega_fn = omega_deriv[[sel]], ln_inc = lnInc[i])
		ggtmp	<- rbind(ggtmp, c(lnInc[i], tmp[j], u))
	}
}
names(ggtmp)	<- c("lnInc", "y", "u")
ggplot(ggtmp, aes(y, u)) + geom_point() + geom_line() + 
		facet_wrap(~lnInc, scales = "free_y")
		
# Compute the 2nd derivative to see if omega''(y) <= 0 
tmp	<- seq(10, 200, 5)
ggtmp	<- data.frame()
for(i in 1:length(lnInc)){
	sel	<- which.min(abs(as.numeric(names(omega_deriv)) - lnInc[i]))
	ggtmp	<-	rbind(ggtmp, cbind(rep(lnInc[i], length(tmp)), tmp, omega_deriv[[sel]](tmp, deriv = 2)))
}
names(ggtmp)	<- c("lnInc", "y", "w_y2")
ggplot(ggtmp, aes(y, w_y2)) + geom_point() + geom_line() + 
		geom_hline(yintercept = 0, col = "red", size = .25) + 
		facet_wrap(~lnInc) + labs(title = "2nd Dervative of omega(y)")

# Compare the method of solving expenditure: FOC vs. optimize
a	<- simExp_fn(lambda, omega_fn_ls = omega_deriv, ln_inc = lnInc, method = "FOC")
b	<- simExp_fn(lambda, omega_fn_ls = omega_deriv, ln_inc = lnInc, method = "Utility")

par(mfrow = c(2, 1))
tmp1	<- sapply(a, function(x) x$root)
tmp2	<- sapply(b, function(x) x$minimum)
plot(tmp1, tmp2, xlim = range(c(tmp1, tmp2)), ylim = range(c(tmp1, tmp2)), xlab = "FOC", ylab = "Utility", 
		main = "Solution to optimal expenditure")
abline(a = 0, b = 1, col ="red")
tmp1	<- log(-sapply(a, function(x) x$obj))
tmp2	<- log(-sapply(b, function(x) x$objective))
plot(tmp1, tmp2, xlim = range(c(tmp1, tmp2)), ylim = range(c(tmp1, tmp2)), xlab = "FOC", ylab = "Utility", 
		main = "Overall utility")
abline(a = 0, b = 1, col ="red")

# Check the monotonicity of y along with lnInc
tmp	<- simExp_fn(lambda, omega_fn_ls = omega_deriv, ln_inc = lnInc, method = "Utility")
tmp	<- sapply(tmp, function(x) x$minimum) 
# sel	<- tmp < c(0, tmp[-length(tmp)]) 
sel	<- rep(TRUE, length(tmp))
tmp.opt	<- tmp[sel]
tmp.inc	<- lnInc[sel]
tmp	<- seq(-20, 20, 1)
ggtmp	<- data.frame()
for(i in 1:length(tmp.inc)){
	sel	<- which.min(abs(as.numeric(names(omega_deriv)) - tmp.inc[i]))
	# u	<- (lambda[1] + lambda[2]*tmp.inc[i]) * omega_deriv[[sel]](tmp.opt[i] + tmp)
	u	<- sapply(tmp.opt[i] + tmp, function(y) -exp_fn(y, lambda, omega_fn = omega_deriv[[sel]], ln_inc = tmp.inc[i]))
	ggtmp	<- rbind(ggtmp, cbind(tmp.inc[i], tmp.opt[i]+tmp, u, 1*(tmp==0)))
}
names(ggtmp)	<- c("lnInc", "y", "u", "sol")
ggplot(ggtmp, aes(y, u)) + geom_point(aes(col=factor(sol))) + geom_line() + 
		facet_wrap(~lnInc, scales = "free") + 
		scale_color_manual(values = c("black", "red"))



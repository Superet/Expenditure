library(ggplot2)
library(reshape2)
library(Rcpp)
library(RcppArmadillo)
library(maxLik)
library(evd)
library(chebpol)
library(foreach)
library(doParallel)
options(error = quote({dump.frames(to.file = TRUE)}))

model_name 	<- "MDCEV_share"
# setwd("~/Documents/Research/Store switching/processed data")
# plot.wd	<- '~/Desktop'
# source("../Exercise/Multiple_discrete_continuous_model/0_Allocation_function.R")
# source("../Exercise/main/share_allocation/ctrfact_sim_functions.r")

setwd("/home/brgordon/ccv103/Exercise/run")
source("0_Allocation_function.R")
source("ctrfact_sim_functions.r")

run_id		<- 1
seg_id		<- 1
make_plot	<- TRUE

ver.date	<- "2015-11-14"
load(paste("estrun_",run_id,"/MDCEV_est_seg",seg_id,"_", ver.date,".rdata",sep=""))

#############
# Functions #
#############
solveExp_fn	<- function(lambda, omega_fn, ln_inc, method = "FOC"){
# Solve for the optimal expediture given parameter lambda and a given omega function, at a single value of log(income)	

# lambda	... A vector of utility paraemter of purchase utility
# omega_fn	... Spline function for a given income level
# ln_inc	...	A scalor of log(income)
#==============================================================================#
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
			init.interval	<- c(10, 200)
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
			}else if (sol$minimum <= init.interval[1] + .1  | sol$minimum >= init.interval[2] - .1){
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

simExp_fn	<- function(lambda, omega_fn_ls, ln_inc, method = "FOC"){
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

###################################
# Simulate omega using quadrature # 
###################################
GH_num_nodes<- 30
y.interval	<- c(.1, 500)
y.nodes.c	<- chebknots(GH_num_nodes, interval = y.interval)[[1]]

# Register parallel computing
mycore 	<- 2
cl		<- makeCluster(mycore, type = "FORK")
registerDoParallel(cl)
cat("Register", mycore, "core parallel computing. \n")

omega_parallel	<- function(eps_draw){
	out		<- matrix(NA, length(lnInc), GH_num_nodes)
	for(j in 1:length(lnInc)){
		tmpX_list1	<- lapply(tmpX_list, function(x) rep(1, GH_num_nodes) %*% t(c(x, x*lnInc[j])))
		tmpsol 		<- incl_value_fn(param_est=tmp_coef, base= beta0_base, X_list=tmpX_list1, y=y.nodes.c, Q=Inf, price=tmp_price, 
							R=R, Ra=R, qz_cons = 0, exp_outside = FALSE, quant_outside = FALSE, eps_draw = rep(1, GH_num_nodes) %*% t(eps_draw) )
		out[j,]		<- tmpsol$omega
	}
	return(out)
}

pct			<- proc.time()
tmp			<- foreach(i = 1:numsim) %dopar% {
	omega_parallel(eps_draw[i,])
}
omega_draw_c	<- array(NA, c(numsim, length(lnInc), GH_num_nodes), dimnames = list(NULL, lnInc, y.nodes.c))
for(i in 1:numsim){
	omega_draw_c[i,,]	<- tmp[[i]]
}

use.time	<- proc.time() - pct
cat("Simulation of omega with random draws finishes with", use.time[3]/60, "min.\n")
cat("--------------------------------------------------------\n")
stopCluster(cl)
cat("Stop clustering. \n")

# Derive chebychev interpolation function
tmp		<- apply(omega_draw_c, c(2, 3), mean, na.rm =T, trim = .025)
omega_deriv_cheb	<- lapply(1:length(lnInc), function(i) chebappx(tmp[i,], interval = y.interval))
names(omega_deriv_cheb) <- lnInc

save.image("test_cheb.rdata")

######################################
# Compare omega from the two methods # 
######################################
lambda	<- coef(sol.top2)
plots	<- list(NULL)
idx		<- 1

# Plot omega over a grid values of y
tmp		<- seq(10, 200, 5)
ggtmp	<- data.frame()
for(i in 1:length(lnInc)){
	sel	<- which.min(abs(as.numeric(names(omega_deriv)) - lnInc[i]))
	for(j in 1:length(tmp)){
		u	<- -exp_fn(tmp[j], lambda, omega_fn = omega_deriv[[sel]], ln_inc = lnInc[i])
		ggtmp	<- rbind(ggtmp, data.frame(lnInc=lnInc[i], y = tmp[j], u = u, method = "spline"))
		u	<- -exp_fn(tmp[j], lambda, omega_fn = omega_deriv_cheb[[sel]], ln_inc = lnInc[i])
		ggtmp	<- rbind(ggtmp, data.frame(lnInc=lnInc[i], y = tmp[j], u = u, method = "cheb"))
	}
}
plots[[idx]]	<- ggplot(ggtmp, aes(y, u, col = method, alpha = .8)) + geom_point() + geom_line() + 
						facet_wrap(~lnInc, scales = "free_y") + 
						labs(title = "Utility function over y")
idx 	<- idx + 1

# Check if optimal solution for y lie in smooth areas. 
tmps	<- simExp_fn(lambda, omega_fn_ls = omega_deriv, ln_inc = lnInc, method = "Utility")
tmpc	<- simExp_fn(lambda, omega_fn_ls = omega_deriv_cheb, ln_inc = lnInc, method = "Utility")
tmp1	<- sapply(tmps, function(x) x$minimum) 
tmp2	<- sapply(tmpc, function(x) x$minimum)
tmp	<- seq(-20, 20, 1)
ggtmp	<- data.frame()
for(i in 1:length(lnInc)){
	sel	<- which.min(abs(as.numeric(names(omega_deriv)) - lnInc[i]))
	# u	<- (lambda[1] + lambda[2]*tmp.inc[i]) * omega_deriv[[sel]](tmp.opt[i] + tmp)
	u	<- sapply(tmp1[i] + tmp, function(y) -exp_fn(y, lambda, omega_fn = omega_deriv[[sel]], ln_inc = lnInc[i]))
	ggtmp	<- rbind(ggtmp, data.frame(lnInc = lnInc[i], y = tmp1[i]+tmp, u = u, sol = 1*(tmp==0), method="spline"))
	u	<- sapply(tmp2[i] + tmp, function(y) -exp_fn(y, lambda, omega_fn = omega_deriv_cheb[[sel]], ln_inc = lnInc[i]))
	ggtmp	<- rbind(ggtmp, data.frame(lnInc = lnInc[i], y = tmp2[i]+tmp, u = u, sol = 1*(tmp==0), method="cheb"))
}
names(ggtmp)	<- c("lnInc", "y", "u", "sol", "method")
plots[[idx]]	<- ggplot(ggtmp, aes(y, u, linetype = method, shape = method)) + geom_point(aes(col=factor(sol))) + geom_line() + 
				facet_wrap(~lnInc, scales = "free") + 
				scale_color_manual(values = c("black", "red"))

pdf(paste(plot.wd,"/graph_test_cheb.pdf", sep =""), width = 10, height = 10)
for(i in 1:idx){
	print(plots[[i]])
}
dev.off()


save.image("test_cheb.rdata")

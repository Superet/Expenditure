#############
# Functions #
#############
mysplfun	<- function(x, y){
	fit	<- smooth.spline(x, y)
	f	<- function(newx, deriv = 0){
		return(predict(fit, newx, deriv = deriv)$y)
	}
	return(f)
}

mytrimfun	<- function(x, alpha = .05){
	sort.x	<- sort(x)
	dropn	<- floor(alpha * length(x))
	drop.idx<- c(1:dropn, (length(x) - dropn + 1):length(x))
	return(sort.x[-drop.idx])
}

expFOC_fn <- function(y, lambda, omega_fn, ln_inc){
# First order condition of optimal expenditure for one single observation

# y			... Expenditure scalor
# lambda	...	A vector of utility parameters of purchase utility
# omega_fn	...	Spline function for a given income level within an income group
# ln_inc	... A scalor of log(income)
#==============================================================================#
	foc	<- (lambda[1] + lambda[2] * ln_inc) * omega_fn(y, deriv = 1) - 1
	return(foc)
}

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

solveExp_fn	<- function(lambda, omega_fn, ln_inc, method = "FOC", lower = NULL, upper = NULL){
# Solve for the optimal expediture given parameter lambda and a given omega function, at a single value of log(income)	

# lambda	... A vector of utility paraemter of purchase utility
# omega_fn	... Spline function for a given income level
# ln_inc	...	A scalor of log(income)
#==============================================================================#
	if(is.null(lower) & is.null(upper)){
		init.interval	<- c(50, 200)
	}else{
		if(is.null(lower)) 	{ lower <- 50 }
		if(is.null(upper))	{ upper	<- 250 }
		init.interval	<- c(lower, upper)
	}
	switch(method, 
		FOC = {
			sol	<- try(uniroot(expFOC_fn, init.interval, lambda = lambda, omega_fn = omega_fn, ln_inc = ln_inc), silent = TRUE)
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
					out	<- sol$root
				}
			}else{
				out	<- sol$root
			}
		}, 
		
		Utility = {			
			sol	<- try(optimize(exp_fn, interval = init.interval, lambda = lambda, omega_fn = omega_fn, ln_inc = ln_inc), silent = TRUE)
			if(class(sol) == "try-error" ){
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
	out		<- rep(NA, length(ln_inc))		
	lvinc	<- unique(ln_inc)				# Levels of discrete income bin
	for(i in 1:length(lvinc)){
		sel <- ln_inc == lvinc[i]
		sel1<- which(names(omega_fn_ls) == lvinc[i])
		out[sel]	<- solveExp_fn(lambda, omega_fn_ls[[sel1]], lvinc[i], method, lower = lower, upper = upper)
	}
	return(out)
}

SimWrapper_fn	<- function(omega_fn_ls, ln_inc, lambda, param_est, base, X_list, price, eps_draw, 
							sim.y = NULL, method = "FOC", use.bound = FALSE){
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
		lower	<- NULL
		if(sum(sel) > 0){
			sel1<- which(abs(as.numeric(names(omega_fn_ls)) - lvinc[i]) < 1e-6)
			if(is.null(sim.y)){
				if(use.bound & i>1) {lower	<- max(y[1:(i-1)]) }
				y[sel]	<- solveExp_fn(lambda, omega_fn_ls[[sel1]], lvinc[i], method = method, lower = lower)
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

SimOmega_fn	<- function(ln_inc, lambda, param_est, base, X_list, price, lnInc_lv, y.nodes, eps_draw, method, interp.method){
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
			tmpsol 		<- incl_value_fn(param_est=param_est, base= beta0_base, X_list=tmpX_list, y=y.nodes, Q=Inf, price=rep(1,numnodes) %*% t(price), 
						R=R, Ra=R, qz_cons = 0, exp_outside = FALSE, quant_outside = FALSE, eps_draw = rep(1, numnodes) %*% t(eps_draw[i,]) )
			omega_new[i,j,] <- tmpsol$omega
		}	
	}	

	# Interpolation functions or smoothing splines
	if(interp.method == "spline"){
		tmp1	<- apply(omega_new, c(2,3), mytrimfun)
		omega.ls.new	<- lapply(1:length(lnInc_lv), function(i) 
									mysplfun(x = rep(y.nodes, each = dim(tmp1)[1]), y = as.vector(tmp1[,i,])))							
	}else{
		tmp1	<- apply(omega_new, c(2, 3), mean, trim = .05)
		omega.ls.new <- lapply(1:length(lnInc_lv), function(i) chebfun(x = y.nodes, y = tmp1[i,], interval = y.interval))
	}	
	# omega.ls.new <- lapply(1:length(lnInc_lv), function(i) splinefun(y.nodes, colMeans(omega_new[,i,], na.rm = T), method = "natural"))
	names(omega.ls.new)	<- lnInc_lv

	# Simulate expenditure and share
	out	<- SimWrapper_fn(omega_fn_ls = omega.ls.new, ln_inc = ln_inc, lambda = lambda, param_est = param_est, 
						base = base, X_list = X_list, price = price, eps_draw = eps_draw, method = method)
	return(out)
}

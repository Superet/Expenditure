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

spl2dfun	<- function(data, fml = omega ~ te(lnInc, y), return.fit = FALSE){
	fit	<- gam(fml, data = data)
	f	<- function(newx, deriv = 0){
		if(deriv == 0){
			return(predict(fit, newdata = newx, type = "link"))
		}else if(deriv == 1){
			# We only compute the marginal derivative along the y dimsion
			eps	<- 1e-7
			newx1 <- newx
			newx1[,"y"] <- newx[,"y"] + eps
			X0	<- predict(fit, newdata = newx, type = "lpmatrix")
			X1	<- predict(fit, newdata = newx1, type = "lpmatrix")
			d	<- (X1 - X0) %*% coef(fit) /eps
			return(d)
		}
	}
	if(return.fit){
		return(fit)
	}else{
		return(f)
	}
}

mytrimfun<- function (x, alpha = 0.05){
	if(alpha == 0){
		return(x)
	}else{
		n	<- length(x)
		lo 	<- floor(n * alpha) + 1
		hi 	<- n + 1 - lo
		out <- sort.int(x, partial = unique(c(lo, hi)), na.last = FALSE)[lo:hi]
		return(out)
	}
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
		init.interval	<- c(50, 250)
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

simExp_fn	<- function(lambda, omega_deriv, ln_inc, method = "FOC", lower = NULL, upper = NULL){
# A function to simulate optimal expenditure for a vector log(income)

# lambda		... A vector of utility parameters for purchase utility
# omega_deriv	...	A function of 2d spline 
# ln_inc		... A vector of log(income) as input data
#==============================================================================#
	out		<- rep(NA, length(ln_inc))		
	lvinc	<- unique(ln_inc)				# Levels of discrete income bin
	for(i in 1:length(lvinc)){
		sel <- ln_inc == lvinc[i]
		f 	<- function(y, ...){
			newx	<- data.frame(lnInc = lvinc[i], y = y)
			omega_deriv(newx, ...)
		}
		out[sel]	<- solveExp_fn(lambda, f, lvinc[i], method, lower = lower, upper = upper)
	}
	return(out)
}

SimWrapper_fn	<- function(omega_deriv, ln_inc, lambda, param_est, base, X_list, price, eps_draw, 
							sim.y = NULL, method = "FOC", use.bound = FALSE, ret.sim = FALSE, par.draw = NULL){
# This function simulates expenditure and share, having known omega_fn_ls. 

# omega_deriv	...	A function of 2d spline 
# ln_inc		...	A vector of log(income) as simulation input data
# lambda		... A vector of utility parameters for purchase utility
# param_est		... A vector of parameters in the conditional allocation model
# base			... Integer index of which retailer has intercept of 0 in the conditional allocation model
# X_list		...	A list of retail attributes
# price			...	A vector of price 
# esp_draw		... A matrix of random draw
# sim.y			... If not NULL, then expenditure is given and only expenditure share is simulated. 
# method		... The method of solving for optimal expenditure. Take value of "FOC" (first-order condition equation)
#					and "Utility" (directly maximizing utility function)
# use.bound		...	Logical variable indicating if use bounds in solving for optimal expenditure
# ret.sim		...	Logical variable indicating if return all the simulation values for each random raw, or 
#					return average. 
# base.ls		... A list of baseline simulations. Each elment is a array of expenditre at each retail format.
# 					If it is not null, then we compute the standard error of the difference between the current
#					simulation and baseline simulation. 
#=========================================================================================================#
	
	numsim	<- nrow(eps_draw)	
	if(is.null(sim.y)){
		y		<- matrix(NA, numsim, length(ln_inc))
	}else{
		y 		<- sim.y
	}
	if(is.null(par.draw)){
		shr.draw	<- matrix(0, numsim, length(param_est))
	}else{
		lambda.draw	<- par.draw[,c("lambda1", "lambda2")]
		shr.draw	<- par.draw[,setdiff(colnames(par.draw), c("lambda1", "lambda2"))]
	}
	e		<- array(NA, c(numsim, length(ln_inc), R), dimnames = list(iter = 1:numsim, lnInc = ln_inc, retailer = 1:R))
	omega	<- matrix(NA, numsim, length(ln_inc))
	omega1	<- matrix(NA, numsim, length(ln_inc))
	lvinc	<- sort(unique(ln_inc))
	
	for(i in 1:length(lvinc)){
		sel <- ln_inc == lvinc[i]
		lower	<- NULL
		if(sum(sel) > 0){			
			# Omega function conditional on income level 
			f 	<- function(y, ...){
				newx	<- data.frame(lnInc = lvinc[i], y = y)
				omega_deriv(newx, ...)
			}
			tmpX_list1	<- lapply(X_list, function(x) rep(1, sum(sel)) %*% t(c(x, x*lvinc[i])))
			
			# Note that expenditure y is independent of eps_draw
			if(is.null(sim.y)){
				if(is.null(par.draw)){
					y[,sel]	<- solveExp_fn(lambda, f, lvinc[i], method = method, lower = lower)		# This is a scalor
				}else{
					y[,sel]	<- apply(lambda.draw, 1, function(x) solveExp_fn(lambda+x, f, lvinc[i], method = method, lower = lower) )
				}
			}

			for(j in 1:numsim){				
				tmpsol 		<- try(incl_value_fn(param_est+shr.draw[j,], base, X_list=tmpX_list1, y=y[j,sel], Q=Inf, price=rep(1, sum(sel)) %*% t(price), 
							R=R, Ra=R, qz_cons = 0, exp_outside = FALSE, quant_outside = FALSE, 
							eps_draw = rep(1, sum(sel)) %*% t(eps_draw[j,]) ), silent = TRUE)
				if(class(tmpsol) != "try-error"){
					omega[j,sel] 	<- tmpsol$omega			
					omega1[j,sel]	<- omega_deriv(newx = data.frame(lnInc = lvinc[i], y = y[j,sel]))
					e[j,sel,]		<- tmpsol$e
				}
			}			
		}
	}

	out	<- cbind(colMeans(y, na.rm=T), colMeans(omega, na.rm = T), colMeans(omega1, na.rm = T), apply(e, c(2,3), mean, na.rm= T))	
	colnames(out)	<- c("Exp", "omega", "spl.omega", fmt_name)
	if(ret.sim){
		out	<- list(Average = out, y = y, Allocation = e)
	}
	return(out)
}

SimOmega_fn	<- function(ln_inc, lambda, param_est, base, X_list, price, lnInc_lv, y.nodes, eps_draw, 
						method, interp.method, ret.sim, alpha = 0.05, par.draw = NULL){
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
			tmpsol 		<- try(incl_value_fn(param_est=param_est, base= beta0_base, X_list=tmpX_list, y=y.nodes, Q=Inf, price=rep(1,numnodes) %*% t(price), 
						R=R, Ra=R, qz_cons = 0, exp_outside = FALSE, quant_outside = FALSE, eps_draw = rep(1, numnodes) %*% t(eps_draw[i,]) ), silent = TRUE)
			if(class(tmpsol) == "try-error"){
				cat("Error at computing omega at eps_draw[",i,"] and lnInc[",j,"].\n")
				next
			}else{
				omega_new[i,j,] <- tmpsol$omega
			}
		}	
	}	

	# Interpolation functions or smoothing splines
	if(interp.method == "spline"){
		tmp1		<- apply(omega_new, c(2,3), mytrimfun, alpha = alpha)
		tmpdat	<- melt(tmp1)
		names(tmpdat)	<- c("iter", "lnInc", "y", "omega")
		omega_new	<- spl2dfun(tmpdat)
		# if(!is.list(tmp1)){
		# }else{
		# 	# In case of NA in omega_new, tmp1 turns out to be list, then I take it carefully
		# 	omega.new	<- vector("list", length(lnInc_lv))
		# 	for(i in 1:length(lnInc_lv)){
		# 		tmp1	<- lapply(1:numnodes, function(j) mytrimfun(omega_new[,i,j], alpha = alpha))		# tmp1 is a list of vector of different dim.
		# 		tmpx	<- unlist(lapply(1:numnodes, function(j) rep(y.nodes[j], length(tmp1[[j]]))))
		# 		omega.new[[i]]	<- mysplfun(x = tmpx, y = unlist(tmp1))
		# 	}
		# }			
	}else{
		tmp1	<- apply(omega_new, c(2, 3), mean, trim = .05)
		omega.new <- lapply(1:length(lnInc_lv), function(i) chebfun(x = y.nodes, y = tmp1[i,], interval = y.interval))
		names(omega.new)	<- lnInc_lv
	}	
	
	# Simulate expenditure and share
	out	<- SimWrapper_fn(omega_deriv = omega.new, ln_inc = ln_inc, lambda = lambda, param_est = param_est, 
						base = base, X_list = X_list, price = price, eps_draw = eps_draw, method = method, 
						ret.sim = ret.sim, par.draw = par.draw)
	return(out)
}

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

omega.lm <- function(data, fml = omega ~ y, return.fit = FALSE){
	fit	<- lm(fml, data = data)
	f	<- function(newx, deriv = 0){
		if(deriv == 0){
			return(predict(fit, newdata = newx))
		}else if(deriv == 1){
			# We only compute the marginal derivative along the y dimsion
			eps	<- 1e-7
			newx1 <- newx
			newx1[,"y"] <- newx[,"y"] + eps
			X0	<- predict(fit, newdata = newx)
			X1	<- predict(fit, newdata = newx1)
			d	<- (X1 - X0)/eps
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

expFOC_fn <- function(y, lambda, omega_fn, Z, inc){
# First order condition of optimal expenditure for one single observation

# y			... Expenditure scalor
# lambda	...	A vector of utility parameters of purchase utility
# omega_fn	...	Spline function for a given income level within an income group
# ln_inc	... A scalor of log(income)
#==============================================================================#
	foc	<- (Z %*% lambda) * omega_fn(y, deriv = 1) - 1
	return(foc)
}

exp_fn <- function(y, lambda, omega_fn, Z, inc){
# Utility function expenditure for one single observation

# y			... Expenditure scalor
# lambda	...	A vector of utility parameters of purchase utility
# omega_fn	...	Spline function for a given income level within an income group
# ln_inc	... A scalor of log(income)
#==============================================================================#
	out	<- as.vector(Z %*% lambda) * omega_fn(y) + (inc - y)
	return(-out)
}

solveExp_fn	<- function(lambda, omega_fn, inc, Z, method = "FOC", lower = NULL, upper = NULL, ...){
# Solve for the optimal expediture given parameter lambda and a given omega function, at a single value of log(income)	

# lambda	... A vector of utility paraemter of purchase utility
# omega_fn	... Spline function for a given income level
# ln_inc	...	A scalor of log(income)
#==============================================================================#
	if(is.null(lower) & is.null(upper)){
		init.interval	<- c(50, 1000)
	}else{
		if(is.null(lower)) 	{ lower <- 100 }
		if(is.null(upper))	{ upper	<- 1000 }
		init.interval	<- c(lower, upper)
	}
	switch(method, 
		FOC = {
			sol	<- try(uniroot(expFOC_fn, init.interval, lambda = lambda, omega_fn = omega_fn, Z = Z, inc = inc), silent = TRUE)
			# The root solution is sensitive to the solution range. 
			# So we search for wider range gradually
			if(class(sol) == "try-error"){
				maxpt	<- seq(1000, 1500, by = 100)
				i 		<- 1
				while(i <= length(maxpt) & class(sol) == "try-error"){
					sol	<- try(uniroot(expFOC_fn, c(1, maxpt[i]), lambda = lambda, omega_fn = omega_fn, Z = Z, inc = inc), silent = TRUE)
					i	<- i + 1
				}
				if(class(sol) == "try-error"){
					out <- 0
					cat("Expenditure solution was not found at inc =", inc, ".\n")
				}else{
					out	<- sol$root
				}
			}else{
				out	<- sol$root
			}
		}, 
		
		Utility = {			
			sol	<- try(optimize(exp_fn, interval = init.interval, lambda = lambda, omega_fn = omega_fn,  Z = Z, inc = inc), silent = TRUE)
			if(class(sol) == "try-error" ){
				interval	<- c(1, 500)
				sol	<- try(optimize(exp_fn, interval = interval, lambda = lambda, omega_fn = omega_fn,  Z = Z, inc = inc), silent = TRUE)
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

simExp_fn	<- function(lambda, omega_deriv, inc, Z, sim.data, method = "FOC", lower = NULL, upper = NULL){
# A function to simulate optimal expenditure for a vector log(income)

# lambda		... A vector of utility parameters for purchase utility
# omega_deriv	...	A function of 2d spline 
# ln_inc		... A vector of log(income) as input data
#==============================================================================#
	out		<- rep(NA, length(inc))		
	lvinc	<- unique(inc)				# Levels of discrete income bin
	for(i in 1:length(lvinc)){
		sel <- inc == inc[i]
		f 	<- function(y, ...){
			newx	<- sim.data
			newx$y	<- y
			omega_deriv(newx, ...)
		}
		out[sel]	<- solveExp_fn(lambda, omega_fn = f, inc = lvinc[i], Z = Z[sel,], method, lower = lower, upper = upper)
	}
	return(out)
}

SimWrapper_fn	<- function(omega_deriv, lambda, shr.par, base, X1, X2, price, demo_mat, eps_draw, ret.idx, demo.x = FALSE,
							sim.y = NULL, method = "FOC", use.bound = FALSE, ret.sim = FALSE, par.draw = NULL, exp.lower = NULL, exp.upper = NULL){
# This function simulates expenditure and share, having known omega_fn_ls. 

# omega_deriv	...	A function of 2d spline 
# ln_inc		...	A vector of log(income) as simulation input data
# lambda		... A vector of utility parameters for purchase utility
# param_est		... A vector of parameters in the conditional allocation model
# base			... Integer index of which retailer has intercept of 0 in the conditional allocation model
# X_list		...	A list of retail attributes, each element is a matrix (length(ln_inc) \times nx)
# price			...	A matrix of price (length(ln_inc) \times R )
# esp_draw		... A matrix of random draw
# sim.y			... If given as a matrix (numsim \times length(ln_inc)), then expenditure is given and only expenditure share is simulated. 
# method		... The method of solving for optimal expenditure. Take value of "FOC" (first-order condition equation)
#					and "Utility" (directly maximizing utility function)
# use.bound		...	Logical variable indicating if use bounds in solving for optimal expenditure
# ret.sim		...	Logical variable indicating if return all the simulation values for each random raw, or 
#					return average. 
#=========================================================================================================#
	
	numsim	<- nrow(eps_draw)	
	N 		<- nrow(demo_mat)
	R		<- length(X1)
	nx		<- ncol(as.matrix(X1[[1]]))
	inc		<- exp(demo_mat[,"ln_inc"])
	if(is.null(sim.y)){
		y		<- matrix(NA, numsim, N)
	}else{
		y 		<- sim.y
	}
	if(is.null(par.draw)){
		shr.draw	<- matrix(0, numsim, length(shr.par))
	}else{
		lambda.draw	<- par.draw[,c("lambda1", "lambda2")]
		shr.draw	<- par.draw[,setdiff(colnames(par.draw), c("lambda1", "lambda2"))]
	}
	e		<- array(NA, c(numsim, N, R), dimnames = list(iter = 1:numsim, obs = 1:N, retailer = 1:R))
	omega	<- matrix(NA, numsim, N)
	
	# Solve for expenditure 
	if(is.null(sim.y)){
		newdat	<- data.frame(demo_mat)
		newdat$price	<- price
		if(demo.x){
			newdat$X<- as.matrix(cbind(demo_mat[,setdiff(colnames(demo_mat), c("1", "Intercept", "ln_inc"))], do.call(cbind, lapply(X1, function(x) x[,ret.idx]))))
		}else{
			newdat$X<- as.matrix(do.call(cbind, lapply(X1, function(x) x[,ret.idx])))
		}
		if(ncol(demo_mat) == 3){
			Z1	<- demo_mat[,1:2]
		}else{
			Z1	<- as.matrix(demo_mat)
		}

		for(i in 1:N){
			f 	<- function(y, ...){
				newdata	<- cbind(y = y, newdat[i,])
				omega_deriv(newdata, ...)
			}
			
			y[,i]	<- solveExp_fn(lambda, omega_fn = f, inc = inc[i], Z = Z1[i,], method = method, lower = exp.lower, upper = exp.upper)		# This is a scalor
		}
	}

	# Solve for allocation 
	for(j in 1:numsim){				
		sel			<- !is.na(y[j,]) & y[j,] > 0 
		if(sum(sel)<N){
			tmpX1		<- lapply(X1, function(x) x[sel,])
			tmpX2		<- lapply(X2, function(x) x[sel,])
			tmpsol 		<- try(incl_value_fn(shr.par, nx, base, X1 = tmpX1, X2=tmpX2, y = y[j,sel], 
							price[sel,], R, eps_draw = eps_draw[j,], returnmax = TRUE), 
						   silent = TRUE)
		}else{
			tmpsol 		<- try(incl_value_fn(shr.par, nx, base, X1 = X1, X2=X2, y = y[j,], 
							price, R, eps_draw = rep(1, N) %*% t(eps_draw[j,]), returnmax = TRUE), 
						   silent = TRUE)
		}
		if(class(tmpsol) != "try-error"){
			omega[j,sel] 	<- tmpsol$omega	
			e[j,sel,]		<- tmpsol$e
		}
	}
	
	# out	<- cbind(colMeans(y, na.rm=T), colMeans(omega, na.rm = T), colMeans(omega1, na.rm = T), apply(e, c(2,3), mean, na.rm= T))	
	out	<- cbind(colMeans(y, na.rm=T), colMeans(omega, na.rm = T), apply(e, c(2,3), mean, na.rm= T))	
	
	if(length(fmt_name) == R){
		# colnames(out)	<- c("Exp", "omega", "spl.omega", fmt_name)
		colnames(out)	<- c("Exp", "omega", fmt_name)
	}else{
		# colnames(out)	<- c("Exp", "omega", "spl.omega", names(X_list))
		colnames(out)	<- c("Exp", "omega", names(X1))
	}
	
	if(ret.sim){
		out	<- list(Average = out, y = y, Allocation = e)
	}
	return(out)
}

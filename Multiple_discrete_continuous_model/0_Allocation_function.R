# Purchae utility function 
uP_fn <- function(e, psi, gamma, price, R, Ra, qz_cons = 0, exp_outside = TRUE, quant_outside = TRUE){
# Negative purchase utility 
# e: Ra-dimentional vector
# psi, gamma is Ra-dimension vector, price is R-dimensional vector.
# exp_outside and quant_outside are indicator whether there is expenditure/quantity outside option. 
	e0 		<- e[1:R]
	e_R1	<- e[(R+1):Ra]
	u0_part1<- log(e0/(gamma[1:R]*price) + 1)
	if(exp_outside){
		if(quant_outside){
			u0_part2	<- c(log(e_R1[1]/gamma[(R+1)]), log(e_R1[2]/gamma[Ra]+qz_cons) )
		}else{
			u0_part2	<- log(e_R1[1]/gamma[(R+1)])
		}
	}else{
		if(quant_outside){
			u0_part2	<- log(e_R1/gamma[Ra]+qz_cons)
		}else{
			u0_part2	<- NULL
		}
	}
	
	u0		<- c(u0_part1, u0_part2 )
	out 	<- sum(psi*gamma*u0 ) 
	return(out)
}

uPGrad_fn <- function(e, psi, gamma, price, R, Ra, qz_cons = 0, exp_outside = TRUE, quant_outside = TRUE){
	e0 		<- e[1:R]
	e_R1	<- e[(R+1):Ra]
	u0_part1<- 1/(e0+price*gamma[1:R])
	if(exp_outside){
		if(quant_outside){
			u0_part2	<- c(1/e_R1[1], 1/(e_R1[2]+qz_cons*gamma[Ra])  )
		}else{
			u0_part2	<- 1/e_R1[1]
		}
	}else{
		if(quant_outside){
			u0_part2	<- 1/(e_R1 + qz_cons*gamma[Ra])
		}else{
			u0_part2	<- NULL
		}
	}
	
	u0		<- c(u0_part1, u0_part2 )
	out 	<- psi*gamma*u0 
	return(out)	
}

Allocation_constr_fn <- function(y, psi, gamma, Q, price, R, Ra, qz_cons = 0, exp_outside = TRUE, quant_outside = TRUE, 
						  inits=NULL, silent=FALSE){
# Solve constrained optimization with "constrOptim". 						  
	
	# Construct the constrained equations: ui %*% theta - ci >= 0
	ui	<- diag(Ra)				# Ra non-negativity constraints. 
	ci 	<- rep(0, Ra)
	if(exp_outside){
		ui	<- rbind(ui, c(rep(-1, R+1),rep(0, Ra-R-1 )) )
		ci	<- c(ci, -y)
		if(quant_outside){
			ci[(length(ci)-1)] <- (-qz_cons)
			ui <- rbind(ui, c(-1/price, 0, -1))
			ci <- c(ci, -Q)
		}
	}else{
		if(quant_outside){
			ci[length(ci)] <- (-qz_cons)
			ui	<- rbind(ui, c(rep(-1,R),0), c(-1/price, -1) )
			ci	<- c(ci, -y, -Q)
		}else{
			ui 	<- rbind(ui, rep(-1,R))
			ci	<- c(ci, -y)
		}
	}
	
	# Set initial values for the optimization
	if(is.null(inits)){
		eps <- 1e-6
		if(exp_outside){
			if(quant_outside){
				tmp	  <- min(y/price, Q) * .99
				inits <- list(c(tmp/(Ra-1)*price, tmp/(Ra-1), Q - tmp - eps), 
						  c(rep(eps, R), y-(R+2)*eps, Q - sum(price)*eps - eps),
						  c(tmp/R * price, y- sum(tmp/R*price) -eps, eps) )	
			}else{
				tmp	  <- min(y/price) * .99
				inits <- list(c(tmp/Ra*price, tmp/Ra), 
						  	c(tmp/R*price, eps) )	
			}
		}else{
			if(quant_outside){
				tmp		<- min(y/price, Q) * .99
				inits	<- list( c(tmp/R * price,  Q-tmp - eps), 
								 c(rep(eps,R), Q - sum(eps/price)- eps) )
			}else{
				tmp	  <- min(y/price) * .99
				inits	<- list( psi/sum(psi) * y * .99, tmp/R*price)
			}
		}
	}
	
	# Select the optimal value.
	if(Q <= 0){
		e <- c(rep(0,R),y)
		if(Ra>R+1){
			e <- c(e, 0)
		}
		return(list(e = e, max = log(y) ) )
	}else{
		sol.list <- vector("list", length(inits))
		for(j in 1:length(inits)){
			sol.list[[j]] <- try(constrOptim(theta = inits[[j]], f=uP_fn, grad=uPGrad_fn, ui=ui, ci=ci, control = list(fnscale = -1), psi=psi,
						 gamma=gamma, price=price, R=R, Ra=Ra, qz_cons = qz_cons, exp_outside = exp_outside, quant_outside = quant_outside), silent=silent)
		}
		
		if(length(inits) == 1){
			sol	<- sol.list[[1]]
			if(class(sol) == "try-error"){
				sol <- list(par = rep(NA, Ra), value=NA, convergence = NA)
			}
		}else{
			sol.list1 <- sol.list[sapply(sol.list, function(x) !inherits(x, "try-error"))]
			sel <- which.max(sapply(sol.list1, function(x) x$value))
			if(length(sel)==0){
				sol <- list(par = rep(NA, Ra), value=NA, convergence = NA)
			}else{
				sol <- sol.list1[[sel]]
				if(sol$convergence != 0){
					cat("Constrained optimization does not converge at value y=",y,", Q=",Q,"\n")
				}
			}
		}		
		return(list(e = sol$par, max = sol$value, convergence = sol$convergence))
	}			
}

Allocation_nlop_fn <- function(y, psi, gamma, Q, price, R, Ra, qz_cons = 0, exp_outside = TRUE, quant_outside = TRUE, 
						  inits=NULL, silent=FALSE, nloptr.opts = list("algorithm"="NLOPT_LD_MMA", "maxeval" = 200, xtol_rel = 1e-8)){
# Solve the constrained optimization with augmented Lagrangian algorithm

	# Define the objective function and inequality function
	obj_fn <- function(e) {-uP_fn(e, psi=psi, gamma=gamma, price=price, R=R, Ra=Ra, qz_cons = qz_cons, exp_outside = exp_outside, quant_outside = quant_outside)}
	grad_fn<- function(e) {-uPGrad_fn(e, psi=psi, gamma=gamma, price=price, R=R, Ra=Ra, qz_cons = qz_cons, exp_outside = exp_outside, quant_outside = quant_outside)}
	if(exp_outside){
		if(quant_outside){
			gin	<- function(e) {c(e, y - sum(e[1:(R+1)]), Q - sum(e[1:R]/price) - e[Ra] )}
			gjac<- function(e) {rbind(diag(Ra), c(rep(-1, R+1),0), c(-1/price, 0, -1)) }
		}
	}else{
		if(!quant_outside){
			gin <- function(e) { return(sum(e) - y) }
			gjac<- function(e){
				out	<- matrix(NA, 1, length(e))
				out[1,]	<- rep(1, length(e))
				return(out)
			}
		}
	}
	
	# Set initial values for the optimization
	if(is.null(inits)){
		eps <- 1e-6
		if(exp_outside){
			if(quant_outside){
				tmp	  <- min(y/price, Q) * .99
				inits <- list(c(tmp/(Ra-1)*price, tmp/(Ra-1), Q - tmp - eps), 
						  c(rep(eps, R), y-(R+2)*eps, Q - sum(price)*eps - eps),
						  c(tmp/R * price, y- sum(tmp/R*price) -eps, eps) )	
			}else{
				tmp	  <- min(y/price) * .99
				inits <- list(c(tmp/Ra*price, tmp/Ra), 
						  	c(tmp/R*price, eps) )	
			}
		}else{
			if(quant_outside){
				tmp		<- min(y/price, Q) * .99
				inits	<- list( c(tmp/R * price,  Q-tmp - eps), 
								 c(rep(eps,R), Q - sum(eps/price)- eps) )
			}else{
				inits	<- list( psi/sum(psi) * y )
			}
		}
	}
	
	
	# Solve for the optimal values. 
	if(Q <= 0){
		e <- c(rep(0,R),y)
		if(Ra>R+1){
			e <- c(e, 0)
		}
		return(list(e = e, max = log(y) ) )
	}else{
		sol.list <- vector("list", length(inits))
		for(j in 1:length(inits)){
			sol.list[[j]] <- try(nloptr(x0=inits[[j]], eval_f = obj_fn, eval_grad_f=grad_fn, 
							lb = rep(0, length(inits[[j]])), ub = rep(y, length(inits[[j]])), 
							eval_g_ineq = gin, eval_jac_g_ineq=gjac, opts = nloptr.opts), 
							silent=silent)
		}
		
		if(length(inits) == 1){
			sol	<- sol.list[[1]]
			if(class(sol) == "try-error"){
				sol <- list(solution = rep(NA, Ra), objective=NA, status = NA)
			}
		}else{
			sol.list1 <- sol.list[sapply(sol.list, function(x) !inherits(x, "try-error"))]
			sel <- which.min(sapply(sol.list1, function(x) x$objective))
			if(length(sel)==0){
				sol <- list(solution = rep(NA, Ra), objective=NA, status = NA)
			}else{
				sol <- sol.list1[[sel]]
				# if(sol$status != 0){
				# 	cat("Constrained optimization does not converge at value y=",y,", Q=",Q,"\n")
				# }
			}
		}
		return(list(e = sol$solution, max = -sol$objective, convergence = sol$status))
	}
}

Allocation_fn <- function(y, psi, gamma, Q, price, R, Ra, qz_cons = 0, exp_outside = TRUE, quant_outside = TRUE, optim_method = c("nloptr"), inits=NULL,...){
	# Check conditions;
	if(exp_outside & quant_outside){
		if(Ra != R+2) stop("Ra should be equal to R+2.")
	}else if(!exp_outside & !quant_outside){
		if(Ra != R) stop("Ra should be equal to R.")
	}else{
		if(Ra != R+1) stop("Ra should be equal to R+1.")
	}
	
	switch(optim_method, 
		constrOptim = { out <- Allocation_constr_fn(y = y, psi=psi, gamma=gamma, Q = Q, price=price, 
				R=R, Ra=Ra, qz_cons = qz_cons, exp_outside = exp_outside, quant_outside = quant_outside, inits=inits, ...) }, 
		nloptr = { out <- Allocation_nlop_fn(y = y, psi=psi, gamma=gamma, Q = Q, price=price, 
				R=R, Ra=Ra, qz_cons = qz_cons, exp_outside = exp_outside, quant_outside = quant_outside, inits=inits, ...)}
	)
	return(out) 
}

param_assignR	<- function(param, nx, R, base = 0){
	beta	<- param[1:nx]
	if(base ==0 ){
		beta0	<- c(0, param[(nx+1):(nx+R-1)])
	}else{
		beta0	<- rep(0, R)
		beta0[-base]	<- param[(nx+1):(nx+R-1)]
	}
	gamma	<- param[(nx+R):(nx+R*2-1)]
	# sigma	<- 1
	sigma	<- exp(param[(nx+R*2)])
	return(list(beta = beta, beta0 = beta0, gamma = gamma, sigma = sigma))
}

incl_value_fn <- function(param_est, nx, base, X_list, y, Q, price, R, Ra, qz_cons = 0, exp_outside = TRUE, quant_outside = TRUE, eps_draw = NULL, inits=NULL,...){
# param_est		...	parameter estimates
# X_array		...	Attributes array, dimension(R*Nobs*nx)
# y				... A Nobs-dimensional expenditure budget. 
# Q				... Q is the scalor of the quantity constraint.
	param_out	<- param_assignR(param_est, nx, R, base)
	beta		<- param_out$beta
	beta0		<- param_out$beta0
	gamma		<- param_out$gamma
	sigma		<- param_out$sigma
	Nobs		<- length(y)
	if(exp_outside){
		psi_R1 	<- param_out$psi_R1
		gamma	<- c(gamma, 1)
	}
	if(quant_outside){
		psi_R2	<- param_out$psi_R2
		gamma	<- c(gamma, 1)
	}
	
	# Compute the parameters for allocation function 
	if(is.null(eps_draw)){
		eps_draw 	<- rep(1, Nobs) %*% t(rgev(R)) * sigma				# matrix(rgev(Nobs*R), Nobs, R)
	}
	xbeta	<- sapply(1:R, function(i) X_list[[i]] %*% beta + beta0[i])
	psi		<- exp(xbeta + eps_draw)
	if(exp_outside)		{psi <- cbind(psi, psi_R1) }
	if(quant_outside) 	{psi <- cbind(psi, psi_R2) }
	
	# Call allocation function 
	omega 	<- rep(NA, Nobs) 
	e		<- matrix(NA, Nobs, Ra)
	for(i in 1:Nobs){
		if(y[i] < 1e-5){
			omega[i]	<- 0 
			e[i,]		<- 0
		}else{
			sol	<- Allocation_fn(y = y[i], psi=psi[i,], gamma=gamma, Q = Q, price=price[i,], 
					R=R, Ra=Ra, qz_cons = qz_cons, exp_outside = exp_outside, quant_outside = quant_outside, inits=inits, ...)
			omega[i]	<- sol$max
			e[i,]		<- sol$e
		}
	}
	
	# Return results 
	return(list(e = e, omega = omega))
}	

cheb.basis	<- function(x, M, interval = c(0, 1), adj = TRUE){
# This function returns the Chebyshev basis function. 
# x			...	a vector of interpolating knots
# M			...	a scalor of Chebyshev degree
# interval	...	interval of knots
# adj		... logical value of whether to adjust x to [0,1]
	n		<- length(x)
	if(adj){
		x.adj	<- 2*(x - interval[1])/(interval[2] - interval[1]) - 1
	}else{
		x.adj	<- x
	}
	
	out	<- matrix(NA, n, M)
	out[,1]	<- 1
	out[,2]	<- x.adj
	for(i in 3:M){
		out[,i]	<- 2 * x.adj * out[,(i-1)] - out[,(i-2)]
	}
	return(out)
}

cheb.1d.basis <- function(x, M, interval = c(0, 1), adj = T){
# This function returns the 1st derivative of Chebyshev basis function. 
# x			...	a vector of interpolating knots
# M			...	a scalor of Chebyshev degree
# interval	...	interval of knots
# adj		... logical value of whether to adjust x to [0,1]
	n		<- length(x)
	if(adj){
		x.adj	<- 2*(x - interval[1])/(interval[2] - interval[1]) - 1
	}else{
		x.adj	<- x
	}
	T		<- matrix(NA, n, M)
	out		<- matrix(NA, n, M)
	T[,1]	<- 1
	T[,2]	<- x.adj
	out[,1]	<- 0
	out[,2]	<- 1
	for(i in 3:M){
		T[,i]	<- 2 * x.adj * T[,(i-1)] - T[,(i-2)]
		out[,i]	<- 2*x.adj*out[,(i-1)] + 2*T[,(i-1)] - out[,(i-2)]
	}
	return(out)	
}

chebfun	<- function(x, y, interval = c(0, 1), ... ){
	df	<- length(x)
	T	<- cheb.basis(x = x, M = df, interval = interval, adj = TRUE)
	alpha	<- solve(t(T) %*% T) %*% t(T) %*% y
	f	<- function(x, deriv = 0, dT = NULL){
		if(deriv == 0){
			T1	<- cheb.basis(x = x, M = df, interval = interval, adj = TRUE)
			out	<-as.vector(T1 %*% alpha)
		}else if(deriv == 1){
			if(is.null(dT)){
				dT	<- cheb.1d.basis(x = x, M = df, interval = interval)
			}
			out	<- dT %*% (alpha * 2 /(interval[2] - interval[1]))
		}	
		return( out )
	}
	return(f)
}

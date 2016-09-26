# Purchae utility function 
uP_fn <- function(e, psi, gamma, price, R){
# Negative purchase utility 
# e: Ra-dimentional vector
# psi, gamma is Ra-dimension vector, price is R-dimensional vector.
# exp_outside and quant_outside are indicator whether there is expenditure/quantity outside option. 
	u0		<- log(e/(gamma*price) + 1)
	if(length(e) == R){
		out 	<- sum(psi*gamma*u0 ) 
	}else{
		out 	<- rowSums(psi*gamma*u0 ) 
	}
	
	return(out)
}

uPGrad_fn <- function(e, psi, gamma, price, R){
	u0		<- 1/(e+price*gamma)
	out 	<- psi*gamma*u0 
	return(out)	
}

Allocation_fn <- function(y, psi, gamma, price, R, returnmax = TRUE){
	if(length(y) == 1){
		bu	<- psi/price
		idx	<- order(bu, decreasing = T)
		sorted.bu	<- bu[idx]
		mu	<- 0
		M 	<- 0
		e	<- rep(0, R)
		while(mu/y < sorted.bu[(M+1)] & M < R){
			M	<- M + 1
			sel	<- idx[1:M]
			mu	<- sum(gamma[sel]*psi[sel]) / (1 + sum(gamma[sel]*price[sel]/y))
		}
		sel 	<- idx[1:M]
		e[sel] 	<- gamma[sel]*(psi[sel]*y/mu - price[sel])
	}else{
		N 	<- length(y)
		bu	<- psi/price
		idx	<- t(apply(bu, 1, order, decreasing = T))
		sorted.bu	<- t(apply(bu, 1, sort, decreasing = T))
		mu	<- rep(0, N)
		M	<- rep(0, N)
		for(r in 1:R){
			sel		<- mu/y < sapply(1:N, function(i) sorted.bu[i,(M[i]+1)])
			if(sum(sel) > 0){
				M[sel]	<- M[sel] + 1
				pe.ind	<- t(sapply(1:N, function(i) 1*(1:R %in% idx[i,(1:M[i])])))
				if(sum(sel) == 1){
					mu[sel]	<- sum(gamma[sel,]*psi[sel,]*pe.ind[sel,]) / (1 + sum(gamma[sel,]*price[sel,]*pe.ind[sel,])/y[sel])
				}else{
					mu[sel]	<- rowSums(gamma[sel,]*psi[sel,]*pe.ind[sel,]) / (1 + rowSums(gamma[sel,]*price[sel,]*pe.ind[sel,])/y[sel])
				}
				
			}
		}
		e	<- gamma*(psi*y/mu - price)*pe.ind
	}
	
	if(returnmax){
		omega	<- uP_fn(e, psi, gamma, price, R)
		return(list(e = e, omega = omega))
	}else{
		return(e)
	}
}

incl_value_fn <- function(param_est, nx, base, X1, X2, y, price, R, eps_draw = NULL, returnmax = TRUE, ...){
# param_est		...	parameter estimates
# X_array		...	Attributes array, dimension(R*Nobs*nx)
# y				... A Nobs-dimensional expenditure budget. 
# Q				... Q is the scalor of the quantity constraint.
	beta		<- param_est[1:nx]
	gamma0		<- param_est[(nx+1):(length(param_est)-1)]
	sigma		<- exp(param_est[length(param_est)])
	Nobs		<- length(y)
	gamma		<- sapply(X2, function(x) exp(x%*%gamma0))
	if(is.null(dim(gamma))){
		gamma	<- matrix(gamma, nrow = 1)
	}
	# Compute the parameters for allocation function 
	if(is.null(eps_draw)){
		eps_draw 	<- rep(1, Nobs) %*% t(rgev(R)) * sigma				# matrix(rgev(Nobs*R), Nobs, R)
	}
	xbeta	<- sapply(1:R, function(i) X1[[i]] %*% beta )
	psi		<- exp(xbeta + eps_draw)
	
	# Call allocation function 
	out		<- Allocation_fn(y, psi, gamma,price, R, returnmax)
	return(out)
}	

marginal_fn <- function(param, base, y, Q, X0_list, X1_list, price0, price1, R, Ra, random = FALSE, 
qz_cons = 0, exp_outside = TRUE, quant_outside = TRUE, eps_draw = NULL, inits=NULL, e0=NULL, single=FALSE, ...){
	param_out	<- param_assign(param, nx, R, base)
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
		
	# Compute the parameters for allocation function with two X
	if(random){
		if(is.null(eps_draw)){
			eps_draw 	<- matrix(rgev(Nobs*R), Nobs, R)
		}
	}else{
		eps_draw 	<- matrix(0, Nobs, R)
	}	
	xbeta0	<- sapply(1:R, function(i) X0_list[[i]] %*% beta + beta0[i])
	psi0	<- exp(xbeta0 + eps_draw)
	xbeta1	<- sapply(1:R, function(i) X1_list[[i]] %*% beta + beta0[i])
	psi1	<- exp(xbeta1 + eps_draw)
	
	# Call allocation function 	
	if(is.null(e0)){
		e0		<- matrix(NA, Nobs, R)
		e1		<- matrix(NA, Nobs, R)
		for(i in 1:Nobs){
			sol0	<- Allocation_fn(y = y[i], psi=psi0[i,], gamma=gamma, Q = Q, price=price0[i,], 
					R=R, Ra=Ra, qz_cons = qz_cons, exp_outside = exp_outside, quant_outside = quant_outside, inits=inits, ...)
			e0[i,]	<- sol0$e[1:R]		
			if(!single){
				sol1	<- Allocation_fn(y = y[i], psi=psi1[i,], gamma=gamma, Q = Q, price=price1[i,], 
				R=R, Ra=Ra, qz_cons = qz_cons, exp_outside = exp_outside, quant_outside = quant_outside, inits=inits, ...)
				e1[i,]		<- sol1$e[1:R]
			}	
		}
	}else{
		e1		<- matrix(NA, Nobs, R)
		for(i in 1:Nobs){
			sol1	<- Allocation_fn(y = y[i], psi=psi1[i,], gamma=gamma, Q = Q, price=price1[i,], 
			R=R, Ra=Ra, qz_cons = qz_cons, exp_outside = exp_outside, quant_outside = quant_outside, inits=inits, ...)
			e1[i,]		<- sol1$e[1:R]
		}
	}
	
	# Compute the average difference.
	s0 	<- e0/y
	s1 	<- e1/y
	elas<- colMeans((s1 - s0)/s0)

	return(list(s0 = s0, s1 = s1, elasticity = elas))
}

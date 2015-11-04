
# Bellman initiation 
Bellman_init <- function(param_list,K,GH_init,cheby_init,verbose=F){
	x	<- cheby_init$Adjusted.base.nodes
	nx 	<- nrow(x)
	colnames(x) <- state.name
	num_cheb	<- cheby_init$Cheb.coefficient.num
	m	<- with(param_list, (tau1-beta*log(K)+(1-beta)*lambda1) / (2*tau2-2*lambda2*(1-beta) ))	
	n	<- with(param_list, 2*tau2/(2*tau2 - 2*lambda2*(1-beta)) )
	Phi.bound <- 10e-8
	cat("m=",m,",n=",n,"\n")
	
	x_prim	<- array(NA,c(nx, K, dim(GH_init$nodes)))
	flow_u	<- array(NA,c(nx, K))
	ET		<- array(NA,c(nx, K, num_cheb))
	n.nodes	<- nrow(GH_init$nodes)
	for(i in 1:nx){
		for(k in 1:K){
			sigma 	<- sigma_Q[k]
			a0 		<- -(m/n + mu_Q_fn(k,x[i,y.idx]))/sigma
			b0 		<- -(m/(n-1) + mu_Q_fn(k,x[i,y.idx]))/sigma
			mu 		<- x[i,I.idx] + mu_Q_fn(k,x[i,y.idx])				
			if(pnorm(a0)<Phi.bound){
				EQ_lower	<- 0 
				EQ2_lower	<- 0
			}else{
				rho 		<- dnorm(a0)/pnorm(a0)
				EQ_lower	<- mu - sigma*rho
				EQ2_lower	<- mu^2 - 2*mu*sigma*rho + sigma^2*(1+a0*rho)
			}
			
			if(pnorm(-b0)<Phi.bound){
				EQ_upper	<- 0
				EQ2_upper	<- 0
			}else{
				rho 		<- dnorm(b0)/pnorm(-b0)
				EQ_upper	<- mu + sigma*rho
				EQ2_upper	<- mu^2 + 2*mu*sigma*rho +sigma^2*(1+b0*rho)
			}
			
			p.mid <- pnorm(b0) - pnorm(a0)
			if(p.mid<Phi.bound){
				EQ_mid		<- 0 
				EQ2_mid		<- 0
				Ec_mid		<- 0
				Ec2_mid		<- 0
				EI_mid		<- 0
				EI2_mid		<- 0 
			}else{
				EQ_mid		<- mu - sigma*(dnorm(b0)-dnorm(a0))/p.mid
				EQ2_mid		<- mu^2 - 2*mu*sigma*(dnorm(b0) - dnorm(a0))/p.mid + sigma^2*(1 - (b0*dnorm(b0) - a0*dnorm(a0))/p.mid)
				Ec_mid		<- m*p.mid + n*EQ_mid
				Ec2_mid		<- m^2*p.mid + 2*m*n*EQ_mid + n^2*EQ2_mid
				EI_mid		<- -m*p.mid - (n-1)*EQ_mid
				EI2_mid		<- m^2*p.mid + 2*m*(n-1)*EQ_mid + (n-1)^2*EQ2_mid
			}
			
			# Expected consumption, inventory and flow utility
			Ec 			<- EQ_upper + Ec_mid
			Ec2			<- EQ2_upper + Ec2_mid
			EI_prim		<- EQ_lower + EI_mid
			EI2_prim	<- EQ2_lower + EI2_mid
			flow_u[i,k] <- with(param_list, lambda1*Ec + lambda2*Ec2 - tau1*EI_prim - tau2*EI2_prim) + omega_fn(k,x[i,y.idx])
			
			# The nodes at which to compute integral of value function
			I_prim		<- cbind(-m-(n-1)*(mu + sqrt(2)*sigma*GH_init$nodes[,I.idx]), 
								 mu + sqrt(2)*sigma*GH_init$nodes[,I.idx]	)
			I_prim		<- apply(I_prim,1,function(x) min(max(x[1],0),x[2]) )
			y_prim		<- x[i,y.idx] + sqrt(2)*param_list$sigma_y * GH_init$nodes[,y.idx]
			x_prim[i,k,,I.idx] <- sapply(I_prim, function(x) min(max(x,cheby_init$lower[I.idx]),cheby_init$upper[I.idx]) )
			x_prim[i,k,,y.idx] <- sapply(y_prim, function(x) min(max(x,cheby_init$lower[y.idx]),cheby_init$upper[y.idx]) )
			
			# Re-normalize the nodes to compute the expected basis function
			x_prim_adj	<- 2*(x_prim[i,k,,] - rep(1,n.nodes)%*%t(cheby_init$lower))/
								(rep(1,n.nodes) %*% t(cheby_init$upper - cheby_init$lower)) - 1
			T 			<- apply(x_prim_adj,1,cheb_polynm_multD,cheby_init$Polynomial.degree)
			ET[i,k,] 	<- apply(T,1,function(x) x %*% GH_init$weights)/pi
			
		}
	}
	
	policy_init <- rep(1,nx)
	theta_init	<- rep(0,num_cheb)
	value_init	<- cheby_init$Basis.polynomial %*% theta_init
	
	if(verbose){
		cat("Valuation function initialization\n")
		cat("\nThe state to be evaluated:\n"); print(x)
		cat("\nFlow utility at the current parameter: \n"); print(flow_u)
		# cat("\nThe states to be evaluated in the GH-quadrature:\n")
		# for(i in 1:nx){
		# 	for(k in 1:K){
		# 		cat("State",i,"k=",k,"\n"); print(x_prim[i,k,,])
		# 	}
		# }
		# cat("\nExpected basis function in the value function: \n")
		# for(k in 1:K){
		# 	cat("k=",k,"\n"); print(ET[,k,]);
		# }
	}
	
	list(param_list = param_list,
		state 		= x, 
		next.state 	= x_prim, 
		policy	 	= policy_init,
		Cheb.theta 	= theta_init,
		value.fn	= value_init,
		flow.utility= flow_u,
		ET			= ET,
		consumption = c(m=m,n=n),
		status		= 1
		)
}


# Bellman operator 
Bellman_operator <- function(param_list,DP_list){
	Ew <- apply(DP_list$ET,c(1,2),function(x) x%*% DP_list$Cheb.theta)
	v  <- DP_list$flow.utility + param_list$beta * Ew
	max_v <- apply(v,1,max)
	Tw <- log(rowSums(exp(v-max_v)))
	return(Tw)
}
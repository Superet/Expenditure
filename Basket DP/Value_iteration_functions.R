
# Bellman initiation 
Bellman_init <- function(param_list,K,GH_init,cheby_init){
	x	<- cheby_init$Adjusted.base.nodes
	nx 	<- nrow(x)
	colnames(x) <- state.name
	num_cheb	<- cheby_init$Cheb.coefficient.num
	
	# Use GH-quadrature	
	sigma_c <- abs(0.5*(1-param_list$beta)/(param_list$gamma*(1-param_list$beta) +2*param_list$tau2 )) 
	x_prim	<- array(NA,c(nx,K,dim(GH_init$nodes)))
	flow.u	<- array(NA,c(nx,K))
	ET		<- array(NA,c(nx,K,num_cheb))
	n.nodes <- nrow(GH_init$nodes)	
	for(i in 1:nx){
		for(k in 1:K){
			stock			<- x[i,1] + Q[k]
			mu_c 			<- .5*(param_list$beta*log(K) - param_list$tau1 - 2*param_list$tau2*stock)/
								(param_list$gamma*(1 - param_list$beta) + param_list$tau2)
			Ec				<- ifelse(mu_c<0, 0, ifelse(mu_c>stock,stock,mu_c))
			Ec2				<- ifelse(mu_c<0, 0, ifelse(mu_c>stock,stock^2,mu_c^2+sigma_c^2))
			x_prim[i,k,,1]	<- stock - (mu_c + sqrt(2) * sigma_c * GH_init$nodes[,1])
			x_prim[i,k,,1]	<- ifelse(x_prim[i,k,,1]<0, 0, ifelse(x_prim[i,k,,1]>stock,stock,x_prim[i,k,,1]))
			x_prim[i,k,,2] 	<- x[i,2] + sqrt(2)*param_list$sigma_y * GH_init$nodes[,2]
			x_prim[i,k,,1] 	<- min(max(x_prim[i,k,,1],cheby_init$lower[1]),cheby_init$upper[1])
			x_prim[i,k,,2] 	<- min(max(x_prim[i,k,,2],cheby_init$lower[2]),cheby_init$upper[2])
			
			# Compute the expected consumption utility			
			flow.u[i,k]		<- -param_list$tau1*stock - param_list$tau2*stock^2 - 
								(param_list$gamma + param_list$tau2*(3+param_list$beta)/(1-param_list$beta))*Ec2 + 
								Ec*param_list$beta*(log(K) - param_list$tau1 - 2*param_list$tau2*stock)/(1-param_list$beta)			
			
			# Re-normalize the nodes for Chebychev approximation and compute integrated basis function
			x_prim_adj		<- 2*(x_prim[i,k,,] - rep(1,n.nodes)%*%t(cheby_init$lower))/
								(rep(1,n.nodes) %*% t(cheby_init$upper - cheby_init$lower)) - 1
			T 				<- apply(x_prim_adj,1,cheb_polynm_multD,cheby_init$Polynomial.degree)
			ET[i,k,] 		<- apply(T,1,function(x) x %*% GH_init$weights)/pi
		}
	}
	
	# # Use Monte-Carlo to compute the integral
	# numsim 	<- 100
	# nu_draw <- rnorm(numsim)
	# y_draw	<- t(sapply(x[,2],rnorm,n=numsim,sd=param_list$sigma_y))
	# x_prim	<- array(NA,c(nx,K,numsim,2))
	# flow.u	<- array(NA,c(nx,K))
	# ET		<- array(NA,c(nx,K,num_cheb))
	# for(i in 1:nx){
	# 	for(k in 1:K){
			# stock			<- x[i,1] + Q[k]
	# 		c_draw 			<- .5*(param_list$beta*log(K) - param_list$tau1 - 2*param_list$tau2*(x[i,1]+Q[k]) - (1-param_list$beta)*nu_draw)/
	# 							(param_list$gamma*(1 - param_list$beta) + param_list$tau2)
	# 		c_draw			<- max(min(c_draw,x[i,1] + Q[k]),0)					
	# 		Ec				<- mean(c_draw)
	# 		Ec2 			<- mean(c_draw^2)			
	# 		flow.u[i,k]		<- -param_list$tau1*(x[i,1]+Q[k]) - param_list$tau2*(x[i,1]+Q[k])^2 - 
	# 							(param_list$gamma + param_list$tau2*(3+param_list$beta)/(1-param_list$beta))*Ec2 + 
	# 							Ec*param_list$beta*(log(K) - param_list$tau1 - 2*param_list$tau2*(x[i,1]+Q[k]))/(1-param_list$beta)	
	# 		
	# 		# The state in the next period
	# 		x_prim[i,k,,1]	<- x[i,1] + Q[k] - c_draw
	# 		x_prim[i,k,,2]	<- y_draw[i,]
	# 		x_prim[i,k,,1] 	<- min(max(x_prim[i,k,,1],cheby_init$lower[1]),cheby_init$upper[1])
	# 		x_prim[i,k,,2] 	<- min(max(x_prim[i,k,,2],cheby_init$lower[2]),cheby_init$upper[2])
	# 												
	# 		# Re-normalize the nodes for Chebychev approximation and compute integrated basis function
	# 		x_prim_adj		<- 2*(x_prim[i,k,,] - rep(1,numsim)%*%t(cheby_init$lower))/
	# 							(rep(1,numsim) %*% t(cheby_init$upper - cheby_init$lower)) - 1
	# 		T 				<- apply(x_prim_adj,1,cheb_polynm_multD,cheby_init$Polynomial.degree)
	# 		ET[i,k,] 		<- rowMeans(T)
	# 	}
	# }
		
	policy_init <- rep(1,nx)
	theta_init	<- rep(0,num_cheb)
	value_init	<- cheby_init$Basis.polynomial %*% theta_init
	
	list(state 		= x, 
		next.state 	= x_prim, 
		policy	 	= policy_init,
		Cheb.theta 	= theta_init,
		value.fn	= value_init,
		flow.utility= flow.u,
		ET			= ET,
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
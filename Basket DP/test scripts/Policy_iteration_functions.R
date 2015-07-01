policy_iteration <- function(DP_list,verbose=F){
# Status: 0(converge), 1(initialization), 2(value function converge), 3(value function not converge)
	norm_dif <- set_list$tol + 1
	i 	 <- 1
	if(DP_list$status==1){
		cat("Policy function starts from the initialization.\n")
	}else{
		cat("Policy function continues from previous iteration.\n")
	}
	
	while((norm_dif>=set_list$tol | DP_list$status!=2) & i<= set_list$policy.iter.max){
		DP_list <- policy_eval(DP_list)
		policy_new <- policy_improve(DP_list)
		norm_dif	<- max(abs(policy_new - DP_list$policy))
		DP_list$policy <- policy_new
		
		if(i%%set_list$display.freq==0){
			cat("Policy function iteration,",i,".\n")
		}
		i <- i + 1
	}
	DP_list$Iteration <- i
	if(DP_list$status==2){
		if(norm_dif<=set_list$tol){
			DP_list$status <- 0
			cat("Policy function iteration converges. \n")
		}else{
			cat("Value function has converged at the current policy, but policy function iteration has not converged. \n")
		}
	}else{
		cat("Value function iteration has not converged at the current policy.\n")
	}
	return(DP_list)
}

policy_eval <- function(DP_list){
# Value function given a policy function 
# Update Cheb.theta and value.fn in the DP_list
	param_list <- DP_list$param_list
	policy <- DP_list$policy
	nx <- nrow(DP_list$state)
	norm_dif <- set_list$tol + 1
	i <- 1
	while(norm_dif>=set_list$tol & i<=set_list$value.iter.max){
		Ew <- rep(NA,nx)
		v <- rep(NA,nx)
		for(j in 1:nx){
			Ew[j] <- DP_list$ET[j,policy[j],] %*% DP_list$Cheb.theta
			v[j] <- DP_list$flow.utility[j,policy[j]] + param_list$beta * Ew[j]
		}
		theta.new <- calculate_cheb_coef(cheby_init,v)
		norm_dif <- max(abs(theta.new - DP_list$Cheb.theta))
		DP_list$Cheb.theta 	<- theta.new
		DP_list$value.fn	<- cheby_init$Basis.polynomial %*% theta.new
		
		if(i%%set_list$display.freq==0){
			cat("	Inner valuation function iteration,",i,".\n")
		}
		i <- i + 1
	}
	if(norm_dif>=set_list$tol){
		DP_list$status <- 3
	}else{
		DP_list$status <- 2
	}
	return(DP_list)
}

policy_improve <- function(DP_list){
# Improve policy function given a value function 
	param_list <- DP_list$param_list
	value_fn <- DP_list$value.fn
	Ew <- apply(DP_list$ET,c(1,2),function(x) x%*% DP_list$Cheb.theta)
	v  <- DP_list$flow.utility + param_list$beta * Ew
	policy.new <- apply(v,1,which.max)
	return(policy.new)
}


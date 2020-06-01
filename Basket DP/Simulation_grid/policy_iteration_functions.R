policy_iteration_fn <- function(DP_list,u_fn,transition_fn,policy_improve_fn,ccp_fn=NULL,print_level=0,plot_state=NULL,plot_path=NULL,plot=FALSE){
	# Status: 0(converge), 1(initialization), 2(policy funciton not converged)

	norm_dif	<- set_list$tol + 1
	i 			<- 1
	if(DP_list$status==1){
		cat("Policy function starts from the initialization.\n")
	}else{
		cat("Policy function continues from previous iteration.\n")
	}
	
	pct	<- proc.time()
	plot_data <- data.frame(NULL)
	while(norm_dif>=set_list$tol & i<= set_list$policy_iter_max){
		Ew 			<- policy_eval_fn(DP_list$state,DP_list$policy,DP_list$value_fn,u_fn,transition_fn,DP_list$param_list)
		value_dif	<- max(abs(Ew - DP_list$value_fn))
		DP_list$value_fn <- Ew
		policy_sol 	<- policy_improve_fn(DP_list,plot_state=plot_state,plot)
		policy_new	<- policy_sol$policy
		norm_dif	<- max(abs(unlist(policy_new) - unlist(DP_list$policy)) )
		DP_list$policy <- policy_new
		# In case of singlularity, we use the maximum value from the policy improve as new value funciton
		if(value_dif < 10e-8){
			DP_list$value_fn	<- policy_sol$value
		}
		
		if(i%%set_list$display.freq==0){
			pct1 <- proc.time() - pct
			cat("Policy function iteration,",i,":", pct1[3]/60,"min, norm_dif =", norm_dif,", value_dif=", value_dif,". \n")
			if(print_level == 1){
				cat("----------------------------------\n")
				cat("The new policy funciton for k and the value functions:\n")
				print( data.frame(DP_list$state,policy=DP_list$policy,value=DP_list$value_fn))
				# dynamic_plot_fn(plot_path,DP_list,iteration=i)
			}
			
			if(!is.null(plot_state)){
				plot_data <- rbind(plot_data,data.frame(Iteration=i,policy_sol$plot_data))
			}
			pct <- proc.time()
		}
		i <- i + 1		
	}
	DP_list$Iteration	<- i - 1
	if(!is.null(ccp_fn)){
		DP_list$CCP 		<- ccp_fn(DP_list)
	}	
	
	if(!is.null(plot_state)){
		plot_data <- rbind(plot_data,data.frame(Iteration=i,policy_sol$plot_data))
	}
	if(norm_dif>set_list$tol){
		DP_list$status <- 2
		cat("Policy function iteration has not converged. \n")
	}else{
		DP_list$status <- 3
		cat("Policy function iteration has converged. \n")
	}
	if(is.null(plot_state)){
		return(list(DP_list=DP_list))
	}else{
		return(list(DP_list=DP_list, plot_data=plot_data))
	}
}

policy_eval_fn <- function(state,policy,W,u_fn,transition_fn,param_list){
	ns	<- nrow(state)
	Pi	<- transition_fn(state,policy)
	u 	<- u_fn(state,policy,param_list)
	B	<- diag(1,nrow(Pi)) - beta * Pi
	if(det(B)<=10e-8){
		cat("In policy function, B matrix is singluar, value function is not updated.\n")
		W	<- W
	}else{
		W	<- solve(B) %*% u
	}
	return(W)
}


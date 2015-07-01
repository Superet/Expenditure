value_iteration_fn <- function(DP_list, Bellman_operator, ccp_fn=NULL,print_level=0){
	norm_dif	<- control_list$tol + 1
	i 			<- 1
	if(DP_list$status==1){
		cat("Value function iteration starts from the initialization.\n")
	}else{
		cat("Value function iteration continues from previous iteration.\n")
	}
	
	pct	<- proc.time()
	while(norm_dif>=control_list$tol & i<= control_list$value_max_iter){
		sol 	<- Bellman_operator(DP_list)
		W_new	<- sol$value
		norm_dif<- max(abs(W_new- DP_list$value_fn))
		DP_list$value_fn	<- W_new
		DP_list$policy		<- sol$policy
		if(i%%control_list$display_freq==0){
			pct1 <- proc.time() - pct
			cat("Value function iteration,",i,":", pct1[3]/60,"min, norm_dif =", norm_dif,". \n")
			if(print_level == 1){
				cat("----------------------------------\n")
				cat("The new policy funciton for k and the value functions:\n")
				print( data.frame(DP_list$state,policy=DP_list$policy,value=DP_list$value_fn))
			}
			pct <- proc.time()
		}
		i <- i + 1		
	}
	DP_list$Iteration <- i-1
	if(!is.null(ccp_fn)){
		DP_list$CCP 		<- ccp_fn(DP_list)
	}	
	if(norm_dif>control_list$tol){
		DP_list$status <- 2
		cat("Value function iteration has not converged. \n")
	}else{
		DP_list$status <- 3
		cat("Value function iteration has converged. \n")
	}
	return(DP_list)
}
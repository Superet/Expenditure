DP_likelihood <- function(param,data){
	if(!any(is.na(names(param)))){
		names(param) <- NULL
	}
	cur_param_list <- list(	lambda1	= param[1],
							lambda2	= param[2],
							tau1 	= param[3],
							tau2 	= param[4],
							beta 	= param_list$beta,
							sigma_y = param_list$sigma_y)
	if(cur_param_list$lambda1<0 | cur_param_list$lambda2>0){
		return(NA)
	}else{
		print(unlist(cur_param_list))
	}
							
	# Solve the dynamic solution for the given parameters
	DP_init <- Bellman_init(cur_param_list,K,GH_init,cheby_init,verbose=F) 
	DP_list <- policy_iteration(DP_init,verbose=F)
	
	p <- rep(NA,t)
	for(i in 1:t){
		sol <- choice_probability_once(data$state[i,],DP_list)
		p[i]<- sol$probability[data$choice[i]]
	}
	print(sum(log(p)))
	return(log(p))
}
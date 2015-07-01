# Flow utility function
singleu_fn <- function(c,x,param_list){
	if(x<min(x_grid) ){
		u <- NA
	}else{
		u	<- with(param_list, lambda*log(c) - tau1*(x-c) -tau2*(x-c)^2)
	}
	return(u)
}

flowu_vec_fn <- function(state,policy,param_list){
	ns	<- nrow(state)
	u	<- rep(NA,ns)
	for(i in 1:ns){
		u[i] <- singleu_fn(c=policy[i],x=state[i,1],param_list)
	}
	return(u)
}

# Construct transition matrix given policy function 
transition_fix_fn <- function(state,policy){
	ns	<- nrow(state)
	Pi	<- matrix(0,ns,ns)
	for(i in 1:ns){
		x_next	<- state[i,1] - policy[i] + Q_kernal$nodes
		sel 	<- sapply(x_next, function(x) ifelse(x<min(x_grid),NA,max(which(x_grid<=x)) ))
		sel1	<- which(!is.na(sel))
		sel 	<- setdiff(sel,NA)
		Pi[i,sel] <- Q_kernal$weight[sel1]
		Pi[i,]	<- Pi[i,]/sum(Pi[i,])
	}
	return(Pi)
}

# Policy improvement function
solve_c_fn <- function(DP_list,plot_state=plot_state,plot){
	state	<- DP_list$state
	W		<- DP_list$value_fn
	ns		<- nrow(state)
	Vc 		<- matrix(NA,ns,length(c_grid))
	for(i in 1:length(c_grid)){
		Vc[,i] 	<- flowu_vec_fn(state,policy=rep(c_grid[i],ns),DP_list$param_list)
		x_next	<- as.vector(state) - c_grid[i] + rep(1,ns) %*% t(Q_kernal$nodes)
		for(j in 1:ns){
			x_next[j,] <- sapply(x_next[j,], function(x) min(max(x,min(x_grid)),max(x_grid)))
			sel <- sapply(x_next[j,],function(x) ifelse(x<min(x_grid),NA,max(which(x_grid<=x)) ))
			Vc[j,i] <- Vc[j,i] + beta * (as.vector(W[sel]) %*% Q_kernal$weight )
		}
	}
	c_interior	<- apply(Vc,1,function(x) max(which(x==max(x,na.rm=T))))	# Or use random if tied
	value		<- apply(Vc,1,max,na.rm=T)
	c_star		<- c_grid[c_interior]
	if(!is.null(plot_state)){
		plot_data <- data.frame(state = plot_state,Vc[plot_state,])
		colnames(plot_data) <- c("state",c_grid)
		plot_data <- melt(plot_data,id="state")
		plot_data$max_pnt <- c_star[plot_data$state]
	}	
	return(list(policy = c_star, value = value, plot_data = plot_data))
}
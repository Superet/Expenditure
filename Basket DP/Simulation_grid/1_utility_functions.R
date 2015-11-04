# Flow utility function
singleu_fn <- function(c,I,param_list){
	if(I+Q-c<min(I_grid) ){
		u <- NA
	}else{
		u	<- with(param_list, lambda*log(c+0.01) - tau1*(I+Q-c) -tau2*(I+Q-c)^2)
	}
	return(u)
}

flowu_vec_fn <- function(state,policy,param_list){
	ns	<- nrow(state)
	u	<- rep(NA,ns)
	for(i in 1:ns){
		u[i] <- singleu_fn(c=policy[i],I=state[i,1],param_list)
	}
	return(u)
}

# Construct transition matrix given policy function 
transition_fix_fn <- function(state,policy){
	ns	<- nrow(state)
	Pi	<- matrix(0,ns,ns)
	for(i in 1:ns){
		I_next	<- state[i,] + Q - policy[i]
		sel 	<- max(which(I_grid<=I_next))
		Pi[i,sel] <- 1
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
		I_next	<- state + Q - c_grid[i]
		sel		<- sapply(I_next,function(x) ifelse(x<min(I_grid),NA,max(which(I_grid<=x)) ))
		Vc[,i] 	<- Vc[,i] + beta*W[sel]
	}
	c_interior	<- apply(Vc,1,function(x) max(which(x==max(x,na.rm=T))))
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

# Policy improvement function with Brent's method
solve_c_brent_fn <- function(DP_list,plot_state=NULL,plot=F){
	state	<- DP_list$state
	W		<- DP_list$value_fn
	ns		<- nrow(state)
	cstar	<- rep(NA,ns)
	value	<- rep(NA,ns)
	
	V_fn <- approxfun(state,W)
	value_fn <- function(c,I){
		I_next <- I+Q-c
		out <- singleu_fn(c,I,DP_list$param_list) + beta*V_fn(I_next)
		return(out)
	}
	for(i in 1:ns){
		sol 		<- optimize(value_fn, interval=c(0,state[i,1]+Q), maximum=T, I=state[i,1])
		cstar[i]	<- sol$maximum
		value[i]	<- sol$objective
	}
	return(list(policy = cstar, value=value))
}

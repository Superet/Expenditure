# Flow utility function
singleu_fn <- function(c,I,Q,param_list){
	u	<- with(param_list, lambda*log(c+1) - tau1*(I+Q-c) -tau2*(I+Q-c)^2)
	return(u)
}

flow_Eu_fn <- function(state,policy,param_list){
	ns	<- nrow(state)
	u	<- matrix(NA,ns,nQ)
	nQ	<- length(Q_kernal$nodes)
	for(i in 1:ns){
		for(j in 1:nQ){
			u[i,j] <- singleu_fn(c=policy$c[i,j],I=state[i,1],Q = Q_kernal$nodes[j],param_list)
		}
	}
	Eu	<- rowSums(u * Q_kernal$weight[policy$k,])
	return(Eu)
}

# Construct transition matrix given policy function 
transition_fn <- function(state,policy){
	ns	<- nrow(state)
	Pi	<- matrix(0,ns,ns)
	for(i in 1:ns){
		I_next	<- state[i,1] +Q_kernal$nodes - policy$c[i,]
		sel		<- sapply(I_next, function(x) ifelse(x<min(I_grid),1,max(which(I_grid<=x)) ))
		if(length(unique(sel))==nQ){
			Pi[i,sel]	<- Q_kernal$weight[policy$k[i],]
		}else{
			for(j in 1:nQ){
				Pi[i,sel[j]] <- Pi[i,sel[j]] + Q_kernal$weight[policy$k[i],j]
			}
		}
	}
	return(Pi)
}

# Policy improvement function
solve_kc_fn <- function(DP_list,plot_state=plot_state,plot){
	state	<- DP_list$state
	W		<- DP_list$value_fn
	ns		<- nrow(state)
	nc		<- length(c_grid)
	vc 		<- array(NA,c(ns,nQ,nc),dimnames = list(1:ns,Q_kernal$nodes,c_grid))
	for(i in 1:ns){
		for(j in 1:nQ){
			for(k in 1:nc){
				I_next 	<- state[i,1] + Q_kernal$nodes[j] - c_grid[k]
				sel		<- ifelse(I_next<min(I_grid),NA,max(which(I_grid<=I_next)) )
				if(!is.na(sel)){
					vc[i,j,k]	<- singleu_fn(c=c_grid[k],I=state[i,1],Q=Q_kernal$nodes[j],param_list) + beta*W[sel]
				}
			}
		}
	}
	c_idx	 	<- apply(vc,c(1,2),function(x) ifelse(all(is.na(x)), 1, 
																	 min(which(x==max(x,na.rm=T))) ))
	Vc			<- apply(vc,c(1,2),function(x) ifelse(all(is.na(x)), NA, max(x,na.rm=T)) )
	# V			<- Vc %*% t(Q_kernal$weight)
	V			<- sapply(1:K,function(i) rowSums(Vc * rep(1,ns)%*%t(Q_kernal$weight[i,]), na.rm=T))
	policy_k	<- apply(V, 1, function(x) ifelse(all(is.na(x)), 1,min(which(x==max(x,na.rm=T)))) )
	value		<- apply(V, 1, function(x) ifelse(all(is.na(x)), min(V,na.rm=T)-1, max(x,na.rm=T)))
	c_star		<- sapply(1:nQ,function(i) c_grid[c_idx[,i]])
	dimnames(c_star) <- list(1:ns,Q_kernal$nodes)
	if(!is.null(plot_state)){
		plot_data 		<- melt(vc[plot_state,,])
		names(plot_data)<- c("state_idx","Q","c","value")
		# tmp				<- melt(c_star)
		# names(tmp)		<- c("state_idx","Q","max_pnt")
		# plot_data		<- merge(plot_data,tmp,by=c("state_idx","Q"),all=T)
	}	
	return(list(policy = list(k=policy_k, c=c_star), value = value, plot_data = plot_data))
}
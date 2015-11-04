# Flow utility function
singleu_fn <- function(c,I,Q,param_list){
# 	u	<- with(param_list, lambda*log(c+0.01) - tau1*(I+Q-c) -tau2*(I+Q-c)^2 - tau3*Q)
	u	<- with(param_list, lambda*log(c+0.01) - tau1*(I+Q-c) - tau3*Q )
# 	u	<- with(param_list, lambda*log(c+0.01) - tau1*I - tau3*Q )
	return(u)
}

flowu_vec_fn <- function(state,policy,param_list){
	ns	<- nrow(state)
	u	<- rep(NA,ns)
	for(i in 1:ns){
		u[i] <- singleu_fn(c=policy$c[i,policy$k[i]],I=state[i,1],Q = Q_grid[policy$k[i]],param_list)
	}
	return(u)
}

# Construct transition matrix given policy function 
transition_fix_fn <- function(state,policy){
	ns	<- nrow(state)
	Pi	<- matrix(0,ns,ns)
	for(i in 1:ns){
		I_next	<- state[i,1] + Q_grid[policy$k[i]] - policy$c[i,policy$k[i]]
		sel 	<- ifelse(I_next<min(I_grid),1,max(which(I_grid<=I_next)) )
		Pi[i,sel]	<- 1
	}
	return(Pi)
}

# Policy improvement function
solve_c_fn <- function(DP_list,plot_state=plot_state,plot){
	state	<- DP_list$state
	W		<- DP_list$value_fn
	ns		<- nrow(state)
	nc		<- length(c_grid)
	Vc 		<- array(NA,c(ns,K,nc),dimnames = list(1:ns,1:K,c_grid))
	for(i in 1:ns){
		for(k in 1:K){
			for(j in 1:nc){
				I_next		<- state[i,1] + Q_grid[k] - c_grid[j]
				sel			<- ifelse(I_next<min(I_grid),NA,max(which(I_grid<=I_next)) )
				if(!is.na(sel)){
					Vc[i,k,j]	<- singleu_fn(c_grid[j],state[i,1],Q=Q_grid[k],DP_list$param_list) + beta*W[sel]
				}
			}
		}
	}
	c_idx	 	<- apply(Vc,c(1,2),function(x) ifelse(all(is.na(x)), 1, 
																	 min(which(x==max(x,na.rm=T))) ))
	V			<- apply(Vc,c(1,2),function(x) ifelse(all(is.na(x)), NA, max(x,na.rm=T)) )
	policy_k	<- apply(V, 1, function(x) min(which(x==max(x,na.rm=T))))
	value		<- apply(V,1,max,na.rm=T)
	c_star		<- sapply(1:K,function(i) c_grid[c_idx[,i]])
	if(!is.null(plot_state)){
		plot_data <- melt(Vc[plot_state,,])
		names(plot_data) <- c("state_idx","k","c","value")
	}	
	return(list(policy = list(k=policy_k, c=c_star), value = value, plot_data = plot_data))
}

solve_brent_fn <- function(DP_list,plot_state=NULL, plot=F){
	state	<- DP_list$state
	W		<- DP_list$value_fn
	ns		<- nrow(state)
	c_star	<- matrix(NA, ns, K)
	value	<- matrix(NA, ns, K)
	
# 	V_fn <- approxfun(state,W, rule = 2)
	V_fn <- splinefun(state,W, method = "natural")
	value_fn <- function(c,I,Q){
		I_next <- I+Q-c
		out <- singleu_fn(c,I,Q,DP_list$param_list) + beta*V_fn(I_next)
		return(out)
	}
	
	for(i in 1:ns){
		for(k in 1:K){
			if(state[i,1]+Q_grid[k]<=0){
				c_star[i,k]	<- 0
				value[i,k]	<- value_fn(0,state[i,1],Q_grid[k])
			}else{
				sol 		<- optimize(value_fn, interval=c(0,state[i,1]+Q_grid[k]), maximum=T, 
										I=state[i,1], Q=Q_grid[k], tol = 0.0000001)
				c_star[i,k]	<- sol$maximum
				value[i,k]	<- sol$objective
			}
		}
	}
	policy_k 	<- apply(value,1,which.max)
	value_max	<- apply(value,1,max)
	return(list(policy = list(k = policy_k, c = c_star), value = value_max))
}


plot_cvalue_fn <- function(DP_list, c_grid, plot_state=NULL, Q_index=NULL, plot = TRUE){
	state	<- DP_list$state
	W		<- DP_list$value_fn
	nc		<- length(c_grid)
	if(is.null(plot_state)) {plot_state <- 1:nrow(state) }
	if(is.null(Q_index))	{Q_index	<- 1:length(Q_grid)}
	ns		<- length(plot_state)
	nQ 		<- length(Q_index)
	
	c_star	<- as.matrix(DP_list$policy$c[plot_state,Q_index])
	dimnames(c_star) <- list(plot_state, Q_grid[Q_index])
	vc		<- array(NA, c(ns, nQ, nc+1), 
				dimnames = list(plot_state, Q_grid[Q_index], c(c_grid, "")) )
	
	# NOTE: extrapolate the value function outside the grid.
	V_fn <- approxfun(state,W)
	value_fn <- function(c,I,Q){
		I_next <- I+Q-c
		out <- singleu_fn(c,I,Q,DP_list$param_list) + beta*V_fn(I_next)
		return(out)
	}
	
	for(i in 1:ns){
		for(j in 1:nQ){
			for(k in 1:nc){
				vc[i,j,k] <- value_fn(c_grid[k], state[plot_state[i],1], Q_grid[Q_index[j]])
			}
			vc[i,j,(nc+1)] <- value_fn(c_star[i,j], state[plot_state[i],1], Q_grid[Q_index[j]])
		}
	}
	
	ggtmp 		<- melt(vc)
	sel 		<- is.na(ggtmp$Var3)
	ggtmp[sel,"Var3"] <- melt(c_star)$value
	ggtmp$max_pnt <- 0
	ggtmp[sel,"max_pnt"] <- 1
	names(ggtmp)[1:3]	<- c("s_index", "Q", "c")
	ggtmp$I		<- state[ggtmp$s_index, 1]
		
	if(plot){
		for(i in 1:nQ){
			quartz()
			print(ggplot(subset(ggtmp, Q==Q_grid[Q_index[i]]), aes(c, value)) + geom_point(aes(col=factor(max_pnt))) + 
						geom_line() + 
						facet_wrap(~I) + 
						scale_color_manual(values=c("black","red")) + 
						guides(color=FALSE) + 
						labs(title = paste("Value function curve at Q=",Q_grid[Q_index[i]], "\n(facet by I)", sep=""))
			)
		}
	}else{
		return(ggtmp)
	}
}

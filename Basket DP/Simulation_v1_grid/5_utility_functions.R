# Flow utility function
singleu_fn <- function(c,I,Q,param_list){
	u	<- with(param_list, lambda*log(c+.01) - tau1*(I+Q-c) - tau3*Q)
	return(u)
}

flow_Eu_fn <- function(state,policy,param_list){
	ns	<- nrow(state)
	u	<- matrix(NA,ns,nQ)
	Eu	<- rep(NA,ns)
	nQ	<- length(Q_kernal$nodes)
	for(i in 1:ns){
		for(j in 1:nQ){
			u[i,j] <- singleu_fn(c=policy$c[i,j],I=state[i,I_idx],Q = Q_kernal$nodes[j],param_list)
		}
		Eu[i] <- sum(u[i,] * Q_kernal$weight[[s_idx[i,y_idx]]][policy$k[i],])
	}
	return(Eu)
}

# Construct transition matrix given policy function 
transition_fn <- function(state,policy){
	ns	<- nrow(state)
	Pi	<- matrix(0,ns,ns)
	for(i in 1:ns){
		I_next	<- state[i,I_idx] + Q_kernal$nodes - policy$c[i,]
		sel		<- sapply(I_next, function(x) ifelse(x<min(I_grid),1,max(which(I_grid<=x)) ))
		for(j in 1:nQ){
			sel1	<- state[,I_idx]==I_grid[sel[j]]
			Pi[i,sel1] <- Pi[i,sel1] + Q_kernal$weight[[s_idx[i,y_idx]]][policy$k[i],j]
		}
		for(j in 1:ny){
			sel1	<- state[,y_idx] == y_grid[j]
			Pi[i,sel1] <- Pi[i,sel1] * y_kernal$weight[s_idx[i,y_idx],j]
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
				I_next 	<- state[i,I_idx] + Q_kernal$nodes[j] - c_grid[k]
				sel		<- ifelse(I_next<min(I_grid),NA,max(which(I_grid<=I_next)) )
				if(!is.na(sel)){
					sel1		<- state[,I_idx] == I_grid[sel]
					vc[i,j,k]	<- singleu_fn(c=c_grid[k],I=state[i,I_idx],Q=Q_kernal$nodes[j],DP_list$param_list) + 
									beta* sum(W[sel1] * y_kernal$weight[s_idx[i,y_idx],])
				}
			}
		}
	}
	c_idx	 	<- apply(vc,c(1,2),function(x) ifelse(all(is.na(x)), 1, 
																	 min(which(x==max(x,na.rm=T))) ))
	Vc			<- apply(vc,c(1,2),function(x) ifelse(all(is.na(x)), NA, max(x,na.rm=T)) )
	# V			<- Vc %*% t(Q_kernal$weight)
	V			<- t(sapply(1:ns,function(i) Q_kernal$weight[[s_idx[i,y_idx]]] %*% Vc[i,] ))
	policy_k	<- apply(V, 1, function(x) ifelse(all(is.na(x)), 1,min(which(x==max(x,na.rm=T)))) )
	value		<- apply(V, 1, function(x) ifelse(all(is.na(x)), min(V,na.rm=T)-1, max(x,na.rm=T)))
	c_star		<- sapply(1:nQ,function(i) c_grid[c_idx[,i]])
	dimnames(c_star) <- list(1:ns,Q_kernal$nodes)
	if(!is.null(plot_state)){
		plot_data 		<- melt(vc[plot_state,,])
		names(plot_data)<- c("state_idx","Q","c","value")
		tmp				<- melt(c_star)
		names(tmp)		<- c("state_idx","Q","max_pnt")
		plot_data		<- merge(plot_data,tmp,by=c("state_idx","Q"),all.x=T)
	}else{
		plot_data		<- NULL
	}	
	return(list(policy = list(k=policy_k, c=c_star), value = value, plot_data = plot_data))
}


# Policy improvement with Brent's method (Bellman operator as well)
solve_c_brent_fn <- function(DP_list, inter_spline = TRUE){
	state	<- DP_list$state
	W		<- DP_list$value_fn
	ns		<- nrow(state)
	nQ		<- length(Q_kernal$nodes)
	c_star	<- matrix(NA, ns, nQ)
	vc		<- matrix(NA, ns, nQ)
	if(!is.null(control_list$inter_spline)){
		inter_spline <- control_list$inter_spline
	}
	
	# NOTE: extrapolate the value function outside the grid.
	if(inter_spline){
		V1_fn <- splinefun(state[state[,"y"]==1,I_idx],W[state[,"y"]==1], method = "natural")
		V2_fn <- splinefun(state[state[,"y"]==2,I_idx],W[state[,"y"]==2], method = "natural")	
	}else{
		V1_fn <- approxfun(state[state[,"y"]==1,I_idx],W[state[,"y"]==1], rule = 2)
		V2_fn <- approxfun(state[state[,"y"]==2,I_idx],W[state[,"y"]==2], rule = 2)	
	}
	value_fn <- function(c,I,Q,sel_y){
		I_next 	<- I+Q-c
		v		<- c(V1_fn(I_next), V2_fn(I_next))
		out 	<- singleu_fn(c,I,Q,DP_list$param_list) + beta*y_kernal$weight[sel_y,] %*% v
		return(out)
	}
	
	policy_k <- rep(NA, ns)
	value	 <- rep(NA, ns)	
	for(i in 1:ns){
		for(j in 1:nQ){
			if(state[i,1]+Q_kernal$nodes[j]<=0){
				c_star[i,j]	<- 0
				vc[i,j]		<- value_fn(0,state[i,1],Q_kernal$nodes[j], sel_y=s_idx[i,y_idx] )
			}else{
				sol 		<- optimize(value_fn, interval=c(0,state[i,1]+Q_kernal$nodes[j]), 
								maximum=T, I=state[i,1], Q=Q_kernal$nodes[j], sel_y=s_idx[i,y_idx],
								tol = control_list$brent_tol)
				c_star[i,j]	<- sol$maximum
				vc[i,j]		<- sol$objective
			}		
		}
		V 		<- Q_kernal$weight[[s_idx[i,y_idx]]] %*% vc[i,] 
		V_max 	<- max(V)
		policy_k[i] <- which.max(V)
		value[i] <- log(sum(exp(V - V_max))) + V_max
	}
	return(list(policy = list(k = policy_k, c = c_star), value = value))
}

plot_cvalue_fn <- function(DP_list, c_grid, plot_state=NULL, Q_index=NULL, control_list, plot = TRUE){
	state	<- DP_list$state
	W		<- DP_list$value_fn
	nc		<- length(c_grid)
	if(is.null(plot_state)) {plot_state <- 1:nrow(state) }
	if(is.null(Q_index))	{Q_index		<- Q_kernal$nodes }
	ns		<- length(plot_state)
	nQ 		<- length(Q_index)
	
	c_star	<- as.matrix(DP_list$policy$c[plot_state,Q_index])
	dimnames(c_star) <- list(plot_state, Q_kernal$nodes[Q_index])
	vc		<- array(NA, c(ns, nQ, nc+1), 
				dimnames = list(plot_state, Q_kernal$nodes[Q_index], c(c_grid, "")) )
	
	# NOTE: extrapolate the value function outside the grid.
	inter_spline <- control_list$inter_spline
	if(inter_spline){
		V1_fn <- splinefun(state[state[,"y"]==1,I_idx],W[state[,"y"]==1], method = "natural")
		V2_fn <- splinefun(state[state[,"y"]==2,I_idx],W[state[,"y"]==2], method = "natural")	
	}else{
		V1_fn <- approxfun(state[state[,"y"]==1,I_idx],W[state[,"y"]==1], rule = 1)
		V2_fn <- approxfun(state[state[,"y"]==2,I_idx],W[state[,"y"]==2], rule = 1)	
	}
	value_fn <- function(c,I,Q,sel_y){
		if(I+Q-c<0){
			return(NA)
		}else{
			I_next 	<- I+Q-c
			v		<- c(V1_fn(I_next), V2_fn(I_next))
			out 	<- singleu_fn(c,I,Q,DP_list$param_list) + beta*y_kernal$weight[sel_y,] %*% v
			return(out)		
		}
	}
	
	for(i in 1:ns){
		for(j in 1:nQ){
			for(k in 1:nc){
				vc[i,j,k] <- value_fn(c_grid[k], state[plot_state[i],I_idx], 
					Q_kernal$nodes[Q_index[j]], sel_y = s_idx[plot_state[i],y_idx])
			}
			vc[i,j,(nc+1)] <- value_fn(c_star[i,j], state[plot_state[i],I_idx], 
					Q_kernal$nodes[Q_index[j]], sel_y = s_idx[plot_state[i],y_idx])
			
		}
	}
	
	ggtmp 		<- melt(vc)
	sel 		<- is.na(ggtmp$Var3)
	ggtmp[sel,"Var3"] <- melt(c_star)$value
	ggtmp$max_pnt <- 0
	ggtmp[sel,"max_pnt"] <- 1
	names(ggtmp)[1:3]	<- c("s_index", "Q", "c")
	ggtmp$I		<- state[ggtmp$s_index, I_idx]
	ggtmp$y		<- state[ggtmp$s_index, y_idx]
	
	
	if(plot){
		for(i in 1:nQ){
			quartz()
			print(ggplot(subset(ggtmp, Q==Q_kernal$nodes[Q_index[i]]), aes(c, value)) + geom_point(aes(col=factor(max_pnt))) + 
						geom_line() + 
						facet_wrap(I ~ y) + 
						scale_color_manual(values=c("black","red")) + 
						guides(color=FALSE) + 
						labs(title = paste("Value function curve at Q=",Q_kernal$nodes[Q_index[i]], "\n(facet by I,y)", sep=""))
			)
		}
	}else{
		return(ggtmp)
	}
}

# ----------------------------------------------------------------------------------------- #
##############
# Simulation # 
##############
CCP_single_fn <- function(state_single, DP_list){
	V	<- rep(NA, nQ)
	sel	<- ifelse(state_single[I_idx]<min(I_grid), 1, max(which(I_grid<=state_single[I_idx])))
	idx	<- which(s_idx[,I_idx]==sel & DP_list$state[,y_idx]==state_single[y_idx])
	for(i in 1:nQ){
		I_next 	<- state_single[I_idx] + Q_kernal$nodes[i] - DP_list$policy$c[idx,i]
		sel		<- ifelse(I_next<min(I_grid),1,max(which(I_grid<=I_next)) )
		sel1	<- DP_list$state[,I_idx] == I_grid[sel]
		V[i]	<- singleu_fn(c=DP_list$policy$c[idx,i],I=DP_list$state[idx,I_idx],Q=Q_kernal$nodes[i],DP_list$param_list) + 
						beta* sum(DP_list$value_fn[sel1] * y_kernal$weight[s_idx[idx,y_idx],])
	}
	EV	<- Q_kernal$weight[[s_idx[idx,y_idx]]] %*% V
	EV_max <- max(EV)
	EV	<- EV - EV_max
	CCP	<- exp(EV)/sum(exp(EV))
	return(CCP)
}

CCP_fullstate_fn <- function(DP_list){
	CCP <- matrix(NA, ns, K)
	for(i in 1:ns){
		CCP[i,]	<- CCP_single_fn(state[i,], DP_list)
	}
	return(CCP)
}

simulate_sequence_fn <- function(init_state, T, DP_list){
	state	<- matrix(NA, T, ncol(DP_list$state))
	colnames(state) <- colnames(DP_list$state)
	Q	 	<- rep(NA, T)
	c		<- rep(NA, T)
	k		<- rep(NA, T)
	ccp		<- matrix(NA, T, K)
	
	state[1,] <- init_state
	for(i in 1:T){	
		sel	<- ifelse(state[i,I_idx]<min(I_grid), 1, max(which(I_grid<=state[i,I_idx])))
		idx	<- which(s_idx[,I_idx]==sel & DP_list$state[,y_idx]==state[i,y_idx])
		ccp[i,]		<- CCP_single_fn(state[i,], DP_list)
		k[i]	<- sample(1:K, 1, prob = ccp[i,])
		Q[i]	<- sample(Q_kernal$nodes, 1, prob = Q_kernal$weight[[s_idx[idx,y_idx]]][k[i],])
		c[i]	<- DP_list$policy$c[idx,which(Q_kernal$nodes==Q[i])]
		
		if(i<T){
			state[(i+1),I_idx] <- state[i,I_idx] + Q[i] - c[i]			
			state[(i+1),y_idx] <- sample(y_kernal$nodes, 1, prob = y_kernal$weight[which(y_kernal$nodes==state[i,y_idx]),])
		}	
	}
	return(cbind(t = 1:T, state, Q = Q, c=c, k=k, ccp=ccp))
}

# ----------------------------------------------------------------------------------------- #
##############
# Estimation # 
##############
ll_fullstate_fn <- function(theta, state, choice_seq, DP_init){
	DP_list 	<- DP_init
	DP_list$param_list 	<- list(lambda 	= theta[1],
								tau1 	= theta[2],
								tau2	= theta[3],
								tau3	= theta[4]
								)
	sol 	<- policy_iteration_fn(DP_list,u_fn=flow_Eu_fn,transition_fn=transition_fn,policy_improve_fn=solve_kc_fn,
								ccp_fn=CCP_fullstate_fn, print_level=0)
	DP_list	<- sol$DP_list
	ccp		<- DP_list$CCP
	idx		<- rep(NA,nrow(state))
	for(i in 1:nrow(state)){
		sel		<- ifelse(state[i,I_idx]<min(I_grid), 1, max(which(I_grid<=state[i,I_idx])))
		idx[i]	<- which(s_idx[,I_idx]==sel & DP_list$state[,y_idx]==state[i,y_idx])
	}
	ll		<- sapply(1:nrow(state), function(i) log(ccp[idx[i],choice_seq[i]]))
	return(ll)
}

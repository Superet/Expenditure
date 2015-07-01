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

# Policy improvement with Brent's method (Bellman operator as well)
solve_c_brent_fn <- function(DP_list, inter_spline = TRUE){
	state	<- DP_list$state
	W		<- DP_list$value_fn
	ns		<- nrow(state)
	nQ		<- length(Q_kernal$nodes)
	c_star	<- array(NA, c(ns, K, nQ))
	vc		<- array(NA, c(ns, K, nQ))
	V		<- matrix(NA, ns, K)
	Qsigma 	<- Q_kernal$sigma
	Qmu		<- Q_kernal$mu
	value	<- rep(NA, ns)
	
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
	
	# Compute the optimal consumption at given discrete choice and quantity realization;
	for(i in 1:ns){
		sel_y <- state[i, y_idx]
		for(k in 1:K){
			if(Qsigma[sel_y,k]==0){
				tmpQ <- 0 
				if(state[i,I_idx]+tmpQ <=0){
					c_star[i,k,]	<- 0
					vc[i,k,]		<- value_fn(0,state[i,1], tmpQ, sel_y=sel_y)
				}else{
					sol 		<- optimize(value_fn, interval=c(0,state[i,1]+tmpQ), 
								maximum=T, I=state[i,1], Q=tmpQ, sel_y=sel_y,
								tol = control_list$brent_tol)
					c_star[i,k,]	<- sol$maximum
					vc[i,k,]		<- sol$objective
				}
				V[i,k]	<- vc[i,k,1]		
			}else{
				for(j in 1:nQ){
					tmpQ <- exp(sqrt(2)*Qsigma[sel_y,k]*Q_kernal$nodes[j] + Qmu[sel_y, k])
					if(state[i,I_idx]+tmpQ <=0){
						c_star[i,k,j]	<- 0
						vc[i,k,j]		<- value_fn(0,state[i,1], tmpQ, sel_y=sel_y)
					}else{
						sol 		<- optimize(value_fn, interval=c(0,state[i,1]+tmpQ), 
									maximum=T, I=state[i,1], Q=tmpQ, sel_y=sel_y,
									tol = control_list$brent_tol)
						c_star[i,k,j]	<- sol$maximum
						vc[i,k,j]		<- sol$objective
					}			
				}
				V[i,k] <- t(vc[i,k,]) %*% Q_kernal$weight / sqrt(pi)
			}
		}
		tmp_vmax <- max(V[i,])
		value[i] <- log(sum(exp(V[i,] - tmp_vmax))) + tmp_vmax
	}
	
	policy_k	<- apply(V, 1, which.max)	
	return(list(policy = list(k = policy_k, c = c_star), value = value))
}

plot_cvalue_fn <- function(DP_list, c_grid, plot_state=NULL, k=1, control_list, plotnodes = NULL,plot = TRUE){
	state	<- DP_list$state
	W		<- DP_list$value_fn
	nc		<- length(c_grid)
	if(is.null(plot_state)) {plot_state <- 1:nrow(state) }
	
	Q_nodes <- Q_kernal$nodes
	mu		<- Q_kernal$mu
	sigma	<- Q_kernal$sigma
	nQ 		<- length(Q_nodes)
	Q_index <- seq((k-1)*nQ+1, k*nQ)
	ns		<- length(plot_state)
	sel_y 	<- unique(state[plot_state,y_idx])
	if(length(sel_y)>1){
		cat("There are more than 1 level in y, choose another plots.\n")
	}
	Q 		<- exp(sqrt(2)*sigma[sel_y,k]*Q_nodes + mu[sel_y, k])
	
	c_star	<- as.matrix(DP_list$policy$c[plot_state,Q_index])
	dimnames(c_star) <- list(plot_state, Q)
	vc		<- array(NA, c(ns, nQ, nc+1), 
				dimnames = list(plot_state, Q, c(c_grid, "")) )
	
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
					Q[j], sel_y = s_idx[plot_state[i],y_idx])
			}
			vc[i,j,(nc+1)] <- value_fn(c_star[i,j], state[plot_state[i],I_idx], 
					Q[j], sel_y = s_idx[plot_state[i],y_idx])
			
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
	
	if(is.null(plotnodes)) {plotnodes <- 1:nQ }
	tmpQ <- unique(ggtmp$Q)[plotnodes]
	if(plot){
		for(i in 1:length(plotnodes)){
			quartz()
			print(ggplot(subset(ggtmp, Q==tmpQ[i]), aes(c, value)) + geom_point(aes(col=factor(max_pnt))) + 
						geom_line() + 
						facet_wrap(I ~ y) + 
						scale_color_manual(values=c("black","red")) + 
						guides(color=FALSE) + 
						labs(title = paste("Value function curve at Q=",round(tmpQ[i],2), "\n(facet by I,y)", sep=""))
			)
		}
	}else{
		return(ggtmp)
	}
}

simulate_seq_fn <- function(init_state, TT, draw_all = FALSE, y_seq, Q_seq, k_seq, DP_list, control_list){
	state 		<- DP_list$state
	Q_kernal 	<- DP_list$Q_kernal
	y_kernal	<- DP_list$y_kernal
	W			<- DP_list$value_fn
	V1			<- V[DP_list$state[,y_idx]==1]
	V2 			<- V[DP_list$state[,y_idx]==2]
	
	# Value function 
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
	
	# The output object
	DataState	<- matrix(NA, TT, ncol(DP_list$state))
	names(DataState) <- names(DP_list$state)
	if(draw_all){
		Q_seq <- rep(NA, TT)
		k_seq <- rep(NA, TT)
	}else{
		DataState[,y_idx] <- y_seq
	}
	c_seq		<- rep(NA, TT)
	ccp			<- matrix(NA, TT, K)
	
	# Solve for the optimal choice at each state
	DataState[1,] <- init_state
	for(i in 1:TT){
		# If the choice k and Q is to be simulated, we solve for k first, and simulate Q conditionally;
		# Otherwise, we directly solve for optimal consumption level that determines the state in the next period. 
		sel_y <- DataState[i,2]
		# Call c++ function that computes CCP for a given state
		tmp_sol <- CCP_fnC(matrix(DataState[i,],nrow=1), I_grid,V1,I_grid, V2,K, DP_list$param, DP_list$beta, 
								Q_kernal, y_kernal,makedraw=T,control_list)
		ccp[i,]	<- tmp_sol$ccp
		if(draw_all){
			k_seq[i]<- tmp_sol$draw
			Q_seq[i]<- rlnorm(1, Q_kernal$mu[sel_y, k_seq[i]], Q_kernal$sigma[sel_y, k_seq[i]])
		}
		tmp_sol <- optimize(value_fn, interval=c(0,DataState[i,1]+Q_seq[i]), 
					maximum=T, I=DataState[i,1], Q=Q_seq[i], sel_y=sel_y,
					tol = control_list$brent_tol)
		c_seq[i]<- tmp_sol$maximum
		
		# Detemine the state of next period
		if(i<TT){
			DataState[(i+1),I_idx] <- DataState[i,I_idx] + Q_seq[i] - c_seq[i]
			if(draw_all){
				DataState[(i+1),y_idx] <- sample(y_kernal$nodes, 1, prob = y_kernal$weight[sel_y,])
			}
		}	
	}
	return(list(data.frame(	t=1:TT, I=DataState[,I_idx], y=DataState[,y_idx], k = k_seq, Q = Q_seq, 
						c = c_seq),
				ccp = ccp) )
}


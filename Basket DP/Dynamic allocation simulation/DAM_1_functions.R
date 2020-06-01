#================================================================================================================================#
# Dynamic allocation model 

# Model: 
# V(I,y,nu) = max_{Q,c}[E_{epsilon}(max_{e_s} \sum_{s=1}^{S+1} uP(e_s,epsilon_s)) + uC(c) - C(I') + nu + beta*E_{I',y',nu'}V(I',y',nu')]
# 		s.t. sum_{s=1}^{S+1} <= y, sum_{s=1}^{S} e_s/p_s + q_{S+2} <= Q, c<=I+Q, Q>=0
# where 
# uP		...	the purchase utility, uP(e_s, epsilon_s) = gamma_s*psi_s*exp(epsilon_s)*log(e_s/(p_s*gamma_s)+1); 
# 			For identification reason, we set epsilon_{S+1} = epsilon_{S+1} = 0
# 			Parameterization of psi: psi_s = exp(X*beta_s)
# uC		... consumption utility. uC(c) = lambda * log(c+0.01)
# C(I')	... inventory cost. C(I) = tau1*I
# 
# Transition:
# I'	= I + Q - c
# y'	= pi(y)		first order Markorvian process. 
# 
# Two-step calculation 
# Step 1: inclusive value of expected purchase utility
# omega(y,Q) = E_{epsilon}(max_{e_s} \sum_{s=1}^{S+1} uP(e_s,epsilon_s))
# 			s.t. sum_{s=1}^{S+1} <= y, sum_{s=1}^{S} e_s/p_s + q_{S+2} <= Q
# 			
# Step 2: Dynamic programming of basket quantity and consumption 
# V(I,y,nu) = max_{Q,c}[omega(y,Q) + uC(c) - C(I') + nu + beta*E_{I',y',nu'}V(I',y',nu')]
# 		s.t. c<=I+Q, Q>=0
# 
# Versions: 
# This is V1: y is discrete, Q is discrete, gamma's are set to be 1. 
#================================================================================================================================#

# The functions contained in this file: 
# uP_fn(e,psi,gamma,price,S)					... (Negative) purchase utility funciton
# uflow_fn(c, I, Q, sel_y, sel_k, omega, param)	... Flow utility funciton
# Allocation_fn(y,Q,psi,gamma,price,S)			... Compute constrained allocation problem 
# Solve_DP_brent_fn(V_old, state, K, Q_grid, omega, param, beta, control_list,DataState)				
#												... Bellman operator
# Bellman_operator(DP_list, DataState)			... a wrapper function of Bellman operator
# plot_cvalue_fn(DP_list, c_grid, plot_state=NULL, Q_index=NULL, k_index = 1, control_list, plot = TRUE)
#												... A function plots the value function along consumption values 
# simulate_seq_fn(init_state, TT, draw_all = FALSE, y_seq, Q_seq, k_seq, DP_list, control_list, alpha, gamma, omega)
#												... A funciton simulates the sequence of Q,e,c given the initial state. 

#================================================================================================================================#

# Purchae utility function 
uP_fn <- function(e, psi, gamma,price, S){
# Negative purchase utility 
# e: S+2 vecotr
# psi, gamma, price is (S+2)-dimension vector 
	e0 		<- e[1:S]
	e_S1	<- e[(S+1):length(e)]
	u0		<- c(log(e0/(gamma[1:S]*price) + 1), log(e_S1) )
	out 	<- sum(psi*gamma*u0 ) 
	return(-out)
}

uPGrad_fn <- function(e, psi, gamma,price, S){
	e0 		<- e[1:S]
	e_S1	<- e[(S+1):length(e)]
	u0		<- c(1/(e0+price*gamma[1:S]), 1/e_S1 )
	out 	<- psi*gamma*u0 
	return(-out)	
}

# Flow utility function
uflow_fn <- function(c, I, Q, sel_y, sel_k, omega, param, nu=0){
	lambda	<- param[1]
	tau1  	<- param[2]
	u		<- lambda*log(c+.01+nu) - tau1*(I+Q-c) + omega[sel_y, sel_k]
	return(u)
}

# Allocation function - solve a constrained optimization
Allocation_fn <- function(y, psi, gamma, Q=NULL, price,S,Sa, inits=NULL, silent=FALSE){
	# We have Sa non-negativity constraints, one budget constraint, and one quantity constraint
	# ui <- rbind(diag(S+2), c(rep(-1, S+1), 0) )
	# inits is a list of initial values
	ui <- rbind(diag(Sa), c(rep(-1, S+1),rep(0, Sa-S-1 )) )
	ci <- c(rep(0,Sa), -y)
	if(is.null(inits)){
		inits <- list(rep(.01, Sa))
	}
		
	if(is.null(Q)){
		sol <- constrOptim(theta = inits, f=uP_fn, grad=NULL, ui=ui, ci=ci,psi=psi, gamma=gamma, price=price, S=S)
		return(list(e = sol$par, max = -sol$value))
	}else{
		if(Q<=0){
			e <- c(rep(0,S),y)
			if(Sa>S+1){
				e <- c(e, 0)
			}
			return(list(e = e, max = log(y) ) )
		}else{
			ui <- rbind(ui, c(-1/price, 0, -1))
			ci <- c(ci, -Q)
			sol.list <- vector("list", length(inits))
			for(j in 1:length(inits)){
				sol.list[[j]] <- try(constrOptim(theta = inits[[j]], f=uP_fn, grad=uPGrad_fn, ui=ui, ci=ci,psi=psi,
									 gamma=gamma, price=price, S=S), silent=silent)
			}
			sol.list1 <- sol.list[sapply(sol.list, function(x) !inherits(x, "try-error"))]
			sel <- which.min(sapply(sol.list1, function(x) x$value))
			if(length(sel)==0){
				sol <- list(par = rep(NA, Sa), value=NA, convergence = NA)
			}else{
				sol <- sol.list1[[sel]]
				if(sol$convergence != 0){
					cat("Constrained optimization does not converge at value y=",y,", Q=",Q,"\n")
				}
			}
			return(list(e = sol$par, max = -sol$value, convergence = sol$convergence))
		}
	}
}

# Bellman operator: use brent's method to solve for optimal consumption
Solve_DP_brent_fn <- function(V_old, state, K, Q_grid, omega, param, beta, control_list,DataState=NULL){
	ns		<- nrow(state)
	nQ		<- length(Q_grid)
	inter_spline <- control_list$inter_spline
	
	# Value function interpolation: cspline or linear
	if(inter_spline){
		V1_fn <- splinefun(state[state[,"y"]==1,I_idx],V_old[state[,"y"]==1], method = "natural")
		V2_fn <- splinefun(state[state[,"y"]==2,I_idx],V_old[state[,"y"]==2], method = "natural")	
	}else{
		# NOTE: extrapolate the value function outside the grid.
		V1_fn <- approxfun(state[state[,"y"]==1,I_idx],V_old[state[,"y"]==1], rule = 2)
		V2_fn <- approxfun(state[state[,"y"]==2,I_idx],V_old[state[,"y"]==2], rule = 2)	
	}
	value_fn <- function(c,I,Q,sel_y,sel_k){
		I_next 	<- I+Q-c
		v		<- c(V1_fn(I_next), V2_fn(I_next))
		out 	<- uflow_fn(c, I, Q, sel_y, sel_k, omega, param) + beta*y_kernal$weight[sel_y,] %*% v
		return(out)
	}
	
	# Change to data state if it is given 
	if(!is.null(DataState)){
		state 	<- DataState
		ns 		<- nrow(state)
	}
	
	# Assign the objects that store results from loops
	c_star	<- array(NA, c(ns, K))
	vc		<- array(NA, c(ns, K))
	V		<- matrix(NA, ns)
	ccp 	<- matrix(NA, ns, K)
	
	# Maximization loop over each state 
	for(i in 1:ns){
		sel_y <- state[i,y_idx]
		for(k in 1:K){
			if(state[i,1]+Q_grid[k]<=0){
				c_star[i,k]	<- 0
				vc[i,k]		<- value_fn(0,state[i,1],Q_grid[k], sel_y, sel_k = k)
			}else if(Q_grid[k]>y_value[sel_y]/min(price)){
				vc[i,k]		<- (-100)
				c_star[i,k]	<- NA			
			}else{
				sol 		<- optimize(value_fn, interval=c(0,state[i,1]+Q_grid[k]), 
								maximum=T, I=state[i,1], Q=Q_grid[k], sel_y=sel_y,sel_k = k,
								tol = control_list$brent_tol)
				c_star[i,k]	<- sol$maximum
				vc[i,k]		<- sol$objective
			}
		}
		tmp_vmax <- max(vc[i,])
		V[i] 	 <- log(sum(exp(vc[i,] - tmp_vmax))) + tmp_vmax
		ev		 <- exp(vc[i,] - tmp_vmax)
		ccp[i,]	 <- (ev)/sum(ev)
	}
	policy_k <- apply(V, 1, which.max)
	return(list(policy = list(k = policy_k, c = c_star), value = V, ccp = ccp))
}


# Bellman operator wrapper
Bellman_operator <- function(DP_list, DataState=NULL){
	state	<- DP_list$state
	K		<- DP_list$K
	Q_grid 	<- DP_list$Q_grid
	omega	<- DP_list$omega
	param	<- DP_list$param
	beta	<- DP_list$beta
	V_old	<- DP_list$value_fn
	out <- Solve_DP_brent_fn(V_old, state, K, Q_grid, omega, param, beta, control_list,DataState)
	return(out)
}

# Check the concavity of value function when solving for c; 
# Plot the value function a long a grid of c.
plot_cvalue_fn <- function(DP_list, c_grid, plot_state=NULL, Q_index=NULL, control_list, plot = TRUE){
	state	<- DP_list$state
	V		<- DP_list$value_fn
	nc		<- length(c_grid)
	if(is.null(plot_state)) {plot_state <- 1:nrow(state) }
	if(is.null(Q_index))	{Q_index		<- 1:K }
	ns		<- length(plot_state)
	nQ 		<- length(Q_index)
	
	c_star	<- as.matrix(DP_list$policy$c[plot_state,Q_index])
	dimnames(c_star) <- list(plot_state, Q_grid)
	vc		<- array(NA, c(ns, nQ, nc+1), 
				dimnames = list(plot_state, Q_index, c(c_grid, "")) )
	
	inter_spline <- control_list$inter_spline
	# Value function interpolation: cspline or linear
	if(inter_spline){
		V1_fn <- splinefun(state[state[,"y"]==1,I_idx],V[state[,"y"]==1], method = "natural")
		V2_fn <- splinefun(state[state[,"y"]==2,I_idx],V[state[,"y"]==2], method = "natural")	
	}else{
		# NOTE: extrapolate the value function outside the grid.
		V1_fn <- approxfun(state[state[,"y"]==1,I_idx],V[state[,"y"]==1], rule = 1)
		V2_fn <- approxfun(state[state[,"y"]==2,I_idx],V[state[,"y"]==2], rule = 1)	
	}
	value_fn <- function(c,I,Q,sel_y,sel_k){
		I_next 	<- I+Q-c
		if(I_next<0){
			return(NA)
		}else{
			v		<- c(V1_fn(I_next), V2_fn(I_next))
			out 	<- uflow_fn(c, I, Q, sel_y, sel_k, omega, param) + beta*y_kernal$weight[sel_y,] %*% v
			return(out)
		}	
	}
		
	for(i in 1:ns){
		for(j in 1:nQ){
			for(k in 1:nc){
				vc[i,j,k] <- value_fn(c_grid[k], state[plot_state[i],I_idx], 
					Q_grid[Q_index[j]], sel_y = s_idx[plot_state[i],y_idx], sel_k=j)
			}
			if(is.na(c_star[i,j])){
				vc[i,j,(nc+1)] <- NA
			}else{
				vc[i,j,(nc+1)] <- value_fn(c_star[i,j], state[plot_state[i],I_idx], 
					Q_grid[Q_index[j]], sel_y = s_idx[plot_state[i],y_idx], sel_k=j)
			}
		}
	}
	
	ggtmp 		<- melt(vc)
	sel 		<- is.na(ggtmp$Var3)
	ggtmp[sel,"Var3"] <- melt(c_star)$value
	ggtmp$max_pnt <- 0
	ggtmp[sel,"max_pnt"] <- 1
	names(ggtmp)[1:3]	<- c("s_index", "Q_index", "c")
	ggtmp$I		<- state[ggtmp$s_index, I_idx]
	ggtmp$y		<- state[ggtmp$s_index, y_idx]
	ggtmp$Q		<- Q_grid[ggtmp$Q_index]
	
	if(plot){
		for(i in 1:nQ){
			quartz()
			print(ggplot(subset(ggtmp, Q_index==Q_index[i]), aes(c, value)) + geom_point(aes(col=factor(max_pnt))) + 
						geom_line() + 
						facet_wrap(I ~ y) + 
						scale_color_manual(values=c("black","red")) + 
						guides(color=FALSE) + 
						labs(title = paste("Value function curve at Q=",
											round(Q_grid[Q_index[i]],2), "\n(facet by I,y)", sep=""))
			)
		}
	}else{
		return(ggtmp)
	}
}


simulate_seq_fn <- function(init_state, TT, draw_all = FALSE, y_seq, Q_seq, k_seq, DP_list, control_list, 
							alpha, gamma, omega, nu=NULL){
	# Extract the solution to dynamic programming
	state 		<- DP_list$state
	Q_grid	 	<- DP_list$Q_grid
	y_kernal	<- DP_list$y_kernal
	V			<- DP_list$value_fn
	inter_spline <- control_list$inter_spline
	if(is.null(nu)){ nu <- rep(0, TT) }
	
	# Value function
	if(inter_spline){
		V1_fn <- splinefun(state[state[,"y"]==1,I_idx],V[state[,"y"]==1], method = "natural")
		V2_fn <- splinefun(state[state[,"y"]==2,I_idx],V[state[,"y"]==2], method = "natural")	
	}else{
		# NOTE: extrapolate the value function outside the grid.
		V1_fn <- approxfun(state[state[,"y"]==1,I_idx],V[state[,"y"]==1], rule = 1)
		V2_fn <- approxfun(state[state[,"y"]==2,I_idx],V[state[,"y"]==2], rule = 1)	
	}
	value_fn <- function(c,I,Q,sel_y,sel_k,nu){
		I_next 	<- I+Q-c
		if(I_next<0){
			return(NA)
		}else{
			v		<- c(V1_fn(I_next), V2_fn(I_next))
			out 	<- uflow_fn(c, I, Q, sel_y, sel_k, omega, param,nu) + beta*y_kernal$weight[sel_y,] %*% v
			return(out)
		}	
	}
	
	# The objects that store output
	DataState	<- matrix(NA, TT, ncol(DP_list$state))
	names(DataState) <- names(DP_list$state)
	if(draw_all){
		k_seq <- rep(NA, TT)
		e_seq <- matrix(NA, TT, S+2)
	}else{
		DataState[,y_idx] <- y_seq
	}
	c_seq		<- rep(NA, TT)
	ccp			<- matrix(NA, TT, K)
	
	# Loop over to find three chocies: basket k, expenditure allocation e, and consumption c 
	DataState[1,] <- init_state
	for(i in 1:TT){
		sel_y	<- DataState[i,2]
		tmp_sol <- Bellman_operator(DP_list, DataState=matrix(DataState[i,], nrow=1))
		ccp[i,]	<- tmp_sol$ccp
		# If the choice k and Q is to be simulated, we solve for k first, and simulate Q conditionally;
		# Otherwise, we directly solve for optimal consumption level that determines the state in the next period.
		if(draw_all){
			k_seq[i]	<- which.max(rmultinom(1, 1, ccp[i,]))
			eps			<- rnorm(S)					# Draw allocation error
			psi			<- c(exp(alpha + eps), rep(1, Sa-S))
			w			<- Allocation_fn(y=y_value[sel_y], psi=psi, gamma=c(gamma,rep(1,Sa-S)), Q=Q_grid[k_seq[i]], price=price,S=S, Sa=Sa)
			e_seq[i,]	<- w$e
		}
		tmp_sol 		<- optimize(value_fn, interval=c(0,DataState[i,1]+Q_grid[k_seq[i]]), 
						maximum=T, I=DataState[i,1], Q=Q_grid[k_seq[i]], sel_y=sel_y,sel_k=k_seq[i],nu=nu[i],
						tol = control_list$brent_tol)
		c_seq[i]		<- tmp_sol$maximum
		
		# Detemine the state of next period
		if(i<TT){
			DataState[(i+1),I_idx] <- DataState[i,I_idx] + Q_grid[k_seq[i]] - c_seq[i]
			if(draw_all){
				DataState[(i+1),y_idx] <- sample(y_kernal$nodes, 1, prob = y_kernal$weight[sel_y,])
			}
		}
	}
	return(list(data.frame(	t=1:TT, I=DataState[,I_idx], y=DataState[,y_idx], k = k_seq, e = e_seq, 
						c = c_seq),
				ccp = ccp) )
}
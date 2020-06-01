# Purchae utility function 
uP_fn <- function(e, psi, gamma,price, S){
# Negative purchase utility 
# e: S+2 vecotr
# psi, gamma, price is (S+2)-dimension vector 
	e0 		<- e[1:S]
	e_S1	<- e[(S+1):length(e)]
	u0		<- c(log(e0/(gamma[1:S]*price) + 1), log(e_S1[1]/gamma[(S+1)]), log(e_S1[2]/gamma[Sa]+qz_cons) )
	out 	<- sum(psi*gamma*u0 ) 
	return(-out)
}

uPGrad_fn <- function(e, psi, gamma,price, S){
	e0 		<- e[1:S]
	e_S1	<- e[(S+1):length(e)]
	u0		<- c(1/(e0+price*gamma[1:S]), 1/e_S1[1], 1/(e_S1[2]+qz_cons*gamma[Sa]) )
	out 	<- psi*gamma*u0 
	return(-out)	
}

# Flow utility function
uflow_fn <- function(c, I, Q, y, sel_k, omega, param, nu=0){
	lambda	<- param[1]
	tau1  	<- param[2]
	tau2	<- param[3]
	u		<- lambda*log(c+.01+nu) - tau1*(I+Q-c) -tau2*1*(Q>0)+ omega(sel_k,y)
	return(u)
}

# Value function 
value_fn <- function(c,I,Q,y,sel_k, y_grid, v1_interp, omega, param, beta, nu=0){
	I_next 	<- I+Q-c
	y_next 	<- exp(sqrt(2)*y_kernal$sigma*y_kernal$nodes + y_kernal$mu + y_kernal$rho*log(y))

	# First interpolate the value function at I_next along y grid;
	v2			<- sapply(1:length(y_grid), function(i) v1_interp[[i]](I_next) ) 
	# Second interpolate the value function at y_next 
	v2_interp	<- splinefun(x=y_grid, y=v2, method="natural")
	v			<- v2_interp(y_next)
	
	# Calculate the RHS of Bellman equation
	out 	<- uflow_fn(c, I, Q, y, sel_k, omega, param) + beta*y_kernal$weight %*% v /sqrt(pi)
	return(out)
}

# Allocation function - solve a constrained optimization
Allocation_fn <- function(y, psi, gamma, Q=NULL, price,S,Sa, inits=NULL, silent=FALSE){
	# We have Sa non-negativity constraints, one budget constraint, and one quantity constraint
	# ui <- rbind(diag(S+2), c(rep(-1, S+1), 0) )
	# inits is a list of initial values
	ui <- rbind(diag(Sa), c(rep(-1, S+1),rep(0, Sa-S-1 )) )
	ci <- c(rep(0,Sa), -y)
	if(is.null(inits)){
		tmp	  <- min(y/price, Q) * .99
		sel   <- which.min(price)
		inits <- list(c(rep(tmp/(Sa-1), Sa-1), Q - tmp), 
					  c(rep(1e-8, S), y-(S+2)*1e-8, Q - sum(price)*1e-8 - 1e-8),
					  c(rep(tmp/S, S), y- sum(tmp/S*price) -1e-8, 1e-8) )
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
	y_grid 	<- unique(state[,y_idx])
	
	# Create spline functions along the first dimension
	v1_interp <- vector("list", length=length(y_grid))
	for(i in 1:length(y_grid)){
		sel <- state[,y_idx]==y_grid[i]
		v1_interp[[i]] <- splinefun(x=state[sel,I_idx], y=V_old[sel], method="natural")
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
		for(k in 1:K){
			if(state[i,1]+Q_grid[k]<=0){
				c_star[i,k]	<- 0
				vc[i,k]		<- value_fn(c=0,I=state[i,1],Q=Q_grid[k],y=state[i,2],sel_k=k, y_grid, v1_interp, omega, param, beta)
			}else if(Q_grid[k]>state[i,2]/min(price)){
				vc[i,k]		<- (-500)
				c_star[i,k]	<- NA			
			}else{
				sol 		<- optimize(value_fn, interval=c(0,state[i,1]+Q_grid[k]), 
								maximum=T, I=state[i,1], Q=Q_grid[k], y=state[i,2],sel_k = k,
								y_grid = y_grid, v1_interp = v1_interp, omega = omega, 
								param = param, beta = beta,tol = control_list$brent_tol)
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
	return(list(policy = list(k = policy_k, c = c_star), value = V, ccp = ccp, vc=vc))
}

# Bellman operator wrapper
Bellman_operator <- function(DP_list, control_list, DataState=NULL){
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

simulate_seq_fn <- function(init_state, TT, draw_all = FALSE, y_seq, Q_seq, k_seq, DP_list, control_list, 
							alpha, alpha_a, gamma, nu=NULL){
	# Extract the solution to dynamic programming
	state 		<- DP_list$state
	Q_grid	 	<- DP_list$Q_grid
	y_kernal	<- DP_list$y_kernal
	V			<- DP_list$value_fn
	y_grid 		<- unique(state[,y_idx])
	psi_a		<- exp(alpha_a)
	if(is.null(nu)){ nu <- rep(0, TT) }
	
	# Create spline functions along the first dimension
	v1_interp <- vector("list", length=length(y_grid))
	for(i in 1:length(y_grid)){
		sel <- state[,y_idx]==y_grid[i]
		v1_interp[[i]] <- splinefun(x=state[sel,I_idx], y=V[sel], method="natural")
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
		tmp_sol <- Bellman_operator(DP_list, control_list, DataState=matrix(DataState[i,], nrow=1))
		ccp[i,]	<- tmp_sol$ccp
		if(draw_all){
			k_seq[i]	<- which.max(rmultinom(1, 1, ccp[i,]))
			eps			<- rnorm(S)					# Draw allocation error
			psi			<- c(exp(alpha + eps), psi_a)
			w			<- Allocation_fn(y=DataState[i,2], psi=psi, gamma=c(gamma,rep(1,Sa-S)), Q=Q_grid[k_seq[i]], price=price,S=S, Sa=Sa)
			e_seq[i,]	<- w$e
		}
		tmp_sol 		<- optimize(value_fn, interval=c(0,DataState[i,1]+Q_grid[k_seq[i]]), 
								maximum=T, I=DataState[i,1], Q=Q_grid[k_seq[i]], y=DataState[i,2],sel_k = k_seq[i],
								y_grid = y_grid, v1_interp = v1_interp, omega = DP_list$omega, 
								param = DP_list$param, beta = DP_list$beta,tol = control_list$brent_tol)
		c_seq[i]		<- tmp_sol$maximum
		
		# Detemine the state of next period
		if(i<TT){
			DataState[(i+1),I_idx] <- DataState[i,I_idx] + Q_grid[k_seq[i]] - c_seq[i]
			if(draw_all){
				DataState[(i+1),y_idx] <- exp(rnorm(1, y_kernal$mu + y_kernal$rho*log(DataState[i,y_idx]), y_kernal$sigma) )
			}
		}
	}
	return(list(data.frame(	t=1:TT, I=DataState[,I_idx], y=DataState[,y_idx], k = k_seq, e = e_seq, 
						c = c_seq),
				ccp = ccp) )
}

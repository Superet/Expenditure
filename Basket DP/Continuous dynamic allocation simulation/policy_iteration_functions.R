flowU	<- function(y, c, I, Inc, Rc, lnZ, param, omega_fn,Q_fn){
	lambda		<- param[1]
	lambda_o	<- param[2]
	tau			<- param[3]
	u 			<- omega_fn(y) + lambda*log(c+.01) + lambda_o*(Inc - y - exp(lnZ)) - tau*(I+Q_fn(y) - c)
	return(u)
}

spl_init <- function(state, V){
# NOTE: 1st and 4th dimension are continuous variables. (Inventory I and exogenous shock Z)
	x1		<- sort(unique(state[,1]))
	x2		<- sort(unique(state[,2]))
	x3		<- sort(unique(state[,3]))
	x4		<- sort(unique(state[,4]))
	
	V_spl_ls<- vector("list", length(x2))
	for(i in 1:length(x2)){
		V_spl_ls[[i]] <- vector("list", length(x3))
		for(j in 1:length(x3)){
			V_spl_ls[[i]][[j]] <- vector("list", length(x4))
			for(k in 1:length(x4)){
				sel	<- state[,2] == x2[i] & state[,3] == x3[j] & state[,4] == x4[k]
				V_spl_ls[[i]][[j]][[k]] <- splinefun(x1, V[sel], method="natural")
			}
		}
	}
	V_spl <- function(s1){
		x2_idx 	<- which(x2 == s1[2])
		x3_idx	<- which(x3 == s1[3])
		y4		<- rep(NA, length(x4))
		for(i in 1:length(x4)){
			y4[i]	<- V_spl_ls[[x2_idx]][[x3_idx]][[i]](s1[1])
		}
		spl4 	<- splinefun(x4, y4, method="natural")
		spl4(s1[4])
	}
	return(V_spl)
}

valuefn	<- function(y, c, I, Inc_index, Rc, lnZ, V_spl, Inc_weight, DP_list){
	param	<- DP_list$param
	state	<- DP_list$state
	beta	<- DP_list$beta
	Q_fn	<- DP_list$Q_fn
	omega_fn<- DP_list$omega
	lnZ_kernel	<- DP_list$lnZ_kernel
	lnZ_nodes	<- lnZ_kernel$nodes
	lnZ_weight	<- lnZ_kernel$weight
	Inc_kernel	<- DP_list$Inc_kernel
	Inc_nodes	<- Inc_kernel$nodes
	mu			<- param[4]
	sigma		<- exp(param[5])
	
	# Interpolate the value function evaluated at the next state
	I_next	<- I + Q_fn(y) - c
	lnZ_next<- sqrt(2)*sigma*lnZ_nodes + mu 
	v1		<- matrix(NA, length(Inc_nodes), length(lnZ_next))
	for(i in 1:length(Inc_nodes)){
		for(j in 1:length(lnZ_next)){
			s_next	<- c(I_next, i, Rc, lnZ_next[j])
			v1[i,j] <- V_spl(s_next)
		}
	}
	ev 	<- Inc_weight %*% v1 %*% lnZ_weight / sqrt(pi)
	V 	<- flowU(y, c, I, Inc_nodes[Inc_index], Rc, lnZ, param, omega_fn,Q_fn) + beta * ev
	return(V)
}

policy_impv	<- function(V_old, DP_list, control_list){
	state 	<- DP_list$state
	Q_fn	<- DP_list$Q_fn
	k		<- DP_list$k
	Inc_nodes	<- DP_list$Inc_kernel$nodes
	Inc_weight 	<- DP_list$Inc_kernel$weight
	n_Inc	<- ncol(Inc_weight)
	ns		<- nrow(state)
	nu		<- 2
	policy	<- matrix(NA, ns, nu); colnames(policy) <- c("y","c")
	V		<- rep(NA, nrow(state))
	multimin_max <- control_list$multimin_max
	Imax	<- max(state[,1])
	
	# Construct spline for the value function
	V_spl 	<- spl_init(state, V_old)
	
	# Solve the maximization problem at each state
	eps 	<- 1e-4
	for(i in 1:ns){
		# value function wrapper 
		my_value <- function(u){
			y <- u[1]
			c <- u[2]
			idx <- state[i,3] * n_Inc + state[i,2]
			-valuefn(y, c, I=state[i,1], Inc_index=state[i,2], Rc=state[i,3], lnZ=state[i,4], V_spl, Inc_weight=Inc_weight[idx,], DP_list)
		}
		
		ui 		<- rbind(diag(nu), -1*diag(nu))
		ui[4,1] <- k
		ui 		<- rbind(ui, c(-k, 1))
		Inc 	<- Inc_nodes[state[i,2]]
		ci 		<- c(rep(0, nu), -Inc + exp(state[i,4]), - state[i,1], state[i,1] - Imax)
		inits 	<- c(Inc - exp(state[i,4]) - eps, (state[i,1]+k*(Inc - exp(state[i,4]) - eps))/2)
		tmp_sol	<- constrOptim(inits, f=my_value, ui=ui, ci=ci, method="Nelder-Mead",outer.iterations = multimin_max)		
		tmp		<- tmp_sol$par 
		
		# Check boundary
		sel <- which(ui %*% tmp - ci < eps)
		if(length(sel)>0){
			if(any(sel <= nu)){
				sel1 	<- sel[sel<=nu]
				tmp[sel1] <- 0
			}else{
				if((nu+1) %in% sel){
					tmp[1]	 <- -ci[(nu+1)]
				}
				if((nu+2) %in% sel ){
					tmp[2] <- state[i,1] + k*tmp[1]
				}
			}
		}
		policy[i,] 	<- tmp
		V[i]		<- -my_value(tmp)
	}
	
	return(list(policy=policy, V=V))
}

policy_eval <- function(policy, V_old, V0, DP_list, control_list, M=NULL){
	state		<- DP_list$state
	Inc_weight	<- DP_list$Inc_kernel$weight

	if(is.null(M)){
		M 		<- control_list$inner_max
	}
	tol		<- control_list$inner_tol
	
	V_new <- V0
	dif <- max(abs(V_new - V_old))
	if(dif > tol){
		k 		<- 0 
		V_old 	<- V0
		while(k <= M & dif > tol){
			# Construct spline for the value function
			V_spl	<- spl_init(state, V_old)

			# Compute new value function 
			for(i in 1:nrow(state)){
				V_new[i] <- valuefn(policy[i,1], policy[i,2], I=state[i,1], Inc_index=state[i,2], Rc = state[i,3], lnZ=state[i,4], 
									V_spl, Inc_weight=Inc_weight[(n_Inc*state[i,3] + state[i,2]),], DP_list)
			}
			
			# Check norm distance
			dif 	<- max(abs(V_new - V_old))
			V_old 	<- V_new
			k 		<- k + 1
			cat("Inner iteration", k, ", norm=", dif,"\n")
		}
	}
	return(V_new)
}

policy_iter <- function(DP_list, control_list, print.level=0){
	tol		<- control_list$tol
	max_iter<- control_list$max_iter
	V_old	<- DP_list$value_fn
	
	# Compute the policy under the initial value functions
	tmp			<- policy_impv(V_old, DP_list)
	policy_old	<- tmp$policy
	V_old		<- tmp$V
	
	dif 	<- tol + 1
	k 		<- 0
	while(k <= max_iter & dif >= tol){
		# Policy improvement 
		tmp			<- policy_impv(V_old, DP_list)
		policy_new	<- tmp$policy
		V0			<- tmp$V
		
		# Policy evaluation 
		V_new 		<- policy_eval(policy_new, V_old, V0=V0, DP_list, control_list)
		dif 		<- max(abs(policy_new - policy_old))
		k			<- k + 1
		V_old		<- V_new
		policy_old	<- policy_new
		
		if(print.level>0){
			cat("Policy iteration", k, ", norm =", dif, "\n")
		}
	}
	
	DP_list$value_fn 	<- policy_eval(policy_new, V_old, V_new, DP_list, control_list, M=control_list$max_iter)
	DP_list$policy		<- policy_new
	DP_list$Iteration	<- k
	if(k > max_iter & dif>= tol){
		cat("Policy iteration does not converge within the maxium iteration.\n")
		DP_list$status <- 2
	}else{
		DP_list$status <- 0
	}
	
	return(DP_list)
}

simulate_seq1_fn<- function(init_state, TT, DP_list, control_list, Inc_draw = NULL, lnZ_draw=NULL){
	state 		<- DP_list$state
	policy		<- DP_list$policy
	k			<- DP_list$k
	Inc_weight	<- DP_list$Inc_kernel$weight
	Inc_nodes	<- DP_list$Inc_kernel$nodes
	Inc_grid	<- sort(unique(state[,2]))
	n_Inc		<- length(Inc_grid)
	Q_fn		<- DP_list$Q_fn 
	Imax		<- max(state[,1])
	mu			<- DP_list$param[4]
	sigma		<- exp(DP_list$param[5])
	
	# Construct policy spline
	y_spl 		<- spl_init(state, policy[,1])
	c_spl 		<- spl_init(state, policy[,2])
	
	DataState 	<- matrix(NA, TT, ncol(state), dimnames=list(t=1:TT, colnames(state)))
	a			<- matrix(NA, TT, ncol(policy), dimnames=list(t=1:TT, c("y","c")))
	astar		<- matrix(NA, TT, ncol(policy), dimnames=list(t=1:TT, c("y","c")))
	if(is.null(Inc_draw)){
		draw_inc <- TRUE
	}else{
		draw_inc <- FALSE
		DataState[,2]	<- Inc_draw
	}
	DataState[,3] <- init_state[3] 
	if(is.null(lnZ_draw)){
		DataState[,4]	<- rnorm(nrow(DataState), mu, sigma)
	}else{
		DataState[,4]	<- lnZ_draw*sigma + mu 
	}
	DataState[1,] <- init_state
	
	for(i in 1:TT){
		# Interpolate policy at state
		tmpy 	<- y_spl(DataState[i,])
		tmpc	<- c_spl(DataState[i,])
		Inc 	<- Inc_nodes[DataState[i,2]]
		upper_y	<- max(Inc - exp(DataState[i,4]), 0)
		a[i,1]	<- max(0, min(tmpy, upper_y))
		q		<- Q_fn(a[i,1])
		lower_c	<- max(0, DataState[i,1] + q - Imax)
		a[i,2]	<- max(lower_c, min(tmpc, DataState[i,1] + q))
		astar[i,]<- c(tmpy, tmpc)
		
		# State transition
		if(i <= TT- 1){
			DataState[(i+1),1]	<- DataState[i,1] + Q_fn(a[i,1]) - a[i,2]
			if(draw_inc){
				DataState[(i+1),2] 	<- sample(Inc_grid, 1, prob=Inc_weight[(n_Inc*DataState[i,3] + DataState[i,2]),])			
			}
		}
	}
	return(data.frame(t=1:TT, DataState, a))
}

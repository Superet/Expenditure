# Check the concavity of value function when solving for c; 
# Plot the value function a long a grid of c.
plot_cvalue_fn <- function(DP_list, c_grid, plot_state=NULL, choice_k=NULL, control_list, plot = TRUE){
	state	<- DP_list$state
	V		<- DP_list$value_fn
	Q_grid	<- DP_list$Q_grid
	TR_grid	<- DP_list$TR_grid
	omega_fn<- DP_list$omega
	K		<- DP_list$K
	param	<- DP_list$param
	beta	<- DP_list$beta
	Inc_weight	<- DP_list$Inc_weight
	y_kernal	<- DP_list$y_kernal
	price	<- DP_list$price
	omega_mat 	<- sapply(1:K, function(x) omega_fn(x, state[,2]))
	
	sol		<- Solve_DP_fnC(state, V, omega_mat, K, param, beta, Q_grid, TR_grid, 
				Inc_weight, y_kernal, price, control_list)
	
	nc		<- length(c_grid)
	if(is.null(plot_state)) {plot_state <- 1:nrow(state) }
	if(is.null(choice_k))	{choice_k		<- 1:K }
	ns		<- length(plot_state)
	nQ 		<- length(choice_k)
	
	c_star	<- as.matrix(DP_list$policy$c[plot_state,choice_k])
	dimnames(c_star) <- list(plot_state, Q_grid[choice_k])
	vc		<- array(NA, c(ns, nQ, nc+1), 
				dimnames = list(plot_state, choice_k, c(c_grid, "")) )
	
	inter_spline <- control_list$inter_spline
	x1		<- sort(unique(state[,1]))
	x2		<- sort(unique(state[,2]))
	x3		<- sort(unique(state[,3]))
	x4		<- sort(unique(state[,4]))
	N1		<- length(x1)
	N2		<- length(x2)
	N3		<- length(x3)
	N4		<- length(x4)
	
	# Value function interpolation: cspline or linear
	spl_list <- lapply(1:(N2*N3*N4), function(i)
					 splinefun(x = x1, y = V[seq(N1*i-N1+1, N1*i)], method="natural") )
	value_fn <- function(c,I,lny,Inc,Rc,Q,k){
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
					Q_grid[choice_k[j]], sel_y = s_idx[plot_state[i],y_idx], sel_k=j)
			}
			if(is.na(c_star[i,j])){
				vc[i,j,(nc+1)] <- NA
			}else{
				vc[i,j,(nc+1)] <- value_fn(c_star[i,j], state[plot_state[i],I_idx], 
					Q_grid[choice_k[j]], sel_y = s_idx[plot_state[i],y_idx], sel_k=j)
			}
		}
	}
	
	ggtmp 		<- melt(vc)
	sel 		<- is.na(ggtmp$Var3)
	ggtmp[sel,"Var3"] <- melt(c_star)$value
	ggtmp$max_pnt <- 0
	ggtmp[sel,"max_pnt"] <- 1
	names(ggtmp)[1:3]	<- c("s_index", "choice_k", "c")
	ggtmp$I		<- state[ggtmp$s_index, I_idx]
	ggtmp$y		<- state[ggtmp$s_index, y_idx]
	ggtmp$Q		<- Q_grid[ggtmp$choice_k]
	
	if(plot){
		for(i in 1:nQ){
			quartz()
			print(ggplot(subset(ggtmp, choice_k==choice_k[i]), aes(c, value)) + geom_point(aes(col=factor(max_pnt))) + 
						geom_line() + 
						facet_wrap(I ~ y) + 
						scale_color_manual(values=c("black","red")) + 
						guides(color=FALSE) + 
						labs(title = paste("Value function curve at Q=",
											round(Q_grid[choice_k[i]],2), "\n(facet by I,y)", sep=""))
			)
		}
	}else{
		return(ggtmp)
	}
}

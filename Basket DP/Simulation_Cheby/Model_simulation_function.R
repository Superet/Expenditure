choice_sequence <- function(I0,y_seq,DP_list){
	m		<- DP_list$consumption["m"]
	n		<- DP_list$consumption["n"]
	t		<- length(y_seq)
	state 	<- cbind(rep(NA,t),y_seq)
	prob	<- matrix(NA,t,K)
	choice	<- rep(NA,t)
	Ec		<- rep(NA,t)
	Q		<- rep(NA,t)
	cons	<- rep(NA,t)
	
	state[1,1] <- I0
	for(i in 1:t){
		opt			<- choice_probability_once(state[i,],DP_list)		
		prob[i,]	<- opt$probability
		choice[i]	<- opt$choice
		Ec[i]		<- opt$consumption
		if(i<t){
			Q[i] 	<- rnorm(1,mean=mu_Q_fn(choice[i],y_seq[i]),sd=sigma_Q[choice[i]])
			cons[i]	<- min(max(m+n*(state[i,I.idx]+ Q[i]),0),state[i,I.idx] + Q[i])
			state[(i+1),I.idx] 	<- state[i,I.idx] + Q[i] - cons[i]
		}
	}
	return(list(state=state,choice=choice,probability=prob,Expected.consumption=Ec,Q=Q,consumption=cons))
}


choice_probability_once <- function(state,DP_list){
# state: a single D-dimentional state
	param_list <- DP_list$param_list
	n.nodes <- nrow(GH_init$nodes)
	m		<- DP_list$consumption["m"]
	n		<- DP_list$consumption["n"]
	x_prim	<- array(NA,c(K,n.nodes,state_dim))
	Ec_star	<- rep(NA,K)
	EI		<- rep(NA,K)
	flow_u	<- rep(NA,K)
	v		<- rep(NA,K)
	for(k in 1:K){
		sigma 	<- sigma_Q[k]
		a0 		<- -(m/n + mu_Q_fn(k,state[y.idx]))/sigma
		b0 		<- -(m/(n-1) + mu_Q_fn(k,state[y.idx]))/sigma
		mu 		<- state[I.idx] + mu_Q_fn(k,state[y.idx])
		if(pnorm(a0)<set_list$Phi.bound){
			EQ_lower	<- 0 
			EQ2_lower	<- 0
		}else{
			rho 		<- dnorm(a0)/pnorm(a0)
			EQ_lower	<- mu - sigma*rho
			EQ2_lower	<- mu^2 - 2*mu*sigma*rho + sigma^2*(1+a0*rho)
		}

		if(pnorm(-b0)<set_list$Phi.bound){
			EQ_upper	<- 0
			EQ2_upper	<- 0
		}else{
			rho 		<- dnorm(b0)/pnorm(-b0)
			EQ_upper	<- mu + sigma*rho
			EQ2_upper	<- mu^2 + 2*mu*sigma*rho +sigma^2*(1+b0*rho)
		}

		p.mid <- pnorm(b0) - pnorm(a0)
		if(p.mid<set_list$Phi.bound){
			EQ_mid		<- 0 
			EQ2_mid		<- 0
			Ec_mid		<- 0
			Ec2_mid		<- 0
			EI_mid		<- 0
			EI2_mid		<- 0 
		}else{
			EQ_mid		<- mu - sigma*(dnorm(b0)-dnorm(a0))/p.mid
			EQ2_mid		<- mu^2 - 2*mu*sigma*(dnorm(b0) - dnorm(a0))/p.mid + sigma^2*(1 - (b0*dnorm(b0) - a0*dnorm(a0))/p.mid)
			Ec_mid		<- m*p.mid + n*EQ_mid
			Ec2_mid		<- m^2*p.mid + 2*m*n*EQ_mid + n^2*EQ2_mid
			EI_mid		<- -m*p.mid - (n-1)*EQ_mid
			EI2_mid		<- m^2*p.mid + 2*m*(n-1)*EQ_mid + (n-1)^2*EQ2_mid
		}

		# Expected consumption, inventory and flow utility
		Ec 			<- EQ_upper + Ec_mid
		Ec2			<- EQ2_upper + Ec2_mid
		EI_prim		<- EQ_lower + EI_mid
		EI2_prim	<- EQ2_lower + EI2_mid
		flow_u[k]	<- with(param_list, lambda1*Ec + lambda2*Ec2 - tau1*EI_prim - tau2*EI2_prim) + omega_fn(k,state[y.idx])
		
		# The nodes at which to compute integral of value function
		I_prim		<- cbind(-m-(n-1)*(mu + sqrt(2)*sigma*GH_init$nodes[,I.idx]), 
							 mu + sqrt(2)*sigma*GH_init$nodes[,I.idx]	)
		I_prim		<- apply(I_prim,1,function(x) min(max(x[1],0),x[2]) )
		y_prim		<- state[y.idx] + sqrt(2)*param_list$sigma_y * GH_init$nodes[,y.idx]
		x_prim[k,,I.idx] <- sapply(I_prim, function(x) min(max(x,cheby_init$lower[I.idx]),cheby_init$upper[I.idx]) )
		x_prim[k,,y.idx] <- sapply(y_prim, function(x) min(max(x,cheby_init$lower[y.idx]),cheby_init$upper[y.idx]) )
		
		# Re-normalize nodes to compute Chebychev approximation of expected value funtion
		x_prim_adj		<- 2*(x_prim[k,,] - rep(1,n.nodes)%*%t(cheby_init$lower))/
							(rep(1,n.nodes) %*% t(cheby_init$upper - cheby_init$lower)) - 1
		T 				<- apply(x_prim_adj,1,cheb_polynm_multD,cheby_init$Polynomial.degree)
		ET				<- apply(T,1,function(x) x %*% GH_init$weights)/pi
		Ew				<- ET %*% DP_list$Cheb.theta
		v[k]			<- flow_u[k] + param_list$beta * Ew	
		Ec_star[k]		<- Ec
		EI[k]			<- EI_prim
	}
	max_v 	<- max(v)
	v		<- v - max_v
	p 		<- exp(v)/sum(exp(v))
	return(list(state=state,expected.next.state=c(EI[which.max(v)],state[2]),choice=which.max(v),probability = p,consumption=Ec_star[which.max(v)]))
}

# Dynamic inventory model version 5
# Discrete basket choice and continous consumption choice

# MODEL
# v(I,y) = max_{k,c} U(c,I,Q) + beta*Ev(I',y'|I,y,c,Q)
# where Q ~ f(Q|k,y)
# U(c,I,Q) = lambda*log(c) - tau1*(I+Q-c) -tau2*(I+Q-c)^2

library(ggplot2)
library(reshape2)
library(scatterplot3d)
library(maxLik)
library(Rcpp)
library(RcppGSL)

# setwd("//tsclient/Resear1/Store switching/Exercise/Basket DP/Simulation_v1_grid")

setwd("~/Documents/Research/Store switching/Exercise/Basket DP/Simulation_v1_grid")
source("5_utility_functions.R")
sourceCpp("5C_Bellman_operator.cpp")

# Set parameters
control_list <- list(	value_max_iter 	= 500,
						tol				= 1e-8,
						display_freq	= 10,
						brent_tol		= 0.00001, 
						brent_max_iter 	= 200,
						inter_spline	= TRUE
				)
beta <- .9
param_list	<- list(lambda 	= 7,
					tau1 	= 0.5,
					tau3	= 2
					)
param <- unlist(param_list)

K		<- 3
step	<- .1					
I_grid	<- seq(0,5,by=step)
y_grid	<- c(1,2)
c_grid	<- seq(0,8,by=step)
state	<- cbind(I=rep(I_grid,length(y_grid)),y=rep(y_grid,each=length(I_grid)))	# state matrix
s_idx	<- cbind(I=rep(1:length(I_grid),length(y_grid)), y=rep(1:length(y_grid),each=length(I_grid)))
ns		<- nrow(state)
ny		<- length(y_grid)
I_idx	<- 1
y_idx	<- 2

y_kernal<- list(nodes = c(1,2), weight = matrix(c(.7,.3,.3,.7),2,2))
Q_grid	<- c(0,.5,1,2,2.5,3)/2
nQ		<- length(Q_grid)				
Q_weight<- array(NA, c(K, ny, nQ))
Q_weight[1,,] <- matrix(c(	1,	0,	0,	0,	0,	0,
							.7,	.2,	.1,	0,	0,	0),ny,nQ,byrow=T,)
Q_weight[2,,] <- matrix(c(	.5,	.2,	.1,	.1,	.05,.05,
							.05,.05,.2,	.3,	.3,	.1),ny,nQ,byrow=T)
Q_weight[3,,] <- matrix(c(	.1,	.2,	.5,	.1,	.05, .05,
							0.05,.05,.1,.1,	.3,	.4),ny,nQ,byrow=T)		
Q_kernal <- list(nodes = Q_grid, weight = list(Q_weight[,1,], Q_weight[,2,] ))

tmp		<- floor(quantile(1:ns, c(1:K)/K))
tmp		<- c(tmp[1], diff(tmp))
policy_k<- rev(unlist(lapply(1:K,function(i) rep(i,tmp[i])) ))

DP_list	<- list(state 		= as.matrix(state), K = K,
				value_fn 	= 0.5*log(state[,I_idx]+1) ,
				policy	 	= list(k = policy_k,c = state[,I_idx] %*% t(rep(1,nQ))),
				param		= param,
				param_list	= param_list,
				beta		= beta,
				state_idx 	= s_idx,
				status 		= 1, 
				Q_kernal	= Q_kernal,
				y_kernal	= y_kernal,
				Iteration 	= 0 )

# Value function iteration with Brent's method in R
system.time(sol <- value_iteration_fnC(DP_list, Bellman_operator = solve_c_brent_fn, print_level=1, control_list))
# DP_list <- sol

# Value function iteration with Brent's method written in cpp
system.time(sol1 <- value_iteration_fnC(DP_list, Bellman_operator = Bellman_operatorC, print_level=1, control_list))
DP_list <- sol1

# Compare C function and R function
plot(sol$value_fn, sol1$value_fn)
abline(0, 1, lty = 2)

tmp1 <- melt(sol$policy$c)
names(tmp1) <- c("state", "Q", "Rfunction")
tmp2 <- melt(sol1$policy$c)
names(tmp2) <- c("state", "Q", "Cppfunction")
ggtmp <- data.frame(tmp1, Cppfunction = tmp2[,3])
ggplot(ggtmp, aes(Rfunction, Cppfunction)) + geom_point() + 
		geom_abline(xintercept = 0, slope = 1, linetype = 2) + 
		facet_wrap( ~ Q)

##############################
# Plot the dynamic solutions #
##############################
# Plot value function 
scatterplot3d(DP_list$state[,I_idx], DP_list$state[,y_idx], DP_list$value_fn, xlab = "I", ylab="y", zlab="value", 
			main="Value funciton")
quartz()
sel <- state[,y_idx]==2
plot(state[sel,I_idx], DP_list$value_fn[sel])

# Plot the probability of policy function k
V		<- DP_list$value_fn
V1		<- V[DP_list$state[,y_idx]==1]
V2 		<- V[DP_list$state[,y_idx]==2]
tmp 	<- CCP_fnC(state, I_grid, V1, I_grid, V2, K, DP_list$param, beta, 
					DP_list$Q_kernal, DP_list$y_kernal, makedraw = FALSE, control_list = control_list)
tmp1 	<- data.frame(state, k1=tmp$ccp[,1], k2=tmp$ccp[,2], k3=tmp$ccp[,3])
ggtmp <- melt(tmp1, id.var=c("I","y"))
names(ggtmp)[3] <- "policy"
quartz()
ggplot(ggtmp, aes(I, value)) + geom_point() + geom_line() + 
		facet_grid(y ~ policy, labeller = "label_both") + 
		labs(y="Probability", title = "Policy function k" )

# Policy function of c
tmp <- data.frame(I=state[,I_idx], y=state[,y_idx], DP_list$policy$c)
names(tmp) <- c("I","y",Q_kernal$nodes)
ggtmp <- melt(tmp,id=c("I","y"))
names(ggtmp) <- c("I","y","Q","c")
ggtmp$Q	<- as.numeric(as.character(ggtmp$Q))
quartz()
ggplot(ggtmp, aes(I,c)) + geom_point() + geom_line() + 
		facet_grid(Q~y, labeller = "label_both") + 
		geom_abline(aes(intercept = Q, slope = 1), linetype=2)

# Check the concavity of value function along c if using value function iteration. 
plot_state <- 88:93
Q_index <- 1
plot_cvalue_fn(DP_list, c_grid, plot_state, Q_index, plot = TRUE, control_list = control_list)

# # Check the concavity of flow utility
# I <- 1
# Q <- 1
# v1 <- with(param_list, lambda*log(c_grid+ 1))
# v2 <- with(param_list, -tau1*(I+Q-c_grid)-tau2*(I+Q-c_grid)^2)
# par(mfrow= c(3,1))
# plot(c_grid,v1)
# plot(c_grid,v2)
# plot(c_grid,v1+ v2)

##############################
# Identification sensitivity # 
##############################
# CCP at baseline parameters
V		<- DP_list$value_fn
V1		<- V[DP_list$state[,y_idx]==1]
V2 		<- V[DP_list$state[,y_idx]==2]
tmp 	<- CCP_fnC(state, I_grid, V1, I_grid, V2, K, DP_list$param, beta, 
					DP_list$Q_kernal, DP_list$y_kernal, makedraw = FALSE, control_list = control_list)
ccp_base<- data.frame(state, k1=tmp$ccp[,1], k2=tmp$ccp[,2], k3=tmp$ccp[,3])

sel_param 	<- "lambda"		# The parameter name that is changed
delta 		<- -8			# The change value
param_list1	<- param_list
param1		<- param
param_list1[[sel_param]] <- param_list[[sel_param]] + delta
param1[sel_param] <- param[sel_param] + delta
DP_list1	<- list(state 		= as.matrix(state), K = K,
				value_fn 	= 0.5*log(state[,I_idx]+1) ,
				policy	 	= list(k = policy_k,c = state[,I_idx] %*% t(rep(1,nQ))),
				param		= param1,
				param_list	= param_list1,
				beta		= beta,
				state_idx 	= s_idx,
				status 		= 1, 
				Q_kernal	= Q_kernal,
				y_kernal	= y_kernal,
				Iteration 	= 0 )
DP_list1 	<- value_iteration_fnC(DP_list1, Bellman_operator = Bellman_operatorC, print_level=1, control_list)
V			<- DP_list1$value_fn
V1			<- V[DP_list1$state[,y_idx]==1]
V2 			<- V[DP_list1$state[,y_idx]==2]
tmp 		<- CCP_fnC(state, I_grid, V1, I_grid, V2, K, DP_list1$param, beta, 
					DP_list1$Q_kernal, DP_list1$y_kernal, makedraw = FALSE, control_list = control_list)
ccp1		<- data.frame(state, k1=tmp$ccp[,1], k2=tmp$ccp[,2], k3=tmp$ccp[,3])

ggtmp 		<- rbind(data.frame(ccp_base, Parameter = "Baseline"), data.frame(ccp1, Parameter = "Experiment"))
ggtmp		<- melt(ggtmp, id.var=c("Parameter","I","y"))
names(ggtmp)<- c("Parameter","I","y","policy","Probability")
quartz()
ggplot(ggtmp, aes(I, Probability, col=Parameter)) + geom_point() + geom_line() + 
		facet_grid(y ~ policy, labeller = "label_both") + 
		labs( title = paste("CCP over states\n Baseline: ", paste(names(param), param, collapse=",", sep="="), 
					"\n Experiment: ", paste(names(param1), param1, collapse=",", sep="="), sep="") )
					
# Compare the consumption difference 
tmp <- data.frame(I=state[,I_idx], y=state[,y_idx], DP_list$policy$c)
names(tmp) <- c("I","y",Q_kernal$nodes)
ggtmp <- melt(tmp,id=c("I","y"))
names(ggtmp) <- c("I","y","Q","c")
ggtmp$Q	<- as.numeric(as.character(ggtmp$Q))
ggtmp$Parameter <- "Baseline"

tmp1 <- data.frame(I=state[,I_idx], y=state[,y_idx], DP_list1$policy$c)
names(tmp1) <- c("I","y",Q_kernal$nodes)
ggtmp1 <- melt(tmp1,id=c("I","y"))
names(ggtmp1) <- c("I","y","Q","c")
ggtmp1$Q	<- as.numeric(as.character(ggtmp1$Q))
ggtmp1$Parameter <- "Experiment"
ggtmp <- rbind(ggtmp, ggtmp1)
quartz()
ggplot(ggtmp, aes(I,c, col=Parameter)) + geom_point() + geom_line() + 
		facet_grid(Q~y, labeller = "label_both") + 
		geom_abline(aes(intercept = Q, slope = 1), linetype=2) + 
		labs( title = paste("Policy function c over states\n Baseline: ", paste(names(param), param, collapse=",", sep="="), 
					"\n Experiment: ", paste(names(param1), param1, collapse=",", sep="="), sep="") )

#####################
# Simulate the data # 
#####################
simulate_sequence_fn <- function(init_state, TT, DP_list){
	state	<- matrix(NA, TT, ncol(DP_list$state))
	colnames(state) <- colnames(DP_list$state)
	Q_kernal <- DP_list$Q_kernal
	V		<- DP_list$value_fn
	V1		<- V[DP_list$state[,y_idx]==1]
	V2 		<- V[DP_list$state[,y_idx]==2]
	Q	 	<- rep(NA, TT)
	c		<- rep(NA, TT)
	k		<- rep(NA, TT)
	ccp		<- matrix(NA, TT, K)
	
	state[1,] <- init_state
	for(i in 1:TT){	
		sol 	<- CCP_fnC(matrix(state[i,],1,2), I_grid, V1, I_grid, V2, K, DP_list$param, beta, 
					DP_list$Q_kernal, DP_list$y_kernal, makedraw = TRUE, control_list = control_list)
		ccp[i,]	<- sol$ccp
		k[i]	<- sol$draw
		Q[i]	<- sample(Q_kernal$nodes, 1, prob = Q_kernal$weight[[state[i,y_idx]]][k[i],])
		c[i]	<- sol$policy$c[which(Q_kernal$nodes==Q[i])]
		
		if(i<TT){
			state[(i+1),I_idx] <- state[i,I_idx] + Q[i] - c[i]			
			state[(i+1),y_idx] <- sample(y_kernal$nodes, 1, prob = y_kernal$weight[which(y_kernal$nodes==state[i,y_idx]),])
		}	
	}
	return(cbind(t = 1:TT, state, Q = Q, c=c, k=k, ccp=ccp))
}

set.seed(99)
TT <- 400
# my_data <- simulate_seq_fn(init_state=state[5,], T, DP_list, I_grid, Q_kernal, y_kernal)    
my_data <- simulate_sequence_fn(init_state=state[80,], TT, DP_list)    

ggtmp <- melt(data.frame(my_data[,1:6]), id="t") 
quartz()
ggplot(ggtmp,aes(t, value)) + geom_point() + geom_line() + 
		facet_grid(variable ~ .)

##############
# Estimation #
##############
theta_init <- param
# sel 	<- theta_init >0	
# theta_init[sel] <- log(theta_init[sel])
# theta_init[!sel] <- log(-theta_init[!sel])
tmp		<- floor(quantile(1:ns, c(1:K)/K))
tmp		<- c(tmp[1], diff(tmp))
policy_k<- rev(unlist(lapply(1:K,function(i) rep(i,tmp[i])) ))
DP_init	<- list(param 		= rep(0, length(theta_init)),
				state 		= as.matrix(state), 
				K 			= K,
				value_fn 	= 0.5*log(state[,I_idx]+1) ,
				policy	 	= list(k = policy_k,c = state[,I_idx] %*% t(rep(1,nQ))),
				state_idx 	= s_idx,
				beta		= beta,
				status 		= 1, 
				Q_kernal	= Q_kernal, 
				y_kernal 	= y_kernal, 
				Iteration	= 0)
system.time(tmp <- ll_fnC(theta_init, DataState = my_data[, c("I","y")], choice_seq=my_data[,"k"], 
						DP_init, Bellman_operator = Bellman_operatorC, control_list = control_list))

pct		<- proc.time()
optsol 	<- maxLik(ll_fnC, DataState = my_data[, c("I","y")], choice_seq=my_data[,"k"], 
						DP_init = DP_init, Bellman_operator = Bellman_operatorC,control_list = control_list,
				start = theta_init,method ="BHHH", print.level=4, fixed = c(F,F,T))
proc.time() - pct



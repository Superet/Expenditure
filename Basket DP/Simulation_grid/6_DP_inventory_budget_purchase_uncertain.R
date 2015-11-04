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
library(glmmML)

# setwd("//tsclient/Resear1/Store switching/Exercise/Basket DP/Simulation_v1_grid")

setwd("~/Documents/Research/Store switching/Exercise/Basket DP/Simulation_v1_grid")
source("6_utility_functions.R")
sourceCpp("6C_Bellman_operator.cpp")

# Set parameters
control_list <- list(	value_max_iter 	= 500,
						tol				= 1e-8,
						display_freq	= 20,
						brent_tol		= 0.00001, 
						brent_max_iter 	= 300,
						inter_spline	= TRUE
				)
beta <- .9
param_list	<- list(lambda = 3,
					tau1 	= 0.4,
					tau3	= 4
					)
param <- unlist(param_list)

K		<- 3
step	<- .5				
I_grid	<- seq(0,10,by=step)
y_grid	<- c(1,2)
c_grid	<- seq(0,8,by=step)
state	<- cbind(I=rep(I_grid,length(y_grid)),y=rep(y_grid,each=length(I_grid)))	# state matrix
s_idx	<- cbind(I=rep(1:length(I_grid),length(y_grid)), y=rep(1:length(y_grid),each=length(I_grid)))
ns		<- nrow(state)
ny		<- length(y_grid)
I_idx	<- 1
y_idx	<- 2
y_kernal<- list(nodes = c(1,2), weight = matrix(c(.8,.2,.2,.8),2,2))

# Gaussian-Hermite quadrature initialization;
nQ <- gh_num_nodes <- 6
Q_weight<- ghq(gh_num_nodes, modified = F)$weights
Q_nodes <- ghq(gh_num_nodes, modified = F)$zeros
Qmu 	<- matrix(c(-50, -.6, -0.3,
					-50, 0, .3),2,K, byrow=T, dimnames = list(y_grid, 1:K))		# The location parameters 
Qsigma 	<- matrix(c(0, .1, .2,
					0, .2, .1),2,K, byrow=T, dimnames = list(y_grid, 1:K))
Q_kernal <- list(nodes = Q_nodes, weight = Q_weight, mu = Qmu , sigma = Qsigma)

# Plot the draws of Q
nsim <- 200
ggtmp <- data.frame(NULL)
for(i in 1:ny){
	for(j in 1:K){
		ggtmp <- rbind(ggtmp, data.frame(y=i, k=j, mu=Qmu[i,j], sigma = Qsigma[i,j], 
										Q = rlnorm(nsim, Qmu[i,j], Qsigma[i,j]) ))
	}
}
ggplot(ggtmp, aes(Q)) + geom_histogram(aes(y=..density..)) + facet_grid(y ~ k)

tmp		<- floor(quantile(1:ns, c(1:K)/K))
tmp		<- c(tmp[1], diff(tmp))
policy_k<- rev(unlist(lapply(1:K,function(i) rep(i,tmp[i])) ))
tmp		<- as.numeric(sapply(1:ns, function(i) singleu_fn(c=0,I=state[i,1],Q=0,param_list)/(1-beta)))

DP_list	<- list(state 		= as.matrix(state), K = K,
# 				value_fn 	= 0.5*log(state[,I_idx]+1) ,
				value_fn	= tmp,
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
# system.time(sol <- value_iteration_fnC(DP_list, Bellman_operator = solve_c_brent_fn, print_level=1, control_list))
# DP_list <- sol

# Value function iteration with Brent's method written in cpp
system.time(sol1 <- value_iteration_fnC(DP_list, Bellman_operator = Bellman_operatorC, print_level=1, control_list))
DP_list <- sol1

# Compare C function and R function
# plot(sol$value_fn, sol1$value_fn)
# abline(0, 1, lty = 2)

##############################
# Plot the dynamic solutions #
##############################
# Plot value function and policy function
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
tmp 				<- data.frame(I=state[,I_idx], y=state[,y_idx], DP_list$policy$c)
names(tmp) 			<- c("I","y",1:ncol(DP_list$policy$c))
ggtmp 				<- melt(tmp,id=c("I","y"))
names(ggtmp) 		<- c("I","y","Q_index","c")
ggtmp$Q_index 		<- as.numeric(as.character(ggtmp$Q_index))
ggtmp$k				<- rep(1:K, each=nQ)[as.numeric(ggtmp$Q_index)]
ggtmp$Qnodes_index	<- with(ggtmp, Q_index - (k-1)*nQ)
ggtmp$Qnodes 		<- Q_nodes[ggtmp$Qnodes_index]
ggtmp$mu 			<- diag(Qmu[ggtmp$y, ggtmp$k])
ggtmp$sigma 		<- diag(Qsigma[ggtmp$y, ggtmp$k])
ggtmp$Q				<- with(ggtmp, exp(sqrt(2)*sigma*Qnodes + mu))

# 3D plot of policy function c along I and Q
sel <- ggtmp$y == 1
quartz()
par(mfrow = c(1,2))
scatterplot3d(ggtmp[sel,"I"], ggtmp[sel,"Q"], ggtmp[sel,"c"], xlab = "I", ylab="Q", zlab="c", 
			main="policy funciton c at y=1")
scatterplot3d(ggtmp[!sel,"I"], ggtmp[!sel,"Q"], ggtmp[!sel,"c"], xlab = "I", ylab="Q", zlab="c", 
			main="policy funciton c at y=2")
			
# 2D plot of policy function at certain points of Q
sel <- 1:4
quartz()
ggplot(subset(ggtmp, Q %in% unique(ggtmp$Q)[sel]), aes(I,c)) + geom_point() + geom_line() + 
		facet_grid(Q~y, labeller = "label_both") + 
		geom_abline(aes(intercept = Q, slope = 1), linetype=2)

# Check the concavity of value function along c if using value function iteration. 
plot_state <- 71:80
plot_cvalue_fn(DP_list, c_grid, plot_state, k=3, plot = TRUE, control_list = control_list, plotnodes = 1:3)

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
delta 		<- -2			# The change value
param_list1	<- param_list
param1		<- param
param_list1[[sel_param]] <- param_list[[sel_param]] + delta
param1[sel_param] <- param[sel_param] + delta
DP_list1	<- DP_list 
DP_list1$param <- param1
DP_list1$param_list <- param_list1
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

# CCP computed at the whole sequence of simulated choices;
set.seed(99)
init_state	<- as.numeric(state[10,])
TT 			<- 500
tmp 		<- simulate_seq_fnC(init_state, TT, draw_all=TRUE, y_seq=rep(NA,TT), Q_seq=rep(NA,TT), 
							k_seq=rep(NA, TT),DP_list, control_list)  
tmp1 		<- simulate_seq_fnC(init_state, TT, draw_all=TRUE, y_seq=rep(NA,TT), Q_seq=rep(NA,TT), 
							k_seq=rep(NA, TT),DP_list1, control_list) 
ggtmp 		<- rbind(data.frame(tmp[[1]],tmp$ccp,Parameter = "Baseline"), data.frame(tmp1[[1]],tmp1$ccp,Parameter = "Experiment"))
ggtmp 		<- melt(ggtmp, id.var=c("Parameter","t","I","y","k"), measure.vars = c("X1","X2","X3"))					
names(ggtmp)<- c("Parameter","t","I","y","k_draw","k","Probability")
quartz()
ggplot(ggtmp, aes(t, Probability)) + geom_point() + geom_line() + 
		facet_grid(Parameter~k, labeller = "label_both") + 
		labs( title = paste("CCP over simulated states\n Baseline: ", paste(names(param), param, collapse=",", sep="="), 
					"\n Experiment: ", paste(names(param1), param1, collapse=",", sep="="), sep="") )

#####################
# Simulate the data # 
#####################
set.seed(99)
TT <- 500
init_state <- as.numeric(state[10,])
# my_data0 <- simulate_seq_fn(init_state,TT,draw_all=TRUE, y_seq=NULL, Q_seq=NULL, k_seq=NULL, DP_list, control_list)
my_data <- simulate_seq_fnC(init_state, TT, draw_all=TRUE, y_seq=rep(NA,TT), Q_seq=rep(NA,TT), 
							k_seq=rep(NA, TT),DP_list, control_list)    

# Plot the simulated data
ggtmp <- melt(data.frame(my_data[[1]]), id="t") 
quartz()
ggplot(ggtmp,aes(t, value)) + geom_point() + geom_line() + 
		facet_grid(variable ~ .)
		
# Compare the imputed states vs. true states
tmp <- simulate_seq_fnC(init_state, TT, draw_all=FALSE, y_seq=my_data[[1]][,"y"], Q_seq=my_data[[1]][,"Q"], 
							k_seq=my_data[[1]][,"k"], DP_list, control_list)  
identical(my_data, tmp)

DP_init 		<- DP_list
DP_init$param 	<- rep(0, length(param))
DataState 		<- as.matrix(my_data[[1]][,c("I","y")])
Q_seq			<- as.vector(my_data[[1]][,"Q"])
choice_seq		<- as.vector(my_data[[1]][,"k"])

#############################
# Likelihood function shape # 
#############################
my_grid <- list(seq(-2, 2, by=.1), seq(-.3, 2, by=0.1), seq(-2,2,by=.1))
ggtmp <- data.frame(NULL)
for(i in 1:length(param)){
	tmp <- NULL
	for(j in 1:length(my_grid[[i]])){
		theta_init <- param
		theta_init[i] <- theta_init[i]+my_grid[[i]][j]
		ll 	<- sum(ll_fnC(theta_init, DataState, choice_seq, Q_seq,  
						DP_init, Bellman_operator = Bellman_operatorC, control_list = control_list))
		tmp <- rbind(tmp, c(theta_init, ll))
	}
	ggtmp <- rbind(ggtmp, data.frame(Variable = names(param)[i], tmp))
}
colnames(ggtmp) <- c("Variable", names(param), "ll")
ggtmp$par <- with(ggtmp, ifelse(Variable == "lambda", lambda, ifelse(Variable=="tau1", tau1, tau3)) )
ggtmp$true <- param[ggtmp$Variable]

ggplot(ggtmp, aes(par, ll)) + geom_point() + geom_line() + 
		facet_wrap(~Variable) + 
		geom_vline(aes(xintercept = true), col="red", linetype= 2) + 
		labs(title = "Likelihood function")

##############
# Estimation #
##############
theta_init 		<- param
theta_init[1] 	<- param[1] + 1
theta_init[2]	<- param[2] + .2


system.time(tmp <- ll_fnC(theta_init, DataState, choice_seq, Q_seq,  
						DP_init, Bellman_operator = Bellman_operatorC, control_list = control_list))

pct		<- proc.time()
optsol 	<- maxLik(ll_fnC, DataState = DataState, choice_seq = choice_seq, Q_seq = Q_seq,
						DP_init = DP_init, Bellman_operator = Bellman_operatorC,control_list = control_list,
				start = theta_init,method ="BFGS", print.level=4, fixed = c(F,F,F))
proc.time() - pct


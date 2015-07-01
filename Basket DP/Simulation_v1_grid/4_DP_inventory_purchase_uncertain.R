# Dynamic inventory model version 3
# Discrete basket choice and continous consumption choice

# MODEL
# v(I) = max_{k,c} U(c,I,Q_k) + beta*Ev(I'|I,c,Q_k)
# where Q_k is a fixed number 
# U(c,I,Q) = lambda*log(c) - tau1*(I+Q-c) -tau2*(I+Q-c)^2

library(ggplot2)
library(reshape2)

setwd("~/Documents/Research/Store switching/Exercise/Basket DP/Simulation_v1_grid")
source("policy_iteration_functions.R")
source("4_utility_functions.R")

# Set parameters
set_list <- list(	value_iter_max 	= 20,
					policy_iter_max = 10,
					tol				= 10e-8,
					display.freq	= 5,
					value_fn_max	= 10e4
				)
beta <- .9
param_list	<- list(lambda = .3,
					tau1 	= 1,
					tau2	= 1.5
					)

step	<- .05					
I_grid	<- seq(.05,5,by=.05)
K		<- 3
c_grid	<- seq(0,7,by=step)
ns		<- length(I_grid)
Q_kernal<- list(nodes = c(0,1,1.5,2), 
				weight = matrix(c(	1,	0,	0,	0,
									0,	.8,	.2,	0,
									.2,	.1,	.2,	.5), byrow=T,K,4))
nQ		<- length(Q_kernal$nodes)

tmp		<- floor(quantile(1:ns, c(1:K)/K))
tmp		<- c(tmp[1], diff(tmp))
policy_k<- rev(unlist(lapply(1:K,function(i) rep(i,tmp[i])) ))

DP_list	<- list(state 		= as.matrix(I_grid), 
				value_fn 	= 0.5*log(I_grid+1),
				policy	 	= list(k = policy_k,c = I_grid %*% t(rep(1,nQ))),
				param_list	= param_list,
				status 		= 1)
plot_state <- 11:20
sol <- policy_iteration_fn(DP_list,u_fn=flow_Eu_fn,transition_fn=transition_fn,policy_improve_fn=solve_kc_fn,print_level=1,plot_state=plot_state)
DP_list <- sol$DP_list

# Plot value function and policy function
par(mfrow = c(2,1))
plot(DP_list$state, DP_list$value_fn, xlab = "I", ylab="V", main="Value funciton")
plot(DP_list$state, DP_list$policy$k, xlab = "I", ylab="k", main="Policy funciton k")

par(mfrow=c(2,2))
for(i in 1:nQ){
	plot(DP_list$state, DP_list$policy$c[,i], xlab = "I", ylab="c", main=paste("Policy function k at Q=",Q_kernal$nodes[i],sep=""))	
}


ggtmp <- sol$plot_data
ggtmp$I <- DP_list$state[ggtmp$state_idx,]
ni <- unique(ggtmp$Iteration)
for(i in 1:length(ni)){
	quartz()
	print(ggplot(subset(ggtmp,Iteration==ni[i]), aes(c, value)) + geom_point() + geom_line() + 
				# facet_grid(.~I,labeller = "label_both") + 
				facet_grid(Q~I) + 
				scale_color_manual(values=c("black","red")) + 
				labs(x = "c", title = paste("Iteration ",ni[i],sep="")) + 
				guides(color=FALSE)
	)
}




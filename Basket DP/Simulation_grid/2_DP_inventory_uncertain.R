# Dynamic inventory model version 2
# Single state of inventory with exogenously uncertain quantity 

# MODEL
# v(I) = max_{c} U(c,I+Q) + beta*Ev(I'|I,c)
# where Q is an exogenous quantity flow, and the expectation is w.r.t. Q. 
# U(c,I) = lambda*log(c) - tau1*(I+Q-c) -tau2*(I+Q-c)^2

# Change state x = I + Q
# W(x) = max_{c} U(c,x) + beta*EW(x'|x,c)

library(ggplot2)
library(reshape2)

setwd("~/Documents/Research/Store switching/Exercise/Basket DP/Simulation_v1_grid")
source("policy_iteration_functions.R")
source("2_utility_functions.R")

# Set parameters
set_list <- list(	value_iter_max 	= 20,
					policy_iter_max = 20,
					tol				= 10e-8,
					display.freq	= 5,
					value_fn_max	= 10e4
				)
beta <- .8
param_list	<- list(lambda = 3,
					tau1 	= .2,
					tau2	= .5
					)
					
I_grid	<- seq(0.1,5,by=.05)
Q_kernal<- list(nodes = c(0, 1, 2), 
				weight = rep(1/3,3)) 
nQ 		<- length(Q_kernal$nodes)
				
x_grid	<- seq(.1,8,by=.05)
c_grid	<- x_grid
ns		<- length(x_grid)

DP_list	<- list(state 		= as.matrix(x_grid), 
				value_fn 	= 3*log(x_grid),
				policy	 	= x_grid,
				param_list	= param_list,
				status 		= 1)
plot_state <- 11:20
sol <- policy_iteration_fn(DP_list,u_fn=flowu_vec_fn,transition_fn=transition_fix_fn,policy_improve_fn=solve_c_fn,print_level=1,plot_state=plot_state)
DP_list <- sol$DP_list

# Plot value function and policy function
par(mfrow = c(2,1))
plot(DP_list$state, DP_list$value_fn, xlab = "x", ylab="V", main="Value funciton")
plot(DP_list$state, DP_list$policy, xlab = "x", ylab="c", main="Policy funciton")

ggtmp <- sol$plot_data
ggtmp$c <- as.numeric(as.character(ggtmp$variable))
ggtmp$x <- DP_list$state[ggtmp$state,]
ggtmp$max <- with(ggtmp,ifelse(max_pnt==variable,1,0))
ni <- unique(ggtmp$Iteration)
for(i in 1:length(ni)){
	quartz()
	print(ggplot(subset(ggtmp,Iteration==ni[i]), aes(c, value)) + geom_point(aes(col=factor(max))) + geom_line() + 
				# facet_grid(.~I,labeller = "label_both") + 
				facet_wrap(~x) + 
				scale_color_manual(values=c("black","red")) + 
				labs(x = "c", title = paste("Iteration ",ni[i],sep="")) + 
				guides(color=FALSE)
	)
}




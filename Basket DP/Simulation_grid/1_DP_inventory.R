# Dynamic inventory model version 1
# Single state of inventory and consumption

# MODEL
# W(I) = max_{c} U(c,I) + beta*EW(I'|I,c)
# U(c,I) = lambda*log(c) - tau1*(I+Q-c) -tau2*(I+Q-c)^2

# NOTE: 
# When I use policy function iteration with grid search for c, policy function c is not smooth. 
# When I use policy function iteration with Brent's method, policy function c should be smooth, but the policy evaluation step still 
# raise the issue of discreteness. 
# Value function iteration solves it. 

library(ggplot2)
library(reshape2)
library(scatterplot3d)

setwd("~/Documents/Research/Store switching/Exercise/Basket DP/Simulation_v1_grid")
source("policy_iteration_functions.R")
source("1_utility_functions.R")
source("value_function_iteration.R")
plot_path <- "~/Desktop/1_graph.pdf"

# Set parameters
set_list <- list(	value_iter_max 	= 100,
					policy_iter_max = 100,
					tol				= 10e-6,
					display.freq	= 5,
					value_fn_max	= 10e4
				)
beta <- .8
param_list	<- list(lambda = 3,
					tau1 	= .2,
					tau2	= 0
					)
Q 		<- 1
c_grid 	<- seq(0,6,by=.1)
I_grid	<- seq(0,5,by=.1)
ns		<- length(I_grid)

DP_list	<- list(state 		= as.matrix(I_grid), 
				value_fn 	= 3*log(I_grid+.1),
				policy	 	= I_grid,
				param_list	= param_list,
				status 		= 1)

# Policy iteration with grids.
plot_state <- 161:170
sol <- policy_iteration_fn(DP_list,u_fn=flowu_vec_fn,transition_fn=transition_fix_fn,policy_improve_fn=solve_c_fn,print_level=1,plot_state=plot_state)
DP_list <- sol$DP_list

# Policy iteration using Brent's method to solve consumption
sol <- policy_iteration_fn(DP_list,u_fn=flowu_vec_fn,transition_fn=transition_fix_fn,policy_improve_fn=solve_c_brent_fn,print_level=1,plot_state=NULL)
DP_list <- sol$DP_list

# Value function iteration
sol <- value_iteration_fn(DP_list,Bellman_operator = solve_c_brent_fn, print_level = 1)
DP_list <- sol

# Plot value function and policy function
par(mfrow = c(2,1))
plot(DP_list$state, DP_list$value_fn, xlab = "I", ylab="V", main="Value funciton")
plot(DP_list$state, DP_list$policy, xlab = "I", ylab="c", main="Policy funciton")

ggtmp <- sol$plot_data
ggtmp$c <- as.numeric(as.character(ggtmp$variable))
ggtmp$I <- DP_list$state[ggtmp$state,]
ggtmp$max <- with(ggtmp,ifelse(max_pnt==variable,1,0))
ni <- unique(ggtmp$Iteration)
for(i in 1:length(ni)){
	quartz()
	print(ggplot(subset(ggtmp,Iteration==ni[i]), aes(c, value)) + geom_point(aes(col=factor(max))) + geom_line() + 
				# facet_grid(.~I,labeller = "label_both") + 
				facet_wrap(~I) + 
				scale_color_manual(values=c("black","red")) + 
				labs(x = "c", title = paste("Iteration ",ni[i],sep="")) + 
				guides(color=FALSE)
	)
}



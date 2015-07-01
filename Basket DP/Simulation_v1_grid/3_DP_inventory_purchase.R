# Dynamic inventory model version 3
# Discrete basket choice and continous consumption choice

# MODEL
# v(I) = max_{k,c} U(c,I,Q_k) + beta*Ev(I'|I,c,Q_k)
# where Q_k is a fixed number 
# U(c,I,Q) = lambda*log(c) - tau1*(I+Q-c) -tau2*(I+Q-c)^2 - tau3*Q

library(ggplot2)
library(reshape2)
library(Rcpp)
library(RcppGSL)

setwd("~/Documents/Research/Store switching/Exercise/Basket DP/Simulation_v1_grid")
# source("policy_iteration_functions.R")
source("3_utility_functions.R")
source("value_function_iteration.R")
sourceCpp("3C_Bellman_operator.cpp")

# Set parameters
set_list <- list(	value_iter_max 	= 500,
					tol				= 1e-10,
					display.freq	= 20,
					value_fn_max	= 10e4
				)
beta <- .9
param_list	<- list(lambda = 10,
					tau1 	= 1,
					tau3	= 3
					)

step	<- .1					
I_grid	<- seq(0,5,by=step)
K		<- 3
Q_grid	<- c(0,1,2)
# K		<- 1
# Q_grid	<- .5
c_grid	<- seq(0,7,by=step)
ns		<- length(I_grid)

tmp		<- floor(quantile(1:ns, c(1:K)/K))
tmp		<- c(tmp[1], diff(tmp))
policy_k<- rev(unlist(lapply(1:K,function(i) rep(i,tmp[i])) ))

# tmp		<- sapply(1:ns, function(i) singleu_fn(c= I_grid[i], I=I_grid[i], Q=Q_grid, param_list ))/(1-beta)
DP_list	<- list(state 		= as.matrix(I_grid), 
				value_fn 	= 3*log(I_grid+1),
# 				value_fn	= tmp,
				policy	 	= list(k = policy_k,c = I_grid %*% t(rep(1,K))),
				param_list	= param_list,
				status 		= 1, 
				beta		= beta, 
				Q_grid		= Q_grid)
				
# Value function iteration with Brent's method
system.time(sol <- value_iteration_fn(DP_list, Bellman_operator = solve_brent_fn, print_level=1))
DP_list <- sol

# Value function iteration with Brent's method written in cpp
system.time(sol1 <- value_iteration_fn(DP_list, Bellman_operator = Bellman_operatorC, print_level=1))
DP_list <- sol1

# Test if the C function matches the R function
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

# Plot value function and policy function
policy_c <- diag(DP_list$policy$c[,DP_list$policy$k])
par(mfrow = c(3,1))
plot(DP_list$state, DP_list$value_fn, xlab = "x", ylab="V", main="Value funciton")
plot(DP_list$state, DP_list$policy$k, xlab = "x", ylab="k", main="Policy funciton k")
plot(DP_list$state, policy_c, xlab = "x", ylab="c", main="Policy funciton c")

# Policy function c at each realization Q
ggtmp <- melt(DP_list$policy$c)
names(ggtmp) <- c("state_idx","Q_index","c")
ggtmp$I <- I_grid[ggtmp$state_idx]
ggtmp$Q	<- Q_grid[ggtmp$Q_index]
quartz()
ggplot(ggtmp, aes(I, c)) + geom_point() + geom_line() + 
		geom_abline(aes(intercept = Q, slope = 1),linetype = 2) + 
		facet_grid(~Q, labeller="label_both") + labs(title = "Policy function c")

# Check whether inner objective function is concave.
plot_state 	<- 1:20
Q_index 	<- 1
plot_cvalue_fn(sol, c_grid,plot_state, Q_index)

plot_state 	<- 6:9
Q_index 	<- 1
c_grid <- seq(0, 5, by=.01)
my_range <- c(.7, 1)

plotdata <- plot_cvalue_fn(sol, c_grid,plot_state, Q_index, plot=F)
tmp <- subset(plotdata, c>= my_range[1] & c<= my_range[2])
tmp1 <- dcast(tmp, c ~ I, value.var="value")
plotdata[plotdata$max_pnt==1,]
apply(tmp1, 2, max, na.rm=T)
my_yrange <- range(tmp$value, na.rm=T)
ggplot(plotdata, aes(c, value, col = factor(I))) + geom_point(aes(shape=factor(max_pnt))) + 
		geom_line() + 
		scale_shape_manual(values = c(1,16)) + 
		xlim(my_range) + ylim(my_yrange)


W <- DP_list$value_fn
V_fn <- approxfun(I_grid,W)
value_fn <- function(c,I,Q){
	I_next <- I+Q-c
	out <- singleu_fn(c,I,Q,DP_list$param_list) + beta*V_fn(I_next)
	return(out)
}
tmpc <- .8
cat("Flow utility function at each state\n")
singleu_fn(c=tmpc, I=I_grid[plot_state[1]], Q=Q_grid, param_list)
singleu_fn(c=tmpc, I=I_grid[plot_state[2]], Q=Q_grid, param_list)
V_fn(I_grid[plot_state[1]] + Q_grid - tmpc)
V_fn(I_grid[plot_state[2]] + Q_grid - tmpc)
# Dynamic discrete choice model solution 

library(rgl)
library(ggplot2)
library(reshape2)
library(maxLik)

setwd("~/Documents/Research/Store switching/Exercise/Basket DP")
setwd("//tsclient/Resear1/Store switching/Exercise/Basket DP")

source('Chebychev interpolation functions.R')
source('Gaussian_Hermite.R')
source('Policy_iteration_functions.R')
source('Simulation/Value_iteration_functions.R')
source('Simulation/Model_simulation_function.R')
source('Simulation/Estimation_functions.R')

### Set parameters ###
set_list <- list(	value.iter.max 	= 20,
					policy.iter.max = 1000,
					tol				= 10e-8,
					display.freq	= 50,
					Phi.bound		= 10e-8
				)
param_list <- list(	lambda1	= 9,
					lambda2	= -1,
					tau1 	= 3,
					tau2 	= -2,
					beta 	= .9,
					sigma_y = 1.5)
K <- 3

# Set purchase utility function and transition function
omega_fn <- function(k,y){
	alpha1 <- c(1,2,3)
	alpha2 <- c(-.2,-.1,-.2)
	if(y < 7 & k==K){
		return(0)
	}else if(y< 4 & k==(K-1)){
		return(0)
	}else{
		return(alpha1[k]*y-alpha2[k]*y^2)
	}
}			

mu_Q_fn <- function(k,y){
	alpha1 <- c(.5,1,1.5)
	if(y < 7 & k==K){
		return(0.1)
	}else if(y< 4 & k==(K-1)){
		return(0.1)
	}else{
		return(alpha1[k]*y)
	}
}
sigma_Q <- c(.1,.5,1)
		
# Initiate Gaussian-Hermite quadruture #
state.name 	<- c('I','y')
I.idx <- 1
y.idx <- 2
state_dim 	<- length(state.name)
GH_num_nodes<- 5
GH_init 	<- Gaussian_Hermite_init(GH_num_nodes,state_dim)

# Initiate Chebychev approximation 
state_lower	<- c(0, 0)
state_upper	<- c(20, 20)
cheby_dgr	<- 3
cheby_num_nodes <- 6
cheby_init 	<- cheb_init(cheby_dgr,cheby_num_nodes,state_dim,state_lower,state_upper)

# Initiate dynamic programming pre-computation 
DP_init <- Bellman_init(param_list,K,GH_init,cheby_init,verbose=TRUE) 
DP_list <- policy_iteration(DP_init,verbose=TRUE)

# Plot value function
plot3d(x=DP_list$state[,I.idx],y=DP_list$state[,y.idx],DP_list$value.fn,xlab="I",ylab="y",zlab="Value function")

#################
# Simulate data # 
#################
# Simulate a single agent problem
set.seed(6)
t 	<- 200
y 	<- 7
for(i in 2:t){
	y.draw <- rnorm(1,y[(i-1)],param_list$sigma_y)
	y <- c(y,min(max(y.draw,state_lower[y.idx]),state_upper[y.idx]))
}
data <- choice_sequence(5,y,DP_list)

# Plot path
ggtmp <- data.frame(t=1:t,choice=data$choice,Inventory=data$state[,I.idx],y=y,Quantity=data$Q,Consumption=data$consumption)
ggtmp <- melt(ggtmp,id=c("t","choice"))
ggplot(ggtmp,aes(t,value)) + geom_point(aes(col=factor(choice))) + geom_line() + 
		facet_wrap(~variable,ncol=1)

# Plot choice probability
plot3d(x=data$state[,I.idx],y=data$state[,y.idx],data$probability[,1],xlab="I",ylab="y",zlab="Pr[1]")
plot3d(x=data$state[,I.idx],y=data$state[,y.idx],data$probability[,2],xlab="I",ylab="y",zlab="Pr[2]")
plot3d(x=data$state[,I.idx],y=data$state[,y.idx],data$probability[,3],xlab="I",ylab="y",zlab="Pr[3]")

##############
# Estiamtion #
##############
param.init <- c(lambda1	= 2,
				lambda2	= -0.3,
				tau1 	= 1,
				tau2 	= -.5)
tmp <- DP_likelihood(param.init,data=data)
sol <- maxLik(DP_likelihood,start=param.init,data=data,method="BHHH")




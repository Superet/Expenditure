# Dynamic allocation model 

# Model: 
# V(I,y,nu) = max_{Q,c}[E_{epsilon}(max_{e_s} \sum_{s=1}^{S+1} uP(e_s,epsilon_s)) + uC(c) - C(I') + nu + beta*E_{I',y',nu'}V(I',y',nu')]
# 		s.t. sum_{s=1}^{S+1} <= y, sum_{s=1}^{S} e_s/p_s + q_{S+2} <= Q, c<=I+Q, Q>=0
# where 
# uP		...	the purchase utility, uP(e_s, epsilon_s) = gamma_s*psi_s*exp(epsilon_s)*log(e_s/(p_s*gamma_s)+1); 
# 			For identification reason, we set epsilon_{S+1} = epsilon_{S+1} = 0
# 			Parameterization of psi: psi_s = exp(X*beta_s)
# uC		... consumption utility. uC(c) = lambda * log(c+0.01)
# C(I')	... inventory cost. C(I) = tau1*I
# 
# Transition:
# I'	= I + Q - c
# y'	= pi(y)		first order Markorvian process. 
# 
# Two-step calculation 
# Step 1: inclusive value of expected purchase utility
# omega(y,Q) = E_{epsilon}(max_{e_s} \sum_{s=1}^{S+1} uP(e_s,epsilon_s))
# 			s.t. sum_{s=1}^{S+1} <= y, sum_{s=1}^{S} e_s/p_s + q_{S+2} <= Q
# 			
# Step 2: Dynamic programming of basket quantity and consumption 
# V(I,y,nu) = max_{Q,c}[omega(y,Q) + uC(c) - C(I') + nu + beta*E_{I',y',nu'}V(I',y',nu')]
# 		s.t. c<=I+Q, Q>=0
# 
# Versions: 
# This is V1: y is discrete, Q is discrete, gamma's are set to be 1. 

library(ggplot2)
library(reshape2)
library(scatterplot3d)
library(maxLik)
library(Rcpp)
library(RcppGSL)
library(glmmML)

setwd("~/Documents/Research/Store switching/Exercise/Basket DP/Dynamic allocation simulation")
source("DAM_1_functions.R")
source("value_function_iteration.R")

##################
# Set parameters # 
##################
control_list <- list(	value_max_iter 	= 500,
						tol				= 1e-8,
						display_freq	= 20,
						brent_tol		= 0.00001, 
						brent_max_iter 	= 300,
						inter_spline	= TRUE
				)

K <- 4					# Number of basket types
S <- 3					# Number of retailers
Sa <- S+2
alpha <- c(	-3, -3, -3)
gamma <- c(1,1,1)
param <- c(lambda = 5, tau1=1)					
beta <- .9

step	<- .5				
I_grid	<- seq(0,5,by=step)
y_grid	<- c(1,2)
y_value <- c(1,2)
c_grid	<- seq(0,8,by=step)
state	<- cbind(I=rep(I_grid,length(y_grid)),y=rep(y_grid,each=length(I_grid)))	# state matrix
s_idx	<- cbind(I=rep(1:length(I_grid),length(y_grid)), y=rep(1:length(y_grid),each=length(I_grid)))
ns		<- nrow(state)
ny		<- length(y_grid)
I_idx	<- 1
y_idx	<- 2
y_kernal<- list(nodes = c(1,2), weight = matrix(c(.8,.2,.2,.8),2,2))
price	<- c(2,1.2,1.5)			# Price 
Q_grid	<- c(0,.5,1,1.5)		# Q have K levels 
nQ		<- K

# Simulate inclusive value 
set.seed(99)
numsim 		<- 50
eps_draw	<- matrix(rnorm(numsim*S), numsim, S)
psi_draw 	<- exp(rep(1,numsim)%*%t(alpha) + eps_draw)
psi_ext 	<- cbind(psi_draw, 1, 1)
gamma_ext	<- c(gamma, 1, 1)
omega_draw	<- array(NA, c(numsim, ny, K))
e_draw 		<- array(NA, c(numsim, ny, K, Sa))
my.silent 	<- TRUE
for(i in 1:numsim){
	for(j in 1:ny){
		for(l in 1:K){
			tmp_inits <- list(rep(.01, Sa), c(rep(.01,S), y_value[j]-sum(price*.01) - .01, .01) )
			sol <- Allocation_fn(y=y_value[j],Q=Q_grid[l],psi=psi_ext[i,], gamma=gamma_ext,
								price=price,S=S, Sa, inits=tmp_inits, silent = my.silent)
			omega_draw[i,j,l] <- sol$max
			e_draw[i,j,l,] <- sol$e
		}
	}
	print(i)
}
sum(is.na(omega_draw))/length(omega_draw)
omega <- apply(omega_draw, c(2,3), mean, na.rm=T)

####################
# Dynamic solution #
####################
tmp		<- floor(quantile(1:ns, c(1:K)/K))
tmp		<- c(tmp[1], diff(tmp))
policy_k<- rev(unlist(lapply(1:K,function(i) rep(i,tmp[i])) ))

DP_list	<- list(state 		= as.matrix(state), K = K,
				value_fn 	= 0.5*log(state[,I_idx]+1) ,
				policy	 	= list(k = policy_k,c = state[,I_idx] %*% t(rep(1,nQ))),
				param		= param,
				beta		= beta,
				state_idx 	= s_idx,
				status 		= 1, 
				Q_grid		= Q_grid,
				y_kernal	= y_kernal,
				omega		= omega,
				Iteration 	= 0 )

sol <- value_iteration_fn(DP_list, Bellman_operator = Bellman_operator)
DP_list <- sol

# Plot value function 
scatterplot3d(DP_list$state[,I_idx], DP_list$state[,y_idx], DP_list$value_fn, xlab = "I", ylab="y", zlab="value", 
			main="Value funciton")
quartz()
sel <- state[,y_idx]==2
plot(state[sel,I_idx], DP_list$value_fn[sel])

# Plot policy function ccp
tmp <- Bellman_operator(DP_list)
tmp1 <- data.frame(state, k1=tmp$ccp[,1], k2=tmp$ccp[,2], k3=tmp$ccp[,3],k4=tmp$ccp[,4])
ggtmp <- melt(tmp1, id.var=c("I","y"))
names(ggtmp)[3] <- "policy"
quartz()
ggplot(ggtmp, aes(I, value)) + geom_point() + geom_line() + 
		facet_grid(y ~ policy, labeller = "label_both") + 
		labs(y="Probability", title = "Policy function k(CCP) " )

# Plot policy function c
tmp 			<- DP_list$policy$c
dimnames(tmp) 	<- list(1:ns, Q_grid)
ggtmp 			<- melt(tmp)
names(ggtmp) 	<- c("state_index","Q","c")
ggtmp$I			<- state[ggtmp$state_index, I_idx]
ggtmp$y			<- state[ggtmp$state_index, y_idx]

quartz()
print(ggplot(ggtmp, aes(I, c)) + geom_point() + geom_line() + 
		facet_grid(Q ~ y, labeller = "label_both") + 
		geom_abline(aes(intercept = Q, slope = 1), linetype=2) + 
		labs(title = paste("Policy function c at y=",i,sep=""))
)

# Check the concavity of value function when solving for continuous consumption;
plot_state <- 10:15
plot_cvalue_fn(DP_list, c_grid, plot_state=plot_state, Q_index=NULL, control_list)

#################
# Simulate data #
#################
init_state <- state[10,]	# Initial state
TT		<- 100				# The length of simulating time series
nu_draw	<- rlnorm(TT)

# No randomness in the consumption function 
my_data <- simulate_seq_fn(init_state, TT, draw_all = TRUE, y_seq=NULL, Q_seq=NULL, k_seq=NULL, DP_list, control_list, 
							alpha, gamma, omega)
# Randomness in the consumption function
my_data <- simulate_seq_fn(init_state, TT, draw_all = TRUE, y_seq=NULL, Q_seq=NULL, k_seq=NULL, DP_list, control_list, 
							alpha, gamma, omega, nu=nu_draw)

# Plot the simulation
ggtmp <- my_data[[1]]
ggtmp1 <- melt(data.frame(ggtmp[,c("t","I","y","k","c")]), id="t") 
quartz()
ggplot(ggtmp1,aes(t, value)) + geom_point() + geom_line() + 
		facet_grid(variable ~ ., scales="free_y")

ggtmp$e.1 <- ggtmp$e.1/y_value[ggtmp$y]
ggtmp$e.2 <- ggtmp$e.2/y_value[ggtmp$y]
ggtmp$e.3 <- ggtmp$e.3/y_value[ggtmp$y]
ggtmp$e.4 <- ggtmp$e.4/y_value[ggtmp$y]
ggtmp2 <- melt(ggtmp[,c("t","e.1","e.2","e.3","e.4")], id="t")
quartz()
ggplot(ggtmp2, aes(t, value)) + geom_point() + geom_line() + 
		facet_grid(.~variable)		

############################
# Identification intuition # 
############################
tmp <- solve_c_brent_fn(DP_list)
ccp_base<- data.frame(state, k1=tmp$ccp[,1], k2=tmp$ccp[,2], k3=tmp$ccp[,3])

sel_param 	<- "lambda"		# The parameter name that is changed
delta 		<- -2			# The change value
param1		<- param
param1[sel_param] <- param[sel_param] + delta
DP_list1	<- DP_list 
DP_list1$param <- param1
DP_list1 	<- value_iteration_fn(DP_list1, Bellman_operator = solve_c_brent_fn)
tmp 		<- solve_c_brent_fn(DP_list1)
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
tmp1 		<- simulate_seq_fn(init_state, TT, draw_all=TRUE, y_seq=rep(NA,TT), Q_seq=rep(NA,TT), 
							k_seq=rep(NA, TT),DP_list1, control_list, alpha, gamma, omega) 
ggtmp 		<- rbind(data.frame(my_data[[1]],my_data$ccp,Parameter = "Baseline"), 
					 data.frame(tmp1[[1]],tmp1$ccp,Parameter = "Experiment"))
ggtmp 		<- melt(ggtmp, id.var=c("Parameter","t","I","y","k"), measure.vars = c("X1","X2","X3"))					
names(ggtmp)<- c("Parameter","t","I","y","k_draw","k","Probability")
quartz()
ggplot(ggtmp, aes(t, Probability)) + geom_point() + geom_line() + 
		facet_grid(Parameter~k, labeller = "label_both") + 
		labs( title = paste("CCP over simulated states\n Baseline: ", paste(names(param), param, collapse=",", sep="="), 
					"\n Experiment: ", paste(names(param1), param1, collapse=",", sep="="), sep="") )


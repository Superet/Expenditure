library(ggplot2)
library(reshape2)
library(scatterplot3d)
library(maxLik)
library(Rcpp)
library(RcppGSL)
library(glmmML)

setwd("~/Documents/Research/Store switching/Exercise/Basket DP/Dynamic allocation simulation")
source("DAM_2_functions.R")
source("value_function_iteration.R")
sourceCpp("DAM_2_functions.cpp")

##################
# Set parameters # 
##################
control_list <- list(	value_max_iter 	= 500,
						tol				= 1e-6,
						display_freq	= 20,
						brent_tol		= 0.00001, 
						brent_max_iter 	= 300,
						inter_spline	= TRUE
				)

K <- 4					# Number of basket types
S <- 3					# Number of retailers
Sa <- S+2
alpha <- c(	-1, -2, -1)
alpha_a <- c(-2, -2) 	# Normalization of the parameters for the outside options
qz_cons <- 0			# The constant parameter in the utility of quantity outside good.
gamma <- c(1,1,1)
param <- c(lambda = 3, tau1=1, tau2=0)					
beta <- .9

step	<- .5				
I_grid	<- seq(0,5,by=step)
c_grid	<- seq(0,8,by=step)
Q_grid	<- c(0,.5,1.2,2)/2		# Q have K levels 
price	<- c(2,3,2)*2				# Price 
nQ		<- K

gh_num_nodes <- 6
y_weight<- ghq(gh_num_nodes, modified = F)$weights
y_nodes <- ghq(gh_num_nodes, modified = F)$zeros
ymu		<- .2
ysigma 	<- .2
y_kernal <- list(nodes = y_nodes, weight = y_weight, mu=ymu, sigma = ysigma, rho = .7)
y_grid 	<- seq(1, 5, length = 10)

# Simulate a sequence of income to check the coverage of y_grid
tmp 	<- rep(NA, 500)
tmp[1]  <- 1
for(i in 2:length(tmp)){ tmp[i] <- exp(rnorm(1, ymu + y_kernal$rho*log(tmp[(i-1)]), ysigma))}
plot(1:500, tmp)
range(tmp)

state	<- cbind(I=rep(I_grid,length(y_grid)),y=rep(y_grid,each=length(I_grid)))	# state matrix
s_idx	<- cbind(I=rep(1:length(I_grid),length(y_grid)), y=rep(1:length(y_grid),each=length(I_grid)))
ns		<- nrow(state)
I_idx	<- 1
y_idx	<- 2

# Simulate inclusive value 
set.seed(99)
numsim 		<- 50
eps_draw	<- matrix(rnorm(numsim*S), numsim, S)
psi_draw 	<- exp(rep(1,numsim)%*%t(alpha) + eps_draw)
psi_ext 	<- cbind(psi_draw, rep(1, numsim) %*% t(exp(alpha_a)))
gamma_ext	<- c(gamma, 1, 1)

# Add more y values to y_grid when computing the inclusive value
tmp_y_grid 	<- y_grid
tmp_ny 		<- length(tmp_y_grid)

omega_draw	<- array(NA, c(numsim, tmp_ny, K))
e_draw 		<- array(NA, c(numsim, tmp_ny, K, Sa))
my.silent 	<- TRUE
for(i in 1:numsim){
	for(j in 1:tmp_ny){
		for(l in 1:K){
			sol <- Allocation_fn(y=tmp_y_grid[j],Q=Q_grid[l],psi=psi_ext[i,], gamma=gamma_ext,
								price=price,S=S, Sa, silent = my.silent)
			omega_draw[i,j,l] <- sol$max
			e_draw[i,j,l,] <- sol$e
		}
	}
	print(i)
}
sum(is.na(omega_draw))/length(omega_draw)
omega <- apply(omega_draw, c(2,3), mean, na.rm=T)
dimnames(omega) <- list(tmp_y_grid, 1:K)

# Spline interpolation of omega function along y
omega_list <- lapply(1:K, function(i) splinefun(x=tmp_y_grid, y=omega[,i], method="natural"))
omega_fn <- function(k, y){
	f  <- omega_list[[k]](y)
	return(f)
}

# Plot the omega function along y
tmp <- seq(min(y_grid), max(y_grid), length=200)
ggtmp <- matrix(NA, length(tmp), K, dimnames = list(tmp, 1:K))
for(i in 1:K){
	ggtmp[,i] <- omega_fn(i, tmp)
}
ggtmp 	<- cbind(melt(ggtmp), pnt = "interp")
ggtmp1 	<- cbind(melt(omega), pnt = "knots")
ggtmp 	<- rbind(ggtmp, ggtmp1)
names(ggtmp) <- c("y","k", "omega","pnt")
ggtmp$k <- factor(ggtmp$k)
ggplot(ggtmp, aes(y, omega, col=k)) + geom_line() +
		geom_point(data = subset(ggtmp, pnt == "knots")) + 
# 		facet_wrap(~k) + 
		labs(title = expression(paste("Cubic spline interpolation of the inclusive value ", omega, "(k,y)"), sep="") ) 

####################
# Dynamic solution #
####################
tmp		<- floor(quantile(1:ns, c(1:K)/K))
tmp		<- c(tmp[1], diff(tmp))
policy_k<- rev(unlist(lapply(1:K,function(i) rep(i,tmp[i])) ))

DP_list	<- list(state 		= as.matrix(state), K = K,
				value_fn 	= 0.5*log(state[,I_idx]+1) ,
# 				value_fn	= sapply(1:ns, function(i) uflow_fn(c=state[i,1]/2, I=state[i,1], Q=0, y=state[i,2], sel_k=1, omega=omega_fn, param))/(1-beta),
				policy	 	= list(k = policy_k,c = state[,I_idx] %*% t(rep(1,nQ))),
				K			= K,
				param		= param,
				beta		= beta,
				state_idx 	= s_idx,
				status 		= 1, 
				Q_grid		= Q_grid,
				y_kernal	= y_kernal,
				omega		= omega_fn,
				price 		= price,
				Iteration 	= 0 )

solC <- value_iteration_fnC(DP_list, Bellman_operatorC, print_level=1, control_list)
DP_list <- solC
# solR <- value_iteration_fnC(DP_list, Bellman_operator, print_level=1, control_list)
# 
# # Check if c++ function and R functions return the same results;
# plot(solC$value_fn, solR$value_fn); abline(a = 0, b= 1, col="red")

# Plot value function 
scatterplot3d(DP_list$state[,I_idx], DP_list$state[,y_idx], DP_list$value_fn, xlab = "I", ylab="y", zlab="value", 
			main="Value funciton", angle = 60)

# 2D value function along I
sel <- unique(floor(length(y_grid)*c(.1, .2, .3, .5, .6,.7,1)))
sel <- y_grid[sel]
ggtmp <- data.frame(state, v=DP_list$value_fn)
quartz()
ggplot(subset(ggtmp, y %in% sel), aes(I, v) ) + geom_point() + geom_line() + 
		facet_wrap(~ y) + 
		labs(title = "Value function in I\n facet by y")

# Plot policy function ccp
tmp <- Bellman_operatorC(DP_list, control_list)
tmp1 <- data.frame(state, k1=tmp$ccp[,1], k2=tmp$ccp[,2], k3=tmp$ccp[,3],k4=tmp$ccp[,4])
ggtmp <- melt(tmp1, id.var=c("I","y"))
names(ggtmp)[3] <- "policy"
quartz()
ggplot(subset(ggtmp, y %in% sel), aes(I, value)) + geom_point() + geom_line() + 
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
print(ggplot(subset(ggtmp, y %in% sel), aes(I, c)) + geom_point() + geom_line() + 
		facet_grid(Q ~ y, labeller = "label_both") + 
		geom_abline(aes(intercept = Q, slope = 1), linetype=2) + 
		labs(title = paste("Policy function c at y grid",sep=""))
)

#################
# Simulate data #
#################
init_state <- state[30,]	# Initial state
TT		<- 200				# The length of simulating time series

# No randomness in the consumption function 
my_data1 <- simulate_seq_fnC(init_state, TT, draw_all=TRUE, y_seq=rep(NA,TT), Q_seq=rep(NA,TT), 
							k_seq=rep(NA, TT),DP_list, control_list)    
my_data <- simulate_seq_fn(init_state, TT, draw_all = TRUE, y_seq=NULL, Q_seq=NULL, k_seq=NULL, DP_list, control_list, 
							alpha, alpha_a, gamma)

# Plot the simulation
ggtmp <- my_data[[1]]
ggtmp$Q <- Q_grid[ggtmp$k]
ggtmp$c_inventory <- with(ggtmp, c/(I+Q))
ggtmp1 <- melt(data.frame(ggtmp[,c("t","I","y","k","c_inventory")]), id="t") 
quartz()
ggplot(ggtmp1,aes(t, value)) + geom_point() + geom_line() + 
		facet_grid(variable ~ ., scales="free_y") + 
		labs(title = "Simulated sequence")

ggtmp$e.1 <- ggtmp$e.1/ggtmp$y
ggtmp$e.2 <- ggtmp$e.2/ggtmp$y
ggtmp$e.3 <- ggtmp$e.3/ggtmp$y
ggtmp$e.4 <- ggtmp$e.4/ggtmp$y
ggtmp2 <- melt(ggtmp[,c("t","e.1","e.2","e.3","e.4")], id="t")
quartz()
ggplot(ggtmp2, aes(t, value)) + geom_point() + geom_line() + 
		facet_grid(.~variable)		

##################
# Identification # 
##################
tmp 	<- Bellman_operatorC(DP_list, control_list)
ccp_base<- data.frame(state, k1=tmp$ccp[,1], k2=tmp$ccp[,2], k3=tmp$ccp[,3], k4=tmp$ccp[,4])

sel_param 	<- "lambda"		# The parameter name that is changed
delta 		<- -2			# The change value
param1		<- param
param1[sel_param] <- param[sel_param] + delta
DP_list1	<- DP_list 
DP_list1$param <- param1
DP_list1 	<- value_iteration_fnC(DP_list1, Bellman_operatorC, print_level=1, control_list)
tmp			<- Bellman_operatorC(DP_list1, control_list)
ccp1		<- data.frame(state, k1=tmp$ccp[,1], k2=tmp$ccp[,2], k3=tmp$ccp[,3], k4=tmp$ccp[,4])

# Plot the CCP functions along states given the two sets of parameters
ggtmp 		<- rbind(data.frame(ccp_base, Parameter = "Baseline"), data.frame(ccp1, Parameter = "Experiment"))
ggtmp		<- melt(ggtmp, id.var=c("Parameter","I","y"))
names(ggtmp)<- c("Parameter","I","y","policy","Probability")
quartz()
ggplot(ggtmp, aes(I, Probability, col=Parameter)) + geom_point() + geom_line() + 
		facet_grid(y ~ policy, labeller = "label_both") + 
		labs( title = paste("CCP over states\n Baseline: ", paste(names(param), param, collapse=",", sep="="), 
					"\n Experiment: ", paste(names(param1), param1, collapse=",", sep="="), sep="") )

# Plot the value function from two sets of parameters. 
ggtmp <- rbind(data.frame(state, V = DP_list$value, Parameter = "Baseline"), 
			   data.frame(state, V = DP_list1$value, Parameter = "Experiment"))
ggplot(ggtmp, aes(I, V, col=Parameter)) + geom_point() + geom_line() + 
		facet_wrap(~y) + 
		labs( title = paste("Value function over states\n Baseline: ", paste(names(param), param, collapse=",", sep="="), 
					"\n Experiment: ", paste(names(param1), param1, collapse=",", sep="="), sep="") )

# CCP computed at the whole sequence of simulated choices;
set.seed(99)
init_state	<- as.numeric(state[30,])
TT 			<- 200
tmp 		<- simulate_seq_fnC(init_state, TT, draw_all=TRUE, y_seq=rep(NA,TT), Q_seq=rep(NA,TT), 
							k_seq=rep(NA, TT),DP_list, control_list)   
tmp1 		<- simulate_seq_fnC(init_state, TT, draw_all=TRUE, y_seq=rep(NA,TT), Q_seq=rep(NA,TT), 
							k_seq=rep(NA, TT),DP_list1, control_list) 
ggtmp 		<- rbind(data.frame(tmp[[1]],tmp$ccp,Parameter = "Baseline"), data.frame(tmp1[[1]],tmp1$ccp,Parameter = "Experiment"))
ggtmp 		<- melt(ggtmp, id.var=c("Parameter","t","I","y","k"), measure.vars = c("X1","X2","X3","X4"))					
names(ggtmp)<- c("Parameter","t","I","y","k_draw","k","Probability")
quartz()
ggplot(ggtmp, aes(t, Probability)) + geom_point() + geom_line() + 
		facet_grid(Parameter~k, labeller = "label_both") + 
		labs( title = paste("CCP over simulated states\n Baseline: ", paste(names(param), param, collapse=",", sep="="), 
					"\n Experiment: ", paste(names(param1), param1, collapse=",", sep="="), sep="") )

############## 					
# Estimation #
############## 
DP_init 		<- DP_list
DP_init$param 	<- rep(0, length(param))
DataState 		<- as.matrix(my_data[[1]][,c("I","y")])
# Q_seq			<- as.vector(my_data[[1]][,"Q"])
Q_seq			<- as.vector(Q_grid[my_data[[1]][,"k"]])
choice_seq		<- as.vector(my_data[[1]][,"k"])

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

#############################
# Likelihood function shape # 
#############################
my_grid <- list(seq(-1, 1, by=.2), seq(-.2, 1, by=0.2), seq(-1,1,by=.2))
ggtmp <- data.frame(NULL)
for(i in 1:length(param)){
	tmp <- NULL
	for(j in 1:length(my_grid[[i]])){
		theta_init <- param
		theta_init[i] <- theta_init[i]+my_grid[[i]][j]
		DP_init$value_fn <- 0.5*log(state[,I_idx]+1) 
		ll 	<- sum(ll_fnC(theta_init, DataState, choice_seq, Q_seq,  
						DP_init, Bellman_operator = Bellman_operatorC, control_list = control_list))
		tmp <- rbind(tmp, c(theta_init, ll))
	}
	ggtmp <- rbind(ggtmp, data.frame(Variable = names(param)[i], tmp))
}
colnames(ggtmp) <- c("Variable", names(param), "ll")
ggtmp$par <- with(ggtmp, ifelse(Variable == "lambda", lambda, ifelse(Variable=="tau1", tau1, tau2)) )
ggtmp$true <- param[ggtmp$Variable]

ggplot(ggtmp, aes(par, ll)) + geom_point() + geom_line() + 
		facet_wrap(~Variable, scales="free") + 
# 		ylim(c(-200,0)) +
		geom_vline(aes(xintercept = true), col="red", linetype= 2) + 
		labs(title = "Likelihood function")


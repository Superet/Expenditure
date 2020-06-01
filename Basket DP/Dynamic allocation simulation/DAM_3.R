library(ggplot2)
library(reshape2)
library(scatterplot3d)
library(maxLik)
library(Rcpp)
library(RcppGSL)
library(glmmML)

# sourceCpp("~/Documents/Research/Store switching/Exercise/Basket DP/0_DAM_functions.cpp")
# run_id	<- "w1t1"

setwd("/home/brgordon/ccv103/processed data")
sourceCpp("../Exercise/Basket DP/0_DAM_functions.cpp")
args <- commandArgs(trailingOnly = TRUE)
print(args)
run_id <- args
make_plot	<- FALSE

control_list <- list(	value_max_iter 	= 500,
						tol				= 1e-6,
						display_freq	= 20,
						brent_tol		= 0.00001, 
						brent_max_iter 	= 300,
						inter_spline	= TRUE
				)
				
#############
# Functions #
#############
matrix_expand <- function(mat, vec){
	n1	<- nrow(mat)
	n2	<- length(vec)
	sel	<- rep(1:n1, n2)
	mat_new <- cbind(mat[sel,], rep(vec, each=n1))
	return(mat_new)
}

my_lag 	<- function(x){
	return(c(NA, x[-length(x)]))
}

FlowU 	<- function(c, I, lny, Inc, Rc, k, Q, TR, param, price_min){
	lambda 	<- param["lambda"]
	tau_b	<- param["tau_b"]
	tau1	<- param["tau1"]
# 	tripcost<- param[(3+k)]
	tau2	<- param["tau2"]
	I_next	<- I + Q - c
	omega	<- omega_fn(k, lny)
	borrow_cost <- 0
	if(exp(lny)/price_min < Q){
		borrow_cost <-  tau_b * log(price_min*Q - exp(lny))
	}
	u 		<- lambda*log(c+0.01) - tau1*I_next - tau2*TR + omega - borrow_cost
	return(u)
}

cal_flow 	<- function(DP_list, c_grid, DataState, choice_k, cstar){
	ns		<- nrow(DataState)
	nc		<- length(c_grid)
	out		<- matrix(NA, ns, nc+1)
	Q_grid 	<- DP_list$Q_grid
	TR_grid	<- DP_list$TR_grid
	param 	<- DP_list$param
	price	<- DP_list$price
	for(i in 1:ns){
		for(j in 1:nc){
			out[i,j] <- FlowU(c_grid[j], DataState[i,1], DataState[i,2], DataState[i,3], 
							DataState[i,4], choice_k, Q_grid[choice_k], TR_grid[choice_k], 
							param, min(price))
		}
		out[i, (nc+1)] <- FlowU(cstar[i], DataState[i,1], DataState[i,2], DataState[i,3], 
							DataState[i,4], choice_k, Q_grid[choice_k], TR_grid[choice_k], 
							param, min(price))
	}
	return(out)
}

#####################
# Simulation set up #
#####################
# Parameter 
param 		<- c(lambda = 1, tau_b = 3, tau1 = .3, tau2 = .5)
if(substr(run_id, 3,4)=="t0"){
	param["tau_b"] <- 0
}
beta		<- .95
K 			<- 4
bskt_typeQ 	<- c(0.1, .5, 0.1, .5)
bskt_typeTR <- c(7,7,9,9)
price 		<- c(250, 300, 350)
n_Inc		<- 4

# Initiate states
my_grid <- list(I = c(0, 1, 2, 3, 4, 5 ), 
				y = c(3, 4.8, 5.2, 5.5, 5.8, 6),  
				Inc	= 1:n_Inc, 
				Rc = c(0, 1))
state	<- as.matrix(my_grid[[1]], ncol=1) 
for(i in 2:length(my_grid)){
	state <- matrix_expand(state, my_grid[[i]])
}
colnames(state) <- c("I", "y", "Inc","Rc")
summary(state)

# Distribution kernal of expenditure 
gh_num_nodes <- 5
y_weight	<- ghq(gh_num_nodes, modified = F)$weights
y_nodes 	<- ghq(gh_num_nodes, modified = F)$zeros
kappa0 		<- c(4, 4.2, 4.8, 4.7)
kappa1		<- kappa0 - .1
rho0		<- .1
rho1		<- .15
sigma0 		<- .15
sigma1 		<- .35
y_kernal	<- list(nodes = y_nodes, weight = y_weight, 
				kappa = cbind(kappa0, kappa1), rho = c(rho0, rho1), sigma = c(sigma0, sigma1))

# Markov probability of income transition
Inc_weight 	<- rbind(diag(n_Inc), diag(n_Inc))
Inc_weight[1,]	<- c(.75, .25,	0,	0)
Inc_weight[2,]	<- c(.1, 	.9,	0,	0)
Inc_weight[3,]	<- c(0,		.1,	.8,	.1)
Inc_weight[4,]	<- c(0,		0,	.2,	.8)
Inc_weight[5,]	<- c(.8, 	.2,	0,	0)
Inc_weight[6,]	<- c(.2, 	.8,	0,	0)
Inc_weight[7,]	<- c(0,		.2,	.7,	.1)
Inc_weight[8,]	<- c(0,		.1,	.2,	.7)


# Omega function
omega_fn 	<- function(k, y){
	if(substr(run_id, 1,2)=="w0"){
		rep(0, length(y))
	}else if (substr(run_id, 1,2)=="w1"){
		-.01*k*y^2 + .6*y
	}else{
		.015*k*y^2 + .4*y
	} 	
}
tmp 		<- seq(0, 7, by=.1)
ggtmp		<- sapply(1:K, function(i) omega_fn(i,tmp))
dimnames(ggtmp) <- list(y=tmp, k=1:K)
ggtmp		<- melt(ggtmp)
if(make_plot){
	ggplot(ggtmp, aes(y, value, col=factor(k))) + geom_point() + geom_line() 
}

####################
# Dynamic solution #
####################
# DP structure 
DP_list	<- list(state 	= state, K = K, 
value_fn= 2 * log(state[,1] + 1), 
# 				value_fn= 0, 
				param  	= param, 
				beta 	= beta,
				Q_grid	= bskt_typeQ, 
				TR_grid = bskt_typeTR, 
				omega	= omega_fn, 
				price	= price,
				y_kernal= y_kernal, 
				Inc_weight = Inc_weight,
				Iteration = 0, 
				status = 1
				)
policy_k_guess <- rep(1, nrow(state))
tmp			<- V_initC(param, DP_list, policy_k = policy_k_guess, policy_c = .5*(state[,1]+bskt_typeQ[policy_k_guess]))
# DP_list$value_fn <- tmp

solC <- value_iteration_fnC(DP_list, Bellman_operatorC, print_level=1, control_list)
DP_list <- solC

# ------------------------------ # 
# Check DP solution 
if(make_plot){
	# Plot value function 
	sel		<- DP_list$state[,3] == 1 & DP_list$state[,4] == 0
	scatterplot3d(DP_list$state[sel,1], DP_list$state[sel,2], DP_list$value_fn[sel], xlab = "I", ylab="y", zlab="value", 
				main="Value funciton", angle = 60)

	# 2D value function along I
	ggtmp <- data.frame(state, v=DP_list$value_fn)
	quartz()
	print(ggplot(ggtmp, aes(I, v, col=factor(Inc), linetype = factor(Rc)) ) + geom_point() + geom_line() + 
			facet_wrap(~ y) + 
			labs(title = "Value function in I\n facet by y")
	)
	
	# Plot policy function ccp
	tmp <- Bellman_operatorC(DP_list, control_list)
	tmp1 <- data.frame(state, tmp$ccp)
	ggtmp <- melt(tmp1, id.var=c("I","y","Inc","Rc"))
	names(ggtmp)[5] <- "policy"
	quartz()
	print( ggplot(ggtmp, aes(I, value, col=factor(Inc), linetype = factor(Rc) )) + geom_point() + geom_line() + 
			facet_grid(y ~ policy, labeller = "label_both") + 
			labs(y="Probability", title = "Policy function k(CCP) " )
	)
	
	# Plot policy function c
	tmp 			<- DP_list$policy$c
	dimnames(tmp) 	<- list(1:nrow(state), 1:K)
	ggtmp 			<- melt(tmp)
	names(ggtmp) 	<- c("state_index","k","c")
	ggtmp$I			<- state[ggtmp$state_index, 1]
	ggtmp$y			<- state[ggtmp$state_index, 2]
	ggtmp$Inc		<- state[ggtmp$state_index, 3]
	ggtmp$Rc		<- state[ggtmp$state_index, 4]
	ggtmp$Q			<- bskt_typeQ[ggtmp$k]

	quartz()
	print(ggplot(ggtmp, aes(I, c, col=factor(Inc), linetype = factor(Rc), alpha=.3 )) + geom_point() + geom_line() + 
			facet_grid(k ~ y, labeller = "label_both") + 
			geom_abline(aes(intercept = Q, slope = 1), linetype=2) + 
			labs(title = paste("Policy function c at y grid",sep=""))
	)
	
	# Shape of objective function along c
	c_grid 		<- seq(.5, 4, .2)
	plot_state 	<- 98:106
	state[plot_state,]
	choice_k 	<- 4
	tmp_cstar	<- DP_list$policy$c[plot_state, choice_k]
	tmp			<- cal_valuefn(DP_list, c_grid, state[plot_state,], choice_k, TRUE, tmp_cstar)
	tmp1 		<- cal_flow(DP_list, c_grid, state[plot_state,], choice_k, tmp_cstar)
	dimnames(tmp) <- list(s_index = plot_state, c = c(c_grid,""))
	dimnames(tmp1) <- list(s_index = plot_state, c = c(c_grid,""))
	ggtmp		<- melt(tmp)
	sel 		<- is.na(ggtmp$c)
	ggtmp[sel,"c"] <- tmp_cstar
	ggtmp$max_pnt <- 0
	ggtmp[sel,"max_pnt"] <- 1
	ggtmp1		<- melt(tmp1)
	sel 		<- is.na(ggtmp1$c)
	ggtmp1[sel,"c"] <- tmp_cstar
	ggtmp1$max_pnt <- 0
	ggtmp1[sel,"max_pnt"] <- 1	
	ggtmp 		<- data.frame(rbind(ggtmp, ggtmp1), v=rep(c("value","flow"), each=nrow(ggtmp1)) )
	ggtmp$I		<- state[ggtmp$s_index, 1]
	ggtmp$lny	<- state[ggtmp$s_index, 2]
	ggtmp$Inc	<- state[ggtmp$s_index, 3]
	ggtmp$Rc	<- state[ggtmp$s_index, 4]

	quartz()
	print(ggplot(ggtmp, aes(c, value, linetype = v, shape = v)) + geom_point(aes(col=factor(max_pnt))) + 
				geom_line() + 
				facet_wrap(I ~ lny) + 
				scale_color_manual(values=c("black","red")) + 
				guides(color=FALSE) + 
				labs(title = paste("Value function curve at Q=",
									round(bskt_typeQ[choice_k],2), "\n(facet by I,y)", sep=""))
	)
}

#################
# Simulate data #
#################
init_state 	<- state[280,]	# Initial state
TT			<- 300				# The length of simulating time series

# No randomness in the consumption function 
set.seed(666) 
my_data <- simulate_seq1_fnC(init_state, TT, draw_all=TRUE, lny_seq=rep(NA,TT), Inc_seq = rep(NA, TT), 
 							Rc_seq=rep(NA, TT), k_seq = rep(NA, TT), DP_list, control_list)
 
# Plot the simulation
if(make_plot){
	ggtmp <- my_data[[1]]
	ggtmp$Q <- bskt_typeQ[ggtmp$k]
	ggtmp$c_inventory <- with(ggtmp, c/(I+Q))
	ggtmp1 <- melt(data.frame(ggtmp[,c("t","I","lny", "Inc", "k","c_inventory")]), id="t") 
	quartz()
	ggplot(ggtmp1,aes(t, value)) + geom_point() + geom_line() + 
			facet_grid(variable ~ ., scales="free_y") + 
			labs(title = "Simulated sequence")
}

# Simulate two households
set.seed(666) 
hh_index 	<- c(1,2,3,4,5,6)
TT_vec		<- c(100, 100, 200, 100, 200, 100)
init_state 	<- as.matrix(state[c(4, 50, 190, 214, 280, 287),]) 
init_state
hhdata		<- data.frame()
for(i in 1:length(hh_index)){
	tmp		<- 	simulate_seq1_fnC(init_state[i,], TT_vec[i], draw_all=TRUE, lny_seq=rep(NA,TT_vec[i]), 
							Inc_seq = rep(NA, TT_vec[i]), Rc_seq=rep(NA, TT_vec[i]), 
							k_seq = rep(NA, TT_vec[i]), DP_list, control_list)
	hhdata	<- rbind(hhdata, tmp[[1]])
}
DataState 	<- cbind(NA, as.matrix(hhdata[,c("lny","Inc","Rc")]))
choice_seq	<- hhdata[,"k"]
cat("Data from", length(hh_index), "people.\n")

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
ggtmp		<- melt(ggtmp, id.var=c("Parameter","I","y","Inc","Rc"))
names(ggtmp)<- c("Parameter","I","y","Inc","Rc","policy","Probability")
selRc		<- 0
if(make_plot){
	for(i in 1:n_Inc){
	quartz()
	print(ggplot(subset(ggtmp, Rc == selRc & Inc==i), aes(I, Probability, col=Parameter)) + geom_point() + geom_line() + 
		facet_grid(y  ~ policy, labeller = "label_both") + 
		labs( title = paste("CCP over Inc=", i, "\n Baseline: ", paste(names(param), param, collapse=",", sep="="), 
					"\n Experiment: ", paste(names(param1), param1, collapse=",", sep="="), sep="") )
	)
	}
}

# Plot the value function from two sets of parameters. 
ggtmp <- rbind(data.frame(state, V = DP_list$value, Parameter = "Baseline"), 
			   data.frame(state, V = DP_list1$value, Parameter = "Experiment"))
selRc		<- 0
if (make_plot){
	quartz()
	ggplot(subset(ggtmp, Rc==selRc), aes(I, V, col=Parameter)) + geom_point() + geom_line() + 
			facet_grid(Inc~y, labeller = "label_both") + 
			labs( title = paste("Value function over states\n Baseline: ", paste(names(param), param, collapse=",", sep="="), 
						"\n Experiment: ", paste(names(param1), param1, collapse=",", sep="="), sep="") )
}

#############################
# Likelihood function shape # 
#############################
my_grid <- list(lambda = seq(-.9, 1, .2), tau_b = seq(-1, 1, .2), tau_1 = seq(-.3, 2, .1), 
				tau_2 = seq(-.5, .5, .1) ) 
ggtmp <- data.frame(NULL)
DP_init	<- DP_list
for(i in 1:length(my_grid)){
	prc <- proc.time()
	tmp <- NULL
	for(j in 1:length(my_grid[[i]])){
		theta_init <- param
		theta_init[i] <- theta_init[i]+my_grid[[i]][j]
		DP_init$value_fn <- 2*log(state[,1]+1) 
		ll 	<- sum(ll_fnC(theta_init, hh_index, TT_vec, init_state, DataState, choice_seq, 
						DP_init, Bellman_operatorC, control_list, policy_k_guess))
		tmp <- rbind(tmp, c(theta_init, ll))
	}
	ggtmp <- rbind(ggtmp, data.frame(Variable = names(param)[i], tmp))
	use.time	<- proc.time() - prc
	cat("Finish", names(param)[i], "with", use.time[3]/60,"min.\n")
}
colnames(ggtmp) <- c("Variable", names(param), "ll")
ggtmp$par		<- NA
for(i in 1:length(my_grid)){
	sel <- ggtmp$Variable == names(param)[i]
	if(length(sel) >0){
		ggtmp[sel,"par"] <- ggtmp[sel,names(param)[i]]
	}
}
ggtmp$true 	<- param[ggtmp$Variable]
ll_comp		<- ggtmp

pdf(paste("run_sim/simulation_DAM_", run_id,".pdf",sep=""), width = 8, height = 8)
ggplot(ggtmp, aes(par, ll)) + geom_point() + geom_line() + 
		facet_wrap(~Variable, scales="free") + 
# 		ylim(c(-200,0)) +
		geom_vline(aes(xintercept = true), col="red", linetype= 2) + 
		labs(title = "Likelihood function")
dev.off()

############## 					
# Estimation #
############## 
DP_init			<- DP_list
DP_init$value_fn <- 2*log(state[,1]+1) 
theta_init 		<- param
theta_init[1] 	<- param[1] + 1
theta_init[3]	<- param[3] + .5
theta_init[4]	<- param[4] + .5
names(theta_init) <- names(param)
cat("Initial values are\n", theta_init,"\n")
myfix 	<- rep(FALSE, length(param))
if(substr(run_id, 3,4)!="t1"){
	sel 	<- names(param) == "tau_b"
	myfix[sel]	<- TRUE
}		
cat("Parameters that are fixed:", myfix,"\n")

system.time(tmp <- ll_fnC(theta_init, hh_index, TT_vec, init_state, DataState, choice_seq, 
						DP_init, Bellman_operatorC, control_list, policy_k_guess))

pct		<- proc.time()
optsol 	<- maxLik(ll_fnC, hh_index = hh_index, TT_vec=TT_vec, init_state=init_state, 
						DataState = DataState, choice_seq = choice_seq, policy_k_guess = policy_k_guess, 
						DP_init = DP_init, Bellman_operator = Bellman_operatorC,control_list = control_list,
				start = theta_init,method ="BFGS", print.level=10, fixed = myfix)
use.time <- proc.time() - pct
cat("Estimation used", use.time[3]/60,"min.\n")
summary(optsol)

save.image(paste("run_sim/simulation_DAM_", run_id,".rdata", sep=""))
cat("\n Having finished simulation.\n")

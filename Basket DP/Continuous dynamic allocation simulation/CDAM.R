library(ggplot2)
library(reshape2)
library(scatterplot3d)
library(glmmML)
library(data.table)
library(Rcpp)
library(RcppGSL)
library(maxLik)

# setwd("~/Documents/Research/Store switching/Exercise/Basket DP/Continuous dynamic allocation simulation")
setwd("/home/brgordon/ccv103/Exercise/Basket DP")

source("policy_iteration_functions.R")
sourceCpp("0_CDAM_functions.cpp")

control_list <- list(	max_iter 		= 50,
						tol				= 1e-3,
						display_freq	= 20,
						inner_max		= 30, 
						inner_tol		= 1e-5,
						NM_sizetol		= 1e-5, 
						NM_startsize 	= .3, 
						NM_iter			= 500, 
						brent_tol 		= 1e-6, 
						brent_iter		= 200, 
						bound_eps		= 1e-4, 
						lnzero			= -30
				)
make_plot <- FALSE

args <- commandArgs(trailingOnly = TRUE)
print(args)
if(length(args)>0){
    for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
    }
}

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

#####################
# Simulation set up #
#####################
# Parameter 
param 		<- c(lambda = 4, lambda_o = 1, tau = .2, mu_z = -2, ln_sigma_z = log(.6))
beta		<- .95
n_Inc		<- 4

# Kernel of exogenous shock Z
gh_num_nodes <- 6
lnZ_weight	<- ghq(gh_num_nodes, modified = F)$weights
lnZ_nodes 	<- ghq(gh_num_nodes, modified = F)$zeros
lnZ_nodes_adj	<- sqrt(2)*exp(param["ln_sigma_z"])*lnZ_nodes + param["mu_z"]
lnZ_kernel	<- list(weight = lnZ_weight, nodes = lnZ_nodes)

# Markov probability of income transition
# Inc_weight		<- diag(n_Inc)
Inc_weight 	<- rbind(diag(n_Inc), diag(n_Inc))
Inc_weight[1,]	<- c(.75, .25,	0,	0)
Inc_weight[2,]	<- c(.15, 	.8,	.05,0)
Inc_weight[3,]	<- c(.1,		.1,	.7,	.1)
Inc_weight[4,]	<- c(0,		.1,	.2,	.7)
Inc_weight[5,]	<- c(.7, 	.1,	.1,	.1)
Inc_weight[6,]	<- c(.2, 	.7,	.05,.05)
Inc_weight[7,]	<- c(.1,	.2,	.6,	.1)
Inc_weight[8,]	<- c(.05,	.15,	.2,	.6)
Inc_nodes		<- 2*(1:n_Inc)
Inc_index		<- 1:n_Inc
Inc_kernel		<- list(weight = Inc_weight, nodes = Inc_nodes)

# Initiate states
my_grid <- list(I = c(0, .5, 1, 2.5, 4, 4.5, 5 ), 
				Inc	= Inc_index, 
				Rc	= c(0,1), 
				lnZ = log(c(.001, .05, .1, .2, .5, 1.2, 5)))	
state	<- as.matrix(my_grid[[1]], ncol=1) 
for(i in 2:length(my_grid)){
	state <- matrix_expand(state, my_grid[[i]])
}
colnames(state) <- c("I","Inc","Rc","lnZ")
summary(state)

# Omega function
omega_fn 	<- function(y){
		-.01*y^2 + .6*y
}
Qk 	<- .1
Q_fn <- function(y){
	Qk*y
}

# Derivative of the two functions 
omega_j	<- function(y){
	-0.02*y + .6
}
Q_j		<- function(y){
	rep(Qk, length(y))
}

####################
# Dynamic solution #
####################
# DP structure 
DP_list	<- list(state 	= state,
value_fn= 2 * log(state[,1] + 1), 
# 				value_fn= 0, 
				param  	= param, 
				beta 	= beta,
				omega	= omega_fn, 
				Q_fn 	= Q_fn, 
				k		= Qk, 
				omega_j	= omega_j,
				Q_j		= Q_j, 
				lnZ_kernel = lnZ_kernel,
				Inc_kernel = Inc_kernel,
				Iteration = 0, 
				status = 1
				)

sol	<- policy_iterC(DP_list, control_list, print_level = 1)
DP_list <- sol

# ------------------------------ # 
# Check DP solution 
if(make_plot){
	# Plot value function 
	sel		<- DP_list$state[,3] == 1 & DP_list$state[,2] == state[1,2]
	scatterplot3d(DP_list$state[sel,1], DP_list$state[sel,4], DP_list$value_fn[sel], xlab = "I", ylab="lnZ", zlab="value", 
				main="Value funciton", angle = 60)

	# 2D value function along I
	ggtmp <- data.frame(state, v=DP_list$value_fn)
	quartz()
	print(ggplot(ggtmp, aes(I, v, col=factor(Inc), linetype = factor(Rc)) ) + geom_point() + geom_line() + 
			facet_wrap(~ lnZ) + 
			labs(title = "Value function in I\n facet by Z")
	)
	
	# Plot policy function y and c
	ggtmp 			<- data.frame(state, y = DP_list$policy[,1], c = DP_list$policy[,2])
	ggtmp$Q			<- Q_fn(ggtmp$y)
	ggtmp$c_rat		<- with(ggtmp, c/(I+Q))
	ggtmp			<- melt(ggtmp, id.var=c("I","Inc","Rc","lnZ","Q"))
	
	quartz()
	print(ggplot(ggtmp, aes(I, value, col=factor(Inc), linetype = factor(Rc), alpha=.3 )) + geom_point() + geom_line() + 
			facet_grid(variable ~ lnZ, scales="free_y") + 
			# geom_abline(aes(intercept = Q, slope = 1), linetype=2) + 
			labs(title = paste("Policy function",sep=""))
	)
	
	# Shape of objective function along c
	c_grid 		<- seq(.5, 4, .2)
	plot_state 	<- 112:120
	state[plot_state,]
	V_spl		<- spl_init(state, DP_list$value_fn)
	tmp			<- matrix(NA, length(plot_state), length(c_grid) + 1, dimnames = list(plot_state, c(c_grid, NA)))
	for(i in 1:length(plot_state)){
		for(j in 1:length(c_grid)){
			sel		<- plot_state[i]
			sel1	<- n_Inc*state[sel,3] + state[sel,2]
			tmp[i,j] <- valuefn(y=DP_list$policy[sel,1], c=c_grid[j], I=state[sel,1], Inc_index=state[sel,2], Rc=state[sel,3], lnZ=state[sel,4], 
								V_spl, Inc_weight[sel1,], DP_list)
		}
		tmp[i, (length(c_grid)+1)] <- valuefn(y=DP_list$policy[sel,1], c=DP_list$policy[sel,2], I=state[sel,1], Inc_index=state[sel,2], 
										Rc=state[sel,3], lnZ=state[sel,4], V_spl, Inc_weight[sel1,], DP_list)
	}
	tmp_cstar	<- DP_list$policy[plot_state, 2]
	ggtmp 		<- melt(tmp)
	names(ggtmp) <- c("s_index","c","value")
	sel 		<- is.na(ggtmp$c)
	ggtmp[sel,"c"] <- tmp_cstar
	ggtmp$max_pnt <- 0
	ggtmp[sel,"max_pnt"] <- 1
	ggtmp$I		<- state[ggtmp$s_index, 1]
	ggtmp$Inc	<- state[ggtmp$s_index, 2]
	ggtmp$Rc	<- state[ggtmp$s_index, 3]
	ggtmp$lnZ	<- state[ggtmp$s_index, 4]
	quartz()
	print(ggplot(ggtmp, aes(c, value, linetype = factor(Inc), shape = factor(Rc))) + geom_point(aes(col=factor(max_pnt))) + 
				geom_line() + 
				facet_grid(I ~ lnZ, labeller = "label_both") + 
				scale_color_manual(values=c("black","red")) + 
				guides(color=FALSE)
	)
}

#################
# Simulate data #
#################
init_state 	<- state[183,]	# Initial state
TT			<- 300				# The length of simulating time series

# No randomness in the consumption function 
set.seed(666) 
my_data <- simulate_seq1_fn(init_state, TT, DP_list, control_list)

# Plot the simulation
if(make_plot){
	ggtmp 	<- my_data
	ggtmp$Q <-  DP_list$k * ggtmp$y
	ggtmp$c_inventory <- with(ggtmp, c/(I+Q))
	ggtmp1 <- melt(data.frame(ggtmp[,c("t","I","Inc","lnZ", "y", "c","c_inventory")]), id="t") 
	quartz()
	ggplot(ggtmp1,aes(t, value)) + geom_point() + geom_line() + 
			facet_grid(variable ~ ., scales="free_y") + 
			labs(title = "Simulated sequence")
}

##################
# Identification #
##################
sel_param	<- "lambda"
delta 		<- 2
param1		<- param
param1[sel_param] <- param[sel_param] + delta
DP_list1	<- DP_list
DP_list1$param	<- param1
DP_list1	<- policy_iterC(DP_list1, control_list, print_level = 1)

# Compare the value function and policy function 
tmp1		<- data.frame(state, value = DP_list$value_fn, policy = DP_list$policy, Parameter="Baseline")
tmp2		<- data.frame(state, value = DP_list1$value_fn, policy = DP_list1$policy, Parameter="Experiment")
ggtmp		<- rbind(tmp1, tmp2)
names(ggtmp)[(ncol(state)+1+1:2)] <- c("policy.y","policy.c")
selRc		<- 1
if(make_plot){
	quartz()
	print(ggplot(subset(ggtmp, Rc==selRc), aes(I, value, col=Parameter)) + 
				geom_point() + geom_line() + 
				facet_grid(Inc ~ lnZ, labeller = "label_both") + 
				labs( title = paste("Value function \nBaseline: ", paste(names(param), param, collapse=",", sep="="), 
							"\n Experiment: ", paste(names(param1), param1, collapse=",", sep="="), sep="") )
	)
	
	quartz()
	print(ggplot(subset(ggtmp, Rc==selRc), aes(I, policy.y, col=Parameter)) + 
				geom_point() + geom_line() + 
				facet_grid(Inc ~ lnZ, labeller = "label_both") + 
				labs( title = paste("Policy function y \nBaseline: ", paste(names(param), param, collapse=",", sep="="), 
							"\n Experiment: ", paste(names(param1), param1, collapse=",", sep="="), sep="") )
	)
	
	quartz()
	print(ggplot(subset(ggtmp, Rc==selRc), aes(I, policy.c, col=Parameter)) + 
				geom_point() + geom_line() + 
				facet_grid(Inc ~ lnZ, labeller = "label_both") + 
				labs( title = paste("Policy function c\nBaseline: ", paste(names(param), param, collapse=",", sep="="), 
							"\n Experiment: ", paste(names(param1), param1, collapse=",", sep="="), sep="") )
	)
}

# Simulate the sequence
TT			<- 200
lnz_draw	<- rnorm(TT)
init_state	<- state[100,]
Inc_draw	<- rep(NA, TT)
Inc_draw[1] <- init_state[2]
for(i in 2:TT){ Inc_draw[i] <- sample(Inc_index, 1, prob = Inc_weight[n_Inc*init_state[3] + Inc_draw[(i-1)],])}
sim1		<- simulate_seq1_fn(init_state, TT, DP_list, control_list, lnZ_draw=lnz_draw, Inc_draw=Inc_draw)
sim2		<- simulate_seq1_fn(init_state, TT, DP_list1, control_list, lnZ_draw=lnz_draw, Inc_draw=Inc_draw)

# Plot all the sequence
ggtmp		<- rbind(data.frame(sim1, Parameter="Baseline"), data.frame(sim2, Parameter="Experiment"))
ggtmp$Q 	<- Q_fn(ggtmp$y)
ggtmp$c_inventory <- with(ggtmp, c/(I+Q))
ggtmp$Inc_lv <- Inc_nodes[ggtmp$Inc]
if(make_plot){
	ggtmp1 <- melt(ggtmp, id=c("Parameter","t") )
	quartz()
	ggplot(ggtmp1,aes(t, value, col=Parameter, alpha = .3)) + geom_point() + geom_line() + 
			facet_grid(variable ~ ., scales="free_y") + 
			labs(title = "Simulated sequence")
}

# Summarize the resulting variables 
tmp		<- data.table(ggtmp)
tmp		<- tmp[,list(sum_c = sum(c), mean_c = mean(c), sd_c = sd(c), 
					 sum_y = sum(y), mean_y = mean(y), sd_y = sd(y), 
					 mean_I = mean(I), mean_c_rat = mean(c_inventory, na.rm=T)), 
				by=list(Parameter)]
ggtmp1	<- data.table(ggtmp)				
ggtmp1	<- ggtmp1[, ':='(y_lag = c(NA,y[-length(y)])), by=list(Parameter)]
tmp1	<- lm(y ~ I + Inc_lv + y_lag, data=subset(ggtmp1, Parameter == "Baseline"))
tmp2	<- lm(y ~ I + Inc_lv + y_lag, data=subset(ggtmp1, Parameter != "Baseline"))
tmp		<- cbind(tmp, rbind(coef(tmp1)[-1], coef(tmp2)[-1]))
tmp

##############
# Estimation #
##############
set.seed(666)
hh_index 	<- 1:6
TT_vec		<- rep(100, length(hh_index))
N			<- sum(TT_vec)
init_state	<- state[c(60, 180, 238, 278, 322, 380),]
choice_seq 	<- matrix(-999, N, 2)
DataState	<- matrix(NA, N, ncol(state))
lnz_draw	<- rnorm(N)
DataState[,4] <- lnz_draw
cnt			<- 1
for(h in 1:length(hh_index)){
	tmpRc			<- init_state[h,3]
	DataState[cnt,] <- init_state[h,]
	for(j in 2:TT_vec[h]){
		sel <- cnt + j - 1
		DataState[sel,2] <- sample(Inc_index, 1, prob = Inc_weight[n_Inc*tmpRc + DataState[(sel-1),2],])
		DataState[sel,3] <- tmpRc
	}
	cnt <- cnt + TT_vec[h]
}

# Simualate data
sim 		<- sim_hhseqC(hh_index, TT_vec, init_state, DataState, choice_seq, DP_list, control_list) 
DataState 	<- sim$DataState
choice_seq 	<- sim$choice_seq
DataState[,4] <- lnz_draw
DataState[,1] <- NA
choice_seq[,2] <- (-999)

# Try GMM objective function
n_M 		<- 4		# Number of moments; 
W			<- diag(n_M)
nz			<- 50
lnz_draw	<- matrix(rnorm(N*nz), N, nz)
# system.time(tmp <- GMM_objC(param, hh_index, TT_vec, init_state, DataState, choice_seq, W, DP_list, lnz_draw, control_list))
DP_init		<- DP_list
DP_init$value_fn <- 2*log(state[,1]+1) 
GMM_wrapper <- function(param){
	GMM_objC(param, hh_index, TT_vec, init_state, DataState, choice_seq, W, DP_init, lnz_draw, control_list)
}

#--------------------#
# GMM function shape #
my_grid	<- list(lambda = seq(-2, 2, .2), tau = seq(-.2, .3, .1), mu_z = seq(-2, 2.5, .2), ln_sigma_z = seq(-.5, 1, .2)) 
ggtmp <- data.frame(NULL)
for(i in 1:length(my_grid)){
	prc <- proc.time()
	tmp <- NULL
	sel	<- names(my_grid)[i]
	for(j in 1:length(my_grid[[i]])){
		theta_init 		<- param
		theta_init[sel] <- theta_init[sel]+my_grid[[i]][j]
		mm				<- GMM_wrapper(theta_init)
		tmp 			<- rbind(tmp, c(theta_init, mm))
	}
	ggtmp 		<- rbind(ggtmp, data.frame(Variable = sel, tmp))
	use.time	<- proc.time() - prc
	cat("Finish", sel, "with", use.time[3]/60,"min.\n")
}
colnames(ggtmp) <- c("Variable", names(param), "m")
ggtmp$par		<- NA
for(i in 1:length(my_grid)){
	sel <- ggtmp$Variable == names(my_grid)[i]
	if(length(sel) >0){
		ggtmp[sel,"par"] <- ggtmp[sel,names(my_grid)[i]]
	}
}
ggtmp$true 	<- param[as.character(ggtmp$Variable)]
mm_comp		<- ggtmp

pdf("graph_mm.pdf", width = 8, height = 8)
ggplot(ggtmp, aes(par, m)) + geom_point() + geom_line() + 
		facet_wrap(~Variable, scales="free") + 
# 		ylim(c(-200,0)) +
		geom_vline(aes(xintercept = true), col="red", linetype= 2) + 
		labs(title = "Objective function of GMM")
dev.off()

save.image(file = "sim_CDAM_mm.rdata")
if(plotmm_id){
	quit(save='n')
}

#--------------------#
# Estimation #
param_init <- param
param_init[1] <- param[1] + 1
param_init[3] <- param[3] + .1
param_init[4] <- param[4] -1
param_init[5] <- param[5] + .2

pct 	<- proc.time()
# sol 	<- optim(param_init, GMM_wrapper, method="BFGS", hessian=TRUE)
sol		<- maxLik(GMM_wrapper, start = param_init, method="BFGS", fixed=2, print.level = 10)
print(summary(sol))
save.image(file = "sim_CDAM.rdata")

use.time <- proc.time() - pct
cat("The optimization finishes with", use.time[3]/60,"min.\n")


cat("This program is done. ")







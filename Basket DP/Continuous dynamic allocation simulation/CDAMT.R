library(ggplot2)
library(reshape2)
library(scatterplot3d)
library(glmmML)
library(data.table)
library(Rcpp)
library(RcppGSL)
library(maxLik)
library(VGAM)

# setwd("~/Documents/Research/Store switching/Exercise/Basket DP/Continuous dynamic allocation simulation")
setwd("/home/brgordon/ccv103/Exercise/Basket DP")
options(error = quote({dump.frames(to.file = TRUE)}))

args <- commandArgs(trailingOnly = TRUE)
print(args)
if(length(args)>0){
    for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
    }
}

# source("policy_iteration_functions.R")
sourceCpp("1_CDAMT_functions.cpp")

control_list <- list(	max_iter 		= 30,
						tol				= 1e-3,
						display_freq	= 20,
						inner_max		= 30, 
						inner_tol		= 1e-5,
						NM_sizetol		= 1e-6, 
						NM_startsize 	= .3, 
						NM_iter			= 500,
						NM_stop			= 12, 
						brent_tol 		= 1e-6, 
						brent_iter		= 200, 
						bound_eps		= 1e-4, 
						lnzero			= -30
				)
make_plot <- FALSE

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

logit_pred <- function(DataState, d){
	sel 	<- d!=0
	mydata	<- data.frame(d = d, I=DataState[,1], Inc=factor(DataState[,2]), Rc=factor(DataState[,3]), lnZ=DataState[,4])
	myfit 	<- try(vglm(d ~ I+Inc+Rc+lnZ+I(I^2), family = "multinomial", mydata[sel,]), TRUE )
	if(class(myfit)!="try-error"){
		phat 	<- try(predict(myfit), TRUE)
	}
	if(class(myfit)=="try-error" | class(phat)=="try-error"){
		p	<- rep(1, length(sel)) %*% t(as.vector(table(d[sel])/length(sel)) )
		phat<- apply(p[sel,-ncol(p)], 2, function(x) x/p[sel,ncol(p)])
	}
	out 	<- matrix(0, length(d), ncol(phat))
	out[sel,] <- phat
	return(out)
}

#####################
# Simulation set up #
#####################
# Parameter 
param 		<- c(lambda = 4, tau1 = .1, tau2 = .05, tau3 = .1, tau4 = .5, lambda_o = 1)
beta		<- .95
n_Inc		<- 4
mu_z		<- -2
sigma_z		<- .6

# Kernel of exogenous shock Z
gh_num_nodes <- 6
lnZ_weight	<- ghq(gh_num_nodes, modified = F)$weights
lnZ_nodes 	<- ghq(gh_num_nodes, modified = F)$zeros
lnZ_nodes_adj	<- sqrt(2)*sigma_z*lnZ_nodes + mu_z
lnZ_kernel	<- list(weight = lnZ_weight, nodes = lnZ_nodes_adj)

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
omega_fn 	<- function(y, d){
		-.015*d*y^2 + .2*d*y
}
Qk 	<- .1
Q_fn <- function(y){
	Qk*y
}

# Derivative of the two functions 
omega_j	<- function(y,d){
	-0.03*d*y + .15*d
}
Q_j		<- function(y){
	rep(Qk, length(y))
}

tmp <- seq(0, 6, .1)
ggtmp <- data.frame(NULL)
for(i in 1:3){
	ggtmp	<- rbind(ggtmp, data.frame(x=tmp, y=omega_fn(tmp, i), d=i))
}
if(make_plot){
	ggplot(ggtmp, aes(x,y,col=factor(d))) + geom_point() + geom_line()
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
				status = 1, 
				D		= 3
				)

system.time(sol	<- policy_iterC(DP_list, control_list, print_level = 1))
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
	ggtmp 			<- data.frame(state, y = DP_list$policy[,2], c = DP_list$policy[,3])
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
	c_grid 		<- seq(.1, 4, .1)
	plot_state 	<- which(state[,1]<2 & state[,2] ==3 & state[,4]==unique(state[,4])[1])
	state[plot_state,]
	tmp			<- plotvalueC(plot_state, c_grid, DP_list, control_list)
	tmp[tmp==0] <- NA
	dimnames(tmp) <- list(plot_state, c(c_grid, NA))
	tmp_cstar	<- DP_list$policy[plot_state, 3]
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
sel 		<- state[,1] == 4 & state[,2]==3 & state[,3] == 1 & state[,4]== unique(state[,4])[1]
init_state 	<- matrix(state[sel,], nrow=1)	# Initial state
TT			<- as.vector(300)				# The length of simulating time series
choice_seq	<- matrix(-999, TT, 3)
DataState	<- matrix(NA, TT, 4)
DataState[,3]<- init_state[3]
DataState[,4]<- rnorm(TT)
tmpRc			<- init_state[3]
DataState[1,] <- init_state
for(j in 2:TT){
	DataState[j,2] <- sample(Inc_index, 1, prob = Inc_weight[n_Inc*tmpRc + DataState[(j-1),2],])
}

# No randomness in the consumption function 
set.seed(666) 
tmp <- sim_hhseqC(1, TT, init_state, DataState, choice_seq, DP_list, control_list, simidx=3)
my_data	<- data.frame(t = 1:TT, tmp$DataState, tmp$choice_seq)
names(my_data)	<- c("t", "I", "Inc", "Rc", "lnZ", "d", "y", "c")

# Plot the simulation
if(make_plot){
	ggtmp 	<- my_data
	ggtmp$Q <-  Q_fn(ggtmp$y)
	ggtmp$c_inventory <- with(ggtmp, c/(I+Q))
	ggtmp1 <- melt(data.frame(ggtmp[,c("t","I","Inc","lnZ", "d","y", "c","c_inventory")]), id="t") 
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
names(ggtmp)[(ncol(state)+1+1:3)] <- c("policy_d", "policy.y","policy.c")
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

##############
# Estimation #
##############
set.seed(666)
hh_index 	<- 1:10
TT_vec		<- rep(100, length(hh_index))
N			<- sum(TT_vec)
init_state	<- state[c(60, 84, 135, 191, 201, 221, 238, 278, 322, 380),]
choice_seq 	<- matrix(-999, N, 3)
DataState	<- matrix(NA, N, ncol(state))
lnz_draw	<- rnorm(N)*sigma_z + mu_z
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
sim 		<- sim_hhseqC(hh_index, TT_vec, init_state, DataState, choice_seq, DP_list, control_list, simidx=3) 
DataState 	<- sim$DataState
choice_seq 	<- sim$choice_seq
DataState[,4] <- lnz_draw
DataState[,1] <- NA
choice_seq[,3] <- (-999)
table(choice_seq[,1])

# Try GMM objective function
n_M 		<- 4		# Number of moments; 
W			<- diag(n_M)
nz			<- 50
lnz_draw	<- matrix(rnorm(N*nz, mu_z, sigma_z), N, nz)
DP_init		<- DP_list
DP_init$value_fn <- 2*log(state[,1]+1) 
Q			<- Q_fn(choice_seq[,2])
omega_mat	<- sapply(1:DP_list$D, function(i) omega_fn(choice_seq[,2], i))
GMM_wrapper <- function(param){
	GMM_objC(param, hh_index, TT_vec, init_state, DataState, choice_seq, W, DP_init, lnz_draw, control_list, logit_pred, Q, omega_mat)
}

#--------------------#
# GMM function shape #
if(run_id=="mm"){
	my_grid	<- list(lambda = seq(-2, 2, .2), tau1 = seq(-.1, .4, .05), tau2 = seq(-.05, .3, .05), tau3 = seq(-.1, .5, .05), 
					tau4 = seq(-.4, .5, .1), lambda_o = seq(-2, 2, .2)) 
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
	print(ggplot(ggtmp, aes(par, m)) + geom_point() + geom_line() + 
			facet_wrap(~Variable, scales="free") + 
	# 		ylim(c(-200,0)) +
			geom_vline(aes(xintercept = true), col="red", linetype= 2) + 
			labs(title = "Objective function of GMM"))
	dev.off()
	
	save.image(file = paste("sim_CDAMT_", run_id,".rdata", sep=""))
	cat("This program is done. ")
}

if(run_id == "est"){
	#--------------------#
	# Estimation #
	param_init <- param
	param_init[1] <- param[1] - 1
	param_init[2] <- param[2] + .05
	param_init[3] <- param[3] + .1
	param_init[4] <- param[4] + .8

	pct 	<- proc.time()
	sol		<- maxLik(GMM_wrapper, start = param_init, method="BFGSR", print.level = 10)
	print(summary(sol))

	save.image(file = "sim_CDAMT.rdata")
	use.time <- proc.time() - pct
	cat("The optimization finishes with", use.time[3]/60,"min.\n")

	cat("This program is done. ")
}



library(ggplot2)
library(reshape2)
library(scatterplot3d)
library(Rcpp)
library(RcppGSL)
library(plm)
library(zoo)
library(data.table)
library(VGAM)

options(error = quote({dump.frames(to.file = TRUE)}))
args <- commandArgs(trailingOnly = TRUE)
print(args)
if(length(args)>0){
    for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
    }
}
# vid_save 	<- vid
# run_id_save <- run_id
vid_save	<- vid	<- 2
run_id_save	<- run_id <- 5
seg_id 		<- as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))

# setwd("~/Documents/Research/Store switching/processed data/Estimation")
# sourceCpp("~/Documents/Research/Store switching/Exercise/main/1_CDAMT_functions.cpp")
# setwd("/home/brgordon/ccv103/Exercise/run")
# setwd("/kellogg/users/marketing/2661703/Exercise/run")
setwd("/sscc/home/c/ccv103/Exercise/run")

sourceCpp("1_CDAMT_functions.cpp")
model_name 	<- "MDCEV_a1b1"
make_plot	<- FALSE

param_mat		<- matrix(0, 6, 12, dimnames = list(c("lambda","lambda_o","tau1","tau2","tau3","tau4"), NULL))
param_mat[,1]	<- c(5,		1,		-3,		-1.2,	2,		-1)
param_mat[,2]	<- c(2.25,	.5,		-1,		-.45,	.5,		-.6)
param_mat[,3]	<- c(1.6,	.3,		-.8,	-.3,	.5,		-.3)
param_mat[,4]	<- c(2,		.3,		-.6,	-.3,	.5,		-.3)
param_mat[,5]	<- c(2,		.5,		-.7,	-.5,	.5,		-.75)
param_mat[,6]	<- c(5,		1,		-1.2,	-.8,	1.2,	-1)
param_mat[,7]	<- c(1.2,	.2,		-.2,	-.16,	.24,	-.1)
param_mat[,8]	<- c(2,		.3,		-.6,	-.26,	.3,		-.3)
param_mat[,9]	<- c(1,		.2,		-.24,	-.16,	.16,	-.3)
param_mat[,10]	<- c(2,		.5,		-.75,	-.4,	.75,	-.5)
param_mat[,11]	<- c(2,		.5,		-1,		-.4,	.75,	-.5)
param_mat[,12]	<- c(5,		1,		-2,		-.8,	1.3,	-1.2)

###################
# Basic functions #
###################
matrix_expand <- function(mat, vec){
	n1	<- nrow(mat)
	n2	<- length(vec)
	sel	<- rep(1:n1, n2)
	mat_new <- cbind(mat[sel,], rep(vec, each=n1))
	return(mat_new)
}

my_forward 	<- function(x){
	return(c(x[-1], NA))
}

recode <- function(old.idx, old.value, new.n){
	old.level 	<- sort(unique(old.idx))
	old.n		<- length(old.level)
	tab			<- sort(table(old.idx), decreasing=T)
	sort.x		<- as.numeric(names(tab))
	
	loop.n 		<- new.n 
	loop.x		<- old.idx
	loop.sx		<- sort.x
	mycut		<- NULL
	for(i in new.n:1){
		qt 		<- quantile(loop.x, c(0:loop.n)/loop.n)
		if(length(unique(qt))==loop.n+1){
			mycut <- c(mycut, qt)
			break
		}
		mycut	<- c(mycut, loop.sx[1])
		loop.x	<- loop.x[loop.x!=loop.sx[1]]
		loop.n	<- loop.n - 1
		loop.sx <- loop.sx[-1]
	}
	new.idx 	<- cut(old.idx, mycut, include.lowest=T, labels=1:new.n)
	new.level	<- tapply(old.value, new.idx, mean, na.rm=T)
	new.value	<- new.level[new.idx]
	cat("Table of original categorical variable vs recoded categorical variable:\n")
	print(table(old.idx, new.idx))
	return(cbind(new.idx, new.value) )
} 

mytransition <- function(x1, x2){
	if(class(x1)!="factor" | class(x2) != "factor"){
		x1 	<- factor(x1)
		x2	<- factor(x2)
	}
	onem	<- diag(length(levels(x1)))
	out		<- table(x1, x2)
	sel		<- which(apply(out, 1, function(x) all(x==0)))
	if(length(sel)>0) { out[sel,] <- onem[sel,] }
	return(out/rowSums(out))
}
				
######################################
# Read data and extract segment data # 
######################################
load(paste("run_",run_id,"/5_est_MDCEV_seg",seg_id,".rdata",sep=""))
rm(list = c("run_id", "vid"))
vid					<- vid_save
run_id				<- run_id_save

# Recode income levels 
sum(is.na(hh_exp$income_real))
hh_exp		<- subset(hh_exp, !is.na(income_real))
n_Inc		<- 8
sel			<-  hh_exp$income_real >= 27
hh_exp[sel,"income_real"] <- 27
tmp			<- recode(hh_exp$income_real, hh_exp$income_midvalue, n_Inc)
hh_exp$income_nodes <- tmp[,1]
Inc_nodes	<- sort(unique(tmp[,2]))/24				# biweekly income
cat("Inc_nodes =", Inc_nodes,"\n")

# # Extract segment
mydata		<- subset(hh_exp, segment == seg_id)
a			<- with(mydata, paste(household_code, biweek, sep="*"))
sel			<- duplicated(a)
cat("Sum(duplicated hh*biweek) =", sum(sel), "\n")
ord			<- order(mydata$household_code, mydata$year, mydata$biweek)
mydata		<- mydata[ord,]

###############################################
# Extract household initial states, DV and IV #
###############################################
# -------------------------------------------# 
# Fill in the no-trip weeks
tmp			<- data.table(mydata)
tmp			<- tmp[,list(minw=min(biweek), maxw=max(biweek)), by=list(household_code)]
tmp			<- tmp[,list(biweek=minw:maxw), by=list(household_code)]
cat("Before filling in non-shoppiing periods, dim(mydata) =", dim(mydata),"\n") 
tmpb		<- with(mydata, paste(household_code, biweek, sep="*"))
mydata		<- merge(data.frame(tmp), mydata, by=c("household_code","biweek"), all.x=T)
cat("After filling in non-shopping periods, dim(mydata)=", dim(mydata),"\n")
tmpa		<- with(mydata, paste(household_code, biweek, sep="*"))
sel			<- !tmpa %in% tmpb
cat("Sum of non-shopping biweeks =", sum(sel), "\n")

# Expenditure duing noshoping trips is 0
mydata[sel,"dol_purchases"] <- 0
mydata[sel,"d"] 			<- 0
sum(is.na(mydata$d))
mydata[is.na(mydata$dol_purchases),"dol_purchases"] <- 0
mydata[is.na(mydata$d),"d"] 						<- 0

# Fill in the missing income and recession 
Rc_t		<- min(hh_exp[hh_exp$recession==1, "biweek"])
tmp			<- min(hh_exp$biweek):max(hh_exp$biweek)
bdate		<- as.Date("2004-01-01", format="%Y-%m-%d") + 14*(tmp-1)		
byear		<- year(bdate)
tmp			<- byear[mydata[sel,"biweek"]]
mydata[sel,"recession"]	<- ifelse(mydata[sel,"biweek"] >= Rc_t, 1, 0)
mydata[sel,"year"]		<- tmp
ord			<- order(mydata$household_code, mydata$year, mydata$biweek)
mydata		<- mydata[ord,]

sel1 		<- is.na(mydata$income_nodes)
tmp1 		<- data.table(mydata[!sel1,])
tmp1		<- tmp1[,list(income_nodes=unique(income_nodes), income_real=unique(income_real)), by=list(household_code, year)]
tmp2		<- merge(mydata[sel1,c("household_code","year")], tmp1, by=c("household_code","year"),all.x=T)
ord			<- order(tmp2$household_code, tmp2$year)
tmp2		<- tmp2[ord,]
mydata[sel1,"income_real"]	<- tmp2$income_real
mydata[sel1,"income_nodes"]	<- tmp2$income_nodes
mydata$income_nodes 		<- factor(mydata$income_nodes)

choice_seq	<- as.matrix(cbind(d=mydata$d, y=mydata$dol_purchases, c=NA))
DataState 	<- as.matrix(cbind(I = rep(NA, nrow(mydata)), Inc = as.numeric(mydata$income_nodes), 
							   Rc=as.numeric(mydata$recession)))
cat("Summary of choice variables:\n"); print(summary(choice_seq)); cat("\n")
cat("Summary of state variables:\n"); print(summary(DataState)); cat("\n")

# -------------------------------------------# 
# Household index and initial stata #
myidx		<- data.table(mydata)
myidx		<- myidx[,list(minw=min(biweek), maxw=max(biweek), stay=max(biweek)-min(biweek)+1, avg_Q = mean(Q,na.rm=T), 
					first_Inc = as.numeric(income_nodes[1]), first_Rc = recession[1]), 
					by=list(household_code)]
myidx		<- myidx[,id:=1:nrow(myidx)]
cat("Check cross-household distribution of stay patten:\n"); summary(myidx); cat("\n")

hh_index 	<- as.vector(myidx$id)
TT_vec		<- as.vector(myidx$stay)
init_state	<- cbind(as.matrix(myidx[,list(avg_Q, first_Inc, first_Rc)]))
cat("Initial state:\n"); summary(init_state); cat("\n")

# Check data consistency
cat("dim(mydata) =", dim(mydata),";")
cat("sum(TT_vec) =", sum(TT_vec),";")
cat("dim(DataState) =", dim(DataState),";")
cat("dim(choice_seq) =",dim(choice_seq),"\n")

########################################################################################################
# Compute DP ingredients: state transition, purchase inclusive value, Quantity function and prediction #
########################################################################################################							
# -------------------------------------------# 
# Model the income transition
# Lag income by 1 year
tmp		<- data.table(mydata)
tmp		<- tmp[, list(income_nodes = unique(income_nodes)), by=list(household_code, recession, year)]
tmp		<- tmp[, income_forward := my_forward(income_nodes), by=list(household_code)]
tmp$income_forward <- factor(tmp$income_forward, levels=1:length(levels(mydata$income_nodes)))

# Transition table: replace 0 diagal with 1 if nobody belongs to some income level.
sel 	<- tmp$recession == 1
Inc_weight <- rbind(mytransition(tmp[!sel,income_nodes], tmp[!sel,income_forward] ), 
				 mytransition(tmp[sel,income_nodes], tmp[sel,income_forward]))
cat("Transition probability of income process:\n"); print(Inc_weight); cat("\n")
Inc_index	<- 1:n_Inc
Inc_kernel	<- list(weight = Inc_weight, nodes = Inc_nodes)

# -------------------------------------------# 
# Inclusive value of purchase utility: omega(y,d)
# omega_dat 	<- melt(omega_nodraw)
omega_dat 	<- omega_draw_mean
omega_list	<- lapply(1:D, function(i) 
			splinefun(omega_dat[omega_dat$d==i,"y"], omega_dat[omega_dat$d==i,"omega"], method="natural"))
omega_fn 	<- function(y,d){
	if(d==0){
		out	<- rep(0, length(y)) 
	}else{
		out <- omega_list[[d]](y) 
	}
	return(as.vector(out))
}

omega_list1	<- lapply(1:D, function(i) 
			splinefun(omega_dat[omega_dat$d==i,"y"], omega_dat[omega_dat$d==i,"omega"], method="natural", derive=1))
domega_fn 	<- function(y,d){
	if(d==0){
		out	<- rep(0, length(y)) 
	}else{
		out <- omega_list[[d]](y, deriv = 1) 
	}
	return(as.vector(out))
}

if(make_plot){
	tmp 		<- seq(10, 1000, by=10)
	ggtmp		<- sapply(1:D, function(i) omega_fn(tmp,i))
	dimnames(ggtmp) <- list(y=tmp, d=1:D)
	ggtmp		<- melt(ggtmp)
	ggplot(ggtmp, aes(y, value, col=factor(d))) + geom_point() + geom_line() 
}

# -------------------------------------------# 
# Quantity function: Q(y)
tmp			<- hh_exp
tmp[is.na(tmp$dol_purchases),"dol_purchases"] <- 0
tmp[is.na(tmp$Q), "Q"]	<- 0
tmp$ly		<- with(tmp, log(dol_purchases+1))
Q_fit 		<- lm(Q ~ ly - 1, data=tmp)
# Q_fit		<- lm(Q ~ ly +factor(d):ly -1, data=tmp)
summary(Q_fit)
Q_coef		<- coef(Q_fit)
Q_fn		<- function(y){
	ly		<- log(y+1)
	out		<- ly * Q_coef
	return(as.vector(out))
}

dQ_fn		<- function(y){
	out		<- Q_coef /(y + 1)
	return(as.vector(out))
}

if(make_plot){
	tmp		<- seq(0, 1000, 10)
	ggtmp	<- data.frame(y = tmp, Q = Q_fn(tmp))
	ggplot(ggtmp, aes(y, Q)) + geom_point() + geom_line()
}

# -------------------------------------------# 
# Multinomial logit prediction of d: P(d|state)
logit_pred <- function(DataState, d){
	sel 	<- d!=0
	mydata	<- data.frame(d = d, I=DataState[,1], Inc=factor(DataState[,2]), Rc=factor(DataState[,3]))
	myfit 	<- try(vglm(d ~ I+Inc+Rc+I(I^2)+Inc:Rc, family = "multinomial", mydata[sel,]), TRUE )
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

###########################
# Euler linear regression # 
###########################
beta	<- .96
tmp			<- mydata
tmp$domega	<- 0
for(d in 1:D){
	sel 	<- tmp$d == d
	tmp[sel,"domega"] <- domega_fn(tmp[sel,"dol_purchases"], d)
}
tmp$dQ		<- dQ_fn(tmp$dol_purchases)
tmp$Q		<- Q_fn(tmp$dol_purchases)
tmp			<- merge(tmp, myidx[,list(household_code,avg_Q)], by = "household_code", all.x=T)
tmp$Inext	<- with(tmp, .1*(avg_Q + Q))
tmp			<- data.table(tmp)
setkeyv(tmp, c("household_code","biweek"))
tmp			<- tmp[,list(domega= domega, domega_next = my_forward(domega), dQ = dQ, dQ_next = my_forward(dQ), 
						Inext = Inext, Inext_next = my_forward(Inext), d = d, d_next = my_forward(d)),
					by = list(household_code)]
tmp			<- tmp[,':='(dv = .1*(beta*dQ/dQ_next-1), x_lambdao = domega - beta*domega_next*dQ/dQ_next, 
						x_tau1 = dQ, x_tau2 = 2*Inext*dQ, x_tau3 = d*dQ)]
summary(lm(dv ~ x_lambdao + x_tau1 + x_tau2 + x_tau3 -1, data = tmp))

#######################
# Set up DP structure #
#######################
# ------------------------------ # 
# State variables 
summary(mydata$Q);  quantile(mydata$Q, c(.01, .025, .975, .99), na.rm=T)

# Initiate states
my_grid <- list(I = c(0, .5, 1, 2, 3, 4, 4.5, 5 ), 
				Inc	= Inc_index, 
				Rc	= c(0,1))
state	<- as.matrix(my_grid[[1]], ncol=1) 
for(i in 2:length(my_grid)){
	state <- matrix_expand(state, my_grid[[i]])
}
colnames(state) <- c("I","Inc","Rc")
cat("Grids of states:\n"); summary(state); cat("\n")

# ------------------------------ # 
# Find DP solution at the intial parameters;
beta	<- .96
param	<- param_mat[,seg_id]
cat("Inital parameters at seg =", seg_id, "are:\n"); print(param); cat("\n")

DP_list	<- list(state 	= state,
				value_fn= 50 * log(state[,1] + 1) , 
				param  	= param, 
				beta 	= beta,
				omega	= list(omega_y = as.numeric(colnames(omega_nodraw)), omega = omega_nodraw), 
				Q_coef	= Q_coef,
				Inc_kernel = Inc_kernel,
				Iteration = 0, 
				status = 1, 
				D		= D
				)
control_list <- list(	max_iter 		= 30,
						tol				= 1e-4,
						display_freq	= 20,
						inner_max		= 50, 
						inner_tol		= 1e-5,
						NM_sizetol		= 1e-6, 
						NM_startsize 	= .5, 
						NM_iter			= 500, 
						NM_stop			= 15,
						brent_tol 		= 1e-6, 
						brent_iter		= 200, 
						bound_eps		= 1e-4, 
						lnzero			= -30
				)				
system.time(tmpsol	<- policy_iterC(DP_list, control_list, print_level = 1))

# ------------------------------ # 
# Check DP solution 
if(make_plot){
	# 2D value function along I
	ggtmp <- data.frame(state, v=tmpsol$value_fn)
	quartz()
	print(ggplot(ggtmp, aes(I, v, linetype = factor(Rc)) ) + geom_point() + geom_line() + 
			facet_wrap(~ Inc) + 
			labs(title = "Value function in I\n facet by Inc")
	)
	
	# Plot policy function y and c
	ggtmp 			<- data.frame(state,d =tmpsol$policy[,1], y = tmpsol$policy[,2], c = tmpsol$policy[,3])
	ggtmp$Q			<- Q_fn(ggtmp$y)
	ggtmp$c_rat		<- with(ggtmp, c/(I+Q))
	ggtmp			<- melt(ggtmp, id.var=c("I","Inc","Rc"))
	
	quartz()
	print(ggplot(ggtmp, aes(I, value, linetype = factor(Rc), alpha=.3 )) + geom_point() + geom_line() + 
			facet_grid(variable ~ Inc, scales="free_y") + 
			# geom_abline(aes(intercept = Q, slope = 1), linetype=2) + 
			labs(title = paste("Policy function",sep=""))
	)
	
	# Plot simulation process
	sel 		<- state[,1] == 4 & state[,2]==3 & state[,3] == 1 
	tmp_init 	<- matrix(state[sel,], nrow=1)	# Initial state
	tmp_TT			<- as.vector(300)				# The length of simulating time series
	tmp_choice	<- matrix(-999, tmp_TT, 3)
	tmp.DataState	<- matrix(NA, tmp_TT, 3)
	tmp.DataState[,3]<- tmp_init[3]
	tmpRc			<- tmp_init[3]
	tmp.DataState[1,] <- tmp_init
	for(j in 2:tmp_TT){
		tmp.DataState[j,2] <- sample(Inc_index, 1, prob = Inc_weight[n_Inc*tmpRc + tmp.DataState[(j-1),2],])
	}

	# No randomness in the consumption function 
	set.seed(666) 
	tmp <- sim_hhseqC(1, tmp_TT, tmp_init, tmp.DataState, tmp_choice, tmpsol, control_list, 0, simidx=3, logit_random=TRUE)
	tmp_data	<- data.frame(t = 1:tmp_TT, tmp$DataState, tmp$choice_seq)
	names(tmp_data)	<- c("t", "I", "Inc", "Rc", "d", "y", "c")

	# Plot the simulation
	ggtmp 	<- tmp_data
	ggtmp$Q <-  Q_fn(ggtmp$y)
	ggtmp$c_inventory <- with(ggtmp, c/(I+Q))
	ggtmp1 <- melt(data.frame(ggtmp[,c("t","I","Inc","d","y", "c","c_inventory")]), id="t") 
	quartz()
	ggplot(ggtmp1,aes(t, value)) + geom_point() + geom_line() + 
			facet_grid(variable ~ ., scales="free_y") + 
			labs(title = "Simulated sequence")
}

##############
# Estimation #
##############
set.seed(666)
n_M 		<- 2 * (D+1)		# Number of moments; 
W			<- diag(n_M)/100
DP_init		<- DP_list
DP_init$value_fn <- 20*log(5*state[,1]+1) 
Q_obs		<- Q_fn(choice_seq[,2])
omega_mat	<- sapply(1:D, function(i) omega_fn(choice_seq[,2], i))
omega_dif	<- apply(omega_mat[,-D], 2, function(x) x - omega_mat[,D])

# Pre-compute the log odds
tmpdata		<- data.table(mydata)
tmpdata		<- tmpdata[, lagy:= c(dol_purchases[1], dol_purchases[-length(dol_purchases)]), by = list(household_code)]
tmpdata		<- tmpdata[,':='(Inc = factor(income_nodes), Rc = factor(recession))]
sel			<- tmpdata$d > 0
lgm			<- try(vglm(d ~ lagy + I(lagy^2) + Inc + Rc + lagy:Rc, data=tmpdata[sel,], family = "multinomial"), TRUE)
if(class(lgm)=="try-error"){
	lgm		<- try(vglm(d ~ lagy + I(lagy^2) + Inc + Rc, data=tmpdata[sel,], family = "multinomial"), TRUE)
}
cat("Summary of multinomial logit of shopping trips:\n"); summary(lgm); cat("\n")
phat		<- matrix(0, nrow(choice_seq), D-1)
phat[sel,]	<- predict(lgm)

# A wrapper function of GMM
par_neg_idx	<- c(3,4)
nest_pred	<- FALSE
GMM_wrapper <- function(param){
	GMM_objC(param, hh_index, TT_vec, init_state, DataState, choice_seq, W, DP_init, control_list, 
			logit_pred, Q_obs, omega_dif, phat, nest_pred = nest_pred)
}
param_init 	<- param
system.time(print(GMM_wrapper(param_init)))

n_sim	<- 500
np		<- length(DP_list$param)
udraw	<- matrix(runif(np*n_sim), np, n_sim, dimnames = list(c("lambda", "lambda_o","tau1","tau2","tau3", "tau4"),1:n_sim))
bound_mat<- matrix(c(	5, 35, 
						.01,  2,
						-4, 2,
					  	-4, -.05,
						.1, 3,
						-6, -.1
					  ), np, 2,
			byrow=T, dimnames = list(c("lambda", "lambda_o", "tau1","tau2","tau3", "tau4"), c("min","max")))
par_draw <- bound_mat[,1] + (bound_mat[,2] - bound_mat[,1]) * udraw

cat("Summary of parameter draws:\n")
summary(t(par_draw)); cat("\n")

# Compute moments over each parameter draws
par_dat 	<- matrix(NA, ncol(par_draw), 21, dimnames = list(1:ncol(par_draw), 
					c("lambda", "lambda_o","tau1","tau2","tau3", "tau4", "status","value","obj",paste("M", 1:n_M, sep=""))))
for(i in 1:ncol(par_draw)){
	cat("Parameter draw", i,"\n"); print(par_draw[,i]); cat("\n")
	tmp		<- GMM_wrapper(par_draw[,i])
	par_dat[i,]	<- c(par_draw[,i], status = tmp$status, value = mean(tmp$value), obj=tmp$obj, tmp$Mbar)
	# cat("obs =", tmp$obj,"\n")
	# cat("Moments are", tmp$Mbar, "\n")
	# cat("-----------------------------------------------------------\n")
}

cat("sum(is.na(obj))=", sum(is.na(par_dat[,"obj"])),"\n")
par_dat	<- par_dat[order(par_dat[,"obj"]),]
sel		<- which(!is.na(par_dat[,"obj"]))[1:5]
cat("The parameters that minimize the objective function are:\n"); 
print(par_dat[sel,]); cat("\n")
cat("The status of DP solutions:\n"); print(table(par_dat[,"status"])); cat("\n")

#-----------------------------------------------------------------#
sel		<- which.min(par_dat[,"obj"])
par_fix <- par_dat[sel,1:6]
names(par_fix)	<- c("lambda", "lambda_o","tau1","tau2","tau3", "tau4")
mygrid 			<- list(lambda = seq(-1, 2, .2), lambda_o = seq(-.1, .2, .05), tau1 = seq(-1, 2, .2), tau2 = seq(-1, 2, .2), 
						tau3 = seq(-1, 1, .2), tau4 = seq(-1, 1, .2))
par_dat1		<- data.frame()
for(i in 1:length(par_fix)){
	for(j in 1:length(mygrid[[i]])){
		tmp						<- par_fix
		tmp[names(mygrid)[i]]	<- tmp[names(mygrid)[i]] + mygrid[[i]][j]
		DP_list$param			<- tmp
		tmpsol					<- policy_iterC(DP_list, control_list)
		tmpsol1					<- GMM_wrapper(tmp)
		par_dat1				<- rbind(par_dat1, data.frame(adjV = names(mygrid)[i], adj_value = tmp[i], 
										 obj = tmpsol1$obj, V = median(tmpsol$value_fn), d = median(tmpsol$policy[,1]), 
										 y = median(tmpsol$policy[,2])))
	}
	cat("Finish", names(mygrid)[i],".\n")
}
sel 			<- which.min(par_dat1[,"obj"])
tmp						<- par_fix
tmp[par_dat1[sel,"adjV"]]	<- tmp[par_dat1[sel,"adjV"]] + par_dat1[sel,"adj_value"]
cat("Parameters at lowest obj are:\n"); print(tmp); cat("\n")

# Plot comparative statics
ggtmp	<- melt(par_dat1[,c("adjV","adj_value","obj", "y")], id.vars=c("adjV","adj_value"))
pdf(paste("graph_comp_seg",seg_id,".pdf",sep=""), width = 8, height = 8)
print(ggplot(ggtmp, aes(adj_value, value)) + geom_point() + geom_line() + 
		facet_wrap(variable ~ adjV, scales="free") + 
		labs(title = paste(paste(names(par_fix), "=", round(par_fix,2), sep=""), collapse=","))
	) 
# Plot policy function y and c
sel 			<- which.min(par_dat1[,"obj"])
tmp						<- par_fix
tmp[par_dat1[sel,"adjV"]]	<- tmp[par_dat1[sel,"adjV"]] + par_dat1[sel,"adj_value"]
DP_list$param	<- tmp
tmpsol			<- policy_iterC(DP_list, control_list)
ggtmp 			<- data.frame(state,d =tmpsol$policy[,1], y = tmpsol$policy[,2], c = tmpsol$policy[,3])
ggtmp$Q			<- Q_fn(ggtmp$y)
ggtmp$c_rat		<- with(ggtmp, c/(I+Q))
ggtmp			<- melt(ggtmp, id.var=c("I","Inc","Rc"))
print(ggplot(ggtmp, aes(I, value, linetype = factor(Rc), alpha=.3 )) + geom_point() + geom_line() + 
		facet_grid(variable ~ Inc, scales="free_y") + 
		# geom_abline(aes(intercept = Q, slope = 1), linetype=2) + 
		labs(title = paste("Policy function",sep=""))
)
dev.off()

save.image(file = paste("run_",run_id,"/test_init_seg",seg_id,"_v",vid,".rdata",sep=""))

cat("File has been saved.\n")
cat("This program is done.\n")

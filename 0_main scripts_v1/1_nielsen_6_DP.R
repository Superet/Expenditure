library(ggplot2)
library(reshape2)
library(scatterplot3d)
library(maxLik)
library(Rcpp)
library(RcppGSL)
library(glmmML)
library(plm)
library(zoo)
library(data.table)
library(doParallel)

# setwd("~/Documents/Research/Store switching/processed data")
# setwd("//tsclient/Resear1/Store switching/processed data")
setwd("/home/brgordon/ccv103/processed data")

sourceCpp("../Exercise/Basket DP/0_DAM_functions.cpp")
model_name 	<- "MDCEV_a1b1"
run_id		<- 1
seg_id		<- 1
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
				
######################################
# Read data and extract segment data # 
######################################
load(paste("run_",run_id,"/result_5_",model_name,".rdata",sep=""))

# Segment households 
hh_exp$income_group <- factor(hh_exp$income_group, levels=paste("Qt",1:4,sep=""))
hh_exp$famsize		<- factor(hh_exp$famsize, levels = c("Single","Two", "Three+"))
hh_exp$income_real	<- factor(hh_exp$income_real)
panelist	<- data.table(hh_exp)
setkeyv(panelist, c("household_code","year","month"))
panelist	<- panelist[,list(first_income = income_group[1], first_famsize = famsize[1]), by=list(household_code)]
table(panelist$first_income, panelist$first_famsize)
seg_index	<- 1:(length(levels(panelist$first_income))*length(levels(panelist$first_famsize)))
names(seg_index) <- paste(rep(levels(panelist$first_income), length(levels(panelist$first_famsize))), 
						  rep(levels(panelist$first_famsize), each=length(levels(panelist$first_income)) ), sep="-")
tmp			<- paste(panelist$first_income, panelist$first_famsize, sep="-")
panelist$segment <- seg_index[tmp]
sel			<- panelist$segment == seg_id
mydata		<- subset(hh_exp, household_code %in% panelist[sel,household_code])
ord			<- order(mydata$household_code, mydata$year, mydata$month)
mydata		<- mydata[ord,]
mydata$Q	<- mydata$food_quant + mydata$nonedible_quant
mydata$ymonth <- paste(mydata$year, mydata$month,sep="-")

####################
# Organize DP data # 
####################
# -------------------------------------------# 
# Household index and initial stata #
myidx		<- data.table(mydata)
myidx		<- myidx[,list(stay = length(month), avg_Q = mean(Q)/my_scale, firstyear = min(year), 
					first_lny = log(dol_purchases[1]), first_Inc = as.numeric(income_group[1]), first_Rc = recession[1]), 
					by=list(household_code)]
myidx		<- myidx[,id:=1:nrow(myidx)]

hh_index 	<- as.vector(myidx$id)
TT_vec		<- as.vector(myidx$stay)
init_state	<- as.matrix(myidx[,list(avg_Q, first_lny, first_Inc, first_Rc)])
choice_seq	<- as.vector(mydata$basket_type)
DataState 	<- as.matrix(cbind(I = rep(NA, nrow(mydata)), lny = log(mydata$dol_purchases), 
							Inc = as.numeric(mydata$income_group), Rc=as.numeric(mydata$recession)))

# -------------------------------------------# 
# Model the expenditure process # 
mydata	<- data.table(mydata)
mydata	<- mydata[, lag_dol_purchases := my_lag(dol_purchases), by = list(household_code)]
mydata	<- mydata[,':='(ln_dol_purchases = log(dol_purchases), ln_lag_dol_purchases = log(lag_dol_purchases))]
mydata	<- data.frame(mydata)

sel		<- mydata$recession==1
fit0	<- plm(ln_dol_purchases ~ ln_lag_dol_purchases + factor(income_group) , data=mydata[!sel,], 
				index = c("household_code","ymonth"), model = "within")
rho0	<- coef(fit0)["ln_lag_dol_purchases"]
kappa0	<- mean(mydata[!sel,"ln_dol_purchases"]) * rep(1, length(levels(mydata$income_group)))
names(kappa0) 	<- paste("factor(income_group)Qt",1:length(levels(mydata$income_group)),sep="")
sel1 	<- names(coef(fit0)[-1])
kappa0[sel1]		<- kappa0[sel1] + coef(fit0)[-1]		
sigma0	<- sum(fit0$residuals^2) / fit0$df.residual

fit1	<- plm(ln_dol_purchases ~ ln_lag_dol_purchases + factor(income_group) , data=mydata[sel,], 
				index = c("household_code","ymonth"), model = "within")
rho1	<- coef(fit1)["ln_lag_dol_purchases"]
kappa1	<- mean(mydata[sel,"ln_dol_purchases"]) * rep(1, length(levels(mydata$income_group)))
names(kappa1) 	<- paste("factor(income_group)Qt",1:length(levels(mydata$income_group)),sep="")
sel1 	<- names(coef(fit1)[-1])
kappa1[sel1]		<- kappa1[sel1] + coef(fit1)[-1]	
sigma1	<- sum(fit1$residuals^2) / fit1$df.residual

gh_num_nodes <- 5
y_weight	<- ghq(gh_num_nodes, modified = F)$weights
y_nodes 	<- ghq(gh_num_nodes, modified = F)$zeros
y_kernal	<- list(nodes = y_nodes, weight = y_weight, 
				kappa = cbind(kappa0, kappa1), rho = c(rho0, rho1), sigma = c(sigma0, sigma1))

# -------------------------------------------# 
# Model the income transition
tmp		<- data.table(mydata)
tmp		<- tmp[, list(income_group = unique(income_group)), by=list(household_code, recession, year)]
tmp		<- tmp[, income_lag := my_lag(income_group), by=list(household_code)]
tmp$income_lag <- factor(tmp$income_lag, levels=1:length(levels(mydata$income_group)))
sel 	<- tmp$recession == 1
tmp3	<- diag(length(levels(mydata$income_group)))
tmp1 	<- table(tmp[!sel,income_lag], tmp[!sel,income_group])	
sel 	<- which(apply(tmp1, 1, function(x) all(x==0)))
if(length(sel)>0){ tmp1[sel,] <- tmp3[sel,] }
tmp2 	<- table(tmp[sel,income_lag], tmp[sel,income_group])
sel 	<- which(apply(tmp2, 1, function(x) all(x==0)))
if(length(sel)>0){ tmp2[sel,] <- tmp3[sel,] }
tmp		<- rbind(tmp1/rowSums(tmp1), tmp2/rowSums(tmp2))
tmp[is.na(tmp)]	<- 0
Inc_weight <- tmp

# Average price
tmp		<- data.table(price_dat)
tmp		<- tmp[,list(price=median(price_paid_norm_index)), by=list(channel_type)]
price 	<- as.vector(tmp$price)

#######################
# Set up DP structure #
#######################
# ------------------------------ # 
# State variables 
summary(mydata$dol_purchases); quantile(mydata$dol_purchases, c(.01, .025, .975, .99))
summary(mydata$Q);  quantile(mydata$Q, c(.01, .025, .975, .99))

my_grid <- list(I = c(.1, .5, 1.5, 2, 3, 4, 5 ), 
				y = c(3, 4.8, 5.2, 5.5, 5.8, 6, 6.5),  
				Inc	= 1:length(levels(mydata$income_group)), 
				Rc = c(0, 1))
state	<- as.matrix(my_grid[[1]], ncol=1) 
for(i in 2:length(my_grid)){
	state <- matrix_expand(state, my_grid[[i]])
}
colnames(state) <- c("I", "y", "Inc","Rc")
summary(state)

# Rescale price and quantity	
my_scale
summary(mydata$Q); 
bskt_typeQ
price
price_rs	<- price*my_scale
bskt_typeQ_rs <- bskt_typeQ/my_scale

# Inclusive value: note to take log of y
tmp 		<- melt(omega_nodraw)
omega_list	<- lapply(1:K, function(i) 
			splinefun(log(tmp[tmp$bskt_type==i,"y"]), tmp[tmp$bskt_type==i,"value"], method="natural"))
omega_fn 	<- function(k, y){
	return(omega_list[[k]](y) )
}
tmp 		<- seq(0, 7, by=.1)
ggtmp		<- sapply(1:K, function(i) omega_fn(i,tmp))
dimnames(ggtmp) <- list(y=tmp, k=1:K)
ggtmp		<- melt(ggtmp)
if(make_plot){
	ggplot(ggtmp, aes(y, value)) + geom_point() + geom_line() + facet_wrap(~k)	
}

# Guess value functions 
param 		<- c(lambda = 3, tau_b = 0, tau1 = 1.5, tau2 = 3, tau3 = 3, tau4 = 70, tau5=4, tau6=4, tau7 = 10, tau8 =5, tau9 = 6, tau10= 15)
beta		<- .95
# DP structure 
DP_list	<- list(state 	= state, K = K, 
				# value_fn= (11 - 2*state[,4]) * log(state[,1] + 1), 
				value_fn= 0 , 
				param  	= param, 
				beta 	= beta,
				Q_grid	= bskt_typeQ_rs, 
				TR_grid = bskt_typeTR, 
				omega	= omega_fn, 
				price	= price_rs,
				y_kernal= y_kernal, 
				Inc_weight = Inc_weight,
				Iteration = 0, 
				status = 1
				)

# tmp			<- do.call(rbind, strsplit(names(bskt_type),"-"))
# tmp1 		<- cut(state[,1], quantile(my_grid$I, c(0:my.break)/my.break), include.lowest=T, labels=tmp[1:my.break,1])
# sel 		<- state[,4] == 0
# tmp2		<- c( as.character(cut(state[sel,2], quantile(my_grid$y, c(0:(my.break-1))/(my.break-1)), include.lowest=T, labels=tmp[seq(my.break+1, K-my.break+1, my.break),2])), 
# 				  as.character(cut(state[!sel,2], quantile(my_grid$y, c(0:(my.break-1))/(my.break-1)), include.lowest=T, labels=tmp[seq(1, K-my.break, my.break),2]))	)
# policy_k_guess <- bskt_type[paste(tmp1,tmp2,sep="-")]
policy_k_guess <- rep(1, nrow(state))
tmp			<- V_initC(param, DP_list, policy_k = policy_k_guess, policy_c = .5*(state[,1]+bskt_typeQ_rs[policy_k_guess]))
DP_list$value_fn <- tmp


pct	<- proc.time()
sol <- value_iteration_fnC(DP_list, Bellman_operatorC, print_level=1, control_list)
print(proc.time() - pct)
DP_list	<- sol

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
	ggtmp$Q			<- bskt_typeQ_rs[ggtmp$k]

	quartz()
	print(ggplot(ggtmp, aes(I, c, col=factor(Inc), linetype = factor(Rc), alpha=.3 )) + geom_point() + geom_line() + 
			facet_grid(k ~ y, labeller = "label_both") + 
			geom_abline(aes(intercept = Q, slope = 1), linetype=2) + 
			labs(title = paste("Policy function c at y grid",sep=""))
	)
}

##############
# Estimation #
##############
# Estimation function wrapper
DP_init	 <- DP_list
DAM_llfn <- function(param){
	return(ll_fnC(param, hh_index, TT_vec, init_state, DataState, choice_seq, DP_init, Bellman_operatorC, control_list, policy_k_guess))
}

theta_init <- list(DP_list$param, 
				c(lambda = 9, tau_b = 0, tau1 = 4, tau2 = 3, tau3 = 3, tau4 = 70, tau5=4, tau6=4, tau7 = 10, tau8 =5, tau9 = 6, tau10= 15), 
				c(lambda = 9, tau_b = 0, tau1 = 3, tau2 = 12, tau3 = 20, tau4 = 80, tau5=10, tau6=15, tau7 = 20, tau8 =10, tau9 = 15, tau10= 15)
				)
print(system.time(tmp <- DAM_llfn(theta_init[[1]])))

# Set up parallel computing 
mylog	<- paste("run_",run_id,"/maxLik_6_",model_name,"_seg",seg_id,"_init",1:length(theta_init),".txt",sep="")
for(i in 1:length(mylog)){ writeLines("", mylog[i]) }
mycore 	<- 3
cl		<- makeCluster(mycore, type = "FORK")
registerDoParallel(cl)

dp_sol_list	<- foreach(i = 1:length(theta_init)) %dopar% {
	sink(mylog[i], append=TRUE)
	cat(paste("Maximizing with initial value ",i," art at",sep=""))
	print(Sys.time())
	maxLik(DAM_llfn, start = theta_init[[i]], method = "BFGS", print = 10, fixed= 2)
	cat(paste("Finish maximizing with initial values ",i," at\n"),sep="")
	print(Sys.time())
	sink()
}
stopCluster(cl)
sapply(dp_sol_list, function(x) x$maximum)
sel		<- which.max(sapply(dp_sol_list, function(x) x$maximum))
dp_sol	<- dp_sol_list[[sel]]

save.image(paste("run_",run_id,"/result_6_",model_name,"_seg",seg_id,"_pst.rdata",sep=""))

cat("File has been saved.\n")



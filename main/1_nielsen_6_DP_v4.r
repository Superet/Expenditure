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
nchain		<- 3
# vid_save 	<- vid
# run_id_save <- run_id
vid_save	<- vid	<- 4
run_id_save	<- run_id <- 5
# arr_idx		<- as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))
seg_id		<- ceiling(arr_idx/nchain)
chain_id	<- arr_idx - (seg_id - 1) * nchain

# setwd("~/Documents/Research/Store switching/processed data/Estimation")
# sourceCpp("~/Documents/Research/Store switching/Exercise/main/1_CDAMT_functions.cpp")
setwd("/home/brgordon/ccv103/Exercise/run")
# setwd("/kellogg/users/marketing/2661703/Exercise/run")
# setwd("/sscc/home/c/ccv103/Exercise/run")

sourceCpp("1_CDAMT_functions.cpp")
model_name 	<- "MDCEV_a1b1"
make_plot	<- FALSE

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
param_mat <- read.csv(paste("run_",run_id,"/inits.csv",sep=""), header = T)

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
mytrim		<- .2
omega_arr	<- apply(omega_draw, c(2,3), mean, na.rm=T, trim = mytrim)
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

if(make_plot){
	tmp 		<- seq(10, 1000, by=10)
	ggtmp		<- sapply(1:D, function(i) omega_fn(tmp,i))
	dimnames(ggtmp) <- list(y=tmp, d=1:D)
	ggtmp		<- melt(ggtmp)
	ggplot(ggtmp, aes(y, value, col=factor(d))) + geom_point() + geom_line() 
}
cat("omega(100,1)/100 =", omega_fn(100,1)/100,"\n")

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
sel		<- param_mat$seg_id == seg_id & param_mat$chain == chain_id
param	<- as.vector(as.matrix(param_mat[sel,-(1:2)])); 
names(param) <- colnames(param_mat)[-c(1:2)]
cat("Inital parameters at seg =", seg_id, "are:\n"); print(param); cat("\n")

DP_list	<- list(state 	= state,
				value_fn= .1*Inc_nodes[state[,2]]/(1-beta) + 20* log(state[,1] + 1), 
				param  	= param, 
				beta 	= beta,
				# omega	= list(omega_y = as.numeric(colnames(omega_nodraw)), omega = omega_nodraw), 
				omega	= list(omega_y = unique(omega_draw_mean$y), omega = omega_arr), 
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
						NM_startsize 	= 1.5, 
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
	print(ggplot(ggtmp, aes(I, value, linetype = factor(Rc))) + geom_point() + geom_line() + 
			facet_grid(variable ~ Inc, scales="free_y") + 
			# geom_abline(aes(intercept = Q, slope = 1), linetype=2) + 
			labs(title = paste("Policy function",sep=""))
	)
	
	# Shape of objective function along c
	c_grid 		<- seq(.1, 4, .1)
	plot_state 	<- which(state[,1]<=2 & state[,2] ==3 & state[,3]==1) 
	state[plot_state,]
	tmp			<- plotvalueC(plot_state, c_grid, tmpsol, control_list)
	tmp[tmp==0] <- NA
	dimnames(tmp) <- list(plot_state, c(c_grid, NA))
	tmp_cstar	<- tmpsol$policy[plot_state, 3]
	ggtmp 		<- melt(tmp)
	names(ggtmp) <- c("s_index","c","value")
	sel 		<- is.na(ggtmp$c)
	ggtmp[sel,"c"] <- tmp_cstar
	ggtmp$max_pnt <- 0
	ggtmp[sel,"max_pnt"] <- 1
	ggtmp$I		<- state[ggtmp$s_index, 1]
	ggtmp$Inc	<- state[ggtmp$s_index, 2]
	ggtmp$Rc	<- state[ggtmp$s_index, 3]
	quartz()
	print(ggplot(ggtmp, aes(c, value, linetype = factor(Inc), shape = factor(Rc))) + geom_point(aes(col=factor(max_pnt))) + 
				geom_line() + 
				facet_grid(I~., labeller = "label_both") + 
				scale_color_manual(values=c("black","red")) + 
				guides(color=FALSE)
	)
}

##############
# Estimation #
##############
set.seed(666)
n_M 		<- 2 * (D+1)		# Number of moments; 
W			<- diag(n_M)/1e4
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

# ----------------------------------------------------------------------------------# 
# A wrapper function of GMM
par_neg_idx	<- 4
par_exp_idx <- c(1,2,4,5)
nest_pred	<- FALSE
GMM_wrapper <- function(param){
	GMM_objC(param, hh_index, TT_vec, init_state, DataState, choice_seq, W, DP_init, control_list, 
			logit_pred, Q_obs, omega_dif, phat, nest_pred = nest_pred)$obj
}

# Prior distribution 
lprior_fun 	<- function(param_conv){
	m		<- 0
	sigma	<- 10
	out		<- dnorm(param_conv, m, sigma, log = TRUE)
	return(sum(out))
}

lobj_fun 	<- function(param_conv){
	param	<- param_conv
	param[par_exp_idx] <- exp(param_conv[par_exp_idx])
	param[par_neg_idx] <- -param[par_neg_idx]
	L		<- - GMM_wrapper(param)
	if(is.na(L)){
		L 	<- -Inf
	}
	J		<- param_conv[par_exp_idx]
	J[which(par_exp_idx==par_neg_idx)] <- - J[which(par_exp_idx==par_neg_idx)]
	return(L + sum(J) + lprior_fun(param_conv))
}

my.mtretrop	<- function(lobj_fun, inits, scale=1, niter, nburn=0, nthin=1,  print_freq = 20){
	theta.old	<- inits
	p.old		<- lobj_fun(inits)
	np			<- length(inits)
	p.ser		<- rep(NA, niter)
	accept		<- rep(NA, niter)
	chain		<- matrix(NA, niter, np)
	cnt			<- 0
	prc			<- proc.time()
	for(i in 1:niter){
		theta.new	<- theta.old + rnorm(np, 0, scale)
		p.new		<- lobj_fun(theta.new)
		alpha		<- min(1, exp(p.new - p.old))
		if(runif(1) < alpha){
			accept[i]	<- 1
			chain[i,]	<- theta.new
			p.ser[i]	<- p.new
			p.old		<- p.new
			theta.old	<- theta.new
		}else{
			accept[i]	<- 0
			chain[i,]	<- theta.old
			p.ser[i]	<- p.old
		}
		if(i%%print_freq == 0){
			sel		<- cnt + 1:print_freq
			use.time<- proc.time() - prc
			cat("Iteration",cnt+1,"-", cnt+print_freq, ", acceptance =", round(mean(accept[sel]), 2), ", using", use.time[3]/60, "min.\n")
			cnt		<- cnt + print_freq
			prc		<- proc.time()
		}	
	}
	keep.idx 	<- setdiff(1:niter, 1:nburn)
	keep.idx	<- keep.idx[seq(1,length(keep.idx), nthin)]
	out.draw	<- chain[keep.idx,]
	out.dev		<- p.ser[keep.idx]
	
	# Compute DIC 
	Dbar 		<- mean(-2*out.dev)
	theta_bar	<- colMeans(out.draw)
	Dtilde		<- lobj_fun(theta_bar)
	pD			<- Dbar - Dtilde
	DIC 		<- pD + Dbar
	
	return(list(draw = out.draw, deviance = out.dev, accept = mean(accept), Dbar = Dbar, pD = pD, DIC = DIC, 
				last.draw = out.draw[nrow(out.draw),]))
}

# Set sampling paraemters 
set.seed(666)
param_init 	<- param
tmp			<- rep(1,length(par_exp_idx))
tmp[which(par_exp_idx==par_neg_idx)] <- -1
param_init[par_exp_idx]	<- log(param_init[par_exp_idx]*tmp)

# Search for scale parameter 
mygrid 		<- c(.02, seq(.05, .25, .05))
niter		<- 1000
scale_ls	<- vector("list", length(mygrid)); names(scale_ls)	<- mygrid; 
for(i in 1:length(mygrid)){
	pct			<- proc.time()
	scale_ls[[i]]<- my.mtretrop(lobj_fun, param_init, scale = mygrid[i], niter, nburn = 0, nthin = 1,  print_freq = 500)
	use.time	<- proc.time() - pct
}

cat("The accpetance rates along the grid are:\n"); sapply(scale_ls, function(x) x$accept); cat("\n")
sel			<- which.min(abs(sapply(scale_ls, function(x) x$accept) - .25 ))
cat("The scale paramter closest to .25 acceptance is", mygrid[sel],"\n")

mcmc.scale 	<- mygrid[sel]
niter		<- 12000
nburn		<- 2000
nthin		<- 40
print_freq	<- 500

pct			<- proc.time()
sol			<- my.mtretrop(lobj_fun, scale_ls[[sel]]$last.draw, scale = mcmc.scale, niter, nburn, nthin,  print_freq)
use.time	<- proc.time() - pct

# Check convergence
cat("Acceptance rate =", sol$accept, "\n")

pdf(paste("run_",run_id,"/graph_acf_seg", seg_id,"_chain", chain_id,".pdf", sep=""), width = 6, height = 8)
par(mfrow = c(3,2))
for(i in 1:length(param_init)){
	acf(sol$draw[,i], main = names(param_init)[i])
}
title(paste("ACF for chain ", chain_id, sep=""), outer = T)
dev.off()

if(chain_id==1){
	save.image(file = paste("run_",run_id,"/6_est_seg",seg_id,"_v", vid,"_chain", chain_id,".rdata",sep=""))
}else{
	save(sol, file = paste("run_",run_id,"/6_est_seg",seg_id,"_v", vid,"_chain", chain_id,".rdata",sep=""))
}

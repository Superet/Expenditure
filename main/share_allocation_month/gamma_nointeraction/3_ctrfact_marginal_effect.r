library(ggplot2)
library(reshape2)
library(maxLik)
library(evd)
library(data.table)
library(doParallel)
library(foreach)
# library(chebpol)
library(nloptr)
library(mgcv)
options(error = quote({dump.frames(to.file = TRUE)}))

# seg_id	<- as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))

args <- commandArgs(trailingOnly = TRUE)
print(args)
if(length(args)>0){
    for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
    }
}
seg_id <- 1
cat("seg_id =", seg_id, ".\n")

model_name 	<- "MDCEV_share"
# setwd("~/Documents/Research/Store switching/processed data/Estimation")
# plot.wd	<- '~/Desktop'
# source("../../Exercise/Multiple_discrete_continuous_model/0_Allocation_function.R")
# source("../../Exercise/main/share_allocation_month/ctrfact_sim_functions.r")

# setwd("/home/brgordon/ccv103/Exercise/run")
# setwd("/kellogg/users/marketing/2661703/Exercise/run")
# setwd("/sscc/home/c/ccv103/Exercise/run")
setwd("U:/Users/ccv103/Documents/Research/Store switching/run")

run_id		<- 12
plot.wd		<- getwd()
make_plot	<- TRUE
ww			<- 10
ar			<- .6

source("0_Allocation_function.R")
source("ctrfact_sim_functions.r")

# Load estimation data 
ver.date	<- "2016-07-19"
cpi.adj		<- TRUE

if(cpi.adj){
	loadf	<- paste("estrun_",run_id,"/MDCEV_cpi_est_seg",seg_id,"_", ver.date,".rdata",sep="")
}else{
	loadf	<- paste("estrun_",run_id,"/MDCEV_est_seg",seg_id,"_", ver.date,".rdata",sep="")
}
loadf
load(loadf)
rm(list = intersect(ls(), c("gamfit", "shr","model_name", "tmpdat")))

# Set simulation parameters
week.price		<- FALSE
interp.method	<- "spline"					# Spline interpolation or Chebyshev interpolation
exp.method		<- "Utility"				# Utility maximization or finding roots for first order condition
trim.alpha		<- 0
numsim			<- 500
draw.par		<- FALSE
sim.omega		<- FALSE
fname			<- paste("nostate_marginal_seg",seg_id,sep="")
if(draw.par){ fname <- paste(fname, "_pardraw", sep="") }
if(sim.omega){fname	<- paste(fname, "_simomega", sep="") }
fname			<- paste(fname, "_sim", numsim, "_", as.character(Sys.Date()), sep="")
cat("Output file name is", fname, ".\n")

#############
# Functions #
#############
incl_value_fn	<- function(param_est, nx, base, X1, X2, y, Q, price, R, Ra, qz_cons = 0, exp_outside = TRUE, quant_outside = TRUE, eps_draw = NULL, inits=NULL,...){
# param_est		...	parameter estimates
# X_array		...	Attributes array, dimension(R*Nobs*nx)
# y				... A Nobs-dimensional expenditure budget. 
# Q				... Q is the scalor of the quantity constraint.
	beta		<- param_est[1:nx]
	gamma0		<- param_est[(nx+1):(length(param_est)-1)]
	sigma		<- exp(param_est[length(param_est)])
	Nobs		<- length(y)
	gamma		<- sapply(X2, function(x) exp(x%*%gamma0))
	if(is.null(dim(gamma))){
		gamma	<- matrix(gamma, nrow = 1)
	}
	if(exp_outside){
		psi_R1 	<- param_out$psi_R1
		gamma	<- cbind(gamma, 1)
	}
	if(quant_outside){
		psi_R2	<- param_out$psi_R2
		gamma	<- cbind(gamma, 1)
	}
	
	# Compute the parameters for allocation function 
	if(is.null(eps_draw)){
		eps_draw 	<- rep(1, Nobs) %*% t(rgev(R)) * sigma				# matrix(rgev(Nobs*R), Nobs, R)
	}
	xbeta	<- sapply(1:R, function(i) X1[[i]] %*% beta)
	psi		<- exp(xbeta + eps_draw)
	if(exp_outside)		{psi <- cbind(psi, psi_R1) }
	if(quant_outside) 	{psi <- cbind(psi, psi_R2) }
	
	# Call allocation function 
	omega 	<- rep(NA, Nobs) 
	e		<- matrix(NA, Nobs, Ra)
	for(i in 1:Nobs){
		if(y[i] < 1e-5){
			omega[i]	<- 0 
			e[i,]		<- 0
		}else{
			sol	<- Allocation_fn(y = y[i], psi=psi[i,], gamma=gamma[i,], Q = Q, price=price[i,], 
					R=R, Ra=Ra, qz_cons = qz_cons, exp_outside = exp_outside, quant_outside = quant_outside, inits=inits, ...)
			omega[i]	<- sol$max
			e[i,]		<- sol$e
		}
	}
	
	# Return results 
	return(list(e = e, omega = omega))
}

SimWrapper_fn	<- function(omega_deriv, ln_inc, rec, lambda, param_est, base, X1, X2, price, eps_draw, 
							sim.y = NULL, method = "FOC", use.bound = FALSE, ret.sim = FALSE, par.draw = NULL, share.state = NULL){
# This function simulates expenditure and share, having known omega_fn_ls. 

# omega_deriv	...	A function of 2d spline 
# ln_inc		...	A vector of log(income) as simulation input data
# lambda		... A vector of utility parameters for purchase utility
# param_est		... A vector of parameters in the conditional allocation model
# base			... Integer index of which retailer has intercept of 0 in the conditional allocation model
# X_list		...	A list of retail attributes, each element is a matrix (length(ln_inc) \times nx)
# price			...	A matrix of price (length(ln_inc) \times R )
# esp_draw		... A matrix of random draw
# sim.y			... If given as a matrix (numsim \times length(ln_inc)), then expenditure is given and only expenditure share is simulated. 
# method		... The method of solving for optimal expenditure. Take value of "FOC" (first-order condition equation)
#					and "Utility" (directly maximizing utility function)
# use.bound		...	Logical variable indicating if use bounds in solving for optimal expenditure
# ret.sim		...	Logical variable indicating if return all the simulation values for each random raw, or 
#					return average. 
#=========================================================================================================#
	
	numsim	<- nrow(eps_draw)	
	R		<- length(X1)
	nx		<- ncol(as.matrix(X1[[1]]))
	if(is.null(sim.y)){
		y		<- matrix(NA, numsim, length(ln_inc))
	}else{
		y 		<- sim.y
	}
	if(is.null(par.draw)){
		shr.draw	<- matrix(0, numsim, length(param_est))
	}else{
		lambda.draw	<- par.draw[,c("lambda1", "lambda2")]
		shr.draw	<- par.draw[,setdiff(colnames(par.draw), c("lambda1", "lambda2"))]
	}
	e		<- array(NA, c(numsim, length(ln_inc), R), dimnames = list(iter = 1:numsim, lnInc = ln_inc, retailer = 1:R))
	omega	<- matrix(NA, numsim, length(ln_inc))
	omega1	<- matrix(NA, numsim, length(ln_inc))
	lvinc	<- sort(unique(ln_inc))
	
	# Set up model matrix
	for(i in 1:length(lvinc)){
		sel <- ln_inc == lvinc[i]
		lower	<- NULL		
		tmpX1	<- lapply(X1, function(x) x[sel,])
		tmpX2	<- lapply(X2, function(x) x[sel,])
		# Omega function conditional on income level 
		f 	<- function(y, ...){
			tmp0	<- data.frame(matrix(price[i,], nrow=1, dimnames = list(NULL, colnames(price))))
			newx	<- data.frame(lnInc = lvinc[i], y = y, tmp0)
			if(!is.null(share.state)){
				Xlag	<- matrix(share.state[i,], nrow = 1)
				newx$Xlag	<- Xlag
			}
			omega_deriv(newx, ...)
		}

		# Note that expenditure y is independent of eps_draw
		if(is.null(sim.y)){
			if(is.null(par.draw)){
				y[,sel]	<- solveExp_fn(lambda, f, lvinc[i], rec = rec[i], method = method, lower = lower)		# This is a scalor
			}else{
				y[,sel]	<- apply(lambda.draw, 1, function(x) solveExp_fn(lambda+x, f, lvinc[i], rec = rec[i], method = method, lower = lower) )
			}
		}

		for(j in 1:numsim){				
			tmpsol 		<- try(incl_value_fn(param_est+shr.draw[j,], nx, base, X1=tmpX1, X2=tmpX2, y=y[j,sel], Q=Inf, 
						price=matrix(price[sel,], nrow = sum(sel)), 
						R=R, Ra=R, qz_cons = 0, exp_outside = FALSE, quant_outside = FALSE, 
						eps_draw = rep(1,sum(sel)) %*% t(eps_draw[j,]) ), silent = TRUE)
			if(class(tmpsol) != "try-error"){
				omega[j,sel] 	<- tmpsol$omega	
				# tmp0			<- data.frame(matrix(price[i,], nrow=1, dimnames = list(NULL, colnames(price))))
				# newx 			<- data.frame(lnInc = lvinc[i], y = y[j,sel], tmp0)		
				# if(!is.null(share.state)){
				# 	newx$Xacs	<- matrix(share.state[sel,], nrow=1)
				# }
				# omega1[j,sel]	<- omega_deriv(newx)
				e[j,sel,]		<- tmpsol$e
			}
		}			
	}

	# out	<- cbind(colMeans(y, na.rm=T), colMeans(omega, na.rm = T), colMeans(omega1, na.rm = T), apply(e, c(2,3), mean, na.rm= T))	
	out	<- cbind(colMeans(y, na.rm=T), colMeans(omega, na.rm = T), apply(e, c(2,3), mean, na.rm= T))	
	
	if(length(fmt_name) == R){
		# colnames(out)	<- c("Exp", "omega", "spl.omega", fmt_name)
		colnames(out)	<- c("Exp", "omega", fmt_name)
	}else{
		# colnames(out)	<- c("Exp", "omega", "spl.omega", names(X1))
		colnames(out)	<- c("Exp", "omega", names(X1))
	}
	
	if(ret.sim){
		out	<- list(Average = out, y = y, Allocation = e)
	}
	return(out)
}

###############################
# Prepare simulation elements # 
###############################
# Data required: parameters, income level, price, retail attributes, and random draws
# For each simulation scenario, if income, price or retail attributes change, we need to re-simulate inclusive values.

# Set simulation parameters
lambda 	<- coef(lsol)	#coef(sol.top2)
shr.par	<- coef(sol)

#-----------------------#
# Construct income data # 
# Take the households' income in 2007 as basis 
selyr	<- 2007
tmp		<- data.table(subset(mydata, year %in% selyr))
tmp		<- tmp[,list(income = unique(income_midvalue)), by = list(household_code, year)]
sim.data<- data.frame(tmp)[,c("household_code","income")]
sim.unq		<- data.frame(income2007 = unique(sim.data[,-1]))

# Counterfactual scenario: income is lower by 10%. 
my.change		<- 0	#.1
# lnInc_08		<- lnInc + log(1 - my.change)
sim.unq$income2008	<- (1 - my.change) * sim.unq$income2007		
sim.unq$Inc07	<- log(sim.unq[,"income2007"])
sim.unq$Inc08	<- log(sim.unq[,"income2008"])
sim.unq			<- sim.unq[order(sim.unq$income2007),]
cat("dim(sim.unq) =", dim(sim.unq), "\n")
sim.unq			<- sim.unq[c(1,nrow(sim.unq)),]

#----------------------------#
# Average price in year 2007
if(week.price){
	tmp 		<- dcast(subset(price_dat, year %in% selyr),	 
						scantrack_market_descr + year + biweek ~ channel_type, value.var = "bsk_price_paid_2004") 
	price.07	<- setNames( colMeans(as.matrix(tmp[,4:(3+R)]), na.rm=T), fmt_name)
}else{
	tmp 		<- dcast(subset(price_dat, year %in% selyr),	 
						scantrack_market_descr + year ~ channel_type, value.var = "bsk_price_paid_2004") 
	price.07	<- setNames( colMeans(as.matrix(tmp[,3:(2+R)]), na.rm=T), fmt_name)
}
cat("The average price level in 2007:\n"); print(price.07);cat("\n")

# Average retail attributes in 2007
selcol	<- c("size_index", "ln_upc_per_mod", "ln_num_module","overall_prvt")
X_list0<- setNames( lapply(fmt_name, function(x) colMeans(as.matrix(subset(fmt_attr, channel_type == x & year%in% selyr)[,selcol]))), 
					 fmt_name)
cat("The average retail attributes in 2007:\n"); print(do.call(rbind, X_list0)); cat("\n")

# Expand X_list and price to match the nobs of income 
price.07	<- rep(1, nrow(sim.unq)) %*% matrix(price.07, nrow = 1)
colnames(price.07)	<- fmt_name
X1.07	<- X2.07	<- setNames(vector("list", R), fmt_name)
for(i in 1:R){
	tmp1	<- matrix(0, nrow(sim.unq), R-1)
	tmp2	<- matrix(0, nrow(sim.unq), R)
	if(i<beta0_base){
		tmp1[,i]	<- 0
	}else if(i>beta0_base){
		tmp1[,(i-1)]	<- 0
	}
	X1.07[[i]]	<- cbind(tmp1, rep(1, nrow(sim.unq)) %*% matrix(X_list0[[i]], nrow = 1))
	X2.07[[i]]	<- cbind(tmp2, rep(1, nrow(sim.unq)) %*% matrix(X_list0[[i]], nrow = 1))
	colnames(X1.07[[i]])	<- c(fmt_name[-beta0_base], selcol)
	colnames(X2.07[[i]])	<- c(fmt_name, selcol)
}

#-------------------#
# Take random draws #
set.seed(666)
eps_draw	<- matrix(rgev(numsim*R, scale = exp(shr.par["ln_sigma"])), numsim, R)
if(draw.par){
	par_se	<- c(sqrt(diag(vcov(sol.top2))), sqrt(diag(vcov(sol))) )
	par_se[is.na(par_se)]	<- 0
	par_draw <- sapply(par_se, function(x) rnorm(numsim, mean = 0, sd = x))
}else{
	par_draw <- NULL
}

##############
# Simulation #
##############
if(interp.method == "spline"){
	y.nodes		<- quantile(mydata$dol, c(0:50)/50)
	y.nodes		<- sort(unique(c(y.nodes , seq(600, 1000, 100)) ))
}else{
	# Set up Chebyshev interpolation
	GH_num_nodes<- 100
	y.interval	<- c(.1, 1000)
	y.nodes		<- chebknots(GH_num_nodes, interval = y.interval)[[1]]
}
numnodes<- length(y.nodes)

# Register parallel computing
mycore 	<- 3
# cl		<- makeCluster(mycore, type = "FORK")
cl		<- makeCluster(mycore, type = "PSOCK")			# Register cores in Windows
registerDoParallel(cl)
cat("Register", mycore, "core parallel computing. \n")

# Simulate expenditure and expenditure share. 
pct				<- proc.time()
sim.base07		<- SimWrapper_fn(omega_deriv, ln_inc = sim.unq$Inc07, rec = rep(0, nrow(sim.unq)), lambda = lambda, param_est = shr.par, base = beta0_base, 
					X1 = X1.07, X2 = X2.07, price = price.07, eps_draw = eps_draw, method = exp.method, ret.sim = TRUE, par.draw = par_draw,
					share.state = NULL)
use.time		<- proc.time() - pct						
cat("2007 Baseline simulation finishes with", use.time[3]/60, "min.\n")

#------------------------#
# Simulation of strategy #
# Assume only grocery stores change retail attributes. 
iv.change <- .1
(sel.retail	<- which(fmt_name == "Discount Store"))

#-------------------------------#
# Simulation of price reduction # 
pct			<- proc.time()
sim.ls07	<- foreach(v = selcol, .errorhandling = "remove", .packages=c("mgcv", "nloptr")) %dopar% {
	out1.avg	<- data.frame()
	out1.alc	<- array(NA, c(numsim, nrow(sim.unq), R), 
						dimnames = list(iter = 1:numsim, lnInc = sim.unq[,"Inc07"], retail = fmt_name))
	X1_new 		<- X1.07
	X2_new		<- X2.07
	if(v %in% c("ln_upc_per_mod", "ln_num_module")){
		X1_new[[sel.retail]][,v]	<- X1_new[[sel.retail]][,v] + log(1 + iv.change)
		X2_new[[sel.retail]][,v]	<- X2_new[[sel.retail]][,v] + log(1 + iv.change)
	}else{
		X1_new[[sel.retail]][,v]	<- X1_new[[sel.retail]][,v] * (1 + iv.change)
		X2_new[[sel.retail]][,v]	<- X2_new[[sel.retail]][,v] * (1 + iv.change)
	}
	
	out1	<- SimWrapper_fn(omega_deriv = omega_deriv, ln_inc = sim.unq[,"Inc07"], rec = rep(0, nrow(sim.unq)),
							lambda = lambda, param_est = shr.par, base = beta0_base, 
							X1 = X1_new, X2 = X2_new, price = price.07, sim.y = sim.base07$y, 
							eps_draw = eps_draw, method = exp.method, ret.sim = TRUE, par.draw = par_draw, share.state = shr.stat[,-beta0_base])
	out1.avg	<- data.frame(out1$Average, lnInc = sim.unq[,"Inc07"], retailer = fmt_name[sel.retail], 
							Var = v, change = iv.change )
	out1$Average<- out1.avg
	return(out1)
}
use.time	<- proc.time() - pct
cat("Parallel simulation of retail price finishes with", use.time[3]/60, "min.\n")

stopCluster(cl)
cat("Stop clustering. \n")

# Check the marginal effects. 
tmp0	<- colMeans(sim.base07$Allocation[,1,] >0 )
tmp1	<- t(sapply(sim.ls07, function(x) colMeans(x$Allocation[,1,]>0)))
tmp		<- tmp1 - rep(1, length(selcol)) %*% t(tmp0)
tmp		<- tmp[,sel.retail]
names(tmp)	<- selcol
cat("Impact on shopping incidence:\n"); print(tmp); cat("\n")

tmp0	<- sim.base07$Allocation[,1,]/sim.base07$y[,1]
tmp0	<- tmp0[,sel.retail]
tmp0	<- mean(tmp0[tmp0>0])
tmp1	<- lapply(sim.ls07, function(x) x$Allocation[,1,]/x$y[,1])
tmp1	<- sapply(tmp1, function(x) mean(x[,sel.retail][x[,sel.retail]>0]))
tmp		<- setNames(tmp1 - tmp0, selcol)
cat("Impact on conditional expenditure:\n"); print(tmp); cat("\n")

tmp0	<- sim.base07$Average[1,fmt_name]/sim.base07$Average[1,"Exp"]
tmp1	<- sapply(sim.ls07, function(x) as.matrix(x$Average[1,gsub(" ", ".", fmt_name)]/x$Average[1,"Exp"]))
tmp		<- setNames(tmp1[sel.retail,] - tmp0[sel.retail], selcol)
cat("Impact on exp. share:\n"); print(tmp); cat("\n")

################
# Save results #
################
rm(list  = intersect(ls(), c("ar", "args", "cl", "ggtmp", "ggtmp1", "ggtmp2", "i", "lastFuncGrad", "lastFuncParam", "make_plot", "mycore", 
			"myfix", "plot.wd", "s1_index", "sel", "selyr", "price", "sol", "sol.top", "sol.top2", "tmp", "tmp_coef", 
			"tmp_price", "tmp1", "tmp2", "use.time", "ver.date", "var1", "ww", "f", 
			"numnodes", "out", "out1", "pct", "tmpd1", "tmpd2", "tmpdat", "u", "W", "y", "y.nodes",
			"Allocation_constr_fn", "Allocation_fn", "Allocation_nlop_fn", "cheb.1d.basis", "cheb.basis", "chebfun", 
			"exp_fn", "expFOC_fn", "incl_value_fn", "mysplfun", "mytrimfun", "param_assignR", "simExp_fn", "SimOmega_fn", 
			"SimWrapper_fn", "solveExp_fn", "spl2dfun", "uP_fn", "uGrad_fn")))

save.image(paste("estrun_",run_id, "/", fname, ".rdata",sep=""))

cat("This program is done.\n")

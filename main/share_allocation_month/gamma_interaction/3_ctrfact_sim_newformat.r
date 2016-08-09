library(ggplot2)
library(reshape2)
library(maxLik)
library(evd)
library(data.table)
library(doParallel)
library(foreach)
library(nloptr)
library(mgcv)
options(error = quote({dump.frames(to.file = TRUE)}))

seg_id	<- as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))
args <- commandArgs(trailingOnly = TRUE)
print(args)
if(length(args)>0){
    for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
    }
}
seg_id 	<- 1
cat("seg_id =", seg_id, ".\n")

# setwd("~/Documents/Research/Store switching/processed data/Estimation")
# plot.wd	<- '~/Desktop'
# source("../../Exercise/Multiple_discrete_continuous_model/0_Allocation_function.R")
# source("../../Exercise/main/share_allocation_month/ctrfact_sim_functions.r")

# setwd("/home/brgordon/ccv103/Exercise/run")
# setwd("/kellogg/users/marketing/2661703/Exercise/run")
# setwd("/sscc/home/c/ccv103/Exercise/run")
# setwd("E:/Users/ccv103/Documents/Research/Store switching/run")
setwd("U:/Users/ccv103/Documents/Research/Store switching/run")

run_id		<- 11

plot.wd		<- getwd()
make_plot	<- TRUE
ww			<- 10
ar			<- .6

source("0_Allocation_function.R")
source("ctrfact_sim_functions.r")

# Load estimation data 
ver.date	<- "2016-07-22"
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
trim.alpha		<- 0.05
numsim			<- 500	#numsim1		#<- 1000
draw.par		<- FALSE
sim.omega		<- FALSE
fname			<- paste("ctrfact_newf_seg",seg_id,sep="")
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
# lambda 	<- coef(sol.top2)
lambda	<- coef(lsol)
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
my.change		<- 0
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
	price0	<- setNames( colMeans(as.matrix(tmp[,4:(3+R)]), na.rm=T), fmt_name)
}else{
	tmp 		<- dcast(subset(price_dat, year %in% selyr),	 
						scantrack_market_descr + year ~ channel_type, value.var = "bsk_price_paid_2004") 
	price0	<- setNames( colMeans(as.matrix(tmp[,3:(2+R)]), na.rm=T), fmt_name)
}
cat("The average price level in 2007:\n"); print(price0);cat("\n")

# Average retail attributes in 2007
selcol	<- c("size_index", "ln_upc_per_mod", "ln_num_module","overall_prvt")
X_list0	<- setNames( lapply(fmt_name, function(x) colMeans(as.matrix(subset(fmt_attr, channel_type == x & year%in% selyr)[,selcol]))), 
					 fmt_name)
cat("The average retail attributes in 2007:\n"); print(do.call(rbind, X_list0)); cat("\n")

# Expand X_list and price to match the nobs of income 
price.07	<- rep(1, nrow(sim.unq)) %*% matrix(price0, nrow = 1)
colnames(price.07)	<- fmt_name
X2.08	<- X2.07	<- X1.08	<- X1.07	<- setNames(vector("list", R), fmt_name)
for(i in 1:R){
	tmp2	<- tmp1	<- matrix(0, nrow(sim.unq), (R-1)*2)
	if(i<beta0_base){
		tmp1[,i]	<- 1
		tmp2[,i]	<- 1
		tmp2[,(R-1+i)]	<- 1
	}else if(i>beta0_base){
		tmp1[,(i-1)]	<- 1
		tmp2[,(i-1)]	<- 1
		tmp2[,(R-1+i-1)]	<- 1
	}
	X1.07[[i]]	<- cbind(tmp1, rep(1, nrow(sim.unq)) %*% matrix(X_list0[[i]], nrow = 1), matrix(0, nrow(sim.unq), length(selcol)))
	X1.08[[i]]	<- cbind(tmp2, rep(1, nrow(sim.unq)) %*% matrix(X_list0[[i]], nrow = 1), rep(1, nrow(sim.unq)) %*% matrix(X_list0[[i]], nrow = 1))
	colnames(X1.07[[i]])	<- colnames(X1.08[[i]])	<- c(fmt_name[-beta0_base], paste("rec*", fmt_name[-beta0_base],sep=""),selcol, paste("rec*", selcol,sep=""))
	
	tmp2	<- tmp1	<- matrix(0, nrow(sim.unq), R*2)
	tmp1[,i]<- 1
	tmp2[,i]<- 1
	tmp2[,(R+i)]	<- 1

	X2.07[[i]]	<- cbind(tmp1, rep(1, nrow(sim.unq)) %*% matrix(X_list0[[i]], nrow = 1), matrix(0, nrow(sim.unq), length(selcol)))
	X2.08[[i]]	<- cbind(tmp2, rep(1, nrow(sim.unq)) %*% matrix(X_list0[[i]], nrow = 1), rep(1, nrow(sim.unq)) %*% matrix(X_list0[[i]], nrow = 1))
	colnames(X2.07[[i]])	<- colnames(X2.08[[i]])	<- c(fmt_name, paste("rec*", fmt_name,sep=""),selcol, paste("rec*", selcol,sep=""))		
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

# Simulate expenditure and expenditure share. 
pct				<- proc.time()
sim.base07		<- SimWrapper_fn(omega_deriv, ln_inc = sim.unq$Inc07, rec = rep(0, nrow(sim.unq)), lambda = lambda, param_est = shr.par, base = beta0_base, 
					X1 = X1.07, X2 = X2.07, price = price.07, eps_draw = eps_draw, method = exp.method, ret.sim = TRUE, par.draw = par_draw,
					share.state = NULL)
use.time		<- proc.time() - pct						
cat("2007 Baseline simulation finishes with", use.time[3]/60, "min.\n")

pct				<- proc.time()
sim.base08		<- SimWrapper_fn(omega_deriv1, ln_inc = sim.unq$Inc08, rec = rep(0, nrow(sim.unq)), lambda = lambda, param_est = shr.par, base = beta0_base, 
					X1 = X1.08, X2 = X2.08, price = price.07, eps_draw = eps_draw, method = exp.method, ret.sim = TRUE, par.draw = par_draw,
					share.state = NULL)
use.time		<- proc.time() - pct						
cat("2008 Baseline simulation finishes with", use.time[3]/60, "min.\n")

# -------------- #
# Counterfactual #
# Introduct a new small format similar to convenience stores
# Scenarios:
# a(2): if intercept inherents the grocery stores, discount stores;
# b(2): if price inherents the intercept
 
ret.a	<- c("Discount Store")
sim.scn	<- c("Itself", "Small")
newf.sim07	<- newf.sim <- setNames(vector("list", length(ret.a)*length(sim.scn)), paste(ret.a, sim.scn, sep="_"))
eps_draw_new	<- cbind(eps_draw, rgev(numsim, scale = exp(shr.par["ln_sigma"])) )

cat(" ---------------------------------------------------------------------- \n")
pct		<- proc.time()
# Set new parameters
sel		<- which(fmt_name == ret.a)
shr.par.new	<- rep(NA, length(shr.par)+4)
beta1	<-  c(ifelse(sel==beta0_base, 0, ifelse(sel<beta0_base, shr.par[paste("beta_", sel, sep="")], shr.par[paste("beta_", sel-1, sep="")])), 
			  ifelse(sel==beta0_base, 0, ifelse(sel<beta0_base, shr.par[paste("beta_", sel+R-1, sep="")], shr.par[paste("beta_", sel-1+R-1, sep="")])) )
beta2	<-  c(shr.par[paste("gamma0_",c(sel, sel+R),sep="")] )
shr.par.new	<- c(shr.par[1:(R-1)], beta1[1], shr.par[(R):(2*(R-1))], beta1[2], shr.par[(2*(R-1)+1):nx], 
				 shr.par[(nx+1):(nx+R)], beta2[1], shr.par[(nx+R+1):(nx+2*R)], beta2[2], shr.par[(nx+2*R+1):length(shr.par)])
cat("The parameter vector for adding a new format from", ret.a, "is:\n"); print(shr.par.new); cat("\n")
	
for(j in 1:length(sim.scn)){
	# Set new attributes
	X2.08.new	<- X2.07.new	<- X1.08.new	<- X1.07.new	<- setNames(vector("list", R+1), c(fmt_name, "new"))

	for(k in 1:R){
		# Covariate matrix for baseline function
		tmp2	<- tmp1	<- matrix(0, nrow(sim.unq), 2*R)
		if(k<beta0_base){
			tmp1[,k]	<- 1
			tmp2[,k]	<- 1
			tmp2[,(R+k)]<- 1
		}else if(k>beta0_base){
			tmp1[,(k-1)]	<- 1
			tmp2[,(k-1)]	<- 1
			tmp2[,(k-1+R-1)]	<- 1
		}
		X1.07.new[[k]]	<- cbind(tmp1, rep(1, nrow(sim.unq)) %*% matrix(X_list0[[k]], nrow = 1), matrix(0, nrow(sim.unq), length(selcol)) )
		X1.08.new[[k]]	<- cbind(tmp2, rep(1, nrow(sim.unq)) %*% matrix(X_list0[[k]], nrow = 1), 
															rep(1, nrow(sim.unq)) %*% matrix(X_list0[[k]], nrow = 1) )

		# Covariate matrix for satiation													
		tmp1	<- tmp2	<- matrix(0, nrow(sim.unq), 2*(R+1))
		tmp1[,k]	<- 1
		tmp2[,c(k, k+R+1)]	<- 1		
		X2.07.new[[k]]	<- cbind(tmp1, rep(1, nrow(sim.unq)) %*% matrix(X_list0[[k]], nrow = 1), matrix(0, nrow(sim.unq), length(selcol)) )
		X2.08.new[[k]]	<- cbind(tmp2, rep(1, nrow(sim.unq)) %*% matrix(X_list0[[k]], nrow = 1), 
															rep(1, nrow(sim.unq)) %*% matrix(X_list0[[k]], nrow = 1) )
	}
	tmp2	<- tmp1	<- matrix(0, nrow(sim.unq), 2*R)
	if(sel!=beta0_base){
		tmp1[,R]	<- 1
		tmp2[,R]	<- 1
		tmp2[,2*R]	<- 1
	}
	if(sim.scn[j] == "Itself"){
		tmp		<- X_list0[[ret.a]]
	}else{
		tmp		<- X_list0[["Convenience Store"]]
		tmp["overall_prvt"]	<- X_list0[[ret.a]]["overall_prvt"]
	}
	X1.07.new[[(R+1)]]	<- cbind(tmp1, rep(1, nrow(sim.unq)) %*% matrix(tmp, nrow = 1), matrix(0, nrow(sim.unq), length(selcol)))
	X1.08.new[[(R+1)]]	<- cbind(tmp2, rep(1, nrow(sim.unq)) %*% matrix(tmp, nrow = 1),
	 											rep(1, nrow(sim.unq)) %*% matrix(tmp, nrow = 1))

	# Covariate matrix for satiation 
	tmp1	<- tmp2	<- matrix(0, nrow(sim.unq), 2*(R+1))
	tmp1[,(R+1)]	<- 1
	tmp2[,c(R+1, 2*(R+1))]	<- 1		
	X2.07.new[[(R+1)]]	<- cbind(tmp1, rep(1, nrow(sim.unq)) %*% matrix(tmp, nrow = 1), matrix(0, nrow(sim.unq), length(selcol)))
	X2.08.new[[(R+1)]]	<- cbind(tmp2, rep(1, nrow(sim.unq)) %*% matrix(tmp, nrow = 1),
	 											rep(1, nrow(sim.unq)) %*% matrix(tmp, nrow = 1))

	# Set prices 
	price.new		<- price.07
	price.new	<- cbind(price.new, new = price.07[,ret.a])
	
	sidx	<- j
	newf.sim07[[sidx]]	<- SimWrapper_fn(omega_deriv, ln_inc = sim.unq$Inc07, rec = rep(0, nrow(sim.unq)), lambda = lambda, 
						param_est = shr.par.new, base = beta0_base, sim.y = sim.base07$y, 
						X1 = X1.07.new, X2 = X2.07.new, price = price.new, eps_draw = eps_draw_new, method = exp.method, ret.sim = TRUE, par.draw = par_draw,
						share.state = NULL)
	newf.sim[[sidx]]	<- 	SimWrapper_fn(omega_deriv1, ln_inc = sim.unq$Inc08, rec = rep(0, nrow(sim.unq)), lambda = lambda, 
							param_est = shr.par.new, base = beta0_base, sim.y = sim.base08$y, 
							X1 = X1.08.new, X2 = X2.08.new, price = price.new, eps_draw = eps_draw_new, method = exp.method, ret.sim = TRUE, par.draw = par_draw,
							share.state = NULL)
}
use.time	<- proc.time() - pct
cat("Counterfactual with", ret.a, "finishes with", use.time[3]/60, "min.\n")

# Compare the difference 
fmt_name1	<- c(fmt_name, "new")
tmp1	<- c(sim.base08$Average[1,fmt_name],0)
tmp2	<- newf.sim[[1]]$Average[1,fmt_name1]
mkt.dif	<- setNames(tmp2/sum(tmp2)-tmp1/sum(tmp1), fmt_name1)
cat("Difference of market share:\n"); print(round(mkt.dif*100,2)); cat("\n")

################
# Save results #
################
rm(list  = intersect(ls(), c("ar", "args", "cl", "ggtmp", "ggtmp1", "ggtmp2", "i", "lastFuncGrad", "lastFuncParam", "make_plot", "mycore", 
			"myfix", "plot.wd", "s1_index", "sel", "selyr", "price", "sol", "sol.top", "sol.top2", "tmp", "tmp_coef", 
			"tmp_price", "tmp1", "tmp2", "use.time", "ver.date", "var1", "ww", "f", "dT","lag_nodes","sidx","price0", "X_ls0",
			"numnodes", "out", "out1", "pct", "tmpd1", "tmpd2", "tmpdat", "u", "W", "y", "y.nodes", "tmp.dif", "inc.tab", 
			"Allocation_constr_fn", "Allocation_fn", "Allocation_nlop_fn", "cheb.1d.basis", "cheb.basis", "chebfun", 
			"exp_fn", "expFOC_fn", "incl_value_fn", "mysplfun", "mytrimfun", "param_assignR", "simExp_fn", "SimOmega_fn", "tmpX",
			"SimWrapper_fn", "solveExp_fn", "spl2dfun", "uP_fn", "uPGrad_fn", "X_list", "beta0_7", "d.psi", "gamfit","j", "k")))

save.image(paste("estrun_",run_id, "/", fname, ".rdata",sep=""))

cat("This program is done.\n")
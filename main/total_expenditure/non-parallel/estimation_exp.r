cat("This program begins to run at", as.character(Sys.time()), ".\n")

library(ggplot2)
library(reshape2)
library(Rcpp)
library(RcppArmadillo)
library(maxLik)
library(evd)
library(data.table)
library(doParallel)
library(foreach)
library(nloptr)
library(mgcv)
options(error = quote({dump.frames(to.file = TRUE)}))

args <- commandArgs(trailingOnly = TRUE)
print(args)
if(length(args)>0){
    for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
    }
}

# seg_id	<- as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))

run_id			<- 1
make_plot		<- TRUE
interp.method	<- "spline"			# "cheb"
trim.alpha		<- 0 #0.05
cpi.adj			<- TRUE
week.price		<- FALSE

# setwd("~/Documents/Research/Store switching/Processed_data/processed data_20160809")
# plot.wd	<- '~/Desktop'
# sourceCpp(paste("../../Exercise/Multiple_discrete_continuous_model/", model_name, ".cpp", sep=""))
# source("../../Exercise/Multiple_discrete_continuous_model/0_Efficient_Allocation_function.R")
# source("../../Exercise/main/ctrfact_sim_functions.r")

# setwd("/home/brgordon/ccv103/Exercise/run")
# setwd("/kellogg/users/marketing/2661703/Exercise/run")
# setwd("/sscc/home/c/ccv103/Exercise/run")
setwd("U:/Users/ccv103/Documents/Research/Store switching/run")

source("0_Efficient_Allocation_function.R")
source("ctrfact_sim_functions.r")

ver.date<- ""
# (fname	<- paste("estrun_",run_id,"/MDCEV_est_", ver.date,".rdata",sep=""))
(fname	<- paste("estrun_",run_id,"/MDCEV_annual_est_2016-08-22.rdata",sep=""))
load(fname)

################################
# Simulate the inclusive value # 
################################
# Register parallel computing
mycore 	<- 3
# cl		<- makeCluster(mycore, type = "FORK")
cl		<- makeCluster(mycore, type = "PSOCK")			# Register cores in Windows
registerDoParallel(cl)
cat("Register", mycore, "core parallel computing. \n")

#-------------------------------------------- # 
# Compute the inclusive value with random draws
shr.par	<- coef(sol)
set.seed(687)
numsim 		<- 500
eps_draw	<- matrix(rgev(numsim*R, scale = exp(shr.par["ln_sigma"])), numsim, R)
Z1		<- as.matrix(cbind(1, Rec = mydata$rec, Z))

# Subset simulation data
sel			<- duplicated(mydata[,c("household_code", "year")])
sht.X1		<- lapply(X1, function(x) x[!sel,])
sht.X2		<- lapply(X2, function(x) x[!sel,])
sht.price	<- price[!sel,]
sht.y		<- y[!sel]
ret.idx		<- 26:30			# The column index in X1 to denote retail attributes

pct			<- proc.time()
omega_draw	<- foreach(i = 1:numsim, .packages= c("nloptr"), .combine = 'rbind', .errorhandling='remove') %dopar% {
	incl_value_fn(param_est = shr.par, nx, base=beta0_base, X1 = sht.X1, X2= sht.X2, y = sht.y, price = sht.price, R, 
					eps_draw = rep(1,length(sht.y)) %*% t(eps_draw[i,]), returnmax = TRUE)$omega
}
use.time	<- proc.time() - pct
cat("Simulation of omega with random draws finishes with", use.time[3]/60, "min.\n")
cat("--------------------------------------------------------\n")
stopCluster(cl)
cat("Stop clustering. \n")
summary(c(omega_draw))

# Fit spline funciton
# tmpdat	<- data.frame(omega = colMeans(omega_draw), y = sht.y, ln_inc = ln_inc[!sel])
tmpdat	<- data.frame(omega = apply(omega_draw/10000, 2, mean, na.rm=T, trim = 0.1), y = sht.y, ln_inc = ln_inc[!sel])
tmpdat$price	<- price[!sel,]
tmp		<- do.call(cbind, lapply(sht.X1, function(x) x[,ret.idx]))
tmpdat$X<- as.matrix(cbind(Z1[!sel,setdiff(colnames(Z1), c("1", "ln_inc"))], tmp))
omega_deriv	<- spl2dfun(tmpdat, fml = omega ~ te(y, ln_inc) + price + X)
rm(list = c("sht.X1", "sht.X2", "sht.price", "sht.y", "tmpdat"))

###################################
# Upper level expenditue decision # 
###################################
# Derivative of omega w.r.t. y
tmpdat	<- data.frame(y = y, ln_inc = ln_inc)
tmpdat$price	<- price
tmp		<- do.call(cbind, lapply(X1, function(x) x[,ret.idx]))
tmpdat$X<- as.matrix(cbind(Z1[,setdiff(colnames(Z1), c("1", "ln_inc"))], tmp))
o.deriv	<- omega_deriv(tmpdat, deriv = 1)
rm(list = "tmpdat")

M_fn	<- function(lambda){
# Moment function: lambda * omega'(y,Inc) - 1 = 0
	m1	<- (Z1 %*% lambda)* o.deriv - 1
	m	<- cbind(m1, m1*ln_inc)
	mbar<- colMeans(m)
	return(list(moment = mbar, omega.derive = o.deriv, m = m))
}

GMM_fn	<- function(lambda, W = NULL){
# We use unit diagnial matrix as the weighting matric in GMM
	mm 	<- M_fn(lambda)
	m	<- mm$moment
	if(is.null(W)){
		W	<- diag(length(m))
	}
	obj	<- - t(m) %*% W %*% m				# negative moment function for maxLik
	return(obj)
}

lsol	<- lm(1/o.deriv ~ Z1 - 1)
summary(lsol)

pct		<- proc.time()
init	<- matrix(c(coef(lsol), rep(0, ncol(Z1))), 2, ncol(Z1), byrow = T)
system.time(tmp <- GMM_fn(init[1,]))

tmp_sol	<- tmp_sol1 <- list(NULL)
for(i in 1:nrow(init)){
	tmp_sol[[i]]	<- maxLik(GMM_fn, start=init[i,], method="BFGS")
	tmp_sol1[[i]]	<- maxLik(GMM_fn, start=init[i,], method="NM")
}

# Choose the solution with max(-MM)
tmp		<- sapply(tmp_sol, function(x) ifelse(is.null(x$maximum), NA, x$maximum))
tmp1	<- sapply(tmp_sol1, function(x) ifelse(is.null(x$maximum), NA, x$maximum))
if(max(tmp1, na.rm =T) > max(tmp, na.rm=T)){
	tmp_sol	<- tmp_sol1
	cat("NM optimization has better results.\n")
}
sel 	<- which(abs(tmp - max(tmp, na.rm=T)) < 1e-4, arr.ind=T)
sel1	<- sapply(tmp_sol[sel], function(x) !any(is.na(summary(x)$estimate[,"Std. error"])) & !any(summary(x)$estimate[,"Std. error"]==Inf) )
sel		<- ifelse(sum(sel1)==1, sel[sel1], ifelse(sum(sel1)==0, sel[1], sel[sel1][1]))
sol.top	<- tmp_sol[[sel]]
if(!sol.top$code %in% c(0,1,2,8)){
	cat("Top level estimation does NOT have normal convergence.\n")
}
use.time <- proc.time() - pct
cat("--------------------------------------------------------\n")
summary(sol.top)
cat("Top level estimation finishes with", use.time[3]/60, "min.\n")

# 2nd step with updated W 
m			<- M_fn(lambda = coef(sol.top))
W			<- solve(var(m$m))
cat("Optimal weighting matrix is:\n"); print(W); cat("\n")
sol.top2	<- maxLik(GMM_fn, start=coef(sol.top), method="BFGS", W = W)
summary(sol.top2)
cat("Top level estimation finishes.\n")

####################
# Save the results #
####################
rm(list= intersect(ls(), c("Allocation_constr_fn","Allocation_fn","Allocation_nlop_fn","incl_value_fn","i","MDCEV_ll_fnC",
		  "MDCEV_LogLike_fnC","MDCEV_wrapper","tmp","tmp1","tmp2","tmp_sol","sel","sel1","sel2","param_assign","tmpX_list", 
		  "use.time", "pct", "uP_fn","uPGrad_fn", "theta_init", "make_plot", "ord","panelist","tmpidx","tmpn","cl", 
		  "tmp_sol1", "GMM_fn", "M_fn", "init", "m", "mycore", "param_assignR", "ggtmp", "ggtmp1", "tmpdat",
		  "interp.method", "trim.alpha", "mysplfun", "mytrimfun", "expFOC_fn", "exp_fn", "solveExp_fn","dT","i", "j", "k",  
		  "simExp_fn", "SimWrapper_fn", "SimOmega_fn", "cheb.1d.basis", "cheb.basis", "chebfun", "omega_parallel", 
		  "lastFuncGrad", "lastFuncParam", "args", "plot.wd", "tmpv1", "tmpv2", "beta.init")))

save.image(file = fname)

cat("This program is done. ")

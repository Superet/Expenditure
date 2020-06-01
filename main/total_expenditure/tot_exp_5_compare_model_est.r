# This script estimates model. The moduels follow as: 
# 1. Segment households based on initial income level and subset data; 
# 2. Prepare estimation data; 
# 3. Estimate the conditional allocation model at different initial values; 
# 4. Simulate inclusive values; 
# 5. Estimate upper level expenditure model. 

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

seg_id	<- as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))
run_id			<- 1
ver.date		<- "2016-09-11"
cat("seg_id =", seg_id,".\n")

# setwd("~/Documents/Research/Store switching/Processed_data")
# plot.wd	<- '~/Desktop'
# sourceCpp(paste("../Exercise/Multiple_discrete_continuous_model/", model_name, ".cpp", sep=""))
# source("../Exercise/Multiple_discrete_continuous_model/0_Efficient_Allocation_function.R")
# source("../Exercise/main/ctrfact_sim_functions.r")

# setwd("/home/brgordon/ccv103/Exercise/run")
# setwd("/kellogg/users/marketing/2661703/Exercise/run")
setwd("/sscc/home/c/ccv103/Exercise/run")
# setwd("U:/Users/ccv103/Documents/Research/Store switching/run")

(fname	<- paste("estrun_",run_id,"/MDCEV_est_seg", seg_id, "_", ver.date,".rdata",sep=""))
load(fname)
rm(list = intersect(ls(), c("sol.top", "sol.top2", "lsol","hh_dist", "hh_dpt", "shr", "omega_deriv", "omega_draw", "o.deriv")))

sourceCpp("MDCEV_exp_outside.cpp")
source("0_Efficient_Allocation_function.R")
source("ctrfact_sim_functions.r")

(fname	<- paste("estrun_",run_id,"/comp_MDCEV_est_seg", seg_id, "_", ver.date,".rdata",sep=""))

##########################
# Fit simultaneous model # 
##########################
sum(mydata$income_midvalue/12 - mydata$dol <= 0 )
sel 	<- mydata$income_mid/12 - mydata$dol > 0 #& (1:nrow(mydata) %in% sample(1:nrow(mydata), 5000))
e.n		<- as.matrix(mydata[sel,paste("DOL_", gsub("\\s", "_", fmt_name), sep="")] )
inc.n	<- as.vector(mydata[sel,"income_midvalue"]/12)
e1.n	<- mydata[sel,"income_mid"]/12 - rowSums(e.n)
pe_idx.n	<- 1 * (cbind(e.n, e1.n) > 0)
M.n		<- rowSums(pe_idx.n)
price.n	<- as.matrix(mydata[sel,paste("PRC_", gsub("\\s", "_", fmt_name), sep="")])
X1.n	<- lapply(X2, function(x) x[sel,])
X2.n	<- lapply(X2, function(x) x[sel,])
nx		<- ncol(X1.n[[1]])
d1.n	<- e1.n + 1

llfunc	<- function(theta){
	MDCEV_ll_fnC(param = theta, nx = nx, e = e.n, e1 = e1.n, p = price.n, X1 = X1.n, X2 = X2.n, 
				M = M.n, pe_idx = pe_idx.n, d1 = d1.n)$ll
}

# Normalize the coefficients for demographics relative to those with grocery 
beta.init	<- c(t(inits[c(1:2),]))
beta.init	<- c(beta.init, rowMeans(inits[-c(1:2),]))
gamma.init	<- c(t(inits2[c(1:2),]))
gamma.init	<- c(gamma.init, rowMeans(inits2[-c(1:2),]))
init		<- c(beta.init, gamma.init/100, -.5) 
names(init)	<- c(colnames(X1.n[[1]]), paste("gamma0_", colnames(X2[[1]]), sep=""), "ln_sigma")
system.time(tmp <- llfunc(init))

pct <- proc.time()
sol.n	<- maxLik(llfunc, start=init, method="BFGS")
use.time <- proc.time() - pct
cat("MDCEV estimation finishes with", use.time[3]/60, "min.\n")
print(summary(sol.n))
cat("--------------------------------------------------------\n")

########################
# Simulate expenditure # 
########################
# Simulate the expenditure # 
incl_value_fn <- function(param_est, nx, X1, X2, y, price, R, eps_draw = NULL, returnmax = TRUE, ...){
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
	Nobs		<- length(y)
	# Compute the parameters for allocation function 
	if(is.null(eps_draw)){
		eps_draw 	<- rep(1, Nobs) %*% t(rgev(R)) * sigma				# matrix(rgev(Nobs*R), Nobs, R)
	}
	xbeta	<- sapply(1:R, function(i) X1[[i]] %*% beta )
	psi		<- exp(xbeta + eps_draw)
	psi		<- cbind(psi, 1)
	gamma	<- cbind(gamma, 1)
	price	<- cbind(price, 1)
	
	# Call allocation function 
	out		<- Allocation_fn(y, psi, gamma,price, R + 1, returnmax)
	return(out)
}

par.n	<- coef(sol.n)
selyr	<- 2007
sel		<- mydata$year == selyr & !duplicated(mydata[,c("household_code", "year")]) 
mydata$ln_inc	<- log(mydata$income_mid)
sim.unq	<- mydata[sel,c("household_code", "ln_inc","rec")]
sim.unq$price	<- as.matrix(mydata[sel,paste("PRC_", gsub("\\s", "_", fmt_name), sep="")])
X1.07	<- lapply(X2, function(x) x[sel,])
X2.07	<- lapply(X2, function(x) x[sel,])

pct		<- proc.time()
ehat1	<- array(NA, c(numsim, nrow(sim.unq), R+1), dimnames = list(NULL, sim.unq$household_code, c(fmt_name, "outside")))
for(i in 1:numsim){
	ehat1[i,,]	<- incl_value_fn(par.n, nx, X1.07, X2 = X2.07, 
					y = exp(sim.unq$ln_inc)/12, price = sim.unq$price, R = R, 
					eps_draw = rep(1, nrow(sim.unq)) %*% t(eps_draw[i,]), returnmax = FALSE)
}
use.time	<- proc.time() - pct
cat("Simulation of simultaneous model finishes with", use.time[3]/60,"min.\n")

##########################
# Fit linear regressions # 
##########################
lm.fit	<- vector("list", length = R)
ehat.lm	<- matrix(NA, nrow(sim.unq), R, dimnames = list(sim.unq$household_code, fmt_name))
for(i in 1:R){
	tmpy		<- mydata[,paste("DOL_", gsub("\\s", "_", fmt_name[i]), sep="")]
	tmpx		<- cbind(rec = mydata[,"rec"], price = price[,i], X2[[i]][,-c(1:(R*2))])
	lm.fit[[i]]	<- lm(tmpy ~ tmpx)

	newx	<- cbind(1, Rec = 0, price = sim.unq$price[,i], X2.07[[i]][,-c(1:(R*2))])
	ehat.lm[,i]	<- newx %*% coef(lm.fit[[i]])
}

rm(list = intersect(ls(), c("Allocation_fn", "attr_mat", "beta.init", "cpi.ajd", "d1.n", "e1.n", "e.n", "exp_fn", 
			"expFOC_fn", "gamma.init", "i", "incl_value_fn", "inc.n", "init", "inits", "inits2", 
			"lastFuncGrad", "lastFuncParam", "llfunc", "max.method", "MDCEV_ll_fnC", "M.n", "model_name", "mysplfun", 
			"mytrimfun", "newx", "omega_draw", "omega.lm", "par.n", "pct", "use.time", "pe_idx.n", "price.n", 
			"sel", "selyr", "simExp_fn", "SimWrapper_fn", "theta.init", "tmp","tmp0","tmp3","tmpfit", "tmpsol","tmpx",
			"tmpy","u", "uP_fn", "uPGrad_fn", "use.time", "ver.date", "W","week.price", "X1.07", "X2.07", "X1.n", "X2.n",
			 "cpi.adj", "ln_inc", "max.dist", "numsim", "pan", "price", "X1", "X2", "rec","s1_index", "shr","solveExp_fn") ))

save.image(file = fname)

cat("This program is done. ")

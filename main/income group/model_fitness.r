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
library(Rcpp)
library(RcppArmadillo)
options(error = quote({dump.frames(to.file = TRUE)}))

seg_id <- 1
cat("seg_id =", seg_id, ".\n")

# setwd("~/Documents/Research/Store switching/Processed_data/processed data_20160809")
# plot.wd	<- '~/Desktop'
# sourceCpp("../../Exercise/Multiple_discrete_continuous_model/MDCEV_exp_outside.cpp")
# source("../../Exercise/Multiple_discrete_continuous_model/0_Efficient_Allocation_function.R")
# source("../../Exercise/main/ctrfact_sim_functions.r")

# setwd("/home/brgordon/ccv103/Exercise/run")
# setwd("/kellogg/users/marketing/2661703/Exercise/run")
# setwd("/sscc/home/c/ccv103/Exercise/run")
setwd("U:/Users/ccv103/Documents/Research/Store switching/run")

sourceCpp("MDCEV_exp_outside.cpp")
source("0_Efficient_Allocation_function.R")
source("ctrfact_sim_functions.r")

run_id 		<- 2
ver.date	<- "2016-08-24"
(loadf		<- paste("estrun_",run_id,"/marginal_seg",seg_id,"_sim500_", ver.date,".rdata",sep=""))
load(loadf)

ehat0 <- SimWrapper_fn(omega_deriv, lambda, shr.par, base = beta0_base, 
					X1 = X1.07, X2=X2.07, price = sim.unq$price, 
					demo_mat = as.matrix(cbind(Intercept = 1, rec = 0, ln_inc = sim.unq$ln_inc)), 
					eps_draw = matrix(0, 1, R), ret.idx = ret.idx, sim.y = sim.base07$y, ret.sim = FALSE)
								
##########################
# Fit linear regressions # 
##########################
lm.fit	<- vector("list", length = R)
selyr	<- 2007
# sel		<- mydata$year == selyr & !duplicated(mydata[,c("household_code", "year")]) # & mydata$household_code < 2250000
sel		<- mydata$year == selyr & !duplicated(mydata[,c("household_code", "year")]) # & mydata$household_code < 2250000

ehat.lm	<- matrix(NA, nrow(sim.unq), R)
for(i in 1:R){
	tmpx	<- as.matrix(cbind(Rec = mydata$rec, X1[[i]][,(ret.idx[1]:nx)]))
	tmpy	<- mydata[,paste("SHR_", gsub("\\s", "_", fmt_name), sep="")[i]]*mydata$dol
	lm.fit[[i]]	<- lm(tmpy ~ tmpx)
	
	sel1	<- attr_mat$channel_type == fmt_name[i]
	tmp		<- as.matrix(attr_mat[sel1,selcol1][sel,])
	newx	<- cbind(1, Rec = 0, tmp, tmp*0)
	ehat.lm[,i]	<- newx %*% coef(lm.fit[[i]])
}

##########################
# Fit simultaneous model # 
##########################
sum(mydata$income_midvalue/12 - mydata$dol <= 0 )
sel 	<- mydata$income_midvalue/12 - mydata$dol > 0
e.n		<- as.matrix(mydata[sel,paste("SHR_", gsub("\\s", "_", fmt_name), sep="")] * mydata[sel,"dol"])
inc.n	<- as.vector(mydata[sel,"income_midvalue"]/12)
e1.n	<- mydata[sel,"income_midvalue"]/12 - rowSums(e.n)
pe_idx.n	<- 1 * (cbind(e.n, e1.n) > 0)
M.n		<- rowSums(pe_idx.n)
price.n	<- as.matrix(mydata[sel,paste("PRC_", gsub("\\s", "_", fmt_name), sep="")])
X1.n	<- lapply(X1, function(x) x[sel,])
X2.n	<- vector("list", length = R)
for(i in 1:R){
	tmp	<- matrix(0, sum(sel), R)
	tmp[,i]	<- 1
	X2.n[[i]]	<- tmp
}
d1.n	<- e1.n + 1

llfunc	<- function(theta){
	MDCEV_ll_fnC(param = theta, nx = nx, e = e.n, e1 = e1.n, p = price.n, X1 = X1.n, X2 = X2.n, 
				M = M.n, pe_idx = pe_idx.n, d1 = d1.n)$ll
}
init	<- c(shr.par[1:nx], 5, 12, 6, 6, 20, 60, -.5)
system.time(tmp <- llfunc(init))

pct <- proc.time()
sol.n	<- maxLik(llfunc, start=init, method="BFGS")
use.time <- proc.time() - pct
cat("MDCEV estimation finishes with", use.time[3]/60, "min.\n")
print(summary(sol.n))
cat("--------------------------------------------------------\n")

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
tmpX1	<- tmpX2	<- vector("list", length = R)
for(i in 1:R){
	tmp	<- X1
}
ehat1	<- incl_value_fn(par.n, nx, X1.07, X2 = lapply(X2.07, function(x) x[,1:R]), 
				y = exp(sim.unq$ln_inc)/12, price = sim.unq$price, R = R, 
				eps_draw = matrix(0, nrow(sim.unq), R), returnmax = FALSE)

###########################
# Calculate model fitness # 
###########################
tmp	<- data.table(subset(mydata, year == selyr & household_code %in% sim.unq$household_code))
tmp	<- tmp[,list(DOL_Convenience_Store = mean(SHR_Convenience_Store*dol), DOL_Discount_Store = mean(SHR_Discount_Store*dol), 
			DOL_Dollar_Store = mean(SHR_Dollar_Store*dol), DOL_Drug_Store = mean(SHR_Drug_Store*dol), 
			DOL_Grocery	= mean(SHR_Grocery*dol), DOL_Warehouse_Club = mean(SHR_Warehouse_Club*dol) ), 
			by = list(household_code)]
identical(as.numeric(tmp$household_code), as.numeric(sim.unq$household_code))			
e.true	<- as.matrix(tmp)[,-1]
esign.true	<- 1*(e.true > 0)
# ehat.ls	<- list(Propose = sim.base07$Average[,fmt_name], Regression = ehat.lm, Simultaenous = ehat1[,1:R])
ehat.ls	<- list(Propose = ehat0[,fmt_name], Regression = ehat.lm, Simultaenous = ehat1[,1:R])
esignhat.ls	<- lapply(ehat.ls, function(x) 1*(x!=0))
sel		<- esign.true == 0
tmp.tab <- cbind(mse = sapply(ehat.ls, function(x) mean((x-e.true)^2)), 
				 accuracy = sapply(esignhat.ls, function(x) sum(x==esign.true)/length(e.true) ), 
				 hit.zero = sapply(esignhat.ls, function(x) sum((x==esign.true)[sel])/sum(sel) ) 	)
cat("Model fitness:\n"); print(tmp.tab); cat("\n")

# Expand the full data 
tmp		<- subset(mydata, year == selyr & household_code %in% sim.unq$household_code)
tmpn	<- setNames(1:nrow(sim.unq), sim.unq$household_code)
sel		<- tmpn[as.character(tmp$household_code)]
e.true	<- as.matrix(tmp[,paste("SHR_", gsub("\\s", "_", fmt_name), sep="")])*tmp$dol
esign.true	<- 1*(e.true > 0)
# ehat.ls	<- list(Propose = sim.base07$Average[,fmt_name], Regression = ehat.lm, Simultaenous = ehat1[,1:R])
ehat.ls	<- list(Propose = ehat0[sel,fmt_name], Regression = ehat.lm[sel,], Simultaenous = ehat1[sel,1:R])
esignhat.ls	<- lapply(ehat.ls, function(x) 1*(x!=0))
sel		<- esign.true == 0
tmp.tab <- cbind(mse = sapply(ehat.ls, function(x) mean((x-e.true)^2)), 
				 accuracy = sapply(esignhat.ls, function(x) sum(x==esign.true)/length(e.true) ), 
				 hit.zero = sapply(esignhat.ls, function(x) sum((x==esign.true)[sel])/sum(sel) ) 	)
cat("Model fitness:\n"); print(tmp.tab); cat("\n")

(fname	<- paste("comp_fit_seg", seg_id, "_", Sys.Date(), ".rdata", sep=""))

save.image(file = paste("estrun_", run_id, "/", fname, sep=""))

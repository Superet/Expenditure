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
library(plm)
# library(chebpol)
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
cat("seg_id =", seg.id, "\.\n")

# setwd("~/Documents/Research/Store switching/processed data")
# plot.wd	<- '~/Desktop'
# sourceCpp(paste("../Exercise/Multiple_discrete_continuous_model/", model_name, ".cpp", sep=""))
# source("../Exercise/Multiple_discrete_continuous_model/0_Allocation_function.R")

# setwd("/home/brgordon/ccv103/Exercise/run")
# setwd("/kellogg/users/marketing/2661703/Exercise/run")
setwd("/sscc/home/c/ccv103/Exercise/run")
model_name 	<- "MDCEV_share"
run_id		<- 4
# seg_id		<- 1
make_plot	<- TRUE
plot.wd		<- paste(getwd(), "/estrun_",run_id, sep="")
ww			<- 10
ar			<- .6
interp.method	<- "spline"			# "cheb"
trim.alpha		<- 0.05
cpi.adj			<- TRUE
nb				<- 100

ver.date	<- "2016-02-26"
if(cpi.adj){
	fname	<- paste("estrun_",run_id,"/MDCEV_cpi_est_seg",seg_id,"_", ver.date,".rdata",sep="")
}else{
	fname	<- paste("estrun_",run_id,"/MDCEV_est_seg",seg_id,"_", ver.date,".rdata",sep="")
}
fname
load(fname)
sourceCpp(paste(model_name, ".cpp", sep=""))
source("0_Allocation_function.R")
source("ctrfact_sim_functions_v2.r")

rm(list = c("gamfit","tmpdat", "shr", "y", "price", "X_list", "ln_inc", "s1_index"))

#############
# Functions #
#############
M_fn	<- function(lambda, omega_deriv, y, ln_inc, dT = NULL){
# Moment function: lambda * omega'(y,Inc) - 1 = 0
	o.deriv	<- rep(NA, length(y))
	if(interp.method == "cheb"){
		for(i in 1:length(lnInc)){
			sel				<- ln_inc == lnInc[i]
			o.deriv[sel]	<- omega_deriv[[i]](y[sel], deriv = 1, dT = dT[sel,])
		}
	}else{
		newx	<- data.frame(lnInc = ln_inc, y = y)
		o.deriv	<- omega_deriv(newx, deriv = 1)
	}
	m1	<- (lambda[1] + lambda[2] * ln_inc)* o.deriv - 1
	m	<- cbind(m1, m1*ln_inc)
	mbar<- colMeans(m)
	return(list(moment = mbar, omega.derive = o.deriv, m = m))
}

GMM_fn	<- function(lambda, omega_deriv, y, ln_inc, dT, W = NULL){
# We use unit diagnial matrix as the weighting matric in GMM
	mm 	<- M_fn(lambda, omega_deriv, y, ln_inc, dT = dT)
	m	<- mm$moment
	if(is.null(W)){
		W	<- diag(length(m))
	}
	obj	<- - t(m) %*% W %*% m				# negative moment function for maxLik
	return(obj)
}

my.bt1	<- function(idxb, coef1, coef2){	
# idxb	...	the index of households in the data
# coef1	...	parameter estimates from the bottom level 
# coef2	... parameter estimates from the top level
	
	# Outcome variables as matrix
	sel		<- paste("SHR_", gsub("\\s", "_", fmt_name), sep="")
	shr		<- as.matrix(mydata[idxb,sel])
	y		<- as.vector(mydata[idxb,"dol"])
	ln_inc	<- mydata[idxb,"ln_inc"]
	sel		<- paste("PRC_", gsub("\\s", "_", fmt_name), sep="")
	price	<- as.matrix(mydata[idxb,sel])

	# Select the index that has the positive expenditure
	tmp 	<- price * 1 *(shr> 0)
	s1_index<- apply(tmp, 1, which.max)

	# Match retailers' attributes
	tmp1	<- unique(fmt_attr$year)
	tmp2	<- unique(fmt_attr$scantrack_market_descr)
	tmp		<- paste(rep(tmp2, each=length(tmp1)), rep(tmp1, length(tmp2)), sep="-")
	tmpn	<- 1:length(tmp)		
	names(tmpn) <- tmp
	sel 	<- with(mydata[idxb,], paste(scantrack_market_descr,year, sep="-"))
	sel1	<- tmpn[sel]
	selcol	<- c("size_index", "ln_upc_per_mod", "ln_num_module","overall_prvt")
	nx 		<- length(selcol) * 2

	X_list 	<- vector("list", length=length(fmt_name))
	for(i in 1:length(fmt_name)){
		sel2		<- fmt_attr$channel_type == fmt_name[i]
		tmp			<- fmt_attr[sel2,selcol]
		tmp1 		<- as.matrix(tmp[sel1,])
		X_list[[i]]	<- cbind(tmp1, tmp1 * ln_inc)
	}
	
	MDCEV_wrapper <- function(param){
		MDCEV_ll_fnC(param, nx, shr, y, s1_index, price, X_list, beta0_base)$ll
	}
	sol.bot <- maxLik(MDCEV_wrapper, start=coef1, method="BFGS", fixed=myfix)
	out		<- coef(sol.bot)
	
	# Upper level 
	sol.top <- maxLik(GMM_fn, start=coef2, method="BFGS", omega_deriv = omega_deriv, y = y, ln_inc = ln_inc, W = W)
	out		<- c(out, coef(sol.top))
	
	return(out)
}

############################
# Bootstrap standard error #
############################
selcol

# Draw bootstrap sample index 
set.seed(666)
unq.hh	<- unique(mydata$household_code)
nh		<- length(unq.hh)
all.idx <- split(1:nrow(mydata), mydata$household_code)
idxb	<- lapply(1:nb, function(i){ unlist(all.idx[as.character(sample(unq.hh, nh, replace = T))]) })
coef1	<- coef(sol)
coef2	<- coef(sol.top2)
cat("The coefficient estimates are:\n"); print(c(coef1, coef2)); cat("\n")

pct		<- proc.time()
# beta.bt	<- sapply(1:nb, function(i) my.bt1(idxb[[i]], coef1, coef2))
beta.bt	<- matrix(NA, length(c(coef1, coef2)), nb)
for(i in 1:nb){
	beta.bt[,i]	<- my.bt1(idxb[[i]], coef1, coef2)
	print(i)
}
use.time	<- proc.time() - pct
cat("Boostrap standard procedure finishes using", use.time[3]/60, "min.\n")

# Compute se
se.bt		<- apply(beta.bt, 1, sd)
tmp.tab		<- cbind(c(coef1, coef2), se.bt)
cat("Estimates and bootstrap se:\n"); print(round(tmp.tab, 4)); cat("\n"); 

####################
# Save the results #
####################
ls()
rm(list=c("hh_exp", "Allocation_constr_fn","Allocation_fn","Allocation_nlop_fn","incl_value_fn","i","j","MDCEV_ll_fnC",
		  "MDCEV_LogLike_fnC","MDCEV_wrapper","tmp","tmp1","tmp2","tmp_sol","sel","sel1","sel2","selcol","param_assign", 
		  "use.time", "pct", "uP_fn","uPGrad_fn", "theta_init", "make_plot", "tmpsol","ord","panelist","tmpidx","tmpn"))

save.image(paste("estrun_",run_id,"/MDCEV_estbt_seg",seg_id, "_nb", nb, "_", Sys.Date(),".rdata",sep=""))

cat("This program is done. ")

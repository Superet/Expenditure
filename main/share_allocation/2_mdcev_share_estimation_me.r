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
library(lbfgs)
library(evd)
library(data.table)
library(doParallel)
library(foreach)
library(plm)
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
# cat("seg_id =", seg.id, "\.\n")

run_id			<- 5
seg_id			<- 1
make_plot		<- TRUE
interp.method	<- "spline"			# "cheb"
trim.alpha		<- 0.05
cpi.adj			<- TRUE

# setwd("~/Documents/Research/Store switching/processed data")
# plot.wd	<- '~/Desktop'
# sourceCpp("~/Documents/Research/Store switching/Exercise/Multiple_discrete_continuous_model/MDCEV_share_ME.cpp")
# source("../Exercise/Multiple_discrete_continuous_model/0_Allocation_function.R")

# setwd("/home/brgordon/ccv103/Exercise/run")
# setwd("/kellogg/users/marketing/2661703/Exercise/run")
setwd("/sscc/home/c/ccv103/Exercise/run")

sourceCpp("MDCEV_share_ME.cpp")
source("0_Allocation_function.R")

if(cpi.adj){
	fname	<- paste("estrun_",run_id,"/MDCEV_cpi_estme_seg",seg_id,"_", as.character(Sys.Date()), sep = "")
}else{
	fname	<- paste("estrun_",run_id,"/MDCEV_estme_seg",seg_id,"_", as.character(Sys.Date()), sep = "")
}
cat("The output file name is", fname, ".\n")

#################################
# Read data and subsetting data #
#################################
load("hh_biweek_exp.rdata")

# Extract 5% random sample
# length(unique(hh_exp$household_code))
# sel			<- sample(unique(hh_exp$household_code), .01*length(unique(hh_exp$household_code)) )
# hh_exp_save	<- hh_exp
# hh_exp		<- subset(hh_exp, household_code %in% sel)

# Retail formats 
fmt_name 	<- as.character(sort(unique(fmt_attr$channel_type)))
R			<- length(fmt_name)

# Convert some date and factor variables
hh_exp$month <- month(as.Date("2004-1-1", format="%Y-%m-%d") + 14*(hh_exp$biweek-1))
hh_exp$famsize		<- factor(hh_exp$famsize, levels = c("Single","Two", "Three+"))

# Compute expenditure share
sel			<- paste("DOL_", gsub("\\s", "_", fmt_name), sep="")
sel1		<- gsub("DOL", "SHR", sel)
for(i in 1:length(fmt_name)){
	hh_exp[,sel1[i]] <- hh_exp[,sel[i]]/hh_exp$dol
}
if(cpi.adj){
	hh_exp$income_midvalue	<- hh_exp$income_midvalue/hh_exp$cpi
}
hh_exp$ln_inc<- log(hh_exp$income_midvalue)

# Subset data 
selcol		<- c("household_code", "biweek", "dol", "year", "month", "income_midvalue", "ln_inc","first_incomeg", "scantrack_market_descr",paste("SHR_", gsub("\\s", "_", fmt_name), sep=""))
mydata		<- subset(hh_exp[,selcol], as.numeric(first_incomeg) == seg_id & dol > .1 )

################################
# Organize data for estimation #
################################
# The data required for estimation: 
# outcome variables: y (expenditure), shr (expenditure share);
# explanatory variables: price, X_list

# Format attributes
fmt_attr 		<- subset(fmt_attr, year > 2003)
fmt_attr$name 	<- paste(fmt_attr$scantrack_market_descr,fmt_attr$year, sep="-")
fmt_attr$ln_num_module 	<- log(fmt_attr$num_module)
fmt_attr$ln_upc_per_mod <- log(fmt_attr$avg_upc_per_mod)

# Match biweekly prices
price_dat 		<- subset(price_dat, year > 2003)
tmp		<- price_dat
tmp$name<- gsub("\\s", "_", tmp$channel_type)
tmp$name<- paste("PRC_", tmp$name, sep="")
price 	<- dcast(tmp, scantrack_market_descr + year ~ name, value.var = "bsk_price_paid_2004")

# Merge price data to the trip data
mydata <- merge(mydata, price, by=c("scantrack_market_descr", "year"), all.x=T)
ord		<- order(mydata$household_code, mydata$biweek)
mydata	<- mydata[ord,]
nh			<- length(unique(mydata[,"household_code"]))
idx			<- setNames(1:nh, unique(mydata[,"household_code"]))
mydata$id	<- idx[as.character(mydata$household_code)]

# Outcome variables as matrix
sel		<- paste("SHR_", gsub("\\s", "_", fmt_name), sep="")
shr		<- as.matrix(mydata[,sel])
y		<- as.vector(mydata$dol)
ln_inc	<- mydata$ln_inc
sel		<- paste("PRC_", gsub("\\s", "_", fmt_name), sep="")
price	<- as.matrix(mydata[,sel])
idx		<- mydata$id
M		<- rowSums(shr > 0)

# Select the index that has the positive expenditure
tmp 	<- price * 1 *(shr> 0)
s1_index<- apply(tmp, 1, which.max)

# Match retailers' attributes
tmp1	<- unique(fmt_attr$year)
tmp2	<- unique(fmt_attr$scantrack_market_descr)
tmp		<- paste(rep(tmp2, each=length(tmp1)), rep(tmp1, length(tmp2)), sep="-")
tmpn	<- setNames(1:length(tmp), tmp)
sel 	<- with(mydata, paste(scantrack_market_descr,year, sep="-"))
sel1	<- tmpn[sel]						# The index that matches the estimation data
selcol	<- c("size_index", "ln_upc_per_mod", "ln_num_module","overall_prvt")		# Covariates 

beta0_base 	<- which(fmt_name == "Grocery")						# The retailer for which the intercept is set 0 
X_list 	<- vector("list", length=length(fmt_name))
for(i in 1:length(fmt_name)){
	sel2		<- fmt_attr$channel_type == fmt_name[i]
	tmp			<- fmt_attr[sel2,selcol]
	tmp1 		<- as.matrix(tmp[sel1,])
	tmp2		<- matrix(0, nrow(shr), R-1)
	if(i < beta0_base){
		tmp2[,i]	<- ln_inc
	}else if(i > beta0_base){
		tmp2[,(i-1)] <- ln_inc
	}
	X_list[[i]]	<- cbind(tmp2, tmp1, tmp1 * ln_inc)
}
nx 		<- length(selcol) * 2 + R-1			# Number of covariates + interaction * ln_income + ln_income

##############
# Estimation #
##############
ll_L1	<- function(param){
	MDCEV_ll_fnC(param, nx = nx, shr = shr, y = y, p = price, X_list = X_list, idx = idx, nh = nh, M = M, base = beta0_base)
}

grad_L1	<- function(param){
	MDCEV_grad_fnC(param, nx = nx, shr = shr, y = y, p = price, X_list = X_list, idx = idx, nh = nh, M = M, base = beta0_base)
}

# Set initial values
theta.init	<- list(c(rep(-.01, R-1), -1, -.1, .5, -1.5, 	.1, .1, -.1, .1, 	-1, -1, -1, -.5, -1,	1.8, 3.3, 2, 1.8, 2.5, 4, -.3, rep(0.1, nh*(R-1))),
					c(rep(-.01, R-1), -.5,-.1, -.3, -.4, 	.05,.1, .03, .05, 	-1, -.8, -.8,-.3,-1,  	1.7, 3.3, 2, 1.8, 2.7, 4.1, -.4, rep(0.1, nh*(R-1))), 
					c(rep(-.01, R-1), -.8,-.4, 1, -6, 		.1, .1, -.1, .6, 	-1.5,-.8,-1.4,-.3,-.3,	1.5, 3.4, 2, 2, 2, 2.9, 4.3,-.4, rep(0.1, nh*(R-1)))
					)
system.time(print(ll_L1(theta.init[[1]])))
system.time(tmp <- grad_L1(theta.init[[1]]))

L1C			<- 1000			# Average number of observations per household
pct			<- proc.time()
est.l1 <- lbfgs::lbfgs(ll_L1, grad_L1, vars = theta.init[[seg_id]], max_iterations = 1000, 
						orthantwise_c=L1C, orthantwise_start = nx+2*R+1, orthantwise_end = length(theta.init[[1]]))
use.time 	<- proc.time() - pct
cat("Estimation with penalty finishes using", use.time[3]/60, "min.\n")

# Extract parameter estimates
shr.par		<- setNames(est.l1$par[1:(nx+2*R)], c(paste("beta_",1:nx, sep=""), paste("beta0_",setdiff(1:R, beta0_base),sep=""), 
							 paste("gamma_",1:R,sep=""), "ln_sigma"))
shr.par[grep("gamma", names(shr.par))]	<- exp(shr.par[grep("gamma", names(shr.par))])
cat("Estimation of fixed effects are\n"); print(shr.par); cat("\n")
est.hes		<- optimHess(est.l1$par, fn = ll_L1, gr = grad_L1)
est.hes.fe	<- est.hes[1:(nx+2*R),1:(nx+2*R)]
par.v		<- chol2inv(chol(est.hes.fe))
shr.se 		<- setNames(sqrt(diag(par.v1)), names(shr.par) )
cat("s.e. of parameters are:\n"); print(shr.se); cat("\n")

shr.re		<- setNames(est.l1$par[(nx+2*R+1):length(theta.init[[1]])], 
					paste("re", rep(1:nh, each = R-1), "_", rep(setdiff(1:R, beta0_base), nh), sep="") )
mean(shr.re == 0 )

####################
# Save the results #
####################
rm(list=intersect(ls(), c("hh_exp", "Allocation_constr_fn","Allocation_fn","Allocation_nlop_fn","incl_value_fn","i","MDCEV_ll_fnC",
		  "MDCEV_LogLike_fnC","tmp","tmp1","tmp2","tmp_sol","sel","sel1","sel2","param_assign","tmpX_list", 
		  "use.time", "pct", "uP_fn","uPGrad_fn", "theta_init", "make_plot", "ord","panelist","tmpidx","tmpn","cl", 
		  "tmp_sol1", "GMM_fn", "M_fn", 
		  "interp.method", "trim.alpha", "mysplfun", "mytrimfun", "expFOC_fn", "exp_fn", "solveExp_fn", 
		  "simExp_fn", "SimWrapper_fn", "SimOmega_fn", "cheb.1d.basis", "cheb.basis", "chebfun", "omega_parallel", 
		  "lastFuncGrad", "lastFuncParam", "args")) )

save.image(file = paste(fname, ".rdata", sep="") )

cat("This program is done. ")
						
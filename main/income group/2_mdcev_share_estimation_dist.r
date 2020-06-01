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
setwd("U:/Users/ccv103/Documents/Research/Store switching/run")
run_id			<- 2
seg_id			<- 1
load(paste("estrun_",run_id,"/MDCEV_est_seg", seg_id, "_2016-08-22.rdata",sep=""))
theta_init		<- coef(sol)
theta_init[c(11,16)]	<- theta_init[c(11,16)] * 100
rm(list = setdiff(ls(), "theta_init"))

model_name 		<- "MDCEV_share_gamma"
make_plot		<- TRUE
interp.method	<- "spline"			# "cheb"
trim.alpha		<- 0 #0.05
cpi.adj			<- TRUE
week.price		<- FALSE
run_id			<- 2
seg_id			<- 1
cat("seg_id =", seg_id, "\n")

# setwd("~/Documents/Research/Store switching/Processed_data/processed data_20160809")
# plot.wd	<- '~/Desktop'
# sourceCpp(paste("../../Exercise/Multiple_discrete_continuous_model/", model_name, ".cpp", sep=""))
# source("../../Exercise/Multiple_discrete_continuous_model/0_Efficient_Allocation_function.R")
# source("../../Exercise/main/ctrfact_sim_functions.r")

# setwd("/home/brgordon/ccv103/Exercise/run")
# setwd("/kellogg/users/marketing/2661703/Exercise/run")
# setwd("/sscc/home/c/ccv103/Exercise/run")

sourceCpp(paste(model_name, ".cpp", sep=""))
source("0_Efficient_Allocation_function.R")
source("ctrfact_sim_functions.r")

(fname	<- paste("estrun_",run_id,"/MDCEV_est_seg", seg_id, "_", Sys.Date(),".rdata",sep=""))

#################################
# Read data and subsetting data #
#################################
load("hh_month_exp.rdata")
load("hh_fmt.rdata")

# Extract 5% random sample
# length(unique(hh_exp$household_code))
# sel			<- sample(unique(hh_exp$household_code), .01*length(unique(hh_exp$household_code)) )
# hh_exp_save	<- hh_exp
# hh_exp		<- subset(hh_exp, household_code %in% sel)

# Retail formats 
fmt_name 	<- as.character(sort(unique(fmt_attr$channel_type)))
R			<- length(fmt_name)

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
hh_exp$rec		<- 1*(hh_exp$year >= 2008)

cat("Table of segments in the expenditure data:\n"); print(table(hh_exp$first_incomeg)); cat("\n")

# Fill in the missing distance for the channels that are not available for households; 
tmp			<- data.table(hh_dist)
tmp			<- tmp[, list(nozip = all(is.na(panelist_zip_code))), by = list(household_code, year)]
tmp			<- tmp[nozip==TRUE,]
unique(tmp$household) %in% hh_exp$household_code

max.dist	<- 100
sum(hh_dist$distance_haver >max.dist|is.na(hh_dist$distance_haver), na.rm=T)/nrow(hh_dist)
sel			<- hh_dist$distance_haver >max.dist|is.na(hh_dist$distance_haver)
hh_dist[sel,"distance_haver"]	<- max.dist

# Subset data 
tmp		<- data.table(hh_dist)
tmp		<- tmp[,list(nchannel = length(channel_type)), by = list(household_code, year)]

selcol		<- c("household_code", "dol", "year", "month", "rec", "income_midvalue", "ln_inc","first_incomeg", "scantrack_market_descr",paste("SHR_", gsub("\\s", "_", fmt_name), sep=""))
mydata		<- subset(hh_exp[,selcol], dol > .1 & as.numeric(first_incomeg) == seg_id)
sel			<- tmp$nchannel == R
sel			<- paste(mydata$household_code, mydata$year, sep="*") %in% paste(tmp[sel,household_code], tmp[sel,year], sep="*")
mydata		<- mydata[sel,]

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
if(week.price){
	tmp		<- price_dat
	tmp$name<- gsub("\\s", "_", tmp$channel_type)
	tmp$name<- paste("PRC_", tmp$name, sep="")
	price 	<- dcast(tmp, scantrack_market_descr + biweek ~ name, value.var = "bsk_price_paid_2004")
	tmpidx	<- 2

	# Impute the missing price for some regioins;
	sum(is.na(price))
	tmp			<- data.table(price_dat)
	tmp			<- tmp[,list(price = mean(bsk_price_paid_2004, na.rm=T)), by=list(biweek, channel_type)]
	tmp			<- data.frame(tmp)
	sel			<- which(is.na(as.matrix(price[,(tmpidx+1):ncol(price)])), arr.ind=T)
	dim(sel)
	fmt_name[unique(sel[,2])]
	for(i in 1:length(unique(sel[,2]))){
		sel1		<- unique(sel[,2])[i]
		sel2		<- sel[sel[,2]==sel1,1]
		tmp1		<- merge(subset(tmp, channel_type==fmt_name[sel1]), price[sel2,1:3], by=c("biweek"))
		tmp2		<- as.vector(tmp1$price)
		names(tmp2)	<- with(tmp1, paste(scantrack_market_descr, biweek,sep="-"))
		price[sel2, (sel1+tmpidx)] <- tmp2[with(price[sel2,], paste(scantrack_market_descr, biweek,sep="-"))]
	}
	
	# Merge price data to the trip data
	mydata <- merge(mydata, price, by=c("scantrack_market_descr", "biweek"), all.x=T)
}else{
	# Merge price data to the trip data
	tmp		<- price_dat
	tmp$name<- gsub("\\s", "_", tmp$channel_type)
	tmp$name<- paste("PRC_", tmp$name, sep="")
	price 	<- dcast(tmp, scantrack_market_descr + year ~ name, value.var = "bsk_price_paid_2004")
	mydata <- merge(mydata, price, by=c("scantrack_market_descr", "year"), all.x=T)
}
ord		<- order(mydata$household_code, mydata$year, mydata$month)
mydata	<- mydata[ord,]

# Outcome variables as matrix
sel		<- paste("SHR_", gsub("\\s", "_", fmt_name), sep="")
shr		<- as.matrix(mydata[,sel])
y		<- as.vector(mydata$dol)
ln_inc	<- mydata$ln_inc
sel		<- paste("PRC_", gsub("\\s", "_", fmt_name), sep="")
price	<- as.matrix(mydata[,sel])

# Select the index that has the positive expenditure
tmp 	<- price * 1 *(shr> 0)
s1_index<- apply(tmp, 1, which.max)

# Match retailers' attributes
tmp1	<- unique(fmt_attr$year)
tmp2	<- unique(fmt_attr$scantrack_market_descr)
tmp		<- paste(rep(tmp2, each=length(tmp1)), rep(tmp1, length(tmp2)), sep="-")
tmpn	<- 1:length(tmp)		
names(tmpn) <- tmp
sel 	<- with(mydata, paste(scantrack_market_descr,year, sep="-"))
sel1	<- tmpn[sel]
selcol	<- c("size_index", "ln_upc_per_mod", "ln_num_module","overall_prvt")
nx 		<- length(selcol)*2 + (R-1)*2 
beta0_base 	<- which(fmt_name == "Grocery")

# Normalize distance by zipcode
tmp		<- data.table(hh_dist[,c("household_code", "panelist_zip_code","year", "channel_type", "distance_haver")])
tmp		<- tmp[, maxd:=max(distance_haver, na.rm=T), by = list(panelist_zip_code)]
tmp		<- tmp[, distance_haver := distance_haver/maxd]
tmp		<- data.frame(tmp)

attr_mat	<- merge(mydata[,c("household_code", "year", "month")], tmp[,c("household_code", "year", "channel_type", "distance_haver")], 
				by = c("household_code", "year"), all.x=T)
attr_mat	<- attr_mat[order(attr_mat$household_code, attr_mat$year, attr_mat$month),]	
cat("dim(mydata) =", dim(mydata), ".\n")
cat("dim(attr_mat) =", dim(attr_mat), ", dim(mydata)*R =",nrow(mydata)*R, ".\n")

X2	<- X1 	<- vector("list", length=length(fmt_name))
for(i in 1:length(fmt_name)){
	sel2		<- fmt_attr$channel_type == fmt_name[i]
	tmp			<- fmt_attr[sel2,selcol]
	tmp1 		<- as.matrix(cbind(distance = attr_mat[attr_mat$channel_type==fmt_name[i],"distance_haver"], tmp[sel1,]))
	tmp2		<- matrix(0, nrow(shr), (R-1)*2)
	if(i < beta0_base){
		tmp2[,i]	<- 1
		tmp2[,(R-1)+i]	<- mydata$rec
	}else if(i > beta0_base){
		tmp2[,(i-1)] <- 1
		tmp2[,((R-1) + (i-1))] <- mydata$rec
	}
	tmp3	<- matrix(0, nrow(shr), R*2)
	tmp3[,i]<- 1
	tmp3[,(R+i)]	<- mydata$rec
	X1[[i]]	<- cbind(tmp2, tmp1, tmp1*mydata$rec)
	X2[[i]]	<- cbind(tmp3, tmp1, tmp1*mydata$rec)
}
nx 		<- ncol(X1[[1]])
cat("nx =", nx,"\n")

##############
# Estimation #
##############
# Set estimation control: fixing parameters
fixgamma 	<- FALSE
fixsigma 	<- FALSE
if(fixgamma){
	myfix 	<- (nx+R):(nx+2*R-1)
}else{
	myfix 	<- NULL
}

if(fixsigma){
	myfix	<- c(myfix, nx+2*R)
}

MDCEV_wrapper <- function(param){
	MDCEV_ll_fnC(param, nx, shr, y, s1_index, price, X1, X2, beta0_base)$ll
}

# # Estimation with multiple initial values. 
# theta_init	<- list(c( -1, -.8, -1, -.4,-1.2,	.2, .03, -.2, -.06, .1, 		.01, .2, .5, .2, 2.3,	0, -.1, .1, -.2,-.6, 
# 					  1.5, 3.5, 2.5, 2.5, 3.5, 4.5,		-.2, -.5, -.5, .2, 0, 1,	.01, 0,.1,-.1,-2.5, 		0,-.1,.1,.1,1.5,	 	-.45), 
# 					c(-1, -.7, -1.2, -.3,-.8,	.2, .05,.4,.1,.5,				.01, .1, .6, .1,2, 		0, -.05, .1, -.15, -.5,
# 					  1.8, 3.5, 2.1, 2.1,3,4, 	-.1, -.7, .1, .1, 0, 1.5,	.01, 0,.1,0,-2.7, 		0, 0, .2,-.3,1.5,	-.5),
# 					c(.7, -.5, 1.4, -1.2,-.3,	-.12,-.38,-.15,.1,	1.8, 	.01, .1, .4, .37, -.3, 	rep(-.1, length(selcol)+1),
# 					  1.6, 4, 2, 2.2, 3.2,5, 	0, .5, .1, .1, -.2, 3,	rep(-.1, length(selcol)*2+2),	-.5))
# for(i in 1:length(theta_init)){ 
# 	names(theta_init[[i]])	<- c(paste("beta_", 1:nx, sep=""), paste("gamma0_", 1:((R+length(selcol)+1)*2), sep=""), "ln_sigma")
# }
cat("Initial values are:\n"); print(theta_init); cat("\n")
system.time(tmp <- MDCEV_wrapper(theta_init) )

pct <- proc.time()
sol		<- maxLik(MDCEV_wrapper, start=theta_init, method="BFGS", fixed=myfix)
use.time <- proc.time() - pct
cat("MDCEV estimation finishes with", use.time[3]/60, "min.\n")
print(summary(sol))
cat("--------------------------------------------------------\n")

save.image(file = fname)

################################
# Simulate the inclusive value # 
################################
# For simulation, we need elements: income level, expenditure, (average) price, (average) retail attributes
# if(cpi.adj){
# 	lnInc	<- sort(unique(hh_exp$ln_inc))
# 	lnInc	<- lnInc[seq(2, length(lnInc), 2)]
# }else{
# 	lnInc	<- sort(unique(hh_exp$ln_inc))
# 	lnInc	<- c(lnInc + log(.8), lnInc + log(.9), lnInc)
# }
cat("Range of expenditure in the dara:", range(mydata$dol), "\n")
cat("Probability of expenditure > 1000: ", sum(mydata$dol > 1000)/nrow(mydata), "\n")

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

# Subset simulation data
sel			<- duplicated(mydata[,c("household_code", "year")])
sht.X1		<- lapply(X1, function(x) x[!sel,])
sht.X2		<- lapply(X2, function(x) x[!sel,])
sht.price	<- price[!sel,]
sht.y		<- y[!sel]

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

# Fit spline funciton
tmpdat	<- data.frame(omega = colMeans(omega_draw), y = sht.y, ln_inc = ln_inc[!sel])
tmpdat$price	<- price[!sel,]
tmp		<- do.call(cbind, lapply(sht.X1, function(x) x[,(2*(R-1)+1):(2*(R-1)+5)]))
tmpdat$X<- as.matrix(tmp)
omega_deriv	<- spl2dfun(tmpdat, fml = omega ~ te(y, ln_inc) + price + X)
rm(list = c("sht.X1", "sht.X2", "sht.price", "sht.y", "tmpdat"))

save.image(file = fname)

###################################
# Upper level expenditue decision # 
###################################
rec		<- mydata$rec

tmpdat	<- data.frame(y = y, ln_inc = ln_inc)
tmpdat$price	<- price
tmp		<- do.call(cbind, lapply(X1, function(x) x[,(2*(R-1)+1):(2*(R-1)+5)]))
tmpdat$X<- as.matrix(tmp)
o.deriv	<- omega_deriv(tmpdat, deriv = 1)
rm(list = "tmpdat")

M_fn	<- function(lambda, omega_deriv, y, ln_inc, dT = NULL){
# Moment function: lambda * omega'(y,Inc) - 1 = 0
	m1	<- (lambda[1] + lambda[2] * rec)* o.deriv - 1
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

lsol	<- lm(1/o.deriv ~ factor(rec))
pct		<- proc.time()
init	<- matrix(c(coef(lsol), 
					.1, -.2, 
					.05, -.1), 
					3, 2, byrow = T)
colnames(init)	<- c("lambda1", "lambda2")
dT		<- NULL
if(interp.method == "cheb"){
	dT <- cheb.1d.basis(y, numnodes, interval = y.interval)				# Derivative of Chebshev basis
}
system.time(tmp <- GMM_fn(init[1,], omega_deriv, y, ln_inc, dT))

tmp_sol	<- tmp_sol1 <- list(NULL)
for(i in 1:nrow(init)){
	tmp_sol[[i]]	<- maxLik(GMM_fn, start=init[i,], method="BFGS", omega_deriv = omega_deriv, y = y, ln_inc = ln_inc, dT = dT)
	tmp_sol1[[i]]	<- maxLik(GMM_fn, start=init[i,], method="NM", omega_deriv = omega_deriv, y = y, ln_inc = ln_inc, dT = dT)
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
m			<- M_fn(lambda = coef(sol.top), omega_deriv = omega_deriv, y = y , ln_inc = ln_inc, dT = dT)
W			<- solve(var(m$m))
cat("Optimal weighting matrix is:\n"); print(W); cat("\n")
sol.top2	<- maxLik(GMM_fn, start=coef(sol.top), method="BFGS", omega_deriv = omega_deriv, y = y, ln_inc = ln_inc, dT = dT, W = W)
summary(sol.top2)
cat("Top level estimation finishes.\n")

####################
# Save the results #
####################
rm(list= intersect(ls(), c("hh_exp", "Allocation_constr_fn","Allocation_fn","Allocation_nlop_fn","incl_value_fn","i","MDCEV_ll_fnC",
		  "MDCEV_LogLike_fnC","MDCEV_wrapper","tmp","tmp1","tmp2","tmp_sol","sel","sel1","sel2","param_assign","tmpX_list", 
		  "use.time", "pct", "uP_fn","uPGrad_fn", "theta_init", "make_plot", "ord","panelist","tmpidx","tmpn","cl", 
		  "tmp_sol1", "GMM_fn", "M_fn", "init", "m", "mycore", "param_assignR", "ggtmp", "ggtmp1", "tmpdat",
		  "interp.method", "trim.alpha", "mysplfun", "mytrimfun", "expFOC_fn", "exp_fn", "solveExp_fn","dT","i", "j", "k",  
		  "simExp_fn", "SimWrapper_fn", "SimOmega_fn", "cheb.1d.basis", "cheb.basis", "chebfun", "omega_parallel", 
		  "lastFuncGrad", "lastFuncParam", "args", "plot.wd", "tmpv1", "tmpv2")))

save.image(file = fname)

cat("This program is done. ")
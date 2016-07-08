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
library(nloptr)
library(mgcv)
library(zipcode)
options(error = quote({dump.frames(to.file = TRUE)}))

args <- commandArgs(trailingOnly = TRUE)
print(args)
if(length(args)>0){
    for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
    }
}

# seg_id	<- as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))

model_name 		<- "MDCEV_share"
run_id			<- 10
seg_id		<- 1
make_plot		<- TRUE
interp.method	<- "spline"			# "cheb"
trim.alpha		<- 0 #0.05
cpi.adj			<- TRUE
week.price		<- FALSE
cat("seg_id =", seg_id, "\n")

# setwd("~/Documents/Research/Store switching/processed data")
# plot.wd	<- '~/Desktop'
# sourceCpp(paste("../Exercise/Multiple_discrete_continuous_model/", model_name, ".cpp", sep=""))
# source("../Exercise/Multiple_discrete_continuous_model/0_Allocation_function.R")
# source("../Exercise/main/share_allocation_month/ctrfact_sim_functions.r")

# setwd("/home/brgordon/ccv103/Exercise/run")
# setwd("/kellogg/users/marketing/2661703/Exercise/run")
# setwd("/sscc/home/c/ccv103/Exercise/run")
setwd("U:/Users/ccv103/Documents/Research/Store switching/run")

sourceCpp(paste(model_name, ".cpp", sep=""))
source("0_Allocation_function.R")
source("ctrfact_sim_functions.r")

if(cpi.adj){
	fname	<- paste("estrun_",run_id,"/MDCEV_cpi_est_seg",seg_id,"_", Sys.Date(),".rdata",sep="")
}else{
	fname	<- paste("estrun_",run_id,"/MDCEV_est_seg",seg_id,"_", Sys.Date(),".rdata",sep="")
}
fname

#################################
# Read data and subsetting data #
#################################
load("hh_month_exp.rdata")
panelists	<- read.csv("2_panelists.csv", header = T)
penetrat	<- read.csv("iri_channel_penetration.csv", header = T)

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
panelists$zip	<- clean.zipcodes(panelists$panelist_zip_code)
panelists$zip3	<- as.numeric(substr(panelists$zip, 1, 3))
length(unique(panelists$zip3))
length(unique(intersect(panelists$zip3, penetrat$zip3)))
hh_exp		<- merge(hh_exp, panelists[, c("household_code", "zip3")], by = "household_code", all.x=T)
sum(is.na(hh_exp$zip3))
if(sum(is.na(hh_exp$zip3))>0){
	hh_exp	<- subset(hh_exp, !is.na(zip3))
}
cat("Table of segments in the expenditure data:\n"); print(table(hh_exp$first_incomeg)); cat("\n")

# Change the name of channel types 
(tmp	<- setNames(fmt_name, c("Convenience", "Mass", "Dollar_variety", "Drug", "Conventional_supermarkets", "Wholesale_clubs")))
levels(penetrat$channel_type)	<- sapply(levels(penetrat$channel_type), function(x) ifelse(x%in% names(tmp), tmp[x], x))
levels(penetrat$channel_type)

# Subset data 
selcol		<- c("household_code", "dol", "year", "month", "income_midvalue", "ln_inc","first_incomeg", "scantrack_market_descr","zip3",
				paste("SHR_", gsub("\\s", "_", fmt_name), sep=""))
mydata		<- subset(hh_exp[,selcol], as.numeric(first_incomeg) == seg_id & dol > .1 )
mydata$rec	<- 1 * (mydata$year >= 2008)

# Drop the households who did not get matched penetration
tmp		<- unique(panelists$zip3)
sel		<- tmp[!tmp %in% penetrat$zip3] 
mydata	<- subset(mydata, !zip3 %in% sel)

# # Lag of expenditure share
# mydata		<- data.table(mydata)
# for(i in fmt_name){
# 	tmpv1	<- as.name(paste("LAG_", gsub(" ", "_", i), sep=""))
# 	tmpv2	<- as.name(paste("SHR_", gsub(" ", "_", i), sep=""))
# 	mydata	<- mydata[, eval(tmpv1):=c(NA, eval(tmpv2)[-length(year)]), by = list(household_code)]
# }
# 
# # Drop the first observation for each household
# dim(mydata)
# mydata	<- mydata[!is.na(LAG_Grocery),]
# dim(mydata)
# mydata	<- data.frame(mydata)

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

# Accessability matrix
tmp			<- dcast(subset(penetrat, channel_type %in% fmt_name), zip3 ~ channel_type, value.var = "penetration")
acs			<- as.matrix(tmp[,-1])
rownames(acs)	<- tmp$zip3
acs			<- acs[as.character(mydata$zip3),]

# Match retailers' attributes
tmp1	<- unique(fmt_attr$year)
tmp2	<- unique(fmt_attr$scantrack_market_descr)
tmp		<- paste(rep(tmp2, each=length(tmp1)), rep(tmp1, length(tmp2)), sep="-")
tmpn	<- 1:length(tmp)		
names(tmpn) <- tmp
sel 	<- with(mydata, paste(scantrack_market_descr,year, sep="-"))
sel1	<- tmpn[sel]
selcol	<- c("size_index", "ln_upc_per_mod", "ln_num_module","overall_prvt")
nx 		<- (length(selcol)+1)*2 + R - 1
beta0_base 	<- which(fmt_name == "Grocery")

X_list 	<- vector("list", length=length(fmt_name))
for(i in 1:length(fmt_name)){
	sel2		<- fmt_attr$channel_type == fmt_name[i]
	tmp			<- fmt_attr[sel2,selcol]
	tmp1 		<- cbind(access = acs[,i], as.matrix(tmp[sel1,]))
	tmp2		<- matrix(0, nrow(shr), R-1)
	if(i < beta0_base){
		tmp2[,i]	<- mydata$rec
	}else if(i > beta0_base){
		tmp2[,(i-1)] <- mydata$rec
	}
	X_list[[i]]	<- cbind(tmp2, tmp1, mydata$rec*tmp1)
}

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
	MDCEV_ll_fnC(param, nx, shr, y, s1_index, price, X_list, beta0_base)$ll
}

# Estimation with multiple initial values. 
if(fixsigma){
	theta_init	<- list(c(rep(0, R-1), -1, -.2, 1, -2, 	.1, .1, -.1, .1,	-2, -1, -1, -.5, -2,	3, 15, 4, 3, 6, 33, 0), 
						c(rep(0, R-1), -1, -.3, 1, -2, 	.1, .0, -.1, .3, 	-2, -1, -1, -.5, -1, 	3, 15, 4, 4, 7, 37, 0), 
						c(rep(0, R-1), -1, -.5, 1.5, -4,	.1, -.1,-.5, .5, 	-3, -1, -2, -1, -.8, 	3, 15, 4, 4, 8, 40, 0  ) )	
}else if(cpi.adj){
	theta_init	<- list(c(-.2, .03, -.2, -.06, .1, 	1, .2, .1, .2, 1.3, 	rep(-.1, length(selcol)+1), 1, -.2, 1, 1, -2,  		6, 34, 8, 8, 15, 60, -.45), 
						c(-.2, -.05,-.2,-.08,.1,	1, .3, .4, .5, 1, 		rep(-.1, length(selcol)+1), 0, -.6, .8, .5, -1.5, 	6, 35, 8, 8, 18, 70, -.5), 
						c(-.3,-.12,-.38,-.15,.1,	1, .1, .4, .37, -.3, 	rep(-.1, length(selcol)+1), -1, .7, -.5, 1.4, -1.2, 5, 37, 8, 9, 20, 80, -.5 ) )
	
}else{
	theta_init	<- list(c(rep(0, R-1), -1, -.1, 1, -1.5, 	.1, .1, -.1, .1, 	-1, -1, -1, -.5, -1,  5, 20, 4, 4, 6, 50, -.3), 
						c(rep(0, R-1), -.6, -.1, .1, -.8, 	.1, .1, -.05, -.1, 	-1.2, -1, -1, -.3, -.5, 6, 25, 6, 5, 15, 60, -.4), 
						c(rep(0, R-1), -1, -.3, 1.2, -6, 	.1, 0, -.1, .6, 	-1.5, -1, -1.5, -.5, -.3, 5, 30, 6, 6, 15, 80, -.4 ),
						c(rep(0, R-1), -1, -.2, 1.1, -3, 	.1, .05, -.1, .3, 	-1.2, -1, -1.2, -.5, -.6, 6, 25, 5, 6, 10, 65, -.35  ) )
}
cat("Initial values are:\n"); print(theta_init); cat("\n")
system.time(tmp <- MDCEV_wrapper(theta_init[[1]]) )

tmp_sol <- vector("list", length(theta_init))
pct <- proc.time()
for(i in 1:length(theta_init)){
	names(theta_init[[i]]) 	<- c(paste("beta_",1:nx, sep=""), paste("beta0_",setdiff(1:R, beta0_base),sep=""), 
								 paste("gamma_",1:R,sep=""), "ln_sigma")
	tmp_sol[[i]]	<- maxLik(MDCEV_wrapper, start=theta_init[[i]], method="BFGS", fixed=myfix)
	print(i)
}
tmp		<- sapply(tmp_sol, function(x) ifelse(is.null(x$maximum), NA, x$maximum))
cat("Maximium values from different initial values:\n"); print(tmp); cat("\n")
sel 	<- which(abs(tmp - max(tmp, na.rm=T)) < 1e-4, arr.ind=T)
sel1	<- sapply(tmp_sol[sel], function(x) !any(is.na(summary(x)$estimate[,"Std. error"])) & !any(summary(x)$estimate[,"Std. error"]==Inf) )
sel		<- ifelse(sum(sel1)==1, sel[sel1], ifelse(sum(sel1)==0, sel[1], sel[sel1][1]))
# sel 	<- which.max(sapply(tmp_sol, function(x) ifelse(is.null(x$maximum), NA, x$maximum)))
sol		<- tmp_sol[[sel]]
if(!sol$code %in% c(0,1,2,8)){
	cat("MDCEV does NOT have normal convergence.\n")
}
use.time <- proc.time() - pct
cat("MDCEV estimation finishes with", use.time[3]/60, "min.\n")
print(summary(sol))
cat("--------------------------------------------------------\n")

save.image(file = fname)

################################
# Simulate the inclusive value # 
################################
# For simulation, we need elements: income level, expenditure, (average) price, (average) retail attributes
if(cpi.adj){
	lnInc	<- sort(unique(hh_exp$ln_inc))
	lnInc	<- lnInc[seq(2, length(lnInc), 2)]
}else{
	lnInc	<- sort(unique(hh_exp$ln_inc))
	lnInc	<- c(lnInc + log(.8), lnInc + log(.9), lnInc)
}
cat("Range of expenditure in the dara:", range(mydata$dol), "\n")
cat("Probability of expenditure > 1000: ", sum(mydata$dol > 1000)/nrow(mydata), "\n")

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

# Average price
if(week.price){
	tmp 		<- dcast(price_dat,	 scantrack_market_descr + year + biweek ~ channel_type, value.var = "bsk_price_paid_2004") 
	tmp_price	<- rep(1,numnodes) %*% t(colMeans(as.matrix(tmp[,4:(3+R)]), na.rm=T) )
}else{
	tmp 		<- dcast(price_dat,	 scantrack_market_descr + year ~ channel_type, value.var = "bsk_price_paid_2004") 
	tmp_price	<- rep(1,numnodes) %*% t(colMeans(as.matrix(tmp[,3:(2+R)]), na.rm=T) )
}

# Average retail attributes
tmpX_list	<- lapply(fmt_name, function(x) colMeans(as.matrix(subset(fmt_attr, channel_type == x)[,selcol])))		

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
numsim 		<- 1000
eps_draw	<- matrix(rgev(numsim*R, scale = exp(shr.par["ln_sigma"])), numsim, R)
acs_nodes	<- acs[sample(1:nrow(acs), numsim),]

omega_parallel	<- function(eps_draw, acs, rec ){
	out		<- matrix(NA, length(lnInc), numnodes)
	for(j in 1:length(lnInc)){
		tmpX_list1	<- vector("list", R)
		for(k in 1:R){
			tmp	<- rep(0, R-1)
			if(k < beta0_base){	tmp[k] <- rec}
			if(k > beta0_base){ tmp[(k-1)]	<- rec}
			tmp1	<- c(acs[k], tmpX_list[[k]])
			tmpX_list1[[k]]	<- rep(1, numnodes) %*% t(c(tmp, tmp1, tmp1*rec))
		}
		tmpsol 		<- incl_value_fn(param_est=shr.par, nx=nx, base= beta0_base, X_list=tmpX_list1, y=y.nodes, Q=Inf, price=tmp_price, 
					R=R, Ra=R, qz_cons = 0, exp_outside = FALSE, quant_outside = FALSE, eps_draw = rep(1, numnodes) %*% t(eps_draw) )
		out[j,]		<- tmpsol$omega
	}
	return(out)
}

pct			<- proc.time()
tmp			<- foreach(i = 1:numsim, .packages= c("nloptr")) %dopar% {
	omega_parallel(eps_draw[i,], acs = acs_nodes[i,], rec = 0)
}
tmp1		<- foreach(i = 1:numsim, .packages= c("nloptr")) %dopar% {
	omega_parallel(eps_draw[i,], acs = acs_nodes[i,], rec = 1)
}

omega_draw1	<-	omega_draw	<- array(NA, c(numsim, length(lnInc), numnodes), dimnames = list(iter = NULL, lnInc = lnInc, y = y.nodes))
for(i in 1:numsim){
	omega_draw[i,,]	<- tmp[[i]]
	omega_draw1[i,,]<- tmp1[[i]]
}

use.time	<- proc.time() - pct
cat("Simulation of omega with random draws finishes with", use.time[3]/60, "min.\n")
cat("--------------------------------------------------------\n")
stopCluster(cl)
cat("Stop clustering. \n")

# NOTE: we trim 5% omega data to smooth things out. 
if(interp.method == "spline"){
	tmp	<- apply(omega_draw, c(2,3), mytrimfun, alpha = trim.alpha)
	tmpdat	<- melt(tmp)
	names(tmpdat) <- c("iter", "lnInc", "y", "omega")
	tmpdat$Xacs	<- acs_nodes[tmpdat$iter,]						# Note that lags are perfectly correlated 
	omega_deriv	<- spl2dfun(tmpdat, fml = omega ~ te(lnInc, y) + Xacs)

	tmp	<- apply(omega_draw1, c(2,3), mytrimfun, alpha = trim.alpha)
	tmpdat	<- melt(tmp)
	names(tmpdat) <- c("iter", "lnInc", "y", "omega")
	tmpdat$Xacs	<- acs_nodes[tmpdat$iter,]						# Note that lags are perfectly correlated 
	omega_deriv1<- spl2dfun(tmpdat, fml = omega ~ te(lnInc, y) + Xacs)
}else{
	omega_deriv	<- lapply(1:length(lnInc), function(i) chebfun(x = y.nodes, y = tmp[i,], interval = y.interval))
	names(omega_deriv)	<- lnInc
}
rm(list = "tmpdat")

# Check the concavity of omega function 
if(make_plot){
	ggtmp	<- melt(apply(omega_draw, c(2, 3), mean, na.rm = T))
	names(ggtmp)	<- c("lnInc", "y", "omega")
	ggtmp$lnInc		<- round(ggtmp$lnInc, 3)
	
	tmp1	<- seq(80, 120, 1)
	tmpdat	<- data.frame(lnInc = rep(lnInc, each = length(tmp1)), y = rep(tmp1, length(lnInc)))
	tmpdat$Xacs <- rep(1, length(tmp1)*length(lnInc)) %*% t(acs_nodes[1,])
	tmp2	<- omega_deriv(tmpdat)
	ggtmp1	<- cbind(tmpdat, omega = tmp2)
	ggtmp1$lnInc	<- round(ggtmp1$lnInc, 3)
	
	pdf(paste("estrun_",run_id,"/omega_seg",seg_id,"_", Sys.Date(),".pdf",sep=""), width = 10, height = 6)
	print(ggplot(ggtmp, aes(y, omega)) + geom_line() + facet_wrap(~lnInc, scales = "free_y") + 
			xlim(c(0, 1000)) + labs(title = "Omega function at interpolating nodes")
		)
	print(ggplot(ggtmp1, aes(y, omega)) + geom_line() + 
			facet_wrap(~lnInc, scales = "free") + labs(title = "Omega function")
		)
	print(vis.gam(gamfit))
	dev.off()
}

rm(list = c("tmpdat"))
save.image(file = fname)

###################################
# Upper level expenditue decision # 
###################################
o.deriv	<- rep(NA, length(y))
if(interp.method == "cheb"){
	for(i in 1:length(lnInc)){
		sel				<- ln_inc == lnInc[i]
		o.deriv[sel]	<- omega_deriv[[i]](y[sel], deriv = 1, dT = dT[sel,])
	}
}else{
	newx	<- data.frame(lnInc = ln_inc, y = y)
	newx$Xacs<- acs
	sel		<- mydata$rec == 0
	o.deriv[sel]	<- omega_deriv(newx[sel,], deriv = 1)
	o.deriv[!sel]	<- omega_deriv1(newx[!sel,], deriv = 1)
}

M_fn	<- function(lambda, o.deriv, y, ln_inc){
# Moment function: lambda * omega'(y,Inc) - 1 = 0
	m1	<- (lambda[1] + lambda[2] * ln_inc)* o.deriv - 1
	m	<- cbind(m1, m1*ln_inc)
	mbar<- colMeans(m)
	return(list(moment = mbar, m = m))
}

GMM_fn	<- function(lambda, o.deriv, y, ln_inc, W = NULL){
# We use unit diagnial matrix as the weighting matric in GMM
	mm 	<- M_fn(lambda, o.deriv, y, ln_inc)
	m	<- mm$moment
	if(is.null(W)){
		W	<- diag(length(m))
	}
	obj	<- - t(m) %*% W %*% m				# negative moment function for maxLik
	return(obj)
}

lsol	<- lm(1/o.deriv ~ ln_inc)
pct		<- proc.time()
init	<- matrix(c(coef(lsol), 
					.1, -.02, 
					-.01, .01), 
					3, 2, byrow = T)
colnames(init)	<- c("lambda1", "lambda2")
dT		<- NULL
if(interp.method == "cheb"){
	dT <- cheb.1d.basis(y, numnodes, interval = y.interval)				# Derivative of Chebshev basis
}
system.time(tmp <- GMM_fn(init[1,], omega_deriv, y, ln_inc, dT))

tmp_sol	<- tmp_sol1 <- list(NULL)
for(i in 1:nrow(init)){
	tmp_sol[[i]]	<- maxLik(GMM_fn, start=init[i,], method="BFGS", o.deriv = o.deriv, y = y, ln_inc = ln_inc)
	tmp_sol1[[i]]	<- maxLik(GMM_fn, start=init[i,], method="NM", o.deriv = o.deriv, y = y, ln_inc = ln_inc)
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
m			<- M_fn(lambda = coef(sol.top), o.deriv = o.deriv, y = y , ln_inc = ln_inc)
W			<- solve(var(m$m))
cat("Optimal weighting matrix is:\n"); print(W); cat("\n")
sol.top2	<- maxLik(GMM_fn, start=coef(sol.top), method="BFGS", o.deriv = o.deriv, y = y, ln_inc = ln_inc, W = W)
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

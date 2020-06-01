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

seg_id	<- as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))
args <- commandArgs(trailingOnly = TRUE)
print(args)
if(length(args)>0){
    for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
    }
}
cat("seg_id =", seg_id, "\n")

setwd("/sscc/home/c/ccv103/Exercise/run")
run_id			<- 2
load(paste("estrun_",run_id,"/MDCEV_est_seg", seg_id, "_2016-08-28.rdata",sep=""))
rm(list = c("sol.top", "sol.top2"))

model_name 		<- "MDCEV_share_gamma"
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

sourceCpp(paste(model_name, ".cpp", sep=""))
source("0_Efficient_Allocation_function.R")
source("ctrfact_sim_functions.r")

(fname	<- paste("estrun_",run_id,"/MDCEV_est_seg", seg_id, "_", Sys.Date(),".rdata",sep=""))

################################
# Simulate the inclusive value # 
################################
selcol	<- unique(c("distance", selcol))
for(i in 1:R){
	colnames(X1[[i]])	<- c(paste("Intercept:",fmt_name[-beta0_base], sep=""), paste("Rec:",fmt_name[-beta0_base], sep=""), selcol, paste("Rec:",selcol,sep=""))
	colnames(X2[[i]])	<- c(paste("Intercept:",fmt_name, sep=""), paste("Rec:",fmt_name, sep=""), selcol, paste("Rec:",selcol,sep=""))
}

cat("Range of expenditure in the dara:", range(mydata$dol), "\n")
cat("Probability of expenditure > 1000: ", sum(mydata$dol > 1000)/nrow(mydata), "\n")

# Register parallel computing
mycore 	<- 3
cl		<- makeCluster(mycore, type = "FORK")
# cl		<- makeCluster(mycore, type = "PSOCK")			# Register cores in Windows
registerDoParallel(cl)
cat("Register", mycore, "core parallel computing. \n")

#-------------------------------------------- # 
# Compute the inclusive value with random draws
# Subset simulation data
sel			<- duplicated(mydata[,c("household_code", "year")])
sht.X1		<- lapply(X1, function(x) x[!sel,])
sht.X2		<- lapply(X2, function(x) x[!sel,])
sht.price	<- price[!sel,]
sht.y		<- y[!sel]
delta		<- .1

pct			<- proc.time()
omega_draw1	<- foreach(i = 1:numsim, .packages= c("nloptr"), .combine = 'rbind', .errorhandling='remove') %dopar% {
	incl_value_fn(param_est = shr.par, nx, base=beta0_base, X1 = sht.X1, X2= sht.X2, y = sht.y + delta, price = sht.price, R, 
					eps_draw = rep(1,length(sht.y)) %*% t(eps_draw[i,]), returnmax = TRUE)$omega
}
use.time	<- proc.time() - pct
cat("Simulation of omega with random draws finishes with", use.time[3]/60, "min.\n")
cat("--------------------------------------------------------\n")
stopCluster(cl)
cat("Stop clustering. \n")

if(trim.alpha == 0){
	tmpdat	<- data.frame(omega = (colMeans(omega_draw1) - colMeans(omega_draw))/delta, y = sht.y, ln_inc = ln_inc[!sel])
}else{
	tmp1	<- apply(omega_draw, 2, mean, trim = trim.alpha, na.rm=T)
	tmp2	<- apply(omega_draw1, 2, mean, trim = trim.alpha, na.rm=T)
	tmpdat	<- data.frame(omega = (tmp2-tmp1)/delta, y = sht.y, ln_inc = ln_inc[!sel])
}
summary(tmpdat$omega)
tmpdat$price	<- sht.price
tmpdat$X		<- as.matrix(do.call(cbind, lapply(sht.X1, function(x) x[,selcol])))	
omega_deriv		<- spl2dfun(tmpdat, fml = omega ~ y + I(y^2) + price + X)

# Plot omega function
(tmp		<- quantile(y, c(.1, .99)))
tmpdat1	<- data.frame(y = seq(tmp[1], tmp[2], length = 100))
tmpdat1$price	<- rep(1, nrow(tmpdat1)) %*% t(colMeans(tmpdat$price))
tmpdat1$X<- rep(1, nrow(tmpdat1)) %*% t(colMeans(tmpdat$X))
ggtmp	<- data.frame(y = tmpdat1$y, d1 = omega_deriv(tmpdat1), d2 = omega_deriv(tmpdat1, deriv = 1))
ggtmp	<- melt(ggtmp, id.vars = "y")

pdf(gsub(".rdata", ".pdf",fname), width = 6, height = 4.5)
ggplot(ggtmp, aes(y, value)) + geom_line() + facet_wrap(~variable, scales = "free")
dev.off()

###################################
# Upper level expenditue decision # 
###################################
rec		<- mydata$rec

tmpdat	<- data.frame(y = y, ln_inc = ln_inc)
tmpdat$price	<- price
tmp		<- do.call(cbind, lapply(X1, function(x) x[,selcol]))
tmpdat$X<- as.matrix(tmp)
o.deriv	<- omega_deriv(tmpdat)
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
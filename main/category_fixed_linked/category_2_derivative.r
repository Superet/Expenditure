
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
cat("seg_id =",seg_id, ".\n")

interp.method	<- "spline"			# "cheb"
trim.alpha		<- 0.05
# seldpt			<- "HBC"	#"NFG" #"DG" 

# setwd("~/Documents/Research/Store switching/Processed_data")
# plot.wd	<- '~/Desktop'
# sourceCpp(paste("../Exercise/Multiple_discrete_continuous_model/", model_name, ".cpp", sep=""))
# source("../Exercise/Multiple_discrete_continuous_model/0_Efficient_Allocation_function.R")
# source("../Exercise/main/ctrfact_sim_functions.r")

# setwd("/home/brgordon/ccv103/Exercise/run")
# setwd("/kellogg/users/marketing/2661703/Exercise/run")
setwd("/sscc/home/c/ccv103/Exercise/run")
# setwd("U:/Users/ccv103/Documents/Research/Store switching/run")

source("0_Efficient_Allocation_function.R")
source("ctrfact_sim_functions.r")

ver.date <- "2016-10-01"
# seldpt	<- "DG"
(fname	<- paste("estrun_4/fixb_est_", seldpt, "_seg", seg_id, "_", ver.date, ".rdata",sep=""))
load(fname)
rm(list = intersect(ls(), c("omega_draw", "omega_deriv", "fname", "lsol","beta.init","gamma.init", 
			"inits", "inits2","o.deriv","W","ver.date","sol.top","sol.top2","theta.init")))
run_id			<- 5			
(fname	<- paste("estrun_",run_id,"/fixb_est_", seldpt, "_seg", seg_id, "_", Sys.Date(),".rdata",sep=""))

# Subset simulation data
sel			<- duplicated(mydata[,c("household_code", "year")]) & mydata$dol <= quantile(mydata$dol, .95)
sht.X1		<- lapply(X1, function(x) x[!sel,])
sht.X2		<- lapply(X2, function(x) x[!sel,])
sht.price	<- price[!sel,]
sht.y		<- y[!sel]
delta		<- .1

pct			<- proc.time()
omega_draw	<- foreach(i = 1:numsim, .packages= c("nloptr"), .combine = 'rbind', .errorhandling='remove') %do% {
	incl_value_fn(param_est = shr.par, nx, base=beta0_base, X1 = sht.X1, X2= sht.X2, y = sht.y, price = sht.price, R, 
					eps_draw = rep(1,length(sht.y)) %*% t(eps_draw[i,]), returnmax = TRUE)$omega
}

omega_draw1	<- foreach(i = 1:numsim, .packages= c("nloptr"), .combine = 'rbind', .errorhandling='remove') %do% {
	incl_value_fn(param_est = shr.par, nx, base=beta0_base, X1 = sht.X1, X2= sht.X2, y = sht.y+delta, price = sht.price, R, 
					eps_draw = rep(1,length(sht.y)) %*% t(eps_draw[i,]), returnmax = TRUE)$omega
}

use.time	<- proc.time() - pct
cat("Simulation of omega with random draws finishes with", use.time[3]/60, "min.\n")
cat("--------------------------------------------------------\n")
print(summary(c(omega_draw)))
# stopCluster(cl)
# cat("Stop clustering. \n")

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
# omega_deriv		<- spl2dfun(tmpdat, fml = omega ~ s(y) + price)
omega_deriv		<- spl2dfun(tmpdat, fml = omega ~ y + I(y^2) + price)

# Plot omega function
tmp		<- dcast(fmt_dpt[fmt_dpt$department == dpt_name[seldpt],c("scantrack_market_descr", "year", "channel_type","unitprice_paid")], 
				scantrack_market_descr + year ~ channel_type, value.var = "unitprice_paid")
tmp1	<- as.vector(colMeans(tmp[,fmt_name], na.rm=T))
tmp		<- attr_mat[,c("household_code", "year", "channel_type", selcol)]
tmp2	<- unlist(lapply(fmt_name, function(x) colMeans(tmp[tmp$channel_type==x, selcol], na.rm=T)))
(tmp		<- quantile(y, c(.1, .99)))
tmpdat	<- data.frame(y = seq(tmp[1], tmp[2], length = 100))
tmpdat$price	<- rep(1, nrow(tmpdat)) %*% t(tmp1)
tmpdat$X<- rep(1, nrow(tmpdat)) %*% t(tmp2)
ggtmp	<- data.frame(y = tmpdat$y, d1 = omega_deriv(tmpdat), d2 = omega_deriv(tmpdat, deriv = 1))
ggtmp	<- melt(ggtmp, id.vars = "y")

pdf(gsub(".rdata", ".pdf",fname), width = 6, height = 6)
ggplot(ggtmp, aes(y, value)) + geom_line() + facet_wrap(~variable, scales = "free")
dev.off()

rm(list = intersect(ls(), c("args","tmpdat", "shr.x", "sht.y", "sht.X1", "sht.X2", "sht.price", "sel", "pct", "use.time", "i", 
			"Allocation_constr_fn","Allocation_fn","Allocation_nlop_fn","incl_value_fn", "omega.lm", "mytrimfun", "ggtmp","simExp_fn",
			"exp_fn", "expFOC_fn", "f", "SimWrapper_fn","solveExp_fn", "spl2dfun", "tmp", "tmp1", "tmp2", "uP_fn","uPGrad_fn", "mysplfun"))
	)

save.image(file = fname)

cat("This program is done. \n")



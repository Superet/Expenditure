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
library(chebpol)
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

# setwd("~/Documents/Research/Store switching/processed data")
# plot.wd	<- '~/Desktop'
# sourceCpp(paste("../Exercise/Multiple_discrete_continuous_model/", model_name, ".cpp", sep=""))
# source("../Exercise/Multiple_discrete_continuous_model/0_Allocation_function.R")

setwd("/home/brgordon/ccv103/Exercise/run")
# setwd("/kellogg/users/marketing/2661703/Exercise/run")
# setwd("/sscc/home/c/ccv103/Exercise/run")
model_name 	<- "MDCEV_share"
run_id		<- 3
# seg_id		<- 1
make_plot	<- TRUE
plot.wd		<- paste(getwd(), "/estrun_",run_id, sep="")
ww			<- 10
ar			<- .6

ver.date	<- "2016-02-16"
load(paste("estrun_",run_id,"/MDCEV_est_seg",seg_id,"_", ver.date,".rdata",sep=""))
sourceCpp(paste(model_name, ".cpp", sep=""))
source("0_Allocation_function.R")
source("ctrfact_sim_functions_v2.r")

rm(list = c("gamfit","model_name", "tmpdat"))

#######################
# Robust otpimization #
#######################
MDCEV_wrapper <- function(param){
	MDCEV_ll_fnC(param, nx, shr, y, s1_index, price, X_list, beta0_base)$ll
}

# Continue maximization with SANN, using the best coef as starting values
pct			<- proc.time()
sol.sann	<- maxLik(MDCEV_wrapper, start=coef(sol), method="SANN", fixed=myfix)
use.time	<- proc.time() - pct
cat("The estimation with SANN finishes with", use.time[3]/60, "min.\n")
print(summary(sol.sann)); cat("\n")
cat("The difference of coefficient estiamtes between sol and SANN =", max(abs(coef(sol) - coef(sol.sann))), ".\n")

# Check log likelihood concavity around the estimates
FindOrder	<-function(x){
	sig.x	<- format(signif(abs(x), 0), scientific = FALSE)
	if(sig.x < 1){
		x1	<- gsub("(.*)(\\.)|([0]*$)","", sig.x)
		d	<- -nchar(x1)
	}else{
		d	<- nchar(as.character(trunc(as.numeric(sig.x)))) - 1
	}
	return(sign(x)*10^d)
}

my.coef		<- coef(sol.sann)
ggtmp		<- data.frame()
my.grid		<- -10:10
if(length(myfix) > 0){
	coef.grid	<- sapply(my.coef[-myfix], function(x) FindOrder(x) * my.grid) + rep(1, length(my.grid)) %*% t(my.coef[-myfix])
}else{
	coef.grid	<- sapply(my.coef, function(x) FindOrder(x) * my.grid) + rep(1, length(my.grid)) %*% t(my.coef)
}

for(i in 1:ncol(coef.grid)){
	for(j in 1:nrow(coef.grid)){
		tmp		<- my.coef
		tmp[colnames(coef.grid)[i]] <- coef.grid[j,i]
		ll		<-  sum(MDCEV_wrapper(tmp) )
		ggtmp	<- rbind(ggtmp, data.frame(Coef = colnames(coef.grid)[i], value = coef.grid[j,i], ll = ll))
	}
}
ggtmp$est	<- my.coef[as.character(ggtmp$Coef)]
ggtmp$est.ind	<- ifelse(ggtmp$value == ggtmp$est, 1, 0)

if(make_plot){
	pdf(paste(plot.wd, "/graph_init_ll_seg", seg_id, ".pdf",sep=""), width = 8, height = 8)
	print(ggplot(ggtmp, aes(value, ll)) + geom_point(aes(color = factor(est.ind))) + 
			geom_line() + 
			scale_color_manual(values = c("black", "red")) +
			facet_wrap(~Coef, scales = "free") + 
			guides(color = FALSE)
		)
	dev.off()
}

# Save results 
rm(list=c("hh_exp", "Allocation_constr_fn","Allocation_fn","Allocation_nlop_fn","incl_value_fn","i","MDCEV_ll_fnC",
		  "MDCEV_LogLike_fnC","MDCEV_wrapper","tmp","tmp1","tmp2","param_assign",
		  "use.time", "pct", "uP_fn","uPGrad_fn", "make_plot", "ord","panelist","tmpidx","tmpn","cl", 
		  "tmp_sol1", "GMM_fn", "M_fn", "init", "mycore", "param_assignR", "ggtmp", "ggtmp1", 
		  "mysplfun", "mytrimfun", "expFOC_fn", "exp_fn", "solveExp_fn", 
		  "simExp_fn", "SimWrapper_fn", "SimOmega_fn", "cheb.1d.basis", "cheb.basis", "chebfun", "omega_parallel", 
		  "lastFuncGrad", "lastFuncParam", "args"))
save.image(paste("estrun_",run_id,"/MDCEV_bottom_test_seg",seg_id, "_", Sys.Date(),".rdata",sep=""))

cat("This program is done. ")



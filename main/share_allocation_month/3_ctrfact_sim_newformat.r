library(ggplot2)
library(reshape2)
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
seg_id 	<- 1
cat("seg_id =", seg_id, ".\n")

# setwd("~/Documents/Research/Store switching/processed data/Estimation")
# plot.wd	<- '~/Desktop'
# source("../../Exercise/Multiple_discrete_continuous_model/0_Allocation_function.R")
# source("../../Exercise/main/share_allocation/ctrfact_sim_functions.r")

# setwd("/home/brgordon/ccv103/Exercise/run")
# setwd("/kellogg/users/marketing/2661703/Exercise/run")
# setwd("/sscc/home/c/ccv103/Exercise/run")
# setwd("E:/Users/ccv103/Documents/Research/Store switching/run")
setwd("U:/Users/ccv103/Documents/Research/Store switching/run")

run_id		<- 9

plot.wd		<- getwd()
make_plot	<- TRUE
ww			<- 10
ar			<- .6

source("0_Allocation_function.R")
source("ctrfact_sim_functions.r")

# Load estimation data 
ver.date	<- "2016-06-20"
cpi.adj		<- TRUE

if(cpi.adj){
	loadf	<- paste("estrun_",run_id,"/MDCEV_cpi_est_seg",seg_id,"_", ver.date,".rdata",sep="")
}else{
	loadf	<- paste("estrun_",run_id,"/MDCEV_est_seg",seg_id,"_", ver.date,".rdata",sep="")
}
loadf
load(loadf)
rm(list = intersect(ls(), c("gamfit", "shr","model_name", "tmpdat")))

# Set simulation parameters
week.price		<- FALSE
interp.method	<- "spline"					# Spline interpolation or Chebyshev interpolation
exp.method		<- "Utility"				# Utility maximization or finding roots for first order condition
trim.alpha		<- 0.05
numsim			<- 1000	#numsim1		#<- 1000
draw.par		<- FALSE
sim.omega		<- FALSE
fname			<- paste("ctrfact_newf_seg",seg_id,sep="")
if(draw.par){ fname <- paste(fname, "_pardraw", sep="") }
if(sim.omega){fname	<- paste(fname, "_simomega", sep="") }
fname			<- paste(fname, "_sim", numsim, "_", as.character(Sys.Date()), sep="")
cat("Output file name is", fname, ".\n")

###############################
# Prepare simulation elements # 
###############################
# Data required: parameters, income level, price, retail attributes, and random draws
# For each simulation scenario, if income, price or retail attributes change, we need to re-simulate inclusive values.

# Set simulation parameters
lambda 	<- coef(sol.top2)
shr.par	<- coef(sol)

#-----------------------#
# Construct income data # 
# Take the households' income in 2007 as basis 
selyr	<- 2007
tmp		<- data.table(subset(mydata, year %in% selyr))
tmp		<- tmp[,list(income = unique(income_midvalue)), by = list(household_code, year)]
sim.data<- data.frame(tmp)[,c("household_code","income")]
sim.unq		<- data.frame(income2007 = unique(sim.data[,-1]))

# Counterfactual scenario: income is lower by 10%. 
my.change		<- .1
lnInc_08		<- lnInc + log(1 - my.change)
sim.unq$income2008	<- (1 - my.change) * sim.unq$income2007		
sim.unq$Inc07	<- log(sim.unq[,"income2007"])
sim.unq$Inc08	<- log(sim.unq[,"income2008"])
sim.unq			<- sim.unq[order(sim.unq$income2007),]
cat("dim(sim.unq) =", dim(sim.unq), "\n")

#----------------------------#
# Average price in year 2007
if(week.price){
	tmp 		<- dcast(subset(price_dat, year %in% selyr),	 
						scantrack_market_descr + year + biweek ~ channel_type, value.var = "bsk_price_paid_2004") 
	price0	<- setNames( colMeans(as.matrix(tmp[,4:(3+R)]), na.rm=T), fmt_name)
}else{
	tmp 		<- dcast(subset(price_dat, year %in% selyr),	 
						scantrack_market_descr + year ~ channel_type, value.var = "bsk_price_paid_2004") 
	price0	<- setNames( colMeans(as.matrix(tmp[,3:(2+R)]), na.rm=T), fmt_name)
}
cat("The average price level in 2007:\n"); print(price0);cat("\n")

# Average retail attributes in 2007
selcol	<- c("size_index", "ln_upc_per_mod", "ln_num_module","overall_prvt")
X_ls0	<- setNames( lapply(fmt_name, function(x) colMeans(as.matrix(subset(fmt_attr, channel_type == x & year%in% selyr)[,selcol]))), 
					 fmt_name)
cat("The average retail attributes in 2007:\n"); print(do.call(rbind, X_ls0)); cat("\n")

# The average shares for each income level 
tmp		<- as.matrix(subset(mydata, year %in% selyr)[,paste("SHR_",gsub(" ", "_", fmt_name), sep="")])
tmp		<- tmp * subset(mydata, year %in% selyr)$dol
tmp1	<- subset(mydata, year %in% selyr)$income_midvalue
# # Share by income level 
# shr.stat<- apply(tmp, 2, function(x) tapply(x, tmp1, sum))
# shr.stat<- shr.stat/rowSums(shr.stat)
# Average share
shr.stat<- colSums(tmp)
shr.stat<- rep(1, nrow(sim.unq)) %*% t(shr.stat/sum(shr.stat))

# Expand X_list and price to match the nobs of income 
price.07	<- rep(1, nrow(sim.unq)) %*% matrix(price0, nrow = 1)
colnames(price.07)	<- fmt_name
X_list08	<- X_list07	<- setNames(vector("list", R), fmt_name)
for(i in 1:R){
	tmp2	<- tmp1	<- matrix(0, nrow(sim.unq), R-1)
	if(i<beta0_base){
		tmp1[,i]	<- sim.unq$Inc07
		tmp2[,i]	<- sim.unq$Inc08
	}else if(i>beta0_base){
		tmp1[,(i-1)]	<- sim.unq$Inc07
		tmp2[,(i-1)]	<- sim.unq$Inc08
	}
	X_list07[[i]]	<- cbind(tmp1, lag = shr.stat[,i], rep(1, nrow(sim.unq)) %*% matrix(X_ls0[[i]], nrow = 1))
	X_list08[[i]]	<- cbind(tmp2, lag = shr.stat[,i], rep(1, nrow(sim.unq)) %*% matrix(X_ls0[[i]], nrow = 1))
	colnames(X_list07[[i]])	<- colnames(X_list08[[i]])	<- c(fmt_name[-beta0_base], "lag",selcol)
}

#-------------------#
# Take random draws #
set.seed(666)
eps_draw	<- matrix(rgev(numsim*R, scale = exp(shr.par["ln_sigma"])), numsim, R)
if(draw.par){
	par_se	<- c(sqrt(diag(vcov(sol.top2))), sqrt(diag(vcov(sol))) )
	par_se[is.na(par_se)]	<- 0
	par_draw <- sapply(par_se, function(x) rnorm(numsim, mean = 0, sd = x))
}else{
	par_draw <- NULL
}

##############
# Simulation #
##############
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

# Simulate expenditure and expenditure share. 
pct				<- proc.time()
sim.base07		<- SimWrapper_fn(omega_deriv, ln_inc = sim.unq$Inc07, lambda = lambda, param_est = shr.par, base = beta0_base, 
					X_list = X_list07, price = price.07, eps_draw = eps_draw, method = exp.method, ret.sim = TRUE, 
					par.draw = par_draw, share.state = shr.stat[,-beta0_base])
use.time		<- proc.time() - pct						
cat("2007 Baseline simulation finishes with", use.time[3]/60, "min.\n")

pct				<- proc.time()
sim.base08		<- SimWrapper_fn(omega_deriv, ln_inc = sim.unq$Inc08, lambda = lambda, param_est = shr.par, base = beta0_base, 
					X_list = X_list08, price = price.07, eps_draw = eps_draw, method = exp.method, ret.sim = TRUE, 
					par.draw = par_draw, share.state = shr.stat[,-beta0_base])
use.time		<- proc.time() - pct						
cat("2008 Baseline simulation finishes with", use.time[3]/60, "min.\n")

# -------------- #
# Counterfactual #
# Introduct a new small format similar to convenience stores
# Scenarios:
# a(2): if intercept inherents the grocery stores, discount stores;
# b(2): if price inherents the intercept
 
ret.a	<- c("Discount Store")
inherent.price	<- c(TRUE, FALSE)
newf.sim07	<- newf.sim <- setNames(vector("list", length(ret.a)*2), paste(rep(ret.a, each =2), "_prc", c(1,0), sep = ""))
eps_draw_new	<- cbind(eps_draw, rgev(numsim, scale = exp(shr.par["ln_sigma"])) )

for(i in 1:length(ret.a)){
	cat(" ---------------------------------------------------------------------- \n")
	pct		<- proc.time()
	# Set new parameters
	sel		<- which(fmt_name == ret.a[i])
	beta1_7	<- ifelse(sel==beta0_base, 0, ifelse(sel<beta0_base, shr.par[paste("beta_", sel, sep="")], shr.par[paste("beta_", sel-1, sep="")]))
	beta0_7	<- ifelse(sel==beta0_base, 0, shr.par[paste("beta0_", sel,sep="")])
	shr.par.new	<- c(shr.par[1:(R-1)], beta1_7, shr.par[R:nx], shr.par[grep("beta0", names(shr.par))], beta0_7, 
					shr.par[grep("gamma", names(shr.par))], shr.par[paste("gamma_", sel, sep="")], shr.par[grep("sigma",names(shr.par))])
	cat("The parameter vector for adding a new format from", ret.a[i], "is:\n"); print(shr.par.new); cat("\n")
	
	for(j in 1:length(inherent.price)){
		# Set new attributes
		X.new08	<- X.new07	<- setNames(vector("list", R+1), c(fmt_name, "new"))
		for(k in 1:R){
				tmp2	<- tmp1	<- matrix(0, nrow(sim.unq), R)
				if(k<beta0_base){
					tmp1[,k]	<- sim.unq$Inc07
					tmp2[,k]	<- sim.unq$Inc08
				}else if(k>beta0_base){
					tmp1[,(k-1)]	<- sim.unq$Inc07
					tmp2[,(k-1)]	<- sim.unq$Inc08
				}
				X.new07[[k]]	<- cbind(tmp1, lag = shr.stat[,k], rep(1, nrow(sim.unq)) %*% matrix(X_ls0[[k]], nrow = 1) )
				X.new08[[k]]	<- cbind(tmp2, lag = shr.stat[,k], rep(1, nrow(sim.unq)) %*% matrix(X_ls0[[k]], nrow = 1) )
				colnames(X.new07[[k]])	<- colnames(X.new08[[k]])	<- c(fmt_name[-beta0_base], "new", "lag", selcol)
			}
		tmp2	<- tmp1	<- matrix(0, nrow(sim.unq), R)
		if(sel!=beta0_base){
			tmp1[,R]	<- sim.unq$Inc07
			tmp2[,R]	<- sim.unq$Inc08
		}
		tmp		<- X_ls0[["Convenience Store"]]
		tmp["overall_prvt"]	<- X_ls0[[ret.a[i]]]["overall_prvt"]
		X.new07[[(R+1)]]	<- cbind(tmp1, lag = 0, rep(1, nrow(sim.unq)) %*% matrix(tmp, nrow = 1))
		X.new08[[(R+1)]]	<- cbind(tmp2, lag = 0, rep(1, nrow(sim.unq)) %*% matrix(tmp, nrow = 1))

		# Set prices 
		price.new		<- price.07
		if(inherent.price[j]){
			price.new	<- cbind(price.new, new = price.07[,ret.a[i]])
		}else{
			price.new	<- cbind(price.new, new = price.07[,"Convenience Store"])
		}
		
		sidx	<- (i-1)*2 + j
		newf.sim07[[sidx]]	<- SimWrapper_fn(omega_deriv = omega_deriv, ln_inc = sim.unq[,"Inc07"], 
					lambda = lambda, param_est = shr.par.new, 
					base = beta0_base, X_list = X.new07, price = price.new, sim.y = sim.base07$y, 
					eps_draw = eps_draw_new, method = exp.method, ret.sim = TRUE, par.draw = par_draw, share.state = shr.stat[,-beta0_base])
		newf.sim[[sidx]]	<- SimWrapper_fn(omega_deriv = omega_deriv, ln_inc = sim.unq[,"Inc08"], 
					lambda = lambda, param_est = shr.par.new, 
					base = beta0_base, X_list = X.new08, price = price.new, sim.y = sim.base08$y, 
					eps_draw = eps_draw_new, method = exp.method, ret.sim = TRUE, par.draw = par_draw, share.state = shr.stat[,-beta0_base])
	}
	use.time	<- proc.time() - pct
	cat("Counterfactual with", ret.a[i], "finishes with", use.time[3]/60, "min.\n")
}

# Compare the difference 
fmt_name1	<- c(fmt_name, "new")
inc.tab	<- table(mydata[mydata$year==2007,c("income_midvalue")])
tmp1	<- sim.base08$Average
tmp2	<- newf.sim[["Discount Store_prc1"]]$Average
mkt1	<- c(sapply(1:R, function(i) sum(tmp1[,fmt_name[i]]*inc.tab)), 0)
mkt2	<- sapply(1:(R+1), function(i) sum(tmp2[,fmt_name1[i]]*inc.tab))
mkt.dif	<- setNames(mkt2/sum(mkt2)-mkt1/sum(mkt1), fmt_name1)
cat("Difference of market share:\n"); print(round(mkt.dif*100,2)); cat("\n")

tmp1	<- cbind(tmp1[, fmt_name]/tmp1[,"Exp"], 0)
tmp2	<- tmp2[, c(fmt_name, "new")]/tmp2[,"Exp"]
tmp		<- tmp2 - tmp1
tmp.dif	<- apply(tmp, 2, function(x) weighted.mean(x, w = inc.tab))
cat("Mean of share difference is:\n"); print(round(tmp.dif*100, 2)); cat("\n")

################
# Save results #
################
rm(list  = intersect(ls(), c("ar", "args", "cl", "ggtmp", "ggtmp1", "ggtmp2", "i", "lastFuncGrad", "lastFuncParam", "make_plot", "mycore", 
			"myfix", "plot.wd", "s1_index", "sel", "selyr", "price", "sol", "sol.top", "sol.top2", "tmp", "tmp_coef", 
			"tmp_price", "tmp1", "tmp2", "use.time", "ver.date", "var1", "ww", "f", "dT","lag_nodes","sidx","price0", "X_ls0",
			"numnodes", "out", "out1", "pct", "tmpd1", "tmpd2", "tmpdat", "u", "W", "y", "y.nodes", "tmp.dif", "inc.tab", 
			"Allocation_constr_fn", "Allocation_fn", "Allocation_nlop_fn", "cheb.1d.basis", "cheb.basis", "chebfun", 
			"exp_fn", "expFOC_fn", "incl_value_fn", "mysplfun", "mytrimfun", "param_assignR", "simExp_fn", "SimOmega_fn", "tmpX",
			"SimWrapper_fn", "solveExp_fn", "spl2dfun", "uP_fn", "uPGrad_fn", "X_list", "beta0_7", "d.psi", "gamfit","j", "k")))

save.image(paste("estrun_",run_id, "/", fname, ".rdata",sep=""))

cat("This program is done.\n")

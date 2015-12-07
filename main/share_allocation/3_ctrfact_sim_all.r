library(ggplot2)
library(reshape2)
library(Rcpp)
library(RcppArmadillo)
library(maxLik)
library(evd)
library(data.table)
library(doParallel)
library(foreach)
library(chebpol)
options(error = quote({dump.frames(to.file = TRUE)}))

args <- commandArgs(trailingOnly = TRUE)
print(args)
if(length(args)>0){
    for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
    }
}

model_name 	<- "MDCEV_share"
# setwd("~/Documents/Research/Store switching/processed data")
# plot.wd	<- '~/Desktop'
# source("../Exercise/Multiple_discrete_continuous_model/0_Allocation_function.R")
# source("../Exercise/main/share_allocation/ctrfact_sim_functions.r")

setwd("/home/brgordon/ccv103/Exercise/run")
# setwd("/kellogg/users/marketing/2661703/Exercise/run")
# setwd("/sscc/home/c/ccv103/Exercise/run")
run_id		<- 1
plot.wd		<- getwd()
make_plot	<- TRUE
ww			<- 6.5
ar			<- .6

source("0_Allocation_function.R")
source("ctrfact_sim_functions.r")

# Load estimation data 
ver.date	<- "2015-11-25"
load(paste("estrun_",run_id,"/MDCEV_est_seg",seg_id,"_", ver.date,".rdata",sep=""))
interp.method	<- "spline"					# Spline interpolation or Chebyshev interpolation
exp.method		<- "Utility"				# Utility maximization or finding roots for first order condition

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
tmp 		<- dcast(subset(price_dat, year %in% selyr),	 
					scantrack_market_descr + year + biweek ~ channel_type, value.var = "bsk_price_paid_2004") 
price.07	<- setNames( colMeans(as.matrix(tmp[,4:(3+R)]), na.rm=T), fmt_name)
cat("The average price level in 2007:\n"); print(price.07);cat("\n")

# Average retail attributes in 2007
selcol	<- c("size_index", "ln_upc_per_mod", "ln_num_module","overall_prvt")
X_list07<- setNames( lapply(fmt_name, function(x) colMeans(as.matrix(subset(fmt_attr, channel_type == x & year%in% selyr)[,selcol]))), 
					 fmt_name)
cat("The average retail attributes in 2007:\n"); print(do.call(rbind, X_list07)); cat("\n")

#-------------------#
# Take random draws #
set.seed(666)
numsim 		<- 400
eps_draw	<- matrix(rgev(numsim*R), numsim, R)

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

# Register parallel computing
mycore 	<- 2
cl		<- makeCluster(mycore, type = "FORK")
registerDoParallel(cl)
cat("Register", mycore, "core parallel computing. \n")

# --------------------#
# Simulate inclusive values 
omega_parallel	<- function(eps_draw, lnInc, y.nodes, X_list, price){
	numnodes<- length(y.nodes)
	out		<- matrix(NA, length(lnInc), numnodes)
	for(j in 1:length(lnInc)){
		tmpX_list1	<- lapply(X_list, function(x) rep(1, numnodes) %*% t(c(x, x*lnInc[j])))
		tmpsol 		<- incl_value_fn(param_est=shr.par, base= beta0_base, X_list=tmpX_list1, y=y.nodes, Q=Inf, price=rep(1, numnodes) %*% t(price), 
							R=R, Ra=R, qz_cons = 0, exp_outside = FALSE, quant_outside = FALSE, eps_draw = rep(1, numnodes) %*% t(eps_draw) )
		out[j,]		<- tmpsol$omega
	}
	return(out)
}

pct		<- proc.time()
tmp		<- foreach(i = 1:numsim) %dopar% {
	omega_parallel(eps_draw = eps_draw[i,], lnInc = lnInc, y.nodes = y.nodes, X_list = X_list07, price = price.07)
}
tmp1	<- foreach(i = 1:numsim) %dopar% {
	omega_parallel(eps_draw = eps_draw[i,], lnInc = lnInc_08, y.nodes = y.nodes, X_list = X_list07, price = price.07)
}
omega.07	<- array(NA, c(numsim, length(lnInc), numnodes), dimnames = list(NULL, lnInc, y.nodes))
omega.08	<- array(NA, c(numsim, length(lnInc_08), numnodes), dimnames = list(NULL, lnInc_08, y.nodes))
for(i in 1:numsim){
	omega.07[i,,]	<- tmp[[i]]
	omega.08[i,,]	<- tmp1[[i]]
}
use.time	<- proc.time() - pct
cat("Simulation of omega with random draws finishes with", use.time[3]/60, "min.\n")

# Interpolation functions
if(interp.method == "spline"){
	# omega.ls.07 <- lapply(1:length(lnInc), function(i) splinefun(y.nodes, omega.07[i,], method = "natural"))
	# omega.ls.08 <- lapply(1:length(lnInc), function(i) splinefun(y.nodes, omega.08[i,], method = "natural"))
	tmp1	<- apply(omega.07, c(2,3), mytrimfun)
	tmp2	<- apply(omega.08, c(2,3), mytrimfun)
	omega.ls.07	<- lapply(1:length(lnInc), function(i) 
								mysplfun(x = rep(y.nodes, each = dim(tmp1)[1]), y = as.vector(tmp1[,i,])))
	omega.ls.08	<- lapply(1:length(lnInc), function(i) 
								mysplfun(x = rep(y.nodes, each = dim(tmp2)[1]), y = as.vector(tmp2[,i,])))								
}else{
	omega.07	<- apply(omega.07, c(2, 3), mean, trim = .05)
	omega.08	<- apply(omega.08, c(2, 3), mean, trim = .05)
	omega.ls.07 <- lapply(1:length(lnInc), function(i) chebfun(x = y.nodes, y = omega.07[i,], interval = y.interval))
	omega.ls.08 <- lapply(1:length(lnInc), function(i) chebfun(x = y.nodes, y = omega.08[i,], interval = y.interval))
}
names(omega.ls.07)	<- lnInc
names(omega.ls.08)	<- lnInc_08

# Simulate expenditure and expenditure share. 
pct				<- proc.time()
sim.base07		<- SimWrapper_fn(omega.ls.07, ln_inc = sim.unq$Inc07, lambda = lambda, param_est = shr.par, 
						base = beta0_base, X_list = X_list07, price = price.07, eps_draw = eps_draw, method = exp.method)
use.time		<- proc.time() - pct						
cat("2007 Baseline simulation finishes with", use.time[3]/60, "min.\n")

pct				<- proc.time()
sim.base08		<- SimWrapper_fn(omega_fn_ls = omega.ls.08, ln_inc = sim.unq$Inc08, lambda = lambda, param_est = shr.par, 
						base = beta0_base, X_list = X_list07, price = price.07, eps_draw = eps_draw, method = exp.method)
use.time		<- proc.time() - pct						
cat("2008 Baseline simulation finishes with", use.time[3]/60, "min.\n")

# Check the baseline simulation 
tmp		<- seq(50, 250, 10)
tmp1	<- sapply(omega.ls.07, function(f) f(tmp))
tmp2	<- sapply(omega.ls.08, function(f) f(tmp))
ggtmp	<- rbind(data.frame(melt(tmp1), year = 2007), data.frame(melt(tmp2), year = 2008))
names(ggtmp)	<- c("yidx", "lnInc", "omega", "year")
ggtmp$y			<- tmp[ggtmp$yidx]
ggtmp$Income	<- exp(ggtmp$lnInc)

tmp1	<- melt(sim.base07[,"Exp"] * sim.base07[,fmt_name])
tmp2	<- melt(sim.base08[,"Exp"] * sim.base08[,fmt_name])
ggtmp1	<- rbind(data.frame(tmp1, year = 2007), data.frame(tmp2, year = 2008))
names(ggtmp1)	<- c("inc_idx", "Retail", "Exp", "year")
ggtmp1$Inc	<- sim.unq[ggtmp1$inc_idx,"income2007"]
sel		<- ggtmp1$year == 2008
ggtmp1[sel,"Inc"]	<- sim.unq[ggtmp1[sel,"inc_idx"], "income2008"]
ggtmp1$Inc1	<- factor(ggtmp1$Inc, levels = sort(ggtmp1$Inc))

if(make_plot){
	pdf(paste(plot.wd, "/estrun_",run_id,"/graph_ctrfact_checkbase_seg", seg_id, ".pdf", sep=""), width = ww, height = ww)
	print(ggplot(ggtmp, aes(y, omega, col = factor(Income))) + geom_line() + 
			facet_wrap(~year) + 
			guides(col = FALSE) + 
			labs(title = "Inclusive value omega(Inc, y)")
		)
	print(ggplot(ggtmp1, aes(Inc1, Exp, fill = Retail)) + geom_bar(stat = "identity", position = "stack") + 
			facet_wrap(~year) + 
			labs(x = "Income", title = "Baseline simulation expenditure allocation") + 
			theme(axis.text.x = element_text(angle = 60))
		)
	dev.off()
}

#------------------------#
# Simulation of strategy #
# Assume only grocery stores change retail attributes. 
iv.change		<- .1			# Percentage increase of independent variable

sim.parallel	<- function(var){
	out			<- data.frame()
	for(i in 1:R){
		sel.retail	<- fmt_name[i]
		X_list_new	<- X_list07
		if(var %in% c("ln_upc_per_mod", "ln_num_module")){
			X_list_new[[sel.retail]][var]	<- X_list_new[[sel.retail]][var] + log(1 + iv.change)
		}else{
			X_list_new[[sel.retail]][var]	<- X_list_new[[sel.retail]][var] * (1 + iv.change)
		}
		out1	<- SimOmega_fn(ln_inc = sim.unq$Inc08, lambda = lambda, param_est = shr.par, 
					base = beta0_base, X_list = X_list_new, price = price.07, 
					lnInc_lv = lnInc_08, y.nodes = y.nodes, eps_draw = eps_draw, method = exp.method, interp.method = interp.method)
		out1 	<- data.frame(out1,lnInc = sim.unq$Inc08, retailer = fmt_name[i], Var = var)
		out		<- rbind(out, out1)
	}
	out$Var	<- as.character(out$Var)
	return(out)
}

pct			<- proc.time()
sim.X.ls	<- foreach(var = selcol, .combine = rbind) %dopar% {
	sim.parallel(var)
}
use.time	<- proc.time() - pct
cat("Parallel simulation of retail attributes finishes with", use.time[3]/60, "min.\n")

#-------------------------------#
# Simulation of price reduction # 
pct			<- proc.time()
sim.prc.ls	<- foreach(sel.retail = fmt_name, .combine = rbind) %dopar% {
	price.new	<- price.07
	price.new[sel.retail]	<- price.07[sel.retail] * (1 - iv.change)
	out		<- SimOmega_fn(ln_inc = sim.unq$Inc08, lambda = lambda, param_est = shr.par, 
				base = beta0_base, X_list = X_list07, price = price.new, 
				lnInc_lv = lnInc_08, y.nodes = y.nodes, eps_draw = eps_draw, method = exp.method, interp.method = interp.method)
	out		<- data.frame(out, lnInc = sim.unq$Inc08, retailer = as.character(sel.retail))
	return(out)
}
use.time	<- proc.time() - pct
cat("Parallel simulation of retail price finishes with", use.time[3]/60, "min.\n")

stopCluster(cl)
cat("Stop clustering. \n")

################
# Save results #
################
save.image(paste("estrun_",run_id,"/ctrfact_seg",seg_id, "_", Sys.Date(), ".rdata",sep=""))

cat("This program is done.\n")


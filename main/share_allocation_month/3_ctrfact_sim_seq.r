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

# seg_id	<- as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))
seg_id 	<- 1
args <- commandArgs(trailingOnly = TRUE)
print(args)
if(length(args)>0){
    for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
    }
}
cat("seg_id =", seg_id, ".\n")

model_name 	<- "MDCEV_share"
# setwd("~/Documents/Research/Store switching/processed data/Estimation")
# plot.wd	<- '~/Desktop'
# source("../../Exercise/Multiple_discrete_continuous_model/0_Allocation_function.R")
# source("../../Exercise/main/share_allocation_month/ctrfact_sim_functions.r")

# setwd("/home/brgordon/ccv103/Exercise/run")
# setwd("/kellogg/users/marketing/2661703/Exercise/run")
# setwd("/sscc/home/c/ccv103/Exercise/run")
setwd("U:/Users/ccv103/Documents/Research/Store switching/run")

run_id		<- 10
plot.wd		<- getwd()
make_plot	<- TRUE
ww			<- 10
ar			<- .6

source("0_Allocation_function.R")
source("ctrfact_sim_functions.r")

# Load estimation data 
ver.date	<- "2016-06-30"
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
trim.alpha		<- 0
numsim			<- 100	
draw.par		<- FALSE
sim.omega		<- FALSE
fname			<- paste("ctrfact_seq_seg",seg_id,sep="")
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
tmp		<- tmp[,list(income = unique(income_midvalue)), by = list(zip3, household_code, year)]
sim.data<- data.frame(tmp)[,c("household_code","income","zip3")]
sim.unq		<- data.frame(unique(sim.data[,-1]))
names(sim.unq)	<- gsub("income", "income2007", names(sim.unq))

# Counterfactual scenario: income is lower by 10%. 
my.change		<- 0 #	.1
lnInc_08		<- lnInc + log(1 - my.change)
sim.unq$income2008	<- (1 - my.change) * sim.unq$income2007		
sim.unq$Inc07	<- log(sim.unq[,"income2007"])
sim.unq$Inc08	<- log(sim.unq[,"income2008"])
sim.unq			<- sim.unq[order(sim.unq$income2007),]
cat("dim(sim.unq) =", dim(sim.unq), "\n")

#----------------------------#
# Accessibility nodes
tmp			<- dcast(subset(penetrat, channel_type %in% fmt_name), zip3 ~ channel_type, value.var = "penetration")
acs			<- as.matrix(tmp[,-1])
rownames(acs)	<- tmp$zip3
acs_nodes		<- acs[as.character(sim.unq$zip3),]

# Average price in year 2007
if(week.price){
	tmp 		<- dcast(subset(price_dat, year %in% selyr),	 
						scantrack_market_descr + year + biweek ~ channel_type, value.var = "bsk_price_paid_2004") 
	price.07	<- setNames( colMeans(as.matrix(tmp[,4:(3+R)]), na.rm=T), fmt_name)
}else{
	tmp 		<- dcast(subset(price_dat, year %in% selyr),	 
						scantrack_market_descr + year ~ channel_type, value.var = "bsk_price_paid_2004") 
	price.07	<- setNames( colMeans(as.matrix(tmp[,3:(2+R)]), na.rm=T), fmt_name)
}
cat("The average price level in 2007:\n"); print(price.07);cat("\n")

# Average retail attributes in 2007
selcol	<- c("size_index", "ln_upc_per_mod", "ln_num_module","overall_prvt")
X_list0<- setNames( lapply(fmt_name, function(x) 
									colMeans(as.matrix(subset(fmt_attr, channel_type == x & year%in% selyr)[,selcol]))), 
					 fmt_name)
cat("The average retail attributes in 2007:\n"); print(do.call(rbind, X_list0)); cat("\n")

# Expand X_list and price to match the nobs of income 
price.07	<- rep(1, nrow(sim.unq)) %*% matrix(price.07, nrow = 1)
colnames(price.07)	<- fmt_name
X_list08	<- X_list07	<- setNames(vector("list", R), fmt_name)
for(i in 1:R){
	tmp2	<- tmp1	<- matrix(0, nrow(sim.unq), R-1)
	if(i<beta0_base){
		tmp1[,i]	<- 0			# Recession = 0 for year = 2007
		tmp2[,i]	<- 1			# Recession = 1 for year = 2008
	}else if(i>beta0_base){
		tmp1[,(i-1)]	<- 0
		tmp2[,(i-1)]	<- 1
	}

	X_list07[[i]]	<- cbind(tmp1, acs_nodes[,i], rep(1, nrow(sim.unq)) %*% t(X_list0[[i]]), matrix(0, nrow(sim.unq), length(selcol)+1))
	X_list08[[i]]	<- cbind(tmp2, acs_nodes[,i], rep(1, nrow(sim.unq)) %*% t(X_list0[[i]]), acs_nodes[,i], rep(1, nrow(sim.unq)) %*% t(X_list0[[i]]))
	colnames(X_list07[[i]])	<- colnames(X_list08[[i]])	<- c(fmt_name[-beta0_base], "access", selcol, paste("rec*",c("access", selcol), sep=""))
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

# Compare the deterministic vs. random components
tmp		<- param_assignR(shr.par, nx, R, beta0_base)
tmp1	<- sapply(1:R, function(i) X_list[[i]] %*% tmp$beta + tmp$beta0[i])
tmp2	<- sapply(1:R, function(i) X_list07[[i]] %*% tmp$beta + tmp$beta0[i])
tmp3	<- sapply(1:R, function(i) X_list08[[i]] %*% tmp$beta + tmp$beta0[i])
cat("The expected xbeta from the actual data:\n"); print((colMeans(tmp1))); cat("\n")
cat("The expected xbeta of the simulation X in year 2007:\n"); print((colMeans(tmp2))); cat("\n")
cat("The expected xbeta of the simulation X in year 2008:\n"); print((colMeans(tmp3))); cat("\n")
cat("Summary stats of random draws:\n"); print(summary(c(eps_draw))); cat("\n")

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
mycore 	<- 4
# cl		<- makeCluster(mycore, type = "FORK")
cl		<- makeCluster(mycore, type = "PSOCK")			# Register cores in Windows
registerDoParallel(cl)
cat("Register", mycore, "core parallel computing. \n")

# Simulate expenditure and expenditure share. 
pct				<- proc.time()
sim.base07		<- foreach(i = 1:nrow(sim.unq), .combine = rbind, .packages=c("nloptr","mgcv")) %dopar%{
	tmpX	<- lapply(X_list07, function(x) matrix(x[i,], nrow = 1))
	out		<- SimWrapper_fn(omega_deriv, ln_inc = sim.unq[i,"Inc07"], lambda = lambda, param_est = shr.par, base = beta0_base, 
						X_list = tmpX, price = matrix(price.07[i,], nrow=1), eps_draw = eps_draw, method = exp.method, ret.sim = TRUE, par.draw = par_draw, 
						share.state = matrix(acs_nodes[i,], nrow = 1))
	out$Average
}
use.time		<- proc.time() - pct						
cat("2007 Baseline simulation finishes with", use.time[3]/60, "min.\n")

pct				<- proc.time()
sim.base08		<- foreach(i = 1:nrow(sim.unq), .combine = rbind, .packages=c("nloptr","mgcv")) %dopar%{
	tmpX	<- lapply(X_list08, function(x) matrix(x[i,], nrow = 1))
	out		<- SimWrapper_fn(omega_deriv1, ln_inc = sim.unq[i,"Inc08"], lambda = lambda, param_est = shr.par, base = beta0_base, 
						X_list = tmpX, price = matrix(price.07[i,], nrow=1), eps_draw = eps_draw, method = exp.method, ret.sim = TRUE, par.draw = par_draw, 
						share.state = matrix(acs_nodes[i,], nrow = 1))
	return(out$Average)
}
use.time		<- proc.time() - pct						
cat("2008 Baseline simulation finishes with", use.time[3]/60, "min.\n")

# Check the baseline simulation 
tmp1	<- sim.base07[,fmt_name]/sim.base07[,"Exp"]
tmp2	<- sim.base08[,fmt_name]/sim.base08[,"Exp"]
tmp		<- tmp2-tmp1
cat("Summary of baseline income elasticity:\n"); print(tmp); cat("\n")

# Plot omega over y
tmp		<- seq(50, 250, 10)
tmpd1	<- data.frame(lnInc = rep(unique(sim.unq$Inc07), each = length(tmp)), y = rep(tmp, length(unique(sim.unq$Inc07))))
tmpd1$Xacs	<- kronecker(rep(1, nrow(tmpd1)), matrix(acs_nodes[1,], nrow=1))
tmp1	<- omega_deriv(tmpd1)
ggtmp	<- cbind(tmpd1[,c("lnInc","y")], omega = tmp1, year = 2007)
ggtmp$Income	<- exp(ggtmp$lnInc)

# Plot the simulation baseline -- expenditure allocation by income 
tmp1	<- melt(data.frame(sim.unq[,c("income2007", "zip3")], sim.base07[,fmt_name]), id.vars = c("income2007", "zip3"))
tmp2	<- melt(data.frame(sim.unq[,c("income2007", "zip3")], sim.base08[,fmt_name]), id.vars = c("income2007", "zip3"))
ggtmp1	<- rbind(data.frame(tmp1, year = 2007), data.frame(tmp2, year = 2008))
names(ggtmp1)	<- c("income2007", "zip3", "Retail", "Exp", "year")
# ggtmp1$Inc	<- round(exp(ggtmp1$inc_idx), 0)
# ggtmp1$Inc1	<- factor(ggtmp1$Inc, levels = sort(unique(ggtmp1$Inc)))
ggtmp1		<- data.table(ggtmp1)
ggtmp1		<- ggtmp1[,list(Exp=sum(Exp)), by = list(income2007, Retail, year)]
ggtmp1		<- ggtmp1[,Exp:=Exp/sum(Exp), by = list(income2007, year)]
ggtmp1		<- ggtmp1[,list(Exp=mean(Exp)), by = list(income2007, Retail, year)]

# Check if utility function at the top level is concave
tmp		<- seq(300, 350, 1)
u		<- matrix(NA, length(tmp), length(lnInc), dimnames = list(y = tmp, lnInc = lnInc))
for(i in 1:length(lnInc)){
	f 	<- function(y, ...){
		newx	<- data.frame(lnInc = lnInc[i], y = y)
		newx$Xacs	<- rep(1, length(y)) %*% matrix(acs_nodes[1,], nrow=1)
		omega_deriv(newx, ...)
	}
	u[,i]<- -exp_fn(tmp, lambda, f, ln_inc = lnInc[i])
}
ggtmp2	<- melt(u, value.name = "utility")

if(make_plot){
	pdf(paste(plot.wd, "/estrun_",run_id,"/graph_checkbase_", fname, ".pdf", sep=""), width = ww/ar, height = ww)
	print(ggplot(ggtmp, aes(y, omega, col = factor(Income))) + geom_line() + 
			# facet_wrap(~year) + 
			guides(col = FALSE) + 
			labs(title = "Inclusive value omega(Inc, y)")
		)
	print(ggplot(ggtmp1, aes(income2007, Exp, fill = Retail)) + geom_bar(stat = "identity", position = "stack") + 
			facet_wrap(~year) + 
			labs(x = "Income", title = "Baseline simulation expenditure allocation") + 
			theme(axis.text.x = element_text(angle = 60))
		)
	print(ggplot(ggtmp2, aes(y, utility)) + geom_point() + geom_line() + facet_wrap(~lnInc, scales = "free"))
	dev.off()
}

save.image(paste("estrun_",run_id, "/", fname, ".rdata",sep=""))

#------------------------#
# Simulation of strategy #
# Assume only grocery stores change retail attributes. 
iv.change	<- -.1

#-------------------------------#
# Simulation of price reduction # 
pct			<- proc.time()
sim.prc.ls08	<- foreach(sel.retail = fmt_name, .errorhandling = "remove", .packages=c("mgcv", "nloptr")) %dopar% {
	out1.avg	<- data.frame()
	out1.alc	<- array(NA, c(numsim, nrow(sim.unq), R), 
						dimnames = list(iter = 1:numsim, lnInc = sim.unq[,"Inc08"], retail = fmt_name))
	
	price.new	<- price.07
	price.new[,sel.retail]	<- price.07[,sel.retail] * (1 + iv.change)
	if(sim.omega){
		out1	<- SimOmega_fn(ln_inc = sim.unq[,"Inc08"], lambda = lambda, param_est = shr.par, 
					base = beta0_base, X_list = X_list08, price = price.new, 
					lnInc_lv = lnInc_08, y.nodes = y.nodes, eps_draw = eps_draw, method = exp.method, 
					interp.method = interp.method, ret.sim = TRUE, alpha = trim.alpha, par.draw = par_draw)
	}else{
		out1	<- SimWrapper_fn(omega_deriv = omega_deriv1, ln_inc = sim.unq[,"Inc08"], lambda = lambda, param_est = shr.par, 
					base = beta0_base, X_list = X_list08, price = price.new, sim.y = rep(1,numsim) %*% t(sim.base08[,"Exp"]), 
					eps_draw = eps_draw, method = exp.method, ret.sim = TRUE, par.draw = par_draw, share.state = acs_nodes)
	}
	out1.avg	<- rbind(out1.avg, data.frame(out1$Average, 
										lnInc = sim.unq[,"Inc08"], retailer = as.character(sel.retail), 
										Var = 'price', change = iv.change ))
	out1.alc	<- out1$Allocation

	out1	<- list(Average = out1.avg, Allocation = out1.alc)
	return(out1)
}
use.time	<- proc.time() - pct
cat("Parallel simulation of retail price finishes with", use.time[3]/60, "min.\n")

stopCluster(cl)
cat("Stop clustering. \n")

################
# Save results #
################
rm(list  = intersect(ls(), c("ar", "args", "cl", "ggtmp", "ggtmp1", "ggtmp2", "i", "lastFuncGrad", "lastFuncParam", "make_plot", "mycore", 
			"myfix", "plot.wd", "s1_index", "sel", "selyr", "price", "sol", "sol.top", "sol.top2", "tmp", "tmp_coef", 
			"tmp_price", "tmp1", "tmp2", "use.time", "ver.date", "var1", "ww", "f", 
			"numnodes", "out", "out1", "pct", "tmpd1", "tmpd2", "tmpdat", "u", "W", "y", "y.nodes",
			"Allocation_constr_fn", "Allocation_fn", "Allocation_nlop_fn", "cheb.1d.basis", "cheb.basis", "chebfun", 
			"exp_fn", "expFOC_fn", "incl_value_fn", "mysplfun", "mytrimfun", "param_assignR", "simExp_fn", "SimOmega_fn", 
			"SimWrapper_fn", "solveExp_fn", "spl2dfun", "uP_fn", "uGrad_fn")))

save.image(paste("estrun_",run_id, "/", fname, ".rdata",sep=""))

cat("This program is done.\n")


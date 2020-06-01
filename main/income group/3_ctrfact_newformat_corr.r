library(ggplot2)
library(reshape2)
library(maxLik)
library(evd)
library(data.table)
library(doParallel)
library(foreach)
# library(chebpol)
library(nloptr)
library(mgcv)
options(error = quote({dump.frames(to.file = TRUE)}))

# seg_id	<- as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))

args <- commandArgs(trailingOnly = TRUE)
print(args)
if(length(args)>0){
    for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
    }
}
# seg_id <- 1
cat("seg_id =", seg_id, ".\n")

model_name 	<- "MDCEV_share"
# setwd("~/Documents/Research/Store switching/Processed_data/processed data_20160809")
# plot.wd	<- '~/Desktop'
# sourceCpp(paste("../../Exercise/Multiple_discrete_continuous_model/", model_name, ".cpp", sep=""))
# source("../../Exercise/Multiple_discrete_continuous_model/0_Efficient_Allocation_function.R")
# source("../../Exercise/main/ctrfact_sim_functions.r")

# setwd("/home/brgordon/ccv103/Exercise/run")
# setwd("/kellogg/users/marketing/2661703/Exercise/run")
# setwd("/sscc/home/c/ccv103/Exercise/run")
setwd("U:/Users/ccv103/Documents/Research/Store switching/run")

run_id		<- 2
plot.wd		<- getwd()
make_plot	<- TRUE
ww			<- 10
ar			<- .6

source("0_Efficient_Allocation_function.R")
source("ctrfact_sim_functions.r")

# Load estimation data 
ver.date	<- "2016-08-28"		#"2016-08-22"
cpi.adj		<- TRUE

loadf	<- paste("estrun_",run_id,"/MDCEV_est_seg",seg_id,"_", ver.date,".rdata",sep="")
loadf
load(loadf)
rm(list = intersect(ls(), c("gamfit", "shr","model_name", "tmpdat")))

# Set simulation parameters
week.price		<- FALSE
interp.method	<- "spline"					# Spline interpolation or Chebyshev interpolation
exp.method		<- "Utility"				# Utility maximization or finding roots for first order condition
trim.alpha		<- 0
numsim			<- 500
draw.par		<- FALSE
sim.omega		<- FALSE
fname			<- paste("newfmt_cor_seg",seg_id,sep="")
fname			<- paste(fname, "_sim", numsim, "_", as.character(Sys.Date()), sep="")
cat("Output file name is", fname, ".\n")

###############################
# Prepare simulation elements # 
###############################
# Data required: parameters, income level, price, retail attributes, and random draws
# For each simulation scenario, if income, price or retail attributes change, we need to re-simulate inclusive values.

# Set simulation parameters
lambda 	<- coef(lsol)	#coef(sol.top2)
shr.par	<- coef(sol)
ret.idx	<- 11:15		# Column index flagging the retail attribugtes in X1
nx0		<- 2			# Number of demograhic var in X1

# Add parameter names 
selcol1			<- c("distance_haver", selcol)
names(shr.par)	<- c(paste(rep(c("Intercept", "Rec"), each = R-1), ":", rep(fmt_name[-beta0_base], nx0), sep=""), 
				    selcol1, paste("Rec:", selcol1, sep=""), 
					paste("gamma0_", rep(c("Intercept", "Rec"), each=R), ":", rep(fmt_name, nx0), sep=""), 
					paste("gamma0_", selcol1, sep=""), paste("gamma0_Rec:", selcol1, sep=""), "ln_sigma"
					)

#-----------------------#
# Construct income data # 
selyr	<- 2007
sel		<- mydata$year == selyr & !duplicated(mydata[,c("household_code", "year")]) # & mydata$household_code < 2250000
sim.unq	<- mydata[sel,c("household_code", "ln_inc")]
sim.unq$price	<- as.matrix(mydata[sel,paste("PRC_", gsub("\\s", "_", fmt_name), sep="")])
tmp		<- do.call(cbind, lapply(X1, function(x) x[sel,ret.idx]))
sim.unq$X<- as.matrix( tmp)
cat("dim(sim.unq) =", dim(sim.unq), "\n")

# attr_mat	<- merge(mydata[,c("household_code", "year", "month", "scantrack_market_descr")], 
# 					hh_dist[,c("household_code", "year", "channel_type", "distance_haver")], 
# 				by = c("household_code", "year"), all.x=T)
tmp		<- data.table(hh_dist[,c("household_code", "panelist_zip_code","year", "channel_type", "distance_haver")])
tmp		<- tmp[, maxd:=max(distance_haver, na.rm=T), by = list(panelist_zip_code)]
tmp		<- tmp[, distance_haver := distance_haver/maxd]
tmp		<- data.frame(tmp)
attr_mat	<- merge(mydata[,c("household_code", "year", "month", "scantrack_market_descr")], 
					tmp[,c("household_code", "year", "channel_type", "distance_haver")], 
					by = c("household_code", "year"), all.x=T)
attr_mat	<- merge(attr_mat, fmt_attr[,c("scantrack_market_descr", "year", "channel_type", selcol)], 
				by = c("scantrack_market_descr", "year", "channel_type"))				
attr_mat	<- attr_mat[order(attr_mat$household_code, attr_mat$year, attr_mat$month),]	
cat("dim(mydata) =", dim(mydata), ".\n")
cat("dim(attr_mat) =", dim(attr_mat), ", dim(mydata)*R =",nrow(mydata)*R, ".\n")
selcol1		<- c("distance_haver", selcol)

# Retail attributes
X1.07	<- lapply(X1, function(x) x[sel,])
X2.07	<- lapply(X2, function(x) x[sel,])
X2.08	<- X1.08 	<- vector("list", length=length(fmt_name))
for(i in 1:length(fmt_name)){
	tmp0<- rep(0, R-1)
	if(i<beta0_base){
		tmp0[i]	<- 1
	}else if(i > beta0_base){
		tmp0[(i-1)]	<- 1
	}
	tmp1<- as.matrix(cbind(1, rec=rep(1, sum(sel))))	
	tmp2<- kronecker(tmp1, t(tmp0))
	colnames(tmp2)	<- paste(rep(c("Intercept", "Rec"), each = R-1), fmt_name[-beta0_base], sep=":")
	tmp3<- rep(0, R)
	tmp3[i]	<- 1
	tmp3<- kronecker(tmp1, t(tmp3))
	colnames(tmp3)	<- paste(rep(c("Intercept", "Rec"), each = R), fmt_name, sep=":")
	sel1	<- attr_mat$channel_type == fmt_name[i]
	tmp	<- as.matrix(attr_mat[sel1,selcol1][sel,])
	X1.08[[i]]	<- cbind(tmp2, tmp, tmp*1)
	X2.08[[i]]	<- cbind(tmp3, tmp, tmp*1)
	colnames(X1.08[[i]])	<- c(colnames(tmp2), colnames(tmp), paste(colnames(tmp), ":Rec", sep=""))
	colnames(X2.08[[i]])	<- c(colnames(tmp3), colnames(tmp), paste(colnames(tmp), ":Rec", sep=""))
	
	colnames(X1.07[[i]])	<- c(colnames(tmp2), colnames(tmp), paste(colnames(tmp), ":Rec", sep=""))
	colnames(X2.07[[i]])	<- c(colnames(tmp3), colnames(tmp), paste(colnames(tmp), ":Rec", sep=""))
}

#-------------------#
# Take random draws #
set.seed(666)
eps_draw	<- matrix(rgev(numsim*R, scale = exp(shr.par["ln_sigma"])), numsim, R)
par_draw <- NULL

##############
# Simulation #
##############
# Register parallel computing
mycore 	<- 3
# cl		<- makeCluster(mycore, type = "FORK")
cl		<- makeCluster(mycore, type = "PSOCK")			# Register cores in Windows
registerDoParallel(cl)
cat("Register", mycore, "core parallel computing. \n")

# ------------------ #
# Basline simulation #
# Simulate expenditure and expenditure share. 
# pct				<- proc.time()
# sim.base07		<- SimWrapper_fn(omega_deriv, lambda, shr.par, base = beta0_base, 
# 								X1 = X1.07, X2=X2.07, price = sim.unq$price, 
# 								demo_mat = as.matrix(cbind(Intercept = 1, rec = 0, ln_inc = sim.unq$ln_inc)), 
# 								eps_draw = eps_draw, ret.idx = ret.idx, sim.y = NULL, ret.sim = TRUE, method = exp.method)
# use.time		<- proc.time() - pct						
# cat("2007 Baseline simulation finishes with", use.time[3]/60, "min.\n")

pct				<- proc.time()
sim.base08		<- SimWrapper_fn(omega_deriv, lambda, shr.par, base = beta0_base, 
								X1 = X1.08, X2=X2.08, price = sim.unq$price, 
								demo_mat = as.matrix(cbind(Intercept = 1, rec = 1, ln_inc = sim.unq$ln_inc)), 
								eps_draw = eps_draw, ret.idx = ret.idx, sim.y = NULL, ret.sim = TRUE, method = exp.method)
use.time		<- proc.time() - pct						
cat("2008 Baseline simulation finishes with", use.time[3]/60, "min.\n")

# -------------- #
# Counterfactual #
# Introduct a new small format similar to convenience stores
# Scenarios:
# a(2): if intercept inherents the grocery stores, discount stores;
# b(2): if price inherents the intercept
 
ret.a	<- c("Discount Store")
sim.scn	<- c("Convenience Store", "Drug Store")
newf.sim07	<- newf.sim <- setNames(vector("list", length(ret.a)*length(sim.scn)), paste(ret.a, sim.scn, sep="_"))

cat(" ---------------------------------------------------------------------- \n")
pct		<- proc.time()

# Set new parameters
sel		<- which(fmt_name == ret.a)
tmp		<- shr.par[1:(nx0*(R-1))]
beta1	<- rbind(matrix(tmp, R-1, nx0), tmp[grep(ret.a, names(tmp))])
tmp		<- shr.par[(nx+1):(nx0*R+nx)]
beta2	<- rbind(matrix(tmp, R, nx0), tmp[grep(ret.a, names(tmp))])
shr.par.new	<- c(c(beta1), shr.par[(nx0*(R-1)+1):nx], c(beta2), shr.par[-c(1:(nx0*R+nx))])
cat("The parameter vector for adding a new format from", ret.a, "is:\n"); print(shr.par.new); cat("\n")

sel		<- mydata$year == selyr & !duplicated(mydata[,c("household_code", "year")]) # & mydata$household_code < 2250000
fmt_name1<- c(fmt_name, "New")
rho.vec	<- seq(-1, 1, .1)
xi.draw	<- rgev(numsim, scale = 1) 

for(j in 1:length(sim.scn)){
	# Set new attributes
	X2.08.new	<- X2.07.new	<- X1.08.new	<- X1.07.new	<- setNames(vector("list", R+1), fmt_name1)

	for(i in 1:R){
		# Covariate matrix for baseline function
		tmp0<- rep(0, R)
		if(i<beta0_base){
			tmp0[i]	<- 1
		}else if(i > beta0_base){
			tmp0[(i-1)]	<- 1
		}
		tmp1<- as.matrix(cbind(1, rec=rep(1, nrow(sim.unq))))	
		tmp2<- kronecker(tmp1, t(tmp0))
		colnames(tmp2)	<- paste(rep(c("Intercept", "Rec"), each = R), fmt_name1[-beta0_base], sep=":")
		tmp3<- rep(0, R+1)
		tmp3[i]	<- 1
		tmp3<- kronecker(tmp1, t(tmp3))
		colnames(tmp3)	<- paste(rep(c("Intercept", "Rec"), each = R+1), fmt_name1, sep=":")
		sel1	<- attr_mat$channel_type == fmt_name[i]
		tmp	<- as.matrix(attr_mat[sel1,selcol1][sel,])
		X1.08.new[[i]]	<- cbind(tmp2, tmp, tmp*1)
		X2.08.new[[i]]	<- cbind(tmp3, tmp, tmp*1)
		colnames(X1.08.new[[i]])	<- c(colnames(tmp2), colnames(tmp), paste(colnames(tmp), ":Rec", sep=""))
		colnames(X2.08.new[[i]])	<- c(colnames(tmp3), colnames(tmp), paste(colnames(tmp), ":Rec", sep=""))
	}
	
	# Set attributes for the new format
	tmp0	<- rep(0, R)
	tmp0[R]	<- 1
	tmp1	<- as.matrix(cbind(1, rec=rep(1, nrow(sim.unq))))	
	tmp2	<- kronecker(tmp1, t(tmp0))
	colnames(tmp2)	<- paste(rep(c("Intercept", "Rec"), each = R), fmt_name1[-beta0_base], sep=":")
	tmp3	<- rep(0, R+1)
	tmp3[(R+1)]	<- 1
	tmp3<- kronecker(tmp1, t(tmp3))
	colnames(tmp3)	<- paste(rep(c("Intercept", "Rec"), each = R+1), fmt_name1, sep=":")
	
	sel1	<- attr_mat$channel_type == sim.scn[j]
	tmp	<- as.matrix(attr_mat[sel1,selcol1][sel,])
	sel1	<- attr_mat$channel_type == ret.a
	tmp[,"overall_prvt"]	<- attr_mat[sel1,"overall_prvt"][sel]
	X1.08.new[[(R+1)]]	<- cbind(tmp2, tmp, tmp*1)
	X2.08.new[[(R+1)]]	<- cbind(tmp3, tmp, tmp*1)
	colnames(X1.08.new[[(R+1)]])	<- c(colnames(tmp2), colnames(tmp), paste(colnames(tmp), ":Rec", sep=""))
	colnames(X2.08.new[[(R+1)]])	<- c(colnames(tmp3), colnames(tmp), paste(colnames(tmp), ":Rec", sep=""))
	
	# Set prices 
	price.new	<- sim.unq$price
	price.new	<- cbind(price.new, New = sim.unq$price[,which(fmt_name==ret.a)])
	
	newf.sim[[j]]	<- foreach(rho = rho.vec, .errorhandling = "remove", .packages=c("mgcv", "nloptr")) %dopar% {
		eps_draw_new<- cbind(eps_draw, rho*eps_draw[,which(fmt_name==ret.a)] + sqrt(1-rho^2)*exp(shr.par["ln_sigma"])*xi.draw)
		out	<- SimWrapper_fn(omega_deriv, lambda, shr.par.new, base = beta0_base, 
							X1 = X1.08.new, X2=X2.08.new, price = price.new, 
							demo_mat = as.matrix(cbind(Intercept = 1, rec = 1, ln_inc = sim.unq$ln_inc)),  
							eps_draw = eps_draw_new, ret.idx = ret.idx+nx0, sim.y = sim.base08$y, ret.sim = FALSE, method = exp.method)
		out	<- cbind(rho = rho, out)	
		return(out)			
	}
}
use.time	<- proc.time() - pct
cat("Counterfactual with", ret.a, "finishes with", use.time[3]/60, "min.\n")
stopCluster(cl)
cat("Stop clustering. \n")

# Compare the difference 
tmp1	<- colSums(cbind(sim.base08$Average[,fmt_name],0), na.rm=T)
sel		<- which(rho.vec == 0)
tmp2	<- colSums(newf.sim[[1]][[sel]][,fmt_name1], na.rm=T)
mkt.dif	<- setNames(tmp2/sum(tmp2)-tmp1/sum(tmp1), fmt_name1)
cat("Difference of market share:\n"); print(round(mkt.dif*100,2)); cat("\n")

tmp2	<- colSums(newf.sim[[2]][[sel]][,fmt_name1], na.rm=T)
mkt.dif	<- setNames(tmp2/sum(tmp2)-tmp1/sum(tmp1), fmt_name1)
cat("Difference of market share in scenario", sim.scn[2], "\n"); print(round(mkt.dif*100,2)); cat("\n")

# Calculat the sum of market share of discount stores and new format
ggtmp	<- data.frame()
for(i in 1:length(newf.sim[[1]])){
	tmp	<- colSums(newf.sim[[1]][[i]][,fmt_name1], na.rm=T)
	tmp	<- tmp/sum(tmp)
	ggtmp	<- rbind(ggtmp, data.frame(rho = newf.sim[[1]][[i]][1,"rho"], share = tmp[ret.a] + tmp["New"]))
}
cat("Net market share over a pricing range is:\n"); print(ggtmp); cat("\n")

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

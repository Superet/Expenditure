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

 
# args <- commandArgs(trailingOnly = TRUE)
# print(args)
# if(length(args)>0){
#     for(i in 1:length(args)){
#       eval(parse(text=args[[i]]))
#     }
# }
# cat("seg_id =", seg_id, ".\n")

# setwd("~/Documents/Research/Store switching/Processed_data/processed data_20160809")
# plot.wd	<- '~/Desktop'
# sourceCpp(paste("../../Exercise/Multiple_discrete_continuous_model/", model_name, ".cpp", sep=""))
# source("../../Exercise/Multiple_discrete_continuous_model/0_Efficient_Allocation_function.R")
# source("../../Exercise/main/ctrfact_sim_functions.r")

# setwd("/home/brgordon/ccv103/Exercise/run")
# setwd("/kellogg/users/marketing/2661703/Exercise/run")
setwd("/sscc/home/c/ccv103/Exercise/run")
# setwd("U:/Users/ccv103/Documents/Research/Store switching/run")

source("0_Efficient_Allocation_function.R")
source("ctrfact_sim_functions.r")

# Load estimation data 
run_id		<- 6
ver.date	<- "2016-10-08"
loadf	<- paste("estrun_",run_id,"/retcat_explog_", ver.date,".rdata",sep="")
loadf
load(loadf)

# Set simulation parameters
exp.method		<- "Utility"				# Utility maximization or finding roots for first order condition
trim.alpha		<- 0
numsim			<- 500
draw.par		<- FALSE
sim.omega		<- FALSE
fname			<- "retcat_marg_effct"
if(draw.par){ fname <- paste(fname, "_pardraw", sep="") }
if(sim.omega){fname	<- paste(fname, "_simomega", sep="") }
fname			<- paste(fname, "_sim", numsim, "_", as.character(Sys.Date()), sep="")
cat("Output file name is", fname, ".\n")

# Simulation function 
SimWrapper_fn 	<- function(par, Inc, X, gamma, eps_draw, price = NULL, ret.sim = FALSE){
	R		<- length(X)
	N		<- length(Inc)
	numsim	<- nrow(eps_draw)
	if(is.null(price)){
		price	<- matrix(1, N, R)
	}
	e		<- array(NA, c(numsim, N, R), dimnames = list(iter = 1:numsim, obs = 1:N, retailer = 1:R))
	omega	<- matrix(NA, numsim, N)
	
	for(i in 1:numsim){
		psi		<- exp(sapply(X, function(x) x%*%par) + rep(1, N) %*% t(eps_draw[i,]))/gamma
		e[i,,]	<- Allocation_fn(Inc, psi, gamma, price, R, returnmax = FALSE)
	}
	out	<- apply(e, c(2,3), mean, na.rm= T)
	colnames(out)	<- names(X)
	if(ret.sim){
		out	<- list(Average = out, Allocation = e)
	}
	return(out)
}

###############################
# Prepare simulation elements # 
###############################
# Data required: parameters, income level, price, retail attributes, and random draws
# For each simulation scenario, if income, price or retail attributes change, we need to re-simulate inclusive values.

# Set simulation parameters
# sel		<- sapply(exp.est.ls, function(x) !is.null(x))
# lambda 	<- sapply(exp.est.ls[sel], coef)
lambda		<- do.call(cbind, exp.est.ls)
ret.idx	<- 11:15		# Column index flagging the retail attribugtes in X1
dpt.name	<- c("DG", "GM","NFG","RF", "HBC", "OTHER")
dpt.lab		<- c("DRY GROCERY","GENERAL MERCHANDISE", "NON-FOOD GROCERY", "REFRIGERATED FROZEN", "HEALTH BEAUTY CARE","PRODUCE OTHER")
nc			<- length(dpt.name)

#-----------------------#
# Construct income data # 
selyr	<- 2007
sel		<- expdata.full$year == selyr & !duplicated(expdata.full[,c("household_code", "year")]) # & expdata.full$household_code < 2250000
sim.unq	<- expdata.full[sel,c("household_code", "income_mid","rec", "segment", "scantrack_market_descr")]
sim.unq$income	<- sim.unq$income_mid/12
sim.unq	<- sim.unq[order(sim.unq$household_code),]
# sim.unq	<- sim.unq[sample(1:nrow(sim.unq), 100),]
cat("dim(sim.unq) =", dim(sim.unq), "\n")

# Organize price and retail attributes for each channel and department 
selcol	<- c("size_index", "ln_upc_per_mod", "ln_num_module", "prvt_overall")

# Normalize distance by zipcode
sel		<- hh_dist$year == selyr
dist	<- data.table(hh_dist[sel,c("household_code", "panelist_zip_code","year", "channel_type", "distance_haver")])
dist	<- dist[, maxd:=max(distance_haver, na.rm=T), by = list(panelist_zip_code)]
dist	<- dist[, distance := distance_haver/maxd]
dist	<- data.frame(dist)
dist	<- merge(sim.unq[,c("household_code","scantrack_market_descr")], dist, by = c("household_code"), all.x=T)
dist	<- dist[order(dist$household_code, dist$channel_type),]
cat("dim(dist) =", dim(dist), ", dim(sim.unq)*R =",nrow(sim.unq)*R, ".\n")
beta0.base	<- which(dpt.name == "OTHER")

sel			<- fmt_dpt$year == selyr
attr_mat	<- merge(sim.unq[,c("household_code","scantrack_market_descr")], 
				fmt_dpt[sel,c("scantrack_market_descr", "channel_type","department","unitprice_paid",selcol)], 
				by = c("scantrack_market_descr"), all.x=T)
attr_mat$department <- factor(attr_mat$department, levels = dpt.lab, labels = dpt.name)				
attr_mat	<- attr_mat[order(attr_mat$household_code, attr_mat$channel_type, attr_mat$department),]	
cat("dim(attr_mat) =", dim(attr_mat), ", dim(sim.unq)*nc*R =",nrow(sim.unq)*nc*R, ".\n")

# Retail attributes for simulation
ord		<- order(sim.unq$household_code)
X1.07	<- X1.08 	<- vector("list", length=length(fmt_name)+1)
names(X1.07)	<- names(X1.08)	<- c(fmt_name, "Outside")
gamma.07	<- gamma.08	<- matrix(1, nrow(sim.unq), R+1)
for(i in 1:length(fmt_name)){
	# Price matrix 
	sel		<- attr_mat$channel_type == fmt_name[i]
	tmp		<- dcast(attr_mat[sel,c("household_code", "department", "unitprice_paid")], 
					household_code ~ department, value.var = "unitprice_paid")
	names(tmp)	<- c("household_code", paste("PRC_", dpt.name, sep=""))
	price	<- as.matrix(tmp[ord,paste("PRC_",dpt.name, sep="")])
	
	# Retail attributes
	sel		<- attr_mat$channel_type == fmt_name[i] 
	tmp		<- do.call(cbind,lapply(dpt.name[-beta0.base], function(x) attr_mat[sel&attr_mat$department==x,selcol]))
	
	# Fit first and second derivative 
	for(j in 1:nseg){
		sel		<- sim.unq$segment == j 
		tmpdat	<- data.frame(y = rep(0, sum(sel)))
		tmpdat$price	<- price[sel,-beta0.base]
		tmpdat$X		<- as.matrix(tmp[sel,])
		tmpdat$rec		<- 0
		gamma.07[sel,i]	<- exp(predict(omega.lsf[[j]][[i]], tmpdat))
		tmpdat$rec		<- 1
		gamma.08[sel,i]	<- exp(predict(omega.lsf[[j]][[i]], tmpdat))
	}
	
	# Model matrix 
	tmp0	<- rep(0, R)
	tmp0[i]	<- 1
	tmp1	<- as.matrix(cbind(1, rec=sim.unq$rec))	
	tmp2	<- kronecker(tmp1, t(tmp0))
	colnames(tmp2)	<- paste(rep(c("Intercept", "Rec"), each = R), fmt_name, sep=":")
	sel		<- dist$channel_type == fmt_name[i]
	X1.07[[i]]	<- as.matrix(cbind(tmp2, dist[sel,"distance"], dist[sel,"distance"]*0))
	X1.08[[i]]	<- as.matrix(cbind(tmp2, dist[sel,"distance"], dist[sel,"distance"]*1))
	colnames(X1.07[[i]])	<- colnames(X1.08[[i]]) <- c(colnames(tmp2), "distance", "distance:Rec")
}
X1.07[[(R+1)]]	<- matrix(0, nrow(sim.unq), ncol(X1.07[[1]]))
X1.08[[(R+1)]]	<- matrix(0, nrow(sim.unq), ncol(X1.08[[1]]))
rm(list = c("price", "tmp0", "tmp1", "tmp2", "tmpdat", "tmpls", "tmp"))

#-------------------#
# Take random draws #
set.seed(666)
eps_draw	<- matrix(rgev(numsim*(R+1), scale = 1), numsim, R+1)

##############
# Simulation #
##############
# Register parallel computing
mycore 	<- 3
cl		<- makeCluster(mycore, type = "FORK")
# cl		<- makeCluster(mycore, type = "PSOCK")			# Register cores in Windows
registerDoParallel(cl)
cat("Register", mycore, "core parallel computing. \n")

# ------------------ #
# Basline simulation #
# Simulate expenditure and expenditure share. 
pct				<- proc.time()
sim.base07		<- foreach(i = 1:nseg, .errorhandling = "remove", .packages=c("mgcv", "nloptr")) %dopar% {
	sel		<- sim.unq$segment == i
	out		<- SimWrapper_fn(lambda[,i], Inc = sim.unq[sel,"income"], X = lapply(X1.07, function(x) x[sel,]), 
							gamma = gamma.07[sel,], eps_draw = eps_draw, price = NULL, ret.sim = TRUE)
	out$Average	<- cbind(segment = i, household_code = sim.unq[sel,"household_code"], out$Average)
	return(out)
}
use.time		<- proc.time() - pct						
cat("2007 Baseline simulation finishes with", use.time[3]/60, "min.\n")

pct				<- proc.time()
sim.base08		<- foreach(i = 1:nseg, .errorhandling = "remove", .packages=c("mgcv", "nloptr")) %dopar% {
	sel		<- sim.unq$segment == i
	out		<- SimWrapper_fn(lambda[,i], Inc = sim.unq[sel,"income"], X = lapply(X1.08, function(x) x[sel,]), 
							gamma = gamma.08[sel,], eps_draw = eps_draw, price = NULL, ret.sim = TRUE)
	out$Average	<- cbind(segment = i, household_code = sim.unq[sel,"household_code"], out$Average)
	return(out)
}
use.time		<- proc.time() - pct
cat("2008 Baseline simulation finishes with", use.time[3]/60, "min.\n")

iv.change <- .1
(sel.retail	<- which(fmt_name == "Discount Store"))

#-------------------------------- #
# Simulation of fixed expenditure #
# Price matrix 
sel		<- attr_mat$channel_type == fmt_name[sel.retail]
tmp		<- dcast(attr_mat[sel,c("household_code", "department", "unitprice_paid")], 
				household_code ~ department, value.var = "unitprice_paid")
names(tmp)	<- c("household_code", paste("PRC_", dpt.name, sep=""))
price	<- as.matrix(tmp[ord,paste("PRC_",dpt.name[-beta0.base], sep="")])

# Retail attributes
sel		<- attr_mat$channel_type == fmt_name[sel.retail] 
attr_sel	<- do.call(cbind,lapply(dpt.name[-beta0.base], function(x) attr_mat[sel&attr_mat$department==x,selcol]))

pct			<- proc.time()
sim.ls07	<- foreach(v = c(selcol, "price"), .errorhandling = "remove", .packages=c("mgcv", "nloptr")) %dopar% {
	out1	<- data.frame()
	for(j in 1:(nc-1)){
		price.new	<- price
		attr.new	<- attr_sel
		if(v == "price"){
			price.new[,j]	<- price.new[,j]*(1 + iv.change)
		}else if(v %in% c("ln_upc_per_mod", "ln_num_module")){
			sel	<- which(selcol == v) + (j-1)*length(selcol)
			attr.new[,sel]	<- attr.new[,sel] + log(1 + iv.change)
		}else{
			sel	<- which(selcol == v) + (j-1)*length(selcol)
			attr.new[,sel]	<- attr.new[,sel] * (1 + iv.change)
		}
		
		gamma.new	<- gamma.07
		for(i in 1:nseg){
			sel		<- sim.unq$segment == i
			tmp		<- data.frame(y = rep(0, sum(sel)))
			tmp$price	<- price.new[sel,]
			tmp$X		<- as.matrix(attr.new[sel,])
			tmp$rec		<- 0
			gamma.new[sel,sel.retail]	<- exp(predict(omega.lsf[[i]][[sel.retail]], newdata = tmp))
			tmp		<- SimWrapper_fn(lambda[,i], Inc = sim.unq[sel,"income"], X = lapply(X1.07, function(x) x[sel,]), 
									gamma = gamma.new[sel,], eps_draw = eps_draw, price = NULL, ret.sim = FALSE)
			tmp		<- data.frame(tmp, retailer = fmt_name[sel.retail], Var = v, department = dpt.name[j], change = iv.change, segment = i, household_code = sim.unq[sel,"household_code"])	
			out1	<- rbind(out1, tmp)
		}
	}
	return(out1)
}
use.time	<- proc.time() - pct
cat("Parallel simulation of retail price finishes with", use.time[3]/60, "min.\n")
names(sim.ls07)	<- sapply(sim.ls07, function(x) as.character(unique(x$Var)))

pct			<- proc.time()
sim.ls08	<- foreach(v = c(selcol, "price"), .errorhandling = "remove", .packages=c("mgcv", "nloptr")) %dopar% {
	out1	<- data.frame()
	for(j in 1:(nc-1)){
		price.new	<- price[,-beta0.base]
		attr.new	<- attr_sel
		if(v == "price"){
			price.new[,j]	<- price.new[,j]*(1 + iv.change)
		}else if(v %in% c("ln_upc_per_mod", "ln_num_module")){
			sel	<- which(selcol == v) + (j-1)*length(selcol)
			attr.new[,sel]	<- attr.new[,sel] + log(1 + iv.change)
		}else{
			sel	<- which(selcol == v) + (j-1)*length(selcol)
			attr.new[,sel]	<- attr.new[,sel] * (1 + iv.change)
		}
		
		gamma.new	<- gamma.08
		for(i in 1:nseg){
			sel		<- sim.unq$segment == i
			tmp		<- data.frame(y = rep(0, sum(sel)))
			tmp$price	<- price.new[sel,]
			tmp$rec		<- 1
			tmp$X		<- as.matrix(attr.new[sel,])
			gamma.new[sel,sel.retail]	<- exp(predict(omega.lsf[[i]][[sel.retail]], newdata = tmp))
			tmp		<- SimWrapper_fn(lambda[,i], Inc = sim.unq[sel,"income"], X = lapply(X1.08, function(x) x[sel,]), 
									gamma = gamma.new[sel,], eps_draw = eps_draw, price = NULL, ret.sim = FALSE)
			tmp		<- data.frame(tmp, retailer = fmt_name[sel.retail], Var = v, department = dpt.name[j], change = iv.change, segment = i, household_code = sim.unq[sel,"household_code"])	
			out1	<- rbind(out1, tmp)
		}
	}
	return(out1)
}
use.time	<- proc.time() - pct
cat("Parallel simulation of retail price finishes with", use.time[3]/60, "min.\n")
names(sim.ls08)	<- sapply(sim.ls08, function(x) as.character(unique(x$Var)))

stopCluster(cl)
cat("Stop clustering. \n")
rm(list = c("X1.07", "X1.08"))

##############################
# Summarize marginal effects # 
##############################
sel.retail	<- "Discount.Store" 
selcol1	<- c(selcol, "price")
base.07	<- data.frame(do.call(rbind, lapply(sim.base07, function(x) x$Average)))
base.08	<- data.frame(do.call(rbind, lapply(sim.base08, function(x) x$Average)))
max(abs(base.07[,"household_code"] - as.numeric(sim.ls07[[1]]$household_code)))

# Mean elasticity across population (by variable by category)
tmp1	<- tmp2	<- setNames(vector("list", length(sim.ls07)), names(sim.ls07))
tmp		<- NULL
for(i in 1:length(sim.ls07)){
	tmp1[[i]]	<- tmp2[[i]]	<- setNames(vector("list", nc-1), dpt.name[-beta0.base])
	for(j in 1:(nc-1)){
		sel	<- sim.ls07[[i]]$department == dpt.name[j]
		tmp1[[i]][[j]]	<- (sim.ls07[[i]][sel,sel.retail] - base.07[,sel.retail])/base.07[,sel.retail]/sim.ls07[[i]][sel,"change"]
		tmp2[[i]][[j]]	<- (sim.ls08[[i]][sel,sel.retail] - base.08[,sel.retail])/base.08[,sel.retail]/sim.ls08[[i]][sel,"change"]
		sel	<- !(is.na(tmp1[[i]][[j]]) | is.na(tmp2[[i]][[j]]) | abs(tmp1[[i]][[j]]) == Inf | abs(tmp2[[i]][[j]]) == Inf)
		tmp0	<- t.test(tmp2[[i]][[j]][sel], tmp1[[i]][[j]][sel])
		tmp		<- rbind(tmp, c(tmp0$estimate, tmp0$p.value))
	}
}
dimnames(tmp)	<- list(paste(rep(names(sim.ls07), each = nc-1), dpt.name[-beta0.base], sep=":"), c("2008", "2007", "pvalue"))
cat("T test of mean elasticity:\n"); print(tmp); cat("\n")

# Overall Revenue
tmp0	<- colSums(base.07[,gsub("\\s", ".",fmt_name)])
tmp1	<- do.call(rbind, sim.ls07)
tmp1	<- split(tmp1, paste(tmp1$Var, tmp1$department, sep="*"))
tmp1	<- lapply(tmp1, function(x) colSums(x[,gsub("\\s", ".",fmt_name)]))
tmp1	<- sapply(tmp1, function(x) (x[sel.retail] - tmp0[sel.retail])/tmp0[sel.retail]/iv.change)

tmp0	<- colSums(base.08[,gsub("\\s", ".",fmt_name)])
tmp2	<- do.call(rbind, sim.ls08)
tmp2	<- split(tmp2, paste(tmp2$Var, tmp2$department, sep="*"))
tmp2	<- lapply(tmp2, function(x) colSums(x[,gsub("\\s", ".",fmt_name)]))
tmp2	<- sapply(tmp2, function(x) (x[sel.retail] - tmp0[sel.retail])/tmp0[sel.retail]/iv.change)
tmp.tab	<- cbind(tmp1, tmp2)
colnames(tmp.tab)	<- c("2007", "2008")
cat("Elasticity on overall revenue:\n"); print(tmp.tab); cat("\n")

# Market share
tmp0	<- colSums(base.07[,gsub("\\s", ".",fmt_name)])
tmp0	<- tmp0/sum(tmp0)
tmp1	<- do.call(rbind, sim.ls07)
tmp1	<- split(tmp1, paste(tmp1$Var, tmp1$department, sep="*"))
tmp1	<- lapply(tmp1, function(x) colSums(x[,gsub("\\s", ".",fmt_name)]))
tmp1	<- lapply(tmp1, function(x) x/sum(x))

tmp02	<- colSums(base.08[,gsub("\\s", ".",fmt_name)])
tmp02	<- tmp02/sum(tmp02)
tmp2	<- do.call(rbind, sim.ls08)
tmp2	<- split(tmp2, paste(tmp2$Var, tmp2$department, sep="*"))
tmp2	<- lapply(tmp2, function(x) colSums(x[,gsub("\\s", ".",fmt_name)]))
tmp2	<- lapply(tmp2, function(x) x/sum(x))
tmp.tab1	<- cbind(sapply(tmp1, function(x) (x[sel.retail] - tmp0[sel.retail])), 
					 sapply(tmp2, function(x) (x[sel.retail] - tmp02[sel.retail])) )
tmp.tab2	<- cbind(sapply(tmp1, function(x) (x[sel.retail] - tmp0[sel.retail])/tmp0[sel.retail]/iv.change), 
					 sapply(tmp2, function(x) (x[sel.retail] - tmp02[sel.retail])/tmp02[sel.retail]/iv.change))
colnames(tmp.tab1)	<- colnames(tmp.tab2) <- c("2007", "2008")
cat("Change in market share:\n"); print(tmp.tab1); cat("\n")
cat("Elasticity on market share:\n"); print(tmp.tab2); cat("\n")

################
# Save results #
################
rm(list  = intersect(ls(), c("ar", "args", "cl", "ggtmp", "ii", "i", "make_plot", "mycore", "cpi.adj","sel1", 
			"myfix", "plot.wd", "sel", "selyr", "price", "sol", "sol.top", "sol.top2", "tmp", "tmp_coef", 
			"tmp_price", "tmp1", "tmp2", "use.time", "ver.date", "var1", "ww", "f", 
			"pct", "tmpd1", "tmpd2", "tmpdat", "out", "out1", "tmp1", "tmp2", "tmp3","tmp.ls","tmp0",
			"Allocation_constr_fn", "Allocation_fn", "Allocation_nlop_fn", "cheb.1d.basis", "cheb.basis", "chebfun", 
			"exp_fn", "expFOC_fn", "incl_value_fn", "mysplfun", "mytrimfun", "param_assignR", "simExp_fn", "SimOmega_fn", 
			"SimWrapper_fn", "solveExp_fn", "spl2dfun", "uP_fn", "uGrad_fn", "attr_sel", "attr_mat", "attr.new", "gamma.new", 
			"j","loadf", "price.new", "uPGrad_fn")))

save.image(paste("estrun_",run_id, "/", fname, ".rdata",sep=""))

cat("This program is done.\n")

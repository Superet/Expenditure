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
fname			<- "retcat_newformat"
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
dist	<- dist[order(dist$household_code),]
cat("dim(dist) =", dim(dist), ", dim(sim.unq)*R =",nrow(sim.unq)*R, ".\n")
beta0.base	<- which(dpt.name == "OTHER")

sel			<- fmt_dpt$year == selyr
attr_mat	<- merge(sim.unq[,c("household_code","scantrack_market_descr")], 
				fmt_dpt[sel,c("scantrack_market_descr", "channel_type","department","unitprice_paid",selcol)], 
				by = c("scantrack_market_descr"), all.x=T)
attr_mat	<- attr_mat[order(attr_mat$household_code,attr_mat$channel_type, attr_mat$department),]	
attr_mat$department <- factor(attr_mat$department, levels = dpt.lab, labels = dpt.name)
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
sim.base08		<- foreach(i = 1:nseg, .errorhandling = "remove", .packages=c("mgcv", "nloptr")) %dopar% {
	sel		<- sim.unq$segment == i
	out		<- SimWrapper_fn(lambda[,i], Inc = sim.unq[sel,"income"], X = lapply(X1.08, function(x) x[sel,]), 
							gamma = gamma.08[sel,], eps_draw = eps_draw, price = NULL, ret.sim = TRUE)
	out$Average	<- cbind(segment = i, household_code = sim.unq[sel,"household_code"], out$Average)
	return(out)
}
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
eps_draw_new	<- cbind(eps_draw, rgev(numsim)) 

cat(" ---------------------------------------------------------------------- \n")
pct		<- proc.time()

# Set new parameters
sel			<- grep(ret.a, rownames(lambda))
lambda.new	<- matrix(0, nrow(lambda)+2, ncol(lambda))
lambda.new[setdiff(1:nrow(lambda.new),c(R+1, 2*R+2)),]	<- lambda
lambda.new[c(R+1, 2*R+2),]	<- lambda[sel,]
cat("The parameter vector for adding a new format from", ret.a, "is:\n"); print(lambda.new); cat("\n")
fmt_name1<- c(fmt_name, "New")
delta	<- seq(0, 1, .1)			# price = (1-delta)*discount + delta*convenience 

# Price level for discount stores 
sel		<- attr_mat$channel_type == ret.a
tmp		<- dcast(attr_mat[sel,c("household_code", "department", "unitprice_paid")], 
				household_code ~ department, value.var = "unitprice_paid")
names(tmp)	<- c("household_code", paste("PRC_", dpt.name, sep=""))
price	<- as.matrix(tmp[,paste("PRC_",dpt.name, sep="")])

pct		<- proc.time()
for(j in 1:length(sim.scn)){
	# Set new attributes
	X.08.new	<- setNames(vector("list", R+2), c(fmt_name1, "Outside"))
	for(i in 1:R){
		tmp		<- cbind(X1.08[[i]], 0, 0)
		colnames(tmp)	<- c(colnames(X1.08[[1]]), "Intercept:New", "Rec:New")
		X.08.new[[i]]	<- tmp[,c(grep("Intercept", colnames(tmp)), grep("Rec:", colnames(tmp)), grep("distance", colnames(tmp)))]
	}

	# Set the model matrix in the upper model for the new format
	tmp0	<- rep(0, R+1)
	tmp0[(R+1)]	<- 1
	tmp1	<- as.matrix(cbind(1, rec=rep(1, nrow(sim.unq))))	
	tmp2	<- kronecker(tmp1, t(tmp0))
	colnames(tmp2)	<- paste(rep(c("Intercept", "Rec"), each = R+1), fmt_name1, sep=":")	
	sel1	<- dist$channel_type == sim.scn[j]
	tmp		<- as.matrix(dist[sel1,"distance"])
	X.08.new[[(R+1)]]	<- cbind(tmp2, distance = tmp, tmp*1)
	colnames(X.08.new[[(R+1)]])	<- c(colnames(tmp2), "distance", "distance:Rec")
	X.08.new[[(R+2)]]	<- matrix(0, nrow(sim.unq), ncol(X.08.new[[1]]))
	
	# Set the retail attributes in the lower model for the format
	sel			<- attr_mat$channel_type == sim.scn[j]
	tmp			<- attr_mat[sel,]
	tmp[,"prvt_overall"]	<- attr_mat[attr_mat$channel_type==ret.a,"prvt_overall"]
	attr.new	<- as.matrix(do.call(cbind,lapply(dpt.name[-beta0.base], function(x) tmp[tmp$department==x,selcol])))
	
	# Price level at the counterfactual format 
	sel		<- attr_mat$channel_type == sim.scn[j]
	tmp		<- dcast(attr_mat[sel,c("household_code", "department", "unitprice_paid")], 
					household_code ~ department, value.var = "unitprice_paid")
	names(tmp)	<- c("household_code", paste("PRC_", dpt.name, sep=""))
	price1	<- as.matrix(tmp[,paste("PRC_",dpt.name, sep="")])
	
	# Set prices 
	newf.sim[[j]]	<- foreach(d = delta, .errorhandling = "remove", .packages=c("mgcv", "nloptr")) %dopar% {
		gamma.new	<- cbind(gamma.08[,1:R], New = 0, Outside = 1)
		price.new	<- (1-d) * price + d*price1
		out1		<- data.frame()
		for(i in 1:nseg){
			sel		<- sim.unq$segment == i
			tmp		<- data.frame(y = rep(0, sum(sel)))
			tmp$price	<- price.new[sel,-beta0.base]
			tmp$X		<- as.matrix(attr.new[sel,])
			tmp$rec		<- 1
			gamma.new[sel,(R+1)]	<- exp(predict(omega.lsf[[i]][[which(fmt_name == ret.a)]], newdata = tmp))
			tmp		<- SimWrapper_fn(lambda.new[,i], Inc = sim.unq[sel,"income"], X = lapply(X.08.new, function(x) x[sel,]), 
									gamma = gamma.new[sel,], eps_draw = eps_draw_new, price = NULL, ret.sim = FALSE)
			tmp		<- data.frame(d = d, household_code = sim.unq[sel,"household_code"], tmp)	
			out1	<- rbind(out1, tmp)
		}
		return(out1)
	}
}
use.time	<- proc.time() - pct
cat("Counterfactual with", ret.a, "finishes with", use.time[3]/60, "min.\n")
stopCluster(cl)
cat("Stop clustering. \n")

# Compare the difference 
tmp1	<- do.call(rbind, lapply(sim.base08, function(x) x$Average[,fmt_name]))
tmp1	<- c(colSums(tmp1, na.rm=T), 0)
tmp2	<- colSums(newf.sim[[1]][[1]][,gsub("\\s", ".",fmt_name1)], na.rm=T)
mkt.dif	<- setNames(tmp2/sum(tmp2)-tmp1/sum(tmp1), fmt_name1)
cat("Difference of market share:\n"); print(round(mkt.dif*100,2)); cat("\n")

tmp2	<- colSums(newf.sim[[2]][[1]][,gsub("\\s", ".",fmt_name1)], na.rm=T)
mkt.dif	<- setNames(tmp2/sum(tmp2)-tmp1/sum(tmp1), fmt_name1)
cat("Difference of market share in scenario", sim.scn[2], "\n"); print(round(mkt.dif*100,2)); cat("\n")

# Calculat the sum of market share of discount stores and new format
ggtmp	<- data.frame()
for(i in 1:length(delta)){
	tmp	<- colSums(newf.sim[[1]][[i]][,gsub("\\s", ".",fmt_name1)], na.rm=T)
	tmp	<- tmp/sum(tmp)
	ggtmp	<- rbind(ggtmp, data.frame(d = newf.sim[[1]][[i]][1,"d"], share = tmp[which(fmt_name==ret.a)] + tmp["New"]))
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
			"SimWrapper_fn", "solveExp_fn", "spl2dfun", "uP_fn", "uPGrad_fn", "X_list", "beta0_7", "d.psi", "gamfit","j", "k",
			"attr_mat", "attr.new", "dist", "price1", "X.08.new", "gamma.new")))

save.image(paste("estrun_",run_id, "/", fname, ".rdata",sep=""))

cat("This program is done.\n")


# # Plot the counterfactual results
# # Overall market share
# bar.col	<- "#4F81BD"
# 
# fmt_name1	<- c(fmt_name, "New")
# fmt_lab1	<- c(fmt_name, "Discount Store Express")
# fmt_lab1[fmt_lab1=="Grocery"]	<- "Grocery Store"
# tmp1	<- do.call(rbind, lapply(sim.base08, function(x) x$Average[,fmt_name]))
# tmp1	<- c(colSums(tmp1, na.rm=T), 0)
# tmp2	<- colSums(newf.sim[[1]][[1]][,gsub("\\s", ".",fmt_name1)], na.rm=T)
# ggtmp	<- data.frame(dif = tmp2/sum(tmp2)-tmp1/sum(tmp1), retailer = fmt_name1)
# ord 	<- order(ggtmp[,1])
# ggtmp$retailer <- factor(ggtmp$retailer, levels =fmt_name1[ord], labels = fmt_lab1[ord])
# text.col	<- rep("black", R+1)
# text.col[c(2,7)]	<- bar.col
# ggtmp$colind	<- ifelse(ggtmp$retailer %in% c("Discount Store Express", "Discount Store"), 1, 0) 
# 
# plots	<- list(NULL)
# plots[[1]]	<- ggplot(ggtmp, aes(retailer, dif)) + 
# 				geom_bar(stat = "identity") + 
# 				scale_y_continuous(labels = scales::percent) +
# 				xlab("") + ylab("Percentage point change")+ coord_flip() + 
# 				theme_bw() + theme(legend.position = "none")
# 
# plots[[2]]	<- ggplot(ggtmp, aes(retailer, dif, fill = factor(colind))) + 
# 				geom_bar(stat = "identity") + 
# 				scale_y_continuous(labels = scales::percent) +
# 				scale_fill_manual(values = c("grey70", bar.col)) + 
# 				xlab("") + ylab("Percentage point change")+ coord_flip() + 
# 				theme_bw() + theme(legend.position = "none", axis.text.y = element_text(colour=text.col))
# 
# plots[[3]]	<- ggplot(ggtmp, aes(retailer, dif, fill = factor(colind), alpha = factor(colind))) + 
# 				geom_bar(stat = "identity", fill = bar.col) + 
# 				scale_y_continuous(labels = scales::percent) +
# 				scale_fill_manual(values = c("grey70", bar.col)) +
# 				scale_alpha_manual(values = c(0, 1), guide='none') + 
# 				xlab("") + ylab("Percentage point change")+ coord_flip() + 
# 				theme_bw() + theme(legend.position = "none", axis.text.y = element_text(colour=text.col))
# 
# pdf(paste(plot.wd, "/graph_newformat.pdf",sep=""), width =4.5, height = 4.5*ar)
# for(i in 1:length(plots)){ print(plots[[i]])}
# dev.off()
# 
# # Price range
# ggtmp	<- do.call(rbind, newf.sim[[1]])
# tmp		<- split(ggtmp, ggtmp$d)
# tmp		<- lapply(tmp, function(x) colSums(x[,gsub("\\s", ".",fmt_name1)]))
# tmp		<- lapply(tmp, function(x) x/sum(x))
# ggtmp	<- sapply(tmp, function(x) x - tmp1)
# ggtmp1	<- ggtmp["Discount.Store",] + ggtmp["New",]
# ggtmp1	<- data.frame(d = as.numeric(names(ggtmp1)), Total_share = ggtmp1)
# 
# pdf(paste(plot.wd, "/graph_newformat_pricing.pdf",sep=""), width =4.5, height = 4.5*ar)
# ggplot(ggtmp1, aes(d, Total_share)) + geom_line() + 
# 		scale_y_continuous(labels = scales::percent) + 
# 		labs(x = expression(delta), y = "Net share gain") + 
# 		theme_bw()
# dev.off()
# 

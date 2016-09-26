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

run_id		<- 1
plot.wd		<- paste(getwd(), "/estrun_", run_id, sep="")
make_plot	<- TRUE
ww			<- 10
ar			<- .6

source("0_Efficient_Allocation_function.R")
source("ctrfact_sim_functions.r")

# Load estimation data 
ver.date	<- "2016-09-11"
loadf	<- paste("estrun_",run_id,"/MDCEV_est_", ver.date,".rdata",sep="")
loadf
load(loadf)

# Set simulation parameters
week.price		<- FALSE
interp.method	<- "spline"					# Spline interpolation or Chebyshev interpolation
exp.method		<- "Utility"				# Utility maximization or finding roots for first order condition
trim.alpha		<- 0
numsim			<- 500
draw.par		<- FALSE
sim.omega		<- FALSE
fname			<- paste("newfmt_sim", numsim, "_", as.character(Sys.Date()), sep="")
cat("Output file name is", fname, ".\n")

###############################
# Prepare simulation elements # 
###############################
# Data required: parameters, income level, price, retail attributes, and random draws
# For each simulation scenario, if income, price or retail attributes change, we need to re-simulate inclusive values.

# Set simulation parameters
lambda 	<- sapply(exp.est.ls, coef)
shr.par	<- sapply(shr.est.ls, coef)
ret.idx	<- 11:15		# Column index flagging the retail attribugtes in X1
nx0		<- 2			# Number of demograhic var in X1


#-----------------------#
# Construct income data # 
selyr	<- 2007
sel		<- mydata.full$year == selyr & !duplicated(mydata.full[,c("household_code", "year")]) 
mydata.full$ln_inc	<- log(mydata.full$income_mid)
sim.unq	<- mydata.full[sel,c("household_code", "ln_inc","rec", "segment")]
sim.unq$price	<- as.matrix(mydata.full[sel,paste("PRC_", gsub("\\s", "_", fmt_name), sep="")])
# sim.unq	<- sim.unq[sample(1:nrow(sim.unq), 100),]

# Normalize distance by zipcode
tmp		<- data.table(hh_dist[,c("household_code", "panelist_zip_code","year", "channel_type", selcol)])
tmp		<- tmp[, maxd:=max(distance_haver, na.rm=T), by = list(panelist_zip_code)]
tmp		<- tmp[, distance_haver := distance_haver/maxd]
tmp		<- data.frame(subset(tmp, year == selyr))
attr_mat	<- merge(sim.unq[,c("household_code", "ln_inc")], tmp[,c("household_code", "channel_type", selcol)], 
				by = c("household_code"), all.x=T)
attr_mat	<- attr_mat[order(attr_mat$household_code),]	
cat("dim(sim.unq) =", dim(sim.unq), ".\n")
cat("dim(attr_mat) =", dim(attr_mat), ", dim(sim.unq)*R =",nrow(sim.unq)*R, ".\n")
beta0_base 	<- which(fmt_name == "Grocery")

X2	<- X1 	<- vector("list", length=length(fmt_name))
for(i in 1:length(fmt_name)){
	tmp0<- rep(0, R-1)
	if(i<beta0_base){
		tmp0[i]	<- 1
	}else if(i > beta0_base){
		tmp0[(i-1)]	<- 1
	}
	tmp1<- as.matrix(cbind(1, rec=sim.unq$rec))	
	tmp2<- kronecker(tmp1, t(tmp0))
	colnames(tmp2)	<- paste(rep(c("Intercept", "Rec"), each = R-1), fmt_name[-beta0_base], sep=":")
	tmp3<- rep(0, R)
	tmp3[i]	<- 1
	tmp3<- kronecker(tmp1, t(tmp3))
	colnames(tmp3)	<- paste(rep(c("Intercept", "Rec"), each = R), fmt_name, sep=":")
	sel	<- attr_mat$channel_type == fmt_name[i]	
	tmp	<- as.matrix(attr_mat[sel,selcol])
	X1[[i]]	<- cbind(tmp2, tmp, tmp*sim.unq$rec)
	X2[[i]]	<- cbind(tmp3, tmp, tmp*sim.unq$rec)
	colnames(X1[[i]])	<- c(colnames(tmp2), colnames(tmp), paste(colnames(tmp), ":Rec", sep=""))
	colnames(X2[[i]])	<- c(colnames(tmp3), colnames(tmp), paste(colnames(tmp), ":Rec", sep=""))
}
nx 		<- ncol(X1[[1]])
cat("nx =", nx,"\n")

sim.unq$X<- as.matrix( do.call(cbind, lapply(X1, function(x) x[,ret.idx])))
cat("dim(sim.unq) =", dim(sim.unq), "\n")

# Retail attributes for simulation
X1.07	<- X1
X2.07	<- X2
X2.08	<- X1.08 	<- vector("list", length=length(fmt_name))
for(i in 1:length(fmt_name)){
	tmp0<- rep(0, R-1)
	if(i<beta0_base){
		tmp0[i]	<- 1
	}else if(i > beta0_base){
		tmp0[(i-1)]	<- 1
	}
	tmp1<- as.matrix(cbind(1, rec=rep(1, nrow(sim.unq))))	
	tmp2<- kronecker(tmp1, t(tmp0))
	colnames(tmp2)	<- paste(rep(c("Intercept", "Rec"), each = R-1), fmt_name[-beta0_base], sep=":")
	tmp3<- rep(0, R)
	tmp3[i]	<- 1
	tmp3<- kronecker(tmp1, t(tmp3))
	colnames(tmp3)	<- paste(rep(c("Intercept", "Rec"), each = R), fmt_name, sep=":")
	sel1	<- attr_mat$channel_type == fmt_name[i]
	tmp	<- as.matrix(attr_mat[sel1,selcol])
	X1.08[[i]]	<- cbind(tmp2, tmp, tmp*1)
	X2.08[[i]]	<- cbind(tmp3, tmp, tmp*1)
	colnames(X1.08[[i]])	<- c(colnames(tmp2), colnames(tmp), paste(colnames(tmp), ":Rec", sep=""))
	colnames(X2.08[[i]])	<- c(colnames(tmp3), colnames(tmp), paste(colnames(tmp), ":Rec", sep=""))
}

#-------------------#
# Take random draws #
set.seed(666)
eps_draw	<- lapply(1:nseg, function(i) matrix(rgev(numsim*R, scale = exp(shr.par["ln_sigma",i])), numsim, R))
par_draw <- NULL

rm(list = intersect(ls(), c("hh_month", "hh_dist", "hh_dpt")))

##############
# Simulation #
##############
# Register parallel computing
mycore 	<- 2
cl		<- makeCluster(mycore, type = "FORK")
# cl		<- makeCluster(mycore, type = "PSOCK")			# Register cores in Windows
registerDoParallel(cl)
cat("Register", mycore, "core parallel computing. \n")

# ------------------ #
# Basline simulation #
# Simulate expenditure and expenditure share. 
# pct				<- proc.time()
# sim.base07		<- foreach(i = 1:nseg, .errorhandling = "remove", .packages=c("mgcv", "nloptr")) %dopar% {
# 	sel		<- sim.unq$segment == i
# 	out		<- SimWrapper_fn(omega.ls[[i]], lambda[,i], shr.par[,i], base = beta0_base, 
# 				X1 = lapply(X1.07, function(x) x[sel,]), X2=lapply(X2.07, function(x) x[sel,]), price = sim.unq$price[sel,], 
# 				demo_mat = as.matrix(cbind(Intercept = 1, rec = 0, ln_inc = sim.unq[sel,"ln_inc"])), 
# 				eps_draw = eps_draw[[i]], ret.idx = ret.idx, sim.y = NULL, ret.sim = TRUE, method = exp.method)
# 	out$Average	<- cbind(segment = i, household_code = sim.unq[sel,"household_code"], out$Average)
# 	return(out)
# }
# use.time		<- proc.time() - pct						
# cat("2007 Baseline simulation finishes with", use.time[3]/60, "min.\n")

pct				<- proc.time()
sim.base08		<- foreach(i = 1:nseg, .errorhandling = "remove", .packages=c("mgcv", "nloptr")) %dopar% {
	sel		<- sim.unq$segment == i
	out		<- SimWrapper_fn(omega.ls[[i]], lambda[,i], shr.par[,i], base = beta0_base, 
				X1 = lapply(X1.08, function(x) x[sel,]), X2=lapply(X2.08, function(x) x[sel,]), price = sim.unq$price[sel,], 
				demo_mat = as.matrix(cbind(Intercept = 1, rec = 1, ln_inc = sim.unq[sel,"ln_inc"])), 
				eps_draw = eps_draw[[i]], ret.idx = ret.idx, sim.y = NULL, ret.sim = TRUE, method = exp.method)
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
eps_draw_new	<- lapply(1:nseg, function(i) cbind(eps_draw[[i]], rgev(numsim, scale = exp(shr.par["ln_sigma",i]))) )

cat(" ---------------------------------------------------------------------- \n")
pct		<- proc.time()

# Set new parameters
shr.par.new	<- NULL
for(i in 1:nseg){
	sel		<- which(fmt_name == ret.a)
	tmp		<- shr.par[1:(nx0*(R-1)),i]
	beta1	<- rbind(matrix(tmp, R-1, nx0), tmp[grep(ret.a, names(tmp))])
	tmp		<- shr.par[(nx+1):(nx0*R+nx),i]
	beta2	<- rbind(matrix(tmp, R, nx0), tmp[grep(ret.a, names(tmp))])
	shr.par.new	<- cbind(shr.par.new, 
						c(c(beta1), shr.par[(nx0*(R-1)+1):nx,i], c(beta2), shr.par[-c(1:(nx0*R+nx)),i]))
}
cat("The parameter vector for adding a new format from", ret.a, "is:\n"); print(shr.par.new); cat("\n")

fmt_name1	<- c(fmt_name, "New")
fmt_lab1	<- c(fmt_name, "Discount Store Express")
delta		<- seq(0, 1, .1)			# price = (1-delta)*discount + delta*convenience 

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
		tmp	<- as.matrix(attr_mat[sel1,selcol])
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
	tmp	<- as.matrix(attr_mat[sel1,selcol])
	sel1	<- attr_mat$channel_type == ret.a
	tmp[,"prvt_overall"]	<- attr_mat[sel1,"prvt_overall"]
	X1.08.new[[(R+1)]]	<- cbind(tmp2, tmp, tmp*1)
	X2.08.new[[(R+1)]]	<- cbind(tmp3, tmp, tmp*1)
	colnames(X1.08.new[[(R+1)]])	<- c(colnames(tmp2), colnames(tmp), paste(colnames(tmp), ":Rec", sep=""))
	colnames(X2.08.new[[(R+1)]])	<- c(colnames(tmp3), colnames(tmp), paste(colnames(tmp), ":Rec", sep=""))
	
	# Set prices 
	newf.sim[[j]]	<- foreach(d = delta, .errorhandling = "remove", .packages=c("mgcv", "nloptr"), .combine = "rbind") %dopar% {
		price.new	<- sim.unq$price
		price.new	<- cbind(price.new, 
			New = (1-d) * sim.unq$price[,which(fmt_name==ret.a)] + d*sim.unq$price[,which(fmt_name=="Convenience Store")])
		
		out1	<- data.frame()	
		for(i in 1:nseg){
			sel		<- sim.unq$segment == i
			out		<- SimWrapper_fn(omega.ls[[i]], lambda[,i], shr.par.new[,i], base = beta0_base, 
						X1 = lapply(X1.08.new, function(x) x[sel,]), X2=lapply(X2.08.new, function(x) x[sel,]), 
						price = price.new[sel,], 
						demo_mat = as.matrix(cbind(Intercept = 1, rec = 1, ln_inc = sim.unq[sel,"ln_inc"])), 
						eps_draw = eps_draw_new[[i]], ret.idx = ret.idx+nx0, sim.y = sim.base08[[i]]$y, ret.sim = FALSE, method = exp.method)
			out		<- cbind(segment = i, household_code = sim.unq[sel,"household_code"], d = d, out)
			out1	<- rbind(out, out1)
		}	
		return(out1)			
	}
}
use.time	<- proc.time() - pct
cat("Counterfactual with", ret.a, "finishes with", use.time[3]/60, "min.\n")
stopCluster(cl)
cat("Stop clustering. \n")

################################
# Summarize the counterfactual #
################################
base.08	<- do.call(rbind, lapply(sim.base08, function(x) x$Average))
bar.col	<- "#4F81BD"

# Compare the difference 
tmp0	<- colSums(cbind(base.08[,fmt_name],0), na.rm=T)
tmp0	<- tmp0/sum(tmp0)
tmp1	<- lapply(newf.sim, function(x) colSums(x[x$d==0,fmt_name1], na.rm=T))
tmp1	<- lapply(tmp1, function(x) x/sum(x))
ggtmp	<- sapply(tmp1, function(x) x-tmp0)
colnames(ggtmp)	<- gsub("Discount Store_", "", colnames(ggtmp))
cat("Difference of market share:\n"); print(round(ggtmp*100,2)); cat("\n")
ord 	<- order(ggtmp[,1])
ggtmp	<- melt(ggtmp)
ggtmp$retailer <- factor(ggtmp$Var1, levels = fmt_name1[ord], labels = fmt_lab1[ord])
text.col	<- rep("black", R+1)
sel			<- grep("Discount", fmt_lab1[ord])
text.col[sel]	<- bar.col
ggtmp$colind	<- ifelse(ggtmp$retailer %in% fmt_lab1[c(2,(R+1))], 1, 0)

plots	<- list(NULL)
plots[[1]]	<- ggplot(subset(ggtmp, Var2 == "Convenience Store"), aes(retailer, value)) + 
				geom_bar(stat = "identity") + 
				scale_y_continuous(labels = scales::percent) +
				xlab("") + ylab("Percentage point change")+ coord_flip() + 
				theme_bw() + theme(legend.position = "none")

plots[[2]]	<- ggplot(subset(ggtmp, Var2 == "Convenience Store"), aes(retailer, value, fill = factor(colind))) + 
				geom_bar(stat = "identity") + 
				scale_y_continuous(labels = scales::percent) +
				scale_fill_manual(values = c("grey70", bar.col)) + 
				xlab("") + ylab("Percentage point change")+ coord_flip() + 
				theme_bw() + theme(legend.position = "none", axis.text.y = element_text(colour=text.col))

plots[[3]]	<- ggplot(subset(ggtmp, Var2 == "Convenience Store"), aes(retailer, value, fill = factor(colind), alpha = factor(colind))) + 
				geom_bar(stat = "identity", fill = bar.col) + 
				scale_y_continuous(labels = scales::percent) +
				scale_fill_manual(values = c("grey70", bar.col)) +
				scale_alpha_manual(values = c(0, 1), guide='none') + 
				xlab("") + ylab("Percentage point change")+ coord_flip() + 
				theme_bw() + theme(legend.position = "none", axis.text.y = element_text(colour=text.col))

pdf(paste(plot.wd, "/graph_newformat.pdf",sep=""), width =4.5, height = 4.5*ar)
for(i in 1:length(plots)){ print(plots[[i]])}
dev.off()

# Price range
tmp		<- split(newf.sim[["Discount Store_Convenience Store"]], newf.sim[["Discount Store_Convenience Store"]][,"d"])
tmp		<- lapply(tmp, function(x) colSums(x[, fmt_name1]))
tmp		<- lapply(tmp, function(x) x/sum(x))
ggtmp	<- sapply(tmp, function(x) x - tmp0)
ggtmp1	<- ggtmp["Discount Store",] + ggtmp["New",]
ggtmp1	<- data.frame(d = as.numeric(names(ggtmp1)), Total_share = ggtmp1)

pdf(paste(plot.wd, "/graph_pricing_", fname, ".pdf",sep=""), width =4.5, height = 4.5*ar)
ggplot(ggtmp1, aes(d, Total_share)) + geom_line() + 
		scale_y_continuous(labels = scales::percent) + 
		labs(x = expression(delta), y = "Net share gain") + 
		theme_bw()
dev.off()

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



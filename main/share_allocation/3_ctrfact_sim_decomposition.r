library(ggplot2)
library(reshape2)
library(maxLik)
library(evd)
library(data.table)
library(doParallel)
library(foreach)
library(nloptr)
library(mgcv)
library(scales)
options(error = quote({dump.frames(to.file = TRUE)}))

seg_id	<- as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))
cat("seg_id =", seg_id, ".\n")

args <- commandArgs(trailingOnly = TRUE)
print(args)
if(length(args)>0){
    for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
    }
}

# setwd("~/Documents/Research/Store switching/processed data")
# plot.wd	<- '~/Desktop'
# source("../Exercise/Multiple_discrete_continuous_model/0_Allocation_function.R")
# source("../Exercise/main/share_allocation/ctrfact_sim_functions_v2.r")

# setwd("/home/brgordon/ccv103/Exercise/run")
# setwd("/kellogg/users/marketing/2661703/Exercise/run")
setwd("/sscc/home/c/ccv103/Exercise/run")
run_id		<- 4
plot.wd		<- getwd()
make_plot	<- TRUE
ww			<- 6.5
ar			<- .6

# Load estimation data 
ver.date	<- "2016-02-26"
cpi.adj		<- TRUE

if(cpi.adj){
	loadf	<- paste("estrun_",run_id,"/MDCEV_cpi_est_seg",seg_id,"_", ver.date,".rdata",sep="")
}else{
	loadf	<- paste("estrun_",run_id,"/MDCEV_est_seg",seg_id,"_", ver.date,".rdata",sep="")
}
loadf
load(loadf)
rm(list = intersect(ls(), c("gamfit", "shr","model_name", "tmpdat")))
source("0_Allocation_function.R")
source("ctrfact_sim_functions_v2.r")

# Set simulation parameters
interp.method	<- "spline"					# Spline interpolation or Chebyshev interpolation
exp.method		<- "Utility"				# Utility maximization or finding roots for first order condition
trim.alpha		<- 0.05
numsim			<- 1000	#numsim1		#<- 1000
draw.par		<- FALSE
sim.omega		<- FALSE
sub.hh			<- TRUE					# Subset households who stay in the panel for all the year. 
fname			<- paste("ctrfact_dcps_seg",seg_id,sep="")
if(draw.par){ fname <- paste(fname, "_pardraw", sep="") }
if(sim.omega){fname	<- paste(fname, "_simomega", sep="") }
fname			<- paste(fname, "_sim", numsim, "_", as.character(Sys.Date()), sep="")
cat("Output file name is", fname, ".\n")

# A randome 5% subsample 
# tmp		<- unique(mydata$household_code)
# sel		<- sample(tmp, ceiling(length(tmp)*.05))
# mydata.save	<- mydata
# mydata	<- subset(mydata, household_code %in% sel)

###############################
# Prepare simulation elements # 
###############################
# Data required: parameters, income level, price, retail attributes, and random draws
# For each simulation scenario, if income, price or retail attributes change, we need to re-simulate inclusive values.

# Set simulation parameters
lambda 	<- coef(sol.top2)
shr.par	<- coef(sol)
rm(list = c("sol", "sol.top", "sol.top2"))

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

# ------------------------------- #
# Set up baseline simulation data #
# NOTE that income has been adjusted by CPI if cpi.adj = TRUE
sim.data	<- mydata[,c("household_code", "year", "scantrack_market_descr", "income_midvalue")]
names(sim.data)		<- c("household_code", "year", "scantrack_market_descr", "income")
if(sub.hh){
	tmp		<- data.table(sim.data)
	tmp		<- tmp[, list(tenure = max(year) - min(year)), by = list(household_code)]
	cat(sum(tmp$tenure == max(tmp$tenure)), "households (", sum(tmp$tenure == max(tmp$tenure))/nrow(tmp), 
			") stay in the panel over the entire periods.\n")
	tmp		<- subset(tmp, tenure == max(tmp$tenure))		
	sim.data<- subset(sim.data, household_code %in% tmp$household_code)
}
sim.data	<- data.table(sim.data)
sim.data	<- unique(sim.data, by = c("household_code", "scantrack_market_descr", "year", "income"))
cat("dim(sim.data) =", dim(sim.data), "\n")

# Reduce to the unique combination of makret-year-income
sim.unq0		<- data.table(sim.data)
sim.unq0		<- unique(sim.unq0, by = c("year", "scantrack_market_descr", "income"))
sim.unq0		<- sim.unq0[,!"household_code", with = FALSE]
setkeyv(sim.unq0, c("scantrack_market_descr", "year", "income"))
cat("dim(sim.unq0) =", dim(sim.unq0), "\n")

# ---------------------------------------------------------------------------------------- #
# Scenario 1: retail attributes at the actual level but fix income at the first year level #
sim.data1	<- data.table(sim.data)
setkeyv(sim.data1, c("household_code", "year"))
sim.data1	<- sim.data1[,income := income[1], by = list(household_code)]		# Set the income at the level of first year

# Reduce to the unique combination of makret-year-income
sim.unq1		<- data.table(sim.data1)
sim.unq1		<- unique(sim.unq1, by = c("year", "scantrack_market_descr", "income"))
sim.unq1		<- sim.unq1[,!"household_code", with = FALSE]
setkeyv(sim.unq1, c("scantrack_market_descr", "year", "income"))
cat("dim(sim.unq1) =", dim(sim.unq1), "\n")

# ------------------------------------------------------------------------------ #
# Scenarior 2: retail attributes fixed at year 1 level but income vary over year #
sim.data2	<- data.table(sim.data)
sim.data2	<- sim.data2[, year0 := year]
sim.data2	<- sim.data2[,year:=2004]		# Set year at 2004

# Reduce to the unique combination of makret-year-income
sim.unq2		<- data.table(sim.data2)
sim.unq2		<- unique(sim.unq2, by = c("year", "scantrack_market_descr", "income"))
sim.unq2		<- sim.unq2[,!"household_code", with = FALSE]
setkeyv(sim.unq2, c("scantrack_market_descr", "year", "income"))
cat("dim(sim.unq2) =", dim(sim.unq2), "\n")

# ---------------------------------- # 
# Combine the unique simulation data # 
sim.unq		<- rbind(sim.unq0, sim.unq1, sim.unq2[,!"year0",with=FALSE])
sim.unq		<- unique(sim.unq, by = c("scantrack_market_descr", "year", "income"))
setkeyv(sim.unq, c("scantrack_market_descr", "year", "income"))
sim.unq		<- sim.unq[,idx:=1:nrow(sim.unq)]
cat("dim(sim.unq) =", dim(sim.unq), ".\n")

# Set up retail attributes 
selcol	<- c("size_index", "ln_upc_per_mod", "ln_num_module","overall_prvt")
X		<- setNames(vector("list", R), fmt_name)
for(i in 1:R){
	tmp	<- merge(sim.unq, subset(fmt_attr, channel_type == fmt_name[i]), by = c("scantrack_market_descr", "year"), all.x = TRUE)
	X[[i]]	<- as.matrix(tmp[,selcol,with = FALSE])
}

# Set up average price index
tmp 		<- dcast(price_dat,	scantrack_market_descr + year ~ channel_type, value.var = "bsk_price_paid_2004", fun.aggregate = mean) 
p			<- merge(sim.unq, tmp, by = c("scantrack_market_descr", "year"), all.x = TRUE)[, fmt_name, with = FALSE]
p			<- as.matrix(p)
cat("Check any missing values in p: sum(is.na(p)) =", sum(is.na(p)), ".\n")

##############
# Simulation #
##############
# Register parallel computing
mycore 	<- 3
cl		<- makeCluster(mycore, type = "FORK")
registerDoParallel(cl)
cat("Register", mycore, "core parallel computing. \n")

# Run simulation #
pct		<- proc.time()
sim	<- foreach(i = 1:numsim) %dopar%{
	SimWrapper_fn(omega_deriv, ln_inc = log(sim.unq$income), lambda = lambda, param_est = shr.par, base = beta0_base, 
						X_list = X, price = p, eps_draw = matrix(eps_draw[i,], nrow = 1), method = exp.method, ret.sim = TRUE, 
						par.draw = par_draw, X.inc.inter = FALSE)$Allocation[1,,]
}
sim.arr	<- array(NA, c(numsim, nrow(sim.unq), R), 
					dimnames = list(iter = 1:numsim, unq.idx = 1:nrow(sim.unq), retailer = fmt_name))
for(i in 1:numsim){
	sim.arr[i,,]	<- sim[[i]]
}
# Take the average over random draws
sim.avg	<- apply(sim.arr, c(2, 3), mean, na.rm = T)
use.time	<- proc.time() - pct
cat("The simulation finishes wtih", use.time[3]/60, "min.\n")

# ------------------------------------------------------------ #
# Scenario 0: retail attributes and income at the actual level # 
# This serves as the baseline
sim0	<- merge(sim.data, sim.unq, by = c("scantrack_market_descr", "year", "income"))
cat("Check dimensions:\n dim(sim.data) =", dim(sim.data), ", sim(sim0) =", dim(sim0), "\n")
tmp		<- sim.avg[sim0$idx,]
rownames(tmp)	<- NULL
sim0	<- cbind(sim0, tmp)

# ---------------------------------------------------------------------------------------- #
# Scenario 1: retail attributes at the actual level but fix income at the first year level #
sim1	<- merge(sim.data1, sim.unq, by = c("scantrack_market_descr", "year", "income"))
cat("Check dimensions:\n dim(sim.data1) =", dim(sim.data1), ", sim(sim1) =", dim(sim1), "\n")
tmp		<- sim.avg[sim1$idx,]
rownames(tmp)	<- NULL
sim1	<- cbind(sim1, tmp)

# ------------------------------------------------------------------------------ #
# Scenarior 2: retail attributes fixed at year 1 level but income vary over year #
sim2	<- merge(sim.data2, sim.unq, by = c("scantrack_market_descr", "year", "income"))
cat("Check dimensions:\n dim(sim.data2) =", dim(sim.data2), ", sim(sim2) =", dim(sim2), "\n")
tmp		<- sim.avg[sim2$idx,]
rownames(tmp)	<- NULL
sim2	<- cbind(sim2, tmp)
# Revert the year variable 
sim2	<- sim2[, year:= year0]
sim2	<- sim2[,!"year0", with = FALSE]

stopCluster(cl)
cat("Stop clustering. \n")

################################
# Compare simulation scenarios # 
################################
ggtmp	<- rbind(data.frame(sim0, scenario = "Baseline"), data.frame(sim1, scenario = "FixIncome"), data.frame(sim2, scenario = "FixAttr"))
sel		<- gsub("\\s", ".", fmt_name)
ggtmp[,sel]	<- ggtmp[,sel]/rowSums(ggtmp[,sel])
ggtmp	<- melt(ggtmp[,setdiff(names(ggtmp), "idx")], id.vars = c("scenario", "household_code", "year", "scantrack_market_descr", "income")) 
ggtmp	<- data.table(ggtmp)
ggtmp	<- ggtmp[, list(value = mean(value)), by = list(scenario, year, scantrack_market_descr, variable)]

pdf(paste(plot.wd,"/estrun_",run_id, "/graph_", fname, ".pdf", sep=""), width = ww, height = ww)
ggplot(ggtmp, aes(year, value, col = scenario)) + 
		geom_smooth(method = "gam", formular = y ~ s(x, bs = "tp")) + 
		facet_wrap(~variable, scales = "free")
dev.off()

rm(list = intersect(ls(), c("Allocation_constr_fn", "Allocation_fn", "Allocation_nlop_fn", "ar", "args", "cheb.1d.basis", "cheb.basis", 
		"chebfun", "cl", "dT", "exp_fn", "expFOC_fn", "ggtmp", "i","incl_value_fn", "loadf", "make_plot", "mycore", "myfix", "mysplfun", 
		"mytrimfun", "param_assignR", "pct", "use.time", "s1_index", "sel","sim.avg", "sim.omega", "simExp_fn", "SimOmega_fn", "SimWrapper_fn", 
		"sol", "sol.top", "sol.top2", "solveExp_fn", "spl2dfun", "tmp", "tmp_price", "trim.alpha", "uP_fn", "uPGrad_fn", "ver.date","W", "y.nodes", "numnodes")))

save.image(file = paste("estrun_",run_id, "/" ,fname, ".rdata", sep=""))

cat("This program is done.\n")

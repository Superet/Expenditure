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

# seg_id	<- as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))

model_name 		<- "MDCEV_share_gamma"
run_id			<- 1
make_plot		<- TRUE
interp.method	<- "spline"			# "cheb"
trim.alpha		<- 0 #0.05
cpi.adj			<- TRUE
week.price		<- FALSE

# setwd("~/Documents/Research/Store switching/Processed_data/processed data_20160809")
# plot.wd	<- '~/Desktop'
# sourceCpp(paste("../../Exercise/Multiple_discrete_continuous_model/", model_name, ".cpp", sep=""))
# source("../../Exercise/Multiple_discrete_continuous_model/0_Efficient_Allocation_function.R")
# source("../../Exercise/main/ctrfact_sim_functions.r")

# setwd("/home/brgordon/ccv103/Exercise/run")
# setwd("/kellogg/users/marketing/2661703/Exercise/run")
# setwd("/sscc/home/c/ccv103/Exercise/run")
setwd("U:/Users/ccv103/Documents/Research/Store switching/run")

sourceCpp(paste(model_name, ".cpp", sep=""))
source("0_Efficient_Allocation_function.R")
source("ctrfact_sim_functions.r")

(fname	<- paste("estrun_",run_id,"/MDCEV_est_", Sys.Date(),".rdata",sep=""))

#################################
# Read data and subsetting data #
#################################
load("hh_month_exp.rdata")
load("hh_fmt.rdata")

# Extract 5% random sample
# length(unique(hh_exp$household_code))
# sel			<- sample(unique(hh_exp$household_code), .005*length(unique(hh_exp$household_code)) )
# hh_exp_save	<- hh_exp
# hh_exp		<- subset(hh_exp, household_code %in% sel)

# Retail formats 
fmt_name 	<- as.character(sort(unique(fmt_attr$channel_type)))
R			<- length(fmt_name)

# Compute expenditure share
sel			<- paste("DOL_", gsub("\\s", "_", fmt_name), sep="")
sel1		<- gsub("DOL", "SHR", sel)
for(i in 1:length(fmt_name)){
	hh_exp[,sel1[i]] <- hh_exp[,sel[i]]/hh_exp$dol
}
if(cpi.adj){
	hh_exp$income_midvalue	<- hh_exp$income_midvalue/hh_exp$cpi
}
hh_exp$ln_inc<- log(hh_exp$income_midvalue)
hh_exp$rec		<- 1*(hh_exp$year >= 2008)
cat("Table of segments in the expenditure data:\n"); print(table(hh_exp$first_incomeg)); cat("\n")

# Fill in the missing distance for the channels that are not available for households; 
tmp			<- data.table(hh_dist)
tmp			<- tmp[, list(nozip = all(is.na(panelist_zip_code))), by = list(household_code, year)]
tmp			<- tmp[nozip==TRUE,]
unique(tmp$household) %in% hh_exp$household_code

max.dist	<- 100
sum(hh_dist$distance_haver >max.dist|is.na(hh_dist$distance_haver), na.rm=T)/nrow(hh_dist)
sel			<- hh_dist$distance_haver >max.dist|is.na(hh_dist$distance_haver)
hh_dist[sel,"distance_haver"]	<- max.dist

# Subset data 
tmp		<- data.table(hh_dist)
tmp		<- tmp[,list(nchannel = length(channel_type)), by = list(household_code, year)]

selcol		<- c("household_code", "dol", "year", "month", "rec", "income_midvalue", "ln_inc","first_incomeg", "scantrack_market_descr",paste("SHR_", gsub("\\s", "_", fmt_name), sep=""))
mydata		<- subset(hh_exp[,selcol], dol > .1 )
sel			<- tmp$nchannel == R
sel			<- paste(mydata$household_code, mydata$year, sep="*") %in% paste(tmp[sel,household_code], tmp[sel,year], sep="*")
mydata		<- mydata[sel,]

# sel			<- duplicated(mydata[,c("household_code", "year")])
# mydata		<- mydata[!sel,]

################################
# Organize data for estimation #
################################
# The data required for estimation: 
# outcome variables: y (expenditure), shr (expenditure share);
# explanatory variables: price, X_list

# Match price
tmp		<- dcast(hh_dist[,c("household_code", "year", "channel_type","unitprice_paid")], 
				household_code + year ~ channel_type, value.var = "unitprice_paid")
names(tmp)	<- c("household_code", "year", paste("PRC_", gsub("\\s", "_", fmt_name), sep=""))
mydata	<- merge(mydata, tmp, by = c("household_code", "year"), all.x=T)
mydata	<- merge(mydata, pan[,c("household_code", "panel_year", "household_size", "age")], 
				by.x = c("household_code", "year"), by.y = c("household_code", "panel_year"),all.x=T)
ord		<- order(mydata$household_code, mydata$year, mydata$month)
mydata	<- mydata[ord,]

# Outcome variables as matrix
sel		<- paste("SHR_", gsub("\\s", "_", fmt_name), sep="")
shr		<- as.matrix(mydata[,sel])
y		<- as.vector(mydata$dol)
ln_inc	<- mydata$ln_inc
sel		<- paste("PRC_", gsub("\\s", "_", fmt_name), sep="")
price	<- as.matrix(mydata[,sel])

# Select the index that has the positive expenditure
tmp 	<- price * 1 *(shr> 0)
s1_index<- apply(tmp, 1, which.max)

# Match retailers' attributes
demo.col<- c("ln_inc", "household_size", "age")
selcol	<- c("distance_haver","size_index", "ln_upc_per_mod", "ln_num_module", "prvt_overall")
hh_dist$ln_num_module	<- log(hh_dist$num_module)
hh_dist$ln_upc_per_mod 	<- log(hh_dist$num_upc_per_mod)
beta0_base 	<- which(fmt_name == "Grocery")

Z	<- mydata[,demo.col]
attr_mat	<- merge(mydata[,c("household_code", "year", "month")], hh_dist[,c("household_code", "year", "channel_type", selcol)], 
				by = c("household_code", "year"), all.x=T)
attr_mat	<- attr_mat[order(attr_mat$household_code, attr_mat$year, attr_mat$month),]	
cat("dim(mydata) =", dim(mydata), ".\n")
cat("dim(attr_mat) =", dim(attr_mat), ", dim(mydata)*R =",nrow(mydata)*R, ".\n")
cat("corr(Z):\n"); print(cor(Z)); cat("\n")
cat("corr(attr_mat):\n"); print(cor(attr_mat[,-(1:4)])); cat("\n")

X2	<- X1 	<- vector("list", length=length(fmt_name))
for(i in 1:length(fmt_name)){
	tmp0<- rep(0, R-1)
	if(i<beta0_base){
		tmp0[i]	<- 1
	}else if(i > beta0_base){
		tmp0[(i-1)]	<- 1
	}
	tmp1<- as.matrix(cbind(1, rec=mydata$rec, Z))	
	tmp2<- kronecker(tmp1, t(tmp0))
	colnames(tmp2)	<- paste(rep(c("Intercept", "Rec", demo.col), each = R-1), fmt_name[-beta0_base], sep=":")
	tmp3<- rep(0, R)
	tmp3[i]	<- 1
	tmp3<- kronecker(tmp1, t(tmp3))
	colnames(tmp3)	<- paste(rep(c("Intercept", "Rec", demo.col), each = R), fmt_name, sep=":")
	sel	<- attr_mat$channel_type == fmt_name[i]	
	tmp	<- as.matrix(attr_mat[sel,selcol])
	X1[[i]]	<- cbind(tmp2, tmp, tmp*mydata$rec)
	X2[[i]]	<- cbind(tmp3, tmp, tmp*mydata$rec)
	colnames(X1[[i]])	<- c(colnames(tmp2), colnames(tmp), paste(colnames(tmp), ":Rec", sep=""))
	colnames(X2[[i]])	<- c(colnames(tmp3), colnames(tmp), paste(colnames(tmp), ":Rec", sep=""))
}
nx 		<- ncol(X1[[1]])
cat("nx =", nx,"\n")

##############
# Estimation #
##############
MDCEV_wrapper <- function(param){
	MDCEV_ll_fnC(param, nx, shr, y, s1_index, price, X1, X2, beta0_base)$ll
}

# Guess initial values
inits	<- NULL
inits2	<- NULL
pct		<- proc.time()
for(i in 1:R){
	tmpy	<- 1 * (shr[,i]>0)
	tmpx	<- cbind(Rec = mydata$rec, Z, X1[[i]][,((nx-length(selcol)*2+1):nx)])
	tmpfit	<- glm(y ~ ., data = data.frame(y = tmpy, tmpx), family=binomial(link='logit') )
	inits	<- cbind(inits, coef(tmpfit))
	
	tmpfit	<- lm(shr[,i]*y ~ as.matrix(tmpx))
	inits2	<- cbind(inits2, coef(tmpfit))
	print(i)
}
use.time	<- proc.time() - pct
cat("Initial logit finishes with", use.time[3]/60, "min.\n")

# Normalize the coefficients for demographics relative to those with grocery 
beta.init	<- inits[(1:(length(demo.col)+2)),] - inits[(1:(length(demo.col)+2)),beta0_base]
beta.init	<- c(t(beta.init[,-beta0_base]))
beta.init	<- c(beta.init, rowMeans(inits[-(1:(length(demo.col)+2)),]))
gamma.init	<- c(t(inits2[(1:(length(demo.col)+2)),]))
gamma.init	<- c(gamma.init, rowMeans(inits2[-(1:(length(demo.col)+2)),]))
theta.init	<- c(beta.init*.4, gamma.init/20, -.5)
names(theta.init)	<- c(colnames(X1[[1]]), paste("gamma0_", colnames(X2[[1]]), sep=""), "ln_sigma")
cat("Initial values are:\n"); print(theta.init); cat("\n")
system.time(tmp <- MDCEV_wrapper(theta.init) )

max.method	<- "BFGS"	#"NM"
pct <- proc.time()
sol		<- maxLik(MDCEV_wrapper, start=theta.init, method=max.method)
use.time <- proc.time() - pct
cat("MDCEV estimation finishes with", use.time[3]/60, "min.\n")
print(summary(sol))
cat("--------------------------------------------------------\n")

rm(list= intersect(ls(), c("Allocation_constr_fn","Allocation_fn","Allocation_nlop_fn","incl_value_fn","i","MDCEV_ll_fnC",
		  "MDCEV_LogLike_fnC","MDCEV_wrapper","tmp","tmp1","tmp2","tmp_sol","sel","sel1","sel2","param_assign","tmpX_list", 
		  "use.time", "pct", "uP_fn","uPGrad_fn", "theta_init", "make_plot", "ord","panelist","tmpidx","tmpn","cl", 
		  "tmp_sol1", "GMM_fn", "M_fn", "init", "m", "mycore", "param_assignR", "ggtmp", "ggtmp1", "tmpdat",
		  "interp.method", "trim.alpha", "mysplfun", "mytrimfun", "expFOC_fn", "exp_fn", "solveExp_fn","dT","i", "j", "k",  
		  "simExp_fn", "SimWrapper_fn", "SimOmega_fn", "cheb.1d.basis", "cheb.basis", "chebfun", "omega_parallel", 
		  "lastFuncGrad", "lastFuncParam", "args", "plot.wd", "tmpv1", "tmpv2", "beta.init", "tmpfit")))

save.image(file = fname)

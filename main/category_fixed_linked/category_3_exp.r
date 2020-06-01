
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
seg_id	<- as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))
cat("seg_id =",seg_id, ".\n")

# setwd("~/Documents/Research/Store switching/Processed_data")
# plot.wd	<- '~/Desktop'
# sourceCpp("~/Documents/Research/Store switching/Exercise/main/category_fixed_linked/category_MDCEV.cpp")
# source("../Exercise/Multiple_discrete_continuous_model/0_Efficient_Allocation_function.R")
# source("../Exercise/main/ctrfact_sim_functions.r")

# setwd("/home/brgordon/ccv103/Exercise/run")
# setwd("/kellogg/users/marketing/2661703/Exercise/run")
setwd("/sscc/home/c/ccv103/Exercise/run")
# setwd("U:/Users/ccv103/Documents/Research/Store switching/run")

# Read estimation of allocation model 
run_id		<- 5
# dpt.name	<- c("DG","HBC")
dpt.name	<- c("DG", "GM","NFG", "RF", "HBC","OTHER")
omega.ls	<- vector("list", length = length(dpt.name))						
tmpls		<- ls()

for(ii in 1:length(dpt.name)){
	fname	<- grep(paste(dpt.name[ii], "_seg", seg_id, sep=""), list.files(paste("estrun_", run_id, sep="")), value = TRUE)
	fname	<- grep(".rdata", fname, value = TRUE)
	print(fname)
	load(paste("estrun_", run_id, "/", fname,sep=""))
	omega.ls[[ii]]	<- omega_deriv
	rm(list = setdiff(ls(), c(tmpls, "ii", "tmpls")))
}

load("hh_month_exp_2016-09-26.rdata")
avg.cons	<- read.csv("avg_consump.csv", header = T)
sourceCpp("category_MDCEV.cpp")
(fname	<- paste("estrun_", run_id, "/fixb_exp_seg", seg_id, "_", Sys.Date(), sep=""))

# Extract 5% random sample
# length(unique(hh_month$household_code))
# sel			<- sample(unique(hh_month$household_code), .05*length(unique(hh_month$household_code)) )
# hh_month_save	<- hh_month
# hh_month		<- subset(hh_month, household_code %in% sel)

#################
# Organize data # 
#################
(fmt_name 	<- as.character(sort(unique(hh_dist$channel_type))))
R			<- length(fmt_name)
dpt.lab		<- setNames(c("DRY GROCERY","GENERAL MERCHANDISE", "NON-FOOD GROCERY", "PRODUCE OTHER", "REFRIGERATED FROZEN", "HEALTH BEAUTY CARE"), 
						c("DG", "GM","NFG", "OTHER","RF", "HBC"))
dpt.lab	<- dpt.lab[dpt.name]
nc		<- length(dpt.lab)

# Subset data 
hh_month$rec		<- 1*(hh_month$year >= 2008)
cat("Table of segments in the expenditure data:\n"); print(table(pan$segment)); cat("\n")

# Fill in the missing distance for the channels that are not available for households; 
tmp			<- data.table(hh_dist)
tmp			<- tmp[, list(nozip = all(is.na(panelist_zip_code))), by = list(household_code, year)]
tmp			<- tmp[nozip==TRUE,]

max.dist	<- 100
sum(hh_dist$distance_haver >max.dist|is.na(hh_dist$distance_haver), na.rm=T)/nrow(hh_dist)
sel			<- hh_dist$distance_haver >max.dist|is.na(hh_dist$distance_haver)
hh_dist[sel,"distance_haver"]	<- max.dist

# ------------#
# Subset data #
# Subset the segment 
tmp1		<- subset(pan, ceiling(as.numeric(segment)/3) == seg_id)
expdata		<- subset(hh_month, household_code %in% tmp1$household_code)					
cat("Before cleaning, dim(expdata) =", dim(expdata), "\n")

# Drop the years for which income is not matched
sel			<- with(expdata, paste(household_code, year, sep="*")) %in%
				with(pan, paste(household_code, panel_year, sep="*"))				
sum(sel)/length(sel)
expdata		<- expdata[sel,]
expdata		<- merge(expdata, pan[,c("household_code", "panel_year", "income_mid","scantrack_market_descr")], 
				by.x = c("household_code", "year"), by.y = c("household_code", "panel_year"))
cat("dim(expdata) =", dim(expdata), "\n")
								
# Filter the households who did not have all the channel
tmp			<- data.table(hh_dist)
tmp			<- tmp[,list(nchannel = length(channel_type)), by = list(household_code, year)]
sel			<- tmp$nchannel == R
sel			<- paste(expdata$household_code, expdata$year, sep="*") %in% paste(tmp[sel,household_code], tmp[sel,year], sep="*")	
expdata		<- expdata[sel,]
cat("After cleanning distance, we have", mean(sel), "observations.\n")

tmp 		<- data.table(fmt_dpt)
tmp 		<- tmp[,list(nchannel = length(channel_type)), by = list(department, scantrack_market_descr, year)]			
tmp			<- subset(tmp, year > 2003 & year < 2011)
sel			<- tmp$nchannel < R
sel			<- !paste(expdata$scantrack_market_descr, expdata$year, sep="*") %in% 
			   paste(tmp[sel,scantrack_market_descr], tmp[sel,year], sep="*")	
expdata		<- expdata[sel,]	
cat("After filtering out the households who do not have complete measurement from all channels, we have", mean(sel), "observations.\n")		

# Only keep relevant columns 
tmp		<- matrix(NA, nrow(expdata), length(dpt.name), dimnames = list(NULL, dpt.name))
for(i in 1:length(dpt.name)){
	sel	<- grep(dpt.name[i], names(expdata), value=TRUE)
	tmp[,i]	<- rowSums(as.matrix(expdata[,sel]))
}
expdata	<- cbind(expdata[,c("household_code", "year", "month", "rec", "scantrack_market_descr","income_mid")], tmp)
expdata	<- merge(expdata, avg.cons, by = "household_code", all.x = T)				
cat("dim(expdata) =", dim(expdata), "\n")

# Drop any observations where income is smaller than total expenditure 
tmp		<- expdata$income_mid/12 - rowSums(expdata[,dpt.name])
sel		<- tmp > 0
expdata	<- expdata[sel,]

# Calculate stock variable = last expenditure - average consumption rate
expdata	<- data.table(expdata)
expdata	<- expdata[,':='(s_DG 	= c(NA, DG[-length(DG)]) - AVG_DG, 		s_GM 	= c(NA, GM[-length(GM)]) - AVG_GM, 
 						 s_NFG 	= c(NA, NFG[-length(NFG)]) - AVG_NFG, s_OTHER 	= c(NA, OTHER[-length(OTHER)]) - AVG_OTHER, 
						 s_RF 	= c(NA, RF[-length(RF)]) - AVG_RF, 		s_HBC 	= c(NA, HBC[-length(HBC)]) - AVG_HBC), 
					by = list(household_code)]
expdata	<- data.frame(expdata)
expdata	<- expdata[, setdiff(names(expdata), grep("AVG",names(expdata)) )]
expdata	<- subset(expdata, !is.na(s_DG))
cat("dim(expdata) =", dim(expdata), "\n")

# ----------------------------------------------------------------# 
# Model matrix
# First and second derivative of purchase utility, stock variable # 
ord		<- order(expdata$household_code, expdata$year, expdata$month)
d1		<- d2 <- matrix(NA, nrow(expdata), nc)
X		<- vector("list", length = length(dpt.name))

for(i in 1:length(dpt.name)){
	# Price matrix 
	sel		<- fmt_dpt$department == dpt.lab[i]
	tmp		<- dcast(fmt_dpt[sel,c("scantrack_market_descr", "year", "channel_type","unitprice_paid")], 
					scantrack_market_descr + year ~ channel_type, value.var = "unitprice_paid")
	names(tmp)	<- c("scantrack_market_descr", "year", paste("PRC_", gsub("\\s", "_", fmt_name), sep=""))
	price	<- merge(expdata[,c("scantrack_market_descr", "year")], tmp, by = c("scantrack_market_descr", "year"), all.x=T)
	price	<- as.matrix(price[ord,paste("PRC_", gsub("\\s", "_", fmt_name), sep="")])
	
	# Fit first and second derivative 
	tmpdat	<- data.frame(y = expdata[, dpt.name[i]])
	tmpdat$price	<- price
	d1[,i]	<- omega.ls[[i]](tmpdat)
	d2[,i]	<- omega.ls[[i]](tmpdat, deriv = 1)
	
	# Model matrix 
	tmp0	<- rep(0, nc)
	tmp0[i]	<- 1
	tmp1	<- as.matrix(cbind(1, rec=expdata$rec))	
	tmp2	<- kronecker(tmp1, t(tmp0))
	colnames(tmp2)	<- paste(rep(c("Intercept", "Rec"), each = nc), dpt.name, sep=":")
	X[[i]]	<- as.matrix(cbind(tmp2, s = expdata[,paste("s_",dpt.name[i], sep="")]))
}
rm(list = c("price", "tmp0", "tmp1", "tmp2", "tmpdat", "tmpls"))
ncol(X[[1]])
cat("Summary of d1:\n"); print(summary(d1)); cat(".\n")
cat(mean(d1 <=0), " d1 <=0\n")
sel	<- d1 <= 0				# d1 must be positive
d1[sel]	<- 1e-5

e	<- as.matrix(expdata[,dpt.name])
y0 	<- expdata[,"income_mid"]/12 - rowSums(e)
e_sgn	<- 1 * (e > 0)
M	<- rowSums(e_sgn)	
cat("Summary of e:\n"); print(summary(cbind(e, y0))); cat("\n")
lnJ	<- -2*log(y0) + log(rowSums(-e_sgn*d1/d2)) + rowSums(e_sgn * log(-d2/d1))
summary(lnJ)

##############
# Estimation #
##############
ll_wrapper	<- function(param){
	MDCEV_ll_fnC(param, e = e, y0 = y0, X = X, M = M, d1 = d1, d2 = d2, lnJ = lnJ)
}

# Take initial values 
fit.lm	<- vector("list", length = length(dpt.name))
for(i in 1:nc){
	tmpy	<- -log(d1[,i]) - log(y0)
	tmpx	<- X[[i]][,c(grep(dpt.name[i], colnames(X[[1]]), value = TRUE),"s")]
	fit.lm[[i]]	<- lm(tmpy ~ tmpx - 1)
}
theta.init	<- sapply(fit.lm, coef)
theta.init	<- c(t(theta.init[-nrow(theta.init),]), mean(theta.init[nrow(theta.init),]))
names(theta.init)	<- colnames(X[[1]])
system.time(tmp <- ll_wrapper(theta.init))

max.method	<- "BFGS"	#"NM"
pct <- proc.time()
sol		<- maxLik(ll_wrapper, start=theta.init, method=max.method)
use.time <- proc.time() - pct
cat("MDCEV estimation finishes with", use.time[3]/60, "min.\n")
print(summary(sol))
cat("--------------------------------------------------------\n")

####################
# Save the results #
####################
rm(list= intersect(ls(), c("hh_month", "i", "ii", "ll_wrapper", "max.dist", "MDCEV_ll_fnC", "ord", "tmpx", "tmpy", "sel", "sel1")))

save.image(file = paste(fname, ".rdata", sep=""))

cat("This program is done. ")


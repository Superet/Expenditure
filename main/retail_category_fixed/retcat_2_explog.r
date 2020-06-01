
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
# seg_id 	<- 1
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
run_id		<- 6
fmt_name	<- c("Convenience Store","Discount Store","Dollar Store","Drug Store", "Grocery", "Warehouse Club")
omega.ls	<- vector("list", length = length(fmt_name))						
tmpls		<- ls()

for(ii in 1:length(fmt_name)){
	fname	<- grep(paste(fmt_name[ii], "_seg", seg_id, sep=""), list.files(paste("estrun_", run_id, sep="")), value = TRUE)
	fname	<- grep(".rdata", fname, value = TRUE)
	print(fname)
	load(paste("estrun_", run_id, "/", fname,sep=""))
	
	# Subset simulation data
	sel			<- duplicated(mydata[,c("household_code", "year")]) & mydata$dol <= quantile(mydata$dol, .95)
	sht.X1		<- lapply(X1, function(x) x[!sel,])
	sht.price	<- price[!sel,]
	sht.y		<- y[!sel]
	
	if(trim.alpha == 0){
		tmpdat	<- data.frame(omega = colMeans(omega_draw), y = sht.y, ln_inc = ln_inc[!sel])
	}else{
		tmpdat	<- data.frame(omega = apply(omega_draw, 2, mean, trim = trim.alpha, na.rm=T), y = sht.y, ln_inc = ln_inc[!sel])
	}
	# Scale omega such that for the maximum expenditure y, omega(y) = log(y+1)
	sel1			<- which.max(tmpdat$y)
	myscale			<- log(tmpdat[sel1,"y"]+1)/tmpdat[sel1,"omega"]
	cat("Omega is scaled by", myscale, ".\n")
	tmpdat$omega	<- tmpdat$omega * myscale	
	tmpdat$price	<- sht.price[,-beta0.base]
	tmpdat$X		<- as.matrix(do.call(cbind, lapply(sht.X1[1:(nc-1)], function(x) x[,selcol])))	
	tmpdat$rec		<- mydata[!sel,"rec"]
	tmpdat$dv		<- log(tmpdat$y) - tmpdat$omega
	omega.fit		<- lm(dv ~ rec + price + X, data = tmpdat)
	cat("Summary of the fitted log(gamma):\n"); print(summary(fitted(omega.fit))); cat("\n")
	
	omega.ls[[ii]]	<- omega.fit
	rm(list = setdiff(ls(), c(tmpls, "ii", "tmpls")))
}
gamma.par	<- sapply(omega.ls, coef)

load("hh_month_exp_2016-09-26.rdata")
avg.cons	<- read.csv("avg_consump.csv", header = T)
sourceCpp("category_MDCEV.cpp")
(fname	<- paste("estrun_", run_id, "/retcat_explog_seg", seg_id, "_", Sys.Date(), sep=""))

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
dpt.name	<- c("DG", "GM","NFG","RF", "HBC", "OTHER")
dpt.lab		<- c("DRY GROCERY","GENERAL MERCHANDISE", "NON-FOOD GROCERY", "REFRIGERATED FROZEN", "HEALTH BEAUTY CARE","PRODUCE OTHER")
nc			<- length(dpt.name)

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

tmp 		<- data.table(subset(fmt_dpt, year > 2003 & year < 2011))
tmp 		<- tmp[,list(n = length(department)), by = list(scantrack_market_descr, channel_type, year)]			
sel			<- tmp$n < nc
sel			<- !paste(expdata$scantrack_market_descr, expdata$year, sep="*") %in% 
			   paste(tmp[sel,scantrack_market_descr], tmp[sel,year], sep="*")	
expdata		<- expdata[sel,]	
cat("After filtering out the households who do not have complete measurement from all channels, we have", mean(sel), "% observations.\n")		

# Only keep relevant columns 
expdata	<- expdata[,c("household_code", "year", "month", "rec", "scantrack_market_descr","income_mid",
						grep("DOL_", names(expdata), value=TRUE))]
# expdata	<- merge(expdata, avg.cons, by = "household_code", all.x = T)				
cat("dim(expdata) =", dim(expdata), "\n")

# Drop any observations where income is smaller than total expenditure 
tmp		<- expdata$income_mid/12 - rowSums(expdata[,grep("DOL_", names(expdata), value=TRUE)])
sel		<- tmp > 0
mean(sel)
expdata	<- expdata[sel,]

# ----------------------------------------------------------------# 
# Model matrix
# First and second derivative of purchase utility, stock variable # 
ord		<- order(expdata$household_code, expdata$year, expdata$month)
gamma	<- matrix(NA, nrow(expdata), R)
X		<- vector("list", length = length(fmt_name))
selcol	<- c("size_index", "ln_upc_per_mod", "ln_num_module", "prvt_overall")
fmt_dpt$ln_num_module	<- log(fmt_dpt$num_module)
fmt_dpt$ln_upc_per_mod 	<- log(fmt_dpt$num_upc_per_mod)
beta0.base	<- which(dpt.name == "OTHER")

# Normalize distance by zipcode
dist	<- data.table(hh_dist[,c("household_code", "panelist_zip_code","year", "channel_type", "distance_haver")])
dist	<- dist[, maxd:=max(distance_haver, na.rm=T), by = list(panelist_zip_code)]
dist	<- dist[, distance := distance_haver/maxd]
dist	<- data.frame(dist)
dist	<- merge(expdata[,c("household_code","scantrack_market_descr", "year","month")], dist, by = c("household_code", "year"), all.x=T)
dist	<- dist[order(dist$household_code, dist$year, dist$month, dist$channel_type),]
cat("dim(dist) =", dim(dist), ", dim(expdata)*R =",nrow(expdata)*R, ".\n")

attr_mat	<- merge(expdata[,c("household_code","scantrack_market_descr", "year","month")], 
				fmt_dpt[,c("scantrack_market_descr", "year", "channel_type","department",selcol)], 
				by = c("scantrack_market_descr", "year"), all.x=T)
attr_mat	<- attr_mat[order(attr_mat$household_code, attr_mat$year, attr_mat$month),]	
cat("dim(attr_mat) =", dim(attr_mat), ", dim(expdata)*nc*R =",nrow(expdata)*nc*R, ".\n")

for(i in 1:length(fmt_name)){
	# Price matrix 
	sel		<- fmt_dpt$channel_type == fmt_name[i]
	tmp		<- fmt_dpt[sel,c("scantrack_market_descr", "year", "department","unitprice_paid")]
	tmp$department	<- factor(tmp$department, levels = dpt.lab, labels = dpt.name)
	tmp		<- dcast(tmp, scantrack_market_descr + year ~ department, value.var = "unitprice_paid")
	names(tmp)	<- c("scantrack_market_descr", "year", paste("PRC_", dpt.name, sep=""))
	price	<- merge(expdata[,c("household_code","scantrack_market_descr", "year","month")], tmp, by = c("scantrack_market_descr", "year"), all.x=T)
	price	<- as.matrix(price[ord,paste("PRC_",dpt.name, sep="")])
	
	# Retail attributes
	sel		<- attr_mat$channel_type == fmt_name[i] 
	tmp		<- do.call(cbind,lapply(dpt.lab[-beta0.base], function(x) attr_mat[sel&attr_mat$department==x,selcol]))
	
	# Fit first and second derivative 
	tmpdat	<- data.frame(y = expdata[,paste("DOL_", gsub("\\s", "_", fmt_name[i]), sep="")])
	tmpdat$price	<- price[,-beta0.base]
	tmpdat$rec		<- expdata$rec
	tmpdat$X		<- as.matrix(tmp)
	gamma[,i]		<- exp(predict(omega.ls[[i]], tmpdat))
	
	# Model matrix 
	tmp0	<- rep(0, R)
	tmp0[i]	<- 1
	tmp1	<- as.matrix(cbind(1, rec=expdata$rec))	
	tmp2	<- kronecker(tmp1, t(tmp0))
	colnames(tmp2)	<- paste(rep(c("Intercept", "Rec"), each = nc), fmt_name, sep=":")
	sel		<- dist$channel_type == fmt_name[i]
	X[[i]]	<- as.matrix(cbind(tmp2, dist[sel,"distance"], dist[sel,"distance"]*expdata$rec))
	colnames(X[[i]])	<- c(colnames(tmp2), "distance", "distance:Rec")
}
rm(list = c("price", "tmp0", "tmp1", "tmp2", "tmpdat", "tmpls","attr_mat", "tmp","dist"))
ncol(X[[1]])

e	<- as.matrix(expdata[,paste("DOL_", gsub("\\s", "_", fmt_name), sep="")])
y0 	<- expdata[,"income_mid"]/12 - rowSums(e)
e_sgn	<- 1 * (e > 0)
M	<- rowSums(e_sgn)	
cat("Summary of e:\n"); print(summary(cbind(e, y0))); cat("\n")
lnJ	<- -log(y0+1) + log(rowSums(e+gamma)) - rowSums(e+gamma)
summary(lnJ)
cat("Summary of gamma:\n"); print(summary(gamma)); cat("\n")

# Plot the conditional expected utility function 
tmp 	<- seq(0, 200, 5)
tmp1	<- colMeans(gamma)
ggtmp	<- data.frame(y = tmp, base = log(tmp + 1), sapply(tmp1, function(x) log(tmp/x + 1)))
colnames(ggtmp)	<- c("y", "base", fmt_name)
ggtmp	<- melt(ggtmp, id.vars = "y")
ggtmp1	<- data.frame(y = tmp, base = log(tmp + 1), sapply(tmp1, function(x) x*log(tmp/x + 1)))
colnames(ggtmp1)	<- c("y", "base", fmt_name)
ggtmp1	<- melt(ggtmp1, id.vars = "y")

# ggplot(ggtmp, aes(y, value)) + geom_line() + facet_wrap(~variable, scales = "free")
# ggplot(ggtmp1, aes(y, value)) + geom_line() + facet_wrap(~variable, scales = "free")

##############
# Estimation #
##############
ll_wrapper	<- function(param){
	MDCEV_ll_fnC(param, e = e, e_sgn = e_sgn, y0 = y0, X = X, M = M, gamma = gamma, lnJ = lnJ)
}

# Take initial values 
fit.lm	<- vector("list", length = length(fmt_name))
for(i in 1:R){
	# tmpy	<- log(e[,i]/gamma[,i] + 1) - log(y0+1)
	tmpy	<- log(e[,i] + gamma[,i]) - log(y0+1)
	tmpx	<- X[[i]][,c(grep(fmt_name[i], colnames(X[[1]]), value = TRUE),"distance","distance:Rec")]
	fit.lm[[i]]	<- lm(tmpy ~ tmpx - 1)
}
theta.init	<- sapply(fit.lm, coef)
# theta.init[1,]	<- theta.init[1,] - c(10, 1, 4, 4,0, 3)
theta.init	<- c(t(theta.init[1:2,]), rowMeans(theta.init[3:4,]))
names(theta.init)	<- colnames(X[[1]])
cat("Initial values are:\n", theta.init, ".\n")
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
rm(list= intersect(ls(), c("hh_month", "i", "ii", "ll_wrapper", "max.dist", "MDCEV_ll_fnC", "ord", 
							"tmpx", "tmpy", "sel", "sel1","tmp", "pct", "use.time","args", "lnJ", "M")))

save.image(file = paste(fname, ".rdata", sep=""))

cat("This program is done. ")


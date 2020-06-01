options(error = quote({dump.frames(to.file = TRUE)}))
options(mc.cores = parallel::detectCores())

library(ggplot2)
library(reshape2)
library(data.table)
library(lme4)
library(foreach)
library(doParallel)
library(MCMCglmm)

# Set working directory
arr_idx		<- as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))
# setwd("~/Documents/Research/Store switching/processed data")
# plot.wd	<- '/home/chaoqunchen/Desktop'
# source("../Exercise/main/outreg function.R")
# setwd("/home/brgordon/ccv103/Exercise/run")
# setwd("/kellogg/users/marketing/2661703/Exercise/run")
setwd("/sscc/home/c/ccv103/Exercise/run")
plot.wd 	<- getwd()
ww			<- 6.5
ww1			<- 12
ar 			<- .6
make_plot 	<- TRUE

# Read data 
load("hh_biweek_exp.rdata")

# Extract subsection data
Nsec		<- 30				# Divide the full data by subsections
subpan		<- ceiling(c(0, c(1:Nsec)/Nsec * length(unique(hh_exp$household_code))))
sel			<- (subpan[arr_idx]+1):subpan[(arr_idx+1)]
hh_exp		<- subset(hh_exp, household_sample %in% sel)
cat("The dimension of sub-data is", dim(hh_exp), "\n")

#############
# Functions # 
#############
my_forward 	<- function(x){
	return(c(x[-1], NA))
}

my_lag		<- function(x){
	return(c(NA, x[-length(x)]))
}

#######################################
# Organize data and create covariates # 
#######################################
# Retail formats 
fmt_name 	<- as.character(sort(unique(fmt_attr$channel_type)))
R			<- length(fmt_name)

# Compute expenditure share
sel			<- paste("DOL_", gsub("\\s", "_", fmt_name), sep="")
sel1		<- gsub("DOL", "SHR", sel)
for(i in 1:length(fmt_name)){
	hh_exp[,sel1[i]] <- hh_exp[,sel[i]]/hh_exp$dol
}

#---------------------------------------#
# Construct regression data 
sel			<- paste("DOL_", gsub("\\s", "_", fmt_name), sep="")
sel1		<- gsub("DOL", "IC", sel)
for(i in 1:length(fmt_name)){
	hh_exp[,sel1[i]] <- 1*(hh_exp[,sel[i]]>0)
}
colMeans(hh_exp[,sel1])

# Segment households based on their initial income level 
panelist	<- data.table(hh_exp)
setkeyv(panelist, c("household_code","year","biweek"))
panelist	<- panelist[,list(income = first_income[1], first_famsize = famsize[1]), by=list(household_code)]
tmp			<- quantile(panelist$income, c(0, .33, .67, 1))
num_grp		<- 3
panelist	<- panelist[, first_incomeg := cut(panelist$income, tmp, labels = paste("T", 1:num_grp, sep=""), include.lowest = T)]
hh_exp		<- merge(hh_exp, data.frame(panelist)[,c("household_code", "first_incomeg")], by = "household_code", all.x=T )
hh_exp$first_incomeg	<- factor(hh_exp$first_incomeg, levels = c("T2", "T1", "T3"))
# hh_exp$first_incomeg	<- cut(hh_exp$first_income, c(0, 17, 23, 30), labels = c("T1", "T2", "T3"))

# Add lagged purchases
hh_exp	<- data.table(hh_exp)
setkeyv(hh_exp, c("household_code", "biweek"))
hh_exp	<- hh_exp[,lag_dol := my_lag(dol), by = list(household_code)]
hh_exp	<- data.frame(hh_exp)
hh_exp	<- hh_exp[!is.na(hh_exp$lag_dol),]			# Drop the first observation for each household
hh_exp$lag_dol 	<- hh_exp$lag_dol/100

# Factorize some variables
hh_exp$ln_income	<- log(hh_exp$income_midvalue/hh_exp$cpi)
hh_exp$year			<- as.factor(hh_exp$year)
hh_exp$recession 	<- factor(hh_exp$recession)
hh_exp$month <- month(as.Date("2004-1-1", format="%Y-%m-%d") + 14*(hh_exp$biweek-1))
hh_exp$month		<- as.factor(hh_exp$month)

# Select only the relevant columns
selcol	<- c("household_code", "biweek", "dol", "year", "month",
			paste("IC_", gsub("\\s", "_", fmt_name), sep=""), "ln_income", "lag_dol", "first_incomeg", "recession")
hh_exp	<- hh_exp[,selcol]

############
# Fit glmm # 
############
# Set prior
# NOTE: prior is p(theta)^(1/Nsec) for each subsection of the full data
# The original prior is G ~ IW(V, nu), now the scaled V ~ IW(V/Nsec, (nu+p+1)/Nsec - p - 1)
myprior	<- list(R = list(V = diag(R)/Nsec, n = (2*R+1)/Nsec - R - 1), 
				G = list(G1= list(V = diag(R)/(R+1), n = R)))		
				
# Make initial values
selcol		<- paste("IC_", gsub("\\s", "_", fmt_name), sep="")
for(i in 1:R){
	fml		<- as.formula(paste(selcol[i], "~ ln_income + lag_dol + month + (1|household_code)", sep=""))
	(tmpfit	<- glmer(fml, family = binomial(link = "logit"), data = hh_exp))
	if(i == 1){
		beta_init 	<- fixef(tmpfit)
		beta0_init	<- as.vector(ranef(tmpfit)$household_code[,1])
		omega_init	<- as.numeric(VarCorr(tmpfit)$household_code)
		z_init		<- as.vector(predict(tmpfit))
	}else{
		beta_init	<- rbind(beta_init, fixef(tmpfit))
		beta0_init	<- rbind(beta0_init, as.vector(ranef(tmpfit)$household_code[,1]))
		omega_init	<- c(omega_init, as.numeric(VarCorr(tmpfit)$household_code))
		z_init		<- cbind(z_init, as.vector(predict(tmpfit)))
	}
}	
mystart	<- list(liab = z_init)			

for(i in 1:R){
	fml		<- as.formula(paste(selcol[i], "~ ln_income*first_incomeg + lag_dol + month + (1|household_code)", sep=""))
	(tmpfit	<- glmer(fml, family = binomial(link = "logit"), data = hh_exp))
	if(i == 1){
		beta_init 	<- fixef(tmpfit)
		beta0_init	<- as.vector(ranef(tmpfit)$household_code[,1])
		omega_init	<- as.numeric(VarCorr(tmpfit)$household_code)
		z_init		<- as.vector(predict(tmpfit))
	}else{
		beta_init	<- rbind(beta_init, fixef(tmpfit))
		beta0_init	<- rbind(beta0_init, as.vector(ranef(tmpfit)$household_code[,1]))
		omega_init	<- c(omega_init, as.numeric(VarCorr(tmpfit)$household_code))
		z_init		<- cbind(z_init, as.vector(predict(tmpfit)))
	}
}
mystart_het	<- list(liab = z_init)			

rm(list = setdiff(ls(), c("hh_exp", "myprior", "arr_idx", "mystart", "mystart_het")))

sim	<- MCMCglmm(cbind(IC_Convenience_Store, IC_Discount_Store, IC_Dollar_Store, IC_Drug_Store, IC_Grocery, IC_Warehouse_Club) ~ 
							trait:ln_income + trait:lag_dol + trait:month + trait -1, 
					random = ~ idh(trait):household_code,
					rcov = ~ corg(trait):units, 
					family = rep("categorical", 6), 
					prior = myprior, 
					data = hh_exp,  
					start = mystart, 
					pr = TRUE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE,  
					nitt = 35000, burnin = 5000, thin= 30)

sim_het	<- MCMCglmm(cbind(IC_Convenience_Store, IC_Discount_Store, IC_Dollar_Store, IC_Drug_Store, IC_Grocery, IC_Warehouse_Club) ~ 
							trait:ln_income*first_incomeg + trait:lag_dol + trait:month + trait -1, 
					random = ~ idh(trait):household_code,
					rcov = ~ corg(trait):units, 
					family = rep("categorical", 6), 
					prior = myprior, 
					data = hh_exp,  
					start = mystart_het, 
					pr = TRUE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE,  
					nitt = 60000, burnin = 10000, thin= 50)					

# Save results 
save(sim, sim_het, file = paste("glmm_sub", arr_idx, ".rdata", sep=""))

cat("This program is done. \n")

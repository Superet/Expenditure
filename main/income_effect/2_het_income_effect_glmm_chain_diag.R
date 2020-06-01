options(error = quote({dump.frames(to.file = TRUE)}))
options(mc.cores = parallel::detectCores())

library(ggplot2)
library(reshape2)
library(data.table)
library(gridExtra)
library(lme4)
library(foreach)
library(doParallel)
library(MCMCglmm)

chain_idx	<- as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))
arr_idx		<- 1
codecase	<- "MCMC"
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
# hh_exp		<- read.csv("2_hh_biweek_exp_merge.csv")
# fmt_attr 	<- read.csv("1_format_year_attr.csv")
# price_dat	<- read.csv("1_format_biweek_price.csv")
load("hh_biweek_exp.rdata")

# Extract subsection data
set.seed(666)
panid		<- sample(unique(hh_exp$household_code))
cat("Number of households =", length(panid), "\n")
Nsec		<- 30				# Divide the full data by subsections
subpan		<- ceiling(c(0, c(1:Nsec)/Nsec * length(panid)))
sel			<- panid[(subpan[arr_idx]+1):subpan[(arr_idx+1)]]
hh_exp		<- subset(hh_exp, household_code %in% sel)
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
mydata 	<- hh_exp
sel			<- paste("DOL_", gsub("\\s", "_", fmt_name), sep="")
sel1		<- gsub("DOL", "IC", sel)
for(i in 1:length(fmt_name)){
	mydata[,sel1[i]] <- 1*(mydata[,sel[i]]>0)
}
colMeans(mydata[,sel1])
mydata$first_incomeg	<- cut(hh_exp$first_income, c(0, 17, 23, 30), labels = c("T1", "T2", "T3"))

# Add lagged purchases
mydata	<- data.table(mydata)
setkeyv(mydata, c("household_code", "biweek"))
mydata	<- mydata[,lag_dol := my_lag(dol), by = list(household_code)]
mydata	<- data.frame(mydata)
mydata	<- mydata[!is.na(mydata$lag_dol),]			# Drop the first observation for each household
mydata$lag_dol 	<- mydata$lag_dol/100

# Factorize some variables
mydata$ln_income	<- log(mydata$income_midvalue/mydata$cpi)
mydata$year			<- as.factor(mydata$year)
mydata$recession 	<- factor(mydata$recession)

# Select only the relevant columns
# Since convenience stores incidence is too low, we drop this variable; 
selcol	<- c("household_code", "biweek", "dol", "year", 
			paste("IC_", gsub("\\s", "_", fmt_name), sep=""), "ln_income", "lag_dol", "first_incomeg", "recession")
mydata	<- mydata[,selcol]

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
	fml		<- as.formula(paste(selcol[i], "~ ln_income*first_incomeg + lag_dol + (1|household_code)", sep=""))
	(tmpfit	<- glmer(fml, family = binomial(link = "logit"), data = mydata))
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
		
rm(list = setdiff(ls(), c("mydata", "myprior", "arr_idx", "codecase","chain_idx", "make_plot", "plot.wd", "mystart")))

if(codecase == "diagnosis"){
	CL		<- makeCluster(2, type = "FORK")
	registerDoParallel(CL)
	# registerDoMC(2)
	# clusterExport(cl = CL, c("mydata", "myprior"))
	
	# Run parallel fitting 
	prt		<- proc.time()
	fitls 	<- foreach(i = 1:4, .packages = "MCMCglmm") %dopar%{
	# fitls 	<- parLapply(CL, 1:2, fun = function(cid) {  
	# require(MCMCglmm)
			MCMCglmm(cbind(IC_Discount_Store, IC_Dollar_Store, IC_Drug_Store, IC_Grocery, IC_Warehouse_Club) ~ 
								trait:ln_income + trait:lag_dol + trait -1, 
						random = ~ idh(trait):household_code,
						rcov = ~ corg(trait):units, 
						family = rep("categorical", 5), 
						prior = myprior, 
						data = mydata, 
						saveX = FALSE, saveZ = FALSE, saveXL = FALSE,  
						nitt = 20000, burnin = 7500, thin= 50)
	}					
	use.time	<- proc.time() - prt
	cat("GLMM sampling takes", use.time[3]/60,"min.\n")
	stopCluster(CL)

	# -------------------------------# 
	# Diagnosis MCMC 
	beta.mcmc	<- do.call("mcmc.list", lapply(fitls, function(x) as.mcmc(x$Sol)))
	vcv.mcmc	<- do.call("mcmc.list", lapply(fitls, function(x) as.mcmc(x$VCV)))
	dev.mcmc  <- do.call("mcmc.list", lapply(fitls, function(x) as.mcmc(x$Deviance)))
	
	if(make_plot){
		pdf(paste(plot.wd, "/graph_trace.pdf", sep=""), width = 8, height = 8)
		plot(beta.mcmc)
		plot(vcv.mcmc)
		plot(dev.mcmc)
		dev.off()
	}
		
	print(gelman.diag(beta.mcmc))
	# print(gelman.diag(vcv.mcmc))

	#---------------------------------# 
	# Display results 
	cat("Posterior distribution of beta:\n")
	print(round(cbind(HPDinterval(mcmc(as.matrix(beta.mcmc))), postior.mode = posterior.mode(mcmc(as.matrix(beta.mcmc))) ), 3) )
	cat("Posterior distribution of covariance matrix:\n")
	print(round(cbind(HPDinterval(mcmc(as.matrix(vcv.mcmc))), postior.mode = posterior.mode(mcmc(as.matrix(vcv.mcmc))) ), 3) )
	
	# Save results 
	save.image(file = paste("glmm_diag_sub", arr_idx, ".rdata", sep=""))
}else{
	sim	<- MCMCglmm(cbind(IC_Convenience_Store, IC_Discount_Store, IC_Dollar_Store, IC_Drug_Store, IC_Grocery, IC_Warehouse_Club) ~ 
								trait:ln_income:first_incomeg + trait:lag_dol + trait -1, 
						random = ~ idh(trait):household_code,
						rcov = ~ corg(trait):units, 
						family = rep("categorical", 6), 
						prior = myprior, 
						data = mydata,  
						start = mystart, 
						pr = TRUE, saveX = FALSE, saveZ = FALSE, saveXL = FALSE,  
						nitt = 30000, burnin = 10000, thin= 80)
	
	# Save results 
	save(sim, file = paste("glmm_sub1_chain", chain_idx, ".rdata", sep=""))
}

cat("This program is done. \n")


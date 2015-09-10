library(ggplot2)
library(reshape2)
library(data.table)
library(gridExtra)
library(scales)
library(r2excel)
library(lme4)
library(parallel)
library(rstan)

options(error = quote({dump.frames(to.file = TRUE)}))
options(mc.cores = parallel::detectCores())

# setwd("~/Documents/Research/Store switching/processed data")
# plot.wd	<- '/home/chaoqunchen/Desktop'
# source("../Exercise/main/outreg function.R")
# setwd("/home/brgordon/ccv103/Exercise/run")
setwd("/kellogg/users/marketing/2661703/Exercise/run")
# setwd("/sscc/home/c/ccv103/Exercise/run")
plot.wd 	<- getwd()
ww			<- 6.5
ww1			<- 12
ar 			<- .6

write2csv	<- TRUE
make_plot	<- TRUE
outfile		<- paste(plot.wd, "/2_income_effect_", gsub("-", "", as.character(Sys.Date())), ".rdata", sep="")
outxls		<- paste(plot.wd, "/2_het_income_effect_", gsub("-", "", as.character(Sys.Date())), ".xlsx", sep="")
mywb		<- createWorkbook()
sht1		<- createSheet(mywb, "Sheet1")
sht2		<- createSheet(mywb, "Regression")
sht3		<- createSheet(mywb, "SUR")

# Read data 
hh_exp		<- read.csv("2_hh_biweek_exp_merge_larger.csv")
fmt_attr 	<- read.csv("1_format_year_attr.csv")
price_dat	<- read.csv("1_format_biweek_price.csv")
codebook	<- read.csv("code_book.csv")
source("outreg function.R")

# Extract 5% random sample
# length(unique(hh_exp$household_code))
# sel			<- sample(unique(hh_exp$household_code), .01*length(unique(hh_exp$household_code)) )
# hh_exp_save	<- hh_exp
# hh_exp		<- subset(hh_exp, household_code %in% sel)

# Stan source code
src_code <- "
data{
	int<lower = 1>	N; 
	int<lower = 1>	R; 
	int<lower = 1>	K; 
	int<lower = 1>	nhh; 
	vector[K]		X[N];
	int 			y[N,R];
	cov_matrix[R]	Sigma_prior;
	int<lower = 1> 	hh_idx[N]; 
}
parameters{
	vector[K]	beta_raw[R];
	vector[nhh]	beta0_raw[R];
	vector<lower = 0>[R]	omega_raw;
	corr_matrix[R] Sigma_raw;
	vector[R]	z[N];
}
transformed parameters{
	vector[K]	beta[R]; 
	vector[nhh]	beta0[R];
	vector[R]	omega;
	corr_matrix[R]	Sigma;
	vector[R]	beta0_mean; 
	
	for(i in 1:R){
		beta0_mean[i]	<- mean(beta0_raw[i]);
		beta[i]		<- beta_raw[i]/sqrt(2);
		beta[i,1]	<- beta[i,1] + beta0_mean[i]/sqrt(2);
		beta0[i]	<- (beta0_raw[i] - beta0_mean[i])/sqrt(2); 
		omega[i]	<- omega_raw[i]/sqrt(2); 
		for(j in 1:R){
			if(i == j){	Sigma[i,j]	<- 1; }
			else{ Sigma[i,j]	<- Sigma_raw[i,j]/sqrt(2); }
		}
	}
}
model{
	matrix[N,R] Mu; 
	
	# Prior distribution
	for(i in 1:R){
		beta_raw[i] ~ normal(0, 10);
	}
	omega_raw ~ uniform(0, 10); 
	for(i in 1:R){
		beta0_raw[i] ~ normal(0, omega_raw[i]); 
	}
	Sigma_raw ~ inv_wishart(R + 1, Sigma_prior);
	
	# Linear response
	for(i in 1:N){
		for(j in 1:R){
			Mu[i,j]	<- dot_product(X[i], beta_raw[j]) + beta0_raw[j,hh_idx[i]] ; 
		}
	}
	
	# The distribution of discrete outcome
	for(i in 1:N){
		z[i] ~ multi_normal(row(Mu, i)', Sigma_raw); 
	}
	for(i in 1:N){
		for(j in 1:R){
			y[i,j] ~ bernoulli(Phi(z[i,j])); 
		}
	}
}
"

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
sel			<- paste("DOLP_", gsub("\\s", "_", fmt_name), sep="")
hh_exp$y 	<- rowSums(hh_exp[,sel])
sel1		<- gsub("DOLP", "SHR", sel)
for(i in 1:length(fmt_name)){
	hh_exp[,sel1[i]] <- hh_exp[,sel[i]]/hh_exp$y
}
sel			<- which(names(hh_exp) == "i")
hh_exp		<- hh_exp[,-sel]

# Conver some date and factor variables
hh_exp$month <- month(as.Date("2004-1-1", format="%Y-%m-%d") + 14*(hh_exp$biweek-1))
hh_exp$famsize		<- factor(hh_exp$famsize, levels = c("Single","Two", "Three+"))
hh_exp$Q			<- hh_exp$food_quant + hh_exp$nonedible_quant
hh_exp[is.na(hh_exp$y),"y"] <- 0
hh_exp[is.na(hh_exp$dol_purchases),"dol_purchases"] <- 0
hh_exp[is.na(hh_exp$num_day),"num_day"] <- 0
cpi 				<- unique(price_dat[,c("year","cpi")])

# Segment households based on their initial income level 
panelist	<- data.table(hh_exp)
setkeyv(panelist, c("household_code","year","biweek"))
# panelist	<- panelist[,list(first_income = income_group[1], income = income_real[1], first_famsize = famsize[1]), by=list(household_code)]
panelist	<- panelist[,list(income = income_real[1], first_famsize = famsize[1]), by=list(household_code)]
tmp			<- quantile(panelist$income, c(0, .33, .67, 1))
num_grp		<- 3
panelist	<- panelist[, first_income := cut(panelist$income, tmp, labels = paste("T", 1:num_grp, sep=""), include.lowest = T)]
hh_exp		<- merge(hh_exp, data.frame(panelist)[,c("household_code", "first_income")], by="household_code", all.x = TRUE)
cat("Table of initial income distribution:\n"); print(table(panelist$first_income)); cat("\n")
cat("Table of segments in the expenditure data:\n"); print(table(hh_exp$first_income)); cat("\n")

#---------------------------------------#
# Construct regression data 
mydata 	<- hh_exp
sel			<- paste("DOLP_", gsub("\\s", "_", fmt_name), sep="")
sel1		<- gsub("DOLP", "IC", sel)
for(i in 1:length(fmt_name)){
	mydata[,sel1[i]] <- 1*(mydata[,sel[i]]>0)
}
mydata$ln_income	<- log(mydata$income_midvalue/mydata$cpi)
mydata$month		<- factor(mydata$month)
mydata$ln_dolpurchases <- log(mydata$dol_purchases)
mydata$recession 	<- factor(mydata$recession)
mydata$first_income	<- factor(mydata$first_income, levels = c("T2", "T1", "T3"))
mydata$year			<- as.factor(mydata$year)

# Add lagged purchases
mydata	<- data.table(mydata)
setkeyv(mydata, c("household_code", "biweek"))
mydata	<- mydata[,lag_dol := my_lag(dol_purchases), by = list(household_code)]
mydata	<- data.frame(mydata)
mydata	<- mydata[!is.na(mydata$lag_dol),]

# Model matrix: drop convenience store
N 	<- nrow(mydata)
y 	<- as.matrix(mydata[, paste("IC_", gsub("\\s", "_", fmt_name[-1]), sep="")])
R	<- ncol(y)
X	<- as.matrix(cbind(1, mydata$ln_income, mydata$lag_dol/100))
K	<- ncol(X)
tmp	<- unique(mydata$household_code)
nhh	<- length(tmp)
hh_idx	<- do.call(c, lapply(1:nhh, function(i) rep(i, sum(mydata$household_code == tmp[i]))))

# Make initial values
selcol		<- paste("IC_", gsub("\\s", "_", fmt_name[-1]), sep="")
for(i in 1:R){
	fml		<- as.formula(paste(selcol[i], "~ ln_income + I(lag_dol/100) + (1|household_code)", sep=""))
	(tmpfit	<- glmer(fml, family = binomial(link = "probit"), data = mydata))
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
 
exp_init	<- list(beta_raw 	= beta_init*sqrt(2), 
		beta0_raw	= beta0_init*sqrt(2), 
		omega_raw	= omega_init*sqrt(2),
		Sigma_raw	= cor(y),
		z			= z_init
		)
my_init		<- function(){return(exp_init)}

# Prior
Sigma_prior <- diag(R)

# Data for stan, parameters to save
mydata		<- list(N = N, R = R, K = K, nhh = nhh, X = X, y = y, Sigma_prior = Sigma_prior, hh_idx = hh_idx)
save_par	<- c("beta", "Sigma","omega", "beta0")
system.time(foo			<- stan(model_code = src_code, data = mydata, pars = save_par, init = my_init, chains = 1, iter = 1))
pct 		<- proc.time()

CL 			<-  makeCluster(4)
clusterExport(cl = CL, c("foo", "src_code", "mydata", "save_par", "my_init", "exp_init")) 
fitls <- parLapply(CL, 1:4, fun = function(cid) {  
  	require(rstan)
	stan(fit = foo, data = mydata, pars = save_par, init = my_init,
					chains = 1, iter = 120, warmup = 20, thin = 4, chain_id = cid)
})
stopCluster(CL)
use.time	<- proc.time() - pct
cat("Stan sampline took", use.time[3]/60,"min.\n")

fit		<- sflist2stanfit(fitls)
print(fit)

save.image("2_het_income_stan.rdata")
cat("This program is done. \n")


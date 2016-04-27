
cat("This program begins to run at", as.character(Sys.time()), ".\n")

library(ggplot2)
library(reshape2)
library(evd)
library(data.table)
library(doParallel)
library(nloptr)
library(rstan)

rstan_options(auto_write = TRUE)
options(error = quote({dump.frames(to.file = TRUE)}))
options(mc.cores = parallel::detectCores())

args <- commandArgs(trailingOnly = TRUE)
print(args)
if(length(args)>0){
    for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
    }
}

# arr_id	<- as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))
# arr_id	<- 22
Nsec		<- 20				# Divide the full data by subsections
seg_id	<- ceiling(arr_id/Nsec)
arr_idx	<- arr_id - (seg_id - 1) * Nsec
cat("seg_id =", seg_id, ", arr_idx =", arr_idx, ".\n")

# setwd("~/Documents/Research/Store switching/processed data")
# plot.wd	<- '~/Desktop'
# source("../Exercise/Multiple_discrete_continuous_model/0_Allocation_function.R")

# setwd("/home/brgordon/ccv103/Exercise/run")
setwd("/kellogg/users/marketing/2661703/Expenditure")
# setwd("/sscc/home/c/ccv103/Exercise/run")
# setwd("E:/Users/ccv103/Documents/Research/Store switching/run")

run_id			<- 6#7
cpi.adj			<- TRUE

# Load estimation data 
load(paste("estrun_",run_id,"/MDCEV_cpiest_shrpar.rdata",sep=""))
shr.par		<- shr.par.mat[,seg_id]

# Set stan simulation parameters 
numchain	<- 4
# numiter		<- 1500
# numburn		<- 750
# numthin		<- 3
numiter		<- 1000
numburn		<- 500
numthin		<- 2

# numchain	<- 2
# numiter		<- 20
# numburn		<- 2
# numthin		<- 2

load("hh_biweek_exp.rdata")

# Extract 5% random sample
# length(unique(hh_exp$household_code))
# sel			<- sample(unique(hh_exp$household_code), .01*length(unique(hh_exp$household_code)) )
# hh_exp_save	<- hh_exp
# hh_exp		<- subset(hh_exp, household_code %in% sel)

####################
# Stan source code #
####################
src_code	<- '
data{
	int<lower = 1> N; 
	int<lower = 1> nh;
	int<lower = 1> R; 
	int<lower = 1> nx;
	int<lower = 1> beta0_base; 
	int	 		   idx[N]; 
	vector[N]		shr[R]; 
	vector[N]		shr_sgn[R]; 
	vector[N]		p[R]; 
	vector[N]		y; 
	matrix[N,nx]	X[R]; 	
	vector[N]	M; 
	vector[N]	lones; 
}
parameters{
	vector[R-1]	mu_raw; 
	vector<lower = 0>[R-1] tau_raw; 
	vector[nx]	beta; 
	vector<lower = 0>[R] gamma; 
	vector[R]	beta0_raw[nh]; 
	real<lower = 0> sigma; 
}
transformed parameters{
	vector[R]	mu; 
	vector[R]	tau; 
	vector[R]	beta0[nh]; 
	vector[N]	beta0_re[R]; 
	
	for(j in 1:R){
		if(j==beta0_base){
			mu[j]	<- 0; 
			tau[j]	<- 0;
		}else{
			if(j < beta0_base){
				mu[j]	<- mu_raw[j];
				tau[j]	<- tau_raw[j];
			}else{
				mu[j]	<- mu_raw[(j-1)];
				tau[j]	<- tau_raw[(j-1)]; 
			}
		}
	}
	
	for(i in 1:nh){
		for(j in 1:R){
			beta0[i,j]	<- mu[j] + tau[j]*beta0_raw[i,j]; 
		}
	}
	
	for(i in 1:R){
		for(j in 1:N){
			beta0_re[i,j]	<- beta0[idx[j], i]; 
		}
	}
}
model{
	// Define working objects; 
	vector[N]	tot; 
	vector[N]	d[R]; 
	vector[N]	V[R]; 
	vector[N]	pd; 		// sum_r d*shr_sgn;
	vector[N]	pld; 		// sum_r log(d)*shr_sgn; 
	vector[N]	sev;		// sum_r exp(V/sigma); 
	vector[N]	pV; 		// sum_r V/sigma*shr_sgn
	pd	<- rep_vector(0, N); 
	pld	<- rep_vector(0, N); 
	sev	<- rep_vector(0, N); 
	pV	<- rep_vector(0, N); 
	
	// Prior; 
	mu_raw	~ normal(0, 100);
	tau_raw	~ uniform(0, 10); 
	beta	~ normal(0, 100); 
	gamma	~ uniform(0, 100); 
	sigma	~ uniform(0, 10); 
	for(i in 1:nh){
		beta0_raw[i] ~ normal(0, 1); 
	}
	
	// Likelihood; 
	for(j in 1:R){
		d[j]	<- shr[j] + gamma[j] * p[j]./y; 
		V[j]	<- X[j] * beta + beta0_re[j] - log(d[j]/gamma[j]); 
		pd		<- pd + d[j] .* shr_sgn[j]; 
		pld		<- pld + log(d[j]) .* shr_sgn[j]; 
		sev		<- sev + exp(V[j]/sigma); 
		pV		<- pV + V[j].* shr_sgn[j]/sigma ;
	}
	tot	<- log(pd) - pld - (M - lones)*log(sigma) + pV - M .* log(sev); 

	increment_log_prob(sum(tot)); 
}
'

###############
# Subset data # 
###############
# Retail formats 
fmt_name 	<- as.character(sort(unique(fmt_attr$channel_type)))
R			<- length(fmt_name)

# Convert some date and factor variables
hh_exp$month <- month(as.Date("2004-1-1", format="%Y-%m-%d") + 14*(hh_exp$biweek-1))
hh_exp$famsize		<- factor(hh_exp$famsize, levels = c("Single","Two", "Three+"))

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

# Subset data by income group
selcol		<- c("household_code", "biweek", "dol", "year", "month", "income_midvalue", "ln_inc","first_incomeg", 
				"scantrack_market_descr",paste("SHR_", gsub("\\s", "_", fmt_name), sep=""))
mydata		<- subset(hh_exp[,selcol], as.numeric(first_incomeg) == seg_id & dol > .1 )
ord			<- order(mydata$household_code, mydata$biweek)
mydata		<- mydata[ord,]

# Extract subsection data for MCMC
subpan		<- ceiling(c(0, c(1:Nsec)/Nsec * length(unique(mydata$household_code))))
sel			<- (subpan[arr_idx]+1):subpan[(arr_idx+1)]
nh			<- length(unique(mydata[,"household_code"]))
idx			<- setNames(1:nh, unique(mydata[,"household_code"]))
mydata$id	<- idx[as.character(mydata$household_code)]
mydata	<- subset(mydata, id %in% sel)
cat("The dimension of sub-data is", dim(mydata), "\n")

# Check income variation 
tmp		<- data.table(mydata)
tmp		<- tmp[, list(income = unique(income_midvalue)), by = list(household_code, year)]
tmp		<- tmp[, ':='(within_dif = income - mean(income), first_dif = income - income[1]), by = list(household_code)]
cat("Varition of within-hh income in the subset:", sd(tmp$within_dif), "\n")

rm(list = c("args", "hh_exp", "ord"))

################################
# Organize data for estimation #
################################
# The data required for estimation: 
# outcome variables: y (expenditure), shr (expenditure share);
# explanatory variables: price, X_list

# Format attributes
fmt_attr 		<- subset(fmt_attr, year > 2003)
fmt_attr$name 	<- paste(fmt_attr$scantrack_market_descr,fmt_attr$year, sep="-")
fmt_attr$ln_num_module 	<- log(fmt_attr$num_module)
fmt_attr$ln_upc_per_mod <- log(fmt_attr$avg_upc_per_mod)

# Match biweekly prices
price_dat 		<- subset(price_dat, year > 2003)
tmp		<- price_dat
tmp$name<- gsub("\\s", "_", tmp$channel_type)
tmp$name<- paste("PRC_", tmp$name, sep="")
price 	<- dcast(tmp, scantrack_market_descr + year ~ name, value.var = "bsk_price_paid_2004")
# Merge price data to the trip data
mydata <- merge(mydata, price, by=c("scantrack_market_descr", "year"), all.x=T)

# Outcome variables as matrix
beta0_base 	<- which(fmt_name == "Grocery")						# The retailer for which the intercept is set 0 
N		<- nrow(mydata)
idx		<- mydata$id - min(mydata$id) + 1
nh		<- length(unique(idx))
sel		<- paste("SHR_", gsub("\\s", "_", fmt_name), sep="")
shr		<- as.matrix(mydata[,sel])
shr_sgn	<- 1 * (shr > 0)
M		<- rowSums(shr_sgn)
y		<- as.vector(mydata$dol)
ln_inc	<- mydata$ln_inc
sel		<- paste("PRC_", gsub("\\s", "_", fmt_name), sep="")
price	<- as.matrix(mydata[,sel])

# Match retailers' attributes
tmp1	<- unique(fmt_attr$year)
tmp2	<- unique(fmt_attr$scantrack_market_descr)
tmp		<- paste(rep(tmp2, each=length(tmp1)), rep(tmp1, length(tmp2)), sep="-")
tmpn	<- 1:length(tmp)		
names(tmpn) <- tmp
sel 	<- with(mydata, paste(scantrack_market_descr,year, sep="-"))
sel1	<- tmpn[sel]
selcol	<- c("size_index", "ln_upc_per_mod", "ln_num_module","overall_prvt")
nx 		<- length(selcol) * 2 + R-1

X	<- array(NA, c(R, N, nx))
for(i in 1:length(fmt_name)){
	sel2		<- fmt_attr$channel_type == fmt_name[i]
	tmp			<- fmt_attr[sel2,selcol]
	tmp1 		<- as.matrix(tmp[sel1,])
	tmp2		<- matrix(0, nrow(shr), R-1)
	if(i < beta0_base){
		tmp2[,i]	<- ln_inc
	}else if(i > beta0_base){
		tmp2[,(i-1)] <- ln_inc
	}
	X[i,,]	<- cbind(tmp2, tmp1, tmp1 * ln_inc)
}

rm(list = c("tmp", "tmp1", "tmp2", "sel", "sel1"))

############
# Run MCMC # 
############
# Stan argument #
data_ls	<- list(N=N, nh=nh, R=R, nx=nx, beta0_base=beta0_base, idx=idx, shr=t(shr), shr_sgn=t(shr_sgn), 
				p=t(price), y=y, X=X, M=M, lones = rep(1, N))
save_par<- c("beta", "gamma", "mu", "tau", "sigma", "beta0")				
init_ls	<- list(mu_raw = shr.par[paste("beta0_", setdiff(1:R, beta0_base), sep="")], tau_raw = rep(1, R-1), 
				beta = c(rep(0, R-1), shr.par[paste("beta_",1:(nx-R+1), sep="")]), gamma = shr.par[paste("gamma_", 1:R, sep="")], 
				beta0_raw = matrix(0, nh, R), sigma = exp(shr.par["ln_sigma"]))
my_init	<- function(){return(init_ls)}

# pct		<- proc.time()
# sim 	<- stan(model_code = src_code, data = data_ls, pars = save_par, init = my_init, 
# 	chains = 1, iter = numiter, warmup = numburn, thin = numthin)
# use.time	<- proc.time() - pct
# cat("The simulation finishes with", use.time[3]/60, "min.\n")

prt		<- proc.time()

CL = makeCluster(2)
clusterExport(cl = CL, c("src_code", "data_ls", "save_par", "my_init", "init_ls", "numiter", "numburn", "numthin")) 
fitls <- parLapply(CL, 1:4, fun = function(cid) {  
  	require(rstan)
	stan(model_code = src_code, data = data_ls, pars = save_par, init = my_init,
		refresh = ifelse(cid == 1, 200, -1), 
		chains = 1, iter = numiter, warmup = numburn, thin = numthin, chain_id = cid)
})
# fitls	<- mclapply(1:4, mc.cores = 4, function(cid) {  
#   	require(rstan)
# 	stan(model_code = src_code, data = data_ls, pars = save_par, init = my_init,
# 		refresh = ifelse(cid == 1, 200, -1), 
# 		chains = 1, iter = numiter, warmup = numburn, thin = numthin, chain_id = cid)
# })
use.time1	<- proc.time() - prt
stopCluster(CL)
cat("Stan sampling takes", use.time1[3]/60,"min.\n")
sim		<- sflist2stanfit(fitls)

# Print the mean fixed effects. 
cat("Summary of beta:\n"); monitor(extract(sim, "beta", permuted = FALSE, inc_warmup = FALSE)); cat("\n")
cat("Summary of mu:\n"); monitor(extract(sim, "mu", permuted = FALSE, inc_warmup = FALSE)); cat("\n")
cat("Summary of tau:\n"); monitor(extract(sim, "tau", permuted = FALSE, inc_warmup = FALSE)); cat("\n")
cat("Summary of gamma:\n"); monitor(extract(sim, "gamma", permuted = FALSE, inc_warmup = FALSE)); cat("\n")

# Save results 
if(cpi.adj){
	fname	<- paste("estrun_",run_id,"/MDCEV_stan_cpi_seg",seg_id,"_sub", arr_idx,
				"_chain",numchain,"_iter", numiter, "_burn", numburn, "_thin", numthin,sep="")
}else{
	fname	<- paste("estrun_",run_id,"/MDCEV_stan_seg",seg_id,"_sub", arr_idx,
				"_chain",numchain,"_iter", numiter, "_burn", numburn, "_thin", numthin,sep="")
}
fname

# Plot the sampling results 
pdf(paste(fname, ".pdf", sep=""), width = 10, height = 10)
for(i in setdiff(save_par, "beta0")){
	print(traceplot(sim, pars = i))
}
dev.off()

save.image(file = paste(fname, ".rdata", sep=""))

cat("This program is done. \n")

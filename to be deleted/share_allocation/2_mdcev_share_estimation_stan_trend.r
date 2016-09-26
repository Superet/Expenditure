
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

arr_id	<- as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))
# arr_id	<- 2
seg_id	<- ceiling(arr_id/3)
arr_idx	<- arr_id - (seg_id - 1) * 3
cat("seg_id =", seg_id, ", arr_idx =", arr_idx, ".\n")

# setwd("~/Documents/Research/Store switching/processed data")
# plot.wd	<- '~/Desktop'
# source("../Exercise/Multiple_discrete_continuous_model/0_Allocation_function.R")

# setwd("/home/brgordon/ccv103/Exercise/run")
# setwd("/kellogg/users/marketing/2661703/Exercise/run")
setwd("/sscc/home/c/ccv103/Exercise/run")

source("0_Allocation_function.R")

# Load estimation data 
run_id		<- 4
make_plot	<- TRUE
ver.date	<- "2016-02-26"
cpi.adj		<- TRUE

if(cpi.adj){
	loadf	<- paste("estrun_",run_id,"/MDCEV_cpi_est_seg",seg_id,"_", ver.date,".rdata",sep="")
}else{
	loadf	<- paste("estrun_",run_id,"/MDCEV_est_seg",seg_id,"_", ver.date,".rdata",sep="")
}
loadf
load(loadf)

rm(list = setdiff(ls(), c("beta0_base", "cpi.adj", "seg_id", "run_id", "R", "fmt_name", "fmt_attr",
			"ln_inc", "mydata", "nx", "price", "X_list", "shr", "shr.par", "X_list", "y", "arr_idx")))

# # Extract subsection data
ord		<- order(mydata$household_code, mydata$biweek)
mydata	<- mydata[ord,]
Nsec		<- 3				# Divide the full data by subsections
subpan		<- ceiling(c(0, c(1:Nsec)/Nsec * length(unique(mydata$household_code))))
sel			<- (subpan[arr_idx]+1):subpan[(arr_idx+1)]
nh			<- length(unique(mydata[,"household_code"]))
idx			<- setNames(1:nh, unique(mydata[,"household_code"]))
mydata$id	<- idx[as.character(mydata$household_code)]
mydata_sub	<- subset(mydata, id %in% sel)
cat("The dimension of sub-data is", dim(mydata_sub), "\n")
# mydata_sub 	<- mydata

# Set stan simulation parameters 
numchain	<- 4
numiter		<- 10000
numburn		<- 5000
numthin		<- 20

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
	int<lower = 1> idx[N]; 
	int<lower = 1> year_idx[N]; 
	int<lower = 1> month_idx[N]; 
	matrix[N,R]		shr; 
	matrix[N,R]		shr_sgn; 
	matrix[N,R]		p; 
	vector[N]		y; 
	matrix[N,nx]	X[R]; 	
	vector[N]	M; 
	vector[R]	ones; 
	vector[N]	lones;
}
parameters{
	vector[R-1]	mu_raw; 
	vector<lower = 0>[R-1] tau_raw; 
	vector[R-1]	mu_year_raw; 
	vector<lower = 0>[R-1] tau_year_raw;
	vector[nx]	beta; 
	vector<lower = 0>[R] gamma; 
	vector[R]	beta0_raw[nh]; 
	vector[R]	beta0_year_raw[nh]; 
	vector[R]	beta0_month_raw[nh]; 
	real<lower = 0> sigma; 
}
transformed parameters{
	vector[R]	mu; 
	vector[R]	tau; 
	vector[R]	beta0[nh]; 
	
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
}
model{
	// Define working objects; 
	vector[N]	J; 
	vector[N]	tot; 
	matrix[N,R]		d; 
	matrix[N,R]		V; 
	
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
	for(i in 1:N){
		for(j in 1:R){
			d[i,j]	<- shr[i,j] + gamma[j] * p[i,j]/y[i];
			V[i,j]	<- X[j,i]*beta + beta0[idx[i],j] - log(d[i,j]/gamma[j]); 
		}
	}
	J	<- log((d .* shr_sgn)*ones) - (log(d) .* shr_sgn)*ones; 
	tot	<- J - (M-lones)*log(sigma) + (V/sigma .* shr_sgn)*ones -M .* log(exp(V/sigma)*ones); 
	increment_log_prob(sum(tot)); 
}
'

############
# Run MCMC # 
############
# Organize data 
# Outcome variables as matrix
N		<- nrow(mydata_sub)
idx		<- mydata_sub$id - min(mydata_sub$id) + 1
nh		<- length(unique(idx))
sel		<- paste("SHR_", gsub("\\s", "_", fmt_name), sep="")
shr		<- as.matrix(mydata_sub[,sel])
shr_sgn	<- 1 * (shr > 0)
M	<- rowSums(shr_sgn)
y		<- as.vector(mydata_sub$dol)
ln_inc	<- mydata_sub$ln_inc
sel		<- paste("PRC_", gsub("\\s", "_", fmt_name), sep="")
price	<- as.matrix(mydata_sub[,sel])

# Match retailers' attributes
tmp1	<- unique(fmt_attr$year)
tmp2	<- unique(fmt_attr$scantrack_market_descr)
tmp		<- paste(rep(tmp2, each=length(tmp1)), rep(tmp1, length(tmp2)), sep="-")
tmpn	<- 1:length(tmp)		
names(tmpn) <- tmp
sel 	<- with(mydata_sub, paste(scantrack_market_descr,year, sep="-"))
sel1	<- tmpn[sel]
selcol	<- c("size_index", "ln_upc_per_mod", "ln_num_module","overall_prvt")
nx 		<- length(selcol) * 2

X	<- array(NA, c(R, N, nx))
for(i in 1:length(fmt_name)){
	sel2		<- fmt_attr$channel_type == fmt_name[i]
	tmp			<- fmt_attr[sel2,selcol]
	tmp1 		<- as.matrix(tmp[sel1,])
	X[i,,]	<- cbind(tmp1, tmp1 * ln_inc)
}

# ------------- #
# Stan argument #
data_ls	<- list(N=N, nh=nh, R=R, nx=nx, beta0_base=beta0_base, idx=idx, shr=shr, shr_sgn=shr_sgn, 
				p=price, y=y, X=X, M=M, ones = rep(1, R), lones = rep(1, N))
save_par<- c("beta", "gamma", "sigma", "mu", "tau","beta0")				
init_ls	<- list(mu_raw = shr.par[paste("beta0_", setdiff(1:R, beta0_base), sep="")], tau_raw = rep(0.1, R-1), 
				beta = shr.par[paste("beta_",1:nx, sep="")], gamma = shr.par[paste("gamma_", 1:R, sep="")], 
				beta0_raw = matrix(0, nh, R), sigma = exp(shr.par["ln_sigma"]))
my_init	<- function(){return(init_ls)}

# pct		<- proc.time()
# sim 	<- stan(model_code = src_code, data = data_ls, pars = save_par, init = my_init, 
# 	chains = 1, iter = numiter, warmup = numburn, thin = numthin)
# use.time	<- proc.time() - pct
# cat("The simulation finishes with", use.time[3]/60, "min.\n")

CL = makeCluster(4)
clusterExport(cl = CL, c("src_code", "data_ls", "save_par", "my_init", "init_ls", "numiter", "numburn", "numthin")) 
prt		<- proc.time()
fitls <- parLapply(CL, 1:4, fun = function(cid) {  
  	require(rstan)
	stan(model_code = src_code, data = data_ls, pars = save_par, init = my_init, 
		chains = 1, iter = numiter, warmup = numburn, thin = numthin, chain_id = cid)
})
use.time1	<- proc.time() - prt
stopCluster(CL)
cat("Stan sampling takes", use.time1[3]/60,"min.\n")
sim		<- sflist2stanfit(fitls)

# Save results 
run_id		<- 5
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
plot(sim)
dev.off()

save.image(file = paste(fname, ".rdata", sep=""))

cat("This program is done. \n")

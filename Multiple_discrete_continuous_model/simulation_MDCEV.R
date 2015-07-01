# MDCEV simulation 

library(reshape2)
library(ggplot2)
library(Rcpp)
library(RcppArmadillo)
library(maxLik)
library(evd)

setwd("~/Documents/Research/Store switching/Exercise/Multiple_discrete_continuous_model")
model_name <- "MDCEV_a1b1"

sourceCpp(paste(model_name, ".cpp", sep=""))
source("0_Allocation_function.R")

# Set paraemters 
R		<- 3 		# Number of alternatives
Ra		<- R		# Number of alternatives + number of outside options
exp_outside <- quant_outside <- FALSE
beta0 	<- c(0, -1, -1)
beta	<- c(.5, -.7)
gamma0 	<- gamma	<- c(1, 1, 1)
sigma 	<- 1
if(substr(model_name,7,8)=="a2"){
	Ra 	<- Ra + 1
	exp_outside <- TRUE
	gamma	<- c(gamma, 1)
}
if(substr(model_name,9,10)=="b2"){
	Ra 	<- Ra + 1
	quant_outside <- TRUE
	gamma <- c(gamma, 1)
}
qz_cons	<- 10
psi_R1	<- 1
psi_R2	<- .01

# Simulate data 
set.seed(666666)
nx 		<- length(beta)
N 		<- 500		# Number of observations
X_arr 	<- array(rnorm(N*R*nx), c(R, N, nx))
X_list  <- lapply(1:R, function(i) X_arr[i,,])
price 	<- matrix(runif(N*R, 2, 4), N, R)
Q		<- runif(N, 1, 20)
y		<- rowSums(price) * Q/R 

par(mfrow=c(3,1))
hist(y, breaks=100)
hist(Q, breaks=100)
hist(as.vector(price), breaks=100)

eps_draw<- matrix(rgumbel(N*R), N, R)
xbeta	<- do.call(cbind, lapply(1:R, function(i) X_list[[i]] %*% beta + beta0[i]))
psi		<- exp(xbeta + eps_draw)
if(substr(model_name,7,8)=="a2"){
	psi	<- cbind(psi, psi_R1)
}
if(substr(model_name,9,10)=="b2"){
	psi	<- cbind(psi, psi_R2)
}
e_mat <- matrix(NA, N, Ra)
for(i in 1:N){
	tmp	<- Allocation_fn(y = y[i], psi = psi[i,], gamma, Q = Q[i], price = price[i,], R, Ra, qz_cons, exp_outside, quant_outside)
	e_mat[i,] <- tmp$e
}
eps <- 1e-4
eR 	<- e_mat[,1:R]
eR[eR<=eps] <- 0
mean(eR==0)

sel	<- apply(eR, 1, function(x) !all(x==0))
eR 	<- eR[sel,]
Q	<- Q[sel]
y	<- y[sel]
price <- price[sel,]
X_list <- lapply(X_list, function(x) x[sel,])
e1_index <- apply(eR, 1, function(x) which(x== min(x[x>0])))

# Estimation 
if(substr(model_name,7,8)=="a1"){
	if(substr(model_name,9,10)=="b2"){
		MDCEV_wrapper <- function(param) {
			MDCEV_LogLike_fnC(param, nx, c_q=qz_cons, Q=Q, e=eR, e1_index = e1_index, 
			p=price, X_list = X_list)
		}
	}else{
		MDCEV_wrapper <- function(param){
			MDCEV_LogLike_fnC(param, nx, e1_index = e1_index, e=eR, p=price, X_list = X_list)
		}
	}
}else{
	if(substr(model_name,9,10)=="b2"){
		MDCEV_wrapper <- function(param) {
			MDCEV_LogLike_fnC(param, nx, c_q=qz_cons, y = y, Q=Q, e=eR, p=price, X_list = X_list)
		}
	}else{
		MDCEV_wrapper <- function(param){
			MDCEV_LogLike_fnC(param, nx, y=y, e=eR, p=price, X_list = X_list)
		}
	}
}

# Estimation with multiple initial values. 
# theta_init	<- list(c(beta, beta0[-1]), 
# 					c(1, -1, 2, -1.5), 		# Close to the true values
# 					c(-3, -3, 4, -1)  
# 					)
# Estimation with multiple initial values (full set of parameters)
theta_init0	<- c(beta, beta0, gamma0, log(psi_R1), log(psi_R2), log(sigma))
theta_init2	<- theta_init1 <- theta_init0
theta_init1[1:(nx+R)] <- theta_init0[1:(nx+R)] + rnorm(nx+R)
theta_init2[1:(nx+R)] <- rep(5, nx+R)
theta_init 	<- list(theta_init0, theta_init1, theta_init2)
system.time(tmp <- MDCEV_wrapper(theta_init[[1]]))

tmp_sol <- vector("list", length(theta_init))
pct <- proc.time()
for(i in 1:length(theta_init)){
	names(theta_init[[i]]) <- c(paste("beta_",1:nx, sep=""), paste("beta0_",1:R,sep=""), paste("gamma",1:R, sep=""),
								"ln_psi1", "ln_psi2", "ln_sigma")	
	tmp_sol[[i]]	<- maxLik(MDCEV_wrapper, start=theta_init[[i]], method="BFGS")
}
sel 	<- which.max(sapply(tmp_sol, function(x) ifelse(is.null(x$maximum), NA, x$maximum)))
sol		<- tmp_sol[[sel]]
lapply(tmp_sol, summary)

# Fix scale parameter
myfix <- length(theta_init0)
tmp_sol <- vector("list", length(theta_init))
for(i in 1:length(theta_init)){
	names(theta_init[[i]]) <- c(paste("beta_",1:nx, sep=""), paste("beta0_",1:R,sep=""), paste("gamma",1:R, sep=""),
								"ln_psi1", "ln_psi2", "ln_sigma")	
	tmp_sol[[i]]	<- maxLik(MDCEV_wrapper, start=theta_init[[i]], method="BFGS", fixed = myfix)
}
lapply(tmp_sol, summary)

# Fix scale parameter and psi for outside options 
myfix <- (length(theta_init0) - 2):length(theta_init0)
tmp_sol <- vector("list", length(theta_init))
for(i in 1:length(theta_init)){
	names(theta_init[[i]]) <- c(paste("beta_",1:nx, sep=""), paste("beta0_",1:R,sep=""), paste("gamma",1:R, sep=""),
								"ln_psi1", "ln_psi2", "ln_sigma")	
	tmp_sol[[i]]	<- maxLik(MDCEV_wrapper, start=theta_init[[i]], method="BFGS", fixed = myfix)
}
lapply(tmp_sol, summary)

# Fix scale parameter and psi for outside options and satiation parameters
myfix <- (length(theta_init0) - (2+R)):length(theta_init0)
tmp_sol <- vector("list", length(theta_init))
for(i in 1:length(theta_init)){
	names(theta_init[[i]]) <- c(paste("beta_",1:nx, sep=""), paste("beta0_",1:R,sep=""), paste("gamma",1:R, sep=""),
								"ln_psi1", "ln_psi2", "ln_sigma")	
	tmp_sol[[i]]	<- maxLik(MDCEV_wrapper, start=theta_init[[i]], method="BFGS", fixed = myfix)
}
lapply(tmp_sol, summary)

# Fix scale parameter and psi for outside options, satiation parameters, and one beta intercept
theta_init[[2]][(nx+1)] <- theta_init[[3]][(nx+1)] <- 0
myfix <- c(nx+1, (length(theta_init0) - (2+R)):length(theta_init0))
tmp_sol <- vector("list", length(theta_init))
for(i in 1:length(theta_init)){
	names(theta_init[[i]]) <- c(paste("beta_",1:nx, sep=""), paste("beta0_",1:R,sep=""), paste("gamma",1:R, sep=""),
								"ln_psi1", "ln_psi2", "ln_sigma")	
	tmp_sol[[i]]	<- maxLik(MDCEV_wrapper, start=theta_init[[i]], method="BFGS", fixed = myfix)
}
lapply(tmp_sol, summary)


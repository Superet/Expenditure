# Simulation of allocation model 

# Model: 
# max_{e_s} \sum_{s=1}^{S} gamma_s*psi_s*log(e_s/p_s + 1) + psi_z*log(e_z) + psi_w*log(e_w)
# s.t. 	sum_{s=1}^{S} e_s + e_z <= y
# 	 	sum_{s=1}^{S} e_s/p_s + e_w <= Q

# Normalization: psi_z = psi_w = 1

library(ggplot2)	
library(reshape2)
library(scatterplot3d)

#####################
# Utility functions #
#####################
uP_fn <- function(e, psi, gamma,price, S){
# Negative purchase utility 
# e: S+2 vecotr
# psi, gamma, price is (S+2)-dimension vector 
	e0 		<- e[1:S]
	e_S1	<- e[(S+1):length(e)]
	u0		<- c(log(e0/(gamma[1:S]*price) + 1), log(e_S1[1]), log(e_S1[2]+qz_cons))
	out 	<- sum(psi*gamma*u0 ) 
	return(-out)
}

uPGrad_fn <- function(e, psi, gamma,price, S){
	e0 		<- e[1:S]
	e_S1	<- e[(S+1):length(e)]
	u0		<- c(1/(e0+price*gamma[1:S]), 1/e_S1[1], 1/(e_S1[2]+qz_cons) )
	out 	<- psi*gamma*u0 
	return(-out)
	
}

# Allocation function - solve a constrained optimization
Allocation_fn <- function(y, psi, gamma, Q=NULL, price,S,Sa, inits=NULL){
	# We have Sa non-negativity constraints, one budget constraint, and one quantity constraint
	# ui <- rbind(diag(S+2), c(rep(-1, S+1), 0) )
	# inits is a list of initial values
	ui <- rbind(diag(Sa), c(rep(-1, S+1),rep(0, Sa-S-1 )) )
	ci <- c(rep(0,S+1), -qz_cons, -y)
	if(is.null(inits)){
		tmp	  <- min(y/price, Q) * .99
		inits <- list(c(rep(tmp/(Sa-1), Sa-1), Q - tmp), c(rep(1e-8, S), y-(S+2)*1e-8, Q - sum(price)*1e-8 - 1e-8))
	}
		
	if(is.null(Q)){
		sol <- constrOptim(theta = inits, f=uP_fn, grad=NULL, ui=ui, ci=ci,psi=psi, gamma=gamma, price=price, S=S)
		return(list(e = sol$par, max = -sol$value))
	}else{
		if(Q<=0){
			e <- c(rep(0,S),y)
			if(Sa>S+1){
				e <- c(e, 0)
			}
			return(list(e = e, max = log(y) ) )
		}else{
			ui <- rbind(ui, c(-1/price, 0, -1))
			ci <- c(ci, -Q)
			sol.list <- vector("list", length(inits))
			for(j in 1:length(inits)){
				sol.list[[j]] <- constrOptim(theta = inits[[j]], f=uP_fn, grad=uPGrad_fn, ui=ui, ci=ci,psi=psi, gamma=gamma, price=price, S=S)
			}
			sel <- which.min(sapply(sol.list, function(x) x$value))
			sol <- sol.list[[sel]]
			if(sol$convergence != 0){
				cat("Constrained optimization does not converge at value y=",y,", Q=",Q,"\n")
			}
			return(list(e = sol$par, max = -sol$value, convergence = sol$convergence))
		}
	}
}

###############################
# Allocation rule exploration # 
###############################
# Set base parameters 
S 	<- 2
Sa	<- S+2
psi <- c(1.5, 2, 1, .00001)
gamma <- rep(1, Sa)
qz_cons <- 2
y 	<- 5
Q 	<- 2
price <- c(2, 3)

sol <- Allocation_fn(y, psi, gamma, Q=Q, price=price, S, Sa)

# Set experiment parameters 
psi1	<- c(5, 10, 1, 1)
gamma1 	<- rep(1, Sa)
sol1 	<- Allocation_fn(y, psi1, gamma1, Q=Q, price=price, S, Sa)

# Compare the allocation 
tmp <- rbind(sol$e, sol1$e)
rownames(tmp) <- c("Baseline", "Experiment")
cat("The expenditure at two parameters: \n"); round(tmp, 4)

############################
# Shape of inclusive value #
############################
y_grid 	<- seq(0.1, 10, by=.2)
Q_grid 	<- seq(.1, 10, by=0.2)
ny		<- length(y_grid)
nQ		<- length(Q_grid)
e_arr	<- array(NA, c(ny, nQ, Sa), dimnames = list(y_grid, Q_grid, 1:Sa))
omega	<- matrix(NA, ny, nQ, dimnames= list(y_grid, Q_grid))
for(i in 1:ny){
	for(j in 1:nQ){
		tmp <- Allocation_fn(y = y_grid[i], psi, gamma, Q=Q_grid[j], price=price, S, Sa)
		e_arr[i,j,] <- tmp$e
		omega[i,j]  <- tmp$max
	}
	print(i)
}

# 3D plots of inclusive values, expenditure allocation. 
my.angle <- 60
ggtmp <- melt(omega)
names(ggtmp) <- c("y","Q","value")
scatterplot3d(x=ggtmp$y, y=ggtmp$Q, z=ggtmp$value, xlab="y", ylab="Q", zlab="omega", main="Utility function(inclusive value)",angle=my.angle)

ggtmp <- melt(e_arr[,,(S+1)])
names(ggtmp) <- c("y","Q","value")
quartz()
scatterplot3d(x=ggtmp$y, y=ggtmp$Q, z=ggtmp$value, xlab="y", ylab="Q", zlab="e3", main="Expenditure of the outside option 1",angle=my.angle)

ggtmp <- melt(e_arr[,,Sa])
names(ggtmp) <- c("y","Q","value")
quartz()
scatterplot3d(x=ggtmp$y, y=ggtmp$Q, z=ggtmp$value, xlab="y", ylab="Q", zlab="q4", main="Quantity of the outside option 2",angle=my.angle)

# 2D plots
sel <- c(1, nQ/4, nQ/2, 3*nQ/4, nQ)
ggtmp <- melt(omega[,sel])
names(ggtmp) <- c("y","Q","value")
quartz()
ggplot(ggtmp, aes(y, value)) + geom_point() + geom_line() + 
		facet_wrap(~Q) + 
		labs(title="Total utility at selected quantity constraints")

ggtmp <- melt(e_arr[,sel,(S+1)])
names(ggtmp) <- c("y","Q","value")
quartz()
ggplot(ggtmp, aes(y, value)) + geom_point() + geom_line() + 
		geom_abline(xintercept=0, slope=1, linetype=2) + 
		facet_wrap(~Q) + 
		labs(title="Expenditure of outside option 1\n at selected quantity constraints")

ggtmp <- melt(e_arr[,sel,Sa])
names(ggtmp) <- c("y","Q","value")
quartz()
ggplot(ggtmp, aes(y, value)) + geom_point() + geom_line() + 
		geom_hline(aes(yintercept=Q), linetype=2) + 
		facet_wrap(~Q) + 
		labs(title="Quantity of outside option 2\n at selected quantity constraints")


##################
# Add randomness # 
##################
set.seed(99)
numsim 	<- 50
y_grid 	<- seq(0.1, 10, by=.2)
Q_grid 	<- seq(.1, 6, by=0.5)
ny		<- length(y_grid)
nQ		<- length(Q_grid)
eps		<- matrix(rnorm(numsim*S), numsim, S)
psi_draw <- rep(1, numsim)%*%t(psi) * cbind(exp(eps), 1, 1)
e_arr_d	<- array(NA, c(numsim, ny, nQ, Sa), dimnames = list(1:numsim, y_grid, Q_grid, 1:Sa))
omega_d	<- array(NA, c(numsim, ny, nQ), dimnames= list(1:numsim, y_grid, Q_grid))
for(i in 1:numsim){
	for(j in 1:ny){
		for(k in 1:nQ){
			tmp <- Allocation_fn(y = y_grid[j], psi=psi_draw[i,], gamma, Q=Q_grid[k], price=price, S, Sa)
			e_arr_d[i,j,k,] <- tmp$e
			omega_d[i,j,k]  <- tmp$max
		}
	}
	print(i)
}
omega <- apply(omega_d, c(2,3), mean) 
e_arr <- apply(e_arr_d, c(2,3,4), mean)

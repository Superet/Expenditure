library(reshape2)
library(ggplot2)

olrp <- function(beta, A,B,Q, R, W=NULL){
	m 	<- max(dim(A))
	rb	<- nrow(B)
	cb	<- ncol(B)
	if(is.null(W)) {W = matrix(0, m, cb) }
# 	if(min(abs(eigen(R)$values)) > 1e-5){
# 		stop("eigen value error");
# 	}else{
		p0 <-  - .01*diag(m)
		dd <- 1
		it <- 1
		maxit <- 1000
		while (dd>1e-6 & it<= maxit){
			f0 	= solve(R+ beta * t(B) %*%p0 %*%B) %*% (beta* t(B)%*%p0%*%A + t(W))
			p1	= beta*t(A) %*% p0 %*%A + Q - (beta * t(A) %*% p0 %*% B+ W) %*% f0
			f1 	= solve(R + beta* t(B) %*% p1 %*% B) %*% (beta* t(B)%*%p1%*%A + t(W))
			dd	= max(abs(f1 - f0));
			it = it +1
			p0 = p1
		}
# 	}
	return(list(f = f1, p = p0))
}

myfn <- function(x, lambda1, lambda2){
	lambda2*x^2+ lambda1*x
}

lambda_w1	= 5;
lambda_w2 	= -.1;
lambda_w3 	= 4; 
lambda_w4	= - .1; 
lambda_w5	= - .5; 
lambda_c2 	= -.2;
lambda_c1	= 6;
lambda_o2 	= 0;
lambda_o1	= 4;
tau_I1		= 5; 
tau_I2 		= -.1;
tau_I3 		= -.5; 
tau_S2		= -.2;
tau_S1		= 5; 
k 			= .5;
rho 		= .6; 
beta= .95;
lmu			= c(-Inf,  1.2);
lsigma		= c(0, .2);
mu 			= exp(lmu + .5*lsigma^2)
Sigma 		= diag((exp(lsigma^2) -1)*exp(2*lmu + lsigma^2));
ns			<- 2; 
nu			<- 3; 

# Plot the utilify functions
x 			<- seq(1, 10, by=.1)
par(mfrow= c(2,2))
plot(x, myfn(x, lambda_w1, lambda_w2))
plot(x, myfn(x, lambda_c1, lambda_c2))
plot(x, myfn(x, lambda_o1, lambda_o2))
plot(x, myfn(x, tau1, tau2))

A 	= matrix(c(1, 0, 0, rho), ns,ns,byrow=T)
B	= matrix(c(k, 0, -1, 0, 0, 0), ns, nu, byrow=T)

R2	<- matrix(0, ns, ns); R2[1,1] <- -tau_I2; R2[2,2] <- lambda_o2; 
R1	<- matrix(c(-tau_I1, lambda_o1), 1, ns)
Q2	<- matrix(c(lambda_w2 + lambda_o2 - tau_I2*k^2, lambda_w5 - tau_I3*k, 	tau_I2*k, 
				lambda_w5 - tau_I3*k,				lambda_w4 - tau_S2, 	tau_I3, 
				tau_I2*k,							tau_I3, 				lambda_c2 - tau_I2), nu, nu, byrow=T)
Q1	<- matrix(c(lambda_w1 - lambda_o1-tau_I1*k, lambda_w3 - tau_S1, lambda_c1 + tau_I1), 1, nu)
W	<- matrix(c(-tau_I2*k,	-tau_I3, 	tau_I2, 
				-lambda_o2, 	0,			0), ns, nu, byrow=T)


sol <- olrp(beta, A, B, R2, Q2, W)
F 	<- sol$f
P 	<- sol$p
S 	<- (R1+ Q1 %*% F + 2*beta*t(mu) %*% t(P) %*% (A+ B%*% F)) %*% solve(diag(ns) - beta*A - beta*B %*% F)
S 	<- matrix(S, nrow = 1)
G	<- -.5 * solve(Q2 + beta*t(B) %*% P %*% B) %*% (2*beta*t(B) %*% P %*% mu + beta*t(B) %*% t(S) + t(Q1) )
d 	<- (.5*(Q1 + beta* S %*% B + 2*beta*t(mu) %*% P %*% B) %*% G + beta* S %*% mu + beta * sum(diag(P %*% Sigma)) )/(1-beta)
round(P, 3) 
round(S, 3)
d

# Simulate the data 
TT 			<- 100
init_state	<- c(I=3, B= 5, Inc=1)
state 		<- matrix(NA, TT, length(init_state))
state[1,]	<- init_state
ubar		<- matrix(NA, TT, 2)
ustar		<- matrix(NA, TT, 2)
zdraw		<- rlnorm(TT, lmu[2], lsigma[2])
eps_inc		<- rlnorm(TT, lmu[3], lsigma[3])
for(i in 1:TT){
	ubar[i,] 	<- F %*% state[i,] + G
	ymax 		<- state[i,2] + state[i,3] - zdraw[i]
	ustar[i,1]	<- max(0, min(ubar[i,1], ymax))
	cmax 		<- state[i,1] + k*ustar[i,1]
	ustar[i,2] 	<- max(0,min(ubar[i,2], cmax))
	if(i< TT){
		state[(i+1), 1] <- state[i,1] + k*ustar[i,1] - ustar[i,2]
		state[(i+1), 2] <- state[i,2] + state[i,3] - zdraw[i] - ustar[i,1]
		state[(i+1), 3] <- rho * state[i,3] + eps_inc[i]
	}
}

# Plot the state
tmp		<- data.frame(t = 1:TT, state, ustar)
names(tmp) <- c("t","I","B","Inc","y","c")
ggtmp	<- melt(tmp, id.var="t")
ggtmp$sol	<- "realized"
tmp		<- melt(data.frame(t=1:TT, y=ubar[,1], c=ubar[,2]), id.var="t")
ggtmp	<- rbind(ggtmp, cbind(tmp, sol="star"))
ggplot(ggtmp, aes(t, value, col=sol)) + geom_point() + geom_line() + 
		facet_grid(variable ~. , scales="free_y")




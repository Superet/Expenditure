# Draw utility function 
# Utility 
# u = psi1*log(x1+1) + psi2*log(x2+1) + log(q+1) + beta*log(x1+x2-q+1)
library(ggplot2)

### Functions ###
u_fn <- function(x1,x2){
	q <- (1-beta+x1+x2)/(1+beta)			# Optimal consumption quanity in the first peirod
	psi1*log(x1+1) + psi2*log(x2+1) + log(q+1) + beta*log(x1+x2-q+1)
}

xstar_fn <- function(y,x.init=c(0.01,0.01)){
	obj_fn <- function(xstar){
		x1 <- xstar[1]
		x2 <- xstar[2]
		-u_fn(x1,x2)
	}
	grad_fn <- function(xstar){
		x1 <- xstar[1]
		x2 <- xstar[2]
		dudx1 <- -(1+beta)/(2+x1+x2) - psi1/(1+x1)
		dudx2 <- -(1+beta)/(2+x1+x2) - psi2/(1+x2)
		as.vector(c(dudx1,dudx2))
	}
	ui <- rbind(-matrix(c(p1,p2),1,2),c(1,0),c(0,1))
	constrOptim(x.init,obj_fn,grad=grad_fn,ui=ui,ci=c(-y,0,0))
}

xstar_limit_fn <- function(y,x.init=c(0.01,0.01)){
	x.interior <- xstar_fn(y)$par
	if(x.interior[2]<x2.limit){
		xstar <- c(y/p1,0)
	}else{
		xstar <- x.interior
	}
	xstar
}

id_fn <- function(x1,ubar){
	solve_fn <- function(x2){
		(u_fn(x1,x2) - ubar)^2
	}	
	grad_fn <- function(x2){
		dudx2 <- (1+beta)/(2+x1+x2) + psi2/(1+x2)
		dobj_x1 <- 2*(u_fn(x1,x2) - ubar)*dudx2
		as.vector(dobj_x1)
	}
	constrOptim(1, solve_fn, grad=grad_fn, ui=matrix(1,1,1),ci=0)$par
}

### Simulation ###
# Set parameters 
psi1 <- 1
psi2 <- 1
p1 <- 1.5
p2 <- 1
beta <- .5

y <- 5
y.new <- 2.4
x1.dis <- c(0,.4*y/p1,.9*y/p1)
x2.dis <- c(0,.5*y/p2)
x.comb <- cbind(rep(x1.dis,length(x2.dis)),rep(x2.dis,each=length(x1.dis)))

x1 <- seq(0.01,y/p1,by=.01)
x2 <- rep(NA,length(x1))
x2.new <- rep(NA,length(x1))
for(i in 1:length(x1)){
	sol <- xstar_fn(y)
	ubar <- -sol$value
	x2[i] <- id_fn(x1[i],ubar)
	
	sol <- xstar_fn(y.new)
	ubar <- -sol$value
	x2.new[i] <- id_fn(x1[i],ubar)	
}
budget <- y/p2 - p1*x1/p2
budget.new <- y.new/p2 - p1*x1/p2


# Plot utility 
pdf("~/Desktop/graph_category_utility.pdf",width=5, height=5)
par(mai=c(.8,.8,.5,.5))
plot(x1,x2,type="l",xlim=c(-.01,y/p1),ylim=c(-.01,y/p2))
lines(x1,budget,col="red")
lines(x1,budget.new,col="red",lty="dashed")
lines(x1,x2.new,lty="dashed")
points(x.comb[,1],x.comb[,2],pch=16)
dev.off()

### Comparative statistics ###
x2.limit <- 2
y.seq <- seq(.1,10,by=.01)
x.opt <- matrix(NA,length(y.seq),2)
for(i in 1:length(y.seq)){
	x.opt[i,] <- xstar_limit_fn(y.seq[i])
}
w.opt <- rep(1,length(y.seq)) %*% t(c(p1,p2)) * x.opt
w.opt <- w.opt/(y.seq %*% t(c(1,1)))

# Comparative statist w.r.t beta
y <- 3
beta.seq <- seq(.01,1,by=.01)
x.opt1 <- matrix(NA,length(beta.seq),2)
for(i in 1:length(beta.seq)){
	beta <- beta.seq[i]
	x.opt1[i,] <- xstar_limit_fn(y)
}
w.opt1 <- rep(1,length(beta.seq)) %*% t(c(p1,p2)) * x.opt1
w.opt1 <- w.opt1/y

par(mfrow=c(2,1))
plot(y.seq,w.opt[,1],type="l",xlab="y",ylab=expression(w[1]))
plot(y.seq,w.opt[,2],type="l",xlab="y",ylab=expression(w[2]))

tmp1 <- data.frame(variable="y",iv=y.seq,w2=w.opt[,2])
tmp2 <- data.frame(variable="beta",iv=beta.seq,w2=w.opt1[,2])
ggtmp <- rbind(tmp1,tmp2)


pdf("~/Desktop/graph_compstat.pdf",width=5,height=5)
ggplot(ggtmp,aes(iv,w2)) + geom_line() + 
		facet_wrap(~variable,scales="free",ncol=1) + 
		labs(x="Value",y=paste(expression(w2),", expenditure share at store 2",sep=""),title="Comparative statistics")
dev.off()



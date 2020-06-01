library(Rcpp)
library(RcppGSL)
# library(benchmark)

# sourceCpp("~/Documents/Research/Store switching/Exercise/Basket DP/cspline3d.cpp")
sourceCpp("cspline3d.cpp")

# The right R function.
mysplR <- function(x, y, x0){
	n <- nrow(x0)
	out <- rep(NA, n)
	x1 	<- unique(x[,1])
	x2	<- unique(x[,2])
	x3	<- unique(x[,3])
	n1	<- length(x1)
	n2	<- length(x2)
	n3	<- length(x3)
	
	for(i in 1:n){
# 		cat(i, "\n")
		x3_index <- which(x3==x0[i,3])
		y2 	<- rep(NA, n2)
		for(j in 1:n2){
			index <- (x3_index-1)*n1*n2 + (j-1)*n1 + seq_len(n1)
			y2[j] <- spline(x1, y[index], method="natural", xout=x0[i,1])$y
# 			cat("--",j)
		}
		out[i] <- spline(x2, y2, method="natural", xout=x0[i,2])$y
	}
	return(out)
}

# Test interpolation 
matrix_expand <- function(mat, vec){
	n1	<- nrow(mat)
	n2	<- length(vec)
	sel	<- rep(1:n1, n2)
	mat_new <- cbind(mat[sel,], rep(vec, each=n1))
	return(mat_new)
}

my_grid <- list(I = seq(0, 500, by=50), 
				y = seq(0, 1000, by=200), 
				Inc	= c(1, 2, 3) )
state	<- as.matrix(my_grid[[1]], ncol=1) 
for(i in 2:length(my_grid)){
	state <- matrix_expand(state, my_grid[[i]])
}
my_scale <- 100			# Scale down the inventory and expenditure
state[,c(1,2)] <- state[,c(1,2)]/my_scale
summary(state)

V <- 2*log(state[,1] + state[,3] +1)
xnew <- c(3.6, 2.5, 1)

my_fn(state, V, xnew)
my_fnS(state, V, xnew)
mysplR(state, V, matrix(xnew, nrow=1))

# microbenchmark(my_fn(state, V, xnew), my_fnS(state, V, xnew), mysplR(state, V, matrix(xnew, nrow=1)))

# A sequence of interpolation 
xnew <- rbind(xnew, xnew)
my_fn1(state, V, xnew)
my_fnS1(state, V, xnew)
mysplR(state, V, xnew)

xnew <- cbind(runif(20, 2, 5), runif(20, 6, 9), sample(my_grid[[3]], 20, replace=T))
my_fn1(state, V, xnew)
my_fnS1(state, V, xnew)
mysplR(state, V, xnew)

# microbenchmark(my_fn(state, V, xnew), my_fnS(state, V, xnew), mysplR(state, V, xnew))

cat("This program is done.\n")

# Check chebyshev interpolation in the package of chebpol
# Chebyshev is good for approximating a smooth function, but bad for approximating an oscillating function 
library(chebpol)
library(microbenchmark)

cheb.basis	<- function(x, M, interval = c(0, 1), adj = TRUE){
	n		<- length(x)
	if(adj){
		x.adj	<- 2*(x - interval[1])/(interval[2] - interval[1]) - 1
	}else{
		x.adj	<- x
	}
	
	out	<- matrix(NA, n, M)
	out[,1]	<- 1
	out[,2]	<- x.adj
	for(i in 3:M){
		out[,i]	<- 2 * x.adj * out[,(i-1)] - out[,(i-2)]
	}
	return(out)
}

cheb.1d.basis <- function(x, M, interval = c(0, 1)){
	n		<- length(x)
	x.adj	<- 2*(x - interval[1])/(interval[2] - interval[1]) - 1
	# T		<- cheb.basis(x = x.adj, M=M, interval = interval, adj = FALSE)
	T		<- matrix(NA, n, M)
	out		<- matrix(NA, n, M)
	T[,1]	<- 1
	T[,2]	<- x.adj
	out[,1]	<- 0
	out[,2]	<- 1
	for(i in 3:M){
		T[,i]	<- 2 * x.adj * T[,(i-1)] - T[,(i-2)]
		out[,i]	<- 2*x.adj*out[,(i-1)] + 2*T[,(i-1)] - out[,(i-2)]
	}
	return(out)	
}

chebfun	<- function(x, y, interval = c(0, 1), ... ){
	df	<- length(x)
	T	<- cheb.basis(x = x, M = df, interval = interval, adj = TRUE)
	alpha	<- solve(t(T) %*% T) %*% t(T) %*% y
	f	<- function(x, deriv = 0, dT = NULL){
		if(deriv == 0){
			T1	<- cheb.basis(x = x, M = df, interval = interval, adj = TRUE)
			out	<-as.vector(T1 %*% alpha)
		}else if(deriv == 1){
			if(is.null(dT)){
				dT	<- cheb.1d.basis(x = x, M = df, interval = interval)
			}
			out	<- dT %*% (alpha * 2 /(interval[2] - interval[1]))
		}	
		return( out )
	}
	return(f)
}

# Example: spline interpolation
n <- 9
x <- 1:n
f 	<- function(x){ sin((x-0.5)*pi) }
f	<- function(x){ - .05 * x^3 + .5* x^2 + 1/x}
f.d	<- function(x){  -.15 * x^2 + x - 1/x^2}
y 	<- f(x)

plot(x, y, main = paste("spline[fun](.) through", n, "points"))
lines(spline(x, y))
lines(spline(x, y, n = 201), col = 2)

# Chebshev interpolation
numnodes	<- 9
x.interval 	<- range(x)
knots	<- chebknots(numnodes, intervals = range(x))[[1]]
knots.a	<- (chebknots(numnodes)[[1]] + 1) * (x.interval[2] - x.interval[1]) * .5 + x.interval[1]
identical(knots, knots.a)
knots.y		<- f(knots)

# Compare chebshev approximation from the chebpol package and my own function
f.a		<- chebfun(knots, knots.y, interval = x.interval)
f.b		<- chebappx(knots.y, interval = x.interval)
f.s		<- splinefun(knots, knots.y)

grid 	<- seq(-1, n + 2, .4)
grid	<- seq(x.interval[1], x.interval[2], .4)
a		<- f.a(grid)
b		<- sapply(grid, f.b)
identical(a, b)
max(abs(a-b))

plot(grid, f(grid), main = "Interpolation f")
lines(grid, a, col = "red")
lines(grid, b, col = "blue")
lines(grid, splinefun(knots, knots.y)(grid), col = 'green')

# Compare speed 
microbenchmark(chebfun(knots, knots.y, interval = x.interval), chebappx(knots.y, interval = x.interval))
microbenchmark(f.a(grid), sapply(grid, f.b), f.s(grid))

# Compute the derivative 
dT	<- cheb.1d.basis(x = grid, M = length(knots), interval = x.interval)
a	<- f.a(grid, deriv = 1)
b	<- f.s(grid, deriv = 1)
a1	<- f.a(grid, deriv = 1, dT = dT)
max(abs(a - b))
identical(a, a1)

# Plot the derivatives
plot(grid, f.d(grid), main = "Derivatives of f")
points(grid, a, col = "red")
points(grid, b, col = "blue")

microbenchmark(f.a(grid, deriv = 1), f.s(grid, deriv = 1))
microbenchmark(f.a(grid, deriv = 1, dT = dT), f.s(grid, deriv = 1))



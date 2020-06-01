# This file compares a few optimizers that solve constrained optimization. 
# constrOptim, npltr, alabama, contrOptim.nl, solnp

library(nloptr)
library(alabama)
library(Rsolnp)
library(microbenchmark)
set.seed(12)

#############
# Example 1 #
#############
# This is the example from auglag
fn <- function(x) (x[1] + 3*x[2] + x[3])^2 + 4 * (x[1] - x[2])^2
 
gr <- function(x) {
	g <- rep(NA, 3)
	g[1] <- 2*(x[1] + 3*x[2] + x[3]) + 8*(x[1] - x[2]) 
	g[2] <- 6*(x[1] + 3*x[2] + x[3]) - 8*(x[1] - x[2]) 
	g[3] <- 2*(x[1] + 3*x[2] + x[3])
	g
}
 
hin <- function(x) {
# For comparison, I remove nonlinear constraints	
	h <- rep(NA, 1)
	# h[1] <- 6*x[2] + 4*x[3] - x[1]^3 - 3
	h[1] <- 6*x[2] + 4*x[3] - 3	
	h[2] <- x[1]
	h[3] <- x[2]
	h[4] <- x[3]
	h
}
 
hin.jac <- function(x) {
	j <- matrix(NA, 4, length(x))
	# j[1, ] <- c(-3*x[1]^2, 6, 4)
	j[1, ] <- c(0, 6, 4)	
	j[2, ] <- c(1, 0, 0)
	j[3, ] <- c(0, 1, 0)
	j[4, ] <- c(0, 0, 1)
	j
}

gin		<- function(x) { -hin(x) }		# Constraints in nloptr are g(x) <= 0
gjac	<- function(x) { -hin.jac(x)}

# Linear contraints in contrOptim
ui 	<- matrix(c(1, 0, 0,
				0, 1, 0, 
				0, 0, 1, 
				0, 6, 4), 4, 3, byrow = T)
ci 	<- c(0,0,0,3)

# Run optimization 
init	<- runif(3)
sol.lt	<- setNames(vector("list", 5), c("constrOptim", "alabama", "contrOptim.nl", "solnp", "npltr"))
sol.lt[[1]]	<- constrOptim(init, fn, gr, ui = ui, ci = ci, method = "BFGS")
sol.lt[[2]]	<- auglag(par=init, fn=fn, gr=gr, hin=hin, hin.jac=hin.jac)
sol.lt[[3]]	<- constrOptim.nl(par=init, fn=fn, gr=gr, hin=hin, hin.jac=hin.jac)
sol.lt[[4]]	<- solnp(pars = init, fun = fn, ineqfun = hin, ineqLB = rep(0,4), ineqUB = rep(Inf, 4))
sol.lt[[5]]	<- nloptr(x0=init, eval_f = fn, eval_grad_f=gr, lb = rep(-Inf,length(init)), ub = rep(Inf, length(init)),
							eval_g_ineq = gin, eval_jac_g_ineq=gjac, 
							opts = list("algorithm"="NLOPT_LD_MMA"))
							
# Summarize results 
val.name	<- c("value", "value", "value", "values", "objective" )
par.name	<- c(rep("par", 3), "pars", "solution")
iter.name	<- c(rep("outer.iterations", 3), "outer.iter", "iterations")
tmp		<- sapply(1:length(sol.lt), function(i) min(sol.lt[[i]][[val.name[i]]]))
tmp1	<- sapply(1:length(sol.lt), function(i) sol.lt[[i]][[par.name[i]]])
tmp2	<- sapply(1:length(sol.lt), function(i) sol.lt[[i]][[iter.name[i]]])
tmp.tab	<- rbind(value = tmp, x = tmp1, iter = tmp2)
round(tmp.tab, 3)

# Compare processing time
microbenchmark(
constrOptim(init, fn, gr, ui = ui, ci = ci, method = "BFGS"), 
auglag(par=init, fn=fn, gr=gr, hin=hin, hin.jac=hin.jac),
constrOptim.nl(par=init, fn=fn, gr=gr, hin=hin, hin.jac=hin.jac), 
solnp(pars = init, fun = fn, ineqfun = hin, ineqLB = rep(0,4), ineqUB = rep(Inf, 4)), 
nloptr(x0=init, eval_f = fn, eval_grad_f=gr, lb = rep(-Inf,length(init)), ub = rep(Inf, length(init)),
							eval_g_ineq = gin, eval_jac_g_ineq=gjac, 
							opts = list("algorithm"="NLOPT_LD_MMA"))
	)



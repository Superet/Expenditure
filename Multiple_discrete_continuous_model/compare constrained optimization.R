# This file compares methods of constrained optimization that solves expenditure allocation 

# Model: 
# U = max_{e} sum_r psi_r*gamma_r*log(e_r/(gamma_r*price_r) +1)
# s.t. e_r >= 0 for all r, and sum_r e_r <= y

# Methods: 
# 1. constrOptim to solve utility maximization subject to e_r>=0 for all r, sum_r e_r <= y
# 2. set sum_r e_r = y, only have constraints e_r >=0,
# 	 note that this method may yield negative expenditure e.
# 3. on top of method 2, use variable transformation e_r = exp(x_r) and solve the problem in an unconstrained optimization, 
#	 note that this method may yield negative expenditure e 
# 4. on top of method 3, add one constraint sum_r e_r <= y
# 5. The same as method 1, but use nloptr. 

# Conclusion: 
# a) Method 2 and Mehtod 3 may yield solutions that do not satisfy constriants, since the two methods do not explicitly account for the constraints. 
# b) constrOptim and nloptr perfrom on par, but nolptr yields more robust solutions.

library(microbenchmark)
library(nloptr)

# Purchae utility function 
uP_fn <- function(e, psi, gamma, price, R, Ra, qz_cons = 0, exp_outside = FALSE, quant_outside = FALSE){
# e: Ra-dimentional vector
# psi, gamma is Ra-dimension vector, price is R-dimensional vector.
# R is length of the effective varaible in the optimization algorithm.
# Ra is the length of variable in the original optimization problem. 
# exp_outside and quant_outside are indicator whether there is expenditure/quantity outside option. 
	e0 		<- e[1:R]
	e_R1	<- e[(R+1):Ra]
	u0_part1<- log(e0/(gamma[1:R]*price) + 1)
	if(exp_outside){
		if(quant_outside){
			u0_part2	<- c(log(e_R1[1]/gamma[(R+1)]), log(e_R1[2]/gamma[Ra]+qz_cons) )
		}else{
			u0_part2	<- log(e_R1[1]/gamma[(R+1)])
		}
	}else{
		if(quant_outside){
			u0_part2	<- log(e_R1/gamma[Ra]+qz_cons)
		}else{
			u0_part2	<- NULL
		}
	}
	
	u0		<- c(u0_part1, u0_part2 )
	out 	<- sum(psi*gamma*u0 ) 
	return(out)
}

# Gradient of utility function w.r.t e
uPGrad_fn <- function(e, psi, gamma, price, R, Ra, qz_cons = 0, exp_outside = FALSE, quant_outside = FALSE){
	e0 		<- e[1:R]
	e_R1	<- e[(R+1):Ra]
	u0_part1<- 1/(e0+price*gamma[1:R])
	if(exp_outside){
		if(quant_outside){
			u0_part2	<- c(1/e_R1[1], 1/(e_R1[2]+qz_cons*gamma[Ra])  )
		}else{
			u0_part2	<- 1/e_R1[1]
		}
	}else{
		if(quant_outside){
			u0_part2	<- 1/(e_R1 + qz_cons*gamma[Ra])
		}else{
			u0_part2	<- NULL
		}
	}
	
	u0		<- c(u0_part1, u0_part2 )
	out 	<- psi*gamma*u0 
	return(out)	
}

# Mothod 1: constrOptim
f1 <- function(y, psi, gamma, price, R, e_init){
	# Construct the constrained equations: ui %*% theta - ci >= 0
	ui	<- diag(R)				# Ra non-negativity constraints. 
	ci 	<- rep(0, R)
	ui 	<- rbind(ui, rep(-1,R))
	ci	<- c(ci, -y)
	
	sol <- constrOptim(theta = e_init, f=uP_fn, grad=uPGrad_fn, ui=ui, ci=ci, control = list(fnscale = -1), psi=psi,
				 gamma=gamma, price=price, R=R, Ra=Ra)
	return(list(e = sol$par, max = sol$value, convergence = sol$convergence))
}

# Method 2: set sum_r e_r = y, then only optimize over R-1 free variable. 
f2	<- function(y, psi, gamma, price, R, e_init, base_idx = 1){
	sub_uP_fn <- function(sube, y, base_idx, psi, gamma, price, R, Ra, qz_cons = 0, exp_outside = FALSE, quant_outside = FALSE){
		e	<- rep(0, length(sube)+1)
		e[base_idx]	<- y - sum(sube)
		e[-base_idx]	<- sube
		return(uP_fn(e, psi, gamma, price, R, Ra, qz_cons, exp_outside, quant_outside))
	}

	sub_uPGrad_fn <- function(sube, y, base_idx, psi, gamma, price, R, Ra, qz_cons = 0, exp_outside = FALSE, quant_outside = FALSE){
		e	<- rep(0, length(sube)+1)
		e[base_idx]	<- y - sum(sube)
		e[-base_idx]	<- sube
		gr	<- uPGrad_fn(e, psi, gamma, price, R, Ra, qz_cons, exp_outside, quant_outside)
		out	<- gr - gr[base_idx]
		out	<- out[-base_idx]
		return(out)
	}
	
	ui	<- diag(R-1)
	ci	<- rep(0, R-1)
	sub_init	<- e_init[-base_idx]
	sol	<- constrOptim(theta = sub_init, f=sub_uP_fn, grad=sub_uPGrad_fn, ui=ui, ci=ci, control = list(fnscale = -1), 
				base_idx = base_idx, y =y, psi=psi, gamma=gamma, price=price, R=R, Ra=Ra)
	sube 	<- sol$par
	e	<- rep(0, length(sube)+1)
	e[base_idx]	<- y - sum(sube)
	e[-base_idx]	<- sube
	return(list(e = e, max = sol$value, convergence = sol$convergence))
}

# Method 3: set sum_r e_r = y and use variable transformation e_r = exp(le_r)
f3 <- function(y, psi, gamma, price, R, e_init, base_idx = 1){
	sub_uP_fn <- function(sub_le, y, base_idx, psi, gamma, price, R, Ra, qz_cons = 0, exp_outside = FALSE, quant_outside = FALSE){
		sube<- exp(sub_le)
		e	<- rep(0, length(sube)+1)
		e[base_idx]	<- y - sum(sube)
		e[-base_idx]	<- sube
		return(uP_fn(e, psi, gamma, price, R, Ra, qz_cons, exp_outside, quant_outside))
	}

	sub_uPGrad_fn <- function(sub_le, y, base_idx, psi, gamma, price, R, Ra, qz_cons = 0, exp_outside = FALSE, quant_outside = FALSE){
		sube<- exp(sub_le)
		e	<- rep(0, length(sube)+1)
		e[base_idx]	<- y - sum(sube)
		e[-base_idx]	<- sube
		gr	<- uPGrad_fn(e, psi, gamma, price, R, Ra, qz_cons, exp_outside, quant_outside)
		out	<- gr - gr[base_idx]
		out	<- out[-base_idx]
		out	<- out*sube							# Gradient needs to be multiplied by jacobian. 
		return(out)
	}
	
	sub_init	<- log(e_init[-base_idx])
	sol	<- optim(sub_init, f=sub_uP_fn, gr=sub_uPGrad_fn, control = list(fnscale = -1), 
				base_idx = base_idx, y =y, psi=psi, gamma=gamma, price=price, R=R, Ra=Ra)
	sube 	<- exp(sol$par)
	e	<- rep(0, length(sube)+1)
	e[base_idx]	<- y - sum(sube)
	e[-base_idx]	<- sube
	return(list(e = e, max = sol$value, convergence = sol$convergence))
}

# Method 4: add one constraint sum_r e_r <= y to method 3
f4 <- function(y, psi, gamma, price, R, e_init, base_idx = 1){
	sub_uP_fn <- function(sub_le, y){
		sube<- exp(sub_le)
		e	<- rep(0, length(sube)+1)
		e[base_idx]	<- y - sum(sube)
		e[-base_idx]	<- sube
		return(-uP_fn(e, psi, gamma, price, R, Ra, qz_cons =0 , exp_outside = FALSE, quant_outside=FALSE))
	}

	sub_uPGrad_fn <- function(sub_le, y){
		sube<- exp(sub_le)
		e	<- rep(0, length(sube)+1)
		e[base_idx]	<- y - sum(sube)
		e[-base_idx]	<- sube
		gr	<- uPGrad_fn(e, psi, gamma, price, R, Ra, qz_cons=0, exp_outside=FALSE, quant_outside=FALSE)
		out	<- gr - gr[base_idx]
		out	<- out[-base_idx]
		out	<- out*sube							# Gradient needs to be multiplied by jacobian. 
		return(-out)
	}
	
	gin	<- function(sub_le, y){
		sube<- exp(sub_le)
		return(sum(sube) - y)
	}
	
	gin.jac	<- function(sub_le, y){
		out	<- matrix(NA, 1, length(sub_le))
		out[1,]	<- exp(sub_le)
		return(out)
	}
	
	sub_init	<- log(e_init[-base_idx])
	sol	<- nloptr(x0=sub_init, eval_f = sub_uP_fn, eval_grad_f=sub_uPGrad_fn, 
					eval_g_ineq = gin, eval_jac_g_ineq=gin.jac, opts = list("algorithm"="NLOPT_LD_MMA"), 
					y =y)

	sube 	<- exp(sol$solution)
	e	<- rep(0, length(sube)+1)
	e[base_idx]	<- y - sum(sube)
	e[-base_idx]	<- sube
	return(list(e = e, max = -sol$objective, convergence = sol$status))
}

# Method 5: The same as method 1, but use nloptr.
f5 <- function(y, psi, gamma, price, R, e_init, opts = list("algorithm"="NLOPT_LD_MMA")){	
	objf<- function(e){
		-uP_fn(e, psi, gamma, price, R, Ra, qz_cons =0 , exp_outside = FALSE, quant_outside=FALSE)
	}
	obj_gr	<- function(e){
		-uPGrad_fn(e, psi, gamma, price, R, Ra, qz_cons=0, exp_outside=FALSE, quant_outside=FALSE)
	}
	gin	<- function(e){
		return(sum(e) - y)
	}
	
	gin.jac	<- function(e){
		out	<- matrix(NA, 1, length(e))
		out[1,]	<- rep(1, length(e))
		return(out)
	}
	
	sol	<- nloptr(x0=e_init, eval_f = objf, eval_grad_f=obj_gr, lb = rep(0, length(e_init)), ub = rep(Inf, length(e_init)), 
					eval_g_ineq = gin, eval_jac_g_ineq=gin.jac, opts = opts)
	return(list(e = sol$solution, max = -sol$objective, convergence = sol$status))
}

# Set parameters
set.seed(111)
y	<- 90
R	<- 6
Ra 	<- R
psi	<- runif(R, 0, 2)
gamma	<- runif(R, 1, 30)
price	<- runif(R, 2, 4)
e_init	<- psi/sum(psi) * y * .99

# Solve the constrained optimization
sol.list <- vector("list", 5)
sol.list[[1]]	<- f1(y, psi, gamma, price, R, e_init)
sol.list[[2]]	<- f2(y, psi, gamma, price, R, e_init)
sol.list[[3]]	<- f3(y, psi, gamma, price, R, e_init)
sol.list[[4]]	<- f4(y, psi, gamma, price, R, e_init)
sol.list[[5]]	<- f5(y, psi, gamma, price, R, e_init)

# List results
ans		<- sapply(sol.list, unlist)
cat("Solutuions are:\n"); print(round(ans, 5)); cat("\n")
which.max(ans["max",])

# Evaluate processing time
microbenchmark(
f1(y, psi, gamma, price, R, e_init),
f2(y, psi, gamma, price, R, e_init),
f3(y, psi, gamma, price, R, e_init), 
f4(y, psi, gamma, price, R, e_init), 
f5(y, psi, gamma, price, R, e_init)
)

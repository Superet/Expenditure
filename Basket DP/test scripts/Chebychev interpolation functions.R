# Chebyshev polynomial interpolation functions

# 1. Chebyshev interpolation intialization
# 2. Chebyshev polinomial calculation
# 3. Chebyshev multi-dimensional polinomial calculation
# 4. Chebyshev interoloation

## Chebyshev polynomial appximation: 
# 		f(x) = \sum_n theta_n*T_n(x)
# with n={0,...,N}, N is the number of degree of polynomial.

## Chebyshev regression Algorithm:
# Step 1: Choose degree N, and number of nodes M. 
# Step 2: Compute interpolation node v_m=-cos(pi*(2*m-1)/2*M), and adjust nodes to the interval [a,b], xi_m=a+(v_m+1)*(b-a)/2.
# Step 3: Construct basis funciton, T(v), M \times (N+1) matrix. 
# 		  For D-dimensional approximation, expand the basis function as a D-fold tensor product.(M^D) \times((N+1)^D)
# Step 4: Calculate Chebyshev coefficients theta_n=(sum_m f(xi_m)*T_n(v_m) ) / (sum_m T_n(v_m)^2) for each n. 


#########################################
# Chebyshev interpolation intialization #
#########################################
cheb_init <- function(N,M,D,lower,upper){
# Initialize Chebychev interpolation
#------------------------------------------------------
# INPUT
# N := Degree of Chebyshev polynomial;
# M	:= Number of interpolation nodes, note M >= N+1; 
# D	:= Dimension of interpolation
# lower := a vector of the lower bound;
# upper := a vector of the upper bound. 
	if( M < N+1){
		cat("Error: Number of nodes should be greater than degree of polynomial + 1. ")
		return(NULL)
	}
	if(D!=length(lower)){
		cat("Error: Length of bounds is not equal to the dimension D.")
		return(NULL)
	}
	
	# Compute base nodes and polynomial basis
	base.nodes <- sapply(1:M, function(xx) - cos(.5*pi*(2*xx-1)/M))
	basis.polynm <- t(sapply(base.nodes, function(xx) cheb_polynm_scalar(xx,N)))
	if(D>1){
		base.nodes <- Matrix_expand(base.nodes,D)
		base.polynm <- t(apply(base.nodes,1,cheb_polynm_multD,N))
	}
	adj.base.nodes <- sapply(1:D,function(i) .5*(base.nodes[,i]+1)*(upper[i]-lower[i]) + lower[i])
	
	# Store the denominator that is used for calculating Chebychev coefficients
	theta.denominator.inv <- solve(t(base.polynm) %*% base.polynm)
	A <- theta.denominator.inv %*% t(base.polynm)
	
	list(Polynomial.degree=N, Num.nodes=M, Dimension=D,
		 lower=lower, upper=upper, 
		 Unadjusted.base.nodes=base.nodes, 
		 Adjusted.base.nodes=adj.base.nodes,
		 Basis.polynomial=base.polynm,
		 Cheb.coefficient.num = ncol(theta.denominator.inv),
		 Cheb.coefficent.denominator.inv=theta.denominator.inv,
		 Cheb.coefficent.multp=A)
}

cheb_polynm_scalar <- function(x,N){
# Chebychev polynomial of degree N for a scalar x with support [0,1]
# T_0(x) = 1, T_1(x) = x, T_n(x) = 2xT_{n-1}(x) - T_{n-2}(x)
#-------------------------------------------------------------
	Tx <- rep(NA,N+1)
	Tx[1] <- 1
	Tx[2] <- x
	if(N>1){
		for(i in 2:N){
			Tx[i+1] <- 2*x*Tx[i] - Tx[i-1]
		}
	}
	Tx
}

cheb_polynm_multD <- function(x.vec,N){
# Chebychev polynomial of degree N for a multi-dimensional vector: tensor product 
# phi(x.vec) = prod_d T_{i_d,d}(x_d), for i_d={1,...,(N+1)}
#-------------------------------------------------------------------------------
# INPUT
# x.vec \in [0,1]^D
# OUTPUT
# phi: (N+1)^D vector
	D <- length(x.vec)
	base.T <- sapply(x.vec,function(xx) cheb_polynm_scalar(xx,N))
	phi <- base.T[,1]
	for(i in 2:D){
		phi <- phi %x% base.T[,i]
	}
	phi
}

calculate_cheb_coef <- function(cheby_init,y){
	cheby_init$Cheb.coefficent.multp %*% y
}

###########################
# Chebychev Interpolation # 
###########################
cheb_interpolation <- function(x.new,fxi,cheb.setting,theta=NULL){
# Interpolation f(x.new)
# f(x) = sum_n theta_n*T_n(x)
#---------------------------------------------------------------
# INPUT
# x.new	:= a scalor or a vector that is to be interpolated
# fxi	:= function values evalued at the adjusted base nodes(xi)
# cheb.setting := a setting list that has the structure of cheb_init
	
	# Compute Chebychev coefficients
	if(all(theta==NULL)){
		theta.dominator <- t(cheb.setting$Basis.polynomial) %*% fxi
		theta <- cheb.setting$Cheb.coefficent.denominator.inv %*% theta.dominator
	}
	
	# Subfunction: Compute Chebychev polynomial at a single node
	single.interp <- function(x){
		delta <- cheb.setting$upper-cheb.setting$lower
		adj.x <- 2*(x-cheb.setting$lower)/delta-1
		if(cheb.setting$Dimension==1){
			polynm <- cheb_polynm_scalar(adj.x,N=cheb.setting$Polynomial.degree)
		}else{
			polynm <- cheb_polynm_multD(adj.x,N=cheb.setting$Polynomial.degree)
		}
		return(t(polynm) %*% theta)
	}
	
	if(cheb.setting$Dimension==1){
		if(!is.vector(x.new)){
			return(list(value=single.interp(x.new),theta=theta))
		}else{
			return(list(value=apply(x.new,1,single.interp ),theta=theta))
		}
	}else{
		if(is.vector(x.new)){
			return(list(value=single.interp(x.new),theta=theta))
		}else{
			return(list(value=apply(x.new,1,single.interp ),theta=theta))
		}
	}
}

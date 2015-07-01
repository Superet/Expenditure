library(glmmML)

Gaussian_Hermite_init <- function(num.nodes, dimension){
	nodes <- ghq(num.nodes, modified = F)$zeros
	weights <- ghq(num.nodes, modified = F)$weights
	
	if(dimension < 2){
		return(list(nodes = nodes, weights = weights))
	}else{
		nodes <- Matrix_expand(nodes,dimension)
		weights <- Matrix_expand(weights, dimension)
		return(list(nodes = nodes, weights = apply(weights, 1, prod)) )
	}
}

Matrix_expand <- function(base.vec,Dfold){
# Expand the base vector into D-fold matrix
#---------------------------------------------
	mat <- base.vec
	for(i in 2:Dfold){
		mat <- cbind( mat%x%rep(1,length(base.vec)), rep(base.vec,length(base.vec)))
	}
	mat
}
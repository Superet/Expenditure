library(MASS)

n 		<- 1000
mu 		<- c(2,3,4)
Sigma 	<- matrix(c(10, 3, 2, 3, 2, 1.5, 2, 1.5, 2), 3, 3)
x		<- mvrnorm(n, mu, Sigma)
beta 	<- c(3, .5, 1.2, 2)
y 		<- cbind(1, x) %*% beta + rnorm(n)
fit1	<- lm(y ~ x)

mypca 	<- prcomp(x, center = TRUE, scale = TRUE)
x.pc	<- mypca$x
fit2	<- lm(y ~ x.pc)

# check score computation
cc		<- matrix(rep(mypca$center, each=nrow(x)), nrow(x), 3)
ss		<- matrix(rep(mypca$scale, each=nrow(x)), nrow(x), 3)
a		<- ((x - cc )/ss) %*% mypca$rotation
identical(a, x.pc)

# check coefficient relationship
coef(fit1)
coef(fit2)
(mypca$rotation/mypca$scale) %*% coef(fit2)[-1]
coef(fit2)[1] - sum((mypca$rotation*mypca$center/mypca$scale) %*% coef(fit2)[-1])

# Check prediction
x.new 	<- mvrnorm(100, mu, Sigma)
yhat1	<- cbind(1, x.new) %*% coef(fit1)

cc		<- matrix(rep(mypca$center, each=nrow(x.new)), nrow(x.new), 3)
ss		<- matrix(rep(mypca$scale, each=nrow(x.new)), nrow(x.new), 3)
x.newpc	<- ((x.new - cc )/ss) %*% mypca$rotation
yhat2	<- cbind(1, x.newpc) %*% coef(fit2)
max(abs(yhat1 - yhat2))

coef2	<- c(coef(fit2)[1] - sum((mypca$rotation*mypca$center/mypca$scale) %*% coef(fit2)[-1]), 
			(mypca$rotation/mypca$scale) %*% coef(fit2)[-1])
yhat3	<- cbind(1, x.new) %*% coef2
max(abs(yhat1 - yhat3))	
	
			
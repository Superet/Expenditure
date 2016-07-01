library(mvtnorm)
library(reshape2)
library(ggplot2)
library(data.table)

# Model 
# y = beta_0 + beta_1*t

set.seed(9)
n 	<- 10					# Number of individuals
T	<- 20					# Number of periods
rho <- -.4					# Correlation between intercept and slope
ab	<- rmvnorm(n, rep(0, 2), matrix(c(1, rho, rho, 1), 2,2))
ab[,2] <- ab[,2] * .1
x 	<- 1:T
y 	<- lapply(1:n, function(i) ab[i,1] + ab[i,2]*x + 0)
mydata <- melt(y)
mydata <- data.frame(id = mydata$L1, t = rep(1:T, n), x = x, y = mydata$value)
mydata	<- data.table(mydata)
mydata	<- mydata[,within:= c(0, y[-1] - y[-length(y)]), by = list(id)]
mydata	<- mydata[,within:=within + mean(y[t==1])]

ggplot(mydata, aes(x, y)) + 
	geom_line(aes(group = id), size = .25) + 
	geom_smooth(method = "lm") + 
	geom_smooth(data = subset(mydata, t>1) , aes(x, within), method = "lm", color = "red") 
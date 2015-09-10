my.singleFE <- function(fml, data, panid){
	# Demean model matrix by group mean
	X		<- model.matrix(fml, data)
	X 		<- X[,-which(colnames(X)=="(Intercept)")]
	y 		<- model.extract(model.frame(fml, data), "response")
	mydata 	<- data.table(cbind(y, X))
	varvec	<- colnames(mydata)
	mydata$id <- data[,panid]
	
	# Substract group mean
	for(i in 1:length(varvec)){
		var		<- as.name(varvec[i])
		mydata 	<- mydata[,eval(var):=eval(var) - mean(eval(var)), by = list(id)]
	}
	mydata	<- data.frame(mydata)
	
	# Fit the linear regression and adjust standard error
	myfit	<- lm(y ~ . - 1 - id, data = mydata)
	coeftab <- data.frame(summary(myfit)$coefficients)
	coeftab$df			<- myfit$df.residual
	coeftab$adj_df		<- coeftab$df - length(unique(data[,panid]))
	coeftab$adj_StdErr	<- with(coeftab, Std..Error*sqrt(df/adj_df))
	coeftab$adj_t		<- with(coeftab, Estimate/adj_StdErr)
	coeftab$adj_p		<- with(coeftab, pt(-abs(adj_t), adj_df)*2)
	
	return(coeftab)
}
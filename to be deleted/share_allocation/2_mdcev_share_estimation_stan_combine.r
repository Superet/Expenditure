library(parallelMCMCcombine)
library(rstan)
library(stargazer)
library(abind)
library(reshape2)
library(ggplot2)

# setwd("/home/brgordon/ccv103/Exercise/run")
# setwd("/kellogg/users/marketing/2661703/Expenditure")
run_id 	<- 6
setwd(paste("/sscc/home/c/ccv103/Exercise/run/estrun_", run_id, sep=""))

filename	<- list.files(getwd())
filename	<- filename[grep("rdata", filename)]
filename	<- filename[grep("seg", filename)]

# Order the filename based on array index
tmp			<- sapply(filename, function(x) gsub("[^0-9]+", "",substr(x, regexpr("sub", x), regexpr("_chain", x)-1)))
filename	<- filename[order(as.numeric(tmp))]
(filels		<- lapply(1:3, function(i) filename[grep(paste("seg",i,sep=""), filename)]))
# sim.ls	<- vector("list", length = length(filels))
# hhidx.ls<- vector("list", length = length(filels))

sim.combine	<- vector("list", length(filels))
beta0.combine	<- vector("list", length(filels))
hhidx.combine	<- vector("list", length(filels))
pop.par		<- c("beta", "gamma", "mu", "tau", "sigma")	
tmp_ls		<- ls()

for(ii in 1:length(filels)){
	simls		<- list(NULL)
	beta0		<- NULL
	hhidx1		<- NULL
	for(jj in 1:length(filels[[ii]])){
		load(filels[[ii]][jj])
		simls[[jj]]	<- do.call(cbind, extract(sim, pop.par))
		tmp1		<- apply(extract(sim, "beta0")[[1]], c(2,3), median)		# Median of posterior beta0
		beta0		<- rbind(beta0,  tmp1)
		hhidx1		<- c(hhidx1, unique(mydata$id))
		rm(list = setdiff(ls(), c(tmp_ls, "ii", "R", "nx", "fmt_name", "tmp_ls", "jj", "simls", "beta0", "hhidx1")))
	}
	beta0.combine[[ii]]	<- beta0
	hhidx.combine[[ii]]	<- hhidx1
	
	# Combine all the parameters
	sim.arr	<- do.call("abind", c(simls, along = 3))
	sim.combine[[ii]]	<- consensusMCindep(subchain=sim.arr, shuff=FALSE)
	colnames(sim.combine[[ii]]) <- c(paste("beta",1:nx, sep=""), paste(rep(pop.par[2:4], each = R), 1:R, sep=""), "sigma")
	print(ii)
}

# Plot the prameters
med_qt	<- function(x){
	out	<- quantile(x, prob = c(.025, .5,.975))
	names(out)	<- c("ymin", "y","ymax")
	return(out)
}

ggtmp	<- melt(sim.combine)
names(ggtmp)	<- c("iter", "var", "value", "IncGrp")
pdf("stan_combine.pdf", width = 12, height = 12)
ggplot(ggtmp, aes(IncGrp, value)) + 
	stat_summary(fun.data = "med_qt", geom = "pointrange") + 
	facet_wrap(~var, scales = "free")
dev.off()

# Posterior median 
par.med	<- lapply(sim.combine, function(x) cbind(Estimate = apply(x, 2, median), se = apply(x, 2, sd) ))
do.call(cbind, par.med)

save(sim.combine, par.med, beta0.combine, hhidx.combine, file = "MDCEV_stan_est.rdata")
	
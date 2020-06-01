cat("This program begins to run at", as.character(Sys.time()), ".\n")

library(ggplot2)
library(reshape2)
library(maxLik)
library(evd)
library(data.table)
library(mgcv)
library(stargazer)
library(lmtest)
library(scatterplot3d)
options(error = quote({dump.frames(to.file = TRUE)}))

run_id			<- 6
ver.date		<- "2016-10-08"
nseg			<- 3

# setwd("~/Documents/Research/Store switching/Processed_data")
# plot.wd	<- '~/Desktop'
# source("../Exercise/Multiple_discrete_continuous_model/0_Efficient_Allocation_function.R")
# source("../Exercise/main/ctrfact_sim_functions.r")

# setwd("/home/brgordon/ccv103/Exercise/run")
# setwd("/kellogg/users/marketing/2661703/Exercise/run")
setwd("/sscc/home/c/ccv103/Exercise/run")
# setwd("U:/Users/ccv103/Documents/Research/Store switching/run")

# Read estimation from all segments # 
expdata.full	<- data.frame()
exp.est.ls		<- vector("list", length = nseg)
omega.lsf		<- vector("list", length = nseg)
gamma.ls		<- vector("list", length = nseg)
tmp.ls			<- ls()

for(ii in 1:nseg){
	fname	<- paste("estrun_",run_id,"/retcat_explog_seg", ii, "_", ver.date,".rdata",sep="")
	cat(fname, "\n")
	load(fname)
	
	expdata.full	<- rbind(expdata.full, cbind(expdata, segment = ii))
	# exp.est.ls[[ii]]	<- sol
	tmp					<- sapply(fit.lm, coef)
	print(tmp)
	exp.est.ls[[ii]]	<- c(t(tmp[1:2,]), rowMeans(tmp[3:4,]))
	omega.lsf[[ii]]		<- omega.ls
	gamma.ls[[ii]]		<- gamma

	rm(list = setdiff(ls(), c(tmp.ls, "tmp.ls", "ii", "run_id", "fmt_dpt", "hh_dist", "hh_month", "pan", "fmt_name", "R", "seldpt", 
								"dpt_name")))
}

# Summarize the parameters
tmp		<- vector("list", length = nseg)
for(i in 1:length(tmp)){
	tmp[[i]]	<- summary(exp.est.ls[[i]])$estimate
	class(tmp[[i]])	<- "coeftest"
}
stargazer(tmp, type = "text")
stargazer(tmp, type = "html", out = paste("estrun_", run_id,"/retcat_explog_", ver.date,".html", sep=""))

# Plot the conditional expected utility function 
tmp 	<- seq(0, 200, 5)
tmp1	<- c(sapply(gamma.ls, colMeans))
ggtmp	<- data.frame(y = tmp, base = log(tmp + 1), sapply(tmp1, function(x) log(tmp/x + 1)))
colnames(ggtmp)	<- c("y", "base", paste("seg", rep(1:3, each = R), ":",fmt_name, sep=""))
ggtmp	<- melt(ggtmp, id.vars = "y")
ggtmp1	<- data.frame(y = tmp, base = log(tmp + 1), sapply(tmp1, function(x) x*log(tmp/x + 1)))
colnames(ggtmp1)	<- c("y", "base", paste("seg", rep(1:3, each = R), ":",fmt_name, sep=""))
ggtmp1	<- melt(ggtmp1, id.vars = "y")

# ggplot(ggtmp, aes(y, value)) + geom_line() + facet_wrap(~variable)
# ggplot(ggtmp1, aes(y, value)) + geom_line() + facet_wrap(~variable)


# Save the estimation 
save(exp.est.ls, omega.lsf, fmt_name, hh_dist, expdata.full, nseg, pan, R, run_id, fmt_dpt, 
	file = paste("estrun_",run_id,"/retcat_explog_", ver.date,".rdata",sep=""))


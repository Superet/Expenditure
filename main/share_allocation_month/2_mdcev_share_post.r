library(ggplot2)
library(reshape2)
library(maxLik)
library(evd)
library(data.table)
library(stargazer)
library(lmtest)
library(mgcv)

setwd("~/Documents/Research/Store switching/processed data/Estimation")
plot.wd		<- "/Users/chaoqunchen/Desktop"
ww			<- 6.5
ar			<- .6
run_id		<- 9

# Read estimation from all segments #
nseg	<- 3
shr.est	<- setNames(vector("list", nseg), c("Low", "Med", "High"))
exp.est	<- setNames(vector("list", nseg), c("Low", "Med", "High"))
gam.ls	<- setNames(vector("list", nseg), c("Low", "Med", "High"))
shr.par.ls	<- setNames(vector("list", nseg), c("Low", "Med", "High"))
omega.est	<- vector("list", nseg)
ver.date	<- "2016-06-20"
cpi.adj		<- TRUE
tmpls		<- ls()

for(ii in 1:nseg){
	if(cpi.adj){
		loadf	<- paste("estrun_",run_id,"/MDCEV_cpi_est_seg",ii,"_", ver.date,".rdata",sep="")
	}else{
		loadf	<- paste("estrun_",run_id,"/MDCEV_est_seg",ii,"_", ver.date,".rdata",sep="")
	}
	load(loadf)
	shr.est[[ii]]	<- sol
	exp.est[[ii]]	<- sol.top2
	omega.est[[ii]]	<- omega_deriv
	gam.ls[[ii]]	<- gamfit
	shr.par.ls[[ii]]	<- shr.par
	print(ii)
	rm(list = setdiff(ls(), c(tmpls, "ii", "tmpls", "R", "fmt_name", "beta0_base")) )
}

# Eestimation coefficient summary
tmpls	<- setNames(vector("list", length = length(shr.est)), names(shr.est))
for(i in 1:length(shr.est)){
	tmp	<- rbind(coeftest(shr.est[[i]]), coeftest(exp.est[[i]]))
	class(tmp)	<- "coeftest"
	tmpls[[i]]	<- tmp
}
col.lab	<- c("SizeIndex", "ln(UPC per cat)", "ln(Num cat)", "Private label")
tmplab	<- c(paste("log(I)*", fmt_name[-beta0_base],sep=""), "s[t-1]", col.lab, fmt_name[-beta0_base], 
			paste("gamma-", fmt_name, sep=""), "ln(sigma)", "lambda1", "lambda2") 

# Check covariate label 
tmplab
stargazer(tmpls, type = "text", column.labels = names(shr.est))

# Adjust standard error
# se		<- list(NULL, NULL, 
# 			c(.11, .3077, .4473, 1.7990, .0096, .0272, .0391, .1575, .0746, .0192, .0434, .0280, .0587, .2714, .5214, .1539, .1095, .3727, 
# 				2.1367, .0065, .0419, .0037))

# Export estimate table 
stargazer(tmpls, type = "html", title = "Multi-stage model coefficent estimates", 
			covariate.labels = tmplab, column.labels = names(shr.est),
			no.space = TRUE, 
			out = paste(plot.wd, "/overall_est.html", sep=""))

stargazer(tmpls, type = "latex", title = "Multi-stage model coefficent estimates", 
			covariate.labels = tmplab, column.labels = names(shr.est),
			no.space = TRUE,
			out = paste(plot.wd, "/overall_est.tex", sep=""))
stargazer(tmpls, type = "latex", title = "Multi-stage model coefficent estimates", 
			covariate.labels = tmplab, column.labels = names(shr.est), 
			no.space = TRUE)

# Combine the inclusive value plot
pdf(paste(plot.wd, "/graph_omega.pdf", sep=""), width = 7.5, height = 3)
par(mfrow=c(1,3))
for(i in 1:nseg){
	print(vis.gam(gam.ls[[i]], 
			ticktype="detailed", theta=-35,color="bw", 
			main = paste(names(gam.ls)[i], "-income", sep=""), xlab = "log(I)", ylab = "y", zlab="omega"))
}
dev.off()

# Save the parameters to a data set
shr.par.mat	<- do.call(cbind, shr.par.ls)
save(shr.par.mat, file = paste(plot.wd, "/MDCEV_cpiest_shrpar.rdata", sep=""))
			
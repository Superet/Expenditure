library(ggplot2)
library(reshape2)
library(maxLik)
library(evd)
library(data.table)
library(stargazer)
library(lmtest)
library(mgcv)

setwd("/kellogg/users/marketing/2661703/Exercise/run")
run_id		<- 2
make_plot	<- FALSE
source("outreg function.R")

# setwd("~/Documents/Research/Store switching/processed data")
# source("../Exercise/main/outreg function.R")
plot.wd		<- "/Users/chaoqunchen/Desktop"
ww			<- 6.5
ar			<- .6

# Read estimation from all segments #
nseg	<- 3
shr.est	<- setNames(vector("list", nseg), c("Low", "Med", "High"))
exp.est	<- setNames(vector("list", nseg), c("Low", "Med", "High"))
gam.ls	<- setNames(vector("list", nseg), c("Low", "Med", "High"))
omega.est	<- vector("list", nseg)
ver.date	<- "2016-02-26"
cpi.adj		<- TRUE
tmpls		<- ls()

for(ii in 1:nseg){
	# if(cpi.adj){
	# 	loadf	<- paste("estrun_",run_id,"/MDCEV_cpi_est_seg",ii,"_", ver.date,".rdata",sep="")
	# }else{
	# 	loadf	<- paste("estrun_",run_id,"/MDCEV_est_seg",ii,"_", ver.date,".rdata",sep="")
	# }
	if(cpi.adj){
		loadf	<- paste("MDCEV_cpi_est_seg",ii,"_", ver.date,".rdata",sep="")
	}else{
		loadf	<- paste("MDCEV_est_seg",ii,"_", ver.date,".rdata",sep="")
	}
	load(loadf)
	shr.est[[ii]]	<- sol
	exp.est[[ii]]	<- sol.top2
	omega.est[[ii]]	<- omega_deriv
	gam.ls[[ii]]	<- gamfit
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
col.lab	<- c("SizeIndex", "ln(UPC per mod)", "ln(Num mod)", "Private label")
tmplab	<- c(col.lab, paste("log(I)*", col.lab, sep=""), fmt_name[-beta0_base], 
			paste("gamma-", fmt_name, sep=""), "ln(sigma)", "lambda1", "lambda2") 

# Check covariate label 
tmplab
stargazer(tmpls, type = "text", column.labels = names(shr.est))

# Adjust standard error
se		<- list(NULL, NULL, 
			c(.11, .3077, .4473, 1.7990, .0096, .0272, .0391, .1575, .0746, .0192, .0434, .0280, .0587, .2714, .5214, .1539, .1095, .3727, 
				2.1367, .0065, .0419, .0037))

# Export estimate table 
stargazer(tmpls, type = "html", title = "Multi-stage model coefficent estimates", 
			covariate.labels = tmplab, column.labels = names(shr.est),
			no.space = TRUE, se = se, 
			out = paste(plot.wd, "/overall_est.html", sep=""))

stargazer(tmpls, type = "latex", title = "Multi-stage model coefficent estimates", 
			covariate.labels = tmplab, column.labels = names(shr.est),
			no.space = TRUE, se = se, 
			out = paste(plot.wd, "/overall_est.tex", sep=""))
stargazer(tmpls, type = "latex", title = "Multi-stage model coefficent estimates", 
			covariate.labels = tmplab, column.labels = names(shr.est), se = se, 
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

			
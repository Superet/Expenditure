library(ggplot2)
library(reshape2)
library(maxLik)
library(evd)
library(data.table)
library(stargazer)
library(lmtest)
library(mgcv)

# setwd("/sscc/home/c/ccv103/Exercise/run/estrun_4")
# plot.wd		<- getwd()
run_id		<- 8
make_plot	<- FALSE

setwd("~/Documents/Research/Store switching/processed data/Estimation/estrun_8")
plot.wd		<- "/Users/chaoqunchen/Desktop"
ww			<- 6.5
ar			<- .6

# Read estimation from all segments #
nseg	<- 3
shr.est	<- setNames(vector("list", nseg), c("Low", "Med", "High"))
exp.est	<- setNames(vector("list", nseg), c("Low", "Med", "High"))
gam.ls	<- setNames(vector("list", nseg), c("Low", "Med", "High"))
se.ls	<- setNames(vector("list", nseg), c("Low", "Med", "High"))
omega.est	<- vector("list", nseg)
ver.date1	<- "2016-05-25"	#"2016-05-14" #"2016-03-06"
nb1			<- 100
cpi.adj		<- TRUE
tmpls		<- ls()

for(jj in 1:nseg){
	loadf	<- paste("MDCEV_estbt_seg",jj, "_nb", nb1, "_", ver.date1,".rdata",sep="")
	load(loadf)
	shr.est[[jj]]	<- sol
	exp.est[[jj]]	<- sol.top2
	omega.est[[jj]]	<- omega_deriv
	se.ls[[jj]]		<- se.bt
	#gam.ls[[jj]]	<- gamfit
	print(jj)
	rm(list = setdiff(ls(), c(tmpls, "jj", "tmpls", "R", "fmt_name", "beta0_base")) )
}

# Eestimation coefficient summary
tmpls	<- setNames(vector("list", length = length(shr.est)), names(shr.est))
for(i in 1:length(shr.est)){
	tmp	<- rbind(coeftest(shr.est[[i]]), coeftest(exp.est[[i]]))
	class(tmp)	<- "coeftest"
	tmpls[[i]]	<- tmp
}
col.lab	<- c("SizeIndex", "ln(UPC per mod)", "ln(No. category)", "Private label")
tmplab	<- c(paste("log(I)*", fmt_name[-beta0_base], sep=""), col.lab, paste("log(I)*", col.lab, sep=""), fmt_name[-beta0_base], 
			paste("gamma-", fmt_name, sep=""), "ln(sigma)", "lambda1", "lambda2") 

# Check covariate label 
tmplab
stargazer(tmpls, type = "text", se = se.ls, column.labels = names(shr.est))

# Export estimate table 
stargazer(tmpls, type = "html", title = "Multi-stage model coefficent estimates", 
			covariate.labels = tmplab, column.labels = names(shr.est),
			no.space = TRUE, se = se.ls, 
			out = paste(plot.wd, "/overall_estbt.html", sep=""))

stargazer(tmpls, type = "latex", title = "Multi-stage model coefficent estimates", 
			covariate.labels = tmplab, column.labels = names(shr.est),
			no.space = TRUE, se = se.ls, 
			out = paste(plot.wd, "/overall_estbt.tex", sep=""))
stargazer(tmpls, type = "latex", title = "Multi-stage model coefficent estimates", 
			covariate.labels = tmplab, column.labels = names(shr.est), se = se.ls, 
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


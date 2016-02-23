library(ggplot2)
library(reshape2)
library(maxLik)
library(evd)
library(data.table)
library(stargazer)
library(lmtest)

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
omega.est	<- vector("list", nseg)
ver.date	<- "2016-01-22"

for(i in 1:nseg){
	# load(paste("estrun_",run_id,"/MDCEV_est_seg",i,".rdata",sep=""))
	load(paste("Estimation/estrun_", run_id, "/MDCEV_est_seg",i,"_", ver.date, ".rdata",sep=""))
	shr.est[[i]]	<- sol
	# exp.est[[i]]	<- sol1
	exp.est[[i]]	<- sol.top2
	omega.est[[i]]	<- omega_deriv
	print(i)
}

tmpls	<- lapply(shr.est, coeftest)
col.lab	<- c("SizeIndex", "ln(UPC per mod)", "ln(Num mod)", "Private label")
tmplab	<- c(col.lab, paste("ln(income)*", col.lab, sep=""), fmt_name[-beta0_base], 
			paste("gamma-", fmt_name, sep=""), "ln(sigma)")
stargazer(tmpls, type = "html", title = ver.date, 
			covariate.labels = tmplab, 
			out = paste(plot.wd, "/MDCEV_est1.html", sep=""))
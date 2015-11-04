library(ggplot2)
library(reshape2)
library(Rcpp)
library(RcppArmadillo)
library(maxLik)
library(evd)
library(data.table)
library(r2excel)

setwd("/kellogg/users/marketing/2661703/Exercise/run")
run_id		<- 1
make_plot	<- FALSE
source("outreg function.R")

# setwd("~/Documents/Research/Store switching/processed data")
# source("../Exercise/main/outreg function.R")
plot.wd		<- "/Users/chaoqunchen/Desktop"
ww			<- 6.5
ar			<- .6

outfile 	<- paste(plot.wd, "/share_estimate.xlsx", sep="")
mywb		<- createWorkbook()
sht1		<- createSheet(mywb, "Estimates")
make_plot	<- TRUE

# Read estimation from all segments #
nseg	<- 3
shr.est	<- vector("list", nseg)
exp.est	<- vector("list", nseg)
omega.est	<- vector("list", nseg)
for(i in 1:nseg){
	# load(paste("estrun_",run_id,"/MDCEV_est_seg",i,".rdata",sep=""))
	load(paste("Estimation/MDCEV_est_seg",i,".rdata",sep=""))
	shr.est[[i]]	<- sol
	exp.est[[i]]	<- sol1
	names(omega_deriv)	<- lnInc
	omega.est[[i]]	<- omega_deriv
}

#########################################
# Organize table of parameter estimates # 
#########################################
tmp.list<- lapply(1:nseg, function(i) rbind(summary(shr.est[[i]])$estimate, summary(exp.est[[i]])$estimate) )
selcol	<- c("size_index", "ln_upc_per_mod", "ln_num_module","overall_prvt")
beta0_base 	<- which(fmt_name == "Grocery")
varname	<- c(selcol, paste(selcol,"*log(income)", sep=""), paste("Intercept: ", fmt_name[-beta0_base], sep = ""), 
			paste("gamma_", fmt_name, sep=""), "ln_sigma", "lambda0", "lambda1")
tmp.list<- lapply(tmp.list, function(x) data.frame(x, Var = varname))	
tmp.tab	<- model_outreg(tmp.list,p.given=TRUE,digits=4,modelname=c("Low", "Med", "High"),latex.superscript=F,
						head.name=c("Estimate", "Std..error", "Pr...t.", "Var"))
print(tmp.tab)

xlsx.addHeader(mywb, sht1, value = "Parameter estimates from three income groups", level = 2)						
xlsx.addTable(mywb, sht1, tmp.tab, row.names = FALSE)

########################
# Plot inclusive value #
########################
exp(lnInc)
selinc	<- c(9000, 42500, 85000)
sel 	<- which(round(exp(lnInc), 0) %in% selinc)
tmp		<- seq(10, 200, 1)
ggtmp 	<- data.frame(NULL)
for(i in 1:nseg){
	for(j in 1:length(selinc)){
		ggtmp	<- rbind(ggtmp, data.frame(IncGrp = i, Income = selinc[j], Exp = tmp, omega = omega.est[[i]][[sel[j]]](tmp) ))
	}
}
ggtmp$IncGrp	<- factor(ggtmp$IncGrp, levels = 1:nseg, labels = c("Low", "Med", "High"))
ggtmp$Income 	<- factor(ggtmp$Income)

if(make_plot){
	pdf(paste(plot.wd, "/graph_omega.pdf", sep=""), width = ww, height = ww*ar)
	plots	<- ggplot(ggtmp, aes(Exp, omega, linetype = Income)) + geom_line() + 
					facet_grid(. ~ IncGrp) + 
					xlab("Grocery expenditure") + ylab(expression(omega(Inc, y)))  
	print(plots)
	dev.off()
}


saveWorkbook(mywb, outfile)


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

setwd("U:/Users/ccv103/Documents/Research/Store switching/run")
plot.wd 	<- get.wd()

ww			<- 6.5
ar			<- .6
run_id		<- 10

# Read estimation from all segments #
nseg	<- 3
shr.est	<- setNames(vector("list", nseg), c("Low", "Med", "High"))
exp.est	<- setNames(vector("list", nseg), c("Low", "Med", "High"))
gam.ls	<- setNames(vector("list", nseg), c("Low", "Med", "High"))
shr.par.ls	<- setNames(vector("list", nseg), c("Low", "Med", "High"))
omega.est	<- vector("list", nseg)
ver.date	<- "2016-06-30"
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
	# gam.ls[[ii]]	<- gamfit
	shr.par.ls[[ii]]	<- shr.par
	print(ii)
	rm(list = setdiff(ls(), c(tmpls, "ii", "tmpls", "R", "fmt_name", "beta0_base")) )
}
plot.wd	<- paste(getwd(), "/estrun_", run_id, sep="")

# Eestimation coefficient summary
tmpls	<- setNames(vector("list", length = length(shr.est)), names(shr.est))
for(i in 1:length(shr.est)){
	tmp	<- rbind(coeftest(shr.est[[i]]), coeftest(exp.est[[i]]))
	class(tmp)	<- "coeftest"
	tmpls[[i]]	<- tmp
}
col.lab	<- c("Access", "SizeIndex", "ln(UPC per cat)", "ln(No. cat)", "Private label")
# tmplab	<- c(paste("log(I)*", fmt_name[-beta0_base],sep=""), "s[t-1]", col.lab, fmt_name[-beta0_base], 
# 			paste("gamma-", fmt_name, sep=""), "ln(sigma)", "lambda1", "lambda2") 
tmplab	<- c(paste("rec*", fmt_name[-beta0_base],sep=""), col.lab, paste("rec*", col.lab, sep=""), fmt_name[-beta0_base], 
			paste("gamma-", fmt_name, sep=""), "ln(sigma)", "lambda1", "lambda2")
names(tmplab)	<- c(paste("beta_",1:15, sep=""), paste("beta0_",setdiff(1:R, beta0_base), sep=""), 
					paste("gamma_",1:R, sep=""), "ln_sigma", "lambda1", "lambda2")
			
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
			
# Export a few parameters
(selpar	<- paste("beta_", 12:15, sep=""))
stargazer(tmpls, type = "html", title = "Multi-stage model coefficent estimates", 
			covariate.labels = tmplab[selpar], column.labels = c("Low-income", "Medium-income", "High-income"),
			no.space = TRUE, keep = selpar, 
			out = paste(plot.wd, "/sub_est.html", sep=""))			

# Plot the parameters
ctv				<- qnorm(.95)
ggtmp			<- as.matrix(do.call(rbind, lapply(tmpls, function(x) x[selpar,])))
ggtmp			<- data.frame(ggtmp, Parameter = tmplab[rownames(ggtmp)], 
					IncGrp = rep(c("Low", "Medium", "High"), each = length(selpar)))
ggtmp$IncGrp	<- factor(ggtmp$IncGrp, levels = c("High", "Medium", "Low"), ordered = TRUE)					
ggtmp$Parameter	<- factor(ggtmp$Parameter, levels = c("rec*ln(UPC per cat)","rec*Private label", "rec*SizeIndex", "rec*ln(No. cat)"))

pdf(paste(plot.wd,"/slideg_est_par.pdf", sep=""), width = ww, height = ww*ar)
ggplot(ggtmp, aes(Parameter, Estimate, color = IncGrp)) + 
		geom_pointrange(aes(y = Estimate, ymin = Estimate - ctv* Std..Error, ymax = Estimate + ctv* Std..Error, color = IncGrp), 
				position = position_dodge(width = 0.50)) + 
		geom_hline(yintercept = 0, size = .25, linetype = 2) + 
		coord_flip() + 
		scale_color_grey(name="Income\ngroup", start = 0, end = .8) +  
		guides(color = guide_legend(reverse = TRUE)) + theme_bw()
dev.off()


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
save(shr.par.mat, file = paste("estrun_",run_id, "/MDCEV_cpiest_shrpar.rdata", sep=""))
			
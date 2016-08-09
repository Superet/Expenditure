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
ver.date	<- "2016-07-07"
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
col.lab	<- c("Package size", "Assortment depth", "Assortment width", "Private label")
# tmplab	<- c(paste("log(I)*", fmt_name[-beta0_base],sep=""), "s[t-1]", col.lab, fmt_name[-beta0_base], 
# 			paste("gamma-", fmt_name, sep=""), "ln(sigma)", "lambda1", "lambda2") 
# tmplab	<- c(paste("rec*", fmt_name[-beta0_base],sep=""), "s_{t-1}", col.lab, paste("rec*", col.lab, sep=""), fmt_name[-beta0_base], 
# 			paste("gamma-", fmt_name, sep=""), "ln(sigma)", "lambda1", "lambda2")
tmplab	<- c(paste("rec*", fmt_name[-beta0_base],sep=""), col.lab, paste("rec*", col.lab, sep=""), fmt_name[-beta0_base], 
			paste("gamma-", fmt_name, sep=""), "ln(sigma)", "lambda1", "lambda2")
names(tmplab)	<- c(paste("beta_",1:13, sep=""), paste("beta0_",setdiff(1:R, beta0_base), sep=""), 
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
			out = paste(plot.wd, "/overall_est_", ver.date, ".html", sep=""))

stargazer(tmpls, type = "latex", title = "Multi-stage model coefficent estimates", 
			covariate.labels = tmplab, column.labels = names(shr.est),
			no.space = TRUE,
			out = paste(plot.wd, "/overall_est_", ver.date, ".tex", sep=""))
stargazer(tmpls, type = "latex", title = "Multi-stage model coefficent estimates", 
			covariate.labels = tmplab, column.labels = names(shr.est), 
			no.space = TRUE)
			
# Export a few parameters
(selpar	<- paste("beta_", 11:14, sep=""))
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

# Plot the parameter comparison
ctv				<- qnorm(.95)
sel0			<- paste("beta_", 7:10,sep="")
sel1			<- paste("beta_", 11:14,sep="")
tmp0			<- as.matrix(do.call(rbind, lapply(tmpls, function(x) x[sel0,])))
tmp1			<- as.matrix(do.call(rbind, lapply(tmpls, function(x) x[sel1,])))
tmp1			<- tmp0 + tmp1
tmpv			<- lapply(shr.est, function(x) vcov(x))
se0				<- as.matrix(do.call(rbind, lapply(tmpv, function(x) sqrt(diag(x[sel0,sel0])))))
se1				<- lapply(tmpv, function(x) sapply(1:length(sel0), function(i) sum(x[c(sel0[i], sel1[i]), c(sel0[i], sel1[i])])))
se1				<- as.matrix(do.call(rbind, lapply(se1, sqrt)))
ggtmp			<- rbind(data.frame(par = rownames(tmp0), Estimate = tmp0[,"Estimate"], se = c(t(se0)), 
									IncGrp = rep(c("Low", "Medium", "High"), each = length(sel0)), rec = 0), 
						data.frame(par = rownames(tmp0), Estimate = tmp1[,"Estimate"], se = c(t(se1)), 
									IncGrp = rep(c("Low", "Medium", "High"), each = length(sel0)), rec = 1)
						)
ggtmp$IncGrp	<- factor(ggtmp$IncGrp, levels = c("Low", "Medium", "High"), ordered = TRUE)					
ggtmp$Parameter	<- factor(ggtmp$par, levels = sel0, labels = tmplab[sel0])
ggtmp$Recession	<- ifelse(ggtmp$rec == 0, "Before 2008", "After 2008")

selpar		<- tmplab[sel0]
ggplot(subset(ggtmp, rec == 0 & Parameter %in% selpar), aes(IncGrp, Estimate)) + 
		geom_pointrange(aes(y = Estimate, ymin = Estimate - ctv* se, ymax = Estimate + ctv* se), 
				position = position_dodge(width = 0.50)) + 
		geom_hline(yintercept = 0, size = .25, linetype = 2) + 
		# coord_flip() + 
		facet_wrap(~Parameter, scales = "free") +
		theme_bw()
	
selpar	<- "Assortment width"	
pdf(paste(plot.wd, "/slideg_est", ver.date, ".pdf", sep=""), width = 4.5, height = 4.5*ar)
ggplot(subset(ggtmp, Parameter %in% selpar), aes(IncGrp, Estimate, color = Recession)) + 
		geom_pointrange(aes(y = Estimate, ymin = Estimate - ctv* se, ymax = Estimate + ctv* se), 
				position = position_dodge(width = 0.50)) + 
		geom_hline(yintercept = 0, size = .25, linetype = 2) + 
		# coord_flip() + 
		scale_color_grey(name="", start = 0, end = .8) +  
		facet_wrap(~Parameter, scales = "free") +
		labs(x = "Income group") + 
		guides(color = guide_legend(reverse = TRUE)) + theme_bw()
dev.off()

####################################
# Combine the inclusive value plot #
####################################
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
			
library(ggplot2)
library(reshape2)
library(Rcpp)
library(RcppArmadillo)
library(maxLik)
library(evd)
library(data.table)
library(scales)

setwd("~/Documents/Research/Store switching/processed data")
plot.wd		<- "/Users/chaoqunchen/Desktop"
ww			<- 6.5
ar			<- .8

outfile 	<- paste(plot.wd, "/counter1.xlsx", sep="")
mywb		<- createWorkbook()
sht1		<- createSheet(mywb, "Counterfactual1")
make_plot	<- TRUE

#############
# Functions #
#############
weight_summary	<- function(x, freq, fun, ...){
	x.full	<- unlist(lapply(1:length(x), function(i) rep(x[i], freq[i])))
	out		<- fun(x.full, ...)
	return(out)
}

mean_qt <- function(x, alpha = .25){
	out	<- c(mean(x), quantile(x, c(alpha, 1-alpha)))
	names(out)	<- c("y", "ymin", "ymax")
	return(out)
}

###########################
# Read counterfactal data # 
###########################
# Read estimation from all segments #
nseg	<- 3
sim.ls	<- setNames(vector("list", length = nseg), c("Low", "Med", "High"))
my.simdata	<- data.frame()
ver.date1	<- "2015-12-12"
tmpls	<- ls()

for(ii in 1:nseg){
	load(paste("Estimation/ctrfact_seg",ii, "_", ver.date1, ".rdata",sep=""))
	inc.tab		<- table(sim.data$income)
	tmp			<- rbind(sim.X.ls, cbind(sim.prc.ls, Var = "price"))
	tmp$income2007	<- exp(tmp$lnInc)/(1 - .1)
	tmp$nhh		<- inc.tab[as.character(tmp$income2007)]
	sel			<- gsub("\\s", ".", fmt_name)
	
	# Convert expenditure share to expenditure dollars for each retailer
	tmp[,sel]	<- tmp[,sel]*tmp$Exp
	sim.base07[,fmt_name]	<- sim.base07[,fmt_name] * sim.base07$Exp
	sim.base08[,fmt_name]	<- sim.base08[,fmt_name] * sim.base08$Exp
	
	sim.ls[[ii]]<- list(	Base07 	= cbind(sim.base07, sim.unq, nhh = as.numeric(inc.tab)), 
							Base08 	= cbind(sim.base08, sim.unq, nhh = as.numeric(inc.tab)), 
							Retail	= tmp)
	my.simdata	<- rbind(my.simdata, data.frame(sim.data, IncGrp = names(sim.ls)[ii]))
	rm(list = setdiff(ls(), c(tmpls, "tmpls", "ii", "fmt_name", "R")))
	print(ii)						
}

#################
# Income effect #
#################
# Average income effect
tmp		<- do.call(rbind, lapply(sim.ls, function(x) x$Base07))
tmp1	<- apply(tmp[,c("Exp", fmt_name)], 2, function(x) weight_summary(x, tmp$nhh, mean))
tmp		<- do.call(rbind, lapply(sim.ls, function(x) x$Base08))
tmp2	<- apply(tmp[,c("Exp", fmt_name)], 2, function(x) weight_summary(x, tmp$nhh, mean))
tmp.tab	<- data.frame(Avg2007 = tmp1, Avg2008 = tmp2, Difference = tmp2-tmp1)
cat("Averge income effect:\n"); print(tmp.tab); cat("\n")

#-----------------------------#
# Heterogeneous income effect #
ggtmp	<- data.frame()
for(i in 1:nseg){
	tmp		<- apply(sim.ls[[i]]$Base07[,c("Exp", fmt_name)], 2, function(x) weight_summary(x, sim.ls[[i]]$Base07$nhh, mean_qt))
	tmp1	<- data.frame(t(tmp), year = 2007, IncGrp = names(sim.ls)[i], Var = colnames(tmp))
	tmp		<- apply(sim.ls[[i]]$Base08[,c("Exp", fmt_name)], 2, function(x) weight_summary(x, sim.ls[[i]]$Base08$nhh, mean_qt))
	tmp2	<- data.frame(t(tmp), year = 2008, IncGrp = names(sim.ls)[i], Var = colnames(tmp))
	rownames(tmp1)	<- rownames(tmp2) <- NULL
	ggtmp	<- rbind(ggtmp, tmp1, tmp2)
}

ggplot(subset(ggtmp, Var == "Exp"), aes(IncGrp, y, col = as.factor(year), fill = as.factor(year))) + 
	geom_bar(stat = "identity", position = position_dodge(width=0.8), width = .5) + 
	geom_errorbar(aes(ymin = ymin, ymax = ymax), position = position_dodge(width=0.8), width = .2)

ggplot(ggtmp, aes(IncGrp, y, col = as.factor(year))) + 
	geom_pointrange(aes(y = y, ymin = ymin, ymax = ymax), position = position_dodge(width=0.8)) + 
	facet_wrap(~Var, scales = "free")

####################
# Retail marketing # 
####################
ggtmp	<- data.frame()	
nvar	<- length(unique(sim.ls[[1]]$Retail$Var))
for(i in 1:R){
	tmpall		<- data.frame()
	for(j in 1:nseg){
		sel		<- sim.ls[[j]]$Retail$retailer == fmt_name[i]
		tmp1	<- sim.ls[[j]]$Retail[sel,]
		sel		<- gsub("\\s", ".", fmt_name[i])
		tmp1[,sel] <- tmp1[,sel] - rep(sim.ls[[j]]$Base08[,fmt_name[i]], nvar)
		tmpall	<- rbind(tmpall, tmp1[,c(sel,"retailer", "Var", "income2007", "nhh")])
		tmp		<- tapply(tmp1[,sel], tmp1$Var, function(x) weight_summary(x, sim.ls[[j]]$Base08[,"nhh"], mean_qt))
		tmp		<- data.frame(do.call(rbind, tmp), Var = names(tmp), retail = fmt_name[i], IncGrp = names(sim.ls)[j])
		ggtmp	<- rbind(ggtmp, tmp)
	}
	
	# Outcome difference across all segments
	tmp	<- split(tmp, tmp$Var)
	tmp	<- sapply(tmp, function(x) weight_summary(x[,sel], x$nhh, mean_qt))
	tmp		<- data.frame(do.call(rbind, tmp), Var = names(tmp), retail = fmt_name[i], IncGrp = "Overall")
	ggtmp	<- rbind(ggtmp, tmp)
}

ggplot(ggtmp, aes(Var, y, col = IncGrp)) + geom_pointrange(aes(y = y, ymin=ymin, ymax=ymax)) + 
		facet_wrap(~retail, scales = "free")





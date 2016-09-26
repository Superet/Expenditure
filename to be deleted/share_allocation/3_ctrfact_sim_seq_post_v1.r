library(ggplot2)
library(reshape2)
library(Rcpp)
library(RcppArmadillo)
library(maxLik)
library(evd)
library(data.table)
library(scales)
library(gridExtra)
library(stargazer)

plot.wd		<- "/Users/chaoqunchen/Desktop"

# setwd("/home/brgordon/ccv103/Exercise/run")
run_id		<- 8
ww			<- 6.5
ar			<- .6
ww1			<- 8.5
ar1			<- .5
setwd(paste("~/Documents/Research/Store switching/processed data/Estimation/estrun_", run_id, sep=""))

# outfile 	<- paste(plot.wd, "/counter1.xlsx", sep="")
# mywb		<- createWorkbook()
# sht1		<- createSheet(mywb, "Counterfactual1")
# make_plot	<- TRUE

# Set plotting parameters
annual.week	<- 26
annual		<- TRUE
draw.par		<- FALSE
sim.omega		<- FALSE

###########################
# Read counterfactal data # 
###########################
# Read estimation from all segments #
nseg	<- 3
segname	<- c("Low", "Med", "High")
sim.ls	<- setNames(vector("list", length = nseg), c("Low", "Med", "High"))
my.simdata	<- data.frame()
sim.df		<- data.frame()
numsim		<- 1000
ver.date1	<- "2016-05-13"  #"2016-05-02" 
tmpls	<- ls()

for(ii in 1:nseg){
	# Load data 
	# load(paste("Estimation/estrun_2/ctrfact_seg",ii, "_", ver.date1, ".rdata",sep=""))
	fname			<- paste("ctrfact_seq_seg",ii, sep="")
	if(draw.par){ fname <- paste(fname, "_pardraw", sep="") }
	if(sim.omega){fname	<- paste(fname, "_simomega", sep="") }
	fname			<- paste(fname, "_sim", numsim, "_", ver.date1, sep="")
	load(paste(fname, ".rdata",sep=""))
	
	# Compute the shopping incidence for baseline simulation
	tmp1	<- apply(sim.base07$Allocation, c(2,3), function(x) mean(x>.01))
	colnames(tmp1)	<- paste("p_", fmt_name, sep="")
	tmp2	<- data.frame(lnInc = as.numeric(rownames(tmp1)), sim.base07$Average, tmp1)
	names(tmp2)	<- gsub("\\.", " ", names(tmp2))
	tmp1	<- apply(sim.base08$Allocation, c(2,3), function(x) mean(x>.01))
	colnames(tmp1)	<- paste("p_", fmt_name, sep="")
	tmp3	<- data.frame(lnInc = as.numeric(rownames(tmp1)), sim.base08$Average, tmp1)
	names(tmp3)	<- gsub("\\.", " ", names(tmp3))
		
	# Compute the shopping incidence for price simualtion in 2007
	sim07			<- do.call(rbind, lapply(sim.prc.ls07, function(x) x$Average))
	names(sim07)	<- gsub("\\.", " ", names(sim07))
	tmp1	<- melt(lapply(sim.prc.ls07, function(x) apply(x$Allocation, c(1,3,4), function(xx) mean(xx>0.01))) )
	tmp1	<- dcast(tmp1, L1+change+lnInc ~ retail, value.var = "value")
	names(tmp1)	<- c("retailer", "change", "lnInc", paste("p_", names(tmp1)[-(1:3)], sep=""))
	tmp1$retailer	<- fmt_name[tmp1$retailer]
	sim07		<- merge(sim07, tmp1, by = c("retailer", "change", "lnInc"))
	
	# Compute the shopping incidence for price simualtion in 2007
	sim08			<- do.call(rbind, lapply(sim.prc.ls08, function(x) x$Average))
	names(sim08)	<- gsub("\\.", " ", names(sim08))
	tmp1	<- melt(lapply(sim.prc.ls08, function(x) apply(x$Allocation, c(1,3,4), function(xx) mean(xx>0.01))) )
	tmp1	<- dcast(tmp1, L1+change+lnInc ~ retail, value.var = "value")
	names(tmp1)	<- c("retailer", "change", "lnInc", paste("p_", names(tmp1)[-(1:3)], sep=""))
	tmp1$retailer	<- fmt_name[tmp1$retailer]
	sim08		<- merge(sim08, tmp1, by = c("retailer", "change", "lnInc"))
				
	# If we want to convert the biweekly simulation to annual level
	if(annual){
		sel	<- c("Exp", fmt_name)
		tmp2[,sel]	<- tmp2[,sel] * annual.week
		tmp3[,sel]	<- tmp3[,sel] * annual.week
		sim07[,sel]	<- sim07[,sel] * annual.week
		sim08[,sel]	<- sim08[,sel] * annual.week
	}
		
	sim.ls[[ii]]<- list(Base07 = tmp2, Base08 = tmp3, price07 = sim07, price08 = sim08)
						
	my.simdata	<- rbind(my.simdata, data.frame(sim.data, IncGrp = names(sim.ls)[ii]))
	
	# Delete unused objects
	rm(list = setdiff(ls(), c(tmpls, "tmpls", "ii", "fmt_name", "R", "iv.change.vec")))
	print(ii)						
}

#############
# Functions #
#############
weight_summary	<- function(x, freq, fun, ...){
	x.full	<- unlist(lapply(1:length(x), function(i) rep(x[i], freq[i])))
	out		<- fun(x.full, ...)
	return(out)
}

mean_qt <- function(x, alpha = .025){
	out	<- c(mean(x), quantile(x, c(alpha, 1-alpha)))
	names(out)	<- c("y", "ymin", "ymax")
	return(out)
}

weight_se	<- function(se, freq){
# This function computes the standard error of x.bar = sum_i (freq_i*x_i)/sum_i freq_i
	wt	<- freq/sum(freq)
	v	<- sum(wt^2 * se^2, na.rm =T)
	return(sqrt(v))
}

mytrimfun	<- function (x, alpha = 0.05){
    if (alpha == 0) {
        return(x)
    }else {
        n <- length(x)
        lo <- floor(n * alpha) + 1
        hi <- n + 1 - lo
		if(sum(is.na(x)) >= (lo-1)){
			out	<- x
		}else{
			out <- sort.int(x, partial = unique(c(lo, hi)), na.last = FALSE)[lo:hi]
		}      
        return(out)
    }
}

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

compare_fn <- function(dat1, dat2, byvar, weight.dat, weight.by, fml, val1 = "value", val2 = "value", agg.first = FALSE, y.name = "y"){
	if(nrow(dat1)!=nrow(dat2)){
		stop("Error: nrow(dat1) not equal to nrow(dat2).")
	}
	if(!identical(dat1[,byvar], dat2[, byvar])){
		stop("Error: two data are not sorted equally.")
	}
	if(agg.first){
		tmp	<- paste("list(", paste(byvar, collapse=","), ")")
		x1	<- data.table(merge(dat1, weight.dat, by.var = weight.by, all.x = T))
		x1	<- x1[,list(x1 = weight_summary(get(val1), nhh, mean)), by = eval(parse(text = tmp))]
		x2	<- data.table(merge(dat2, weight.dat, by.var = weight.by, all.x = T))
		x2	<- x2[,list(x2 = weight_summary(get(val2), nhh, mean)), by = eval(parse(text = tmp))]
		out	<- merge(x1, x2, by = byvar)
		out	<- out[, y:=eval(parse(text = fml))]
		out	<- data.frame(out)
		out	<- out[,c(byvar, "y")]
		if(y.name != "y"){
			names(out)	<- gsub("y", y.name, names(out))
		}
	}else{
		x1	<- dat1[,val1]
		x2	<- dat2[,val2]
		y	<- eval(parse(text = fml))
		dat.sum 	<- data.frame(dat1[,unique(c(byvar, weight.by))], y = y)
		# names(dat.sum)	<- c(byvar, "y")
		dat.sum		<- merge(dat.sum, weight.dat, by.var = weight.by, all.x=T)
		dat.sum		<- data.table(dat.sum)
		tmp			<- paste("list(", paste(byvar, collapse=","), ")")
		out			<- dat.sum[,list(y = weight_summary(y, nhh, mean), y.sd = weight_summary(y, nhh, sd)), 
								by = eval(parse(text = tmp))]
		out			<- data.frame(out)
		if(y.name != "y"){
			names(out)	<- gsub("y", y.name, names(out))
		}
	}
	return(out)
}

#################
# Income effect #
#################
iv.change 	<- -.1
p.idx		<- 0
plots		<- list(NULL)
test.alpha	<- .95						# Significance level
crtv		<- qnorm(test.alpha)		# Critical value
ret.ord		<- c("Dollar Store", "Grocery", "Discount Store", "Drug Store", "Convenience Store", "Warehouse Club")
dodge.width	<- .5						# The width parameter in position_dodge

# Income frequency
tab.inc	<- data.table(my.simdata)
tab.inc	<- tab.inc[, list(nhh = length(unique(household_code))), by = list(IncGrp, income)]
setnames(tab.inc, "income", "income2007")
setkeyv(tab.inc, c("IncGrp", "income2007"))
tab.inc$income2007	<- round(tab.inc$income2007, 2)

# Merge the weights into the baseline simulation
base07	<- do.call(rbind, lapply(1:nseg, function(i) cbind(IncGrp=names(sim.ls)[i], sim.ls[[i]]$Base07)))
base07$income2007	<- round(exp(base07$lnInc), 2)
base07	<- merge(base07, tab.inc, by = c("IncGrp", "income2007"), all.x = T)
base08	<- do.call(rbind, lapply(1:nseg, function(i) cbind(IncGrp=names(sim.ls)[i], sim.ls[[i]]$Base08)))
base08$income2007	<- round(exp(base08$lnInc)/(1+iv.change), 2)
base08	<- merge(base08, tab.inc, by = c("IncGrp", "income2007"), all.x = T)

# Check if data is ordered 
base07	<- base07[order(base07$IncGrp, base07$income2007),]
base08	<- base08[order(base08$IncGrp, base08$income2007),]
identical(base07$income2007, base08$income2007)

# ---------------------------------------- #
# Income effect on the overall expenditure #
sel	<- "Exp"
tmp	<- data.table(base07[,c("IncGrp", "nhh")], base2007 = base07[,sel], base2008 = base08[,sel], 
				Dif = base08[,sel] - base07[,sel], Elasticity = (base08[,sel]-base07[,sel])/base07[,sel]/iv.change, 
				PTC = (base08[,sel]-base07[,sel])/(base07[,sel]*iv.change))
tmp1<- tmp[, list(base2007 = weight_summary(base2007, nhh, mean), base2008 = weight_summary(base2008, nhh, mean), 
				Dif = weight_summary(Dif, nhh, mean), Dif.sd = weight_summary(Dif, nhh, sd), 
				Elasticity = weight_summary(Elasticity, nhh, mean), Elasticity.sd = weight_summary(Elasticity, nhh, sd), 
				PTC = weight_summary(PTC, nhh, mean), PTC.sd = weight_summary(PTC, nhh, sd))]
tmp2<- tmp[, list(base2007 = weight_summary(base2007, nhh, mean), base2008 = weight_summary(base2008, nhh, mean), 
				Dif = weight_summary(Dif, nhh, mean), Dif.sd = weight_summary(Dif, nhh, sd), 
				Elasticity = weight_summary(Elasticity, nhh, mean), Elasticity.sd = weight_summary(Elasticity, nhh, sd), 
				PTC = weight_summary(PTC, nhh, mean), PTC.sd = weight_summary(PTC, nhh, sd)), by = list(IncGrp)]
cat("Income effect on overall expenditure:\n"); print(tmp1); cat(".\n")
cat("Income effect on overall expenditure by income group:\n"); print(tmp2); cat(".\n")
test.alpha	<- .9
cat("The", test.alpha*100, "% CI. of elasticity for each income group is:\n"); 
	cbind(tmp2$Elasticity-qnorm(test.alpha)*tmp2$Elasticity.sd, tmp2$Elasticity+qnorm(test.alpha)*tmp2$Elasticity.sd); cat("\n")
	
# -----------------------------------#
# Income effect on expenditure share #
sel 	<- fmt_name
tmp		<- base08[,sel]/base08[,"Exp"] - base07[,sel]/base07[,"Exp"]		# Change in expenditure share
ggtmp	<- melt(data.frame(base07[,c("IncGrp", "nhh")], tmp), id.var = c("IncGrp", "nhh"))
tmp1	<- setNames(fmt_name, gsub(" ", "\\.", fmt_name))
ggtmp$retailer	<- tmp1[as.character(ggtmp$variable)]
ggtmp$IncGrp	<- factor(ggtmp$IncGrp, levels = segname)
ggtmp$IncGrp	<- factor(ggtmp$IncGrp, levels = rev(levels(ggtmp$IncGrp)), ordered = TRUE)
ggtmp	<- data.table(ggtmp)	

ggtmp1	<- ggtmp[, list(change = weight_summary(value, nhh, mean), change.sd = weight_summary(value, nhh, sd)), by = list(retailer)]
ggtmp2	<- ggtmp[, list(change = weight_summary(value, nhh, mean), change.sd = weight_summary(value, nhh, sd)), by = list(IncGrp, retailer)]
# (ret.ord	<- ggtmp1[order(ggtmp1$change),retailer])
ggtmp1$retailer	<- factor(ggtmp1$retailer, levels = ret.ord)
ggtmp2$retailer	<- factor(ggtmp2$retailer, levels = ret.ord)
mylim	<- range(c(ggtmp2$change-crtv*ggtmp2$change.sd, ggtmp2$change+crtv*ggtmp2$change.sd))

# Plot the change in expenditure share of 10% income drop
plots	<- list(NULL)
plots[[1]]	<- ggplot(ggtmp1, aes(retailer, change)) + 
			geom_pointrange(aes(y = change, ymin = change-crtv*change.sd, ymax = change + crtv * change.sd)) + 
			geom_hline(yintercept = 0, size = .25, linetype = 2) + 
			scale_y_continuous(labels = percent, limits = mylim) + 
			xlab("Retail format") + ylab("") + coord_flip() + theme_bw()
plots[[2]]	<- ggplot(ggtmp2, aes(retailer, change, col = IncGrp)) + 
			geom_pointrange(aes(ymin = change-crtv*change.sd, ymax = change + crtv * change.sd, col = IncGrp), 
					position=position_dodge(width=dodge.width)) +
			geom_hline(yintercept = 0, size = .25, linetype = 2) + 		
			scale_color_grey(name="Income\ngroup", start = 0, end = .8) +  
			scale_y_continuous(labels = percent, limits = mylim) + 
			xlab("Retail format") + ylab("") + coord_flip() + 
			guides(color = guide_legend(reverse = TRUE)) + theme_bw()
plots[[3]]	<- get_legend(plots[[2]])
plots[[2]]	<- plots[[2]] + theme(legend.position="none")	

pdf(paste(plot.wd, "/graph_ctr_income_effect_share.pdf", sep=""), width = ww1, height = ww1*ar1)
grid.arrange(plots[[1]], plots[[2]], plots[[3]], nrow = 1, widths = c(.45, .45, .1), 
	# top = "Change of expenditure share for 10% income drop", bottom = "Change of expenditure share")
	bottom = "Change of expenditure share")
dev.off()
	 				
# ------------------------------------#
# Income effect on shopping incidence #
sel 	<- paste("p_", fmt_name, sep="")
tmp		<- base08[,sel] - base07[,sel]		# Change in shopping probability
ggtmp	<- melt(data.frame(base07[,c("IncGrp", "nhh")], tmp), id.var = c("IncGrp", "nhh"))
tmp1	<- setNames(fmt_name, gsub(" ", "\\.", sel))
ggtmp$retailer	<- tmp1[as.character(ggtmp$variable)]
ggtmp$retailer	<- factor(ggtmp$retailer, levels = ret.ord)
ggtmp$IncGrp	<- factor(ggtmp$IncGrp, levels = segname)
ggtmp$IncGrp	<- factor(ggtmp$IncGrp, levels = rev(levels(ggtmp$IncGrp)), ordered = TRUE)
ggtmp	<- data.table(ggtmp)	

ggtmp1	<- ggtmp[, list(change = weight_summary(value, nhh, mean), change.sd = weight_summary(value, nhh, sd)), by = list(retailer)]
ggtmp2	<- ggtmp[, list(change = weight_summary(value, nhh, mean), change.sd = weight_summary(value, nhh, sd)), by = list(IncGrp, retailer)]
cat("On average, the probability of shopping incidence at retail formats:\n"); print(ggtmp1); cat("\n")

# Plot the change of shopping incidence 
mylim	<- range(c(ggtmp2$change-crtv*ggtmp2$change.sd, ggtmp2$change+crtv*ggtmp2$change.sd))
plots	<- list(NULL)
plots[[1]]	<- ggplot(ggtmp1, aes(retailer, change)) + 
			geom_pointrange(aes(y = change, ymin = change-crtv*change.sd, ymax = change + crtv * change.sd)) + 
			geom_hline(yintercept = 0, size = .25, linetype = 2) + 
			ylim(mylim) + 
			xlab("Retail format") + ylab("") + coord_flip() + theme_bw()
plots[[2]]	<- ggplot(ggtmp2, aes(retailer, change, col = IncGrp)) + 
			geom_pointrange(aes(ymin = change-crtv*change.sd, ymax = change + crtv * change.sd, col = IncGrp), 
					position=position_dodge(width=0.6)) +
			geom_hline(yintercept = 0, size = .25, linetype = 2) + 
			scale_color_grey(name="Income\ngroup", start = 0, end = .6) +  
			ylim(mylim) + 
			xlab("Retail format") + ylab("") + coord_flip() + 
			guides(color = guide_legend(reverse = TRUE)) + theme_bw()
plots[[3]]	<- get_legend(plots[[2]])
plots[[2]]	<- plots[[2]] + theme(legend.position="none")	

pdf(paste(plot.wd, "/graph_ctr_income_effect_incidence.pdf", sep=""), width = ww1, height = ww1*ar1)
grid.arrange(plots[[1]], plots[[2]], plots[[3]], nrow = 1, widths = c(.45, .45, .1), 
	 bottom = "Change in incidence probability")
dev.off()

#############################################
# Organize summary table for counterfactual # 
#############################################
# Break-even margin is the original margin m0 such that r1/r0 > m0/m1, for delta = p1/p0
# m0 > (1-delta)*(r1/r0)/(r1/r0 - delta)

delta 	<- .9

# --------------#
# Income effect #
# Expenditure ratio
tmp07	<- melt(cbind(base07[,c("IncGrp", "income2007", fmt_name)]), id.vars = c("IncGrp", "income2007"), value.name = "exp")
tmp08	<- melt(cbind(base08[,c("IncGrp", "income2007", fmt_name)]), id.vars = c("IncGrp", "income2007"), value.name = "exp")
identical(tmp07[,c("IncGrp", "income2007", "variable")], tmp08[,c("IncGrp", "income2007", "variable")])
tmp1	<- compare_fn(tmp07, tmp08, byvar = "variable", weight.dat = tab.inc, weight.by = c("IncGrp", "income2007"), 
						fml = "x2/x1", val1 = "exp", val2 = "exp", agg.first = FALSE)

# Share difference
dat1	<- melt(cbind(base07[,c("IncGrp", "income2007")], base07[,fmt_name]/base07[, "Exp"]), id.vars = c("IncGrp", "income2007"), value.name = "share")
dat2	<- melt(cbind(base08[,c("IncGrp", "income2007")], base08[,fmt_name]/base08[, "Exp"]), id.vars = c("IncGrp", "income2007"), value.name = "share")
identical(dat1[,c("IncGrp", "income2007", "variable")], dat2[,c("IncGrp", "income2007", "variable")])
tmp2	<- compare_fn(dat1, dat2, byvar = "variable", weight.dat = tab.inc, weight.by = c("IncGrp", "income2007"), 
						fml = "x2-x1", val1 = "share", val2 = "share", agg.first = FALSE)

# Shopping incidence
sel		<- paste("p_", fmt_name, sep="")
dat1	<- melt(cbind(base07[,c("IncGrp", "income2007")], base07[,sel]), id.vars = c("IncGrp", "income2007"), value.name = "p")
dat2	<- melt(cbind(base08[,c("IncGrp", "income2007")], base08[,sel]), id.vars = c("IncGrp", "income2007"), value.name = "p")
identical(dat1[,c("IncGrp", "income2007", "variable")], dat2[,c("IncGrp", "income2007", "variable")])
tmp3	<- compare_fn(dat1, dat2, byvar = "variable", weight.dat = tab.inc, weight.by = c("IncGrp", "income2007"), 
						fml = "x2-x1", val1 = "p", val2 = "p", agg.first = FALSE)

sum.tab	<- data.frame(Retail = tmp1$variable, ExpRatio = tmp1$y, ShareDifference = tmp2$y, Incidence = tmp3$y)
sum.tab	<- sum.tab[order(sum.tab$Retail),]
print(sum.tab)
stargazer(sum.tab, summary=FALSE, rownames = FALSE)

# ------------------# 
# Pure price effect # 
tmp1	<- do.call(rbind, lapply(1:length(sim.ls), function(i) cbind(IncGrp = segname[i], sim.ls[[i]]$price07[, c("retailer", "change", "lnInc", fmt_name)])) )
tmp1$income2007	<- round(exp(tmp1$lnInc), 2)
tmp1$income.level	<- 2007
tmp2	<- do.call(rbind, lapply(1:length(sim.ls), function(i) cbind(IncGrp = segname[i], sim.ls[[i]]$price08[, c("retailer", "change", "lnInc", fmt_name)])) )
tmp2$income2007	<- round(exp(tmp2$lnInc)/(1+iv.change), 2)
tmp2$income.level	<- 2008
tmpsim	<- rbind(tmp1, tmp2)

# We consider only own price elasticity
tmpsim	<- subset(tmpsim, change == max(tmpsim$change))		# compute effect at price change 10%
tmpsim	<- melt(tmpsim, id.vars = setdiff(names(tmpsim), fmt_name), value.name = "exp")
tmpsim1	<- tmpsim
tmpsim	<- subset(tmpsim, retailer == variable)
tmpsim	<- tmpsim[order(tmpsim$income.level, tmpsim$variable, tmpsim$IncGrp, tmpsim$income2007),]
tmpbase	<- rbind(cbind(tmp07, income.level = 2007), cbind(tmp08, income.level = 2008))
rownames(tmpbase) 	<- rownames(tmpsim)	<- NULL
identical(tmpbase[,c("IncGrp", "income2007", "variable")], tmpsim[,c("IncGrp", "income2007", "variable")])

# Compute price elasticity at income of 2007
tmp1	<- compare_fn(tmpbase[tmpbase$income.level==2007,], tmpsim[tmpsim$income.level==2007,], byvar = c("income.level", "variable"), 
					weight.dat = tab.inc, weight.by = c("IncGrp", "income2007"), 
					fml = "x2/x1", val1 = "exp", val2 = "exp", agg.first = FALSE)
tmp1$PriceElst	<- with(tmp1, (y-1)/abs(iv.change))					
tmp1$bmargin	<- with(tmp1, (1-delta)*y/(y-delta))
tmp1	<- tmp1[,c("variable", "PriceElst", "bmargin")]
colnames(tmp1)	<- c("Retail", "Price Elasticity07", "Break-even margin07")
tmp1	<- tmp1[order(tmp1$Retail),]

# Compute price elasticity at income of 2008
tmp2	<- compare_fn(tmpbase[tmpbase$income.level==2008,], tmpsim[tmpsim$income.level==2008,], byvar = c("income.level", "variable"), 
					weight.dat = tab.inc, weight.by = c("IncGrp", "income2007"), 
					fml = "x2/x1", val1 = "exp", val2 = "exp", agg.first = FALSE)
tmp2$PriceElst	<- with(tmp2, (y-1)/abs(iv.change))					
tmp2$bmargin	<- with(tmp2, (1-delta)*y/(y-delta))
tmp2	<- tmp2[,c("variable", "PriceElst", "bmargin")]
colnames(tmp2)	<- c("Retail", "Price Elasticity08", "Break-even margin08")
tmp2	<- tmp2[order(tmp2$Retail),]

# Combine income effect and price effect
sum.tab	<- cbind(sum.tab, tmp1[,-1], tmp2[,-1])
sum.tab$ShareDifference <- paste(round(sum.tab$ShareDifference*100, 2),"%", sep="")
rownames(sum.tab)	<- NULL
cat("Summary of the counterfactual results:\n"); print(sum.tab); cat("\n")

# -----------------------------------------------------#
# Test if price elasticity increases with lower income #
selcol	<- c("IncGrp", "income2007", "variable")
identical(tmpbase[,c("income.level", "IncGrp", "income2007", "variable")], tmpsim[,c("income.level", "IncGrp", "income2007", "variable")])
tmp		<- data.frame(tmpbase[,c("income.level", "IncGrp", "income2007", "variable")], els = (tmpsim$exp/tmpbase$exp - 1)/abs(iv.change))
sel		<- tmp$income.level == 2007
tmp1	<- compare_fn(tmp[sel, c(selcol, "els")], tmp[!sel,c(selcol, "els")], 
					byvar = c("variable"), weight.dat = tab.inc, weight.by = c("IncGrp", "income2007"), 
					fml = "x2 - x1", val1 = "els", val2 = "els", agg.first = FALSE)
sel		<- abs(tmp1$y/tmp1$y.sd) >	qnorm(test.alpha)				
cat("Testing if price elasticity increases with lower income :\n"); print(tmp1); cat("\n")
cat("Number of store with significant price elasticity difference is", sum(sel), ".\n")

# --------------------------------------------------------# 
# Examine if price cuts can compensate expenditure losses #
tmp1	<- compare_fn(tmpbase[tmpbase$income.level==2007,], tmpsim[tmpsim$income.level==2008,], 
					byvar = c("variable"), weight.dat = tab.inc, weight.by = c("IncGrp", "income2007"), 
					fml = "x2/x1", val1 = "exp", val2 = "exp", agg.first = FALSE)
tmp1	<- tmp1[order(tmp1$variable),]

# Export the table 
sum.tab	<- cbind(sum.tab, BothIncomePrice = tmp1$y)
rownames(sum.tab) <- NULL
cat("Summary of the counterfactual simulations:\n"); print(sum.tab); cat("\n")
selcol	<- c("Retail", "Price Elasticity07", "Break-even margin07", "Price Elasticity08", "Break-even margin08", "BothIncomePrice")
stargazer(sum.tab[,selcol], summary = FALSE, rownames = FALSE, digits = 2)

# ---------------------- #
# Cross price elasticity #
# Use income.level == 2008 and baseline 2008
crs.els	<- matrix(NA, R, R, dimnames = list(Retailer = fmt_name, Impact = fmt_name))
tmp		<- subset(tmpsim1, income.level == 2008)
tmp		<- tmp[order(tmp$retailer, tmp$variable, tmp$IncGrp, tmp$income2007),]
for(i in 1:R){
	tmp1<- compare_fn(tmp08, tmp[tmp$retailer == fmt_name[i],],
						byvar = c("variable"), weight.dat = tab.inc, weight.by = c("IncGrp", "income2007"), 
						fml = "(x2/x1-1)/(delta-1)", val1 = "exp", val2 = "exp", agg.first = FALSE)
	tmp1<- tmp1[order(tmp1$variable),]					
	crs.els[i,]	<- tmp1$y
}
cat("Price elasticity matrix is:\n"); print(round(crs.els, 2)); cat("\n")
stargazer(crs.els, summary = FALSE, digits = 2)

########################
# Effect of price cuts # 
########################
fmt_label	<- setNames(gsub("\\s", "\n", fmt_name), fmt_name)
ggtmp	<- do.call(rbind, lapply(1:length(sim.ls), function(i) cbind(IncGrp = segname[i], sim.ls[[i]]$price[, c("retailer", "change", "lnInc", fmt_name)])) )
ggtmp$income2007	<- round(exp(ggtmp$lnInc)/(1+iv.change), 2)
ggtmp	<- merge(ggtmp, tab.inc, by = c("IncGrp", "income2007"), all.x=T)
ggtmp	<- ggtmp[order(ggtmp$retailer, ggtmp$change, ggtmp$IncGrp, ggtmp$income2007),]

# Compute the expenditure ratio
# The ratio of expenditure under price cuts and base07
tmpn	<- length(unique(ggtmp$retailer)) * length(unique(ggtmp$change))
identical(tmpn*nrow(base07), nrow(ggtmp))
tmp		<- kronecker(rep(1, tmpn), as.matrix(base07[,fmt_name]))
ggtmp1	<- ggtmp
ggtmp1[,fmt_name]	<- ggtmp1[,fmt_name]/tmp
# Add expenditure ratio at price change = 0
tmp		<- base08[, fmt_name]/base07[,fmt_name]
tmp		<- data.frame(base08[,c("IncGrp","income2007","lnInc", "nhh")], tmp, change = 0)
names(tmp)	<- gsub("\\.", " ", names(tmp))
tmp1	<- do.call(rbind, lapply(fmt_name, function(x) cbind(tmp, retailer = x)))
ggtmp1	<- rbind(ggtmp1, tmp1[, names(ggtmp1)])
ggtmp1$base	<- "Base2007"

tmp		<- kronecker(rep(1, tmpn), as.matrix(base08[,fmt_name]))
ggtmp2	<- ggtmp
ggtmp2[,fmt_name]	<- ggtmp2[,fmt_name]/tmp
ggtmp2$base	<- "Base2008"

# Summerzie the average effect: own-price effect
ggtmp3	<- melt(rbind(ggtmp1, ggtmp2), id.vars = setdiff(names(ggtmp1), fmt_name), value.name = "ratio")
ggtmp3	<- data.table(ggtmp3)
setnames(ggtmp3, "variable", "ImpactRetailer")
ggtmp1	<- ggtmp3[retailer == ImpactRetailer,list(ratio = weight_summary(ratio, nhh, mean), ratio.sd = weight_summary(ratio, nhh, sd) ), 
				by = list(base, retailer, change)] 
ggtmp2	<- ggtmp3[retailer == ImpactRetailer,list(ratio = weight_summary(ratio, nhh, mean), ratio.sd = weight_summary(ratio, nhh, sd) ), 
				by = list(base, retailer, change, IncGrp)]
				
# Compute price elasticity
tmp		<- subset(ggtmp1, base == "Base2008" & change == max(ggtmp1$change))
els		<- setNames(as.numeric((tmp$ratio - 1)/tmp$change), tmp$retailer)
els		<- round(els, 1)

# Order retailers
ret.ord			<- c("Warehouse Club", "Convenience Store", "Drug Store", "Dollar Store", "Discount Store", "Grocery")
ggtmp1$retailer <- factor(ggtmp1$retailer, levels = ret.ord, labels = paste(ret.ord," (",els[ret.ord],")", sep=""))
ggtmp2$retailer <- factor(ggtmp2$retailer, levels = ret.ord, labels = fmt_label[ret.ord])

# Plot the effect
plots	<- list(NULL)
plots[[1]]	<- ggplot(ggtmp1, aes(-change, ratio)) + geom_line() + geom_point() + 
				geom_hline(yintercept = 1, size = .25, linetype = 2) +
				facet_grid(retailer ~ base)
plots[[2]]	<- ggplot(ggtmp2, aes(-change, ratio, color = IncGrp)) + geom_line() + geom_point() + 
				facet_grid(retailer ~ base)

pdf(paste(plot.wd, "/graph_ctrfact_price.pdf", sep=""), width = ww, height = ww*ar)				
p1 <- 	ggplot(subset(ggtmp1, base=="Base2008"), aes(-change, ratio)) + geom_line() + geom_point() + 
			geom_hline(yintercept = 1, size = .25, linetype = 2) +
			facet_wrap(~ retailer) + 
			scale_x_continuous(labels = percent) +
			labs(x = "Price change (%)", y = "Expenditure ratio", title = "Effect of price cuts relative to baseline 2008")
p2 <- 	ggplot(subset(ggtmp1, base=="Base2007"), aes(-change, ratio)) + geom_line() + geom_point() + 
			geom_hline(yintercept = 1, size = .25, linetype = 2) +
			facet_wrap(~ retailer) + 
			scale_x_continuous(labels = percent) +
			labs(x = "Price change (%)", y = "Expenditure ratio", title = "Effect of price cuts relative to baseline 2007")
print(p1); 
print(p2)			
dev.off()	

# Overall profit
ggtmp1	<- data.table(melt(ggtmp, id.vars = setdiff(names(ggtmp), fmt_name)))			
ggtmp1	<- ggtmp1[retailer==variable,]
ggtmp1	<- ggtmp1[, list(revenue = sum(value*nhh)), by = list(retailer, change)]
tmp1	<- colSums(base07[,fmt_name]*base07$nhh)
tmp2	<- colSums(base08[,fmt_name]*base07$nhh)
ggtmp1$Base2007	<- tmp1[as.character(ggtmp1$retailer)]
ggtmp1$Base2008	<- tmp2[as.character(ggtmp1$retailer)]
ggtmp1	<- melt(ggtmp1, id.vars = c("retailer", "change", "revenue"))
ggtmp1	<- ggtmp1[,ratio:=revenue/value]
setnames(ggtmp1, "variable", "base")

# Add margin comparison
# Suppose the old margin is m0, then m0/m1 = (1 + delta)*m0/(delta + m0), where delta = (p1-p0)/p0
p <- ggplot(subset(ggtmp1, base == "Base2007"), aes(-change, ratio)) + geom_line() + geom_point() + 
			geom_hline(yintercept = 1, size = .25, linetype = 2) +
			# facet_grid(retailer ~ .)
			facet_wrap(~ retailer) + 
			scale_x_continuous(labels = percent) +
			xlab("Price change (%)") + ylab("Expenditure ratio")
			
m0		<- .5
delta 	<- seq(-.1, 0, .01)
tmpr	<- (1 + delta) * m0 / (delta + m0)
tmp		<- data.frame(change = delta, ratio = tmpr)
mylim	<- range(ggtmp1$ratio)
p1		<- p	+ geom_line(data = tmp, aes(change, ratio), linetype = 3) + ylim(mylim)
print(p1)

pdf(paste(plot.wd, "/graph_ctrfact_price_margin.pdf", sep=""), width = ww, height = ww*ar)				
print(p1)
dev.off()



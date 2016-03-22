library(ggplot2)
library(reshape2)
library(Rcpp)
library(RcppArmadillo)
library(maxLik)
library(evd)
library(data.table)
library(scales)

# setwd("~/Documents/Research/Store switching/processed data")
# plot.wd		<- "/Users/chaoqunchen/Desktop"

setwd("/home/brgordon/ccv103/Exercise/run")
run_id		<- 4
ww			<- 6.5
ar			<- .8

# outfile 	<- paste(plot.wd, "/counter1.xlsx", sep="")
# mywb		<- createWorkbook()
# sht1		<- createSheet(mywb, "Counterfactual1")
# make_plot	<- TRUE

# Set plotting parameters
out.stat	<- "level"		# Output statistics -- "level" or "elasticity"
ci.method	<- "error"	# Confidence interval method -- "crossHH" or "error"
annual.week	<- 26
annual		<- TRUE
draw.par		<- FALSE
sim.omega		<- FALSE

###########################
# Read counterfactal data # 
###########################
# Read estimation from all segments #
nseg	<- 3
sim.ls	<- setNames(vector("list", length = nseg), c("Low", "Med", "High"))
my.simdata	<- data.frame()
sim.df		<- data.frame()
numsim		<- 1000
ver.date1	<- "2016-03-01" 
tmpls	<- ls()

for(ii in 1:nseg){
	# Load data 
	# load(paste("Estimation/estrun_2/ctrfact_seg",ii, "_", ver.date1, ".rdata",sep=""))
	fname			<- paste("ctrfact_seq_seg",ii, sep="")
	if(draw.par){ fname <- paste(fname, "_pardraw", sep="") }
	if(sim.omega){fname	<- paste(fname, "_simomega", sep="") }
	fname			<- paste(fname, "_sim", numsim, "_", ver.date1, sep="")
	load(paste(fname, ".rdata",sep=""))
	# load(paste("estrun_",run_id, "/", fname, ".rdata",sep=""))
	# 
	# # NOTE: the structure of sim.X.ls -- list[4 intervention variables][6 retailers][Average, Allocation]
	# # Stack the baseline average and average simulation results of each intervention
	# tmp			<- do.call(rbind, lapply(sim.prc.ls, function(x) x$Average))
	# tmp			<- rbind(sim.X.df, cbind(tmp, Var = "price"))
	# tmp$income2007	<- exp(tmp$lnInc)/(1 - .1)
	# tmp1		<- data.frame(sim.base07$Average, lnInc = sim.unq$Inc07, retailer = "", Var = "Base07", income2007 = sim.unq$income2007)
	# tmp2		<- data.frame(sim.base08$Average, lnInc = sim.unq$Inc08, retailer = "", Var = "Base08", income2007 = sim.unq$income2007)
	# tmp			<- rbind(tmp, tmp1, tmp2)
	# 
	# # Compute number of households in each income level 
	# inc.tab		<- table(sim.data$income)
	# tmp$nhh		<- inc.tab[as.character(tmp$income2007)]
	# sim.df		<- rbind(sim.df, cbind(tmp, IncGrp = names(sim.ls)[ii]))
	
	# Combine counterfactual simulation together
	tmp			<- setNames(vector("list", length(sim.X.ls) +1), c(names(sim.X.ls), "price"))			# A list of interventions
	for(j in 1:length(sim.X.ls)){
		tmp1	<- sapply(sim.X.ls[[j]], function(x) levels(x$Average$retailer))
		tmp[[j]]<- setNames(lapply(sim.X.ls[[j]], function(x) x$Allocation), tmp1)
	}
	tmp1	<- sapply(sim.prc.ls, function(x) levels(x$Average$retailer))
	tmp[[length(sim.X.ls) +1]]	<- setNames(lapply(sim.prc.ls, function(x) x$Allocation), tmp1)

	# If we want to convert the biweekly simulation to annual level
	if(annual){
		tmp	<- lapply(tmp, function(l) {lapply(l, function(ll) ll*annual.week)})
		sim.base07$Allocation	<- sim.base07$Allocation * annual.week
		sim.base08$Allocation	<- sim.base08$Allocation * annual.week
	}
	# dimnames(sim.base07$Allocation)	<- arr.name
	# dimnames(sim.base08$Allocation)	<- arr.name
	sim.ls[[ii]]<- c(list(Base07 = sim.base07$Allocation, Base08 = sim.base08$Allocation), tmp)
						
	my.simdata	<- rbind(my.simdata, data.frame(sim.data, IncGrp = names(sim.ls)[ii]))
	
	# Delete unused objects
	rm(list = setdiff(ls(), c(tmpls, "tmpls", "ii", "fmt_name", "R", "iv.change.vec")))
	print(ii)						
}
trim.alpha		<- .05
plot.wd			<- paste(getwd(), "/estrun_", run_id, sep="")

# Set small simulation value to zero 
zero.max		<- .01	# 1 cent
for(i in 1:nseg){
	for(j in 1:length(sim.ls[[1]])){
		if(j > 2){
			for(k in 1:R){
				sim.ls[[i]][[j]][[k]]<- 1*(sim.ls[[i]][[j]][[k]] >= zero.max) * sim.ls[[i]][[j]][[k]]
			}
		}else{
			sim.ls[[i]][[j]]	<- 1*(sim.ls[[i]][[j]] >= zero.max) * sim.ls[[i]][[j]]
		}
	}
}

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

#################
# Income effect #
#################
iv.change 	<- .1
p.idx			<- 0
plots			<- list(NULL)

# Income frequency
tab.inc	<- data.table(my.simdata)
tab.inc	<- tab.inc[, list(nhh = length(unique(household_code))), by = list(IncGrp, income)]
setnames(tab.inc, "income", "income2007")
setkeyv(tab.inc, c("IncGrp", "income2007"))
tab.inc$income2007	<- round(tab.inc$income2007, 2)

# --------------------- # 
# Average income effect #
# Stack all the simlation outcome from all iterations
sim.dt	<- data.frame(NULL)
for(i in 1:nseg){
	tmp1	<- cbind(melt(sim.ls[[i]]$Base07, value.name = "Base07"), 
					Base08 = melt(sim.ls[[i]]$Base08)$value, 
					Difference = melt(sim.ls[[i]]$Base08 - sim.ls[[i]]$Base07)$value, 
					Elasticity = melt(-1/iv.change*(sim.ls[[i]]$Base08 - sim.ls[[i]]$Base07)/sim.ls[[i]]$Base07)$value)
	tmp1$income2007	<- round(exp(tmp1$lnInc),2)				
	tmp1$PTC	<- tmp1$Difference/(-iv.change*tmp1$income2007)
	sim.dt	<- rbind(sim.dt, cbind(IncGrp = names(sim.ls)[i], tmp1))
}
sim.dt$Elasticity	<- ifelse(sim.dt$Elasticity == -Inf, NA, sim.dt$Elasticity)
sim.dt$retailer		<- fmt_name[sim.dt$retailer]

# Compute expected value for each Income group - income level by averaging over iterations
sim.dt	<- data.table(sim.dt)
agg.dt	<- sim.dt[, list(	Base07 = mean(Base07, trim = trim.alpha, na.rm = TRUE), 
							Base08 = mean(Base08, trim = trim.alpha, na.rm = TRUE), 
							Difference = mean(Difference, trim = trim.alpha, na.rm = TRUE), 
							Diff.se = sd(mytrimfun(Difference, trim.alpha), na.rm = T), 
							PTC		= mean(PTC, trim = trim.alpha, na.rm = TRUE), 
							PTC.se	= sd(mytrimfun(PTC, trim.alpha), na.rm = T),
							Elasticity	= mean(Elasticity, trim = trim.alpha, na.rm = TRUE), 
							Elast.se	= sd(mytrimfun(Elasticity, trim.alpha), na.rm = TRUE)
							), 
					by = list(IncGrp, income2007, retailer)]	

# Share measure
shr.dt	<- sim.dt[,':='(Base07 = Base07/sum(Base07)*100, Base08 = Base08/sum(Base08)*100), 
					by = list(IncGrp, income2007, iter)]
shr.dt	<- shr.dt[,':='(Difference = Base08 - Base07, 
						PTC = (Base08 - Base07)/(-iv.change*income2007), 
						Elasticity = -1/iv.change*(Base08-Base07)/Base07)]
shr.dt$Elasticity	<- ifelse(shr.dt$Elasticity == -Inf, NA, shr.dt$Elasticity)									
shr.dt	<- shr.dt[, list(	Base07 = mean(Base07, trim = trim.alpha, na.rm = TRUE), 
							Base08 = mean(Base08, trim = trim.alpha, na.rm = TRUE), 
							Difference = mean(Difference, trim = trim.alpha, na.rm = TRUE), 
							Diff.se = sd(mytrimfun(Difference, trim.alpha), na.rm = T), 
							PTC		= mean(PTC, trim = trim.alpha, na.rm = TRUE), 
							PTC.se	= sd(mytrimfun(PTC, trim.alpha), na.rm = T),
							Elasticity	= mean(Elasticity, trim = trim.alpha, na.rm = TRUE), 
							Elast.se	= sd(mytrimfun(Elasticity, trim.alpha), na.rm = TRUE)
							), 
					by = list(IncGrp, income2007, retailer)]

# Merge in income frequency
dim(agg.dt)
agg.dt	<- merge(agg.dt, tab.inc, by = c("IncGrp", "income2007"))
dim(agg.dt)
shr.dt	<- merge(shr.dt, tab.inc, by = c("IncGrp", "income2007"))

# Compute overall income effect on expenditure 
# NOTE: expenditure does not vary across random draws, so we report across household variation 
tmp.tab	<- agg.dt[, list(Base07 = sum(Base07), Base08 = sum(Base08)), by = list(IncGrp, income2007, nhh)]
tmp.tab	<- tmp.tab[, list(	Base07 = weight_summary(Base07, nhh, mean), 
						  	Base08 = weight_summary(Base08, nhh, mean), 
							Difference = weight_summary(Base08 - Base07, nhh, mean), 
							Diff.sd	= weight_summary(Base08 - Base07, nhh, sd), 
							PTC		= weight_summary((Base08 - Base07)/(-iv.change*income2007), nhh, mean), 
							PTC.sd	= weight_summary((Base08 - Base07)/(-iv.change*income2007), nhh, sd), 
							Elasticity = weight_summary(-1/iv.change*(Base08 - Base07)/Base07, nhh, mean), 
							Elast.sd = weight_summary(-1/iv.change*(Base08 - Base07)/Base07, nhh, sd))]
cat("Overall average income effect on expenditure (across household variation):\n"); print(tmp.tab); cat("\n")

# Compute overall income effect on expenditure at retailers by taking weighted average
tmp.tab	<- agg.dt[, list(Base07 = weight_summary(Base07, nhh, mean), 
						Base08 = weight_summary(Base08, nhh, mean), 
						Difference = weight_summary(Difference, nhh, mean), 
						Diff.se	= weight_se(Diff.se, nhh), 
						PTC = weight_summary(PTC, nhh, mean), 
						PTC.se = weight_se(PTC.se, nhh),
						Elasticity = weight_summary(Elasticity, nhh, mean), 
						Elast.se = weight_se(Elast.se, nhh)), 
					by = list(retailer)]
tmp.tab	<- data.frame(tmp.tab)					
cat("Averge income effect on expenditure at retailers:\n"); print(cbind(tmp.tab[,1], round(tmp.tab[,-1], 3))); cat("\n")

tmp.tab	<- shr.dt[, list(Base07 = weight_summary(Base07, nhh, mean), 
						Base08 = weight_summary(Base08, nhh, mean), 
						Difference = weight_summary(Difference, nhh, mean), 
						Diff.se	= weight_se(Diff.se, nhh), 
						PTC = weight_summary(PTC, nhh, mean), 
						PTC.se = weight_se(PTC.se, nhh),
						Elasticity = weight_summary(Elasticity, nhh, mean), 
						Elast.se = weight_se(Elast.se, nhh)), 
					by = list(retailer)]
tmp.tab	<- data.frame(tmp.tab)					
cat("Averge income effect on expenditure share at retailers:\n"); print(cbind(tmp.tab[,1], round(tmp.tab[,-1], 4))); cat("\n")

# Plot the share difference 
ord		<- tmp.tab[order(tmp.tab$Difference), "retailer"]
tmp.tab$retailer <- factor(tmp.tab$retailer, levels = ord)
p		<- ggplot(tmp.tab, aes(x = retailer, y = Difference)) + 
		geom_bar(stat = "identity", position = position_dodge(width=0.9)) + 
		geom_errorbar(aes(ymin = Difference-1.96*Diff.se, ymax = Difference + 1.96 * Diff.se), position=position_dodge(width=0.9), width=0.25) + 
		labs(ylad = "Share difference") + 
		coord_flip()
p.idx	<- p.idx + 1
plots[[p.idx]]	<- p

#-----------------------------#
# Heterogeneous income effect #
tmp.tab	<- agg.dt[, list(Base07 = sum(Base07), Base08 = sum(Base08)), by = list(IncGrp, income2007, nhh)]
tmp.tab	<- tmp.tab[, list(	Base07 = weight_summary(Base07, nhh, mean), 
						  	Base08 = weight_summary(Base07, nhh, mean), 
							Difference = weight_summary(Base08 - Base07, nhh, mean), 
							Diff.sd	= weight_summary(Base08 - Base07, nhh, sd), 
							PTC		= weight_summary((Base08 - Base07)/(-iv.change*income2007), nhh, mean), 
							PTC.sd	= weight_summary((Base08 - Base07)/(-iv.change*income2007), nhh, sd),
							Elasticity = weight_summary(-1/iv.change*(Base08 - Base07)/Base07, nhh, mean), 
							Elast.sd = weight_summary(-1/iv.change*(Base08 - Base07)/Base07, nhh, sd)), 
						by = list(IncGrp)]
cat("Heterogeneous average income effect on expenditure (across household variation):\n"); print(tmp.tab); cat("\n")

# Plot across household variation of income effect on overall expenditure
ggtmp	<- agg.dt[, list(Base07 = sum(Base07), Base08 = sum(Base08)), by = list(IncGrp, income2007, nhh)]
ggtmp	<- melt(ggtmp, id.vars = c("IncGrp", "income2007", "nhh"))
p.idx	<- p.idx + 1
plots[[p.idx]]	<- ggplot(ggtmp, aes(IncGrp, value, fill = variable)) + geom_boxplot(aes(weight = nhh), outlier.shape = NA) + 
		labs(x = "Income Group", y = "Overall expenditure", 
		title = "Cross household variation of overall expenditure")

# Heterogeneous income effect on expenditure at retailers
tmp.tab	<- agg.dt[, list(Base07 = weight_summary(Base07, nhh, mean), 
						Base08 = weight_summary(Base07, nhh, mean), 
						Difference = weight_summary(Difference, nhh, mean), 
						Diff.se	= weight_se(Diff.se, nhh), 
						PTC = weight_summary(PTC, nhh, mean), 
						PTC.se = weight_se(PTC.se, nhh),
						Elasticity = weight_summary(Elasticity, nhh, mean), 
						Elast.se = weight_se(Elast.se, nhh)), 
					by = list(IncGrp, retailer)]
tmp.tab	<- data.frame(tmp.tab)					
sel 	<- sapply(1:ncol(tmp.tab), function(i) is.numeric(tmp.tab[,i]))					
cat("Heterogeneous income effect on expenditure at retailers:\n"); 
print(cbind(tmp.tab[,!sel], round(tmp.tab[,sel], 4))); cat("\n")
sel		<- with(tmp.tab, Difference/Diff.se > 1.96)
cat("Number of significant level difference =", sum(sel), ".\n")
sel		<- with(tmp.tab, Elasticity/Elast.se > 1.96)
cat("Number of significant Elasticity =", sum(sel), ".\n")

# Plot expenditure at reatilers by income group
p.idx	<- p.idx + 1
if(ci.method == "error"){
	tmp		<- sim.dt[,list(	Base07 = mean(mytrimfun(Base07, trim.alpha), na.rm = T), 
								Base08 = mean(mytrimfun(Base07, trim.alpha), na.rm = T), 
								Base07.se = sd(mytrimfun(Base07, trim.alpha), na.rm = T), 
								Base08.se = sd(mytrimfun(Base08, trim.alpha), na.rm = T)),  
						by = list(IncGrp, income2007, retailer)]
	tmp		<- merge(tmp, tab.inc, by = c("IncGrp", "income2007"))
	tmp		<- tmp[, list(	Base07 = weight_summary(Base07, nhh, mean), 
							Base08 = weight_summary(Base07, nhh, mean),
							Base07.se = weight_se(Base07.se, nhh), 
							Base08.se = weight_se(Base08.se, nhh)), 
					by = list(IncGrp, retailer)]
	ggtmp	<- cbind(melt(tmp[,list(IncGrp, retailer, Base07, Base08)], id.vars = c("IncGrp", "retailer")), 
					se = melt(tmp[,list(IncGrp, retailer, Base07.se, Base08.se)], id.vars = c("IncGrp", "retailer"))$value)
	plots[[p.idx]]	<- ggplot(ggtmp, aes(IncGrp, color = variable)) + 
				geom_pointrange(aes(y = value, ymin = value - 1.96*se, ymax = value + 1.96*se), position = position_dodge(width=0.8)) + 
				facet_wrap(~retailer, scales = "free") +
				labs(title = "Impact on expenditure at retail format of 10% decrease of income")
	# print(p)
}else if(ci.method == "crossHH"){
	tmp		<- sim.dt[,list(	Base07 = mean(mytrimfun(Base07, trim.alpha), na.rm = T), 
								Base08 = mean(mytrimfun(Base07, trim.alpha), na.rm = T)),  
						by = list(IncGrp, income2007, retailer)]
	ggtmp	<- merge(tmp, tab.inc, by = c("IncGrp", "income2007"))
	ggtmp	<- melt(ggtmp, id.vars = c("IncGrp", "income2007", "retailer", "nhh"))
	plots[[p.idx]]		<- ggplot(ggtmp, aes(IncGrp, value, fill = variable)) + 
				geom_boxplot(aes(weight = nhh), outlier.shape = NA) + 
				facet_wrap(~retailer, scales = "free") +
				labs(title = "Impact on expenditure at retail format of 10% decrease of income")
	# print(p)
}

# Plot elasticity of expenditure at retailers by income group
# ggplot(tmp.tab, aes(retailer, Elasticity, col = IncGrp)) + 
# 		geom_pointrange(aes(y = Elasticity, ymin = Elasticity - 1.96*Elast.se, ymax = Elasticity + 1.96*Elast.se), 
# 		position = position_dodge(width=0.8))  
		
# Heterogeneous income effects on expenditure share		
tmp.tab	<- shr.dt[, list(Base07 = weight_summary(Base07, nhh, mean), 
						Base08 = weight_summary(Base07, nhh, mean), 
						Difference = weight_summary(Difference, nhh, mean), 
						Diff.se	= weight_se(Diff.se, nhh), 
						Elasticity = weight_summary(Elasticity, nhh, mean), 
						Elast.se = weight_se(Elast.se, nhh)), 
					by = list(IncGrp, retailer)]
tmp.tab	<- data.frame(tmp.tab)					
sel 	<- sapply(1:ncol(tmp.tab), function(i) is.numeric(tmp.tab[,i]))					
cat("Heterogeneous income effect on expenditure at retailers:\n"); 
print(cbind(tmp.tab[,!sel], round(tmp.tab[,sel], 4))); cat("\n")
sel		<- with(tmp.tab, Difference/Diff.se > 1.96)
cat("Number of significant level difference =", sum(sel), ".\n")
sel		<- with(tmp.tab, Elasticity/Elast.se > 1.96)
cat("Number of significant Elasticity =", sum(sel), ".\n")

# Plot the share change by income group 
incg.col<- c("red","grey30", "grey50")			# Color of low, median, high income bar
tmp.tab$retailer	<- factor(tmp.tab$retailer, levels = ord)
p		<- ggplot(tmp.tab, aes(retailer, Difference, fill = IncGrp, col = IncGrp)) + 
			geom_bar(stat = "identity", position = position_dodge(width=0.9)) + 
			geom_errorbar(aes(ymin = Difference-1.96*Diff.se, ymax = Difference + 1.96 * Diff.se, col = IncGrp), position=position_dodge(width=0.9), width=0.25) + 
			scale_fill_manual(values = rev(incg.col)) + 
			scale_color_manual(values = rev(incg.col)) + 
			xlab("Channel") + ylab("Share difference")+ coord_flip() + 
			guides(fill = guide_legend(reverse = TRUE), color = guide_legend(reverse = TRUE))
p.idx	<- p.idx + 1
plots[[p.idx]]	<- p

pdf(paste(plot.wd, "/graph_ctrfact_seq_sim", numsim, "_inc_",ver.date1, ".pdf", sep=""), width = ww, height = ww*ar)
for(i in 1:length(plots)){
	print(plots[[i]])
}
dev.off()

####################
# Retail marketing # 
####################
ww		<- 9
# Construct plotting data
ggtmp	<- data.frame()
nx		<- 3:length(sim.ls[[1]])

# Structure of sim.ls: list[3 income group][2 baseline + 5 variables][6 retailers] -- array[numsim, income level, outcome retailer]
for(i in nx){
	for(j in 1:R){
		for(k in 1:nseg){
			tmp		<- array(0, dim(sim.ls[[k]][[i]][[j]]), dimnames = dimnames(sim.ls[[k]][[i]][[j]]))
			# For each change value, compute the difference of expenditure from the baseline
			for(l in 1:length(iv.change.vec)){
				tmp[l,,,]	<- sim.ls[[k]][[i]][[j]][l,,,]	- sim.ls[[k]]$Base08
			}		
			tmp		<- tmp[,,,names(sim.ls[[k]][[i]])[j]]					# Only focus on own effect
			tmp1	<- apply(tmp, c(1,3), mean, trim = trim.alpha, na.rm=T)	
			tmp2	<- apply(tmp, c(1,3), function(x) sd(mytrimfun(x, alpha = trim.alpha), na.rm=T))
			tmp1	<- melt(tmp1, value.name = "Difference")
			tmp2	<- melt(tmp2, value.name = "se")
			tmp		<- data.frame(Var = names(sim.ls[[k]])[i], IncGrp = names(sim.ls)[k], retailer = names(sim.ls[[k]][[i]])[j], 
								  tmp1, se = tmp2[,"se"])
			ggtmp	<- rbind(ggtmp, tmp)
		}		
	}
}
unique(ggtmp$lnInc)
ggtmp$income2007	<- round(exp(ggtmp$lnInc) / .9, 2)
ggtmp	<- merge(ggtmp, tab.inc, by = c("IncGrp", "income2007"))
ggtmp	<- data.table(ggtmp)
ggtmp$retailer <- gsub("\\s", "\n", as.character(ggtmp$retailer))

# Combine together overall and income groups
ggtmp1	<- ggtmp[, list(Difference = weight_summary(Difference, nhh, mean), 
						Diff.se = weight_se(se, nhh)), 
					by = list(IncGrp, Var, retailer, change)]
tmp		<- ggtmp[, 	list(Difference = weight_summary(Difference, nhh, mean), 
						Diff.se = weight_se(se, nhh)), 
					by = list(Var, retailer, change)]
ggtmp1	<- rbind(ggtmp1, cbind(IncGrp = "Overall", tmp))
ggtmp1$group	<- ifelse(ggtmp1$IncGrp == "Overall", "Overall", "Income group")
ggtmp1$group	<- factor(ggtmp1$group, levels = c("Overall", "Income group"))
var.unq		<- levels(ggtmp1$Var)
tmplab			<- c("size index", "No. UPC per category", "No. categories", "private label", "price")
cbind(var.unq, tmplab)

plots	<- list(NULL)	
my.color	<- c("#66c2a5", "#fc8d62", "#8da0cb", "black")
for(i in 1:length(var.unq)){
	plots[[i]]<- 	ggplot(subset(ggtmp1, Var %in% var.unq[i]), aes(change, Difference, col = IncGrp)) + 	
					geom_pointrange(aes(y = Difference, ymin=Difference - 1.96*Diff.se, ymax=Difference + 1.96*Diff.se), 
									position=position_dodge(width=0.005)) + 
					geom_line() + 
					geom_hline(yintercept = 0, linetype = 2, size = .25) + 
					facet_grid(retailer~group, scales = "free") + 
					scale_color_manual(values = my.color) + 
					labs(x = "Percentage change", y = "Expenditure difference ($)", 
						title = paste("The effect of ", tmplab[i], sep=""))
}	

pdf(paste(plot.wd, "/graph_ctrfact_seq_sim", numsim, "_", ver.date1,".pdf", sep=""), width = ww, height = ww)
for(i in 1:length(plots)){
	print(plots[[i]])
}
dev.off()

#  Slides buildup
pdf(paste(plot.wd, "/graph_ctrfact_seq_sim", numsim, "_slide_", ver.date1,".pdf", sep=""), width = 6.5, height = 6.5)
print(ggplot(subset(ggtmp1, Var == "price" & IncGrp == "Overall" & retailer == "Grocery"), aes(change, Difference)) + 	
				geom_pointrange(aes(y = Difference, ymin=Difference - 1.96*Diff.se, ymax=Difference + 1.96*Diff.se)) + 
				geom_line() + 
				geom_hline(yintercept = 0, linetype = 2, size = .25) + 
				facet_grid(retailer~group, scales = "free") + 
				labs(x = "Percentage change", y = "Expenditure difference ($)"))
print(ggplot(subset(ggtmp1, Var == "price" & IncGrp == "Overall"), aes(change, Difference)) + 	
				geom_pointrange(aes(y = Difference, ymin=Difference - 1.96*Diff.se, ymax=Difference + 1.96*Diff.se)) + 
				geom_line() + 
				geom_hline(yintercept = 0, linetype = 2, size = .25) + 
				facet_grid(retailer~group, scales = "free") + 
				labs(x = "Percentage change", y = "Expenditure difference ($)")	)				
dev.off()

# --------------------- #
# Plot share elasticity #
# Construct plotting data
ggtmp	<- data.frame()
nx		<- 3:length(sim.ls[[1]])

# Structure of sim.ls: list[3 income group][2 baseline + 5 variables][6 retailers] -- array[numsim, income level, outcome retailer]
for(i in nx){
	for(j in 1:R){
		for(k in 1:nseg){
			tmp		<- array(0, dim(sim.ls[[k]][[i]][[j]]), dimnames = dimnames(sim.ls[[k]][[i]][[j]]))
			# For each change value, compute the difference of expenditure from the baseline
			for(l in 1:length(iv.change.vec)){
				tmp1		<- apply(sim.ls[[k]][[i]][[j]][l,,,], c(1,2), function(x) x/sum(x))
				tmp0		<- apply(sim.ls[[k]]$Base08, c(1,2), function(x) x/sum(x))
				tmp1		<- aperm(tmp1, c(3, 2, 1))
				tmp0		<- aperm(tmp0, c(3, 2, 1))
				tmp[l,,,]	<- tmp1	- tmp0
			}		
			
			tmp		<- tmp[,,,names(sim.ls[[k]][[i]])[j]]					# Only focus on own effect
			tmp1	<- apply(tmp, c(1,3), mean, trim = trim.alpha, na.rm=T)	
			tmp2	<- apply(tmp, c(1,3), function(x) sd(mytrimfun(x, alpha = trim.alpha), na.rm=T))
			tmp1	<- melt(tmp1, value.name = "Difference")
			tmp2	<- melt(tmp2, value.name = "se")
			tmp		<- data.frame(Var = names(sim.ls[[k]])[i], IncGrp = names(sim.ls)[k], retailer = names(sim.ls[[k]][[i]])[j], 
								  tmp1, se = tmp2[,"se"])
			ggtmp	<- rbind(ggtmp, tmp)
		}		
		cat("j=",j,"\n")
	}
	print(i)
}
unique(ggtmp$lnInc)
ggtmp$income2007	<- round(exp(ggtmp$lnInc) / .9, 2)
dim(ggtmp)
ggtmp	<- merge(ggtmp, tab.inc, by = c("IncGrp", "income2007"))
dim(ggtmp)
ggtmp	<- data.table(ggtmp)
ggtmp$retailer <- gsub("\\s", "\n", as.character(ggtmp$retailer))
 
# Combine together overall and income groups
ggtmp1	<- ggtmp[, list(Difference = weight_summary(Difference, nhh, mean), 
						Diff.se = weight_se(se, nhh)), 
					by = list(IncGrp, Var, retailer, change)]
tmp		<- ggtmp[, 	list(Difference = weight_summary(Difference, nhh, mean), 
						Diff.se = weight_se(se, nhh)), 
					by = list(Var, retailer, change)]
ggtmp1	<- rbind(ggtmp1, cbind(IncGrp = "Overall", tmp))
ggtmp1$group	<- ifelse(ggtmp1$IncGrp == "Overall", "Overall", "Income group")
var.unq		<- levels(ggtmp1$Var)
tmplab			<- c("size index", "No. UPC per category", "No. categories", "private label", "price")
cbind(var.unq, tmplab)

plots	<- list(NULL)	
my.color	<- c("#66c2a5", "#fc8d62", "#8da0cb", "black")
for(i in 1:length(var.unq)){
	plots[[i]]<- ggplot(subset(ggtmp1, Var %in% var.unq[i]), aes(change, Difference, col = IncGrp)) + 	
					geom_pointrange(aes(y = Difference, ymin=Difference - 1.96*Diff.se, ymax=Difference + 1.96*Diff.se)) + 
					geom_line() + 
					geom_hline(yintercept = 0, linetype = 2, size = .25) + 
					scale_y_continuous(labels=percent) + 
					facet_grid(retailer~group, scales = "free") + 
					scale_color_manual(values = my.color) + 
					labs(x = "Percentage change", y = "Share difference", 
						title = paste("The effect of ", tmplab[i], sep=""))
}	

pdf(paste(plot.wd, "/graph_ctrfact_seq_sim", numsim, "_share_", ver.date1,".pdf", sep=""), width = ww, height = ww)
for(i in 1:length(plots)){
	print(plots[[i]])
}
dev.off()


cat("This program is done. \n")

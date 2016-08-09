library(ggplot2)
library(reshape2)
library(data.table)
library(mgcv)
library(gridExtra)
library(scales)
options(error = quote({dump.frames(to.file = TRUE)}))

args <- commandArgs(trailingOnly = TRUE)
print(args)
if(length(args)>0){
    for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
    }
}

# setwd("~/Documents/Research/Store switching/processed data/Estimation")
# plot.wd	<- '~/Desktop'
# source("../Exercise/Multiple_discrete_continuous_model/0_Allocation_function.R")
# source("../Exercise/main/share_allocation/ctrfact_sim_functions.r")

# setwd("/home/brgordon/ccv103/Exercise/run")
# setwd("/kellogg/users/marketing/2661703/Exercise/run")
# setwd("/sscc/home/c/ccv103/Exercise/run")
setwd("U:/Users/ccv103/Documents/Research/Store switching/run")

# setwd("~/Documents/Research/Store switching/processed data/Estimation")
run_id		<- 11

###########################
# Read counterfactal data # 
###########################
# Read estimation from all segments #
nseg	<- 3
seg.name	<- c("Low", "Medium", "High")
sim.ls	<- setNames(vector("list", length = nseg), seg.name)
simavg.ls	<- setNames(vector("list", length = nseg), seg.name)
my.simdata	<- data.frame()
numsim		<- 500		#1000
ver.date1	<- "2016-07-24"		#"2016-06-25" 
cpi.adj		<- TRUE
annual		<- TRUE
annual.week	<- 26
tmpls	<- ls()

for(ii in 1:nseg){
	fname			<- paste("ctrfact_newf_seg",ii,sep="")
	loadf			<- paste("estrun_",run_id, "/",fname, "_sim", numsim, "_", ver.date1, ".rdata", sep="")
	load(loadf)
	
	# Convert biweekly simulations to annual simulations
	if(annual){
		sim.base07$Allocation	<- sim.base07$Allocation * annual.week
		sim.base08$Allocation	<- sim.base08$Allocation * annual.week		
		sel	<- c("Exp", fmt_name)
		sim.base07$Average[,sel]	<- sim.base07$Average[,sel] * annual.week
		sim.base08$Average[,sel]	<- sim.base08$Average[,sel] * annual.week
		sim.base07$y				<- sim.base07$y	* annual.week
		sim.base08$y				<- sim.base08$y	* annual.week
		for(j in 1:length(newf.sim)){
			newf.sim[[j]]$y 		<- newf.sim[[j]]$y * annual.week
			newf.sim[[j]]$Allocation 		<- newf.sim[[j]]$Allocation * annual.week
			newf.sim[[j]]$Average[,sel]		<- newf.sim[[j]]$Average[,sel] * annual.week
			newf.sim07[[j]]$y 		<- newf.sim07[[j]]$y * annual.week
			newf.sim07[[j]]$Allocation 		<- newf.sim07[[j]]$Allocation * annual.week
			newf.sim07[[j]]$Average[,sel]		<- newf.sim07[[j]]$Average[,sel] * annual.week
		}
	}
	
	# Combine all the draws
	tmp1	<- lapply(newf.sim07, function(x) x$Allocation)
	names(tmp1)	<- paste(names(tmp1), "_07", sep="")
	tmp2	<- lapply(newf.sim, function(x) x$Allocation)
	sim.ls[[ii]]<- c(list(Base07 = sim.base07$Allocation, Base08 = sim.base08$Allocation), tmp1, tmp2)
	
	# Combine the averge level: share and incidence
	tmp1	<- cbind(apply(sim.base07$Allocation, c(2,3), mean), 0, 
					 apply(1*(sim.base07$Allocation>0), c(2,3), mean), 0)
	tmp		<- exp(as.numeric(rownames(tmp1)))
	tmp2	<- cbind(apply(sim.base08$Allocation, c(2,3), mean), 0, 
					 apply(1*(sim.base08$Allocation>0), c(2,3), mean), 0)
	colnames(tmp1)	<- colnames(tmp2)	<- c(fmt_name, "New", paste("Trip_", fmt_name, sep=""), "Trip_New")
	tmp3	<- lapply(newf.sim07, function(x){ 
				out <- cbind(tmp, apply(x$Allocation, c(2,3), mean), apply(1*(x$Allocation>0), c(2,3), mean)); 
				colnames(out) <- c("income2007", fmt_name, "New", paste("Trip_", fmt_name, sep=""), "Trip_New"); return(out)} )
	names(tmp3)	<- paste(names(tmp3), "_07", sep="")
	tmp4	<- lapply(newf.sim, function(x){ 
				out <- cbind(tmp, apply(x$Allocation, c(2,3), mean), apply(1*(x$Allocation>0), c(2,3), mean)); 
				colnames(out) <- c("income2007", fmt_name, "New", paste("Trip_", fmt_name, sep=""), "Trip_New"); return(out)} )
	simavg.ls[[ii]]	<- c(list(Base07 = cbind(income2007 = tmp, tmp1), 
							  Base08 = cbind(income2007 = tmp, tmp2)), tmp3, tmp4)
	my.simdata	<- rbind(my.simdata, data.frame(sim.data, IncGrp = names(sim.ls)[ii]))
	
	# Delete unused objects
	rm(list = setdiff(ls(), c(tmpls, "tmpls", "ii", "fmt_name", "R", "ret.a", "ret.b")))
	print(ii)
}

# Size plotting parameteres 
make_plot	<- TRUE
ww			<- 10
ar			<- .5
pnt.size	<- .15		# Size in pointrange
position.width <- .5
plot.wd		<- "~/Desktop"

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

get_legend<-function(myggplot){
  # tmp <- ggplot_gtable(ggplot_build(myggplot))
	tmp	<- ggplotGrob(myggplot)
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

#################
# Income effect #
#################
# Income frequency
tab.inc	<- data.table(my.simdata)
tab.inc	<- tab.inc[, list(nhh = length(unique(household_code))), by = list(IncGrp, income)]
setnames(tab.inc, "income", "income2007")
setkeyv(tab.inc, c("IncGrp", "income2007"))
tab.inc$income2007	<- round(tab.inc$income2007, 2)

# Initiate names
newf.name	<- setdiff(names(simavg.ls[[1]]), c("Base07", "Base08"))
fmt_name1	<- c(fmt_name, "New")
my.change	<- -.1

# ---------------------------- #
# Income effect on expenditure # 
tmp	<- do.call(rbind, lapply(simavg.ls, function(x) cbind(income2007 = x$Base07[,"income2007"], 
											Exp2007 = rowSums(x$Base07[, fmt_name1]), 
											Exp2008 = rowSums(x$Base08[, fmt_name1]))) )
tmp	<- data.table(tmp)
tmp	<- tmp[, IncGrp:=unlist(lapply(1:nseg, function(i) rep(names(simavg.ls)[i], nrow(simavg.ls[[i]][[1]] )) ))]
tmp	<- tmp[, income2007:= round(income2007, 2)]
tmp	<- tmp[, ':='(Diff = Exp2008 - Exp2007, Elasticity = (Exp2008-Exp2007)/Exp2007/my.change, MPC = (Exp2008-Exp2007)/(my.change*income2007))]
tmp	<- merge(tmp, tab.inc, by = c("IncGrp", "income2007"), all.x = T)
tmp1	<- tmp[, list(Exp2007 = weight_summary(Exp2007, nhh, mean), Exp2008 = weight_summary(Exp2008, nhh, mean), 
					  Difference = weight_summary(Diff, nhh, mean), Diff_sd = weight_summary(Diff, nhh, sd), 
					  Elasticity = weight_summary(Elasticity, nhh, mean), Elast_sd = weight_summary(Elasticity, nhh, sd), 
					  MPC = weight_summary(MPC, nhh, mean), MPC_sd = weight_summary(MPC, nhh, sd))]
tmp2	<- 	tmp[, list(Exp2007 = weight_summary(Exp2007, nhh, mean), Exp2008 = weight_summary(Exp2008, nhh, mean), 
					  Difference = weight_summary(Diff, nhh, mean), Diff_sd = weight_summary(Diff, nhh, sd), 
					  Elasticity = weight_summary(Elasticity, nhh, mean), Elast_sd = weight_summary(Elasticity, nhh, sd), 
					  MPC = weight_summary(MPC, nhh, mean), MPC_sd = weight_summary(MPC, nhh, sd)), by = list(IncGrp)]
cat("Average income effect on expenditure:\n"); print(tmp1); cat("\n")
cat("Income effect on expenditure by income group:\n"); print(tmp2); cat("\n")

# ----------------------------------- #
# Income effect on shopping incidence #
sel		<- paste("Trip", fmt_name, sep="_")
tmp.dt	<- data.table(do.call(rbind, lapply(simavg.ls, function(x) x$Base08[,sel] - x$Base07[,sel])))
tmp.dt$IncGrp	<- unlist(lapply(1:nseg, function(i) rep(names(simavg.ls)[i], nrow(simavg.ls[[i]][[1]] )) ))
tmp.dt$income2007	<- round(unlist(lapply(1:nseg, function(i) simavg.ls[[i]][[1]][,"income2007"])), 2)
dim(tmp.dt)
tmp.dt		<- merge(tmp.dt, tab.inc, by = c("IncGrp", "income2007"), all.x = T)
dim(tmp.dt)
ggtmp	<- data.table(melt(tmp.dt, id.vars = c("IncGrp", "income2007", "nhh")))
tmp		<- setNames(fmt_name1, paste("Trip", fmt_name, sep="_"))
ggtmp	<- ggtmp[, retailer:= tmp[as.character(ggtmp$variable)]]
ggtmp1	<- ggtmp[,list(Difference = weight_summary(value, nhh, mean), sd = weight_summary(value, nhh, sd)), 
				 by = list(retailer)]
ggtmp2	<- ggtmp[,list(Difference = weight_summary(value, nhh, mean), sd = weight_summary(value, nhh, sd)), 
				 by = list(IncGrp, retailer)]
fmt_ord	<- ggtmp1[order(ggtmp1$Difference), retailer]
ggtmp1$retailer <- factor(ggtmp1$retailer, levels = fmt_ord)
ggtmp2$retailer <- factor(ggtmp2$retailer, levels = fmt_ord)

# Plot the incidence change
plots	<- list(NULL)
plots[[1]]	<- ggplot(ggtmp1, aes(retailer, Difference)) + 
				geom_pointrange(aes(y = Difference, ymin = Difference - 1.96*sd, ymax = Difference + 1.96*sd)) +
				xlab("Retail format") + ylab("Change in incidence probability")+ coord_flip() 
plots[[2]]	<- ggplot(ggtmp2, aes(retailer, Difference, col = IncGrp)) + 
				geom_pointrange(aes(y = Difference, ymin = Difference - 1.96*sd, ymax = Difference + 1.96*sd), position=position_dodge(width=0.3)) +
				scale_color_grey(name="Income\ngroup", start = 0, end = .6) +  
				xlab("Retail format") + ylab("Change in incidence probability")+ coord_flip() + 
				guides(color = guide_legend(reverse = TRUE))
plots[[3]]	<- get_legend(plots[[2]])
plots[[2]]	<- plots[[2]] + theme(legend.position="none")	

grid.arrange(plots[[1]], plots[[2]], plots[[3]], nrow = 1, widths = c(.45, .45, .1))
				
# ---------------------------------- #
# Income effect on expenditure share #
tmp.dt	<- data.table(do.call(rbind, lapply(simavg.ls, function(x) 
						x$Base08[,fmt_name]/rowSums(x$Base08[,fmt_name]) - x$Base07[,fmt_name]/rowSums(x$Base07[,fmt_name]) )))
tmp.dt$IncGrp	<- unlist(lapply(1:nseg, function(i) rep(names(simavg.ls)[i], nrow(simavg.ls[[i]][[1]] )) ))
tmp.dt$income2007	<- round(unlist(lapply(1:nseg, function(i) simavg.ls[[i]][[1]][,"income2007"])), 2)
dim(tmp.dt)
tmp.dt		<- merge(tmp.dt, tab.inc, by = c("IncGrp", "income2007"), all.x = T)
dim(tmp.dt)
ggtmp	<- data.table(melt(tmp.dt, id.vars = c("IncGrp", "income2007", "nhh")))
setnames(ggtmp, "variable", "retailer")
ggtmp1	<- ggtmp[,list(Difference = weight_summary(value, nhh, mean), sd = weight_summary(value, nhh, sd)), 
				 by = list(retailer)]
ggtmp2	<- ggtmp[,list(Difference = weight_summary(value, nhh, mean), sd = weight_summary(value, nhh, sd)), 
				 by = list(IncGrp, retailer)]
ggtmp1$retailer <- factor(ggtmp1$retailer, levels = fmt_ord)
ggtmp2$retailer <- factor(ggtmp2$retailer, levels = fmt_ord)

# Plot the share change
plots	<- list(NULL)
plots[[1]]	<- ggplot(ggtmp1, aes(retailer, Difference)) + 
				geom_pointrange(aes(y = Difference, ymin = Difference - 1.96*sd, ymax = Difference + 1.96*sd)) +
				geom_hline(yintercept = 0, size = .1, linetype = 2) + 
				scale_y_continuous(labels = scales::percent) + 
				xlab("Retail format") + ylab("Expenditure share change")+ coord_flip() 
plots[[2]]	<- ggplot(ggtmp2, aes(retailer, Difference, col = IncGrp)) + 
				geom_pointrange(aes(y = Difference, ymin = Difference - 1.96*sd, ymax = Difference + 1.96*sd), position=position_dodge(width=0.3)) +
				geom_hline(yintercept = 0, size = .1, linetype = 2) + 
				scale_color_grey(name="Income\ngroup", start = 0, end = .6) +  
				scale_y_continuous(labels = scales::percent) +
				xlab("Retail format") + ylab("Expenditure share change")+ coord_flip() + 
				guides(color = guide_legend(reverse = TRUE))
plots[[3]]	<- get_legend(plots[[2]])
plots[[2]]	<- plots[[2]] + theme(legend.position="none")	

grid.arrange(plots[[1]], plots[[2]], plots[[3]], nrow = 1, widths = c(.45, .45, .1))

#############################################
# Simulation results of adding a new format #
#############################################
# Evaluate the effect by focusing on the 2008 income level.
sel.scn		<- "Discount Store_Small"	#"Discount Store_prc1"	# Counterfactual scenarios we foucs on 
sel.base	<- "Base08"				# Fix baseline
scn.name	<- c("Discount Stores")

# Compute the changes in expenditure at retail formats in the simulation 
ggtmp		<- data.frame()
for(i in 1:length(sel.scn)){
	tmp		<- data.table(do.call(rbind,lapply(simavg.ls, function(x) x[[sel.scn[i]]][,fmt_name1] - x[[sel.base]][,fmt_name1])))
	tmp$IncGrp	<- unlist(lapply(1:nseg, function(i) rep(names(simavg.ls)[i], nrow(simavg.ls[[i]][[1]] )) ))
	tmp$income2007	<- round(unlist(lapply(1:nseg, function(i) simavg.ls[[i]][[1]][,"income2007"])), 2)
	tmp$scenario	<- scn.name[i]
	ggtmp	<- rbind(ggtmp, tmp)
}

ggtmp	<- melt(ggtmp, id.vars = c("scenario", "IncGrp", "income2007"))
setnames(ggtmp, "variable", "retailer")
dim(ggtmp)
ggtmp	<- merge(ggtmp, tab.inc, by = c("IncGrp", "income2007"), all.x = T)
dim(ggtmp)
ggtmp1	<- ggtmp[,list(Difference = weight_summary(value, nhh, mean), sd = weight_summary(value, nhh, sd)), 
				 by = list(retailer)]
ggtmp2	<- ggtmp[,list(Difference = weight_summary(value, nhh, mean), sd = weight_summary(value, nhh, sd)), 
				 by = list(IncGrp, retailer)]				
ggtmp1$retailer <- factor(ggtmp1$retailer, levels = fmt_name1)
ggtmp2$retailer <- factor(ggtmp2$retailer, levels = fmt_name1)

# Plot the share change
plots	<- list(NULL)
plots[[1]]	<- ggplot(ggtmp1, aes(retailer, Difference)) + 
				geom_pointrange(aes(y = Difference, ymin = Difference - 1.96*sd, ymax = Difference + 1.96*sd), size = pnt.size) +
				xlab("Retail format") + ylab("Expenditure change")+ coord_flip() 
plots[[2]]	<- ggplot(ggtmp2, aes(retailer, Difference, col = IncGrp)) + 
				geom_pointrange(aes(y = Difference, ymin = Difference - 1.96*sd, ymax = Difference + 1.96*sd), 
							position=position_dodge(width=position.width), size = pnt.size) +
				scale_color_grey(name="Income\ngroup", start = 0, end = .6) +  
				xlab("Retail format") + ylab("Expenditure change")+ coord_flip() + 
				guides(color = guide_legend(reverse = TRUE))
plots[[3]]	<- get_legend(plots[[2]])
plots[[2]]	<- plots[[2]] + theme(legend.position="none")	

pdf(paste(plot.wd, "/graph_newf_dol.pdf", sep=""), width = ww, height = ww*ar)
grid.arrange(plots[[1]], plots[[2]], plots[[3]], nrow = 1, widths = c(.45, .45, .1))
dev.off()

# ------------------------------------------------------------ #
# Compute the changes in expenditure share at 7 retail formats #
ggtmp		<- data.frame()
for(i in 1:length(sel.scn)){
	tmp		<- data.table(do.call(rbind,lapply(simavg.ls, function(x) 
						x[[sel.scn[i]]][,fmt_name1]/rowSums(x[[sel.scn[i]]][,fmt_name1]) - x[[sel.base]][,fmt_name1]/rowSums(x[[sel.base]][,fmt_name1]))))
	tmp$IncGrp	<- unlist(lapply(1:nseg, function(i) rep(names(simavg.ls)[i], nrow(simavg.ls[[i]][[1]] )) ))
	tmp$income2007	<- round(unlist(lapply(1:nseg, function(i) simavg.ls[[i]][[1]][,"income2007"])), 2)
	tmp$scenario	<- scn.name[i]
	ggtmp	<- rbind(ggtmp, tmp)
}

ggtmp	<- melt(ggtmp, id.vars = c("scenario", "IncGrp", "income2007"))
setnames(ggtmp, "variable", "retailer")
dim(ggtmp)
ggtmp	<- merge(ggtmp, tab.inc, by = c("IncGrp", "income2007"), all.x = T)
dim(ggtmp)
ggtmp1	<- ggtmp[,list(Difference = weight_summary(value, nhh, mean), 
						ymin = weight_summary(value, nhh, quantile, prob = .025), 
						ymax = weight_summary(value, nhh, quantile, prob = .975)),
				 by = list(scenario, retailer)]
ggtmp2	<- ggtmp[,list(Difference = weight_summary(value, nhh, mean),
						ymin = weight_summary(value, nhh, quantile, prob = .025), 
						ymax = weight_summary(value, nhh, quantile, prob = .975)), 
				 by = list(scenario, IncGrp, retailer)]								
ggtmp1$retailer <- factor(ggtmp1$retailer, levels = ggtmp1[order(Difference),retailer])
ggtmp2$retailer <- factor(ggtmp2$retailer, levels = ggtmp1[order(Difference),retailer])

# Plot the share change
plots	<- list(NULL)
plots[[1]]	<- ggplot(ggtmp1, aes(retailer, Difference)) + 
				geom_pointrange(aes(y = Difference, ymin = ymin, ymax = ymax), size = pnt.size) +
				facet_grid(scenario ~ ., switch = "y") + 
				scale_y_continuous(labels = scales::percent) +
				xlab("Retail format") + ylab("Change in expenditure share")+ coord_flip() 
plots[[2]]	<- ggplot(ggtmp2, aes(retailer, Difference, col = IncGrp)) + 
				geom_pointrange(aes(y = Difference, ymin = ymin, ymax = ymax), 
						position=position_dodge(width=position.width), size = pnt.size) +
				facet_grid(scenario ~ ., switch = "y") + 
				scale_color_grey(name="Income\ngroup", start = 0, end = .6) +  
				scale_y_continuous(labels = scales::percent) +
				xlab("Retail format") + ylab("Change in expenditure share")+ coord_flip() + 
				guides(color = guide_legend(reverse = TRUE))
plots[[3]]	<- get_legend(plots[[2]])
plots[[2]]	<- plots[[2]] + theme(legend.position="none")	

pdf(paste(plot.wd, "/graph_newf_shr.pdf", sep=""), width = ww, height = ww*ar)
grid.arrange(plots[[1]], plots[[2]], plots[[3]], nrow = 1, widths = c(.45, .45, .1))
dev.off()

pdf(paste(plot.wd, "/graph_newf_shr_discount.pdf", sep=""), width = 6.5, height = 6.5*.6)
# ggplot(subset(ggtmp1, scenario=="Discount Stores"), aes(retailer, Difference)) + 
# 				geom_pointrange(aes(y = Difference, ymin = ymin, ymax = ymax), size = pnt.size) +
# 				geom_hline(yintercept = 0, size = .25, linetype = 2) + 
# 				scale_y_continuous(labels = scales::percent) +
# 				xlab("Retail format") + ylab("Change in expenditure share")+ coord_flip() + 
#				theme_bw()

ggplot(subset(ggtmp1, scenario=="Discount Stores"), aes(retailer, Difference)) + 
				# geom_pointrange(aes(y = Difference, ymin = ymin, ymax = ymax), size = pnt.size) +
				geom_bar(stat = "identity") + 
				geom_hline(yintercept = 0, size = .25, linetype = 2) + 
				scale_y_continuous(labels = scales::percent) +
				xlab("Retail format") + ylab("Change in expenditure share")+ coord_flip() + 
				theme_bw()				
dev.off()

#############################################
# Interaction between income and downsizing # 
#############################################
# Evaluate the effect by focusing on the 2008 income level.
sel.scn0	<- c("Discount Store_prc1_07")			# Counterfactual scenarios we foucs on in 2007
sel.scn1	<- c("Discount Store_prc1")			# Counterfactual scenarios we foucs on 
sel.base0	<- "Base07"
sel.base1	<- "Base08"				# Fix baseline
scn.name	<- c("Discount Stores", "Grocery")

# Compute the changes in expenditure share at retail formats in the simulation 
inter.only	<- TRUE
ggtmp		<- data.frame()
for(i in 1:length(sel.scn0)){
	# The effect of downsizing at 2007 income
	tmp0	<- data.table(do.call(rbind,lapply(simavg.ls, function(x) 
				x[[sel.scn0[i]]][,fmt_name1]/rowSums(x[[sel.scn0[i]]][,fmt_name1]) - x[[sel.base0]][,fmt_name1]/rowSums(x[[sel.base0]][,fmt_name1]))))
	# The effect of downsizing at 2008 income 
	tmp1	<- data.table(do.call(rbind,lapply(simavg.ls, function(x) 
				x[[sel.scn1[i]]][,fmt_name1]/rowSums(x[[sel.scn1[i]]][,fmt_name1]) - x[[sel.base1]][,fmt_name1]/rowSums(x[[sel.base1]][,fmt_name1]))))			
	if(inter.only){
		tmp	<- tmp1 - tmp0
		tmp$IncGrp	<- unlist(lapply(1:nseg, function(i) rep(names(simavg.ls)[i], nrow(simavg.ls[[i]][[1]] )) ))
		tmp$income2007	<- round(unlist(lapply(1:nseg, function(i) simavg.ls[[i]][[1]][,"income2007"])), 2)
		tmp$scenario	<- scn.name[i]
		tmp$condition	<- "Interaction"
		ggtmp	<- rbind(ggtmp, tmp)
	}else{
		tmp0$IncGrp	<- tmp1$IncGrp	<- unlist(lapply(1:nseg, function(i) rep(names(simavg.ls)[i], nrow(simavg.ls[[i]][[1]] )) ))
		tmp0$income2007	<- tmp1$income2007	<- round(unlist(lapply(1:nseg, function(i) simavg.ls[[i]][[1]][,"income2007"])), 2)
		tmp0$scenario	<- tmp1$scenario 	<- scn.name[i]
		tmp0$condition	<- "2007"
		tmp1$condition	<- "2008"
		ggtmp	<- rbind(ggtmp, tmp0, tmp1)
	}
}

ggtmp	<- melt(ggtmp, id.vars = c("condition", "scenario", "IncGrp", "income2007"))
setnames(ggtmp, "variable", "retailer")
dim(ggtmp)
ggtmp	<- merge(ggtmp, tab.inc, by = c("IncGrp", "income2007"), all.x = T)
dim(ggtmp)
ggtmp1	<- ggtmp[,list(Difference = weight_summary(value, nhh, mean), 
						ymin = weight_summary(value, nhh, quantile, prob = .025), 
						ymax = weight_summary(value, nhh, quantile, prob = .975)),
				 by = list(condition, scenario, retailer)]
ggtmp2	<- ggtmp[,list(Difference = weight_summary(value, nhh, mean),
						ymin = weight_summary(value, nhh, quantile, prob = .025), 
						ymax = weight_summary(value, nhh, quantile, prob = .975)), 
				 by = list(condition, scenario, IncGrp, retailer)]				
ggtmp1$retailer <- factor(ggtmp1$retailer, levels = fmt_name1)
ggtmp2$retailer <- factor(ggtmp2$retailer, levels = fmt_name1)


# Plot the share change
plots	<- list(NULL)
plots[[1]]	<- ggplot(ggtmp1, aes(retailer, Difference, shape = condition)) + 
				geom_pointrange(aes(y = Difference, ymin = ymin, ymax = ymax), size = pnt.size, position = "dodge") +
				facet_grid(scenario ~ ., switch = "y") + 
				scale_y_continuous(labels = scales::percent) +
				xlab("Retail format") + ylab("Change in expenditure share")+ coord_flip() 
plots[[2]]	<- ggplot(ggtmp2, aes(retailer, Difference, col = IncGrp, shape = condition)) + 
				geom_pointrange(aes(y = Difference, ymin = ymin, ymax = ymax), 
						position=position_dodge(width=position.width), size = pnt.size) +
				facet_grid(scenario ~ ., switch = "y") + 
				scale_color_grey(name="Income\ngroup", start = 0, end = .6) +  
				scale_y_continuous(labels = scales::percent) +
				xlab("Retail format") + ylab("Change in expenditure share")+ coord_flip() + 
				guides(color = guide_legend(reverse = TRUE))
plots[[3]]	<- get_legend(plots[[2]])
plots[[2]]	<- plots[[2]] + theme(legend.position="none")	

pdf(paste(plot.wd, "/graph_newf_shr_interaction.pdf", sep=""), width = ww, height = ww*ar)
grid.arrange(plots[[1]], plots[[2]], plots[[3]], nrow = 1, widths = c(.45, .45, .1))
dev.off()

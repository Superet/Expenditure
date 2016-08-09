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

# plot.wd		<- "/Users/chaoqunchen/Desktop"

# setwd("/home/brgordon/ccv103/Exercise/run")
run_id		<- 11
ww			<- 6.5
ar			<- .6
ww1			<- 8.5
ar1			<- .5
# setwd(paste("~/Documents/Research/Store switching/processed data/Estimation/estrun_", run_id, sep=""))
setwd(paste("U:/Users/ccv103/Documents/Research/Store switching/run/estrun_", run_id, sep=""))

# Set plotting parameters
annual.month	<- 12
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
numsim		<- 500
ver.date1	<- "2016-07-24"  
tmpls	<- ls()

for(ii in 1:nseg){
	# Load data 
	fname			<- paste("marginal_seg",ii, sep="")
	if(draw.par){ fname <- paste(fname, "_pardraw", sep="") }
	if(sim.omega){fname	<- paste(fname, "_simomega", sep="") }
	fname			<- paste(fname, "_sim", numsim, "_", ver.date1, sep="")
	load(paste(fname, ".rdata",sep=""))
	
	tmp1	<- sim.base07$Average[1,]
	tmp2	<- sim.base08$Average[1,]
	sim07			<- do.call(rbind, lapply(sim.ls07, function(x) x$Average[1,]))
	names(sim07)	<- gsub("\\.", " ", names(sim07))
	sim08			<- do.call(rbind, lapply(sim.ls08, function(x) x$Average[1,]))
	names(sim08)	<- gsub("\\.", " ", names(sim08))
	
	# If we want to convert the biweekly simulation to annual level
	if(annual){
		sel	<- c("Exp", fmt_name)
		tmp1[sel]	<- tmp1[sel] * annual.month
		tmp2[sel]	<- tmp2[sel] * annual.month
		sim07[,sel]	<- sim07[,sel] * annual.month
		sim08[,sel]	<- sim08[,sel] * annual.month
	}
		
	sim.ls[[ii]]<- list(Base07 = tmp1, Base08 = tmp2, sim07 = sim07, sim08 = sim08)
						
	my.simdata	<- rbind(my.simdata, data.frame(sim.data, IncGrp = names(sim.ls)[ii]))
	
	# Delete unused objects
	rm(list = setdiff(ls(), c(tmpls, "tmpls", "ii", "fmt_name", "R", "iv.change.vec")))
	print(ii)						
}

# Segment weight
iv.change <- .1
selcol1	<- selcol	<- c("size_index", "ln_upc_per_mod", "ln_num_module","overall_prvt","price")
wt.data	<- table(my.simdata$IncGrp)
# names(wt.data)	<- c("IncGrp", "nhh")

# Combine the simulation from 3 segments
base07	<- do.call(rbind, lapply(sim.ls, function(x) x$Base07))
base08	<- do.call(rbind, lapply(sim.ls, function(x) x$Base08))
sim07	<- do.call(rbind, lapply(1:nseg, function(i) cbind(sim.ls[[i]]$sim07, IncGrp = names(sim.ls)[i])))
sim08	<- do.call(rbind, lapply(1:nseg, function(i) cbind(sim.ls[[i]]$sim08, IncGrp = names(sim.ls)[i])))
tmp07	<- base07[,fmt_name] * (wt.data %*% matrix(1, R, nrow =1))
tmp08	<- base08[,fmt_name] * (wt.data %*% matrix(1, R, nrow =1))

ggtmp1	<- matrix(NA, 2, length(selcol1), dimnames = list(c("Before", "After"), selcol1))
ggtmp2	<- matrix(NA, 2*nseg, length(selcol1), 
				dimnames = list(paste(rep(c("Before", "After"),each = nseg), segname, sep="-") , selcol1))
sel.retail	<- "Discount Store"
dv.shr	<- TRUE

for(i in 1:length(selcol1)){
	sel <- sim07$Var == selcol1[i]
	tmp1<- sim07[sel,fmt_name] * (wt.data %*% matrix(1, R, nrow =1))
	tmp2<- sim08[sel,fmt_name] * (wt.data %*% matrix(1, R, nrow =1))
	
	# Overall marginal effect 
	# The dependent variable of elasticity
	if(dv.shr){
		tmp10	<- colSums(tmp07)/sum(tmp07)
		tmp20	<- colSums(tmp08)/sum(tmp08)
		tmp11	<- colSums(tmp1)/sum(tmp1)
		tmp21	<- colSums(tmp2)/sum(tmp2)
		ggtmp1[,i]	<- c( (tmp11[sel.retail] - tmp10[sel.retail])/tmp10[sel.retail]/iv.change , 	# in 2007
				 		  (tmp21[sel.retail] - tmp20[sel.retail])/tmp20[sel.retail]/iv.change )		# in 2008
	}else{
		ggtmp1[,i]<- c((sum(tmp1[,sel.retail]) - sum(tmp07[,sel.retail]))/sum(tmp07[,sel.retail])/iv.change , 	# in 2007
				 		(sum(tmp2[,sel.retail]) - sum(tmp08[,sel.retail]))/sum(tmp08[,sel.retail])/iv.change )		# in 2008
	}
	
	# Marginal effect by segment
	if(dv.shr){
		tmp1<- (tmp1/rowSums(tmp1) - tmp07/rowSums(tmp07))/(tmp07/rowSums(tmp07))/iv.change
		tmp2<- (tmp2/rowSums(tmp2) - tmp08/rowSums(tmp08))/(tmp08/rowSums(tmp08))/iv.change
		ggtmp2[,i]	<- c(tmp1[,sel.retail], tmp2[,sel.retail])
	}else{
		tmp1<- (tmp1 - tmp07)/tmp07/iv.change
		tmp2<- (tmp2 - tmp08)/tmp08/iv.change
		ggtmp2[,i]	<- c(tmp1[,sel.retail], tmp2[,sel.retail])		
	}
}

cat("Overall marginal effect:\n"); print(ggtmp1); cat("\n")
cat("Marginal effect by segment:\n"); print(ggtmp2); cat("\n")


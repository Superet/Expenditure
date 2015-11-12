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

# Read estimation from all segments #
nseg	<- 3
sim.ls	<- setNames(vector("list", length = nseg), c("Low", "Med", "High"))
my.simdata	<- data.frame()

for(ii in 1:nseg){
	load(paste("Estimation/ctrfact_seg",ii,".rdata",sep=""))
	tmpn		<- setNames(1:nrow(sim.unq), sim.unq$income2007)
	sel			<- tmpn[as.character(sim.data$income)]
	sim.ls[[ii]]<- list(	Base07 	= sim.base07[sel,], 
							Base08 	= sim.base08[sel,], 
						    Price 	= lapply(sim.prc.ls, function(m) m[sel,]))
	tmp			<- lapply(sim.X.ls, function(l) lapply(l, function(m) m[sel,]))	
	sim.ls[[ii]]<- c(sim.lis[[ii]], tmp) 				
	my.simdata	<- rbind(my.simdata, data.frame(sim.data, IncGrp = names(sim.ls)[ii]))
	rm(list = c("sim.base07", "sim.base08", "sim.fixy08", "sim.fixp08"))						
}

######################
# Plot income effect #
######################




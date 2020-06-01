library(ggplot2)
library(reshape2)
library(maxLik)
library(evd)
library(data.table)
library(doParallel)
library(foreach)
# library(chebpol)
library(nloptr)
library(mgcv)
options(error = quote({dump.frames(to.file = TRUE)}))

# setwd("~/Documents/Research/Store switching/Processed_data/processed data_20160809")
# plot.wd	<- '~/Desktop'
# sourceCpp(paste("../../Exercise/Multiple_discrete_continuous_model/", model_name, ".cpp", sep=""))
# source("../../Exercise/Multiple_discrete_continuous_model/0_Efficient_Allocation_function.R")
# source("../../Exercise/main/ctrfact_sim_functions.r")

# setwd("/home/brgordon/ccv103/Exercise/run")
# setwd("/kellogg/users/marketing/2661703/Exercise/run")
# setwd("/sscc/home/c/ccv103/Exercise/run")
setwd("U:/Users/ccv103/Documents/Research/Store switching/run")

run_id		<- 2
plot.wd		<- getwd()
make_plot	<- TRUE
ww			<- 10
ar			<- .6

# Set plotting parameters
annual.month	<- 12
annual		<- TRUE
draw.par		<- FALSE
sim.omega		<- FALSE

# Load estimation data 
ver.date	<- "2016-08-28" # "2016-08-22"
cpi.adj		<- TRUE
setwd(paste("estrun_", run_id, sep=""))

###########################
# Read counterfactal data # 
###########################
# Read estimation from all segments #
nseg	<- 3
segname	<- c("Low", "Med", "High")
sim.ls	<- setNames(vector("list", length = nseg), c("Low", "Med", "High"))
my.simdata	<- data.frame()
numsim		<- 500
ver.date1	<- "2016-08-29"  
base.07	<- base.08	<- sim07	<- sim08	<- data.frame()
tmpls	<- ls()

for(ii in 1:nseg){
	# Load data 
	fname			<- paste("marginal_seg",ii, sep="")
	fname			<- paste(fname, "_sim", numsim, "_", ver.date1, sep="")
	load(paste(fname, ".rdata",sep=""))
	
	tmp1	<- sim.base07$Average
	tmp2	<- sim.base08$Average
	# tmp3	<- do.call(rbind, lapply(sim.2ls07, function(x) x$Average))
	# tmp4	<- do.call(rbind, lapply(sim.2ls08, function(x) x$Average))
	tmp3	<- do.call(rbind, lapply(sim.ls07, function(x) x$Average))
	tmp4	<- do.call(rbind, lapply(sim.ls08, function(x) x$Average))
	names(tmp3)	<- gsub("\\.", " ", names(tmp3))
	names(tmp4)	<- gsub("\\.", " ", names(tmp4))
	
	# If we want to convert the biweekly simulation to annual level
	if(annual){
		sel	<- c("Exp", fmt_name)
		tmp1[,sel]	<- tmp1[,sel] * annual.month
		tmp2[,sel]	<- tmp2[,sel] * annual.month
		tmp3[,sel]	<- tmp3[,sel] * annual.month
		tmp4[,sel]	<- tmp4[,sel] * annual.month
	}
	base.07	<- rbind(base.07, tmp1)	
	base.08	<- rbind(base.08, tmp2)	
	sim07	<- rbind(sim07, tmp3)
	sim08	<- rbind(sim08, tmp4)		
	
						
	my.simdata	<- rbind(my.simdata, data.frame(sim.unq[,c("household_code", "ln_inc", "price")], IncGrp = segname[ii]))
	
	# Delete unused objects
	rm(list = setdiff(ls(), c(tmpls, "tmpls", "ii", "fmt_name", "R", "iv.change")))
	print(ii)						
}

# Mean elasticity
sel.retail	<- "Discount Store"

tmp1	<- split(sim07, sim07$Var)
tmp1	<- lapply(tmp1, function(x) (x[,sel.retail] - base.07[,sel.retail])/base.07[,sel.retail]/x[,"change"])
tmp2	<- split(sim08, sim08$Var)
tmp2	<- lapply(tmp2, function(x) (x[,sel.retail] - base.08[,sel.retail])/base.08[,sel.retail]/x[,"change"])
tmp		<- lapply(1:length(tmp1), function(i) t.test(tmp2[[i]], tmp1[[i]]))
tmp		<- do.call(rbind, lapply(tmp, function(x) c(x$estimate, x$p.value)))
colnames(tmp)	<- c("2008", "2007", "pvalue")
rownames(tmp)	<- names(tmp1)
cat("T test of mean elasticity:\n"); print(tmp); cat("\n")

# Overall Revenue
tmp0	<- colSums(base.07[,fmt_name])
tmp1	<- split(sim07, sim07$Var)
tmp1	<- lapply(tmp1, function(x) colSums(x[,fmt_name]))
tmp1	<- sapply(tmp1, function(x) (x[sel.retail] - tmp0[sel.retail])/tmp0[sel.retail]/iv.change)
tmp0	<- colSums(base.08[,fmt_name])
tmp2	<- split(sim08, sim08$Var)
tmp2	<- lapply(tmp2, function(x) colSums(x[,fmt_name]))
tmp2	<- sapply(tmp2, function(x) (x[sel.retail] - tmp0[sel.retail])/tmp0[sel.retail]/iv.change)
tmp.tab	<- cbind(tmp1, tmp2)
colnames(tmp.tab)	<- c("2007", "2008")
cat("Elasticity on overall revenue:\n"); print(tmp.tab); cat("\n")

# Market share
tmp0	<- colSums(base.07[,fmt_name])
tmp0	<- tmp0/sum(tmp0)
tmp1	<- split(sim07, sim07$Var)
tmp1	<- lapply(tmp1, function(x) colSums(x[,fmt_name]))
tmp1	<- lapply(tmp1, function(x) x/sum(x))
tmp1	<- sapply(tmp1, function(x) (x[sel.retail] - tmp0[sel.retail])/tmp0[sel.retail]/iv.change)
tmp0	<- colSums(base.08[,fmt_name])
tmp0	<- tmp0/sum(tmp0)
tmp2	<- split(sim08, sim08$Var)
tmp2	<- lapply(tmp2, function(x) colSums(x[,fmt_name]))
tmp2	<- lapply(tmp2, function(x) x/sum(x))
tmp2	<- sapply(tmp2, function(x) (x[sel.retail] - tmp0[sel.retail])/tmp0[sel.retail]/iv.change)
tmp.tab1	<- cbind(tmp1, tmp2)
colnames(tmp.tab1)	<- c("2007", "2008")
cat("Elasticity on overall revenue:\n"); print(tmp.tab); cat("\n")

# Distance marginal effect
setwd("U:/Users/ccv103/Documents/Research/Store switching/run")
load("hh_fmt.rdata")

max.dist	<- 100
sum(hh_dist$distance_haver >max.dist|is.na(hh_dist$distance_haver), na.rm=T)/nrow(hh_dist)
sel			<- hh_dist$distance_haver >max.dist|is.na(hh_dist$distance_haver)
hh_dist[sel,"distance_haver"]	<- max.dist

# Normalize distance by zipcode
tmp		<- data.table(hh_dist[,c("household_code", "panelist_zip_code","year", "channel_type", "distance_haver")])
tmp		<- tmp[, maxd:=max(distance_haver, na.rm=T), by = list(panelist_zip_code)]
tmp		<- tmp[, distance_haver := distance_haver/maxd]
tmp		<- data.frame(tmp)

tmp1	<- merge(my.simdata, subset(tmp, year == 2007 & channel_type == sel.retail), by = "household_code", all.x=T)
ggtmp	<- cbind(tmp1[,c("household_code", "distance_haver", "maxd")], 
				old = base.07[,sel.retail], new = subset(sim07, Var == "distance_haver")[,sel.retail]
				)
ggtmp$d_d	<- with(ggtmp, -.1*distance_haver*maxd)
ggtmp$d_exp	<- with(ggtmp, new - old)
summary(ggtmp$d_d)
summary(ggtmp$d_exp/ggtmp$d_d)
summary(with(ggtmp, d_exp/old/(d_d/(maxd*distance_haver)) ))

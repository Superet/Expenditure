library(ggplot2)
library(reshape2)
library(maxLik)
library(evd)
library(data.table)
library(doParallel)
library(foreach)
library(nloptr)
library(mgcv)
library(scales)
options(error = quote({dump.frames(to.file = TRUE)}))


# setwd("~/Documents/Research/Store switching/processed data")
# plot.wd	<- '~/Desktop'
# source("../Exercise/Multiple_discrete_continuous_model/0_Allocation_function.R")
# source("../Exercise/main/share_allocation/ctrfact_sim_functions_v2.r")

# setwd("/home/brgordon/ccv103/Exercise/run")
# setwd("/kellogg/users/marketing/2661703/Exercise/run")
setwd("/sscc/home/c/ccv103/Exercise/run")
run_id		<- 4
plot.wd		<- getwd()
make_plot	<- TRUE
ww			<- 6.5
ar			<- .6

# outfile 	<- paste(plot.wd, "/counter1.xlsx", sep="")
# mywb		<- createWorkbook()
# sht1		<- createSheet(mywb, "Counterfactual1")
# make_plot	<- TRUE

# Set plotting parameters
draw.par		<- FALSE
sim.omega		<- FALSE

###########################
# Read counterfactal data # 
###########################
# Read estimation from all segments #
nseg	<- 3
sim.ls	<- setNames(vector("list", length = nseg), c("Low", "Med", "High"))
sim.dat	<- data.frame()
numsim		<- 1000
ver.date1	<- "2016-03-09" #"2016-01-25"
tmpls	<- ls()

for(ii in 1:nseg){
	fname			<- paste("ctrfact_dcps_seg",ii,sep="")
	if(draw.par){ fname <- paste(fname, "_pardraw", sep="") }
	if(sim.omega){fname	<- paste(fname, "_simomega", sep="") }
	(fname			<- paste(fname, "_sim", numsim, "_", ver.date1, ".rdata", sep=""))
	load(paste("estrun_", run_id, "/",fname, sep=""))
	sim.ls[[ii]]	<- rbind(data.frame(sim0, scenario = "Baseline"), 
							 data.frame(sim1, scenario = "FixIncome"), 
							 data.frame(sim2, scenario = "FixAttr"))
	sim.dat 		<- rbind(sim.dat, data.frame(sim.data, IncGrp = names(sim.ls)[ii]))
	rm(list = setdiff(ls(), c(tmpls, "tmpls", "R", "fmt_name", "ii")))
	print(ii)
}

# Load the complete data 
load("hh_biweek_exp.rdata")

# ------------------------------------------------------ #
# Segment households based on their initial income level #
panelist	<- data.table(hh_exp)
setkeyv(panelist, c("household_code","year","biweek"))
panelist	<- panelist[,list(income = first_income[1], first_famsize = famsize[1]), by=list(household_code)]
tmp			<- quantile(panelist$income, c(0, .33, .67, 1))
num_grp		<- 3
# panelist	<- panelist[, first_incomeg := cut(panelist$income, tmp, labels = paste("T", 1:num_grp, sep=""), include.lowest = T)]
panelist	<- panelist[, first_incomeg := cut(panelist$income, tmp, labels = c("Low", "Med", "High"), include.lowest = T)]
hh_exp		<- merge(hh_exp, data.frame(panelist)[,c("household_code", "first_incomeg")], by = "household_code", all.x=T )
cat("Table of initial income distribution:\n"); print(table(panelist$first_incomeg)); cat("\n")
cat("Table of segments in the expenditure data:\n"); print(table(hh_exp$first_incomeg)); cat("\n")

#############################
# Plot the simulation trend # 
#############################
# Compute the expenditure share using the actual data # 
ggtmp	<- data.table(hh_exp[,c("biweek", "dol", paste("DOL_", gsub("\\s", "_", fmt_name), sep=""))])
ggtmp	<- ggtmp[,list(dol=sum(dol), 
					   DOL_Convenience_Store = sum(DOL_Convenience_Store), DOL_Discount_Store = sum(DOL_Discount_Store), 
					   DOL_Dollar_Store = sum(DOL_Dollar_Store), DOL_Drug_Store = sum(DOL_Drug_Store), 
					   DOL_Grocery = sum(DOL_Grocery), DOL_Warehouse_Club = sum(DOL_Warehouse_Club)), 
					by = list(biweek)]
for(i in paste("DOL_", gsub("\\s", "_", fmt_name), sep="")){
	ggtmp	<- ggtmp[,eval(as.name(i)):= eval(as.name(i))/dol]
}
ggtmp	<- melt(ggtmp[,!"dol", with = FALSE], id.vars=c("biweek"))
ggtmp$variable <- factor(ggtmp$variable, levels = paste("DOL_", gsub("\\s", "_", fmt_name), sep=""), labels = fmt_name)
ggtmp$biweek	<- as.Date("2003-12-30") + ggtmp$biweek*14



	print(ggplot(subset(ggtmp, variable != "Overall expenditure"), aes(biweek, value)) + 
			stat_smooth(method = "gam", formular = y ~ s(x, sp = .05), col = "black") + 
			facet_wrap(~ variable, ncol = 2, scales = "free_y") + 
			scale_y_continuous(labels=percent) + 
			labs(x = "Biweek", y = "Expenditure share") + 
			theme_bw()
		)
		


# --------------------------------------------------------- #
# Compute the expenditure share in counterfactual scenarios #
tmp		<- lapply(1:nseg, function(i) cbind(sim.ls[[i]][sim.ls[[i]]$scenario %in% c("Baseline", "FixIncome"),], IncGrp = names(sim.ls)[i]))
ggtmp1	<- data.table(do.call(rbind, tmp)) 
ggtmp1	<- ggtmp1[,list(Convenience_Store = sum(Convenience.Store), Discount_Store = sum(Discount.Store), 
   					Dollar_Store = sum(Dollar.Store), Drug_Store = sum(Drug.Store), 
   					Grocery = sum(Grocery), Warehouse_Club = sum(Warehouse.Club)),
				by = list(scenario, year)]	
ggtmp1	<- ggtmp1[,dol:=Convenience_Store + Discount_Store+ Dollar_Store + Drug_Store + Grocery + Warehouse_Club]				
for(i in gsub("\\s", "_", fmt_name)){
	ggtmp1	<- ggtmp1[,eval(as.name(i)):= eval(as.name(i))/dol]
}									
ggtmp1	<- melt(ggtmp1[,!"dol",with =FALSE], id.vars = c("scenario", "year"))

# Combine together the actual annual market share
tmp	<- subset(hh_exp, household_code %in% sim.dat$household_code)
tmp	<- data.table(tmp[,c("year", "dol", paste("DOL_", gsub("\\s", "_", fmt_name), sep=""))])
tmp	<- tmp[,list(dol=sum(dol), 
					   Convenience_Store = sum(DOL_Convenience_Store), Discount_Store = sum(DOL_Discount_Store), 
					   Dollar_Store = sum(DOL_Dollar_Store), Drug_Store = sum(DOL_Drug_Store), 
					   Grocery = sum(DOL_Grocery), Warehouse_Club = sum(DOL_Warehouse_Club)), 
					by = list(year)]
for(i in gsub("\\s", "_", fmt_name)){
	tmp	<- tmp[,eval(as.name(i)):= eval(as.name(i))/dol]
}
tmp	<- melt(tmp[,!"dol", with = FALSE], id.vars=c("year"))
tmp$scenario	 <- "Actual"
setcolorder(tmp, names(ggtmp1))
ggtmp2	<- rbind(ggtmp1, tmp)

ggplot(ggtmp2, aes(year, value, col = scenario)) + geom_line() + 
		facet_wrap(~ variable, ncol = 2, scales = "free_y") + 
		scale_y_continuous(labels=percent) + 
		labs(x = "Biweek", y = "Expenditure share") + 
		theme_bw()


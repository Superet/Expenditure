library(ggplot2)
library(reshape2)
library(data.table)
library(plm)
library(gridExtra)
library(scales)
library(systemfit)
library(stargazer)
library(gridExtra)
library(mvProbit)
library(mgcv)

setwd("~/Documents/Research/Store switching/processed data")
plot.wd	<- '/Users/chaoqunchen/Desktop'
source('~/Documents/Research/Store switching/Exercise/main/income_effect/plm_vcovHC.R')

week.price	<- FALSE
cpi.adj		<- TRUE
write2csv	<- FALSE
make_plot	<- TRUE
fname		<- "expenditure_reg_sub"
if(cpi.adj) { fname <- paste(fname, "_cpi", sep="")}
if(week.price){ fname <- paste(fname, "_wkprc", sep="")}

if(week.price){
	load("hh_biweek_exp_20150812.rdata")
}else{
	load("hh_biweek_exp.rdata")
}
codebook	<- read.csv("code_book.csv")
pan_org		<- read.csv("2_panelists_20160606.csv")

#################
# Organize data # 
#################
# Re-code demographic data
sel			<- grep("panel_year", names(pan_org))
names(pan_org)[sel]	<- "year"
age.code	<- setNames(c(0, 22, 27, 32, 37, 42, 47, 52, 60, 65), 
						0:9)
pan_org$female_head_age	<- age.code[as.character(pan_org$female_head_age)]						
pan_org$male_head_age	<- age.code[as.character(pan_org$male_head_age)]
pan_org$age <- rowMeans(cbind(ifelse(pan_org$female_head_age==0, NA, pan_org$female_head_age), 
						ifelse(pan_org$male_head_age==0, NA, pan_org$male_head_age)),na.rm=T)
pan_org$retire	<- ifelse((pan_org$female_head_occupation==12| pan_org$male_head_occupation ==12), 1, 0)						

# ----------------------------------------------- #
# Collaspse biweekly data to household-year panel #
# Retail formats 
fmt_name 	<- as.character(sort(unique(fmt_attr$channel_type)))
R			<- length(fmt_name)
hh_exp$first_incomeg	<- as.character(hh_exp$first_incomeg)
hh_exp$first_incomeg	<- factor(hh_exp$first_incomeg, levels = paste("T", 1:3, sep=""), labels = c("Low", "Med", "High"))
cat("Table of segments in the expenditure data:\n"); print(table(hh_exp$first_incomeg)); cat("\n")

pan_yr		<- data.table(hh_exp)
pan_yr		<- pan_yr[,list(Income = unique(income_midvalue), dol = sum(dol), 
							DOL_Grocery = sum(DOL_Grocery), DOL_Discount_Store = sum(DOL_Discount_Store), 
							DOL_Warehouse_Club = sum(DOL_Warehouse_Club), DOL_Dollar_Store = sum(DOL_Dollar_Store), 
							DOL_Drug_Store = sum(DOL_Drug_Store), DOL_Convenience_Store = sum(DOL_Convenience_Store)), 
						by = list(first_incomeg, household_code, year, cpi)]
setkeyv(pan_yr, c("household_code", "year"))
pan_yr		<- pan_yr[, ':='(ln_income = log(Income), recession = 1*(year >= 2008))]
pan_yr    	<- pan_yr[,':='(first_income = Income[1], year_id= year - year[1] + 1, tenure = length(year), move = c(0, 1*(diff(Income)!=0)))
                    , by = list(household_code)]
pan_yr		<- pan_yr[,move_all := sum(move), by = list(household_code)]
# pan_yr 	<- pan_yr[, last.income := Income[length(Income)], by = list(household_code)]
# pan_yr$last_incomeg	<- cut(pan_yr$last.income, c(0, 37500, 65000, 200000), labels = c("Low", "Med", "High"))

# Compute expenditure share
pan_yr		<- data.frame(pan_yr)
sel			<- paste("DOL_", gsub("\\s", "_", fmt_name), sep="")
sel1		<- gsub("DOL", "SHR", sel)
for(i in 1:length(fmt_name)){
	pan_yr[,sel1[i]] <- pan_yr[,sel[i]]/pan_yr$dol*100			# Percentage scale
}

# Merge in demographics 
dim(pan_yr)
pan_yr		<- merge(pan_yr, pan_org[,c("household_code", "year", "age", "household_size", "NumChild", "retire")],
					by = c("household_code","year"),all.x =T)					
dim(pan_yr)					
pan_yr		<- data.table(pan_yr)
setkeyv(pan_yr, c("first_incomeg", "household_code","year"))

# ------------------------------- #
# Collapse panelist characristics #
panelist	<- pan_yr[, list(tenure = unique(tenure), large.init.change = max(abs(Income-Income[1])/Income[1] ), 
							delta.lastincome = Income[length(Income)] - Income[(length(Income)-1)], 
							growth = mean(diff(Income)/Income[-length(Income)], na.rm=TRUE), 
							age = age[1], retire = 1*any(retire==1) ), 
						by = list(first_incomeg, household_code)]
summary(panelist)

# -------------------------------------------- #
# Difference data relative to the initial year #
selcol		<- c("dol", sel1)
pan_yrd		<- data.table(pan_yr)
inc.scale	<- 10000
large.shock	<- 30000				# Threshold for a large income change

for(i in selcol){
	v		<- as.name(paste("LAG_", i,sep=""))
	pan_yrd	<- pan_yrd[, eval(v):=c(NA, eval(as.name(i))[-length(year)] ), by = list(household_code) ]
	pan_yrd	<- pan_yrd[, eval(as.name(i)):=eval(as.name(i)) - eval(as.name(i))[1], by = list(household_code)]
}
if(cpi.adj){
	pan_yrd		<- pan_yrd[, ':='(d_Income = Income - Income[1]), by = list(household_code)]
}else{
	pan_yrd		<- pan_yrd[, ':='(d_Income = Income/cpi - Income[1]/cpi[1]), by = list(household_code)]
}
pan_yrd		<- pan_yrd[, ':='(d_Income = d_Income/inc.scale, first_income=first_income/inc.scale)]
pan_yrd		<- pan_yrd[, ':='(large = factor(1*(abs(d_Income)> large.shock/inc.scale)))]
# pan_yrd		<- pan_yrd[, ':='(large = factor(1*(any(abs(d_Income)>large.shock/inc.scale)))), by = list(household_code)]
pan_yrd		<- pan_yrd[, ':='(retire_age = 1*(age[1]>=60),
 							 ADP_Convenience_Store 	= 1*(DOL_Convenience_Store[1]>0), 
							ADP_Discount_Store 		= 1*(DOL_Discount_Store[1]>0), 
							ADP_Dollar_Store 		= 1*(DOL_Dollar_Store[1]>0),
							ADP_Drug_Store 			= 1*(DOL_Drug_Store[1]>0), 
							ADP_Grocery 			= 1*(DOL_Grocery[1]>0), 
							ADP_Warehouse_Club 		= 1*(DOL_Warehouse_Club[1]>0) ), 
						by = list(household_code)]
pan_yrd		<- pan_yrd[,':='(year = factor(year), retire_age = factor(retire_age), retire = factor(retire),
							ADP_Convenience_Store = factor(ADP_Convenience_Store), ADP_Discount_Store = factor(ADP_Discount_Store), 
							ADP_Dollar_Store = factor(ADP_Dollar_Store), ADP_Drug_Store = factor(ADP_Drug_Store), 
							ADP_Grocery = factor(ADP_Grocery), ADP_Warehouse_Club = factor(ADP_Warehouse_Club))]						
table(pan_yrd$large)

###################
# Run regressions # 
###################
dv.vec	<- c("dol", paste("SHR_", gsub("\\s", "_", fmt_name), sep=""))
dv.lab	<- c("Expenditure ($)", paste("Expenditure share at ", fmt_name, " (%)", sep=""))
sel 	<- pan_yrd$year_id > 1
# sel		<- rep(TRUE, nrow(pan_yrd))
# fml.ls	<- list(IncGrp 		= y ~ d_Income:first_incomeg + year + retire + household_size, 
# 				ShockSize	= y ~ d_Income:first_incomeg:large + year + retire + household_size, 
# 				RetireAge	= y ~ d_Income:first_incomeg:retire_age + year + retire + household_size, 
# 				Adoption	= y ~ d_Income:first_incomeg:ADP + year + retire + household_size)
fml.ls	<- list(Inc			= y ~ d_Income + year + retire + household_size,
				IncGrp 		= y ~ d_Income:first_incomeg + year + retire + household_size, 
				YearGrp		= y ~ d_Income:first_incomeg + year:first_incomeg + retire + household_size, 
				SplIncome	= y ~ s(first_income, by = d_Income) + year + retire + household_size )
				# SplYear		= y ~ s(first_income, by = year) + d_Income + retire + household_size )
				

fit.ls	<- setNames(vector("list", length(fml.ls)), names(fml.ls))				
for(i in 1:length(fml.ls)){
	fit.ls[[i]]	<- setNames(vector("list", length(dv.vec)), dv.vec)
	for(j in 1:length(dv.vec)){
		if(names(fml.ls)[i] == "Adoption"){
			if(j >1){
				fml	<- as.formula(paste(dv.vec[j], "~", gsub("ADP", gsub("SHR", "ADP", dv.vec[j]), as.character(fml.ls[[i]]))[3], "+LAG_",dv.vec[j], sep=""))
			}else{
				# fml	<- as.formula(substitute(y ~ x+ LAG_y, list(y = as.name(dv.vec[j]), x = terms(fml.ls[[1]])[[3]])) )
				fml		<- as.formula(paste(dv.vec[j], "~", as.character(fml.ls[[1]])[3], "+LAG_",dv.vec[j], sep=""))
			}
		}else{
			# fml	<- as.formula(substitute(y ~ x, list(y = as.name(dv.vec[j]), x = terms(fml.ls[[i]])[[3]])) )
			fml		<- as.formula(paste(dv.vec[j], "~", as.character(fml.ls[[i]])[3], "+LAG_",dv.vec[j], sep=""))
		}
		if(length(grep("Spl", names(fml.ls)[i])) == 0){
			fit.ls[[i]][[j]]	<- lm(fml, data = pan_yrd[sel,])
		}else{
			fit.ls[[i]][[j]]	<- gam(fml, data = pan_yrd[sel,])
			if(fit.ls[[i]][[j]]$sp < .1){
				fit.ls[[i]][[j]]	<- gam(fml, data = pan_yrd[sel,], sp = .1)
			}
		}
	}
	print(i)
}

# Print out results 
for(i in 1:length(fml.ls)){
	if(length(grep("Spl", names(fml.ls)[i])) > 0){
		next
	}
	cat("-------------------------------------------------------------------\n")
	cat("Results from model", names(fml.ls)[i], "\n")
	stargazer(fit.ls[[i]], type = "text", order = c("d_Income","year"), omit.stat = c("f","rsq", "adj.rsq","ser"))
}

# Print out results 
for(i in 1:length(fml.ls)){
	if(length(grep("Spl", names(fml.ls)[i])) > 0){
		next
	}
	cat("-------------------------------------------------------------------\n")
	cat("Results from model", names(fml.ls)[i], "\n")
	stargazer(fit.ls[[i]], type = "html", order = c("d_Income", "year"), omit.stat = c("f","rsq", "adj.rsq","ser"), 
			# keep = "d_Income", 
			column.labels = c("Exp", fmt_name), dep.var.labels.include = FALSE,
			out = paste("~/Desktop/reg_", names(fml.ls)[i], ".html", sep=""))
}

# Plot spline
i <- 4
pdf(paste("~/Desktop/graph_", names(fml.ls)[i], ".pdf", sep=""), width = 8, height = 8)
par(mfrow = c(4,2))
for(j in 1:length(dv.vec)){
	print(plot(fit.ls[[i]][[j]], xlab = "Initial income ($10,000)", ylab = expression(paste("Estiamte for ", Delta, "Income ($10,000)", sep="")), 
			main = paste("Equation for ", dv.lab[j]))
			)
}
dev.off()

i <- 5
pdf(paste("~/Desktop/graph_", names(fml.ls)[i], ".pdf", sep=""), width = 8, height = 8)
par(mfrow = c(4,2))
for(j in 1:length(dv.vec)){
	print(plot(fit.ls[[i]][[j]], xlab = "Initial income ($10,000)", 
			main = paste("Equation for ", dv.lab[j]))
			)
}
dev.off()


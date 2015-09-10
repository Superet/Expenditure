library(r2excel)
library(data.table)

setwd("~/Documents/Research/Store switching/processed data/Estimation")
source('~/Documents/Research/Store switching/Exercise/main/outreg function.R')
res 	<- read.csv("substitution_test_09sep2015.csv")
res$Pvalue	<- as.numeric(gsub("<", "",as.character(res$Probt)))

# Initiate results excel
tabfile	<- paste("/Users/chaoqunchen/Desktop/result_3_subsitution_",gsub("-","", Sys.Date()),".xlsx", sep="")
mywb 	<- createWorkbook()
sht1	<- createSheet(mywb, "Pre-recession")
sht2 	<- createSheet(mywb, "Cost")

# The distribution of outcome variables by income groups prior to the recession 
tmp 	<- subset(res, test == "Pre-recession" & substr(Parameter,1, 12) == "first_income")
tmp1	<- c("PCTL_DOL", "PCTL_Discount_Store", "PCTL_Grocery", "PCTL_Warehouse_Club", "PCTL_UPRICE", 'PCTL_SIZE')
tmp$Dependent <- factor(as.character(tmp$Dependent), levels=tmp1) 
model_list <- split(tmp, tmp$Dependent)
tmp.tab	<- model_outreg(model_list, p.given = TRUE, head.name = c("Estimate","StdErr","Pvalue","Parameter"), digits = 4)

xlsx.addHeader(mywb, sht1, value = "The position of income groups along the disbution prior to the recession.", level = 2)
xlsx.addHeader(mywb, sht1, value = "Regression of outcome variables against income groups, famsize, and time controls", level = 3)
xlsx.addTable(mywb, sht1, tmp.tab, row.names = FALSE)

# The cost of moving up and down along the distribution by segment 
# Part of the parameter names are trimmed off, so fill those missing parts first
tmp 	<- subset(res, test != "Pre-recession")
sel		<- grep("first_inc", tmp$Parameter)
tmp		<- data.table(tmp[sel,])
# tmp		<- tmp[,Parameter1 := paste(Parameter, " Qt",4:1, sep=""), by = list(test, Dependent)]
setkeyv(tmp, c("test","Dependent"))
tmp		<- data.frame(tmp)

tmp1	<- c("Cost-expenditure", "Cost-expenditurePCT","Cost-Discount_Store", "Cost-Grocery", "Cost-Warehouse_Club", "Cost-SIZE")
tmp_ls 	<- vector("list", length(tmp1))
names(tmp_ls)	<- tmp1
for(i in 1:length(tmp1)){
	tmp2		<- subset(tmp, test == tmp1[i])
	tmp2$Dependent <- as.character(tmp2$Dependent)
	model_list 	<- split(tmp2, tmp2$Dependent)
	tmp_ls[[i]]	<- model_outreg(model_list, p.given = TRUE, head.name = c("Estimate","StdErr","Pvalue","Parameter"), digits = 4)
}

for(i in 1:length(tmp1)){
	xlsx.addHeader(mywb, sht2, value = tmp1[i], level = 3)
	xlsx.addTable(mywb, sht2, tmp_ls[[i]], row.names=F)
	xlsx.addLineBreak(sht2, 1)
}

saveWorkbook(mywb, tabfile)

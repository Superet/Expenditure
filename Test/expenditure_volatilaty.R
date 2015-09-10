library(ggplot2)
library(reshape2)
library(data.table)
library(plm)
library(gridExtra)
library(scales)

options(error = quote({dump.frames(to.file = TRUE)}))

# setwd("~/Documents/Research/Store switching/processed data")
# plot.wd	<- '~/Desktop'
# source("../Exercise/main/outreg function.R")
# setwd("/home/brgordon/ccv103/Exercise/run")
setwd("/kellogg/users/marketing/2661703/Exercise/run")
# setwd("/sscc/home/c/ccv103/Exercise/run")
plot.wd 	<- getwd()
ww			<- 6.5
ww1			<- 8
ar 			<- .6

# Read data 
hh_exp		<- read.csv("2_hh_biweek_exp_merge.csv")
fmt_attr 	<- read.csv("1_format_year_attr.csv")
price_dat	<- read.csv("1_format_biweek_price.csv")
source("outreg function.R")

# Constract regression data
mydata		<- data.table(hh_exp)
mydata		<- mydata[,list(sdol = sd(dol), income_real = unique(income_real), 
						ln_income = log(income_midvalue[1]), first_income = first_income[1], famsize = famsize[1]), 
						by = list(household_code, year)]
table(mydata$first_income)
(tmp		<- quantile(mydata$first_income, c(0:3)/3))
mydata$first_incomeg	<- cut(mydata$first_income, tmp, labels = paste("T",1:3, sep=""))
mydata$year	<- factor(mydata$year)

# Test if spending volatility changes to income change
summary(lm(sdol ~ ln_income, data = mydata))

summary(myfit1 	<- plm(sdol ~ ln_income + famsize, data = mydata, index = c("household_code","year"), model="within"))
myfit2 	<- plm(sdol ~ ln_income*first_incomeg, data = mydata, index = c("household_code","year"), model="within")
myfit3	<- plm(sdol ~ year:first_incomeg, data = mydata, index = c("household_code","year"), model="within")

# Test if households stock more when income decreases.
panelist	<- data.table(hh_exp)
panelist	<- panelist[,list(first_income = first_income[1]), by = list(household_code)]
(tmp		<- quantile(panelist$first_income, c(0:3)/3))
panelist$first_incomeg	<- cut(panelist$first_income, tmp, labels = paste("T",1:3, sep=""), include.lowest= T)
panelist	<- data.frame(panelist[,list(household_code, first_incomeg)])

mydata		<- merge(hh_exp, panelist, by = "household_code", all.x=T)
mydata$ln_income 	<- log(mydata$income_midvalue)
mydata$year	<- factor(mydata$year)

summary(myfit1 	<- plm(storable_quant ~ ln_income + famsize, data = mydata, index = c("household_code","biweek"), model="within") )
summary(myfit2	<- plm(storable_quant ~ ln_income:first_incomeg + famsize, data = mydata, index = c("household_code","biweek"), model="within") )
summary(myfit3 	<- plm(storable_quant ~ year + famsize, data = mydata, index = c("household_code","biweek"), model="within") )

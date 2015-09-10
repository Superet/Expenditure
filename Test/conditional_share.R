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
hh_exp		<- read.csv("2_hh_biweek_exp_merge_larger.csv")
fmt_attr 	<- read.csv("1_format_year_attr.csv")
price_dat	<- read.csv("1_format_biweek_price.csv")
source("outreg function.R")

fmt_name 	<- as.character(sort(unique(fmt_attr$channel_type)))
R			<- length(fmt_name)
hh_exp$num_module	<- hh_exp$num_food_module + hh_exp$num_noned_module

# Compute expenditure share
con_var		<- c("num_day", "num_trip", "num_module")
ggtmp		<- subset(hh_exp, dol_purchases > 0 & !is.na(dol_purchases))
sel			<- paste("DOLP_", gsub("\\s", "_", fmt_name), sep="")
sel1		<- gsub("DOLP", "SHR", sel)
for(i in 1:length(fmt_name)){
	ggtmp[,sel1[i]] <- ggtmp[,sel[i]]/ggtmp$dol_purchases
}
ggtmp		<- ggtmp[,c("household_code", "biweek", "dol_purchases", con_var, sel1)]
ggtmp		<- subset(ggtmp, !is.na(num_day))

# Bin the conditional variables 
hist(hh_exp$dol_purchases)
summary(hh_exp$dol_purchases)
seq_y		<- seq(0, 500, 20)
ggtmp$ybin	<- cut(ggtmp$dol_purchases, seq_y, labels = seq_y[-1])
sel			<- is.na(ggtmp$ybin)
ggtmp[sel,"ybin"]	<- max(seq_y)

table(ggtmp$num_day)
quantile(ggtmp$num_day, c(1:6)/6)
seq1			<- c(0, 2, 3, 4, 6, 8, 14)
ggtmp$day_bin	<- cut(ggtmp$num_day, seq1, labels = seq1[-1])

table(ggtmp$num_trip)
quantile(ggtmp$num_trip, c(1:6)/6)
seq1			<- c(0, 2, 3, 4, 5, 7, 9, 50)
ggtmp$trip_bin	<- cut(ggtmp$num_trip, seq1, labels = seq1[-1])

summary(ggtmp$num_module)
(seq1	<- c(0, quantile(ggtmp$num_module, c(1:6)/6)))
ggtmp$module_bin 	<- cut(ggtmp$num_module, seq1, labels = seq1[-1])

# Loop over the conditional variables
selcol		<- paste("SHR_", gsub("\\s", "_", fmt_name), sep="")
convar1		<- c("day_bin", "trip_bin", "module_bin")
plots		<- list(NULL)
for(i in 1:length(con_var)){
	ggtmp1	<- data.table(melt(ggtmp[,c("ybin",convar1[i], selcol)], id.var = c("ybin", convar1[i])))
	setnames(ggtmp1, convar1[i], "d")
	ggtmp1	<- ggtmp1[,list(value = mean(value)), by = list(ybin, d, variable)]
	ggtmp1$variable <- factor(ggtmp1$variable, levels = selcol, labels = fmt_name)
	quartz()
	plots[[i]]	<- ggplot(ggtmp1, aes(ybin, value, fill = variable )) + geom_bar(stat = "identity") + 
					facet_wrap( ~ d) + 
					labs(title = paste("Expenditure share conditional on expenditure and ", con_var[i], sep=""))
	print(plots[[i]])				
}

# Run ANOVA regressions 
selcol		<- paste("SHR_", gsub("\\s", "_", fmt_name), sep="")
tmptab		<- data.frame()
for(i in 1:length(con_var)){
	for(j in 1:length(selcol)){
		cat("---------------------------------\n")
		cat("ANOVA for", selcol[j], "against", con_var[i], "\n")
		fml		<- as.formula(paste(selcol[j], "~ dol_purchases + ", con_var[i], sep="" ))
		tmpfit	<- lm(fml, data = ggtmp)
		print(tmp1	<- anova(tmpfit))
		tmptab	<- rbind(tmptab, data.frame(CondVar = con_var[i], retailer = selcol[j], Var = rownames(tmp1), tmp1))
	}
}

tmp		<- data.table(tmptab)
tmp		<- tmp[,prop := Sum.Sq/sum(Sum.Sq), by = list(CondVar, retailer)]

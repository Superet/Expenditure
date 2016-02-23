library(ggplot2)
library(reshape2)
library(data.table)
library(plm)
library(gridExtra)
library(scales)
library(r2excel)
library(stargazer)

options(error = quote({dump.frames(to.file = TRUE)}))

# setwd("~/Documents/Research/Store switching/processed data")
# plot.wd	<- '/Users/chaoqunchen/Desktop'
# setwd("/home/brgordon/ccv103/Exercise/run")
setwd("/kellogg/users/marketing/2661703/Expenditure")
# setwd("/sscc/home/c/ccv103/Exercise/run")
plot.wd 	<- getwd()
ww			<- 6.5
ww1			<- 12
ar 			<- .6

write2csv	<- TRUE
make_plot	<- TRUE
outxls		<- paste(plot.wd, "/expenditure_eda", gsub("-", "", as.character(Sys.Date())), ".xlsx", sep="")
mywb		<- createWorkbook()
sht1		<- createSheet(mywb, "Sheet1")

load("hh_biweek_exp.rdata")
codebook	<- read.csv("code_book.csv")

# Extract 5% random sample
# length(unique(hh_exp$household_code))
# sel			<- sample(unique(hh_exp$household_code), .01*length(unique(hh_exp$household_code)) )
# hh_exp_save	<- hh_exp
# hh_exp		<- subset(hh_exp, household_code %in% sel)

#################
# Organize data # 
#################
# Retail formats 
fmt_name 	<- as.character(sort(unique(fmt_attr$channel_type)))
R			<- length(fmt_name)

# Conver some date and factor variables
hh_exp$month <- month(as.Date("2004-1-1", format="%Y-%m-%d") + 14*(hh_exp$biweek-1))
hh_exp$famsize		<- factor(hh_exp$famsize, levels = c("Single","Two", "Three+"))

# Compute expenditure share
sel			<- paste("DOL_", gsub("\\s", "_", fmt_name), sep="")
sel1		<- gsub("DOL", "SHR", sel)
for(i in 1:length(fmt_name)){
	hh_exp[,sel1[i]] <- hh_exp[,sel[i]]/hh_exp$dol
}

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

# ------------------- #
# Household-year data #
pan_yr		<- data.table(hh_exp)
pan_yr		<- pan_yr[,list(Income = unique(income_midvalue), Inc = unique(income_real)), by = list(first_incomeg, household_code, year, cpi)]
setkeyv(pan_yr, c("household_code", "year"))
pan_yr		<- pan_yr[, ':='(ln_income = log(Income), cpi_adj_income = Income/cpi, recession = 1*(year >= 2008))]
pan_yr    	<- pan_yr[,':='(year_id= year - year[1] + 1, tenure = length(year), move = c(0, 1*(diff(Inc)!=0)))
                          	# move = 1*(!all(Inc==Inc[1]))) 
                    , by = list(household_code)]
pan_yr		<- pan_yr[,move_all := sum(move), by = list(household_code)]

################################################
# Variation of independent variables -- income # 
################################################
# ----------------------------------------------------------- #
# Uniqueness of the data -- within-household income variation #
# We plot the number of households with income changes over their tenure. 

# Convert income code into income description
ggtmp 	<- dcast(pan_yr, first_incomeg + tenure + move_all + household_code ~ year_id, value.var = "Inc")
names(ggtmp)[-(1:4)] <- paste("yr", names(ggtmp)[-(1:4)], sep="")
selcol 	<- grep("yr", names(ggtmp))
sel		<- codebook$code_value <= 27					# 27 is the highest value in early years
for(i in selcol){
	ggtmp[,i]	<- factor(ggtmp[,i], levels = codebook[sel,"code_value"], labels = codebook[sel,"description"])
}

# Table household tenure
tmp   <- rowSums(!is.na(ggtmp[,-(1:4)]))
tmp.tab <- table(tmp)
tmp.tab	<- t(rbind(tmp.tab, tmp.tab/sum(tmp.tab)*100))
colnames(tmp.tab)	<- c("Frequency", "Proportion")
tmp.tab	<- data.frame(Tenure = rownames(tmp.tab), tmp.tab)
cat("Table of data tenure of households:\n"); tmp.tab; cat("\n")
tmpn    <- ncol(ggtmp) - 4

# How many people never change income 
sel		<- ggtmp$tenure > 1
tmp.tab1 <- table(ggtmp[sel,"move_all"])
tmp.tab1	<- t(rbind(tmp.tab1, tmp.tab1/sum(tmp.tab1)*100))
colnames(tmp.tab1)	<- c("Frequency", "Proportion")
tmp.tab1	<- data.frame(Num.Changes = rownames(tmp.tab1), tmp.tab1)
cat("Table of number of income changes for households with tenure greater than 1:\n"); tmp.tab1; cat("\n")

# Export tenure
if(write2csv){
	xlsx.addHeader(mywb, sht1, value = "Table of data tenure of households", level = 2)
	xlsx.addTable(mywb, sht1, tmp.tab, row.names=FALSE)
	xlsx.addHeader(mywb, sht1, value = "Table of number of income changes for households with tenure greater than 1", level = 2)
	xlsx.addTable(mywb, sht1, tmp.tab1, row.names=FALSE)
}

# Plot of income transition from year 1
# We use 1st year income as x-axix, 
# We then count the number/proportion of households within a pair of income level
tmp     <- lapply(2:tmpn, function(i) cbind(melt(table(ggtmp$yr1, ggtmp[,paste("yr",i,sep="")])), year = i))
ggtmp1  <- do.call("rbind", tmp)
names(ggtmp1)	<- c("x", "y", "Freq", "year")
ggtmp1        	<- data.table(ggtmp1)
ggtmp1			<- ggtmp1[,':='(Prop= Freq/sum(Freq)*100, TotalChange = sum(Freq*1*(x!=y))/sum(Freq)*100), by = list(year)]
ggtmp1			<- ggtmp1[, year.lab:= paste(year, "-th year(", round(TotalChange, 0), "%)", sep="")]

plots			<- ggplot(ggtmp1, aes(x, y, size = Prop, alpha = .8)) + 
				  geom_point(position= "jitter", pch = 21) + 
				  geom_abline(xintercept = 0, slope = 1, size = .25, linetype = 2) + 
				  scale_size_area(max_size = 15) + 
				  facet_wrap( ~ year.lab) + 
				  # guides(alpha = FALSE, size = guide_legend(title = "Proportion(%)"), col = guide_legend(title = "n-th year")) + 
				  guides(alpha = FALSE, size = guide_legend(title = "Proportion(%)")) + 
				  labs(x = "1st Year Income", y = "nth Year Income", title = "Income transition from first year") + 
				  theme_bw() + 
				  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

if(make_plot){
	pdf(paste(plot.wd,"/graph_income_1styear_base.pdf",sep=""), width = ww1, height = ww1*ar)
	print(plots)
	dev.off()
}

# ------------ #
# Income trend #
# We plot the average income over time by group (from raw data)
ggtmp	<- pan_yr
ggtmp	<- ggtmp[, cpi_income := Income/cpi]
ggtmp	<- melt(data.frame(ggtmp[, list(first_incomeg, household_code, year, Income, cpi_income)]), 
			id.vars = c("first_incomeg", "household_code", "year"))
ggtmp$variable	<- factor(ggtmp$variable, levels = c("Income", "cpi_income"), labels = c("Income", "CPI ajdusted income"))

plots	<- list(NULL)			
# First plot the overall income trend from the raw data
plots[[1]]	<- ggplot(ggtmp, aes(year, value/1000)) + geom_smooth(method = "gam", formula = y ~ s(x, bs = "cr", k = 6)) + 
				facet_wrap(~variable) + 
				labs(x = "Year", y = "Income ($1000)", title = "Average income trend from the raw data")
# Then plot the income trend by income group from the raw data
plots[[2]]	<- ggplot(ggtmp, aes(year, value/1000, col = first_incomeg)) + 
				geom_smooth(method = "gam", formula = y ~ s(x, bs = "cr", k = 6)) + 
				facet_wrap(~variable) + 
				guides(color = guide_legend(title = "Income group")) + 
				labs(x = "Year", y = "Income ($1000)", title = "Income trend by income groups from the raw data")

# We plot the within-hosuehold income trend by group 
plots1	<- list(NULL)			
# First plot the overall income trend from the raw data
ggtmp	<- pan_yr
ggtmp	<- ggtmp[,':='(grand_avg = mean(Income), grand_cpi_avg = mean(Income/cpi))]
ggtmp	<- ggtmp[,':='(	within_income = Income - mean(Income) + grand_avg, 
						within_cpi_income = Income/cpi - mean(Income/cpi) + grand_cpi_avg), 
				by = list(household_code)]
ggtmp	<- melt(data.frame(ggtmp[, list(first_incomeg, household_code, year, within_income, within_cpi_income)]), 
			id.vars = c("first_incomeg", "household_code", "year"))		
ggtmp$variable	<- factor(ggtmp$variable, levels = c("within_income", "within_cpi_income"), labels = c("Income", "CPI ajdusted income"))
			 		
plots1[[1]]	<- ggplot(ggtmp, aes(year, value/1000)) + geom_smooth(method = "gam", formula = y ~ s(x, bs = "cr", k = 6)) + 
				facet_wrap(~variable) + 
				labs(x = "Year", y = "Income ($1000)", title = "Average within-household income trend")

# Then plot the income trend by income group from the raw data
ggtmp	<- pan_yr
ggtmp	<- ggtmp[,':='(grand_avg = mean(Income), grand_cpi_avg = mean(Income/cpi)), by = list(first_incomeg)]
ggtmp	<- ggtmp[,':='(	within_income = Income - mean(Income) + grand_avg, 
						within_cpi_income = Income/cpi - mean(Income/cpi) + grand_cpi_avg), 
				by = list(household_code)]
ggtmp	<- melt(data.frame(ggtmp[, list(first_incomeg, household_code, year, within_income, within_cpi_income)]), 
			id.vars = c("first_incomeg", "household_code", "year"))				
ggtmp$variable	<- factor(ggtmp$variable, levels = c("within_income", "within_cpi_income"), labels = c("Income", "CPI ajdusted income"))
			
plots1[[2]]	<- ggplot(ggtmp, aes(year, value/1000, col = first_incomeg)) + 
				geom_smooth(method = "gam", formula = y ~ s(x, bs = "cr", k = 6)) + 
				facet_wrap(~variable) + 
				guides(color = guide_legend(title = "Income group")) + 
				labs(x = "Year", y = "Income ($1000)", title = "Within-household income trend by income group")

if(make_plot){
	pdf(paste(plot.wd, "/graph_income_trend.pdf", sep=""), width = ww, height = ww*ar)
	print(plots[[1]])
	print(plots[[2]])
	print(plots1[[1]])
	print(plots1[[2]])	
	dev.off()
}

######################
# Channel difference # 
######################
# Distribution of retail attributes # 
# Add market-year price data 
tmp1	<- data.table(price_dat)
tmp1	<- tmp1[year > 2003, list(price = mean(bsk_price_paid_2004, na.rm=T)), by = list(scantrack_market_descr, year, channel_type)]
selcol 	<- c("size_index", "overall_prvt", "num_module", "avg_upc_per_mod")
tmpdat	<- merge(fmt_attr[,c("scantrack_market_descr", "year", "channel_type", selcol)], tmp1, 
				by = c("scantrack_market_descr", "year", "channel_type"))
ggtmp	<- melt(tmpdat, id.vars = c("scantrack_market_descr", "year", "channel_type"))				
tmplab	<- c("Size index", "Total proportion private label","Total num module","Num UPC per module", "Price index")								
ggtmp$variable <- factor(ggtmp$variable, levels= c(selcol,"price"), labels= tmplab)

plots 	<- ggplot(ggtmp, aes(channel_type, value)) + 
	stat_summary(fun.y=mean, fun.ymin=function(x) quantile(x, .25), fun.ymax = function(x) quantile(x, .75), geom ="pointrange", size=.5) + 
	stat_summary(fun.y=mean, fun.ymin=function(x) quantile(x, .025), fun.ymax = function(x) quantile(x, .975), geom ="pointrange", size=.25) + 
	facet_wrap(~variable, scales="free_y") + 
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
	labs(x = "Channel")

if(make_plot){
	pdf(paste(plot.wd, "/graph_retail_distribution.pdf", sep=""), width = ww, height = ww)
	print(plots)
	dev.off()
}

# ------------------------------------------------------------------- #
# Variance decomposition -- the importance of channel differentiation #
# ANOVA decomposition
var.vec 	<- c(selcol, "price")
anv.ls		<- matrix(NA, length(var.vec), 4, dimnames = list(tmplab, c("Channel", "scnatrack market", "Year", "Residuals")))				
for(i in 1:length(var.vec)){
	cat("--------------------------------\n")
	fml		<- as.formula(paste(var.vec[i], "~ channel_type + scantrack_market_descr + as.factor(year)", sep=""))
	tmp1	<- anova(lm(fml, data = tmpdat))
	print(tmp1)
	tmp.tab		<- tmp1[,"Sum Sq"]/sum(tmp1[,"Sum Sq"])
	anv.ls[i,]	<- tmp.tab*100
}				
cat("ANOVA variance decomposition of market-year retail attributes:\n"); print(round(anv.ls, 1)); cat("\n")

if(write2csv){
	xlsx.addHeader(mywb, sht1, value = "ANOVA variance decomposition of market-year retail attributes", level = 2)
	xlsx.addTable(mywb, sht1, anv.ls, row.names=TRUE)
}

###################################################
# Variation of dependent variables -- Expenditure # 
################################################### 
# We show expenditure and expenditure allocation vary by income level and by income group 
selcol	<- paste("DOL_", gsub("\\s", "_", fmt_name), sep="")
ggtmp	<- melt(hh_exp[, c("first_incomeg", "household_code", "income_midvalue", selcol)], 
				id.vars = c("first_incomeg", "household_code", "income_midvalue")) 				
ggtmp	<- data.table(ggtmp)
ggtmp	<- ggtmp[, list(value = mean(value)), by = list(first_incomeg, income_midvalue, variable)]
ggtmp	<- ggtmp[,':='(variable = gsub("_", "\\s", gsub("DOL_", "", variable)), 
					income_midvalue = factor(income_midvalue))]
plots	<- ggplot(ggtmp, aes(income_midvalue, value, fill = variable)) + geom_bar(stat = "identity") + 
				facet_wrap(~first_incomeg) + 
				theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
				guides(fill = guide_legend(title = "Channel")) +
				labs(x = "Income", y = "Expenditure")

if(make_plot){
	pdf(paste(plot.wd, "/graph_exp_composition.pdf", sep=""), width = ww1, height = ww1*ar)
	print(plots)
	dev.off()
}

# ------------------------------------------------------------------------------------ #
# We plot the trend of expenditure and expenditure share by income group from raw data # 
# Raw data pattern of expenditure share
selcol		<- paste("SHR_", gsub("\\s", "_", fmt_name), sep="")
ggtmp 		<- melt(hh_exp[,c("first_incomeg","biweek", "dol", selcol)], id.vars=c("first_incomeg","biweek"))
ggtmp$variable <- factor(ggtmp$variable, levels = c("dol", selcol), labels = c("Overall expenditure", fmt_name))
sp.vec		<- c(-1, .05, .1, .5)
if(make_plot){
	pdf(paste(plot.wd,"/graph_raw_share_gam.pdf",sep=""), width = ww, height = ww*.8)
	print(ggplot(subset(ggtmp, variable == "Overall expenditure"), aes(biweek, value, linetype = first_incomeg)) + 
			stat_smooth(method = "gam") + 
			labs(x = "Biweek", y = "Expenditure") + 
			guides(linetype = guide_legend(title = "Segment")) + 
			theme_bw()
		)	
	for(i in 1:length(sp.vec)){
		print(ggplot(subset(ggtmp, variable != "Overall expenditure"), aes(biweek, value, linetype = first_incomeg)) + 
				stat_smooth(method = "gam", formula = y ~ s(x, sp = sp.vec[i])) + 
				facet_wrap(~ variable, ncol = 2, scales = "free_y") + 
				scale_y_continuous(labels=percent) + 
				labs(x = "Biweek", y = "Expenditure share") + 
				guides(linetype = guide_legend(title = "Segment")) + 
				theme_bw()
			)
	}
	dev.off()
}

saveWorkbook(mywb, outxls)
cat("This program is done.\n")

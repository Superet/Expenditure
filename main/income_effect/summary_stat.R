library(ggplot2)
library(reshape2)
library(data.table)
library(gridExtra)
library(scales)
library(stargazer)

options(error = quote({dump.frames(to.file = TRUE)}))

# setwd("~/Documents/Research/Store switching/processed data")
# plot.wd	<- '/Users/chaoqunchen/Desktop'
# setwd("/home/brgordon/ccv103/Exercise/run")
setwd("/kellogg/users/marketing/2661703/Expenditure")
# setwd("/sscc/home/c/ccv103/Exercise/run")
# setwd("U:/Users/ccv103/Documents/Research/Store switching/run")

plot.wd 	<- getwd()
ww			<- 6.5
ww1			<- 10
ar 			<- .6

make_plot	<- TRUE
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
hh_exp$first_incomeg	<- as.character(hh_exp$first_incomeg)
hh_exp$first_incomeg	<- factor(hh_exp$first_incomeg, levels = paste("T", 1:3, sep=""), labels = c("Low", "Med", "High"))

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
panelist	<- panelist[,list(income = first_income[1], first_famsize = famsize[1]), by=list(first_incomeg, household_code)]
cat("Table of initial income distribution:\n"); print(table(panelist$first_incomeg)); cat("\n")
cat("Table of segments in the expenditure data:\n"); print(table(hh_exp$first_incomeg)); cat("\n")

# ------------------- #
# Household-year data #
pan_yr		<- data.table(hh_exp)
pan_yr		<- pan_yr[,list(Income = unique(income_midvalue), Inc = unique(income_real), scantrack_market_descr = unique(scantrack_market_descr)), 
						by = list(first_incomeg, household_code, year, cpi)]
setkeyv(pan_yr, c("household_code", "year"))
pan_yr		<- pan_yr[, ':='(ln_income = log(Income), cpi_adj_income = Income/cpi, recession = 1*(year >= 2008))]
pan_yr    	<- pan_yr[,':='(year_id= year - year[1] + 1, tenure = length(year), move = c(0, 1*(diff(Inc)!=0)))
                          	# move = 1*(!all(Inc==Inc[1]))) 
                    , by = list(household_code)]
pan_yr		<- pan_yr[,move_all := sum(move), by = list(household_code)]

##########################
# Table of summary stats #
##########################
# Table of tenure
tmp		<- unique(pan_yr, by = "household_code")
tmp1 	<- table(tmp$tenure)
tmp.tab	<- cbind(tmp1, prop.table(tmp1))
tmp.tab[,2]		<- round(tmp.tab[,2]*100, 2)
dimnames(tmp.tab)	<- list(Tenure = 1:nrow(tmp.tab), c("Num. HH", "Proportion (%)"))
cat("Household tenure table:\n"); print(tmp.tab); cat("\n")
stargazer(tmp.tab)

# Summary table of demogrpahics
selcol	<- c("household_size", "condo", "employed", "NumChild")
tmp		<- data.table(hh_exp[,c("household_code", selcol)])
tmp		<- unique(tmp, by = "household_code")
tmp1	<- t(apply(as.matrix(tmp)[,-1], 2, function(x) c(summary(x), SD = sd(x))))
tmptab	<- c("Household size", "Condo", "Employed", "Num. child")
rownames(tmp1)	<- tmptab
tmp1	<- tmp1[,c("Mean", "SD", "Min.", "Median", "Max.")]
cat("Summary stat of household demographics are:\n"); print(tmp1)
stargazer(tmp1)

# Summary stats of biweekly shopping behavior
tmp		<- hh_exp[,c("dol", "num_day", "num_trip")]
tmp$nchannel	<- rowSums(hh_exp[,paste("DOL_", gsub("\\s", "_", fmt_name), sep="")]>0)
tmp1	<- t(apply(tmp, 2, function(x) c(summary(x), SD = sd(x))))
tmp1	<- tmp1[,c("Mean", "SD", "Min.", "Median", "Max.")]
rownames(tmp1)	<- c("Expenditure ($)", "Num. shopping days", "Num. trips", "Num. retail formats")
cat("Summary stats of biweekly shopping pattern:\n"); print(tmp1); cat("\n")
stargazer(tmp1, summary = FALSE, digits = 2)

# Summary stats of biweekly shopping behavior by income group
tmp		<- hh_exp[,c("dol", "num_trip",paste("SHR_", gsub(" ", "_", fmt_name), sep=""))]
tmp$nchannel	<- rowSums(hh_exp[,paste("DOL_", gsub("\\s", "_", fmt_name), sep="")]>0)
tmp		<- tmp[, c("dol", "num_trip", "nchannel", paste("SHR_", gsub(" ", "_", fmt_name), sep=""))]
tmp1	<- t(apply(tmp, 2, function(x) c(mean(x, na.rm=T), tapply(x, hh_exp$first_incomeg, mean, na.rm=T))))
dimnames(tmp1)	<- list(c("Expenditure ($)", "Num. trips", "Num. retail formats", fmt_name), 
						c("Mean", "Low-income", "Med-income", "High-income"))
sel		<- rownames(tmp1) %in% fmt_name
tmp1[!sel,]	<- round(tmp1[!sel,], 1)
tmp1[sel,]	<- paste(round(tmp1[sel,]*100,1), "%", sep="")
cat("Summary stats of biweekly shopping pattern:\n"); print(tmp1); cat("\n")
stargazer(tmp1, summary = FALSE)

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


# Plot of income transition from year 1
# We use 1st year income as x-axix, 
# We then count the number/proportion of households within a pair of income level
tmp		<- unique(pan_yr, by = "household_code")
cat("The number of households who had at least one income change:", length(tmp[tenure>1 & move_all>0,household_code]), 
	"or", length(tmp[tenure>1 & move_all>0,household_code])/length(tmp[tenure>1,household_code]), ".\n" )

tmp		<- pan_yr
setkeyv(tmp, c("household_code", "year"))
# tmp		<- tmp[, income_update:=Income-c(NA, Income[-length(Income)]), by = list(household_code)]
tmp		<- tmp[, income_update:=Income-Income[1], by = list(household_code)]
mylim	<- quantile(tmp$income_update, c(.025, .975), na.rm=T)
plots	<- ggplot(subset(tmp, year_id>1), aes(income_update)) + geom_histogram(aes(y=..density..)) + 
				xlim(mylim) + 
				labs(x = "Income change relative to first year")

if(make_plot){
	pdf(paste(plot.wd,"/graph_income_1styear_base.pdf",sep=""), width = ww, height = ww*ar)
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
ggtmp$variable	<- factor(ggtmp$variable, levels = c("Income", "cpi_income"), labels = c("Income", "CPI ajdusted\n income"))

# Construct within-hh income trend 
ggtmp1	<- pan_yr
ggtmp1	<- ggtmp1[,':='(grand_avg = mean(Income[year==2004]), grand_cpi_avg = mean((Income/cpi)[year==2004]))]
ggtmp1	<- ggtmp1[,':='(	within_income = Income - Income[1] + grand_avg, 
						within_cpi_income = Income/cpi - (Income/cpi)[1] + grand_cpi_avg), 
				by = list(household_code)]
ggtmp1	<- melt(data.frame(ggtmp1[, list(first_incomeg, household_code, year, within_income, within_cpi_income)]), 
			id.vars = c("first_incomeg", "household_code", "year"))		
ggtmp1$variable	<- factor(ggtmp1$variable, levels = c("within_income", "within_cpi_income"), labels = c("Income", "CPI ajdusted\n income"))
ggtmp <- rbind(cbind(ggtmp, Variation = "Total"), cbind(ggtmp1, Variation = "Within-household"))

plots	<- ggplot(subset(ggtmp, variable != "Income"), aes(year, value/1000, linetype = Variation)) + 
				geom_smooth(method = "gam", formula = y ~ s(x, bs = "cr", k = 6), col = "black") + 
				guides(color = guide_legend(title = "Income group")) + 
				labs(x = "Year", y = "Income ($1000)")

if(make_plot){
	pdf(paste(plot.wd, "/graph_income_trend.pdf", sep=""), width = ww, height = ww*ar)
	print(plots)
	dev.off()
}

######################
# Channel difference # 
######################
# Distribution of retail attributes # 
# Number of households in the geographic market
tmp		<- data.table(hh_exp)[,list(household_code, year, scantrack_market_descr)]
tmp		<- unique(tmp, by = "household_code")
tmp1	<- table(tmp$scantrack_market_descr)

# Add market-year price data 
tmp1	<- data.table(price_dat)
tmp1	<- tmp1[year > 2003, list(price = mean(bsk_price_paid_2004, na.rm=T)), by = list(scantrack_market_descr, year, channel_type)]
selcol 	<- c("size_index", "overall_prvt", "num_module", "avg_upc_per_mod")
tmpdat	<- merge(fmt_attr[,c("scantrack_market_descr", "year", "channel_type", selcol)], tmp1, 
				by = c("scantrack_market_descr", "year", "channel_type"))
ggtmp	<- melt(tmpdat, id.vars = c("scantrack_market_descr", "year", "channel_type"))				
tmplab	<- c("Price index", "Total num categories","Num UPC per category",  "Size index", "Total proportion private label")								
ggtmp$variable <- factor(ggtmp$variable, levels= c("price", "num_module", "avg_upc_per_mod", "size_index", "overall_prvt"), labels= tmplab)
ggtmp1	<- data.table(ggtmp)
ggtmp1	<- ggtmp1[, list(middle = mean(value), lower = quantile(value, .25), upper = quantile(value, .75), 
						ymin = quantile(value, .05), ymax = quantile(value, .95)), 
					by = list(variable, channel_type)]
sel		<- order(ggtmp1[variable=="Price index",middle])
ret.ord	<- as.character(ggtmp1[variable=="Price index", channel_type][sel])
ggtmp1$channel_type	<- factor(ggtmp1$channel_type, levels = ret.ord)

# Boxplot
plots 	<- ggplot(ggtmp1, aes(x = channel_type)) + 
	geom_boxplot(aes(middle = middle, ymax = ymax, ymin = ymin, lower = lower, upper = upper), stat = "identity", size = .4) + 
	facet_wrap(~variable, scales="free_y", nrow = 1) + 
	theme(axis.text.x = element_text(angle = 75, hjust = 1)) + 
	labs(x = "Retail formats")
	
if(make_plot){
	pdf(paste(plot.wd, "/graph_retail_distribution.pdf", sep=""), width = ww1, height = ww1*ar)
	print(plots)
	dev.off()
}

# ------------------------------------------------------------------- #
# Variance decomposition -- the importance of channel differentiation #
# ANOVA decomposition
var.vec 	<- c("price", "num_module", "avg_upc_per_mod", "size_index", "overall_prvt")
tmplab
anv.ls		<- matrix(NA, length(var.vec), 4, dimnames = list(tmplab, c("Retail format(%)", "scnatrack market(%)", "Year(%)", "Residuals(%)")))				
for(i in 1:length(var.vec)){
	cat("--------------------------------\n")
	fml		<- as.formula(paste(var.vec[i], "~ channel_type + scantrack_market_descr + as.factor(year)", sep=""))
	tmp1	<- anova(lm(fml, data = tmpdat))
	print(tmp1)
	tmp.tab		<- tmp1[,"Sum Sq"]/sum(tmp1[,"Sum Sq"])
	anv.ls[i,]	<- tmp.tab*100
}				
cat("ANOVA variance decomposition of market-year retail attributes:\n"); print(round(anv.ls, 1)); cat("\n")
stargazer(anv.ls, summary= FALSE, digits = 1)

# ----------------- #
# Price time series #
# Price trend from regressions 
tmpfit 	<- lm(bsk_price_tag_2004 ~ factor(year) + channel_type + scantrack_market_descr, data = price_dat )
tmpfit1	<- lm(bsk_price_paid_2004 ~ factor(year) + channel_type + scantrack_market_descr, data = price_dat )
tmp1	<- summary(tmpfit)$coefficients
tmp2	<- summary(tmpfit1)$coefficients
sel		<- c(1, grep("year", rownames(tmp1)))
ggtmp	<- rbind(data.frame(tmp1[sel,], var = "Regular price"), data.frame(tmp2[sel,], var = "Price paid"))
names(ggtmp)	<- c("Estimate", "se", "t", "p", "var")
ggtmp$year	<- rep(2004:2012, rep = 2)
ggtmp	<- data.table(ggtmp)
ggtmp	<- ggtmp[,':='(base = Estimate[year==2004], base.se = se[year==2004]), by = list(var)]
ggtmp	<- ggtmp[,':='(Estimate = ifelse(year==2004, Estimate, Estimate + base), 
						se = ifelse(year == 2004, se, sqrt(se^2 + base.se^2)))]
plots 	<- ggplot(subset(ggtmp, var == "Price paid"), aes(year, Estimate)) + 
				geom_ribbon(aes(ymin = Estimate - 1.96*se, ymax = Estimate + 1.96*se), fill = "grey70") + 
				geom_line() 

if(make_plot){
	pdf(paste(plot.wd, "/graph_price_trend.pdf", sep=""), width = ww, height = ww*ar)
	print(plots)
	dev.off()
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
tmp		<- sort(unique(ggtmp$income_midvalue))
ggtmp	<- ggtmp[,':='(variable = gsub("_", "\\s", gsub("DOL_", "", variable)), 
					income_midvalue = factor(income_midvalue, levels = tmp))]
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

# -------------------------------------------------------------------- #
# We plot the trend of expenditure and expenditure share from raw data # 
# Raw data pattern of expenditure share
selcol		<- paste("SHR_", gsub("\\s", "_", fmt_name), sep="")
ggtmp 		<- melt(hh_exp[,c("biweek", "dol", selcol)], id.vars=c("biweek"))
ggtmp$variable <- factor(ggtmp$variable, levels = c("dol", selcol), labels = c("Overall expenditure", fmt_name))
ggtmp$biweek	<- as.Date("2003-12-30") + ggtmp$biweek*14
sp.vec		<- c(-1, .05, .1, .5)
if(make_plot){
	pdf(paste(plot.wd,"/graph_raw_share_gam.pdf",sep=""), width = ww, height = ww*.8)	
	for(i in 1:length(sp.vec)){
		print(ggplot(subset(ggtmp, variable == "Overall expenditure"), aes(biweek, value)) + 
				stat_smooth(method = "gam", formula = y ~ s(x, sp=sp.vec[i]), col = "black") + 
				labs(x = "Biweek", y = "Expenditure ($)") + 
				guides(linetype = guide_legend(title = "Income\ngroup")) + 
				theme_bw()
			)
		print(ggplot(subset(ggtmp, variable != "Overall expenditure"), aes(biweek, value)) + 
				stat_smooth(method = "gam", formula = y ~ s(x, sp = sp.vec[i]), col = "black") + 
				facet_wrap(~ variable, ncol = 2, scales = "free_y") + 
				scale_y_continuous(labels=percent) + 
				labs(x = "Biweek", y = "Expenditure share") + 
				guides(linetype = guide_legend(title = "Income\ngroup")) + 
				theme_bw()
			)
	}
	dev.off()
}

# --------------------------------------------------------------# 
# Plot the trend of expenditure share of aggregate market share # 
ggtmp	<- data.table(hh_exp[,c("first_incomeg","biweek", "dol", paste("DOL_", gsub("\\s", "_", fmt_name), sep=""))])
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
tmp		<- ggtmp[, list(value=value[1]), by =variable]
ret.ord	<- as.character(tmp[order(tmp$value, decreasing = T), variable])
ggtmp$variable	<- factor(ggtmp$variable, levels = ret.ord)

if(make_plot){
	pdf(paste(plot.wd,"/graph_raw_aggshare_gam.pdf",sep=""), width = ww, height = ww*.8)
	for(i in 1:length(sp.vec)){
		print(ggplot(subset(ggtmp, variable != "Overall expenditure"), aes(biweek, value)) + 
				stat_smooth(method = "gam", formula = y ~ s(x, sp = sp.vec[i]), col = "black") + 
				facet_wrap(~ variable, ncol = 2, scales = "free_y") + 
				scale_y_continuous(labels=percent) + 
				labs(x = "Biweek", y = "Expenditure share") + 
				guides(linetype = guide_legend(title = "Income\ngroup")) + 
				theme_bw()
			)
	}
	dev.off()
}

cat("This program is done.\n")

library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(Rcpp)
library(RcppArmadillo)
library(maxLik)
library(data.table)
library(plm)
library(zoo)

setwd("~/Documents/Research/Store switching/processed data")
plot.wd <- "~/Desktop"

setwd("//tsclient/Resear1/Store switching/processed data")
plot.wd <- "E:/Users/ccv103/Desktop/temp_results"

ww <- 6.5
ww1 <- 8
ar <- .6

##################
# Plot functions #
##################
trend_plot <- function(data, fml, groupvar, base, panel_index, subs_start, subs_end, ymonth_trans = TRUE, alpha=.05, title="", pdfname=NULL){
	group_level	<- levels(data[,groupvar])
	sel 	<- data[,groupvar] == base
	myfit1 	<- plm(fml, data=data[sel,], index=panel_index, model="within")
	myfit2 	<- plm(fml, data=data[!sel,], index=panel_index, model="within")
	ggtmp	<- rbind(data.frame(summary(myfit1)$coefficients, mkt_type =base), 
					 data.frame(summary(myfit2)$coefficients, mkt_type =setdiff(group_level, base)) )
	
	z 		<- qnorm(1 - alpha/2)
	if(ymonth_trans){
		sel			<- substr(rownames(ggtmp), 8,13) == "ymonth"
		ggtmp		<- ggtmp[sel,]
		ggtmp$time	<- substr(rownames(ggtmp), subs_start, subs_end)
		ggtmp$time	<- gsub("[[:punct:]]","", ggtmp$time)
		ggtmp$time	<- as.Date(as.yearmon(ggtmp$time))
	}else{
		ggtmp$time 	<- as.numeric(substr(rownames(ggtmp), subs_start, subs_end))
	}	
	ggtmp$up	<- with(ggtmp, Estimate + z*Std..Error)
	ggtmp$low	<- with(ggtmp, Estimate - z*Std..Error)
	plots <- ggplot(ggtmp, aes(time, Estimate, col=mkt_type, fill=mkt_type)) + geom_point() + geom_line() + 
				geom_ribbon(aes(ymin=low, ymax=up), alpha=.3) + 
				labs(title=title)
	if(is.null(pdfname)){
		print(plots)
	}else{
		pdf(file=pdfname, width=ww, height=ww*ar)
		print(plots)
		dev.off()
	}
	return(ggtmp)
}

#########################################
# Read data and construct modeling data # 
#########################################
fmt_attr 	<- read.csv("1_format_year_attr.csv")
price_dat	<- read.csv("1_format_month_price.csv")
cpi			<- read.csv("CPI.csv")

# panelist	<- read.csv("2_panelists.csv")
# hh_exp		<- read.csv("2_hh_month_exp_merge.csv")
panelist	<- read.csv("2_all_panelists.csv")
hh_exp		<- read.csv("2_all_hh_month_exp_merge.csv")

# Exclude obs with 0 DOL
fmt_name 	<- unique(fmt_attr$channel_type)
sum(hh_exp$DOL<=0) 
sum(hh_exp$net_dol<=0)
sum(is.na(hh_exp$dol_purchases))

hh_exp 		<- subset(hh_exp, DOL>0 & net_dol>0)
selc 		<- c("food_quant","nonedible_quant", "num_food_module","num_noned_module",
				"num_storable_module","num_nonstr_module", "storable_quant", "nonstr_quant")
sel1 		<- is.na(hh_exp[,selc[1]])
hh_exp[sel1,selc] <- 0
hh_exp$ymonth <- as.yearmon(paste(hh_exp$year, hh_exp$month, sep="-"))

#NOTE: drop duplicated ones
tmp <- paste(hh_exp$household_code, hh_exp$ymonth)
sel <- duplicated(tmp)
sum(sel)
hh_exp		<- hh_exp[!sel,]

# Format attributes
fmt_attr 	<- subset(fmt_attr, year > 2003)
fmt_attr$name <- paste(fmt_attr$market,fmt_attr$year, sep="-")
fmt_attr$ln_num_module <- log(fmt_attr$num_module)
fmt_attr$ln_upc_per_mod <- log(fmt_attr$avg_upc_per_mod)

# Subset panelist data such that only the households in expenditure data show up
panelist$income_year <- panelist$panel_year - 2
sel			<- unique(hh_exp$household_code)
panelist 	<- subset(panelist, household_code %in% sel)


sel 	<- sample(panelist$household_code, 1000)
panelist <- subset(panelist, household_code %in% sel)
hh_exp	<- subset(hh_exp, household_code %in% sel)

#####################
# Format attributes #
#####################
# Retail map
tmp		<- data.table(price_dat)
tmp 	<- tmp[,list(price=mean(price_paid_index)), by=list(market, year, channel_type)]
ggtmp 	<- merge(fmt_attr, data.frame(tmp), by=c("market", "year", "channel_type"))
myxlim	<- range(ggtmp$size_index)
myylim	<- range(ggtmp$price)

pdf(paste(plot.wd,"/graph_retail_map.pdf",sep=""), width=ww1, height=ww1*ar)
for(i in 2004:2010){
	print(	ggplot(subset(ggtmp, year==i), aes(size_index, price, col=channel_type, shape=market)) + 
			geom_point(size = 3) + 
			xlim(myxlim) + ylim(myylim) + 
			xlab("Size index") + ylab("Average price for a UPC") + 
			labs(title=paste("Retail map in ", i, sep="")) + 
			scale_shape(name="Market") 
	)
}
dev.off()

# Stat summary of other attributes across market and years
# selcol	<- c("size_index", "num_module", "num_brands", "avg_brands_per_mod", "num_upc", "avg_upc_per_mod", "avg_prvt_per_mod", "overall_prvt")
# tmp		<- c("Size index","Total num module","Total num brands","Num brands per module", "Num UPC","Num UPC per module",
# 		"Proportion private label per module", "Total proportion private label")
selcol 	<- c("size_index", "overall_prvt", "num_module", "avg_upc_per_mod")
tmp		<- c("Size index", "Total proportion private label","Total num module","Num UPC per module")
ggtmp <- melt(fmt_attr[,c("market", "year", "channel_type", selcol)], 
				id.var=c("market", "year", "channel_type"))
ggtmp$variable <- factor(ggtmp$variable, levels= selcol, labels= tmp)
ggplot(ggtmp, aes(channel_type, value)) + 
	stat_summary(fun.y=mean, fun.ymin=function(x) quantile(x, .25), fun.ymax = function(x) quantile(x, .75), geom ="pointrange", size=.5) + 
	stat_summary(fun.y=mean, fun.ymin=function(x) quantile(x, .025), fun.ymax = function(x) quantile(x, .975), geom ="pointrange", size=.25) + 
	facet_wrap(~variable, scales="free_y") + 
	theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Series of monthly price 
selmkt 	<- "Boston"
ggtmp 	<- subset(price_dat, market==selmkt)
ggtmp$ymonth <- as.Date(paste(ggtmp$year, ggtmp$month, "01",sep="-"), format="%Y-%m-%d")
ggplot(ggtmp, aes(ymonth, price_paid_index, col=channel_type)) + geom_point() + geom_line() + 
		labs(title = paste("Price series in ",selmkt, sep=""))

##################################
# Trend of household expenditure # 
##################################
my_alpha<- .05
myz		<- qnorm(1-my_alpha/2)
ymonth_start	<- 15			# The starting position when extract ymonth-factor in the regression
ymonth_end		<- 22

# Household income between two camps
my_fml	<- as.formula(income_midvalue ~ factor(income_year))
pdfname	<- paste(plot.wd,"/graph_trend_income.pdf",sep="")
ggtmp 	<- trend_plot(data = panelist, fml = my_fml, groupvar = "mkt_type", base = "Badcity", panel_index = c("household_code","income_year"), 
				subs_start = 20, subs_end=23, ymonth_trans=FALSE, title = "Trend of reported income", pdfname = pdfname)

# --------------------------------------------------------- # 			
# Trend of mothly expenditure without separating seasonality
my_fml	<- as.formula(DOL ~ factor(ymonth))
pdfname	<- paste(plot.wd,"/graph_trend_totalexp0.pdf",sep="")
ggtmp 	<- trend_plot(data = hh_exp, fml = my_fml, groupvar = "mkt_type", base = "Badcity", panel_index = c("household_code","ymonth"), 
				subs_start = ymonth_start, subs_end=ymonth_end, title = "Trend of monthly expenditure", pdfname = pdfname)

# Trend of mothly expenditure after controling seasonality		
sel		<- hh_exp$mkt_type == "Badcity"
tmp		<- model.matrix(~ factor(month) + factor(ymonth), data = hh_exp)
selcol	<- substr(colnames(tmp), 19, 22) != "2004" & colnames(tmp) != "(Intercept)"
tmpdata	<- data.frame(household_code = hh_exp[,"household_code"], ymonth = hh_exp[,"ymonth"], 
					DOL = hh_exp[,"DOL"], tmp[,selcol])
my_fml	<- as.formula(paste('DOL ~', paste(colnames(tmpdata)[-c(1:3)], collapse="+")))
tmpdata	<- cbind(tmpdata, mkt_type =hh_exp[,"mkt_type"])
pdfname	<- paste(plot.wd,"/graph_trend_totalexp1.pdf",sep="")
ggtmp 	<- trend_plot(data = tmpdata, fml = my_fml, groupvar = "mkt_type", base = "Badcity", panel_index = c("household_code","ymonth"), 
				subs_start = ymonth_start, subs_end=ymonth_end, pdfname = pdfname, 
				title = "Trend of monthly expenditure after controlling for seasonality")

# --------------------------------------------------------- # 			
# Trend of mothly expenditure without separating seasonality
my_fml	<- as.formula(dol_purchases ~ factor(ymonth))
pdfname	<- paste(plot.wd,"/graph_trend_netdol0.pdf",sep="")
ggtmp 	<- trend_plot(data = hh_exp, fml = my_fml, groupvar = "mkt_type", base = "Badcity", panel_index = c("household_code","ymonth"), 
				subs_start = ymonth_start, subs_end=ymonth_end, pdfname = pdfname, 
				title = "Trend of monthly expenditure")
	
# Trend of mothly expenditure after controling seasonality		
sel		<- hh_exp$mkt_type == "Badcity"
tmp		<- model.matrix(~ factor(month) + factor(ymonth), data = hh_exp)
selcol	<- substr(colnames(tmp), 19, 22) != "2004" & colnames(tmp) != "(Intercept)"
tmpdata	<- data.frame(household_code = hh_exp[,"household_code"], ymonth = hh_exp[,"ymonth"], 
					dol_purchases = hh_exp[,"dol_purchases"], tmp[,selcol])
my_fml	<- as.formula(paste('dol_purchases ~', paste(colnames(tmpdata)[-c(1:3)], collapse="+")))
tmpdata	<- cbind(tmpdata, mkt_type =hh_exp[,"mkt_type"])
pdfname	<- paste(plot.wd,"/graph_trend_netdol1.pdf",sep="")
ggtmp 	<- trend_plot(data = tmpdata, fml = my_fml, groupvar = "mkt_type", base = "Badcity", panel_index = c("household_code","ymonth"), 
				subs_start = ymonth_start, subs_end=ymonth_end, pdfname = pdfname, 
				title = "Trend of monthly expenditure after controlling for seasonality")

# --------------------------------------------------------- # 
# Trend of mothly expenditure without separating seasonality
my_fml	<- as.formula(net_dol ~ factor(ymonth))
pdfname	<- paste(plot.wd,"/graph_trend_dolpur0.pdf",sep="")
ggtmp 	<- trend_plot(data = hh_exp, fml = my_fml, groupvar = "mkt_type", base = "Badcity", panel_index = c("household_code","ymonth"), 
				subs_start = ymonth_start, subs_end=ymonth_end, pdfname=pdfname, 
				title = "Trend of monthly expenditure")

# Trend of mothly expenditure after controling seasonality		
sel		<- hh_exp$mkt_type == "Badcity"
tmp		<- model.matrix(~ factor(month) + factor(ymonth), data = hh_exp)
selcol	<- substr(colnames(tmp), 19, 22) != "2004" & colnames(tmp) != "(Intercept)"
tmpdata	<- data.frame(household_code = hh_exp[,"household_code"], ymonth = hh_exp[,"ymonth"], 
					net_dol = hh_exp[,"net_dol"], tmp[,selcol])
my_fml	<- as.formula(paste('net_dol ~', paste(colnames(tmpdata)[-c(1:3)], collapse="+")))
tmpdata	<- cbind(tmpdata, mkt_type =hh_exp[,"mkt_type"])
pdfname	<- paste(plot.wd,"/graph_trend_dolpur1.pdf",sep="")
ggtmp 	<- trend_plot(data = tmpdata, fml = my_fml, groupvar = "mkt_type", base = "Badcity", panel_index = c("household_code","ymonth"), 
				subs_start = ymonth_start, subs_end=ymonth_end, pdfname = pdfname, 
				title = "Trend of monthly expenditure after controlling for seasonality")

# ------------------------------------------------------------------------------- # 
# Trend of expenditure share at discount stores without controling for seasonality
selcol	<- paste("DOLP_", gsub("\\s", "_",fmt_name), sep="")
tmp		<- hh_exp[,selcol]/hh_exp$dol_purchases
tmpdata	<- cbind(hh_exp[,c("mkt_type","household_code","ymonth","month")], tmp)
my_fml	<- as.formula(DOLP_Discount_Store ~ factor(ymonth))
pdfname	<- paste(plot.wd,"/graph_trend_share_discount0.pdf",sep="")
ggtmp 	<- trend_plot(data = tmpdata, fml = my_fml, groupvar = "mkt_type", base = "Badcity", panel_index = c("household_code","ymonth"), 
				subs_start = ymonth_start, subs_end=ymonth_end, pdfname = pdfname, 
				title="Expenditure share at discount stores")
				
# Trend of expenditure share at discount stores after controling for seasonality
tmp		<- model.matrix(~ factor(month) + factor(ymonth), data = tmpdata)
selcol	<- substr(colnames(tmp), 19, 22) != "2004" & colnames(tmp) != "(Intercept)"
tmpdata	<- data.frame(tmpdata, tmp[,selcol])
my_fml	<- as.formula(paste('DOLP_Discount_Store ~', paste(colnames(tmpdata)[-c(1:11)], collapse="+")))
pdfname	<- paste(plot.wd,"/graph_trend_share_discount1.pdf",sep="")
ggtmp 	<- trend_plot(data = tmpdata, fml = my_fml, groupvar = "mkt_type", base = "Badcity", panel_index = c("household_code","ymonth"), 
				subs_start = ymonth_start, subs_end=ymonth_end, pdfname = pdfname,
				title="Expenditure share at discount stores after controlling for seasonality")

# ------------------------------------------------------------------------------- # 
# Trend of expenditure share at grocery without controling for seasonality
selcol	<- paste("DOLP_", gsub("\\s", "_",fmt_name), sep="")
tmp		<- hh_exp[,selcol]/hh_exp$dol_purchases
tmpdata	<- cbind(hh_exp[,c("mkt_type","household_code","ymonth","month")], tmp)
my_fml	<- as.formula(DOLP_Grocery ~ factor(ymonth))
pdfname	<- paste(plot.wd,"/graph_trend_share_grocery0.pdf",sep="")
ggtmp 	<- trend_plot(data = tmpdata, fml = my_fml, groupvar = "mkt_type", base = "Badcity", panel_index = c("household_code","ymonth"), 
				subs_start = ymonth_start, subs_end=ymonth_end, pdfname = pdfname, 
				title="Expenditure share at grocery")

# Trend of expenditure share at discount stores after controling for seasonality
tmp		<- model.matrix(~ factor(month) + factor(ymonth), data = tmpdata)
selcol	<- substr(colnames(tmp), 19, 22) != "2004" & colnames(tmp) != "(Intercept)"
tmpdata	<- data.frame(tmpdata, tmp[,selcol])
my_fml	<- as.formula(paste('DOLP_Grocery ~', paste(colnames(tmpdata)[-c(1:11)], collapse="+")))
pdfname	<- paste(plot.wd,"/graph_trend_share_grocery1.pdf",sep="")
ggtmp 	<- trend_plot(data = tmpdata, fml = my_fml, groupvar = "mkt_type", base = "Badcity", panel_index = c("household_code","ymonth"), 
				subs_start = ymonth_start, subs_end=ymonth_end, pdfname = pdfname, 
				title="Expenditure share at grocery after controlling for seasonality")

########################################
# Trend of price paid within household # 
########################################
# The proportion of coupons out of the total undiscounted prices. 
tmp		<- model.matrix( ~ factor(month) + factor(ymonth), data=hh_exp)
selcol	<- substr(colnames(tmp), 19, 22) != "2004" & colnames(tmp) != "(Intercept)"
tmpdata	<- data.frame(hh_exp[,c("mkt_type","household_code","ymonth","month","pct_coupon")], tmp[,selcol])
my_fml	<- as.formula(paste('pct_coupon ~', paste(colnames(tmpdata)[-c(1:5)], collapse="+")))
my_fml
pdfname	<- paste(plot.wd,"/graph_trend_coupon1.pdf",sep="")
ggtmp 	<- trend_plot(data = tmpdata, fml = my_fml, groupvar = "mkt_type", base = "Badcity", panel_index = c("household_code","ymonth"), 
				subs_start = ymonth_start, subs_end=ymonth_end, pdfname = pdfname, 
				title="Trend of coupon usage after controlling for seasonality")

# Unit price adjusted with CPI
tmpbase	<- cpi[cpi$Year==2004,"Annual"]
tmpdata	<- data.frame(hh_exp[,c("mkt_type","household_code","ymonth","year","month")], unit_price = with(hh_exp, dol_purchases/(storable_quant+ nonstr_quant)))
tmpdata	<- merge(tmpdata, cpi[,c("Year","Annual")], by.x="year",by.y="Year")
tmpdata$unit_price_adj	<- tmpdata$unit_price/(tmpdata$Annual/tmpbase)
tmp		<- model.matrix( ~ factor(month) + factor(ymonth), data=tmpdata)
selcol	<- substr(colnames(tmp), 19, 22) != "2004" & colnames(tmp) != "(Intercept)"
tmpdata	<- data.frame(tmpdata, tmp[,selcol])
my_fml	<- as.formula(paste('unit_price_adj ~', paste(colnames(tmpdata)[-c(1:8)], collapse="+")))
my_fml
pdfname	<- paste(plot.wd,"/graph_trend_unitprice.pdf",sep="")
ggtmp 	<- trend_plot(data = tmpdata, fml = my_fml, groupvar = "mkt_type", base = "Badcity", panel_index = c("household_code","ymonth"), 
				subs_start = ymonth_start, subs_end=ymonth_end, pdfname = pdfname, 
				title="Trend of Unit price paid by households \nadjusted with 2004 CPI")
												
###############
# Basket type # 
###############
hh_exp$Q	<- hh_exp$storable_quant + hh_exp$nonstr_quant
hh_exp$n_mod<- with(hh_exp, num_food_module + num_noned_module)

# Check distribution of basket quantiy
ggtmp <- heat_map(hh_exp$storable_quant, tmp$nonstr_quant, x.cut, y.cut, breaks=3, T)
quartz()
ggtmp1 <- heat_map(hh_exp$food_quant, tmp$nonedible_quant, x.cut, y.cut, breaks=3, T)

# Correlation betweem basket metrics
cor(hh_exp$Q, hh_exp$n_mod)
cor(hh_exp$food_quant, hh_exp$nonedible_quant)
cor(hh_exp$num_food_module, hh_exp$num_noned_module)
cor(hh_exp$storable_quant, hh_exp$nonstr_quant)

# # Classify bakset type
# my.break	<- 3
# tmp 	 	<- cut(hh_exp$storable_quant, quantile(hh_exp$storable_quant, c(0:my.break)/my.break, na.rm=T), include.lowest=T)
# tmp1 	 	<- cut(hh_exp$nonstr_quant, quantile(hh_exp$nonstr_quant, c(0:my.break)/my.break, na.rm=T), include.lowest=T)
# tmp2 <- table(tmp, tmp1)
# names(dimnames(tmp2)) <- c("Storable_quant","NonStorable_quant")
# bskt_type	<- 1:my.break^2
# names(bskt_type) <- paste(rep(levels(tmp), my.break), rep(levels(tmp1), each=my.break), sep="-")
# tmp2 		<- paste(tmp, tmp1, sep="-")
# hh_exp$basket_type	<- bskt_type[tmp2]
# table(hh_exp$basket_type)

#----------------------------------------------------------------------------- #
# Within household change of basket quantity without controlling for seasonality
my_fml	<- as.formula(Q ~ factor(ymonth))
pdfname	<- paste(plot.wd,"/graph_trend_quantity0.pdf",sep="")
ggtmp 	<- trend_plot(data = hh_exp, fml = my_fml, groupvar = "mkt_type", base = "Badcity", panel_index = c("household_code","ymonth"), 
				subs_start = ymonth_start, subs_end=ymonth_end, pdfname = pdfname, 
				title="Trend of basket quantity")

# Within household change of basket quantity controlling for seasonality
tmp		<- model.matrix( ~ factor(month) + factor(ymonth), data=hh_exp)
selcol	<- substr(colnames(tmp), 19, 22) != "2004" & colnames(tmp) != "(Intercept)"
tmpdata	<- data.frame(hh_exp[,c("mkt_type","household_code","ymonth","month","Q")], tmp[,selcol])
my_fml	<- as.formula(paste('Q ~', paste(colnames(tmpdata)[-c(1:5)], collapse="+")))
my_fml
pdfname	<- paste(plot.wd,"/graph_trend_quantity1.pdf",sep="")
ggtmp 	<- trend_plot(data = tmpdata, fml = my_fml, groupvar = "mkt_type", base = "Badcity", panel_index = c("household_code","ymonth"), 
				subs_start = ymonth_start, subs_end=ymonth_end, pdfname = pdfname, 
				title="Trend of basket quantity after controlling for seasonality")

###########################
# Trend of trip frequency # 
###########################
hist(hh_exp$num_day)
hist(hh_exp$num_trip)

#----------------------------------------------------------------------------- #
# Within household change of shopping days without controlling for seasonality
my_fml	<- as.formula(num_day ~ factor(ymonth))
pdfname	<- paste(plot.wd,"/graph_trend_days0.pdf",sep="")
ggtmp 	<- trend_plot(data = hh_exp, fml = my_fml, groupvar = "mkt_type", base = "Badcity", panel_index = c("household_code","ymonth"), 
				subs_start = ymonth_start, subs_end=ymonth_end, pdfname = pdfname, 
				title="Trend of shopping frequency (number of shopping days)")

# Within household change of shopping days controlling for seasonality
tmp		<- model.matrix( ~ factor(month) + factor(ymonth), data=hh_exp)
selcol	<- substr(colnames(tmp), 19, 22) != "2004" & colnames(tmp) != "(Intercept)"
tmpdata	<- data.frame(hh_exp[,c("mkt_type","household_code","ymonth","month","num_day")], tmp[,selcol])
my_fml	<- as.formula(paste('num_day ~', paste(colnames(tmpdata)[-c(1:5)], collapse="+")))
my_fml
pdfname	<- paste(plot.wd,"/graph_trend_days1.pdf",sep="")
ggtmp 	<- trend_plot(data = tmpdata, fml = my_fml, groupvar = "mkt_type", base = "Badcity", panel_index = c("household_code","ymonth"), 
				subs_start = ymonth_start, subs_end=ymonth_end, pdfname = pdfname, 
				title="Trend of shopping frequency (number of shopping days)\n after controlling for seasonality")

#----------------------------------------------------------------------------- #
# Within household change of number of trips without controlling for seasonality
my_fml	<- as.formula(num_trip ~ factor(ymonth))
pdfname	<- paste(plot.wd,"/graph_trend_trips0.pdf",sep="")
ggtmp 	<- trend_plot(data = hh_exp, fml = my_fml, groupvar = "mkt_type", base = "Badcity", panel_index = c("household_code","ymonth"), 
				subs_start = ymonth_start, subs_end=ymonth_end, pdfname = pdfname, 
				title="Trend of shopping frequency (number of trips)")

# Within household change of number of trips controlling for seasonality
tmp		<- model.matrix( ~ factor(month) + factor(ymonth), data=hh_exp)
selcol	<- substr(colnames(tmp), 19, 22) != "2004" & colnames(tmp) != "(Intercept)"
tmpdata	<- data.frame(hh_exp[,c("mkt_type","household_code","ymonth","month","num_trip")], tmp[,selcol])
my_fml	<- as.formula(paste('num_trip ~', paste(colnames(tmpdata)[-c(1:5)], collapse="+")))
my_fml
pdfname	<- paste(plot.wd,"/graph_trend_trips1.pdf",sep="")
ggtmp 	<- trend_plot(data = tmpdata, fml = my_fml, groupvar = "mkt_type", base = "Badcity", panel_index = c("household_code","ymonth"), 
				subs_start = ymonth_start, subs_end=ymonth_end, pdfname = pdfname, 
				title="Trend of shopping frequency (number of trips) \n after controlling for seasonality")

#######################################
# Use the regression results form SAS/STATA #
#######################################
#------------------------------------- #
# Organize data from SAS
myfit	<- read.csv("trend_fit.csv")

my_alpha<- .05
myz		<- qnorm(1-my_alpha/2)
myfit 	<- subset(myfit, !is.na(StdErr))
myfit$ymonth_date 	<- as.character(myfit$ymonth_date)
myfit$ymonth_date	<- as.Date(myfit$ymonth_date, format = "%d%b%Y")
myfit$up			<- with(myfit, Estimate + myz*StdErr)
myfit$low			<- with(myfit, Estimate - myz*StdErr)
myfit$mkt_type		<- factor(myfit$mkt_type, levels=c("","Badcity","Goodcity"), labels =c("", "Weak city","Strong city"))

tmp		<- setdiff(unique(as.character(myfit$Dependent)), "income_midvalue"); tmp
tmp1	<- c("monthly expenditure", "expenditure share at convenience",
			"expenditure share at discount", "expenditure share at dollar", 
			"expenditure share at drug", "expenditure share at grocery",
			"expenditure share at health food", "expenditure share at warehouse", 
			"basket quantity", "monthly expenditure", "expenditure/income", "number of modules","nonstorable quantity", 
			"number of shopping days","number of trips", "coupon usage", "storable quantity",
			"unit price"
			)

#------------------------------------- #
# Organize data from STATA estimation 
my_alpha<- .05
myz		<- qnorm(1-my_alpha/2)

tmp1	<- read.csv("temp_results/stata_badcity.csv")
tmp1$stat <- rep(c("Estimate","StdErr"), nrow(tmp1)/2)
sel		<- seq(2, nrow(tmp1), by=2)
tmp1[sel,"VARIABLES"] <- tmp1[(sel-1), "VARIABLES"]
tmp1	<- melt(tmp1, id.var=c("VARIABLES","stat"))
tmp1	<- dcast(tmp1, variable + VARIABLES ~ stat, value.var="value")
names(tmp1) <- c("Dependent","Parameter","Estimate","StdErr")
tmp1$mkt_type <- "Badcity"

tmp2	<- read.csv("temp_results/stata_goodcity.csv")
tmp2$stat <- rep(c("Estimate","StdErr"), nrow(tmp2)/2)
sel		<- seq(2, nrow(tmp2), by=2)
tmp2[sel,"VARIABLES"] <- tmp2[(sel-1), "VARIABLES"]
tmp2	<- melt(tmp2, id.var=c("VARIABLES","stat"))
tmp2	<- dcast(tmp2, variable + VARIABLES ~ stat, value.var="value")
names(tmp2) <- c("Dependent","Parameter","Estimate","StdErr")
tmp2$mkt_type <- "Goodcity"

myfit 	<- rbind(tmp1, tmp2)
myfit	<- subset(myfit, Parameter!="Constant")
myfit$ymonth_num <- as.numeric(gsub("[^0-9]","",myfit$Parameter))
myfit$up			<- with(myfit, Estimate + myz*StdErr)
myfit$low			<- with(myfit, Estimate - myz*StdErr)

# Transform the year-month
my_ymonth	<- data.frame(ymonth_num = 540:645, year=2004, month=1)
for(i in 2:nrow(my_ymonth)){
	tmp	<- my_ymonth[(i-1),"month"] + 1
	if(tmp>12){
		my_ymonth[i,"year"] <- my_ymonth[(i-1),"year"] + 1
		my_ymonth[i,"month"]<- 1
	}else{
		my_ymonth[i,"year"] <- my_ymonth[(i-1),"year"] 
		my_ymonth[i,"month"]<- tmp
	}
}
my_ymonth$ymonth_date <- with(my_ymonth, as.Date(paste(year, month, "01", sep="-"), format="%Y-%m-%d"))
myfit	<- merge(myfit, my_ymonth[,c("ymonth_num","ymonth_date")], by="ymonth_num")
tmp <- sort(unique(myfit$Dependent)); tmp
tmp1 <- c("monthly expenditure", "net expenditure", "purchase dollar",
			"expenditure share at convenience","expenditure share at dollar",
			"expenditure share at drug", "expenditure share at discount",  
			 "expenditure share at grocery",
			"expenditure share at health food", "expenditure share at warehouse", 
			"basket quantity", "number of modules","storable quantity","nonstorable quantity", 
			"number of shopping days","number of trips","coupon usage", "unit price"
			)

#------------------------------------- #
# Make the plot 
pdf(paste(plot.wd,"/graph_sas_results.pdf",sep=""), width=ww, height = ww*ar)
for(i in 1:length(tmp)){
	ggtmp <- subset(myfit, Dependent == tmp[i] & mkt_type!="")
	plots <- ggplot(ggtmp, aes(ymonth_date, Estimate)) + geom_point(aes(col=mkt_type), size=1) + geom_line(aes(col=mkt_type)) + 
				geom_ribbon(aes(ymin=low, ymax=up,fill=mkt_type), alpha=.3) + 
				labs(x="Month", title=paste("Trend of ",tmp1[i], sep=""), col="", fill="") 
	print(plots)
}
dev.off()

# The xpenditure share from the full sample 
pdf(paste(plot.wd,"/graph_sas_results1.pdf",sep=""), width=ww, height = ww*ar)
	i <- 3
	ggtmp <- subset(myfit, Dependent == tmp[i] & mkt_type=="")
	plots <- ggplot(ggtmp, aes(ymonth_date, Estimate)) + geom_point(size=.5) + geom_line() + 
				labs(x="Month", title=paste("Trend of ",tmp1[i], sep="")) + geom_smooth()
	print(plots)
dev.off()

# The comparison of income trend;
pdf(paste(plot.wd,"/graph_sas_income.pdf",sep=""), width=ww, height = ww*ar)
ggtmp <- subset(myfit, Dependent == "income_midvalue")
plots <- ggplot(ggtmp, aes(panel_year, Estimate)) + geom_point(aes(col=mkt_type), size=1) + geom_line(aes(col=mkt_type)) + 
			geom_ribbon(aes(ymin=low, ymax=up,fill=mkt_type), alpha=.3) + 
			labs(x = "Year", title="Trend of income", col="", fill="") 
print(plots)
dev.off()

####################
# Category example # 
####################
catdata	<- read.csv("sas_category_example.csv")
unique(catdata$module)
summary(catdata)
quantile(catdata$price, .95)
quantile(catdata$size, .95)
tmp1 	<- c(seq(0, quantile(catdata$price, .95), by=1), Inf)
tmp2	<- c(seq(0, quantile(catdata$size, .95), by= 5), Inf)
ggtmp	<- data.table(subset(catdata, year==2012))
ggtmp	<- ggtmp[,':='(price = cut(price, tmp1, include.lowest=T, labels=tmp1[-length(tmp1)]+ 1), 
						size = cut(size, tmp2, include.lowest=T, labels=tmp2[-length(tmp2)]+ 5) )]
ggtmp	<- ggtmp[,list(n = length(upcv)), by=list(channel_type, price, size)]
pdf(paste(plot.wd, "/graph_example.pdf",sep=""), width=ww, height = ww*ar)
ggplot(ggtmp, aes(size, price, size=n, color=channel_type)) + geom_point(position="jitter", alpha= .6)+ 
		scale_size_area(max_size = 21) + 
		labs(size = "Number of UPC", x = "Size (OZ)", y = "Price", title = paste("Cereal products across channel types")) + 
		theme_classic()
dev.off()




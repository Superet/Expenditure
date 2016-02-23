# Compare retail attributes from two datasets: HMS and RMS
library(ggplot2)
library(reshape2)
library(data.table)

setwd("~/Documents/Research/Store switching/processed data")

#############
# Read data # 
#############
# HMS data
hms.price	<- read.csv("1_format_biweek_price.csv", header = T)
hms.attr	<- read.csv("1_format_year_attr.csv", header = T)

# RMS data
my.wd	<- "~/Documents/Research/Store switching/rms_attributes/annual_files/channel_scantrack"
rms		<- data.frame()
for(i in 2006:2011){
	tmp	<- read.csv(paste(my.wd, "/channel_attr_", i, ".csv", sep=""), header = T)
	rms	<- rbind(rms, tmp)
}
# Remove the observations with missing market 
rms		<- subset(rms, scantrack_market_descr != "")
rms$week_end	<- as.Date(as.character(rms$week_end), format = "%Y-%m-%d")

########################
# Compare weekly price #
########################
# Plot price distribution from two dataset 
# NOTE: Price data from HMS is compuated at biweekly level, while price data from RMS is computed as weekly level. 
summary(hms.price$bsk_price_paid)
summary(rms$bskt_price)
tmp1	<- data.frame(hms.price[,c("scantrack_market_descr", "channel_type", "bsk_price_paid")], data = "HMS")
tmp2	<- data.frame(rms[,c("scantrack_market_descr", "channel_type", "bskt_price")], data = "RMS")
names(tmp1)	<- gsub("bsk_price_paid", "bskt_price", names(tmp1))
ggtmp	<- rbind(tmp1, tmp2)
p 	<- ggplot(ggtmp, aes(bskt_price, fill = data, alpha = .7)) + 
		geom_histogram(aes(y = ..density..), position = "identity") + 
		guides(alpha = FALSE) + 
		labs(x = "Bakset price index", title = "Histogram of basekt price index across scantrack market-channel-week")
print(p)
p	+ xlim(c(0, 15))
p 	+ xlim(c(0, 15)) + facet_wrap(~ channel_type)

# ANOVA decomposition 
(anv.hms	<- anova(lm(bsk_price_paid ~ scantrack_market_descr + as.factor(biweek) + channel_type, data = hms.price)))
(anv.rms	<- anova(lm(bskt_price ~ scantrack_market_descr + as.factor(week_end) + channel_type, data = rms)))
tmp.tab		<- cbind(HMS = anv.hms[,"Sum Sq"]/sum(anv.hms[,"Sum Sq"]), 
					 RMS = anv.rms[,"Sum Sq"]/sum(anv.rms[,"Sum Sq"]))
rownames(tmp.tab)	<- c("scnatrack market", "Time", "Channel", "Residuals")
cat("ANOVA decomposition of price from data sets:\n"); print(tmp.tab); cat("\n")

###################################
# Compare other retail attributes # 
###################################
# NOTE: RMS does not have data of channel convenience stores, dollar stores, mass, wholesale clubs
unique(rms$channel_type)
unique(hms.attr$channel_type)

# Aggregate from weekly data to annual level
rms.yr	<- data.table(rms)
rms.yr	<- rms.yr[, year := year(week_end)]
rms.yr	<- rms.yr[,list(size_index = mean(size_index), num_module = max(num_module), 
						num_brands = max(num_brand), avg_brands_per_mod = max(num_brand_per_mod), 
						num_upc	= max(num_upc), avg_upc_per_mod = max(num_upc_per_mod), 
						avg_prvt_per_mod = max(prvt_per_mod), overall_prvt = max(overall_prvt)), 
				by = list(scantrack_market_descr, year, channel_type)]

# Plot the distribution
mean_qt	<- function(x){
	out	<- c(mean(x, na.rm=T), quantile(x, c(.25, .75), na.rm=T))
	names(out)	<- c("y", "ymin", "ymax")
	return(out)
}

ggtmp	<- rbind(data.frame(rms.yr, data = "RMS"), data.frame(hms.attr, data = "HMS") )
ggtmp	<- melt(ggtmp, id.vars = c("scantrack_market_descr", "year", "channel_type", "data"))
ggplot(ggtmp, aes(channel_type, value, col = data)) + 
	stat_summary(fun.data = mean_qt, geom = "pointrange", position = position_dodge(width=.3)) + 
	facet_wrap(~ variable, scales = "free", ncol = 2) + 
	theme(axis.text.x = element_text(angle = 45)) + 
	labs(x = "Channel type", title = "Distribution of retail attributes across scantrack market-year")
	
# ANOVA decomposition
var.vec 	<- c("size_index", "num_module", "num_brands", "avg_brands_per_mod", "num_upc", 
				"avg_upc_per_mod", "avg_prvt_per_mod", "overall_prvt")
anv.ls		<- setNames(vector("list", length = length(var.vec)), var.vec)				
for(i in 1:length(var.vec)){
	cat("--------------------------------\n")
	cat("Variable", var.vec[i], "\n")
	fml		<- as.formula(paste(var.vec[i], "~ scantrack_market_descr + as.factor(year) + channel_type", sep=""))
	tmp1	<- anova(lm(fml, data = hms.attr))
	tmp2	<- anova(lm(fml, data = rms.yr))
	print(tmp1)
	print(tmp2)
	tmp.tab		<- cbind(HMS = tmp1[,"Sum Sq"]/sum(tmp1[,"Sum Sq"]), 
						 RMS = tmp2[,"Sum Sq"]/sum(tmp2[,"Sum Sq"]))
	rownames(tmp.tab)	<- c("scnatrack market", "Time", "Channel", "Residuals")
	tmp.tab	<- round(tmp.tab * 100, 2)
	anv.ls[[i]]	<- tmp.tab
}				

tmp.tab	<- do.call(cbind, anv.ls)
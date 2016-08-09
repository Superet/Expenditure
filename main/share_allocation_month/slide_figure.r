library(ggplot2)
library(reshape2)
library(data.table)
library(gridExtra)
library(scales)
library(stargazer)

options(error = quote({dump.frames(to.file = TRUE)}))

setwd("~/Documents/Research/Store switching/processed data")
plot.wd	<- '/Users/chaoqunchen/Desktop'
# setwd("/home/brgordon/ccv103/Exercise/run")
# setwd("/kellogg/users/marketing/2661703/Expenditure")
# setwd("/sscc/home/c/ccv103/Exercise/run")
# setwd("U:/Users/ccv103/Documents/Research/Store switching/run")

# plot.wd 	<- getwd()
ww			<- 6
ar 			<- .6

make_plot	<- TRUE
load("hh_month_exp.rdata")
codebook	<- read.csv("code_book.csv")

col.plt	<- c("#8dd3c7", "#fb8072",  "#bebada",  "#80b1d3", "#fdb462", "#b3de69")
fnt.base	<- 18
# col.plt	<- c("#7fc97f", "#beaed4", "#fdc086", "#ffff99", "#386cb0", "#f0027f")

#################
# Organize data # 
#################
# Retail formats 
fmt_name 	<- as.character(sort(unique(fmt_attr$channel_type)))
R			<- length(fmt_name)

# Conver some date and factor variables
hh_exp$first_incomeg	<- as.character(hh_exp$first_incomeg)
hh_exp$first_incomeg	<- factor(hh_exp$first_incomeg, levels = paste("T", 1:3, sep=""), labels = c("Low", "Medium", "High"))


# Compute expenditure share
sel			<- paste("DOL_", gsub("\\s", "_", fmt_name), sep="")
sel1		<- gsub("DOL", "SHR", sel)
for(i in 1:length(fmt_name)){
	hh_exp[,sel1[i]] <- hh_exp[,sel[i]]/hh_exp$dol
}

# ------------------------------------------------------ #
# Segment households based on their initial income level #
panelist	<- data.table(hh_exp)
setkeyv(panelist, c("household_code","year","month"))
panelist	<- panelist[,list(income = first_income[1], first_famsize = household_size[1]), by=list(first_incomeg, household_code)]
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

##############
# Make plots # 
##############
ret.ord	<- c("Grocery", "Discount Store", "Warehouse Club", "Drug Store", "Dollar Store","Convenience Store")
bar.col	<- "#4F81BD"

# Number of retail formats 
ggtmp1	<- rowSums(as.matrix(hh_exp[,paste("DOL_", gsub(" ", "_", fmt_name), sep="")])>0)
tmp		<- table(ggtmp1)
tmp		<- tmp/sum(tmp)
cat("Frequency of number of retail formats:\n"); print(round(tmp*100, 2)); cat("\n")

ggtmp2	<- melt(tmp)
names(ggtmp2)	<- c("n", "value")
ggtmp2	<- subset(ggtmp2, n > 0)

pdf(paste(plot.wd, "/slideg_dist_nformat.pdf", sep=""), width = 3.5, height = 3.5)
ggplot(ggtmp2, aes(n, value)) + geom_bar(stat = "identity", fill = bar.col) + 
		scale_y_continuous(labels = scales::percent) +
		scale_x_continuous(breaks = 1:R) + 
		labs(x = "Number of different retail formats\n visited in a month", y = "Fraction of observations") + 
		theme_bw() 
		# theme_classic(base_size = fnt.base)
dev.off()		

# Overall market share
tmp		<- colSums(hh_exp[,paste("DOL_", gsub(" ", "_", fmt_name), sep="")])
tmp		<- tmp/sum(tmp)
tmp1	<- colMeans(hh_exp[,paste("SHR_", gsub(" ", "_", fmt_name), sep="")], na.rm=T)
cat("Aggregation of market share:\n"); print(tmp); cat("\n")
cat("Average of expenditure share:\n"); print(tmp1); cat("\n")

ggtmp	<- data.frame(x = 1, share = tmp, retail = fmt_name)
ggtmp$retail	<- factor(ggtmp$retail, levels = rev(ret.ord))
ggtmp	<- ggtmp[order(ggtmp$retail,decreasing = T),]
ggtmp$y	<- cumsum(ggtmp$share) - ggtmp$share/2
ggtmp[ggtmp$retail == "Convenience Store","y"]	<- 1.04 #1.01
ggtmp[ggtmp$retail == "Dollar Store","y"]		<- ggtmp[ggtmp$retail == "Dollar Store","y"] + 0.01
# ggtmp$lab	<- paste(ggtmp$retail, " (", round(ggtmp$share*100), "%", ")", sep="")
ggtmp$lab	<- paste(round(ggtmp$share*100), "%", sep="")

pdf(paste(plot.wd, "/slideg_mktshare.pdf", sep=""), width = 2.5, height = 4.5)
ggplot(ggtmp,aes(x, share, fill = retail)) + geom_bar(stat = "identity", width = 1) + 
		geom_text(aes(x = x, y = y, label=lab), size = 4) + 
		geom_text(aes(x= 1.55, y = y, label = retail), size = 4, hjust = 0) + 
		scale_fill_manual(values = rev(col.plt), name = "Retail\nformat") + 
		# scale_color_manual(values = rev(col.plt)) + 
		labs(x = "", y = "Expenditure share") + 
		theme_bw() + 
		xlim(c(.5, 3)) + 
		theme(legend.position = "none", axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank(), 
			axis.title.y=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), 
			panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
			    panel.background = element_blank()) 
dev.off()

# Make an unanooted pie chart
ggtmp$share	<- c(.32, .2, .15, .08, .08, .17)
ggtmp$y	<- cumsum(ggtmp$share) - ggtmp$share/2
ggtmp$lab	<- gsub(" ", "\n", ret.ord)
ggtmp[ggtmp$lab == "Grocery","lab"]	<- "Grocery \n"
ggtmp$lab1	<- paste("s[ht", 1:R, "]", sep="")
ggtmp$retail<- as.character(ret.ord)
ggtmp$retail<- factor(ggtmp$retail, levels = ret.ord)

pdf(paste(plot.wd, "/slideg_share_eg.pdf",sep=""), width = 6, height = 1.8)
ggplot(ggtmp, aes(x, share, fill = retail)) + geom_bar(stat = "identity", width = .2) + 
	geom_text(aes(x = x, y = y, label=lab1), size = 5, parse = TRUE) + 
	geom_text(aes(x = .83, y = y, label = lab), size = 4, vjust = .2) +
	# coord_polar(theta = "y")				# Pie 
	coord_flip() + 
	scale_fill_manual(values = col.plt) + 
	scale_color_manual(values = col.plt) +
	theme_bw() +
	theme(legend.position = "none",
		axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank(), 
		axis.title.y=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), 
		panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank())
dev.off()		

# ---------------------- #
# Shares by income group #
tmp1	<- apply(hh_exp[,paste("DOL_", gsub(" ", "_", fmt_name), sep="")], 2, function(x) tapply(x, hh_exp$first_incomeg, sum))
tmp1	<- tmp1/rowSums(tmp1)
# tmp		<- colSums(as.matrix(hh_exp[,paste("DOL_", gsub(" ", "_", fmt_name), sep="")]))
# tmp		<- tmp/sum(tmp)
# tmp1		<- rbind(tmp1, Average = tmp)
ggtmp	<- melt(tmp1)
names(ggtmp)	<- c("IncGrp", "retail", "share")
tmp		<- tapply(hh_exp$dol, hh_exp$first_incomeg, mean)
ggtmp$Exp	<- tmp[as.character(ggtmp$IncGrp)]
ggtmp$dol	<- ggtmp$share*ggtmp$Exp
ggtmp$retail	<- factor(ggtmp$retail, levels = paste("DOL_", gsub(" ", "_", ret.ord), sep=""), labels = ret.ord)
ggtmp$retail	<- factor(ggtmp$retail, levels = rev(levels(ggtmp$retail)), ordered = TRUE)
ggtmp		<- ggtmp[order(ggtmp$retail, decreasing =T),]
ggtmp$lab	<- paste(round(ggtmp$share*100), "%", sep="")
# sel			<- ggtmp$IncGrp == "Low"
# ggtmp[sel, "lab"]	<- paste(ggtmp[sel,"retail"], " - ", ggtmp[sel,"lab"], sep="")
ggtmp$y		<- c(do.call(rbind, tapply(ggtmp$dol,ggtmp$IncGrp, function(x) cumsum(x) - x/2)))
ggtmp[ggtmp$retail == "Convenience Store","y"]	<- ggtmp[ggtmp$retail == "Convenience Store","Exp"] + 20	#-10
ggtmp[ggtmp$retail == "Dollar Store","y"]	<- ggtmp[ggtmp$retail == "Dollar Store","y"] + 4	#-10
ggtmp[ggtmp$retail == "Drug Store","y"]	<- ggtmp[ggtmp$retail == "Drug Store","y"] -2	#-10

ggtmp$IncGrp1	<- as.numeric(ggtmp$IncGrp)

pdf(paste(plot.wd, "/slideg_share_incgrp.pdf", sep=""), width = 6, height = 6*.8)
ggplot(ggtmp, aes(IncGrp1, dol, fill = retail)) + geom_bar(stat = "identity") + 
		geom_text(aes(x = IncGrp1, y = y, label=lab), size = 4) +
		geom_text(aes(x = IncGrp1, y = max(Exp) + 50, label = paste("$", round(Exp,0), sep="")), size = 5) +  
		# geom_text(data = subset(ggtmp, IncGrp == "High"), aes(x = 3.5, y = y, label = retail, color = retail), 
		# 			hjust = 0, size = 4, fontface = "bold") + 
		geom_text(data = subset(ggtmp, IncGrp == "High"), aes(x = 3.5, y = y, label = retail), 
					hjust = 0, size = 4) +
		scale_fill_manual(values = rev(col.plt), name = "Retail\nformat") + 
		scale_color_manual(values = rev(col.plt)) + 
		scale_x_continuous(limits = c(.5, 4.5), breaks = c(1,2,3), labels = c("Low", "Medium", "High")) + 
		labs(x = "Income group", y = "") + 
		theme(legend.position = "none",
			axis.title.y=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(),
			panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
			    panel.background = element_blank())	
dev.off()	

# ------------------------------------------------------------- #
# Plot the trend of expenditure share of aggregate market share # 
ggtmp	<- data.table(hh_exp[,c("first_incomeg","year", "month", "dol", paste("DOL_", gsub("\\s", "_", fmt_name), sep=""))])
ggtmp	<- ggtmp[,list(dol=sum(dol), 
					   DOL_Convenience_Store = sum(DOL_Convenience_Store), DOL_Discount_Store = sum(DOL_Discount_Store), 
					   DOL_Dollar_Store = sum(DOL_Dollar_Store), DOL_Drug_Store = sum(DOL_Drug_Store), 
					   DOL_Grocery = sum(DOL_Grocery), DOL_Warehouse_Club = sum(DOL_Warehouse_Club)), 
					by = list(year, month)]
for(i in paste("DOL_", gsub("\\s", "_", fmt_name), sep="")){
	ggtmp	<- ggtmp[,eval(as.name(i)):= eval(as.name(i))/dol]
}
ggtmp	<- melt(ggtmp[,!"dol", with = FALSE], id.vars=c("year", "month"))
setkeyv(ggtmp, c("variable","year","month"))
ggtmp	<- ggtmp[, share:=mean(value), by = list(variable)]
# ggtmp	<- ggtmp[, value:=value/value[1], by = list(variable)]			# Growth relative to the first month
ggtmp$variable <- factor(ggtmp$variable, levels = paste("DOL_", gsub("\\s", "_", fmt_name), sep=""), labels = fmt_name)
ggtmp$ymonth	<- as.Date(paste(ggtmp$year, "-", ggtmp$month, "-01", sep=""), format = "%Y-%m-%d")
tmp		<- ggtmp[, list(value=value[1]), by =variable]
ret.ord	<- as.character(tmp[order(tmp$value, decreasing = T), variable])
ggtmp$variable	<- factor(ggtmp$variable, levels = ret.ord)
my.sp		<- .1

pdf(paste(plot.wd,"/slideg_raw_aggshare_gam.pdf",sep=""), width = ww, height = ww*.6)
ggplot(ggtmp, aes(ymonth, value)) + 
		stat_smooth(method = "gam", formula = y ~ s(x, sp = my.sp), col = "black", se = FALSE) + 
		facet_wrap(~ variable, ncol = 3, scales = "free_y") + 
		scale_y_continuous(labels=percent) + 
		labs(x = "", y = "Market share") + 
		guides(linetype = guide_legend(title = "Income\ngroup")) + 
		theme_classic( base_size = fnt.base) + 
		theme_bw()		
dev.off()

ggplot(ggtmp,aes(ymonth, value, color = variable)) + geom_line() + 
		stat_smooth(method = "gam", formula = y ~ s(x, sp = my.sp),  se = FALSE) 

# Plot only the trend for discount stores
pdf(paste(plot.wd,"/slideg_raw_aggshare_gam_dis.pdf",sep=""), width = ww, height = ww*ar)
ggplot(subset(ggtmp, variable == "Discount Store"), aes(ymonth, value)) + 
		stat_smooth(method = "gam", formula = y ~ s(x, sp = my.sp), col = "black", se = FALSE) + 
		scale_y_continuous(labels=percent) + 
		labs(x = "", y = "Market share") + 
		guides(linetype = guide_legend(title = "Income\ngroup")) + 
		theme_bw() + 
dev.off()

# ---------------------- #
# Plot retail attributes #
# Add market-year price data 
tmp1	<- data.table(price_dat)
tmp1	<- tmp1[year > 2003, list(price = mean(bsk_price_paid_2004, na.rm=T)), by = list(scantrack_market_descr, year, channel_type)]
selcol 	<- c("size_index", "overall_prvt", "num_module", "avg_upc_per_mod")
tmpdat	<- merge(fmt_attr[,c("scantrack_market_descr", "year", "channel_type", selcol)], tmp1, 
				by = c("scantrack_market_descr", "year", "channel_type"))
ggtmp	<- melt(tmpdat, id.vars = c("scantrack_market_descr", "year", "channel_type"))				
# tmplab	<- c("Price", "No. categories","No. UPC\n per category",  "Package size", "Private label\n penetration")								
tmplab	<- c("Price", "Assortment\nwidth", "Assortment\ndepth", "Package\nsize", "Private\n label")
ggtmp$variable <- factor(ggtmp$variable, levels= c("price", "num_module", "avg_upc_per_mod", "size_index", "overall_prvt"), labels= tmplab)
ggtmp1	<- data.table(ggtmp)
ggtmp1	<- ggtmp1[, ':='(m = mean(value[channel_type == "Grocery"])), by = list(variable) ]	# Normalize the values relative to grocery
# ggtmp1	<- ggtmp1[, value:= value/m]	
ggtmp1	<- ggtmp1[, list(middle = mean(value), lower = quantile(value, .25), upper = quantile(value, .75), 
						ymin = quantile(value, .05), ymax = quantile(value, .95)), 
					by = list(variable, channel_type)]	
ggtmp1$channel_type	<- factor(ggtmp1$channel_type, levels = ret.ord, labels = gsub(" ", "\n", ret.ord))
my.dodge	<- .4

pdf(paste(plot.wd, "/slideg_attribute.pdf", sep=""), width = 5, height = 5*.8)
ggplot(subset(ggtmp1, variable == "Assortment\nwidth"), aes(x = channel_type)) + 
		geom_boxplot(aes(middle = middle, ymax = ymax, ymin = ymin, lower = lower, upper = upper), 
				stat = "identity",position = position_dodge(width = .5), notchwidth = .1)  + 
		# geom_pointrange(aes(y = middle, ymax = ymax, ymin = ymin), size = .3, position = position_dodge(width = my.dodge)) + 
		# geom_linerange(aes(ymax = upper, ymin = lower, color = channel_type), size = 1, position = position_dodge(width = my.dodge)) + 
		# scale_color_manual(values = col.plt, name = "Retial format") + 
		labs(x = "", y = "Assortment width\n(No. categories)") + 
		theme_bw() + 
		theme( panel.grid.minor = element_blank(), panel.background = element_blank(), legend.position = "bottom") #+ 
		# guides(color = guide_legend(nrow = 3, byrow = TRUE)) 							
dev.off()	





	

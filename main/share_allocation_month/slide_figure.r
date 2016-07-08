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

tmp		<- colSums(hh_exp[,paste("DOL_", gsub(" ", "_", fmt_name), sep="")])
tmp		<- tmp/sum(tmp)
tmp1	<- colMeans(hh_exp[,paste("SHR_", gsub(" ", "_", fmt_name), sep="")], na.rm=T)
cat("Aggregation of market share:\n"); print(tmp); cat("\n")
cat("Average of expenditure share:\n"); print(tmp1); cat("\n")

ggtmp	<- data.frame(x = 1, share = tmp, retail = fmt_name)
ggtmp$retail	<- factor(ggtmp$retail, levels = ret.ord)
ggtmp	<- ggtmp[order(ggtmp$retail,decreasing = T),]
ggtmp$y	<- cumsum(ggtmp$share) - ggtmp$share/2
ggtmp[ggtmp$retail == "Convenience Store","y"]	<- -.01 #1.01
ggtmp$lab	<- paste(ggtmp$retail, " - ", round(ggtmp$share*100), "%", sep="")

pdf(paste(plot.wd, "/slideg_mktshare.pdf", sep=""), width = 5*.4, height = 5)
ggplot(ggtmp,aes(x, share, fill = retail)) + geom_bar(stat = "identity") + 
		geom_text(aes(x = x, y = y, label=lab)) + 
		scale_fill_manual(values = col.plt) + 
		theme_bw() + 
		theme(legend.position = "none",
			axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank(), 
			axis.title.y=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), 
			panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
			    panel.background = element_blank()) 
dev.off()

# Make an unanooted pie chart
ggtmp$share	<- c(.35, .2, .15, .08, .1, .12)
ggtmp$y	<- cumsum(ggtmp$share) - ggtmp$share/2
ggtmp$lab	<- gsub(" ", "\n", ret.ord)
ggtmp[ggtmp$retail == "Grocery","lab"]	<- "Grocery \n"
ggtmp$lab1	<- paste("s[ht", 1:R, "]", sep="")

pdf(paste(plot.wd, "/slideg_share_eg.pdf",sep=""), width = 5, height = 1.5)
ggplot(ggtmp, aes(x, share, fill = retail)) + geom_bar(stat = "identity", width = .3) + 
	geom_text(aes(x = x, y = y, label=lab1), parse = TRUE) + 
	geom_text(aes(x = .75, y = y, label = lab, color = retail), size = 3, vjust = .2) +
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

# Shares by income group
tmp		<- apply(hh_exp[,paste("DOL_", gsub(" ", "_", fmt_name), sep="")], 2, function(x) tapply(x, hh_exp$first_incomeg, sum))
tmp		<- tmp/rowSums(tmp)
ggtmp	<- melt(tmp)
names(ggtmp)	<- c("IncGrp", "retail", "share")
tmp		<- tapply(hh_exp$dol, hh_exp$first_incomeg, mean)
ggtmp$Exp	<- tmp[as.character(ggtmp$IncGrp)]
ggtmp$dol	<- ggtmp$share*ggtmp$Exp
ggtmp$retail	<- factor(ggtmp$retail, levels = paste("DOL_", gsub(" ", "_", ret.ord), sep=""), labels = ret.ord)
ggtmp		<- ggtmp[order(ggtmp$retail, decreasing =T),]
ggtmp$lab	<- paste(round(ggtmp$share*100), "%", sep="")
sel			<- ggtmp$IncGrp == "Low"
ggtmp[sel, "lab"]	<- paste(ggtmp[sel,"retail"], " - ", ggtmp[sel,"lab"], sep="")
ggtmp$y		<- c(do.call(rbind, tapply(ggtmp$dol,ggtmp$IncGrp, function(x) cumsum(x) - x/2)))
ggtmp[ggtmp$retail == "Convenience Store","y"]	<- -5

pdf(paste(plot.wd, "/slideg_share_incgrp.pdf", sep=""), width = 6, height = 6)
ggplot(ggtmp, aes(IncGrp, dol, fill = retail)) + geom_bar(stat = "identity") + 
		geom_text(aes(x = IncGrp, y = y, label=lab)) + 
		scale_fill_manual(values = col.plt) + 
		labs(x = "Income group", y = "Monthly expenditure ($)") + 
		theme(legend.position = "none",
			# axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank(), 
			# axis.title.y=element_blank(),axis.text.y=element_blank(), axis.ticks.y=element_blank(), 
			panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
			    panel.background = element_blank())		
dev.off()	

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
ggtmp$variable <- factor(ggtmp$variable, levels = paste("DOL_", gsub("\\s", "_", fmt_name), sep=""), labels = fmt_name)
ggtmp$ymonth	<- as.Date(paste(ggtmp$year, "-", ggtmp$month, "-01", sep=""), format = "%Y-%m-%d")
tmp		<- ggtmp[, list(value=value[1]), by =variable]
ret.ord	<- as.character(tmp[order(tmp$value, decreasing = T), variable])
ggtmp$variable	<- factor(ggtmp$variable, levels = ret.ord)
my.sp		<- .05

pdf(paste(plot.wd,"/slideg_raw_aggshare_gam.pdf",sep=""), width = ww, height = ww*.8)
ggplot(subset(ggtmp, variable != "Overall expenditure"), aes(ymonth, value)) + 
		stat_smooth(method = "gam", formula = y ~ s(x, sp = my.sp), col = "black", se = FALSE) + 
		facet_wrap(~ variable, ncol = 2, scales = "free_y") + 
		scale_y_continuous(labels=percent) + 
		labs(x = "", y = "Market share") + 
		guides(linetype = guide_legend(title = "Income\ngroup")) + 
		theme_bw()
dev.off()



		
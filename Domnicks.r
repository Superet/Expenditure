library(ggplot2)
library(data.table)
library(reshape2)
library(scales)

setwd("~/Documents/Research/Store switching/Processed_data")
plot.wd	<- "~/Desktop"
trips 	<- read.csv("domnick_trips.csv", header = T)
panelists	<- read.csv("domnick_panelists.csv", header = T)
dominick	<- 61						# NOTE: retailer_code = 61 for Dominick's 

# Difference in epxneidutre share across channels between 2013 and 2014. 
fmt.name	<- c("Grocery", "Discount Store", "Warehouse Club", "Drug Store", "Dollar Store", "Convenience Store")
ggtmp	<- data.table(trips)
ggtmp	<- ggtmp[,list(dol = sum(total_spent)), by = list(household_code, domnick_ever, year, channel_type, retailer_code)]

nhh		<- length(unique(ggtmp$household_code))
nhh1	<- length(unique(ggtmp[domnick_ever==1, household_code]))
ggtmp1	<- ggtmp[,list(dol = sum(dol)/nhh), by = list(domnick_ever, year, channel_type)]
ggtmp1	<- ggtmp1[, ':='(exp = sum(dol), share = dol/sum(dol)), by = list(domnick_ever, year)]
ggtmp1$domnick_ever	<- factor(ggtmp1$domnick_ever, levels = c(1, 0), labels = c("Dominick's customers", "Non-Dominick's customers"))
ggtmp1$channel_type <- factor(ggtmp1$channel_type, levels = fmt.name)
ggtmp1$share.lab	<- paste(round(ggtmp1$share*100, 1), "%", sep="")

pdf(paste(plot.wd,"/fg_domnk_share_channel.pdf",sep=""), width = 7.5, height = 7.5*.6)
ggplot(ggtmp1, aes(channel_type, share, fill = factor(year))) + geom_bar(stat = "identity", position = position_dodge(width=0.9)) + 
		geom_text(position = position_dodge(width=0.9), aes(x = channel_type, label = share.lab, group = year), vjust = -.5) + 
		facet_grid( domnick_ever ~ . ) + 
		scale_y_continuous(labels = scales::percent, limits = c(0, max(ggtmp1$share) + .05)) + 
		guides(fill = guide_legend(title = "Year")) + 
		theme(axis.text.x = element_text(angle = 30)) + 
		theme_bw() + 
		labs(x = "Retail format", y = "Annual expenditure share", title = "Share difference before and after Dominick's closure")
dev.off()
		
# Difference in epxneidutre share within grocery stores between 2013 and 2014. 
ggtmp2	<- ggtmp[,list(dol = sum(dol)/nhh), by = list(domnick_ever, year, channel_type, retailer_code)]
ggtmp2	<- ggtmp2[, ':='(exp = sum(dol), share = dol/sum(dol)), by = list(domnick_ever, year)]
ggtmp2	<- ggtmp2[channel_type == "Grocery",]
ggtmp2$domnick_ever	<- factor(ggtmp2$domnick_ever, levels = c(1, 0), labels = c("Dominick's customers", "Non-Dominick's customers"))
ggtmp2$share.lab	<- paste(round(ggtmp2$share*100, 1), "%", sep="")
tmp		<- setNames(ggtmp2[year == 2013 & domnick_ever == "Dominick's customers", share], 
					ggtmp2[year == 2013 & domnick_ever == "Dominick's customers", retailer_code])
tmp1 	<- names(tmp)[order(tmp, decreasing = T)][1:10]
ggtmp2	<- subset(ggtmp2, retailer_code %in% as.numeric(tmp1))			# Restrict to top 10 retailers 
ggtmp2$retailer_code	<- factor(ggtmp2$retailer_code, levels = as.numeric(tmp1), labels = paste("Top ", 1:length(tmp1), sep=""))

pdf(paste(plot.wd,"/fg_domnk_share_grocery.pdf",sep=""), width = 7.5, height = 7.5*.6)
ggplot(ggtmp2, aes(retailer_code, share, fill = factor(year))) + geom_bar(stat = "identity", position = position_dodge(width=0.9)) + 
		geom_text(position = position_dodge(width=0.9), aes(x = retailer_code, label = share.lab, group = year), vjust = -.3, size = 2.5) + 
		facet_grid( domnick_ever ~ . ) + 
		scale_y_continuous(labels = scales::percent, limits = c(0, max(ggtmp2$share) + .02)) + 
		guides(fill = guide_legend(title = "Year")) + 
		theme_bw() + 
		labs(x = "Retailers", y = "Annual expenditure share", title = "Share difference with Grocery \nbefore and after Dominick's closure")
dev.off()		


# Restrict to the customers who live in downtown and north
sub.zip	<- c(60601, 60602, 60603, 60604, 60605, 60607, 60608, 60611, 60654, 60610, 60614, 60657, 60613, 60640, 60660, 60626, 
			60606, 60612, 60661, 60622, 60647, 60618, 60625)
sum(panelists$panelist_zip_code %in% sub.zip)
tmp		<- panelists[panelists$panelist_zip_code %in% sub.zip,"household_code"]
ggtmp	<- data.table(subset(trips, household_code %in% tmp))
ggtmp	<- ggtmp[,list(dol = sum(total_spent)), by = list(household_code, domnick_ever, year, channel_type, retailer_code)]

library(r2excel)
library(data.table)
library(ggplot2)
library(gridExtra)

setwd("~/Documents/Research/Store switching/processed data/Estimation")
source('~/Documents/Research/Store switching/Exercise/main/outreg function.R')
res 	<- read.csv("substitution_test_16sep2015.csv", stringsAsFactors = FALSE)
res$Pvalue	<- as.numeric(gsub("<", "",as.character(res$Probt)))

# Initiate results excel
tabfile	<- paste("/Users/chaoqunchen/Desktop/result_3_subsitution_",gsub("-","", Sys.Date()),".xlsx", sep="")
mywb 	<- createWorkbook()
sht1	<- createSheet(mywb, "Pre-recession")
sht2 	<- createSheet(mywb, "Cost")
sht3	<- createSheet(mywb, "In-recession")
make_plot	<- TRUE

# Plot setting
plot.wd	<- "~/Desktop"
ww		<- 6
ar 		<- .6
incg.col<- c("red","grey30", "grey50")			# Color of low, median, high income bar

# Function to change factor levels and label
changelab	<- function(x, ord, lab){
	y	<- rep(NA, length(x))
	for(i in 1:length(ord)){
		sel		<- grep(ord[i], x)
		y[sel]	<- lab[i]
	}
	y	<- factor(y, levels = lab)
	return(y)
}

#################################################################################
# The distribution of outcome variables by income groups prior to the recession # 
#################################################################################
tmp 	<- subset(res, test == "Pre-recession" & substr(Parameter,1, 12) == "first_income")
tmp1	<- c("PCTL_DOL", "PCTL_Discount_Store", "PCTL_Grocery", "PCTL_Warehouse_Club", "PCTL_UPRICE", 'PCTL_SIZE')
tmp$Dependent <- factor(as.character(tmp$Dependent), levels=tmp1) 
model_list <- split(tmp, tmp$Dependent)
tmp.tab	<- model_outreg(model_list, p.given = TRUE, head.name = c("Estimate","StdErr","Pvalue","Parameter"), digits = 4)

xlsx.addHeader(mywb, sht1, value = "The position of income groups along the disbution prior to the recession.", level = 2)
xlsx.addHeader(mywb, sht1, value = "Regression of outcome variables against income groups, famsize, and time controls", level = 3)
xlsx.addTable(mywb, sht1, tmp.tab, row.names = FALSE)

if(make_plot){
	tmp$IncGrp	<- changelab(tmp$Parameter, paste("T",1:3, sep=""), c("Low", "Med", "High"))
	tmp$Dependent<- factor(tmp$Dependent, levels = rev(levels(tmp$Dependent)))
	pdf(paste(plot.wd, "/graph_prior.pdf", sep=""), width = ww, height = ww*ar)
	plots 	<- ggplot(subset(tmp, IncGrp != "Low"), aes(Dependent, Estimate, fill = IncGrp)) + 
					geom_bar(stat = "identity", position = position_dodge(width=0.9)) + 
					geom_errorbar(aes(ymin = Estimate-1.96*StdErr, ymax = Estimate + 1.96 * StdErr), position=position_dodge(width=0.9), width=0.25) + 
					scale_fill_manual(values = incg.col[-1]) + 
					xlab("Dependent variable") + coord_flip()
	print(plots)
	dev.off()
}

####################################################################
# The cost of moving up and down along the distribution by segment #
####################################################################
# Part of the parameter names are trimmed off, so fill those missing parts first
tmp 	<- subset(res, substr(test,1,4) == "Cost")
sel		<- grep("first_inc", tmp$Parameter)
tmp		<- data.table(tmp[sel,])
setkeyv(tmp, c("test","Dependent"))
tmp		<- data.frame(tmp)
dv.lv 	<- c("num_mod", "quant", "unit_price", "dol_rt", "dol")
dv.lab	<- c("Num. categories", "Basket quantity", "Unit price", "Exp/Income", "Exp")
tmp$Dependent	<- factor(tmp$Dependent, levels = dv.lv, labels = dv.lab)

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

#------------------------------------------- # 
# Plot parameters
plots	<- list(NULL)
ggtmp	<- subset(tmp, Dependent != "Exp")
ggtmp$IncGrp	<- changelab(ggtmp$Parameter, paste("T",1:3, sep=""), c("Low", "Med", "High"))

# Plot the cost of moving down/up expenditure 
plots[[1]]	<- ggplot(subset(ggtmp, test == "Cost-expenditurePCT"), aes(IncGrp, Estimate, fill = IncGrp)) + 
				geom_bar(stat = "identity", position = position_dodge(width=0.9)) + 
				geom_errorbar(aes(ymin = Estimate-1.96*StdErr, ymax = Estimate + 1.96 * StdErr), position=position_dodge(width=0.9), width=0.25) + 
				scale_fill_manual(values = incg.col) + 
				xlab("") + 
				# theme(legend.position="none") + 
				facet_wrap(~ Dependent, scales = "free")
				
# Plot the cost parameters of retail substitution
ggtmp$Retail 	<- substr(ggtmp$test, 6, nchar(as.character(ggtmp$test)))
# ggtmp$Retail 	<- factor(ggtmp$Retail, levels = c("Grocery", "Discount_Store", "Warehouse_Club"))
plots[[2]]	<- ggplot(subset(ggtmp, test %in% c("Cost-Discount_Store", "Cost-Grocery", "Cost-Warehouse_Club")), 
						aes(Retail, Estimate, fill = IncGrp)) + 
				geom_bar(stat = "identity", position = position_dodge(width=0.9)) + 
				geom_errorbar(aes(ymin = Estimate-1.96*StdErr, ymax = Estimate + 1.96 * StdErr), position=position_dodge(width=0.9), width=0.25) +
				facet_wrap(~Dependent, scales = "free") + 
				theme(axis.text.x = element_text(angle = 45, hjust = .8)) + xlab("") + 
				scale_fill_manual(values = c("red","grey30", "grey50"))

# Plot the cost of 
plots[[3]]	<- ggplot(subset(ggtmp, test == "Cost-SIZE" & Dependent != "Unit price"), aes(IncGrp, Estimate, fill = IncGrp)) + 
				geom_bar(stat = "identity", position = position_dodge(width=0.9)) + 
				geom_errorbar(aes(ymin = Estimate-1.96*StdErr, ymax = Estimate + 1.96 * StdErr), position=position_dodge(width=0.9), width=0.25) + 
				scale_fill_manual(values = incg.col) + 
				xlab("") + 
				facet_wrap(~ Dependent, scales = "free")

if(make_plot){
	pdf(paste(plot.wd, "/graph_subs_cost.pdf",sep=""), width = ww, height = ww*ar)
	print(plots[[1]])
	print(plots[[2]])
	print(plots[[3]])
	dev.off()
}

######################################################
# Within-household substitution during the recession # 
######################################################
# The parameter of recession
tmp 	<- subset(res, substr(test,1,4) == "Post")
sel		<- grep("first_inc", tmp$Parameter)
tmp		<- data.table(tmp[sel,])
# tmp		<- tmp[,Parameter1 := paste(Parameter, " Qt",4:1, sep=""), by = list(test, Dependent)]
setkeyv(tmp, c("test","Dependent"))
tmp		<- data.frame(tmp)

tmp1	<- c("Post-Expenditure", "Post-Discount_Store", "Post-Dollar_Store", "Post-Grocery", "Post-Warehouse_Club", 
			 "Post-Basket_Unit_price", "Post-Purchase_Unit_price", "Post-Basket_Size", "Post-Purchase_Size")
tmp$test <- factor(as.character(tmp$test), levels=tmp1) 
model_list <- split(tmp, tmp$test)
tmp.tab	<- model_outreg(model_list, p.given = TRUE, head.name = c("Estimate","StdErr","Pvalue","Parameter"), digits = 4)

xlsx.addHeader(mywb, sht3, value = "Within-household substitutuion in the recession", level = 2)
xlsx.addHeader(mywb, sht3, value = "Regression of outcome variables against income groups*recession, household FE, month FE", level = 3)
xlsx.addTable(mywb, sht3, tmp.tab, row.names = FALSE)

#--------------------# 
# Plot the estimates # 
ggtmp			<- subset(tmp, !test %in% c("Post-Purchase_Unit_price", "Post-Purchase_Size") )
sel				<- ggtmp$Dependent %in% c("SHR_Grocery", "SHR_Discount_Store", "SHR_Dollar_Store", "SHR_Warehouse_Club")
ggtmp[sel,"Estimate"]	<- ggtmp[sel,"Estimate"] * 100					# Conver the unit to percentage for expenditure share
ggtmp$DV		<- factor(as.character(ggtmp$Dependent), levels = c("dol", "SHR_Grocery", "SHR_Discount_Store", 
						"SHR_Dollar_Store", "SHR_Warehouse_Club", "PCTL_UPRICE", "PCTL_SIZE"), 
						labels = c("Expenditure", 
								"Grocery", "Discount Store", "Dollar Store", "Warehouse Club", 
								"PCTL(Price)", "PCTL(Size)"))
ggtmp$IncGrp	<- changelab(ggtmp$Parameter, paste("T",1:3, sep=""), c("Low", "Med", "High"))


plots	<- list(NULL)
plots[[1]]	<- ggplot(subset(ggtmp, test == "Post-Expenditure"), aes(DV, Estimate, fill = IncGrp)) + 
					geom_bar(stat = "identity", position = position_dodge(width=0.9)) + 
					geom_errorbar(aes(ymin = Estimate-1.96*StdErr, ymax = Estimate + 1.96 * StdErr), position=position_dodge(width=0.9), width=0.25) +
					scale_fill_manual(values = incg.col) +
					xlab("") + ylab("Estimate($)")					
plots[[2]]	<- 	ggplot(subset(ggtmp, test %in% c("Post-Discount_Store", "Post-Dollar_Store", "Post-Grocery", "Post-Warehouse_Club")), 
						aes(DV, Estimate, fill = IncGrp)) + 
						geom_bar(stat = "identity", position = position_dodge(width=0.9)) + 
						geom_errorbar(aes(ymin = Estimate-1.96*StdErr, ymax = Estimate + 1.96 * StdErr), position=position_dodge(width=0.9), width=0.25) +
					theme(axis.text.x = element_text(angle = 45, hjust = .8)) + xlab("") + 
					xlab("") + ylab("Estiamte(%)") + 
					scale_fill_manual(values = incg.col)
plots[[3]]	<- 	ggplot(subset(ggtmp, test %in% c("Post-Basket_Unit_price", "Post-Purchase_Unit_price", "Post-Basket_Size", "Post-Purchase_Size")), 
					aes(DV, Estimate, fill = IncGrp)) + 
					geom_bar(stat = "identity", position = position_dodge(width=0.9)) + 
					geom_errorbar(aes(ymin = Estimate-1.96*StdErr, ymax = Estimate + 1.96 * StdErr), position=position_dodge(width=0.9), width=0.25) +
					theme(axis.text.x = element_text(angle = 45, hjust = .8)) + xlab("") + 
					scale_fill_manual(values = incg.col)
grid.arrange(plots[[1]], plots[[2]] + theme(legend.position = "none"), plots[[3]] + theme(legend.position = "none") , 
				ncol=1, nrow=3)
				
#--------------------------------# 				
# The parameter of income change # 				
tmp 	<- subset(res, substr(test,1,8) == "Post-Inc")
sel		<- grep("first_inc", tmp$Parameter)
tmp		<- tmp[sel,]
ord		<- order(tmp$test, tmp$Dependent)
tmp		<- tmp[ord,]

unique(tmp$test)
tmp1	<- c("Post-Inc-Expenditure", "Post-Inc-Discount_Store", "Post-Inc-Dollar_Store", "Post-Inc-Grocery", "Post-Inc-Warehouse_Club", 
			 "Post-Inc-Basket_Unit_price", "Post-Inc-Basket_Size")
tmp$test <- factor(as.character(tmp$test), levels=tmp1) 
model_list <- split(tmp, tmp$test)
tmp.tab	<- model_outreg(model_list, p.given = TRUE, head.name = c("Estimate","StdErr","Pvalue","Parameter"), digits = 4)

xlsx.addHeader(mywb, sht3, value = "Within-household substitutuion in the recession", level = 2)
xlsx.addHeader(mywb, sht3, value = "Regression of outcome variables against income groups*recession, household FE, month FE", level = 3)
xlsx.addTable(mywb, sht3, tmp.tab, row.names = FALSE)

#--------------------# 
# Plot the estimates # 
ggtmp			<- tmp
sel				<- ggtmp$Dependent %in% c("SHR_Grocery", "SHR_Discount_Store", "SHR_Dollar_Store", "SHR_Warehouse_Club")
ggtmp[sel,"Estimate"]	<- ggtmp[sel,"Estimate"] * 100					# Conver the unit to percentage for expenditure share
ggtmp[sel,"StdErr"]		<- ggtmp[sel,"StdErr"] * 100
ggtmp$DV		<- factor(as.character(ggtmp$Dependent), levels = c("dol", "SHR_Grocery", "SHR_Discount_Store", 
						"SHR_Dollar_Store", "SHR_Warehouse_Club", "PCTL_UPRICE", "PCTL_SIZE"), 
						labels = c("Expenditure", 
								"Grocery", "Discount Store", "Dollar Store", "Warehouse Club", 
								"PCTL(Price)", "PCTL(Size)"))
ggtmp$IncGrp	<- changelab(ggtmp$Parameter, paste("T",1:3, sep=""), c("Low", "Med", "High"))


plots	<- list(NULL)
plots[[1]]	<- ggplot(subset(ggtmp, test == "Post-Inc-Expenditure"), aes(DV, Estimate, fill = IncGrp)) + 
					geom_bar(stat = "identity", position = position_dodge(width=0.9)) + 
					geom_errorbar(aes(ymin = Estimate-1.96*StdErr, ymax = Estimate + 1.96 * StdErr), position=position_dodge(width=0.9), width=0.25) +
					scale_fill_manual(values = incg.col) +
					xlab("") + ylab("Estimate($)")					
plots[[2]]	<- 	ggplot(subset(ggtmp, test %in% c("Post-Inc-Discount_Store", "Post-Inc-Dollar_Store", "Post-Inc-Grocery", "Post-Inc-Warehouse_Club")), 
						aes(DV, Estimate, fill = IncGrp)) + 
						geom_bar(stat = "identity", position = position_dodge(width=0.9)) + 
						geom_errorbar(aes(ymin = Estimate-1.96*StdErr, ymax = Estimate + 1.96 * StdErr), position=position_dodge(width=0.9), width=0.25) +
					theme(axis.text.x = element_text(angle = 45, hjust = .8)) + xlab("") + 
					xlab("") + ylab("Estiamte(%)") + 
					scale_fill_manual(values = incg.col)
plots[[3]]	<- 	ggplot(subset(ggtmp, test %in% c("Post-Inc-Basket_Unit_price", "Post-Inc-Purchase_Unit_price", "Post-Inc-Basket_Size", "Post-Inc-Purchase_Size")), 
					aes(DV, Estimate, fill = IncGrp)) + 
					geom_bar(stat = "identity", position = position_dodge(width=0.9)) + 
					geom_errorbar(aes(ymin = Estimate-1.96*StdErr, ymax = Estimate + 1.96 * StdErr), position=position_dodge(width=0.9), width=0.25) +
					theme(axis.text.x = element_text(angle = 45, hjust = .8)) + xlab("") + 
					scale_fill_manual(values = incg.col)
if(make_plot){
	pdf(paste(plot.wd,"/graph_subs_income.pdf",sep=""), width = ww, height = ww*ar)
	# print(grid.arrange(plots[[1]], plots[[2]] + theme(legend.position = "none"), plots[[3]] + theme(legend.position = "none") , 
	# 				ncol=1, nrow=3))
	print(plots[[1]])
	print(plots[[2]])
	print(plots[[3]])
	dev.off()
}					



saveWorkbook(mywb, tabfile)

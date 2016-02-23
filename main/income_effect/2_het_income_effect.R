library(ggplot2)
library(reshape2)
library(data.table)
library(plm)
library(gridExtra)
library(scales)
library(r2excel)
library(systemfit)

options(error = quote({dump.frames(to.file = TRUE)}))

# setwd("~/Documents/Research/Store switching/processed data")
# plot.wd	<- '/home/chaoqunchen/Desktop'
# source("../Exercise/main/outreg function.R")
# setwd("/home/brgordon/ccv103/Exercise/run")
setwd("/kellogg/users/marketing/2661703/Exercise/run")
# setwd("/sscc/home/c/ccv103/Exercise/run")
plot.wd 	<- getwd()
ww			<- 6.5
ww1			<- 12
ar 			<- .6

write2csv	<- TRUE
make_plot	<- TRUE
outfile		<- paste(plot.wd, "/2_cpiincome_effect_", gsub("-", "", as.character(Sys.Date())), ".rdata", sep="")
outxls		<- paste(plot.wd, "/2_het_cpiincome_effect_", gsub("-", "", as.character(Sys.Date())), ".xlsx", sep="")
mywb		<- createWorkbook()
sht1		<- createSheet(mywb, "Sheet1")
sht2		<- createSheet(mywb, "Regression")
sht3		<- createSheet(mywb, "SUR")

# Read data 
# hh_exp		<- read.csv("2_hh_biweek_exp_merge_larger.csv")
# fmt_attr 	<- read.csv("1_format_year_attr.csv")
# price_dat	<- read.csv("1_format_biweek_price.csv")
load("hh_biweek_exp.rdata")
codebook	<- read.csv("code_book.csv")
source("outreg function.R")

# Extract 5% random sample
# length(unique(hh_exp$household_code))
# sel			<- sample(unique(hh_exp$household_code), .01*length(unique(hh_exp$household_code)) )
# hh_exp_save	<- hh_exp
# hh_exp		<- subset(hh_exp, household_code %in% sel)

#############
# Functions # 
#############
my_forward 	<- function(x){
	return(c(x[-1], NA))
}

my_lag		<- function(x){
	return(c(NA, x[-length(x)]))
}

mytransition <- function(x1, x2){
	if(class(x1)!="factor" | class(x2) != "factor"){
		x1 	<- factor(x1)
		x2	<- factor(x2)
	}
	onem	<- diag(length(levels(x1)))
	out		<- table(x1, x2)
	sel		<- which(apply(out, 1, function(x) all(x==0)))
	if(length(sel)>0) { out[sel,] <- onem[sel,] }
	return(out/rowSums(out))
}

#######################################
# Organize data and create covariates # 
#######################################
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

# Segment households based on their initial income level 
panelist	<- data.table(hh_exp)
setkeyv(panelist, c("household_code","year","biweek"))
panelist	<- panelist[,list(income = first_income[1], first_famsize = famsize[1]), by=list(household_code)]
tmp			<- quantile(panelist$income, c(0, .33, .67, 1))
num_grp		<- 3
panelist	<- panelist[, first_incomeg := cut(panelist$income, tmp, labels = paste("T", 1:num_grp, sep=""), include.lowest = T)]
hh_exp		<- merge(hh_exp, data.frame(panelist)[,c("household_code", "first_incomeg")], by = "household_code", all.x=T )
hh_exp$first_incomeg	<- factor(hh_exp$first_incomeg, levels = c("T2", "T1", "T3"))
# hh_exp$first_incomeg	<- cut(hh_exp$first_income, c(0, 17, 23, 30), labels = c("T1", "T2", "T3"))
cat("Table of initial income distribution:\n"); print(table(panelist$first_incomeg)); cat("\n")
cat("Table of segments in the expenditure data:\n"); print(table(hh_exp$first_incomeg)); cat("\n")

###############################
# Characterize income changes # 
###############################
# The within-household trend of cpi-adjusted income 
pan_yr		<- data.table(hh_exp)
pan_yr		<- pan_yr[,list(Income = unique(income_midvalue), Inc = unique(income_real)), by = list(famsize, first_incomeg, first_income, household_code, year)]
# pan_yr		<- merge(pan_yr, cpi, by="year", all.x=T)
# pan_yr$Income	 <- pan_yr$Income / pan_yr$cpi
pan_yr$ln_income <- log(pan_yr$Income)
pan_yr$recession <- 1 * (pan_yr$year >= 2008)
setkeyv(pan_yr, c("household_code", "year"))
pan_yr    <- pan_yr[,':='(year_id= year - year[1] + 1, tenure = length(year), move = c(0, 1*(diff(Inc)!=0)))
                          	# move = 1*(!all(Inc==Inc[1]))) 
                    , by = list(household_code)]
pan_yr		<- pan_yr[,move_all := sum(move), by = list(household_code)]
pan_yr		<- data.frame(pan_yr)

# Within-household regression: income ~ year/recession
pan_yr$year <- factor(pan_yr$year)
summary(lm(Income ~ year*first_incomeg, data = pan_yr))
myfit 		<- plm(Income ~ first_incomeg*year, data=pan_yr, index = c("household_code","year"), model="within")
tmp			<- summary(myfit)
cat("Income trend regression:\n"); print(tmp); cat("\n")
myfit1		<- plm(Income ~ first_incomeg*recession, data=pan_yr, index = c("household_code","year"), model="within")
tmp1		<- summary(myfit1)
cat("Income trend regression:\n"); print(tmp1); cat("\n")

# Export regression coefficients 
if(write2csv){
	model_list 	<- list(year=data.frame(tmp$coefficients, Var = rownames(tmp$coefficients)), 
						recession = data.frame(tmp1$coefficients, Var = rownames(tmp1$coefficients)))
	tmp.tab		<- model_outreg(model_list, p.given = TRUE, head.name = c("Estimate","Std..Error", "Pr...t.."), digits = 4)
	xlsx.addHeader(mywb, sht1, value = "Fix-effect regressions of household income", level = 2)
	xlsx.addTable(mywb, sht1, tmp.tab, row.names=FALSE)
}

# Plot the income trend 
tmp			<- unique(pan_yr[,c("first_incomeg","year")])
tmp$year	<- factor(tmp$year)
tmpx		<- model.matrix(~first_incomeg*year, data = tmp)[,-(1:3)]
identical(colnames(tmpx), names(coef(myfit)))
ggtmp		<- tmpx %*% coef(myfit)
ggtmp		<- cbind(tmp, ggtmp)
sel 		<- pan_yr$year == 2004
tmp			<- tapply(pan_yr[sel,"Income"], pan_yr[sel,"first_incomeg"], mean)
ggtmp$Income<- (ggtmp$ggtmp + tmp[ggtmp$first_incomeg])/1000
ggtmp$year	<- as.numeric(as.character(ggtmp$year))

plots		<- list(NULL)
plots[[1]]	<- ggplot(pan_yr, aes(as.numeric(year), y = Income/1000, linetype = first_incomeg)) + geom_smooth(method=lm, formula = y ~ x + I(x^2) + I(x^3)) + 
					guides(linetype = guide_legend(title = "Segment")) + 
					labs(x = "Year", y = "Income($1000)", title = "Smoothing line of raw data")
plots[[2]]	<- ggplot(ggtmp, aes(year, Income, linetype = first_incomeg)) + geom_line() + 
			    		labs(x = "Year", y = "Income($1000)", title = "Smoothing line from within-model prediction") + 
						guides(linetype = guide_legend(title = "Segment"))
if(make_plot){
	pdf(paste(plot.wd,"/graph_income.pdf",sep=""), width = ww, height = ww*ar)
	print(plots[[1]])
	print(plots[[2]])
	dev.off()	
}

#------------------------------------------------------#
# Income transition within households; 
ggtmp <- dcast(pan_yr, first_incomeg + tenure + move_all + household_code ~ year_id, value.var = "Inc")
names(ggtmp)[-(1:4)] <- paste("yr", names(ggtmp)[-(1:4)], sep="")
selcol 	<- grep("yr", names(ggtmp))
sel		<- codebook$code_value <= 27					# 27 is the highest value in early years
for(i in selcol){
	ggtmp[,i]	<- factor(ggtmp[,i], levels = codebook[sel,"code_value"], labels = codebook[sel,"description"])
}

# Check household tenure
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
tmp     <- lapply(2:tmpn, function(i) cbind(melt(table(ggtmp$yr1, ggtmp[,paste("yr",i,sep="")])), year = i))
ggtmp1  <- do.call("rbind", tmp)
names(ggtmp1) <- c("x", "y", "Freq", "year")
ggtmp1        <- data.table(ggtmp1)
ggtmp1        <- ggtmp1[, Prop:= Freq/sum(Freq)*100, by = list(year)]
plots		<- list(NULL)
plots[[1]]	<- ggplot(ggtmp1, aes(x, y, size = Prop, col= factor(year), alpha = .8)) + 
				  geom_point(position= "jitter", pch = 21) + 
				  geom_abline(xintercept = 0, slope = 1, size = .25, linetype = 2) + 
				  scale_size_area(max_size = 15) + 
				  facet_wrap( ~ year) + 
				  guides(alpha = FALSE, size = guide_legend(title = "Proportion(%)"), col = guide_legend(title = "n")) + 
				  labs(x = "1st Year Income", y = "nth Year Income", title = "Income transition from 1st year") + 
				  theme_bw() + 
				  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
plots[[2]] 	<- ggplot(ggtmp1, aes(x, y, size = Prop, col= factor(year), alpha = .8)) + 
		        geom_point(position= "jitter", pch = 21) + 
		        geom_abline(xintercept = 0, slope = 1, size = .25, linetype = 2) + 
		        scale_size_area(max_size = 15) + 
		        guides(alpha = FALSE, size = guide_legend(title = "Proportion(%)"), col = guide_legend(title = "n")) + 
		        labs(x = "1st Year Income", y = "nth Year Income", title = "Income transition from 1st year") + 
			  	theme_bw() + 
			  	theme(axis.text.x = element_text(angle = 45, hjust = 1))
if(make_plot){
	pdf(paste(plot.wd,"/graph_income_1styear_base.pdf",sep=""), width = ww1, height = ww1*ar)
	print(plots[[1]])
	print(plots[[2]])
	dev.off()
}

# Year to year transition
tmp     <- lapply(2:tmpn, function(i) cbind(melt(table(ggtmp[,paste("yr",i-1,sep="")], ggtmp[,paste("yr",i,sep="")])), year = i-1))
ggtmp2  <- do.call("rbind", tmp)
names(ggtmp2) <- c("x", "y", "Freq", "year")
ggtmp2  <- data.table(ggtmp2)
ggtmp2  <- ggtmp2[, Prop:=Freq/sum(Freq)*100, by = list(year)]
plots 	<- list(NULL)
plots[[1]]	<- ggplot(ggtmp2, aes(x, y, size = Prop, col= factor(year), alpha = .8)) + 
					geom_point(position= "jitter", pch = 21) + 
					scale_size_area(max_size = 15) + 
					facet_wrap(~ year) + 
					guides(alpha = FALSE, size = guide_legend(title = "Proportion(%)"), col = guide_legend(title = "n")) + 
					labs(x = "Income in focal year n", y = "Income in the following year", title = "Income transition from year to year")
plots[[2]]	<- ggplot(ggtmp2, aes(x, y, size = Prop, col= factor(year), alpha = .8)) + 
					geom_point(position= "jitter", pch = 21) + 
					scale_size_area(max_size = 15) + 
					guides(alpha = FALSE, size = guide_legend(title = "Proportion(%)"), col = guide_legend(title = "n")) + 
					labs(x = "Income in focal year n", y = "Income in the following year", title = "Income transition from year to year")

if(make_plot){
	pdf(paste(plot.wd,"/graph_income_cont_transition.pdf",sep=""), width = ww1, height = ww1*ar)
	print(plots[[1]])
	print(plots[[2]])
	dev.off()
}

# Plot of income transition from year 1, by famsize
ggtmp 	<- dcast(pan_yr, famsize + first_incomeg + tenure + move + household_code ~ year_id, value.var = "Inc")
tmp1	<- 5
names(ggtmp)[-(1:tmp1)] <- paste("yr", names(ggtmp)[-(1:tmp1)], sep="")
tmp   	<- rowSums(!is.na(ggtmp[,-(1:tmp1)]))
cat("Table of data tenure of households:\n"); table(tmp); cat("\n")
cat("Table of data tenure by family size\n"); table(ggtmp$famsize, ggtmp$tenure); cat("\n")
tmpn    <- ncol(ggtmp) - tmp1
tmp1	<- unique(ggtmp$famsize)

plots	<- list(NULL)
for(i in 1:length(tmp1)){
	sel		<- ggtmp$famsize == tmp1[i]
	tmp     <- lapply(2:tmpn, function(j) cbind(melt(table(ggtmp[sel,"yr1"], ggtmp[sel,paste("yr",j,sep="")])), year = j ))
	ggtmp1  <- do.call("rbind", tmp)
	names(ggtmp1) <- c("x", "y", "Freq", "year")
	ggtmp1        <- data.table(ggtmp1)
	ggtmp1        <- ggtmp1[, Prop:= Freq/sum(Freq)*100, by = list(year)]
	plots[[i]]	<- ggplot(ggtmp1, aes(x, y, size = Prop, col= factor(year), alpha = .8)) + 
					  geom_point(position= "jitter", pch = 21) + 
					  geom_abline(xintercept = 0, slope = 1, size = .25, linetype = 2) + 
					  scale_size_area(max_size = 12, breaks = c(0, 2, 4, 8, 12)) + 
					  facet_wrap( ~ year) + 
					  guides(alpha = FALSE, size = guide_legend(title = "Proportion(%)"), col = guide_legend(title = "n")) + 
					  labs(x = "1st Year Income", y = "nth Year Income", title = "Income transition from 1st year")
}

if(make_plot){
	pdf(paste(plot.wd,"/graph_income_1styear_famsize.pdf",sep=""), width = ww1, height = ww1*ar)
	for(i in 1:length(plots)){
		print(plots[[i]])
	}
	dev.off()
}

#------------------------------------------------------#
# Expectation and variance
# Data of current income and next-year income
tmp		<- data.table(pan_yr)
sel		<- tmp$Inc > 27
tmp[sel,"Inc"]	<- 27
tmp		<- tmp[, income_forward := my_forward(Inc), by=list(first_incomeg, household_code)]
tmp1	<- sort(unique(tmp$Inc))
tmp$Inc <- factor(tmp$Inc, levels = tmp1)
tmp$income_forward <- factor(tmp$income_forward, levels = tmp1)
tmp$seg	<- paste(tmp$first_incomeg, tmp$year, sep="-")

# Compute the transition matrix of income for each segment
tmp		<- data.frame(tmp)
sel		<- tmp$year != 2010
tmp1	<- split(tmp[sel,], tmp[sel,"seg"])
tmp2	<- lapply(tmp1, function(x) mytransition(x$Inc, x$income_forward))
sel		<- codebook$code_value <= 27
tmp.mat	<- rep(1, length(unique(tmp$Inc))) %*% t(as.vector(codebook[sel,"mid_range"]))

# Compute expected future income and uncertainty
tmp3	<- matrix(NA, length(tmp2), 2, dimnames = list(names(tmp1), c("Expectation", "SD")))
for(i in 1:length(tmp2)){
	# sel			<- ifelse(i%%6==0, 1, i%%6)
	# tmp.mat1	<- tmp.mat/cpi[sel,"cpi"]
	tmp.mat1	<- tmp.mat
	ee			<- rowSums(tmp2[[i]]*tmp.mat1)
	vv			<- rowSums(tmp2[[i]]*tmp.mat1^2) - ee^2
	tmp.wt		<- table(tmp1[[i]]$Inc)/nrow(tmp1[[i]])
	tmp3[i,1]	<- sum(tmp.wt * ee)
	tmp3[i,2]	<- sqrt(sum(tmp.wt^2 * vv))
}
tmps	<- do.call(rbind, strsplit(rownames(tmp3),"-"))
tmp		<- split(data.frame(tmp3), tmps[,1])
tmp.tab	<- do.call(cbind, tmp)

if(write2csv){
	xlsx.addHeader(mywb, sht1, value = "Expectation and sd of future income:", level = 2)
	xlsx.addTable(mywb, sht1, tmp.tab, row.names = F)
}

######################################
# Plot share trend from the raw data #
######################################
# Aggregate pattern of market share
ggtmp		<- data.table(hh_exp)
ggtmp		<- ggtmp[,list(DOL_Convenience_Store = sum(DOL_Convenience_Store), DOL_Discount_Store = sum(DOL_Discount_Store), 
							DOL_Dollar_Store = sum(DOL_Dollar_Store), DOL_Drug_Store = sum(DOL_Drug_Store), 
							DOL_Grocery = sum(DOL_Grocery), DOL_Warehouse_Club = sum(DOL_Warehouse_Club), y = sum(dol)), 
					by = list(first_incomeg, biweek)]
ggtmp		<- data.frame(ggtmp)
ggtmp[,3:(2+R)]	<- ggtmp[,3:(2+R)]/ggtmp$y
ggtmp		<- melt(ggtmp[,-ncol(ggtmp)], id.vars=c("first_incomeg", "biweek"))
ggtmp		<- subset(ggtmp, !is.na(variable))
tmp			<- tapply(ggtmp$value, ggtmp$variable, mean)
tmp			<- names(tmp)[order(tmp,decreasing=T)]
tmp1		<- gsub("_", " ", gsub("DOL_","",tmp))
ggtmp$variable	<- factor(as.character(ggtmp$variable), levels=tmp, labels=tmp1)
Rc_t		<- min(hh_exp[hh_exp$recession == 1, "biweek"])

if(make_plot){
	pdf(paste(plot.wd, "/graph_raw_share.pdf",sep=""), width=ww1, height = ww1*ar)
	print(
		ggplot(ggtmp, aes(biweek, value, fill = variable, order = as.numeric(variable))) + geom_bar(stat="identity") + 
				facet_wrap( ~ first_incomeg) + 
				labs(x = "Biweek", y = "Expenditure share") + 
				scale_y_continuous(labels=percent) + 
				scale_fill_grey() + 
				guides(fill = guide_legend(reverse = TRUE, title="Retail format")) + 
				theme_bw()
	)
	dev.off()
}

# Raw data pattern of expenditure share
tmp_dv		<- paste("SHR_", gsub("\\s", "_", fmt_name), sep="")
ggtmp 		<- melt(hh_exp[,c("first_incomeg","biweek",tmp_dv)], id.vars=c("first_incomeg","biweek"))
ggtmp$variable <- factor(ggtmp$variable, levels = tmp_dv, labels = fmt_name)
if(make_plot){
	pdf(paste(plot.wd,"/graph_raw_share_gam.pdf",sep=""), width = ww, height = ww*.8)
	print(ggplot(ggtmp, aes(biweek, value, linetype = first_incomeg)) + stat_smooth() + 
			facet_wrap(~ variable, ncol = 2, scales = "free_y") + 
			scale_y_continuous(labels=percent) + 
			labs(x = "Biweek", y = "Expenditure share") + 
			guides(linetype = guide_legend(title = "Segment")) + 
			theme_bw()
		)
	dev.off()
}

#########################################
# Regressions to explore income effects # 
#########################################
#---------------------------------------#
# Construct regression data 
mydata 	<- hh_exp
# Not adjust income by cpi 
# mydata		<- merge(mydata, cpi, by="year", all.x=T)

# Add price 
tmp		<- dcast(price_dat, scantrack_market_descr+biweek ~ channel_type, value.var = "bsk_price_paid")
colnames(tmp)	<- c("scantrack_market_descr", "biweek", paste("PRC_", gsub(" ", "_", fmt_name), sep=""))
mydata 	<- merge(mydata, tmp, by = c("scantrack_market_descr", "biweek"), all.x = T)

# Shopping incidence
tmp_dv		<- paste("SHR_", gsub("\\s", "_", fmt_name), sep="")
for(i in tmp_dv){
	mydata[,i] <- mydata[,i]*100
}
sel			<- paste("DOL_", gsub("\\s", "_", fmt_name), sep="")
sel1		<- gsub("DOL", "IC", sel)
for(i in 1:length(fmt_name)){
	mydata[,sel1[i]] <- 1*(mydata[,sel[i]]>0)
}

mydata$ln_income	<- log(mydata$income_midvalue/mydata$cpi)
sum(mydata$dol == 0 )/nrow(mydata)
mydata				<- subset(mydata, dol > 0)
# mydata$ln_income	<- log(mydata$income_midvalue)
mydata$month		<- factor(mydata$month)
mydata$ln_dol 		<- log(mydata$dol)
mydata$recession 	<- factor(mydata$recession)
mydata$year			<- as.factor(mydata$year)

# Add lagged purchases
mydata	<- data.table(mydata)
mydata	<- mydata[,lag_dol := my_lag(dol), by = list(household_code)]
mydata	<- data.frame(mydata)
mydata	<- subset(mydata, !is.na(lag_dol))
mydata$ln_incomeT1 <- 1*(mydata$first_incomeg == "T1")*mydata$ln_income
mydata$ln_incomeT3 <- 1*(mydata$first_incomeg == "T3")*mydata$ln_income

#---------------------------------------#
# Set up DV and IV, regression models 
dv_vec	<- c("ln_dol", paste("SHR_", gsub("\\s", "_", fmt_name), sep=""), paste("IC_", gsub("\\s", "_", fmt_name), sep="") )
myfml	<- list(Homogeneity		= as.formula("y ~ ln_income + month"), 
				Heterogeneity	= as.formula("y ~ ln_income*first_incomeg + month"), 
				Year			= as.formula("y ~ ln_income*first_incomeg + year + month"), 
				Price			= as.formula("y ~ ln_income*first_incomeg + 
											PRC_Convenience_Store + PRC_Discount_Store + PRC_Dollar_Store + 
											PRC_Drug_Store + PRC_Grocery + PRC_Warehouse_Club + month")
				)

# Run regressions
regrs	<- data.frame()
for(i in 1:length(dv_vec)){
	prc 		<- proc.time()
	for(j in 1:length(myfml)){
		rm(list = "myfit")
		tmp 	<- as.formula(substitute(y ~ x, list(y = as.name(dv_vec[i]), x = terms(myfml[[j]])[[3]])) )
		if(dv_vec[i] != "ln_dol"){
			tmp		<- update(tmp, . ~ . + lag_dol)
		}
		myfit	<- plm(tmp, data = mydata, index = c("household_code","biweek"), model="within")
		tmp2	<- summary(myfit)
		tmp3	<- data.frame(model = names(myfml)[j], DV = dv_vec[i], tmp2$coefficients)
		tmp3$Var<- rownames(tmp3)
		rownames(tmp3) <- NULL
		regrs	<- rbind(regrs, tmp3)
	}
	use.time 	<- proc.time() - prc
	cat("Regressions for", dv_vec[i], "finishes, using", use.time[3]/60, "min.\n")
}

# Print regression results 
tmp			<- subset(regrs, substr(regrs$Var, 1, 5) != "month" & substr(regrs$Var, 1, 4) != "year")
tmpls		<- vector("list", length(myfml))
names(tmpls)<- names(myfml)
for(i in 1:length(myfml)){
	tmp1	<- subset(tmp, model == names(myfml)[i])
	model_list	<- split(tmp1, tmp1$DV)
	tmp.tab	<- model_outreg(model_list, p.given = TRUE, head.name = c("Estimate","Std..Error","Pr...t..","Var"), digits = 4)
	tmpls[[i]]	<- tmp.tab
}

if(write2csv){
	for(i in 1:length(tmpls)){
		xlsx.addHeader(mywb, sht2, value = names(tmpls)[i], level = 3)
		xlsx.addTable(mywb, sht2, tmpls[[i]], row.names = F)
	}
}

#####################################
# SUR regression with fixed effects #
#####################################
my.SURFE <- function(syseq, data, panid, ...){
# This function fits system equations SUR with fixed effects. 
# It requires that independent variables are the same. 
	mfr		<- model.frame(syseq[[1]], data = data, weights = eval(as.name(panid)))
	X		<- model.matrix(syseq[[1]], data = data)
	X 		<- X[,-which(colnames(X)=="(Intercept)")]
	y		<- sapply(syseq, function(x) model.extract(model.frame(x, data), "response") )
	mydata	<- data.table(cbind(y, X))
	varvec	<- colnames(mydata)
	mydata$id	<- model.extract(mfr, "weights")
	
	# Substract group mean
	for(i in 1:length(varvec)){
		var		<- as.name(varvec[i])
		mydata 	<- mydata[,eval(var):=eval(var) - mean(eval(var), na.rm=T), by = list(id)]
	}
	mydata	<- data.frame(mydata)
	
	# Correct the system of formula
	syseq_new	<- syseq
	Xname		<- setdiff(colnames(mydata), c("id", colnames(y)))
	for(i in 1:length(syseq)){
		syseq_new[[i]]	<- as.formula( paste(colnames(y)[i], "~", paste(Xname, collapse = "+"), "-1", sep="") )
	}
	
	# Fit SUR
	myfit 	<- systemfit(syseq_new, data = mydata, method = "SUR", control = systemfit.control( ... ))
	print(summary(myfit))
	coeftab <- data.frame(summary(myfit)$coefficients)
	coeftab$df			<- myfit$df.residual/length(syseq)
	coeftab$adj_df		<- coeftab$df - length(unique(data[,panid]))
	coeftab$adj_StdErr	<- with(coeftab, Std..Error*sqrt(df/adj_df))
	coeftab$adj_t		<- with(coeftab, Estimate/adj_StdErr)
	coeftab$adj_p		<- with(coeftab, pt(-abs(adj_t), adj_df)*2)
	
	return(coeftab)
}

# Fit SUR with fixed effects 
dv_mat	<- cbind(paste("SHR_", gsub("\\s", "_", fmt_name), sep=""), paste("IC_", gsub("\\s", "_", fmt_name), sep=""))
dv_mat	<- dv_mat[-1,]
myfml	<- list(Homogeneity		= as.formula("y ~ ln_income + month + lag_dol"), 
				Heterogeneity	= as.formula("y ~ ln_income + ln_incomeT1 + ln_incomeT3 + month + lag_dol"), 
				Year 			= as.formula("y ~ ln_income + ln_incomeT1 + ln_incomeT3 + year + month + lag_dol"), 
				Price 			= as.formula("y ~ ln_income + ln_incomeT1 + ln_incomeT3 + lag_dol + month + 
													PRC_Convenience_Store + PRC_Discount_Store + PRC_Dollar_Store + 
													PRC_Drug_Store + PRC_Grocery + PRC_Warehouse_Club")
				)
				
regrs1 	<- data.frame()
panid	<- "household_code"
fmln	<- c("Homogeneity", "T1", "T2", "T3")
for(i in 1:ncol(dv_mat)){
	# for(j in 1:4){
	# 	prc 		<- proc.time()
	# 	syseq	<- lapply(dv_mat[,i], function(yn) 
	# 				as.formula(substitute(y ~ x, list(y = as.name(yn), x = terms(myfml[[1]])[[3]])) ) )
	# 	names(syseq)	<- gsub("_", ".", dv_mat[,i])
	# 	if(j == 1){
	# 		tmp		<- try(my.SURFE(syseq, mydata, panid = "household_code", maxiter = 2, residCovWeighted = TRUE, model = FALSE))
	# 	}else{
	# 		tmp		<- try(my.SURFE(syseq, subset(mydata, first_incomeg == fmln[j]), panid = "household_code"))
	# 	}
	# 
	# 	if(class(tmp) != "try-error"){
	# 		tmp1 	<- sapply(rownames(tmp), function(x) strsplit(x, "_")[[1]][1])
	# 		tmp2	<- sapply(1:length(tmp1), function(i) gsub(paste(tmp1[i], "_", sep=""), "", rownames(tmp)[i]))
	# 		tmp3	<- cbind(model = fmln[j], DV = tmp1, Var = tmp2, tmp)
	# 		rownames(tmp3) <- NULL
	# 		regrs1	<- rbind(regrs1, tmp3)
	# 		use.time 	<- proc.time() - prc
	# 		cat("Regressions for", ifelse(i==1, "share", "incidence"), "with", ifelse(j==1, "homogenous", "heterogenous"),
	# 			"responses finishes, using", use.time[3]/60, "min.\n")
	# 	}else{
	# 		cat("Regressions for", ifelse(i==1, "share", "incidence"), "with", ifelse(j==1, "homogenous", "heterogenous"),
	# 			"responses fails.\n")
	# 	}
	# }
	for(j in 1:length(myfml)){
		rm(list = "tmp")
		prc 		<- proc.time()
		syseq	<- lapply(dv_mat[,i], function(yn) 
					as.formula(substitute(y ~ x, list(y = as.name(yn), x = terms(myfml[[j]])[[3]])) ) )
		names(syseq)	<- gsub("_", ".", dv_mat[,i])
		tmp		<- try(my.SURFE(syseq, mydata, panid = "household_code", maxiter = 2, residCovWeighted = TRUE, model = FALSE))
		if(class(tmp) != "try-error"){
			tmp1 	<- sapply(rownames(tmp), function(x) strsplit(x, "_")[[1]][1])
			tmp2	<- sapply(1:length(tmp1), function(i) gsub(paste(tmp1[i], "_", sep=""), "", rownames(tmp)[i]))
			tmp3	<- cbind(model = names(myfml)[j], DV = tmp1, Var = tmp2, tmp)
			rownames(tmp3) <- NULL
			regrs1	<- rbind(regrs1, tmp3)
			use.time 	<- proc.time() - prc
			cat("Regressions for", ifelse(i==1, "share", "incidence"), "with", names(myfml)[j],
				"responses finishes, using", use.time[3]/60, "min.\n")
		}else{
			cat("Regressions for", ifelse(i==1, "share", "incidence"), "with", names(myfml)[j],
				"responses fails.\n")
		}
	}
}

# Print regression results 
tmp			<- subset(regrs1, substr(regrs1$Var, 1, 5) != "month" & substr(regrs1$Var, 1, 4) != "year")
tmpls		<- vector("list", length(myfml))
names(tmpls)<- names(myfml)
for(i in 1:length(myfml)){
	tmp1	<- subset(tmp, model == names(myfml)[i])
	if(nrow(tmp1) !=0 ){
		model_list	<- split(tmp1, tmp1$DV)
		tmp.tab	<- model_outreg(model_list, p.given = TRUE, head.name = c("Estimate","adj_StdErr","adj_p","Var"), digits = 4)
		tmpls[[i]]	<- tmp.tab
	}
}
sel 	<- sapply(tmpls, is.null)
tmpls	<- tmpls[!sel]

if(write2csv){
	for(i in 1:length(tmpls)){
		xlsx.addHeader(mywb, sht3, value = names(tmpls)[i], level = 3)
		xlsx.addTable(mywb, sht3, tmpls[[i]], row.names = F)
	}
}

saveWorkbook(mywb, outxls)

save.image(file = outfile)
cat("This program is done.")


#################
# Post plotting #
#################
plot.wd	<- "~/Desktop"
ww		<- 4
ar		<- .8
incg.col<- c("red","grey30", "grey50")			# Color of low, median, high income bar

# Function to change factor levels and label
changelab	<- function(x, ord, lab){
	y	<- rep(NA, length(x))
	for(i in 1:length(ord)){
		sel		<- x == ord[i]
		y[sel]	<- lab[i]
	}
	y	<- factor(y, levels = lab)
	return(y)
}

# Homogenous income effect for expenditure in SUR regressions
ggtmp1	<- subset(regrs1, model == "Homogeneity" & substr(DV,1,3) == "SHR" & Var == "ln_income")
ggtmp1$Retail	<- factor(as.character(ggtmp1$DV), levels = paste("SHR.", gsub("\\s", ".", fmt_name[-1]), sep=""), 
						labels = fmt_name[-1])
ord		<- order(ggtmp1$Estimate)
ggtmp1$Retail 	<- factor(ggtmp1$Retail, levels = ggtmp1[ord,"Retail"])
ord		<- levels(ggtmp1$Retail)
						
if(make_plot){
	pdf(paste(plot.wd, "/graph_income_shr_avg.pdf", sep=""), width = ww*ar, height = ww)
	plots	<- ggplot(ggtmp1, aes(Retail, Estimate)) + 
				geom_bar(stat = "identity", width = .9) + 
				geom_errorbar(aes(ymin = Estimate-1.96*Std..Error, ymax = Estimate + 1.96 * Std..Error), width=0.25) + 
				# geom_pointrange(aes(ymin = Estimate-1.96*Std..Error, ymax = Estimate + 1.96 * Std..Error)) +
				xlab("Dependent variable") + coord_flip()
	print(plots)
	dev.off()
}

# Heterogenous income effect for expenditure in SUR regressions
ggtmp1	<- subset(regrs1, model == "Heterogeneity" & substr(DV,1,3) == "SHR" & substr(Var,1,9) == "ln_income")
ggtmp1$Retail	<- factor(as.character(ggtmp1$DV), levels = paste("SHR.", gsub("\\s", ".", fmt_name[-1]), sep=""), 
						labels = fmt_name[-1])						
ggtmp1$IncGrp	<- changelab(ggtmp1$Var, paste("ln_income",c("T1","","T3"),sep=""), c("Low", "Med", "High"))
ggtmp1$IncGrp	<- factor(ggtmp1$IncGrp, levels = rev(levels(ggtmp1$IncGrp)), ordered = TRUE)
ggtmp1$Retail 	<- factor(ggtmp1$Retail, levels = ord)

# Adjust the difference to estiamte
sel		<- ggtmp1$Var == "ln_income"
tmp		<- ggtmp1[sel,"Estimate"]
names(tmp)	<- ggtmp1[sel,"DV"]
tmp1	<- tmp[as.character(ggtmp1$DV)]
tmp1[sel]	<- 0
ggtmp1$Estimate1	<- ggtmp1$Estimate + tmp1
tmp		<- ggtmp1[sel,"Std..Error"]
names(tmp)	<- ggtmp1[sel,"DV"]
tmp1	<- tmp[as.character(ggtmp1$DV)]
tmp1[sel]	<- 0
ggtmp1$StdErr	<- sqrt(ggtmp1$Std..Error^2 + tmp1^2)
if(make_plot){
	pdf(paste(plot.wd, "/graph_income_shr_het.pdf", sep=""), width = ww, height = ww)
	plots	<- ggplot(ggtmp1, aes(Retail, Estimate1, fill = IncGrp, col = IncGrp)) + 
				geom_bar(stat = "identity", position = position_dodge(width=0.9)) + 
				geom_errorbar(aes(ymin = Estimate1-1.96*StdErr, ymax = Estimate1 + 1.96 * StdErr, col = IncGrp), position=position_dodge(width=0.9), width=0.25) + 
				# geom_pointrange(aes(ymin = Estimate1-1.96*StdErr, ymax = Estimate1 + 1.96 * StdErr), position=position_dodge(width=0.9)) +
				scale_fill_manual(values = rev(incg.col)) + 
				scale_color_manual(values = rev(incg.col)) + 
				xlab("Dependent variable") + ylab("Estimate")+ coord_flip() + 
				guides(fill = guide_legend(reverse = TRUE), color = guide_legend(reverse = TRUE))
	print(plots)
	dev.off()
}



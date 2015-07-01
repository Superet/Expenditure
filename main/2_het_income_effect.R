library(ggplot2)
library(reshape2)
library(data.table)
library(plm)
library(gridExtra)
library(scales)
library(r2excel)

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

write2csv	<- TRUE
make_plot	<- TRUE
outfile		<- paste(plot.wd, "/2_income_effect.rdata", sep="")
outxls		<- paste(plot.wd, "/2_het_income_effect.xlsx", sep="")
mywb		<- createWorkbook()
sht1		<- createSheet(mywb, "Sheet1")
sht2		<- createSheet(mywb, "Regression")
sht3		<- createSheet(mywb, "Regression w price")

# Read data 
hh_exp		<- read.csv("2_hh_biweek_exp_merge_larger.csv")
fmt_attr 	<- read.csv("1_format_year_attr.csv")
price_dat	<- read.csv("1_format_biweek_price.csv")
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

recode <- function(old.idx, old.value, new.n){
	old.level 	<- sort(unique(old.idx))
	old.n		<- length(old.level)
	tab			<- sort(table(old.idx), decreasing=T)
	sort.x		<- as.numeric(names(tab))
	
	loop.n 		<- new.n 
	loop.x		<- old.idx
	loop.sx		<- sort.x
	mycut		<- NULL
	for(i in new.n:1){
		qt 		<- quantile(loop.x, c(0:loop.n)/loop.n)
		if(length(unique(qt))==loop.n+1){
			mycut <- c(mycut, qt)
			break
		}
		mycut	<- c(mycut, loop.sx[1])
		loop.x	<- loop.x[loop.x!=loop.sx[1]]
		loop.n	<- loop.n - 1
		loop.sx <- loop.sx[-1]
	}
	new.idx 	<- cut(old.idx, mycut, include.lowest=T, labels=1:new.n)
	new.level	<- tapply(old.value, new.idx, mean, na.rm=T)
	new.value	<- new.level[new.idx]
	cat("Table of original categorical variable vs recoded categorical variable:\n")
	print(table(old.idx, new.idx))
	return(cbind(new.idx, new.value) )
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

# Compute expenditure share
sel			<- paste("DOLP_", gsub("\\s", "_", fmt_name), sep="")
hh_exp$y 	<- rowSums(hh_exp[,sel])
hh_exp$income_real <- factor(hh_exp$income_real)
sel1		<- gsub("DOLP", "SHR", sel)
for(i in 1:length(fmt_name)){
	hh_exp[,sel1[i]] <- hh_exp[,sel[i]]/hh_exp$y
}

# Conver some date and factor variables
hh_exp$month <- month(as.Date("2004-1-1", format="%Y-%m-%d") + 14*(hh_exp$biweek-1))
hh_exp$income_group <- factor(hh_exp$income_group, levels=paste("Qt",1:4,sep=""))
hh_exp$famsize		<- factor(hh_exp$famsize, levels = c("Single","Two", "Three+"))
hh_exp$Q			<- hh_exp$food_quant + hh_exp$nonedible_quant
mycut 				<- c(0,1,2,3,5,14) + .5
hh_exp$d			<- as.numeric(cut(hh_exp$num_day, mycut, labels=1:5))
hh_exp[is.na(hh_exp$d),"d"] <- 0
hh_exp[is.na(hh_exp$y),"y"] <- 0
hh_exp[is.na(hh_exp$dol_purchases),"dol_purchases"] <- 0
hh_exp[is.na(hh_exp$num_day),"num_day"] <- 0
cpi 				<- unique(price_dat[,c("year","cpi")])

# Recode inocme 
hh_exp$income_real	<- as.numeric(as.character(hh_exp$income_real))
n_Inc		<- 8
sel			<-  hh_exp$income_real >= 27
hh_exp[sel,"income_real"] <- 27
tmp			<- recode(hh_exp$income_real, hh_exp$income_midvalue, n_Inc)
Inc_nodes	<- sort(unique(tmp[,2]))
hh_exp$income_nodes <- tmp[,1]

# Segment households based on their initial income level 
panelist	<- data.table(hh_exp)
setkeyv(panelist, c("household_code","year","biweek"))
panelist	<- panelist[,list(first_income = income_group[1], first_famsize = famsize[1]), by=list(household_code)]
hh_exp		<- merge(hh_exp, data.frame(panelist)[,c("household_code", "first_income")], by="household_code", all.x = TRUE)

###############################
# Characterize income changes # 
###############################
# The within-household trend of cpi-adjusted income 
pan_yr		<- data.table(hh_exp)
pan_yr		<- pan_yr[,list(Income = unique(income_midvalue)), by = list(first_income, household_code, year)]
pan_yr		<- merge(pan_yr, cpi, by="year", all.x=T)
pan_yr$Income	 <- pan_yr$Income / pan_yr$cpi
pan_yr$ln_income <- log(pan_yr$Income)
pan_yr$recession <- 1 * (pan_yr$year >= 2008)
pan_yr		<- data.frame(pan_yr)

# Within-household regression: income ~ year/recession
pan_yr$year <- factor(pan_yr$year)
summary(lm(Income ~ year*first_income, data = pan_yr))
myfit 		<- plm(Income ~ first_income*year, data=pan_yr, index = c("household_code","year"), model="within")
tmp			<- summary(myfit)
cat("Income trend regression:\n"); print(tmp); cat("\n")
myfit1		<- plm(Income ~ first_income*recession, data=pan_yr, index = c("household_code","year"), model="within")
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
tmp			<- unique(pan_yr[,c("first_income","year")])
tmp$year	<- factor(tmp$year)
tmpx		<- model.matrix(~first_income*year, data = tmp)[,-(1:4)]
identical(colnames(tmpx), names(coef(myfit)))
ggtmp		<- tmpx %*% coef(myfit)
ggtmp		<- cbind(tmp, ggtmp)
sel 		<- pan_yr$year == 2004
tmp			<- tapply(pan_yr[sel,"Income"], pan_yr[sel,"first_income"], mean)
ggtmp$Income<- (ggtmp$ggtmp + tmp[ggtmp$first_income])/1000
ggtmp$year	<- as.numeric(as.character(ggtmp$year))

plots		<- list(NULL)
plots[[1]]	<- ggplot(pan_yr, aes(as.numeric(year), Income, linetype = first_income)) + geom_smooth(method=lm, formula = y ~ x + I(x^2) + I(x^3)) + 
					guides(linetype = guide_legend(title = "Segment")) + 
					labs(tilte = "Smoothing line of raw data")
plots[[2]]	<- ggplot(ggtmp, aes(year, Income, linetype = first_income)) + geom_line() + 
			    		labs(x = "Year", y = "Income($1000)", title = "Smoothing line from within-model prediction") + 
						guides(linetype = guide_legend(title = "Segment"))
if(make_plot){
	pdf(paste(plot.wd,"/graph_income.pdf",sep=""), width = ww, height = ww*ar)
	print(plots[[1]])
	print(plots[[2]])
	dev.off()	
}

#------------------------------------------------------#
# Expectation and variance
# Data of current income and next-year income
tmp		<- data.table(hh_exp)
tmp.lab	<- c("Under 25k", "25k-29,999", "30k-39,999","40k-49,999","50k-59,999","60k-69,999","70k-99,999","100k+")
tmp$income_nodes	<- factor(tmp$income_nodes, levels=1:n_Inc, labels=tmp.lab)
tmp		<- tmp[, list(income_nodes = unique(income_nodes)), by=list(first_income, household_code, year)]
setkeyv(tmp, c("first_income","household_code","year"))
tmp		<- tmp[, income_forward := my_forward(income_nodes), by=list(first_income, household_code)]
tmp$income_forward	<- factor(tmp$income_forward, levels = 1:n_Inc, labels=tmp.lab)
tmp$seg	<- paste(tmp$first_income, tmp$year, sep="-")

# Compute the transition matrix of income for each segment
tmp		<- data.frame(tmp)
sel		<- tmp$year != 2010
tmp1	<- split(tmp[sel,], tmp[sel,"seg"])
tmp2	<- lapply(tmp1, function(x) mytransition(x$income_nodes, x$income_forward))
tmp.mat	<- rep(1, n_Inc) %*% t(Inc_nodes)

# Compute expected future income and uncertainty
tmp3	<- matrix(NA, length(tmp2), 2, dimnames = list(names(tmp1), c("Expectation", "SD")))
for(i in 1:length(tmp2)){
	sel			<- ifelse(i%%6==0, 1, i%%6)
	tmp.mat1	<- tmp.mat/cpi[sel,"cpi"]
	ee			<- rowSums(tmp2[[i]]*tmp.mat1)
	vv			<- rowSums(tmp2[[i]]*tmp.mat1^2) - ee^2
	tmp.wt		<- table(tmp1[[i]]$income_nodes)/nrow(tmp1[[i]])
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

#########################################
# Regressions to explore income effects # 
#########################################
#---------------------------------------#
# Construct regression data 
mydata 	<- hh_exp
tmp_dv		<- paste("SHR_", gsub("\\s", "_", fmt_name), sep="")
for(i in tmp_dv){
	mydata[,i] <- mydata[,i]*100
}
mydata$ln_income	<- log(mydata$income_midvalue)
mydata$month		<- factor(mydata$month)
mydata$ln_dolpurchases <- log(mydata$dol_purchases)
mydata$recession 	<- factor(mydata$recession)

# Add income uncertainty
tmp		<- data.table(hh_exp)
tmp$income_nodes	<- factor(tmp$income_nodes)
tmp		<- tmp[, list(income_nodes = unique(income_nodes)), by=list(first_income, household_code, year)]
tmp$recession <- ifelse(tmp$year >= 2008, 1, 0)
tmp		<- tmp[, income_forward := my_forward(income_nodes), by=list(household_code)]
tmp		<- data.frame(tmp)
tmp$income_forward <- as.factor(tmp$income_forward)

# Transition table: replace 0 diagal with 1 if nobody belongs to some income level.
tmp1	<- sort(unique(hh_exp$first_income))
tmpn	<- max(hh_exp$income_nodes)
trans.dat<- data.frame()
for(i in 1:length(tmp1)){
	for(j in 2004:2009){
		sel		<- tmp$first_income == tmp1[i] & tmp$year ==j
		tmp.tab	<- mytransition(tmp[sel,"income_nodes"], tmp[sel,"income_forward"])
		for(k in as.numeric(rownames(tmp.tab))){
			ss	<- ifelse(k==1, 1, k-1)
			ee	<- ifelse(k==tmpn, tmpn, k+1)
			trans.dat <- rbind(trans.dat, data.frame(first_income = tmp1[i], year = j, income_nodes = k, 
								p_low = sum(tmp.tab[as.character(k),1:ss]), p_high = sum(tmp.tab[as.character(k), ee:tmpn])))
		}
	}
}
mydata	<- merge(mydata, trans.dat, by = c("first_income","year","income_nodes"), all.x=T)


# Add the direction of income changes: increase or decrease compared with previous year
tmp		<- data.table(pan_yr)
setkeyv(tmp, c("household_code","year"))
tmp		<- tmp[, last_income:=c(NA, Income[-length(Income)]), by = list(household_code)]
tmp		<- tmp[, income_increase := 1*(Income - last_income>0)]

mydata	<- merge(mydata, data.frame(tmp)[,c("household_code","year","income_increase")], by=c("household_code","year"), all.x=T)
mydata$income_increase	<- factor(mydata$income_increase)
table(mydata$income_increase)
mydata$year	<- as.factor(mydata$year)

#---------------------------------------#
# Set up DV and IV, regression models 
dv_vec	<- c("ln_dolpurchases", paste("SHR_", gsub("\\s", "_", fmt_name), sep="") )
myfml	<- list(Homogeneity		= as.formula("y ~ ln_income + year + month"), 
				Heterogeneity	= as.formula("y ~ ln_income*first_income + year + month"),
				Assymetry		= as.formula("y ~ ln_income*income_increase + year + month"),
				Het_Asy			= as.formula("y ~ ln_income*first_income*income_increase + year + month"),
				Expectation		= as.formula("y ~ ln_income + p_low + p_high + year + month"),
				Het_expc		= as.formula("y ~ ln_income + p_low*first_income + p_high*first_income + year + month")
				)

# Run regressions
regrs	<- data.frame()
for(i in 1:length(dv_vec)){
	prc 		<- proc.time()
	for(j in 1:length(myfml)){
		tmp 	<- as.formula(substitute(y ~ x, list(y = as.name(dv_vec[i]), x = terms(myfml[[j]])[[3]])) )
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
tmpls		<- vector("list", length(dv_vec))
names(tmpls)<- dv_vec
for(i in 1:length(dv_vec)){
	tmp1	<- subset(tmp, DV == dv_vec[i])
	model_list	<- split(tmp1, tmp1$model)
	tmp.tab	<- model_outreg(model_list, p.given = TRUE, head.name = c("Estimate","Std..Error","Pr...t..","Var"), digits = 4)
	tmpls[[i]]	<- tmp.tab
}

if(write2csv){
	for(i in 1:length(tmpls)){
		xlsx.addHeader(mywb, sht2, value = names(tmpls)[i], level = 3)
		xlsx.addTable(mywb, sht2, tmpls[[i]], row.names = F)
	}
}

#################################
# Add prices to the regressions # 
#################################
# Append price data to the regression data 
tmp		<- data.table(price_dat)
tmp		<- tmp[,list(bsk_price_paid_2004 = mean(bsk_price_paid_2004)), by = list(scantrack_market_descr, biweek, channel_type)] 
tmp		<- dcast(data.frame(tmp), scantrack_market_descr + biweek ~ channel_type, value.var = "bsk_price_paid_2004")
names(tmp)[-(1:2)]	<- paste("PRC_", names(tmp)[-(1:2)], sep="")
names(tmp)	<- gsub("\\s", "_", names(tmp))
mydata	<- merge(mydata, tmp, by = c("scantrack_market_descr","biweek"), all.x=T)

# Re estimate the regressions with price controls 
myfml	<- list(Homogeneity		= as.formula("y ~ ln_income + year + month + PRC_Discount_Store + PRC_Grocery + PRC_Warehouse_Club"), 
				Heterogeneity	= as.formula("y ~ ln_income*first_income + year + month + PRC_Discount_Store + PRC_Grocery + PRC_Warehouse_Club"),
				Assymetry		= as.formula("y ~ ln_income*income_increase + year + month + PRC_Discount_Store + PRC_Grocery + PRC_Warehouse_Club"),
				Het_Asy			= as.formula("y ~ ln_income*first_income*income_increase + year + month + PRC_Discount_Store + PRC_Grocery + PRC_Warehouse_Club"),
				Expectation		= as.formula("y ~ ln_income + p_low + p_high + year + month + PRC_Discount_Store + PRC_Grocery + PRC_Warehouse_Club"),
				Het_expc		= as.formula("y ~ ln_income + p_low*first_income + p_high*first_income + year + month + PRC_Discount_Store + PRC_Grocery + PRC_Warehouse_Club")
				)

# Run regressions
regrs1	<- data.frame()
for(i in 1:length(dv_vec)){
	prc 		<- proc.time()
	for(j in 1:length(myfml)){
		tmp 	<- as.formula(substitute(y ~ x, list(y = as.name(dv_vec[i]), x = terms(myfml[[j]])[[3]])) )
		myfit	<- plm(tmp, data = mydata, index = c("household_code","biweek"), model="within")
		tmp2	<- summary(myfit)
		tmp3	<- data.frame(model = names(myfml)[j], DV = dv_vec[i], tmp2$coefficients)
		tmp3$Var<- rownames(tmp3)
		rownames(tmp3) <- NULL
		regrs1	<- rbind(regrs1, tmp3)
	}
	use.time 	<- proc.time() - prc
	cat("Regressions for", dv_vec[i], "finishes, using", use.time[3]/60, "min.\n")
}

# Print regression results 
tmp			<- subset(regrs1, substr(regrs1$Var, 1, 5) != "month" & substr(regrs1$Var, 1, 4) != "year")
tmpls		<- vector("list", length(dv_vec))
names(tmpls)<- dv_vec
for(i in 1:length(dv_vec)){
	tmp1	<- subset(tmp, DV == dv_vec[i])
	model_list	<- split(tmp1, tmp1$model)
	tmp.tab	<- model_outreg(model_list, p.given = TRUE, head.name = c("Estimate","Std..Error","Pr...t..","Var"), digits = 4)
	tmpls[[i]]	<- tmp.tab
}

if(write2csv){
	for(i in 1:length(tmpls)){
		xlsx.addHeader(mywb, sht3, value = names(tmpls)[i], level = 3)
		xlsx.addTable(mywb, sht3, tmpls[[i]], row.names = F)
	}
}

saveWorkbook(mywb, outxls)

save.image(file = outfile)
cat("This program is done.")

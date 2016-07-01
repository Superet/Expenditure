library(ggplot2)
library(reshape2)
library(data.table)
library(plm)
library(gridExtra)
library(scales)
library(systemfit)
library(stargazer)
library(gridExtra)

options(error = quote({dump.frames(to.file = TRUE)}))

# setwd("~/Documents/Research/Store switching/processed data")
# plot.wd	<- '/Users/chaoqunchen/Desktop'
# source('~/Documents/Research/Store switching/Exercise/main/income_effect/plm_vcovHC.R')

# setwd("/home/brgordon/ccv103/Exercise/run")
setwd("/kellogg/users/marketing/2661703/Expenditure")
# setwd("/sscc/home/c/ccv103/Exercise/run")
plot.wd 	<- paste(getwd(), "/results", sep="")
ww			<- 6.5
ww1			<- 8.5
ar 			<- .8	#.6

week.price	<- FALSE
cpi.adj		<- TRUE
write2csv	<- FALSE
make_plot	<- TRUE
fname		<- "expenditure_reg_exp"
if(cpi.adj) { fname <- paste(fname, "_cpi", sep="")}
if(week.price){ fname <- paste(fname, "_wkprc", sep="")}
# outxls		<- paste(plot.wd, "/", fname, "_", gsub("-", "", as.character(Sys.Date())), ".xlsx", sep="")
# mywb		<- createWorkbook()
# sht1		<- createSheet(mywb, "Regression")
# sht2		<- createSheet(mywb, "SUR")

if(week.price){
	load("hh_biweek_exp_20150812.rdata")
}else{
	load("hh_biweek_exp.rdata")
}
codebook	<- read.csv("code_book.csv")
source("plm_vcovHC.R")

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

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

#################
# Organize data # 
#################
# Retail formats 
fmt_name 	<- as.character(sort(unique(fmt_attr$channel_type)))
R			<- length(fmt_name)

# Conver some date and factor variables
if(week.price){
	# Segment households based on their initial income level 
	panelist	<- data.table(hh_exp)
	setkeyv(panelist, c("household_code","year","biweek"))
	panelist	<- panelist[,list(income = first_income[1], first_famsize = famsize[1]), by=list(household_code)]
	tmp			<- quantile(panelist$income, c(0, .33, .67, 1))
	num_grp		<- 3
	panelist	<- panelist[, first_incomeg := cut(panelist$income, tmp, labels = paste("T", 1:num_grp, sep=""), include.lowest = T)]
	hh_exp		<- merge(hh_exp, data.frame(panelist)[,c("household_code", "first_incomeg")], by = "household_code", all.x=T )
	hh_exp$first_incomeg	<- factor(hh_exp$first_incomeg, levels = c("T1", "T2", "T3"))
	cat("Table of initial income distribution:\n"); print(table(panelist$first_incomeg)); cat("\n")
	cat("Table of segments in the expenditure data:\n"); print(table(hh_exp$first_incomeg)); cat("\n")
}
hh_exp$month <- month(as.Date("2004-1-1", format="%Y-%m-%d") + 14*(hh_exp$biweek-1))
hh_exp$famsize		<- factor(hh_exp$famsize, levels = c("Single","Two", "Three+"))
hh_exp$condo		<- factor(hh_exp$condo, levels = c(0, 1))
hh_exp$employed		<- factor(hh_exp$employed, levels = c(0, 1))
hh_exp$first_incomeg	<- as.character(hh_exp$first_incomeg)
hh_exp$first_incomeg	<- factor(hh_exp$first_incomeg, levels = paste("T", 1:3, sep=""), labels = c("Low", "Med", "High"))
cat("Table of segments in the expenditure data:\n"); print(table(hh_exp$first_incomeg)); cat("\n")

# Compute expenditure share
sel			<- paste("DOL_", gsub("\\s", "_", fmt_name), sep="")
sel1		<- gsub("DOL", "SHR", sel)
for(i in 1:length(fmt_name)){
	hh_exp[,sel1[i]] <- hh_exp[,sel[i]]/hh_exp$dol
}

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

# --------------- #
# Regression data #
mydata 	<- hh_exp
if(cpi.adj){
	mydata$income_midvalue	<- mydata$income_midvalue/mydata$cpi
	mydata$stone_price		<- mydata$stone_price/mydata$cpi
}	
mydata$ln_income	<- log(mydata$income_midvalue)
sum(mydata$dol == 0 )/nrow(mydata)
annual.week	<- 26
mydata$week_income	<- mydata$income_midvalue/annual.week
mydata$month		<- factor(mydata$month)
mydata$ln_dol 		<- log(mydata$dol)
mydata$recession 	<- factor(mydata$recession)

# Add price 
if(week.price){
	tmp		<- dcast(price_dat, scantrack_market_descr+biweek ~ channel_type, value.var = "bsk_price_paid_2004")
	colnames(tmp)	<- c("scantrack_market_descr", "biweek", paste("PRC_", gsub(" ", "_", fmt_name), sep=""))
	mydata 	<- merge(mydata, tmp, by = c("scantrack_market_descr", "biweek"), all.x = T)
}else{
	tmp		<- dcast(price_dat, scantrack_market_descr+year ~ channel_type, value.var = "bsk_price_paid_2004")
	colnames(tmp)	<- c("scantrack_market_descr", "year", paste("PRC_", gsub(" ", "_", fmt_name), sep=""))
	dim(mydata)
	mydata 	<- merge(mydata, tmp, by = c("scantrack_market_descr", "year"), all.x = T)
}
dim(mydata)
mydata$year			<- as.factor(mydata$year)

# Shopping incidence
tmp_dv		<- paste("SHR_", gsub("\\s", "_", fmt_name), sep="")
# for(i in tmp_dv){
# 	mydata[,i] <- mydata[,i]*100
# }
sel			<- paste("DOL_", gsub("\\s", "_", fmt_name), sep="")
sel1		<- gsub("DOL", "IC", sel)
for(i in 1:length(fmt_name)){
	mydata[,sel1[i]] <- 1*(mydata[,sel[i]]>0)
}

# Add lagged purchases
mydata	<- data.table(mydata)
setkeyv(mydata, c("household_code", "biweek"))
mydata	<- mydata[,lag_dol := my_lag(dol), by = list(household_code)]
mydata	<- data.frame(mydata)

# mydata	<- subset(mydata, dol > 0)
# mydata	<- subset(mydata, !is.na(lag_dol))
mydata$ln_income_low <- 1*(mydata$first_incomeg == "Low")*mydata$ln_income
mydata$ln_income_med <- 1*(mydata$first_incomeg == "Med")*mydata$ln_income
mydata$ln_income_high <- 1*(mydata$first_incomeg == "High")*mydata$ln_income
mydata$hh_year			<- paste(mydata$household_code, mydata$year, sep = "*")

# Only keep relevant columns
mydata	<- mydata[, c("household_code", "biweek","dol", "ln_dol",  
					paste("DOL_", gsub("\\s", "_", fmt_name), sep=""), 
					paste("SHR_", gsub("\\s", "_", fmt_name), sep=""), 
					paste("PRC_", gsub("\\s", "_", fmt_name), sep=""),
					paste("IC_", gsub("\\s", "_", fmt_name), sep=""), 
					"first_incomeg", "ln_income", "ln_income_low", "ln_income_med", "ln_income_high", "week_income",
					"year", "month", "lag_dol", "hh_year",
					"stone_price", "household_size", "condo", "employed", "NumChild")]

cat("Number of observations with 0 expenditure =", sum(mydata$dol==0), ", or, ", round(sum(mydata$dol==0)/nrow(mydata), 4), ".\n")

############################################
# Single-equation fixed effect regressions #
############################################
#---------------------------------------#
# Set up DV and IV, regression models 
dv_vec	<- "ln_dol"
myfml	<- list(Homogeneity		= as.formula("y ~ ln_income + month"), 
				HomoYear		= as.formula("y ~ ln_income + year + month"), 
				HomoPrice		= as.formula("y ~ ln_income + year + month + stone_price"),
				Heterogeneity	= as.formula("y ~ ln_income*first_incomeg + year + month"),
				HeteroYear		= as.formula("y ~ ln_income*first_incomeg + month"),
				HeteroPrice		= as.formula("y ~ ln_income*first_incomeg + year + month + stone_price"))

# Run regressions
regrs	<- data.frame()
prc 	<- proc.time()
sel		<- mydata$dol > 0 		# Only focus on positive expenditure 
reg.ls	<- setNames(vector("list", length(myfml)), names(myfml))
for(j in 1:length(myfml)){
	rm(list = intersect(ls(), "myfit"))
	tmp 	<- as.formula(substitute(y ~ x, list(y = as.name(dv_vec), x = terms(myfml[[j]])[[3]])) )
	if(!dv_vec %in% c("ln_dol", "dol")){
		tmp		<- update(tmp, . ~ . + lag_dol)
	}
	myfit	<- plm(tmp, data = mydata[sel,], index = c("household_code","biweek"), model="within")
	reg.ls[[j]]	<- myfit
	tmp2	<- summary(myfit)
	tmp3	<- data.frame(model = names(myfml)[j], DV = dv_vec, tmp2$coefficients)
	tmp3$Var<- rownames(tmp3)
	rownames(tmp3) <- NULL
	regrs	<- rbind(regrs, tmp3)
}

# Compute cluster-robust se
tmp.coef	<- lapply(reg.ls, coeftest)
unc.se		<- lapply(tmp.coef, function(x) x[, "Std. Error"])
cls.se1		<- lapply(reg.ls, function(x) sqrt(diag(vcovHC.plm.new(x, method = "arellano", type = "HC1", cluster = "define", clustervec = mydata[sel,"household_code"]))))
cls.se2		<- lapply(reg.ls, function(x) sqrt(diag(vcovHC.plm.new(x, method = "arellano", type = "HC1", cluster = "define", clustervec = mydata[sel,"hh_year"]))))

# Put together the se 
tmpidx		<- lapply(unc.se, function(x) grep("ln_income", names(x)))
tmp			<- lapply(1:length(tmpidx), function(i) cbind(Uncorrected = unc.se[[i]][tmpidx[[i]]], 
														HHyr.cluster.se = cls.se2[[i]][tmpidx[[i]]],
														HH.cluster.se = cls.se1[[i]][tmpidx[[i]]] ))
names(tmp)	<- names(reg.ls)
cat("Compare cluster robust standard error for regressions of log(expenditure) ~ log(income):\n"); print(tmp); cat("\n")

# Export to HTML
stargazer(reg.ls, type = "html", title = "Single-equation fixed effects regression", 
		 align = TRUE, no.space = TRUE, 
	     se = cls.se1, omit.stat = "f",
		 omit = c("year", "month"), order=c("ln_income"),
		 covariate.labels = c("log(I)", "log(I)*Med", "log(I)*High", "Stone price"), 
		 add.lines = list(c("Year FE", "N","Y","Y","N","Y","Y"), c("Month FE","Y","Y","Y","Y","Y","Y"), c("HH FE", "Y","Y","Y","Y","Y","Y")),
		 notes = c("se clustered over households"),
		 out = paste(plot.wd, "/", fname, "_se_", dv_vec, ".html", sep=""))

# Export to Latex
stargazer(reg.ls, type = "latex", title = "Single-equation fixed effects regression", 
		 align = TRUE, no.space = TRUE, 
		 se = cls.se1, omit.stat = "f",
		 omit = c("year", "month"), order=c("ln_income"), 
		 covariate.labels = c("log(I)", "log(I)*Med", "log(I)*High", "Stone price"), 
		 add.lines = list(c("Year FE", "N","Y","Y","N","Y","Y"), c("Month FE","Y","Y","Y","Y","Y","Y"), c("HH FE", "Y","Y","Y","Y","Y","Y")),
		 out = paste(plot.wd, "/", fname, "_se_", dv_vec, ".tex", sep=""))		

use.time 	<- proc.time() - prc
cat("Regressions for log(dol) finishes, using", use.time[3]/60, "min.\n")

#############################################################################
# Single-equation fixed effect regressions measureing propensity to consume #
#############################################################################
#---------------------------------------#
# Set up DV and IV, regression models 
dv_vec	<- "dol"
myfml	<- list(Homogeneity		= as.formula("y ~ week_income + month"), 
				HomoYear		= as.formula("y ~ week_income + year + month"), 
				HomoPrice		= as.formula("y ~ week_income + year + month + stone_price"), 
				Heterogeneity	= as.formula("y ~ week_income*first_incomeg + month"),
				HeteroYear		= as.formula("y ~ week_income*first_incomeg + year + month"),
				HeteroPrice		= as.formula("y ~ week_income*first_incomeg + year + month + stone_price")
				)

# Run regressions
regrs.ptc	<- data.frame()
prc 		<- proc.time()
sel			<- mydata$dol > 0
reg.linear.ls	<- setNames(vector("list", length(myfml)), names(myfml))
for(j in 1:length(myfml)){
	rm(list = "myfit")
	tmp 	<- as.formula(substitute(y ~ x, list(y = as.name(dv_vec), x = terms(myfml[[j]])[[3]])) )
	if(!dv_vec %in% c("ln_dol", "dol")){
		tmp		<- update(tmp, . ~ . + lag_dol)
	}
	myfit	<- plm(tmp, data = mydata[sel,], index = c("household_code","biweek"), model="within")
	reg.linear.ls[[j]]	<- myfit
	tmp2	<- summary(myfit)
	tmp3	<- data.frame(model = names(myfml)[j], DV = dv_vec, tmp2$coefficients)
	tmp3$Var<- rownames(tmp3)
	rownames(tmp3) <- NULL
	regrs.ptc	<- rbind(regrs.ptc, tmp3)
}

# Compute cluster-robust se
tmp.coef	<- lapply(reg.linear.ls, coeftest)
unc.se		<- lapply(tmp.coef, function(x) x[, "Std. Error"])
cls.se1		<- lapply(reg.linear.ls, function(x) sqrt(diag(vcovHC.plm.new(x, method = "arellano", type = "HC1", cluster = "define", clustervec = mydata[sel,"household_code"]))))
cls.se2		<- lapply(reg.linear.ls, function(x) sqrt(diag(vcovHC.plm.new(x, method = "arellano", type = "HC1", cluster = "define", clustervec = mydata[sel,"hh_year"]))))

# Put together the se 
tmpidx		<- lapply(unc.se, function(x) grep("week_income", names(x)))
tmp			<- lapply(1:length(tmpidx), function(i) cbind(Uncorrected = unc.se[[i]][tmpidx[[i]]], 
														HHyr.cluster.se = cls.se2[[i]][tmpidx[[i]]],
														HH.cluster.se = cls.se1[[i]][tmpidx[[i]]]))
names(tmp)	<- names(reg.linear.ls)
cat("Compare cluster robust standard error for regressions of expenditure dol ~ week_income:\n"); print(tmp); cat("\n")

# Export estimation table
stargazer(reg.linear.ls, type = "html", title = "Single-equation fixed effects regression of expenditure on weekly income", 
		 align = TRUE, no.space = TRUE, se = cls.se1, omit.stat = "f",
		 omit = c("year", "month"), order=c("week_income"), 
		 covariate.labels = c("I", "I*Med", "I*High", "Stone price"), 
		 add.lines = list(c("Year FE","N","Y","Y","N","Y","Y"), c("Month FE","Y","Y","Y","Y","Y","Y"), c("HH FE", "Y","Y","Y","Y","Y","Y")), 
		 notes = c("se clustered over households"),
		 out = paste(plot.wd, "/", fname, "_septc_", dv_vec, ".html", sep=""))

stargazer(reg.linear.ls, type = "latex", title = "Single-equation fixed effects regression of expenditure on weekly income", 
		 align = TRUE, no.space = TRUE, se = cls.se1, omit.stat = "f",
		 omit = c("year", "month"), order=c("week_income"), 
		 covariate.labels = c("I", "I*Med", "I*High", "Stone price"), 
		 add.lines = list(c("Year FE","N","Y","Y","N","Y","Y"), c("Month FE","Y","Y","Y","Y","Y","Y"), c("HH FE", "Y","Y","Y","Y","Y","Y")), 
		 out = paste(plot.wd, "/", fname, "_septc_", dv_vec, ".tex", sep=""))		

use.time 	<- proc.time() - prc
cat("Regressions for", dv_vec, "finishes, using", use.time[3]/60, "min.\n")

# Save results 
rm(list = intersect(ls(), c("tmp", "tmp1", "tmp2", "tmp3", "tmp4", "use.time", "make_plot", "prc","tmp.coef", "tmp_dv", 
							"i", "j", "sel", "sel1", "myfml", "tmpidx")))

save.image(file = paste(plot.wd, "/", fname, "_", as.character(Sys.Date()), ".rdata", sep=""))

cat("This program is done.")
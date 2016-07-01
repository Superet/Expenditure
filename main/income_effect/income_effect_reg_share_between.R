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
# setwd("/kellogg/users/marketing/2661703/Expenditure")
# setwd("/sscc/home/c/ccv103/Exercise/run")
setwd("U:/Users/ccv103/Documents/Research/Store switching/run")
plot.wd 	<- paste(getwd(), "/results", sep="")
ww			<- 6.5
ww1			<- 8.5
ar 			<- .5	

week.price	<- FALSE
cpi.adj		<- TRUE
write2csv	<- FALSE
make_plot	<- TRUE
fname		<- "expenditure_reg_share_btw"
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
mydata$famsize_two	<- 1*(mydata$famsize == "Two")
mydata$famsize_three	<- 1*(mydata$famsize == "Three+")
mydata$hh_year			<- paste(mydata$household_code, mydata$year, sep = "*")

# Only keep relevant columns
mydata	<- mydata[, c("household_code", "biweek","dol", "ln_dol",  
					paste("DOL_", gsub("\\s", "_", fmt_name), sep=""), 
					paste("SHR_", gsub("\\s", "_", fmt_name), sep=""), 
					paste("PRC_", gsub("\\s", "_", fmt_name), sep=""),
					paste("IC_", gsub("\\s", "_", fmt_name), sep=""), 
					"first_incomeg", "ln_income", "ln_income_low", "ln_income_med", "ln_income_high", "week_income",
					"year", "month", "lag_dol", "hh_year",
					"stone_price", "famsize", "condo", "employed", "NumChild", "famsize_two", "famsize_three")]

cat("Number of observations with 0 expenditure =", sum(mydata$dol==0), ", or, ", round(sum(mydata$dol==0)/nrow(mydata), 4), ".\n")

#####################################
# SUR regression with fixed effects #
#####################################
my.SURFE <- function(syseq, data, panid, model = 'within',...){
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
	if(model == "within"){
		for(i in 1:length(varvec)){
			var		<- as.name(varvec[i])
			mydata 	<- mydata[,eval(var):=eval(var) - mean(eval(var), na.rm=T), by = list(id)]
		}
	}else if(model == "between"){
		for(i in 1:length(varvec)){
			var		<- as.name(varvec[i])
			mydata 	<- mydata[,eval(var):=mean(eval(var), na.rm=T), by = list(id)]
		}
		mydata	<- unique(mydata, by = "id")
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
	df		<- myfit$df.residual/length(syseq)
	coeftab$adj_df	<- switch(model, 
							within = df - length(unique(data[,panid])), 
							between = df, 
							pooling = df)

	# Compute cluster-robust standard error 
	# NOTE: clustered standard error has correctted degree of freedom explicity
	cls.cov		<- NA
	if(model != "between"){
		cls.cov		<- lapply(myfit$eq, function(x) vcovHC.sur(x, method="arellano", type = "HC1", clustervec = mydata$id))
		cls.se		<- lapply(cls.cov, function(x) sqrt(diag(x)))
		coeftab$cls_se		<- 	unlist(cls.se)
	}else{
		coeftab$cls_se		<- coeftab[,"Std..Error"]		# For between estimator, clustered se. = se.
	}
	
	return(list(coeftab = coeftab, fit = myfit, cls.cov = cls.cov))
}

# Fit SUR with fixed effects 
dv_mat	<- paste("SHR_", gsub("\\s", "_", fmt_name), sep="")
dv_mat	<- dv_mat[-1]
myfml	<- list(Homogeneity		= as.formula("y ~ ln_income + month + lag_dol+ famsize_two + famsize_three + NumChild"), 
				HomoPrice		= as.formula("y ~ ln_income + month + year + lag_dol + famsize_two + famsize_three + NumChild +
													PRC_Convenience_Store + PRC_Discount_Store + PRC_Dollar_Store + 
													PRC_Drug_Store + PRC_Grocery + PRC_Warehouse_Club"), 
				Heterogeneity 	= as.formula("y ~ ln_income + ln_income_med + ln_income_high + month + lag_dol + famsize_two + famsize_three + NumChild"), 
				HeterPrice 		= as.formula("y ~ ln_income + ln_income_med + ln_income_high + lag_dol + year + month + famsize_two + famsize_three + NumChild +
													PRC_Convenience_Store + PRC_Discount_Store + PRC_Dollar_Store + 
													PRC_Drug_Store + PRC_Grocery + PRC_Warehouse_Club")								
				)
modn	<- c("pooling", "between", "within")
				
regrs1 	<- data.frame()
regrs1.v<- regrs1.clsv	<- setNames(vector("list", length = length(myfml)*length(modn)), 
									paste(rep(names(myfml), each = length(modn)), modn, sep="*"))
panid	<- "household_code"
for(j in 1:length(myfml)){
	rm(list = "tmp")
	prc 		<- proc.time()
	syseq	<- lapply(dv_mat, function(yn) 
				as.formula(substitute(y ~ x, list(y = as.name(yn), x = terms(myfml[[j]])[[3]])) ) )
	names(syseq)	<- gsub("_", ".", dv_mat)
	
	for(k in 1:length(modn)){
		tmp		<- try(my.SURFE(syseq, mydata, panid = "household_code", model = modn[k],maxiter = 2, residCovWeighted = TRUE))
		# Save the coefficient estimates 
		tmpidx	<- (j-1)*length(modn) + k
		regrs1.v[[tmpidx]]	<- vcov(tmp$fit)								# Covariance of se
		if(modn[k] == "between"){
			regrs1.clsv[[tmpidx]]<- regrs1.v[[tmpidx]]			
		}else{
			regrs1.clsv[[tmpidx]]<- setNames(tmp$cls.cov, dv_mat)				# Covariance of clustered se.
		}
		tmp		<- tmp$coeftab
		tmp1 	<- sapply(rownames(tmp), function(x) strsplit(x, "_")[[1]][1])
		tmp2	<- sapply(1:length(tmp1), function(i) gsub(paste(tmp1[i], "_", sep=""), "", rownames(tmp)[i]))
		tmp3	<- cbind(model = names(myfml)[j], DV = tmp1, method = modn[k], Var = tmp2, tmp)
		rownames(tmp3) <- NULL
		regrs1	<- rbind(regrs1, tmp3)
	}
}

# Organize the parameter table 
tmp		<- regrs1[grep("ln_income", regrs1$Var),]
tmp1	<- setNames(fmt_name, paste("SHR.", gsub("\\s", "\\.", fmt_name), sep=""))
tmp$DV1	<- tmp1[as.character(tmp$DV)]
tmp1	<- setNames(c("log(I)", "log(I)*Med", "log(I)*High"),c("ln_income", "ln_income_med", "ln_income_high"))
tmp$Var1<- tmp1[as.character(tmp$Var)]
tmp$model1	<- factor(paste(tmp$model, tmp$method, sep="*"), levels = paste(rep(names(myfml), each = length(modn)), modn, sep="*"))

tmp1	<- split(tmp, tmp$model1)
for(i in 1:length(tmp1)){
	tmp			<- cbind(tmp1[[i]]$Estimate, tmp1[[i]]$cls_se, tmp1[[i]]$Estimate/tmp1[[i]]$cls_se)
	tmp			<- cbind(tmp, pt(-abs(tmp[,3]), tmp1[[i]]$adj_df)*2)
	dimnames(tmp)	<- list(paste(tmp1[[i]]$DV1, tmp1[[i]]$Var1, sep=":"), c("Estimate", "Std. Error", "t value", "Pr(>|t|)")) 
	class(tmp)	<- "coeftest"
	tmp1[[i]]	<- tmp
}

stargazer(tmp1, type = "html", title = "Summary of SUR regressions", dep.var.labels = "Expenditure share",
		 align = TRUE, no.space = TRUE, column.labels = rep(modn, length(myfml)), 
	 	 notes = c("se clustered over households"), 
		add.lines = list(c("Year FE", rep(c("N", "Y", "N", "Y"), each = length(modn))), c("Month FE",rep("Y",length(tmp1))), 
						 c("HH FE", rep(c("N", "N", "Y"), length(myfml))), c("Price", rep(c("N", "Y", "N", "Y"), each = length(modn)))), 
		 out = paste(plot.wd, "/", fname, "_sur_share_sum.html", sep=""))

stargazer(tmp1, type = "latex", title = "Summary of SUR regressions", dep.var.labels = "Expenditure share",
		 align = TRUE, no.space = TRUE, column.labels = rep(modn, length(myfml)), 
	 	 notes = c("se clustered over households"), 
		add.lines = list(c("Year FE", rep(c("N", "Y", "N", "Y"), each = length(modn))), c("Month FE",rep("Y",length(tmp1))), 
						 c("HH FE", rep(c("N", "N", "Y"), length(myfml))), c("Price", rep(c("N", "Y", "N", "Y"), each = length(modn)))), 
		 out = paste(plot.wd, "/", fname, "_sur_share_sum.tex", sep=""))

# Save results 
rm(list = intersect(ls(), c("tmp", "tmp1", "tmp2", "tmp3", "tmp4", "use.time", "make_plot", "prc","tmp.coef", "tmp_dv", 
							"i", "j", "sel", "sel1", "myfml", "tmpidx")))

save.image(file = paste(plot.wd, "/", fname, "_", as.character(Sys.Date()), ".rdata", sep=""))

cat("This program is done.")
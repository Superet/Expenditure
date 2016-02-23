library(ggplot2)
library(reshape2)
library(data.table)
library(plm)
library(gridExtra)
library(scales)
library(r2excel)
library(systemfit)
library(stargazer)

options(error = quote({dump.frames(to.file = TRUE)}))

# setwd("~/Documents/Research/Store switching/processed data")
# plot.wd	<- '/Users/chaoqunchen/Desktop'
# source('~/Documents/Research/Store switching/Exercise/main/outreg function.R')

# setwd("/home/brgordon/ccv103/Exercise/run")
setwd("/kellogg/users/marketing/2661703/Expenditure")
# setwd("/sscc/home/c/ccv103/Exercise/run")
plot.wd 	<- paste(getwd(), "/results", sep="")
ww			<- 6.5
ww1			<- 12
ar 			<- .8	#.6

cpi.adj		<- FALSE
write2csv	<- TRUE
make_plot	<- TRUE
fname		<- "expenditure_reg"
if(cpi.adj) { fname <- paste(fname, "_cpi", sep="")}
outxls		<- paste(plot.wd, "/", fname, "_", gsub("-", "", as.character(Sys.Date())), ".xlsx", sep="")
mywb		<- createWorkbook()
sht1		<- createSheet(mywb, "Regression")
sht2		<- createSheet(mywb, "SUR")

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

#################
# Organize data # 
#################
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

# ------------------------------------------------------ #
# Segment households based on their initial income level #
panelist	<- data.table(hh_exp)
setkeyv(panelist, c("household_code","year","biweek"))
panelist	<- panelist[,list(income = first_income[1], first_famsize = famsize[1]), by=list(household_code)]
tmp			<- quantile(panelist$income, c(0, .33, .67, 1))
num_grp		<- 3
# panelist	<- panelist[, first_incomeg := cut(panelist$income, tmp, labels = paste("T", 1:num_grp, sep=""), include.lowest = T)]
panelist	<- panelist[, first_incomeg := cut(panelist$income, tmp, labels = c("Low", "Med", "High"), include.lowest = T)]
hh_exp		<- merge(hh_exp, data.frame(panelist)[,c("household_code", "first_incomeg")], by = "household_code", all.x=T )
cat("Table of initial income distribution:\n"); print(table(panelist$first_incomeg)); cat("\n")
cat("Table of segments in the expenditure data:\n"); print(table(hh_exp$first_incomeg)); cat("\n")

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
	mydata$ln_income	<- log(mydata$income_midvalue/mydata$cpi)
}else{
	mydata$ln_income	<- log(mydata$income_midvalue)
}
sum(mydata$dol == 0 )/nrow(mydata)
mydata$month		<- factor(mydata$month)
mydata$ln_dol 		<- log(mydata$dol)
mydata$recession 	<- factor(mydata$recession)
mydata$year			<- as.factor(mydata$year)

# Add price 
tmp		<- dcast(price_dat, scantrack_market_descr+biweek ~ channel_type, value.var = "bsk_price_paid_2004")
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

# Only keep relevant columns
mydata	<- mydata[, c("household_code", "biweek","dol", "ln_dol",  
					paste("DOL_", gsub("\\s", "_", fmt_name), sep=""), 
					paste("SHR_", gsub("\\s", "_", fmt_name), sep=""), 
					paste("PRC_", gsub("\\s", "_", fmt_name), sep=""),
					paste("IC_", gsub("\\s", "_", fmt_name), sep=""), 
					"first_incomeg", "ln_income", "ln_income_low", "ln_income_med", "ln_income_high",
					"year", "month", "lag_dol")]

############################################
# Single-equation fixed effect regressions #
############################################
#---------------------------------------#
# Set up DV and IV, regression models 
dv_vec	<- c("dol", "ln_dol", paste("SHR_", gsub("\\s", "_", fmt_name), sep=""), paste("IC_", gsub("\\s", "_", fmt_name), sep="") )
myfml	<- list(Homogeneity		= as.formula("y ~ ln_income + year + month"), 
				Heterogeneity	= as.formula("y ~ ln_income*first_incomeg + month"), 
				Year			= as.formula("y ~ ln_income*first_incomeg + year + month"), 
				Price			= as.formula("y ~ ln_income*first_incomeg + 
											PRC_Convenience_Store + PRC_Discount_Store + PRC_Dollar_Store + 
											PRC_Drug_Store + PRC_Grocery + PRC_Warehouse_Club + month + year")
				)

# Run regressions
regrs	<- data.frame()
for(i in 1:length(dv_vec)){
	prc 		<- proc.time()
	if(dv_vec[i] != "dol"){
		sel	<- mydata$dol > 0
	}else{
		if(dv_vec[i] != "ln_dol"){
			sel <- !is.na(mydata$lag_dol)
		}else{
			sel	<- rep(TRUE, nrow(mydata))
		}
	}
	reg.ls	<- setNames(vector("list", 4), names(myfml))
	for(j in 1:length(myfml)){
		rm(list = "myfit")
		tmp 	<- as.formula(substitute(y ~ x, list(y = as.name(dv_vec[i]), x = terms(myfml[[j]])[[3]])) )
		if(!dv_vec[i] %in% c("ln_dol", "dol")){
			tmp		<- update(tmp, . ~ . + lag_dol)
		}
		myfit	<- plm(tmp, data = mydata[sel,], index = c("household_code","biweek"), model="within")
		reg.ls[[j]]	<- myfit
		tmp2	<- summary(myfit)
		tmp3	<- data.frame(model = names(myfml)[j], DV = dv_vec[i], tmp2$coefficients)
		tmp3$Var<- rownames(tmp3)
		rownames(tmp3) <- NULL
		regrs	<- rbind(regrs, tmp3)
	}
	stargazer(reg.ls, type = "html", title = "Single-equation fixed effects regression", 
			 align = TRUE, no.space = TRUE, 
			 omit = c("year", "month"), order=c("ln_income"), 
			 add.lines = list(c("Year FE","Y","Y","Y","Y"), c("Month FE","Y","Y","Y","Y")), 
			 out = paste(plot.wd, "/", fname, "_se_", dv_vec[i], ".html", sep=""))
			
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
		xlsx.addHeader(mywb, sht1, value = names(tmpls)[i], level = 3)
		xlsx.addTable(mywb, sht1, tmpls[[i]], row.names = F)
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
	
	return(list(coeftab = coeftab, fit = myfit))
}

# Fit SUR with fixed effects 
dv_mat	<- cbind(paste("SHR_", gsub("\\s", "_", fmt_name), sep=""), paste("IC_", gsub("\\s", "_", fmt_name), sep=""))
dv_mat	<- dv_mat[-1,]
colnames(dv_mat) <- c("share", "incidence")
myfml	<- list(Homogeneity		= as.formula("y ~ ln_income + month + year + lag_dol"), 
				Heterogeneity	= as.formula("y ~ ln_income + ln_income_med + ln_income_high + month + lag_dol"), 
				Year 			= as.formula("y ~ ln_income + ln_income_med + ln_income_high + year + month + lag_dol"), 
				Price 			= as.formula("y ~ ln_income + ln_income_med + ln_income_high + lag_dol + year + month + 
													PRC_Convenience_Store + PRC_Discount_Store + PRC_Dollar_Store + 
													PRC_Drug_Store + PRC_Grocery + PRC_Warehouse_Club")
				)
				
regrs1 	<- data.frame()
panid	<- "household_code"
fmln	<- c("Homogeneity", "Low", "Med", "High")
for(i in 1:ncol(dv_mat)){
	for(j in 1:length(myfml)){
		rm(list = "tmp")
		prc 		<- proc.time()
		syseq	<- lapply(dv_mat[,i], function(yn) 
					as.formula(substitute(y ~ x, list(y = as.name(yn), x = terms(myfml[[j]])[[3]])) ) )
		names(syseq)	<- gsub("_", ".", dv_mat[,i])
		tmp		<- try(my.SURFE(syseq, mydata, panid = "household_code", maxiter = 2, residCovWeighted = TRUE, model = FALSE))
		if(class(tmp) != "try-error"){
			# Export the results to html
			tmp1 <- coeftest(tmp$fit)
			tmp2 <- sapply(strsplit(rownames(tmp1), "_"), function(x) x[1])
			tmp3 <- unique(tmp2)
			tmp4 <- lapply(tmp3, function(i) 
					{out <- tmp1[tmp2==i,]; class(out) <- "coeftest"; rownames(out) <- gsub(paste(i, "_", sep=""),"", rownames(out)); return(out)})
			names(tmp4) <- tmp3
			stargazer(tmp4, type = "html", title = "SUR fixed effects regression", dep.var.labels = colnames(dv_mat)[i],
					 align = TRUE, no.space = TRUE, column.labels = gsub("\\.", " ", gsub("SHR.", "",tmp3)), 
					 omit = c("year", "month"), order=c("ln_income"), 
					 out = paste(plot.wd, "/", fname, "_sur_", colnames(dv_mat)[i], "_", names(myfml)[j], ".html", sep=""))
					
			tmp		<- tmp$coeftab
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
	# NOTE: the codes above sometimes don't work. In that case, we run fixed SUR with subset data of income group
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
		xlsx.addHeader(mywb, sht2, value = paste("Model - ", names(tmpls)[i], ":", as.character(myfml[[i]][3]), sep=""), level = 3)
		xlsx.addTable(mywb, sht2, tmpls[[i]], row.names = F)
		# xlsx.addLineBreak(sht2, 1)
	}
}

saveWorkbook(mywb, outxls)

##################
# Plot estimates #
##################
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
	pdf(paste(plot.wd, "/graph_", fname, "_income_shr_avg.pdf", sep=""), width = ww*ar, height = ww)
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
ggtmp1$IncGrp	<- changelab(ggtmp1$Var, paste("ln_income",c("","_med","_high"),sep=""), c("Low", "Med", "High"))
ggtmp1$IncGrp	<- factor(ggtmp1$IncGrp, levels = rev(levels(ggtmp1$IncGrp)), ordered = TRUE)
ggtmp1$Retail 	<- factor(ggtmp1$Retail, levels = ord)

# For median and high income group, the estiamtes are difference relative to low-income group
# So, we add the baseline (the estimate for low-income group) to the difference
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
	pdf(paste(plot.wd, "/graph_", fname, "_income_shr_het.pdf", sep=""), width = ww, height = ww)
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


save.image(file = paste(plot.wd, "/", fname, "_", as.character(Sys.Date()), ".rdata", sep=""))

cat("This program is done.")
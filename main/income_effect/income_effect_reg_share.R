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
fname		<- "expenditure_reg_share"
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

	# Compute cluster-robust standard error 
	# NOTE: clustered standard error has correctted degree of freedom explicity
	cls.cov		<- lapply(myfit$eq, function(x) vcovHC.sur(x, method="arellano", type = "HC1", clustervec = mydata$id))
	# cls.se		<- lapply(myfit$eq, function(x) sqrt(diag(vcovHC.sur(x, method="arellano", type = "HC1", clustervec = mydata$id))))
	cls.se		<- lapply(cls.cov, function(x) sqrt(diag(x)))
	coeftab$cls_se		<- 	unlist(cls.se)
	return(list(coeftab = coeftab, fit = myfit, cls.cov = cls.cov))
}

# Fit SUR with fixed effects 
dv_mat	<- paste("SHR_", gsub("\\s", "_", fmt_name), sep="")
dv_mat	<- dv_mat[-1]
myfml	<- list(Homogeneity		= as.formula("y ~ ln_income + month + lag_dol"), 
				HomoYear		= as.formula("y ~ ln_income + month + year + lag_dol"), 
				HomoPrice		= as.formula("y ~ ln_income + month + year + lag_dol + 
													PRC_Convenience_Store + PRC_Discount_Store + PRC_Dollar_Store + 
													PRC_Drug_Store + PRC_Grocery + PRC_Warehouse_Club"), 
				Heterogeneity 	= as.formula("y ~ ln_income + ln_income_med + ln_income_high + month + lag_dol"), 
				HeterYear 		= as.formula("y ~ ln_income + ln_income_med + ln_income_high + year + month + lag_dol"), 
				HeterPrice 		= as.formula("y ~ ln_income + ln_income_med + ln_income_high + lag_dol + year + month + 
													PRC_Convenience_Store + PRC_Discount_Store + PRC_Dollar_Store + 
													PRC_Drug_Store + PRC_Grocery + PRC_Warehouse_Club")								
				)
				
regrs1 	<- data.frame()
regrs1.v<- regrs1.clsv	<- setNames(vector("list", length = length(myfml)), names(myfml))
panid	<- "household_code"
for(j in 1:length(myfml)){
	rm(list = "tmp")
	prc 		<- proc.time()
	syseq	<- lapply(dv_mat, function(yn) 
				as.formula(substitute(y ~ x, list(y = as.name(yn), x = terms(myfml[[j]])[[3]])) ) )
	names(syseq)	<- gsub("_", ".", dv_mat)
	tmp		<- try(my.SURFE(syseq, mydata, panid = "household_code", maxiter = 2, residCovWeighted = TRUE))
	if(class(tmp) != "try-error"){	
		# Export the results to html
		tmp1 <- coeftest(tmp$fit)
		tmp2 <- sapply(strsplit(rownames(tmp1), "_"), function(x) x[1])
		tmp3 <- unique(tmp2)
		tmp4 <- lapply(tmp3, function(i) 
				{out <- tmp1[tmp2==i,]; class(out) <- "coeftest"; rownames(out) <- gsub(paste(i, "_", sep=""),"", rownames(out)); return(out)})
		names(tmp4) <- tmp3
		cls.se	<- lapply(tmp3, function(i) tmp$coeftab$cls_se[tmp2==i])
		stargazer(tmp4, type = "html", title = "SUR fixed effects regression", dep.var.labels = colnames(dv_mat)[i],
				 align = TRUE, no.space = TRUE, column.labels = gsub("\\.", " ", gsub("SHR.", "",tmp3)), 
				 omit = c("year", "month"), order=c("ln_income"), se = cls.se, 
			 	 notes = c("se clustered over households"), 
				 out = paste(plot.wd, "/", fname, "_sur_share", "_", names(myfml)[j], ".html", sep=""))

		stargazer(tmp4, type = "latex", title = "SUR fixed effects regression", dep.var.labels = colnames(dv_mat)[i],
				 align = TRUE, no.space = TRUE, column.labels = gsub("\\.", " ", gsub("SHR.", "",tmp3)), 
				 omit = c("year", "month"), order=c("ln_income"), digits = 4, se = cls.se, 
				 out = paste(plot.wd, "/", fname, "_sur_share", "_", names(myfml)[j], ".tex", sep=""))				
		
		# Save the coefficient estimates 
		regrs1.v[[j]]	<- vcov(tmp$fit)								# Covariance of se
		regrs1.clsv[[j]]<- setNames(tmp$cls.cov, dv_mat)				# Covariance of clustered se.
		tmp		<- tmp$coeftab
		tmp1 	<- sapply(rownames(tmp), function(x) strsplit(x, "_")[[1]][1])
		tmp2	<- sapply(1:length(tmp1), function(i) gsub(paste(tmp1[i], "_", sep=""), "", rownames(tmp)[i]))
		tmp3	<- cbind(model = names(myfml)[j], DV = tmp1, Var = tmp2, tmp)
		rownames(tmp3) <- NULL
		regrs1	<- rbind(regrs1, tmp3)
		use.time 	<- proc.time() - prc
		cat("Regressions for expenditure share with", names(myfml)[j],
			"responses finishes, using", use.time[3]/60, "min.\n")
	}else{
		cat("Regressions for expenditure share with", names(myfml)[j],
			"responses fails.\n")
	}
}

# Organize the parameter table 
tmp		<- regrs1[grep("ln_income", regrs1$Var),]
tmp1	<- setNames(fmt_name, paste("SHR.", gsub("\\s", "\\.", fmt_name), sep=""))
tmp$DV1	<- tmp1[as.character(tmp$DV)]
tmp1	<- setNames(c("log(I)", "log(I)*Med", "log(I)*High"),c("ln_income", "ln_income_med", "ln_income_high"))
tmp$Var1<- tmp1[as.character(tmp$Var)]

tmp1	<- split(tmp, tmp$model)
for(i in 1:length(tmp1)){
	tmp			<- cbind(tmp1[[i]]$Estimate, tmp1[[i]]$cls_se, tmp1[[i]]$Estimate/tmp1[[i]]$cls_se)
	tmp			<- cbind(tmp, pt(-abs(tmp[,3]), tmp1[[i]]$adj_df)*2)
	dimnames(tmp)	<- list(paste(tmp1[[i]]$DV1, tmp1[[i]]$Var1, sep="*"), c("Estimate", "Std. Error", "t value", "Pr(>|t|)")) 
	class(tmp)	<- "coeftest"
	tmp1[[i]]	<- tmp
}

stargazer(tmp1, type = "html", title = "Summary of SUR regressions", dep.var.labels = "Expenditure share",
		 align = TRUE, no.space = TRUE, column.labels = names(tmp1), 
	 	 notes = c("se clustered over households"), 
		add.lines = list(c("Year FE","N","Y","Y","N","Y","Y"), c("Month FE","Y","Y","Y","Y","Y","Y"), 
						 c("HH FE", "Y","Y","Y","Y","Y","Y"), c("Price", "N","N","Y","N","N","Y")), 
		 out = paste(plot.wd, "/", fname, "_sur_share_sum.html", sep=""))

stargazer(tmp1, type = "latex", title = "Summary of SUR regressions", dep.var.labels = "Expenditure share",
		 align = TRUE, no.space = TRUE, column.labels = names(tmp1), 
	 	 notes = c("se clustered over households"), 
		add.lines = list(c("Year FE","N","Y","Y","N","Y","Y"), c("Month FE","Y","Y","Y","Y","Y","Y"), 
						 c("HH FE", "Y","Y","Y","Y","Y","Y"), c("Price", "N","N","Y","N","N","Y")), 
		 out = paste(plot.wd, "/", fname, "_sur_share_sum.tex", sep=""))

##################
# Plot estimates #
##################
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

test.alpha	<- .9						# Significance level
crtv		<- qnorm(test.alpha)		# Critical value
sel.model	<- c("Homogeneity", "Heterogeneity")		# Select model in the plot
ret.ord		<- c("Dollar Store", "Grocery", "Discount Store", "Drug Store", "Warehouse Club")
dodge.width	<- .5						# The width parameter in position_dodge

# Homogenous income effect for expenditure in SUR regressions
ggtmp1	<- subset(regrs1, model == sel.model[1] & substr(DV,1,3) == "SHR" & Var == "ln_income")
ggtmp1$Retail	<- factor(as.character(ggtmp1$DV), levels = paste("SHR.", gsub("\\s", ".", fmt_name[-1]), sep=""), 
						labels = fmt_name[-1])
# ord		<- order(ggtmp1$Estimate)
# ret.ord	<- ggtmp1[ord,"Retail"]
ggtmp1$Retail 	<- factor(ggtmp1$Retail, levels = ret.ord)

# Heterogenous income effect for expenditure in SUR regressions
ggtmp2	<- subset(regrs1, model == sel.model[2] & substr(DV,1,3) == "SHR" & substr(Var,1,9) == "ln_income")
ggtmp2$Retail	<- factor(as.character(ggtmp2$DV), levels = paste("SHR.", gsub("\\s", ".", fmt_name[-1]), sep=""), 
						labels = fmt_name[-1])						
ggtmp2$IncGrp	<- changelab(ggtmp2$Var, paste("ln_income",c("","_med","_high"),sep=""), c("Low", "Med", "High"))
ggtmp2$IncGrp	<- factor(ggtmp2$IncGrp, levels = rev(levels(ggtmp2$IncGrp)), ordered = TRUE)
ggtmp2$Retail 	<- factor(ggtmp2$Retail, levels = ret.ord)

# For median and high income group, the estiamtes are difference relative to low-income group
# So, we add the baseline (the estimate for low-income group) to the difference
sel		<- ggtmp2$Var == "ln_income"
tmp		<- ggtmp2[sel,"Estimate"]
names(tmp)	<- ggtmp2[sel,"DV"]
tmp1	<- tmp[as.character(ggtmp2$DV)]
tmp1[sel]	<- 0
ggtmp2$Estimate1	<- ggtmp2$Estimate + tmp1

# Compute standard error of income effect for each group
# For median income group, var(b) = sum(cov[ln_income, ln_income_med])
sel		<- setNames(list(c("ln_income"), c("ln_income", "ln_income_med"), c("ln_income", "ln_income_high")), 
					c("ln_income", "ln_income_med","ln_income_high"))
tmp		<- sapply(regrs1.clsv[["Heterogeneity"]], function(x) sapply(sel, function(i) sqrt(sum(x[i,i]))))				
tmp		<- melt(tmp, value.name = c("Cluster_se"))
tmp$Var2<- gsub("\\_", "\\.", tmp$Var2)
names(tmp)[1:2]	<- c("Var", "DV")
ggtmp2	<- merge(ggtmp2, tmp, by = c("Var", "DV"), all.x = T)

tmp1	<- gsub("\\_", "\\.", dv_mat)
tmp		<- sapply(tmp1, function(i) sapply(sel, function(j) sum(regrs1.v[["Heterogeneity"]][paste(i,"_", j,sep=""), paste(i,"_", j,sep="")])))
tmp		<- melt(tmp, value.name = "StdErr")
names(tmp)[1:2]	<- c("Var", "DV")
ggtmp2	<- merge(ggtmp2, tmp, by = c("Var", "DV"), all.x = T)

# Find the range of estimates
mylim		<- range(c(with(ggtmp2, Estimate1-crtv*Cluster_se), with(ggtmp2, Estimate1 + crtv * Cluster_se)))
cat("y axis range is", mylim, "\n")					

plots	<- list(NULL)
plots[[1]]	<- ggplot(ggtmp1, aes(Retail, Estimate)) + 
			geom_pointrange(aes(y = Estimate, ymin = Estimate-crtv*cls_se, ymax = Estimate + crtv * cls_se)) + 
			ylim(mylim) + 
			xlab("Retail format") + coord_flip() + theme_bw()
plots[[2]]	<- ggplot(ggtmp2, aes(Retail, Estimate1, col = IncGrp)) + 
			geom_pointrange(aes(ymin = Estimate1-crtv*Cluster_se, ymax = Estimate1 + crtv * Cluster_se, col = IncGrp), 
					position=position_dodge(width=dodge.width)) +
			scale_color_grey(name="Income\ngroup", start = 0, end = .8) +  
			ylim(mylim) + 
			xlab("Retail format") + ylab("Estimate")+ coord_flip() +
			guides(fill = guide_legend(reverse = TRUE), color = guide_legend(reverse = TRUE)) + 
			theme_bw()
plots[[3]]	<- get_legend(plots[[2]])
plots[[2]]	<- plots[[2]] + theme(legend.position="none")	
					
if(make_plot){
	pdf(paste(plot.wd, "/graph_", fname, "_income_shr_est.pdf", sep=""), width = ww1, height = ww1*ar)
	grid.arrange(plots[[1]], plots[[2]], plots[[3]], nrow = 1, widths = c(.45, .45, .1))
	dev.off()
}

# Convert the estimates to the impacts of 10% income drop 
myscale	<- -.1
mylim		<- range(c(with(ggtmp2, Estimate1*myscale-crtv*Cluster_se*abs(myscale)), 
					   with(ggtmp2, Estimate1*myscale + crtv * Cluster_se*abs(myscale))))
plots	<- list(NULL)
plots[[1]]	<- ggplot(ggtmp1, aes(Retail, Estimate*myscale)) + 
			geom_pointrange(aes(y = Estimate*myscale, ymin = Estimate*myscale-crtv*cls_se*abs(myscale), 
							ymax = Estimate*myscale + crtv*cls_se*abs(myscale))) + 
			geom_hline(yintercept = 0, size = .25, linetype = 2) + 				
			scale_y_continuous(labels = scales::percent, limits = mylim) + 
			xlab("Retail format") + ylab("")+ coord_flip() + theme_bw()
plots[[2]]	<- ggplot(ggtmp2, aes(Retail, Estimate1*myscale, col = IncGrp)) + 
			geom_pointrange(aes(ymin = Estimate1*myscale-crtv*Cluster_se*abs(myscale), 
								ymax = Estimate1*myscale + crtv * Cluster_se*abs(myscale), col = IncGrp), 
					position=position_dodge(width=dodge.width)) +
			geom_hline(yintercept = 0, size = .25, linetype = 2) + 		
			scale_color_grey(name="Income\ngroup", start = 0, end = .8) +  
			scale_y_continuous(labels = scales::percent, limits = mylim) +
			xlab("Retail format") + ylab("")+ coord_flip() + 
			guides(color = guide_legend(reverse = TRUE)) + theme_bw()
plots[[3]]	<- get_legend(plots[[2]])
plots[[2]]	<- plots[[2]] + theme(legend.position="none")	

if(make_plot){
	pdf(paste(plot.wd, "/graph_", fname, "_income_shr.pdf", sep=""), width = ww1, height = ww1*ar)
	grid.arrange(plots[[1]], plots[[2]], plots[[3]], nrow = 1, widths = c(.45, .45, .1), 
		# top = "Change of expenditure share for a 10% income drop", bottom = "Change in expenditure share")
		bottom = "Change in expenditure share")
	dev.off()
}

# -------------------------------------- # 
# Visualize the different standard error # 
# Homogenous income effect for expenditure in SUR regressions
ggtmp1	<- subset(regrs1,substr(model,1,4) == "Homo" & substr(DV,1,3) == "SHR" & Var == "ln_income")
ggtmp1$Retail	<- factor(as.character(ggtmp1$DV), levels = paste("SHR.", gsub("\\s", ".", fmt_name[-1]), sep=""), 
						labels = fmt_name[-1])
ggtmp1$Retail 	<- factor(ggtmp1$Retail, levels = ret.ord)

# Heterogenous income effect for expenditure in SUR regressions
ggtmp2	<- subset(regrs1, substr(model,1,5) == "Heter" & substr(DV,1,3) == "SHR" & substr(Var,1,9) == "ln_income")
ggtmp2$Retail	<- factor(as.character(ggtmp2$DV), levels = paste("SHR.", gsub("\\s", ".", fmt_name[-1]), sep=""), 
						labels = fmt_name[-1])						
ggtmp2$IncGrp	<- changelab(ggtmp2$Var, paste("ln_income",c("","_med","_high"),sep=""), c("Low", "Med", "High"))
ggtmp2$IncGrp	<- factor(ggtmp2$IncGrp, levels = rev(levels(ggtmp2$IncGrp)), ordered = TRUE)
ggtmp2$Retail 	<- factor(ggtmp2$Retail, levels = ret.ord)

# For median and high income group, the estiamtes are difference relative to low-income group
# So, we add the baseline (the estimate for low-income group) to the difference
sel		<- ggtmp2$Var == "ln_income"
tmp		<- ggtmp2[sel,"Estimate"]
names(tmp)	<- ggtmp2[sel,"DV"]
tmp1	<- tmp[as.character(ggtmp2$DV)]
tmp1[sel]	<- 0
ggtmp2$Estimate1	<- ggtmp2$Estimate + tmp1

# For median income group, var(b) = sum(cov[ln_income, ln_income_med])
tmp1	<- tmp2	<- data.frame()
sel		<- setNames(list(c("ln_income"), c("ln_income", "ln_income_med"), c("ln_income", "ln_income_high")), 
					c("ln_income", "ln_income_med","ln_income_high"))
for(i in unique(ggtmp2$model)){
	# Clustered standard error
	tmp		<- sapply(regrs1.clsv[[i]], function(x) sapply(sel, function(j) sqrt(sum(x[j,j]))))				
	tmp		<- melt(tmp, value.name = c("Cluster_se"))
	tmp$Var2<- gsub("\\_", "\\.", tmp$Var2)
	names(tmp)[1:2]	<- c("Var", "DV")
	tmp$model	<- i
	tmp1	<- rbind(tmp1, tmp)
	
	# Standard error
	sel1	<- gsub("\\_", "\\.", dv_mat)
	tmp		<- sapply(sel1, function(k) sapply(sel, function(j) sum(regrs1.v[[i]][paste(k,"_", j,sep=""), paste(k,"_", j,sep="")])))
	tmp		<- melt(tmp, value.name = "StdErr")
	tmp$model	<- i
	names(tmp)[1:2]	<- c("Var", "DV")
	tmp2	<- rbind(tmp2, tmp)	
}
ggtmp2	<- merge(ggtmp2, tmp1, by = c("model", "Var", "DV"), all.x = T)
ggtmp2	<- merge(ggtmp2, tmp2, by = c("model", "Var", "DV"), all.x = T)

plots	<- list(NULL)				
plots[[1]]	<- ggplot(ggtmp1, aes(Retail, Estimate)) + 
			geom_pointrange(aes(y = Estimate, ymin = Estimate-crtv*Std..Error, ymax = Estimate + crtv * Std..Error), size = .5) + 
			geom_errorbar(aes(ymin = Estimate-crtv*cls_se, ymax = Estimate + crtv * cls_se), size = .25, width=0.25) + 
			facet_wrap(~model) + 
			xlab("Dependent variable") + coord_flip()
plots[[2]]	<- ggplot(ggtmp2, aes(Retail, Estimate1, fill = IncGrp, col = IncGrp)) + 
			geom_pointrange(aes(y = Estimate1, ymin = Estimate1-crtv*StdErr, ymax = Estimate1 + crtv *StdErr), size = .5, position = position_dodge(width=0.5)) + 
			geom_errorbar(aes(ymin = Estimate1-crtv*Cluster_se, ymax = Estimate1 + crtv * Cluster_se), size = .25, width = .25, position = position_dodge(width=0.5)) + 
			scale_fill_grey() + 
			scale_color_grey() +  
			facet_wrap(~model) +
			xlab("Retail format") + ylab("Estimate")+ coord_flip() + 
			guides(fill = guide_legend(reverse = TRUE), color = guide_legend(reverse = TRUE))


if(make_plot){
	pdf(paste(plot.wd, "/graph_", fname, "_income_shr_compare.pdf", sep=""), width = ww1, height = ww1*.8)
	for(i in 1:length(plots)){ print(plots[[i]]) }
	dev.off()
}		

# Save results 
rm(list = intersect(ls(), c("tmp", "tmp1", "tmp2", "tmp3", "tmp4", "use.time", "make_plot", "prc","tmp.coef", "tmp_dv", 
							"i", "j", "sel", "sel1", "myfml", "tmpidx")))

save.image(file = paste(plot.wd, "/", fname, "_", as.character(Sys.Date()), ".rdata", sep=""))

cat("This program is done.")
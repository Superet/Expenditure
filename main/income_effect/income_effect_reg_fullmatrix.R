# This file estimate the income effects in full matrix. 
# beta_ij: the changes if a household transit from income level i to income level j. 

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
ar 			<- .8	#.6

week.price	<- FALSE
cpi.adj		<- TRUE
write2csv	<- FALSE
make_plot	<- TRUE
fname		<- "expenditure_reg_fmatrix"
if(cpi.adj) { fname <- paste(fname, "_cpi", sep="")}
if(week.price){ fname <- paste(fname, "_wkprc", sep="")}

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

test.alpha	<- .9						# Significance level
crtv		<- qnorm(test.alpha)		# Critical value
ret.ord		<- c("Dollar Store", "Grocery", "Discount Store", "Drug Store", "Warehouse Club")
dodge.width	<- .5						# The width parameter in position_dodge

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

time_diff <- function(x1, t1, x2, t2){
# This function returns a vector with the length of x2, with elements being the differnce from x1 at the same period
# All arguments need to be sorted by time
	out	<- rep(NA, length(x2))
	t.in<- intersect(t1, t2)
	sel1<- t1 %in% t.in
	sel2<- t2 %in% t.in
	out[sel2]	<- x2[sel2] - x1[sel1]
	return(out)
}

time_diff_init	<- function(x, t, blk, fill.miss = TRUE){
	out			<- rep(NA, length(x))
	blk.unq		<- sort(unique(blk))
	sel			<- blk == blk.unq[1]
	x.init		<- x[sel]
	t.init		<- t[sel]
	i			<- 2
	while(i <= length(blk.unq)){
		sel2	<- blk == blk.unq[i]
		sel.mis	<- !t[sel2] %in% t.init
		a		<- time_diff(x.init, t.init, x[sel2], t[sel2])
		out[sel2][!is.na(a)]	<- a[!is.na(a)]
		if(any(sel.mis) & fill.miss){
			t.init	<- c(t.init, t[sel2][sel.mis])
			x.init	<- c(x.init, x[sel2][sel.mis])
			ord		<- order(t.init)
			t.init	<- t.init[ord]
			x.init	<- x.init[ord]
		}
		i		<- i + 1
	}
	return(out)
}

time_diff_period	<- function(x, t, blk){
	out			<- rep(NA, length(x))
	blk.unq		<- sort(unique(blk))
	sel			<- blk == blk.unq[1]
	i			<- 1
	while(i < length(blk.unq)){
		sel1	<- blk == blk.unq[i]
		sel2	<- blk == blk.unq[(i+1)]
		a		<- time_diff(x[sel1], t[sel1], x[sel2], t[sel2])
		out[sel2][!is.na(a)]	<- a[!is.na(a)]
		i		<- i + 1
	}
	return(out)
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
pan_yr    	<- pan_yr[,':='(year_id= year - year[1] + 1, tenure = length(year), move = c(0, 1*(diff(Inc)!=0)), 
							sgn_init = c(NA,sign(Inc[-1]-Inc[1])), sgn_period = c(NA, sign(diff(Inc))))
                          	# move = 1*(!all(Inc==Inc[1]))) 
                    , by = list(household_code)]
pan_yr		<- pan_yr[,move_all := sum(move), by = list(household_code)]
pan_yr		<- pan_yr[, ':='(last_inc = c(NA, Inc[-length(Inc)]), init_inc = Inc[1]), by = list(household_code)]
pan_yr$Inc	<- factor(pan_yr$Inc)
pan_yr$last_inc	<- factor(pan_yr$last_inc)
pan_yr$init_inc	<- factor(pan_yr$init_inc)
panelist	<- pan_yr[, list(tenure = unique(tenure), p.change = unique(move_all/tenure), 
							large.init.change = max(abs(Income-Income[1])/Income[1] )), 
						by = list(first_incomeg, household_code)]
summary(panelist)
table(pan_yr$init_inc, pan_yr$Inc)
table(pan_yr$last_inc, pan_yr$Inc)

# --------------- #
# Regression data #
# Drop households who stay only for one year
mydata 	<- subset(hh_exp, household_code %in% panelist[tenure>1, household_code])
# mydata	<- hh_exp
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
# mydata$year			<- as.factor(mydata$year)

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

# Caculate biweek within a year
tmp			<- as.vector(by(mydata$biweek, mydata$year, max))
tmp			<- setNames(c(0, tmp[-length(tmp)]), 2004:2010)
mydata$biweek.wyr	<- mydata$biweek - tmp[as.character(mydata$year)]

# Data of differences relative to the first year
selcol		<- c("ln_dol", paste("SHR_", gsub("\\s", "_", fmt_name), sep=""))
delta.init	<- data.table(mydata[, c("household_code", "biweek", "biweek.wyr", 'year',selcol)])
setkeyv(delta.init, c("household_code", "biweek"))
for(v in selcol){
	delta.init	<- delta.init[, eval(as.name(v)):=time_diff_init(eval(as.name(v)), biweek.wyr, year), by = list(household_code)]
}
tmp		<- pan_yr[, list(first_incomeg, household_code, year, init_inc, Inc, year_id)]
setnames(tmp, "init_inc", "base_inc")
dim(delta.init)
delta.init	<- merge(delta.init, tmp[year_id>1,], 
					by = c("household_code", "year"))
dim(delta.init)
sum(is.na(delta.init))		
delta.init	<- delta.init[!is.na(delta.init$ln_dol) & abs(delta.init$ln_dol) != Inf,]

# Data of different from period to period
delta.period	<- data.table(mydata[, c("household_code", "biweek", "biweek.wyr", 'year',selcol)])
setkeyv(delta.period, c("household_code", "biweek"))
for(v in selcol){
	delta.period	<- delta.period[, eval(as.name(v)):=time_diff_period(eval(as.name(v)), biweek.wyr, year), by = list(household_code)]
}
tmp		<- pan_yr[, list(first_incomeg, household_code, year, last_inc, Inc, year_id)]
setnames(tmp, "last_inc", "base_inc")
dim(delta.period)
delta.period	<- merge(delta.period, tmp[year_id>1,], 
					by = c("household_code", "year"))
dim(delta.period)
sum(is.na(delta.period))		
delta.period	<- delta.period[!is.na(delta.period$ln_dol) & abs(delta.period$ln_dol) != Inf,]

############################################
# Single-equation fixed effect regressions #
############################################
#---------------------------------------#
# Set up DV and IV, regression models 
dv_vec	<- "ln_dol"
dat.vec	<- c("delta.init", "delta.period")
myfml	<- y ~ base_inc+Inc - 1

# Run regressions
regrs	<- data.frame()
prc 	<- proc.time()
reg.ls	<- setNames(vector("list", length(dat.vec)), dat.vec)

for(i in 1:length(dat.vec)){
	rm(list = intersect(ls(), "myfit"))
	tmp 	<- as.formula(substitute(y ~ x, list(y = as.name(dv_vec), x = terms(myfml)[[3]])) )
	myfit	<- plm(tmp, data = eval(parse(text = dat.vec[i])), index = c("household_code","biweek"), model="pooling")
	reg.ls[[i]]	<- myfit
	tmp2	<- summary(myfit)
	tmp3	<- data.frame(model = dat.vec[i], DV = dv_vec, tmp2$coefficients)
	tmp3$Var<- rownames(tmp3)
	rownames(tmp3) <- NULL
	regrs	<- rbind(regrs, tmp3)
}					

# Compute cluster-robust se
tmp.coef	<- lapply(reg.ls, coeftest)
unc.se		<- lapply(tmp.coef, function(x) x[, "Std. Error"])
cls.v1		<- lapply(1:length(dat.vec), function(i) vcovHC.plm.new(reg.ls[[i]], method = "arellano", type = "HC1", cluster = "define", 
						clustervec = get("household_code", get(dat.vec[i])) ))
cls.se1		<- lapply(1:length(dat.vec), function(i) sqrt(diag(cls.v1[[i]] )))
# cls.se2		<- lapply(reg.ls, function(x) sqrt(diag(vcovHC.plm.new(x, method = "arellano", type = "HC1", cluster = "define", clustervec = mydata[sel,"hh_year"]))))

# Put together the se 
tmp			<- lapply(1:length(tmpidx), function(i) cbind(Uncorrected = unc.se[[i]], 
														HH.cluster.se = cls.se1[[i]] ))
names(tmp)	<- names(reg.ls)
cat("Compare cluster robust standard error for regressions of log(expenditure) ~ log(income):\n"); print(tmp); cat("\n")

# Plot income elasticity
ggtmp		<- regrs
ggtmp$cls_se	<- unlist(cls.se1)
ggtmp$value		<- as.numeric(gsub("\\D", "", ggtmp$Var)) 
sel			<- grep("base_inc", ggtmp$Var)
tmp1		<- split(ggtmp[sel,], ggtmp[sel,"model"])
tmp2		<- split(ggtmp[-sel,], ggtmp[-sel,"model"])
out			<- list(NULL)
ggtmp		<- data.frame()

for(i in 1:length(reg.ls)){
	tmp		<- setNames(codebook$description, codebook$code_value)
	out1	<- setNames(vector("list", nrow(tmp2[[1]])), tmp[as.character(tmp2[[1]]$value)])
	for(j in 1:length(out1)){
		a	<- tmp1[[i]]
		b	<- tmp2[[i]][j,]
		sel.a	<- grep("base", rownames(cls.v1[[i]]))
		sel.b	<- grep(b$Var, colnames(cls.v1[[i]]))
		tmpv	<- diag(cls.v1[[i]][sel.a, sel.a]) + cls.v1[[i]][sel.b,sel.b] + 2* cls.v1[[i]][sel.a, sel.b] 
		tmp			<- cbind(a$Estimate + b$Estimate, sqrt(tmpv))
		tmp			<- cbind(tmp, tmp[,1]/tmp[,2], pt(-abs(tmp[,1]/tmp[,2]), reg.ls[[i]]$df.residual)*2)
		colnames(tmp)	<- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
		ggtmp1		<- data.frame(model = dat.vec[i], tmp, Base = rownames(tmp), Change = names(out1)[j])
		rownames(ggtmp1)	<- NULL
		ggtmp		<- rbind(ggtmp, ggtmp1)
		class(tmp)	<- "coeftest"
		out1[[j]]	<- tmp
	}
	out[[i]]		<- out1
}

# Export to HTML
tmp		<- setNames(as.character(codebook$description), codebook$code_value)
tmplab	<- tmp[as.character(tmp1[[1]]$value)]
for(i in 1:length(out)){
	stargazer(out[[i]], type = "html", title = "Single-equation fixed effects regression", 
		     column.labels = names(out1),
			 align = TRUE, no.space = TRUE, 
			 covariate.labels = tmplab, 
			 notes = c("se clustered over households"),
			 out = paste(plot.wd, "/", fname, "_se_", dv_vec, "_", dat.vec[i], ".html", sep=""))
}

# Plot bubble plot of elasticity
ggtmp$bvalue	<- as.numeric(gsub("\\D", "", ggtmp$Base)) 
ggtmp$Base1		<- tmp[as.character(ggtmp$bvalue)]
ggtmp$sgn_est	<- factor(sign(ggtmp$Estimate))

pdf(paste(plot.wd, "/", fname, "_elasticity.pdf", sep=""), width = ww1, height = ww1*ar)
for(i in 1:length(dat.vec)){
	print(ggplot(subset(ggtmp, model == dat.vec[i]), aes(Base1, Change, size = abs(Estimate), col = sgn_est)) + 
			geom_point(shape = 21) + 
			scale_size_continuous(range = c(1, 10)) + 
			theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
			labs(title = dat.vec[i]) + 
			theme_bw()
		)
}
dev.off()

##############################
# Fit SUR with fixed effects #
##############################
dv_mat	<- paste("SHR_", gsub("\\s", "_", fmt_name), sep="")
dv_mat	<- dv_mat[-1]
				
regrs1 	<- data.frame()
regrs1.v<- regrs1.clsv	<- setNames(vector("list", length = length(dat.vec)), dat.vec)
panid	<- "household_code"

for(i in 1:length(dat.vec)){
	rm(list = "tmp")
	prc 		<- proc.time()
	syseq	<- lapply(dv_mat, function(yn) 
				as.formula(substitute(y ~ x, list(y = as.name(yn), x = terms(myfml)[[3]])) ) )
	names(syseq)	<- gsub("_", ".", dv_mat)
	
	myfit 	<- systemfit(syseq, data = eval(parse(text = dat.vec[i])), method = "SUR", maxiter = 2, residCovWeighted = TRUE)
	print(summary(myfit))
	coeftab <- data.frame(summary(myfit)$coefficients)
	coeftab$df			<- myfit$df.residual/length(syseq)
	
	# Compute cluster-robust standard error 
	# NOTE: clustered standard error has correctted degree of freedom explicity
	cls.cov		<- lapply(myfit$eq, function(x) vcovHC.sur(x, method="arellano", type = "HC1", 
						clustervec = get("household_code", eval(parse(text = dat.vec[i]))) ))
	cls.se		<- lapply(cls.cov, function(x) sqrt(diag(x)))
	coeftab$cls_se		<- 	unlist(cls.se)
		
	# Save the coefficient estimates 
	regrs1.v[[i]]	<- vcov(myfit)								# Covariance of se
	regrs1.clsv[[i]]<- setNames(cls.cov, dv_mat)				# Covariance of clustered se.
	tmp1 	<- sapply(rownames(coeftab), function(x) strsplit(x, "_")[[1]][1])
	tmp2	<- sapply(1:length(tmp1), function(i) gsub(paste(tmp1[i], "_", sep=""), "", rownames(coeftab)[i]))
	tmp3	<- cbind(model = dat.vec[i], DV = tmp1, Var = tmp2, coeftab)
	rownames(tmp3) <- NULL
	regrs1	<- rbind(regrs1, tmp3)
	use.time 	<- proc.time() - prc
	cat("Regressions for expenditure share with", dat.vec[i],
		"responses finishes, using", use.time[3]/60, "min.\n")
}

# Plot the output
ggtmp		<- regrs1
ggtmp$value		<- as.numeric(gsub("\\D", "", ggtmp$Var)) 
sel			<- grep("base_inc", ggtmp$Var)
tmp1		<- split(ggtmp[sel,], ggtmp[sel,"model"])
tmp2		<- split(ggtmp[-sel,], ggtmp[-sel,"model"])
tmp0		<- unique(regrs1$DV)
ggtmp		<- data.frame()

for(i in 1:length(dat.vec)){
	for(j in 1:length(unique(tmp2[[1]]$Var))){
		for(k in 1:length(dv_mat)){
			a	<- tmp1[[i]][tmp1[[i]]$DV == tmp0[k], ]
			b	<- tmp2[[i]][tmp2[[i]]$DV == tmp0[k],][j,]
			sel.a	<- grep("base", rownames(regrs1.clsv[[i]][[k]]))
			sel.b	<- grep(b$Var, colnames(regrs1.clsv[[i]][[k]]))
			tmpv	<- diag(regrs1.clsv[[i]][[k]][sel.a, sel.a]) + regrs1.clsv[[i]][[k]][sel.b,sel.b] + 2* regrs1.clsv[[i]][[k]][sel.a, sel.b] 
			tmp			<- cbind(a$Estimate + b$Estimate, sqrt(tmpv))
			tmp			<- cbind(tmp, tmp[,1]/tmp[,2], pt(-abs(tmp[,1]/tmp[,2]), reg.ls[[i]]$df.residual)*2)
			colnames(tmp)	<- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
			ggtmp1		<- data.frame(model = dat.vec[i], Retail = dv_mat[k], tmp, Base = a$value, Change = b[,"value"])
			rownames(ggtmp1)	<- NULL
			ggtmp		<- rbind(ggtmp, ggtmp1)
		}
	}
}

tmp		<- setNames(as.character(codebook$description), codebook$code_value)
ggtmp$Base1		<- factor(tmp[as.character(ggtmp$Base)], levels = tmp)
ggtmp$Change1	<- factor(tmp[as.character(ggtmp$Change)], levels = tmp)
ggtmp$sgn_est	<- factor(sign(ggtmp$Estimate))

pdf(paste(plot.wd, "/graph_", fname, "_sur_share.pdf", sep=""), width = ww1, height = ww1*ar)
for(i in 1:length(dat.vec)){
	for(j in 1:length(dv_mat)){
		print(ggplot(subset(ggtmp, model == dat.vec[i] & Retail == dv_mat[j]), aes(Base1, Change1, size = abs(Estimate), col = sgn_est)) + 
				geom_point() + 
				scale_size_area(breaks=c(0, .005, .010, .015, .020, .025, .030, .035, .040), max_size=5)+ 
				theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
				theme_bw() + 
				labs(title = paste("Change of expenditure share (max=1) at\n", fmt_name[(j+1)], " in data ", dat.vec[i],sep=""))	
			)
	}
}
dev.off()

save.image(file = paste(plot.wd, "/",fname, ".rdata", sep=""))


library(ggplot2)
library(reshape2)
library(data.table)
library(plm)
library(gridExtra)
library(scales)
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
tmpcsv		<- paste(plot.wd,"/4_nielsen_regtab_larger.csv", sep="")
if(write2csv){
	f		<- file(tmpcsv, "w")
	writeLines("Output from 1_nielsen_4_sumstat.r\n", f)
	close(f)
}

# Read data 
hh_exp		<- read.csv("2_hh_biweek_exp_merge_larger.csv")
fmt_attr 	<- read.csv("1_format_year_attr.csv")
price_dat	<- read.csv("1_format_biweek_price.csv")
source("outreg function.R")

# load("4_sumstat.rdata")

#############
# Functions # 
#############
within_loess <- function(dv, groupv, data){
	myfml	<- as.formula(paste(dv, "~ factor(month)", sep=""))
	ggtmp 	<- data.frame()
	groupl	<- sort(unique(data[,groupv]))
	for(i in 1:length(groupl)){
		sel		<- data[,groupv] == groupl[i]
		myfit 	<- plm(myfml, data=data[sel,], index = c("household_code","biweek"), model="within")
		res		<- myfit$residuals
		sel1	<- !is.na(data[,dv]) & sel
		sel2	<- data$biweek == 1 & sel
		ggtmp	<- rbind(ggtmp, data.frame(biweek = data[sel1,"biweek"], res = res, mean = mean(hh_exp[sel2,dv], na.rm=T), 
							  Segment = groupl[i]))
	}
	ggtmp$DV	<- dv
	ggtmp$y		<- ggtmp$mean + ggtmp$res
	return(ggtmp)
}

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

############################
# Construct data variables # 
############################
fmt_name 	<- as.character(sort(unique(fmt_attr$channel_type)))
R			<- length(fmt_name)
sel			<- paste("DOLP_", gsub("\\s", "_", fmt_name), sep="")
hh_exp$y 	<- rowSums(hh_exp[,sel])
hh_exp$income_real <- factor(hh_exp$income_real)
sel1		<- gsub("DOLP", "SHR", sel)
for(i in 1:length(fmt_name)){
	hh_exp[,sel1[i]] <- hh_exp[,sel[i]]/hh_exp$y
}
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

#####################
# Format attributes # 
#####################
# Stat summary of other attributes across market and years
# selcol	<- c("size_index", "num_module", "num_brands", "avg_brands_per_mod", "num_upc", "avg_upc_per_mod", "avg_prvt_per_mod", "overall_prvt")
# tmp		<- c("Size index","Total num module","Total num brands","Num brands per module", "Num UPC","Num UPC per module",
# 		"Proportion private label per module", "Total proportion private label")
selcol 	<- c("size_index", "overall_prvt", "num_module", "avg_upc_per_mod")
tmp		<- c("Size index", "Total proportion private label","Total num module","Num UPC per module")
ggtmp <- melt(fmt_attr[,c("scantrack_market_descr", "year", "channel_type", selcol)], 
				id.var=c("scantrack_market_descr", "year", "channel_type"))
ggtmp$variable <- factor(ggtmp$variable, levels= selcol, labels= tmp)

if(make_plot){
	pdf(paste(plot.wd,"/graph_retail_distribution.pdf", sep=""), width = ww, height = ww)
	print( ggplot(ggtmp, aes(channel_type, value)) + 
		stat_summary(fun.y=mean, fun.ymin=function(x) quantile(x, .25), fun.ymax = function(x) quantile(x, .75), geom ="pointrange", size=.5) + 
		stat_summary(fun.y=mean, fun.ymin=function(x) quantile(x, .025), fun.ymax = function(x) quantile(x, .975), geom ="pointrange", size=.25) + 
		facet_wrap(~variable, scales="free_y") + 
		theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
		xlab("Channel type")
	)
	dev.off()
}

# Variance decomposition of format attributes
tmp.tab		<- data.frame()
for(i in 1:length(selcol)){
	myfml	<- paste(selcol[i], "~ scantrack_market_descr + factor(year) + channel_type")
	myfit 	<- anova(lm(myfml, data=fmt_attr))
	tmp.tab	<- rbind(tmp.tab, cbind(Var = selcol[i], myfit))
}
cat("Variance decomposition of format attributes:\n"); print(tmp.tab); cat("\n")

if(write2csv){
	tmp.tab$component 	<- gsub("[0-9]","", rownames(tmp.tab))
	names(tmp.tab)		<- gsub("\\s","_", names(tmp.tab))
	tmp1		<- dcast(tmp.tab[,c("Var","Sum_Sq", "component")], Var ~ component, value.var="Sum_Sq")
	tmp2		<- dcast(tmp.tab[,c("Var","Mean_Sq", "component")], Var ~ component, value.var="Mean_Sq")
	tmp.tab1 	<- rbind(data.frame(tmp1, value="Sum_sq"), data.frame(tmp2, value="Mean_Sq"))
	tmp.tab1	<- tmp.tab1[order(tmp.tab1$Var),]
	tmp.tab1[,2:5]	<- 	tmp.tab1[,2:5]/rowSums(	tmp.tab1[,2:5])
	f 			<- file(tmpcsv, "at")
	writeLines("-----------------------------------------------------", f)
	writeLines("Variance decomposition of format attributes:\n", f)
	write.csv(tmp.tab1, f)
	close(f)
}

# Series of monthly price 
if(make_plot){
	selmkt 	<- "Boston"
	ggtmp 	<- subset(price_dat, scantrack_market_descr==selmkt)
	ggplot(ggtmp, aes(biweek, price_tag_norm_2004, col=channel_type)) + geom_point() + geom_line() + 
			labs(title = paste("Price series in ",selmkt, sep=""))	
}

##########################
# Simply summary statics # 
##########################
# Segment households 
panelist	<- data.table(hh_exp)
setkeyv(panelist, c("household_code","year","biweek"))
panelist	<- panelist[,list(first_income = income_group[1], first_famsize = famsize[1]), by=list(household_code)]
table(panelist$first_income, panelist$first_famsize)
seg_index	<- 1:(length(levels(panelist$first_income))*length(levels(panelist$first_famsize)))
names(seg_index) <- paste(rep(levels(panelist$first_income), length(levels(panelist$first_famsize))), 
						  rep(levels(panelist$first_famsize), each=length(levels(panelist$first_income)) ), sep="-")
tmp			<- paste(panelist$first_income, panelist$first_famsize, sep="-")
panelist$segment <- seg_index[tmp]
hh_exp		<- merge(hh_exp, data.frame(panelist)[,c("household_code", "segment","first_income")], by="household_code", all.x = TRUE)

selcol		<- c("household_code", "biweek", "num_day", "y", paste("SHR_", gsub("\\s", "_", fmt_name), sep=""))
tmp			<- merge(hh_exp[,selcol], data.frame(panelist)[,c("household_code","segment")], by = "household_code", all.x=T)
tmp			<- data.table(tmp)

tmp1		<- tmp[, list(Nhh = length(unique(household_code)), Nobs = length(biweek), mean_y = mean(y, na.rm=T),  mean_day = mean(num_day, na.rm=T), 
						  SHR_Convenience_Store = mean(SHR_Convenience_Store, na.rm=T), SHR_Discount_Store = mean(SHR_Discount_Store, na.rm=T), 
						  SHR_Dollar_Store = mean(SHR_Dollar_Store, na.rm=T), SHR_Drug_Store = mean(SHR_Drug_Store, na.rm=T), 
						  SHR_Grocery = mean(SHR_Grocery, na.rm=T), SHR_Warehouse_Club = mean(SHR_Warehouse_Club, na.rm=T)), by=list(segment)]
tmp2		<- tmp[, list(Nhh = length(unique(household_code)), Nobs = length(biweek), sd_y = sd(y, na.rm=T),  sd_day = sd(num_day, na.rm=T), 
						  SHR_Convenience_Store = sd(SHR_Convenience_Store, na.rm=T), SHR_Discount_Store = sd(SHR_Discount_Store, na.rm=T), 
						  SHR_Dollar_Store = sd(SHR_Dollar_Store, na.rm=T), SHR_Drug_Store = sd(SHR_Drug_Store, na.rm=T), 
						  SHR_Grocery = sd(SHR_Grocery, na.rm=T), SHR_Warehouse_Club = sd(SHR_Warehouse_Club, na.rm=T)), by=list(segment)]

tmp.tab 	<- matrix(NA, nrow(tmp1)*2, ncol(tmp1), dimnames = list(NULL, names(tmp1)))
sel			<- seq(1, 2*nrow(tmp1), 2)
tmp.tab[sel,] <- as.matrix(tmp1)
tmp.tab[(sel+1),] <- as.matrix(tmp2)
tmp.tab[(sel+1),c(2,3)] <- NA

# Format the table 
mydigits 	<- 2
selcol 		<- paste("SHR_", gsub("\\s", "_", fmt_name), sep="")
tmp.df		<- data.frame(tmp.tab)
for(i in 1:nrow(tmp.tab)){
	if(i%%2==1){
		tmp.df[i,selcol]	<- paste(round(tmp.tab[i,selcol]*100,mydigits), "%", sep="")
	}else{
		tmp.df[i,4:5]	<- paste("((", round(tmp.tab[i,4:5], mydigits), "))", sep="")
		tmp.df[i,selcol]	<- paste("((", round(tmp.tab[i,selcol]*100,mydigits), "%))", sep="")
	}
}

tmp.df[,"seg"]	<- names(seg_index)[tmp.df[,"segment"]]
tmp.df$FAMSIZE	<- lapply(strsplit(tmp.df$seg,"-"), function(x) x[2])
tmp.df$FAMSIZE <- factor(tmp.df$FAMSIZE, c("Single","Two", "Three+"))
ord			<- order(tmp.df$FAMSIZE, tmp.df$seg)
tmp.df		<- tmp.df[ord,]
cat("Summary statistics of data by segment:\n");print(tmp.df); cat("\n")

if(write2csv){	
	f 		<- file(tmpcsv, "at")
	writeLines("-----------------------------------------------------", f)
	writeLines("Summary statistics of anaylis data\n", f)
	write.csv(tmp.df,f)
	close(f)
}

##################
# Income pattern # 
##################
pan_yr		<- data.table(hh_exp)
pan_yr		<- pan_yr[,list(Income = unique(income_midvalue)), by = list(first_income, household_code, year)]
pan_yr		<- merge(pan_yr, cpi, by="year", all.x=T)
pan_yr$Income	 <- pan_yr$Income / pan_yr$cpi
pan_yr$ln_income <- log(pan_yr$Income)
pan_yr$recession <- 1 * (pan_yr$year >= 2008)
pan_yr		<- data.frame(pan_yr)

plots		<- list(NULL)
plots[[1]]	<- ggplot(pan_yr, aes(year, Income, linetype = first_income)) + geom_smooth(method=lm, formula = y ~ x + I(x^2) + I(x^3)) + 
					guides(linetype = guide_legend(title = "Segment"))

pan_yr$year <- factor(pan_yr$year)
summary(lm(Income ~ year*first_income, data = pan_yr))
myfit 		<- plm(Income ~ first_income*year, data=pan_yr, index = c("household_code","year"), model="within")
tmp			<- summary(myfit)
cat("Income trend regression:\n"); print(tmp); cat("\n")

myfit1		<- plm(Income ~ first_income*recession, data=pan_yr, index = c("household_code","year"), model="within")
tmp1		<- summary(myfit1)
cat("Income trend regression:\n"); print(tmp1); cat("\n")

if(write2csv){
	model_list 	<- list(year=data.frame(tmp$coefficients, Var = rownames(tmp$coefficients)), 
						recession = data.frame(tmp1$coefficients, Var = rownames(tmp1$coefficients)))
	tmp.tab		<- model_outreg(model_list, p.given = TRUE, head.name = c("Estimate","Std..Error", "Pr...t.."), digits = 4)
	f			<- file(tmpcsv, "at")
	writeLines("Fix-effect regressions of household income:", f)
	write.csv(tmp.tab, f)
	close(f)
}

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
plots[[2]]	<- ggplot(ggtmp, aes(year, Income, linetype = first_income)) + geom_line() + 
			    		labs(x = "Year", y = "Income($1000)") + 
						guides(linetype = guide_legend(title = "Segment"))

if(make_plot){
	pdf(paste(plot.wd,"/graph_income.pdf",sep=""), width = ww, height = ww*ar)
	print(plots[[1]])
	dev.off()

	pdf(paste(plot.wd,"/graph_income_within.pdf",sep=""), width = ww, height = ww*ar)
	print(plots[[2]])
	dev.off()	
}

#------------------------------------------------------#
# Income transition 
hh_exp$income_real	<- as.numeric(as.character(hh_exp$income_real))
n_Inc		<- 8
sel			<-  hh_exp$income_real >= 27
hh_exp[sel,"income_real"] <- 27
tmp			<- recode(hh_exp$income_real, hh_exp$income_midvalue, n_Inc)
Inc_nodes	<- sort(unique(tmp[,2]))
hh_exp$income_nodes <- tmp[,1]

# Lag income by 1 year
tmp		<- data.table(hh_exp)
tmp.lab	<- c("Under 25k", "25k-29,999", "30k-39,999","40k-49,999","50k-59,999","60k-69,999","70k-99,999","100k+")
tmp$income_nodes	<- factor(tmp$income_nodes, levels=1:n_Inc, labels=tmp.lab)
tmp		<- tmp[, list(income_nodes = unique(income_nodes)), by=list(first_income, household_code, recession, year)]
setkeyv(tmp, c("first_income","household_code","year"))
tmp		<- tmp[, income_forward := my_forward(income_nodes), by=list(first_income, household_code)]

# Transition table: replace 0 diagal with 1 if nobody belongs to some income level.
sel 	<- tmp$recession == 1
tmp.tab	<- rbind(mytransition(tmp[!sel,income_nodes], tmp[!sel,income_forward] ), 
				 mytransition(tmp[sel,income_nodes], tmp[sel,income_forward]))

ggtmp	<- rbind(data.frame(melt(tmp.tab[1:n_Inc,]), recession = 0), 
				 data.frame(melt(tmp.tab[-(1:n_Inc),]), recession = 1))				 
# ggplot(ggtmp, aes(Var2, value, fill = factor(recession))) + geom_bar(stat = "identity", position=position_dodge(width = .5)) + 
# 		facet_grid(Var1 ~ .) + labs(x = "Future income")
				
tmpn	<- nrow(tmp.tab)/2
ord		<- c(rbind(1:tmpn, 1:tmpn + tmpn))
tmp.tab	<- tmp.tab[ord,]
cat("Transition probability of income process:\n"); print(tmp.tab); cat("\n")

sel1	<- tmp$year == 2007
sel2	<- tmp$year == 2008
tmp.tab1<- rbind(mytransition(tmp[sel1,income_nodes], tmp[sel1,income_forward]), 
				 mytransition(tmp[sel2,income_nodes], tmp[sel2,income_forward]))
tmp.tab1<- tmp.tab1[ord,]
cat("Transition probability of income process (2007 vs 2008):\n"); print(tmp.tab); cat("\n")

if(write2csv){
	f	<- file(tmpcsv, "at")
	writeLines("Transition matrix of income:", f)
	write.csv(tmp.tab, f)
	writeLines("Transition matrix of income in 2007 vs 2008.", f)
	write.csv(tmp.tab1, f)
	close(f)
}

# Expectation and variance
tmp		<- data.table(hh_exp)
tmp.lab	<- c("Under 25k", "25k-29,999", "30k-39,999","40k-49,999","50k-59,999","60k-69,999","70k-99,999","100k+")
tmp$income_nodes	<- factor(tmp$income_nodes, levels=1:n_Inc, labels=tmp.lab)
tmp		<- tmp[, list(income_nodes = unique(income_nodes)), by=list(first_income, household_code, year)]
setkeyv(tmp, c("first_income","household_code","year"))
tmp		<- tmp[, income_forward := my_forward(income_nodes), by=list(first_income, household_code)]
tmp$income_forward	<- factor(tmp$income_forward, levels = 1:n_Inc, labels=tmp.lab)
tmp$seg	<- paste(tmp$first_income, tmp$year, sep="-")

tmp		<- data.frame(tmp)
sel		<- tmp$year != 2010
tmp1	<- split(tmp[sel,], tmp[sel,"seg"])
tmp2	<- lapply(tmp1, function(x) mytransition(x$income_nodes, x$income_forward))
tmp.mat	<- rep(1, n_Inc) %*% t(Inc_nodes)

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
	f			<- file(tmpcsv, "at")
	writeLines("Expectation and sd of future income:", f)
	write.csv(tmp.tab, f)
	close(f)
}

#############################
# Expenditure share pattern # 
#############################
# Raw data pattern
ggtmp		<- data.table(hh_exp)
ggtmp		<- subset(ggtmp, !is.na(first_income))
ggtmp		<- ggtmp[,list(DOLP_Convenience_Store = sum(DOLP_Convenience_Store), DOLP_Discount_Store = sum(DOLP_Discount_Store), 
							DOLP_Dollar_Store = sum(DOLP_Dollar_Store), DOLP_Drug_Store = sum(DOLP_Drug_Store), 
							DOLP_Grocery = sum(DOLP_Grocery), DOLP_Warehouse_Club = sum(DOLP_Warehouse_Club), y = sum(y)), 
					by = list(first_income, biweek)]
ggtmp		<- data.frame(ggtmp)
ggtmp[,3:(2+R)]	<- ggtmp[,3:(2+R)]/ggtmp$y
ggtmp		<- melt(ggtmp[,-ncol(ggtmp)], id.vars=c("first_income", "biweek"))
ggtmp		<- subset(ggtmp, !is.na(variable))
tmp			<- tapply(ggtmp$value, ggtmp$variable, mean)
tmp			<- names(tmp)[order(tmp,decreasing=T)]
tmp1		<- gsub("_", " ", gsub("DOLP_","",tmp))
ggtmp$variable	<- factor(as.character(ggtmp$variable), levels=tmp, labels=tmp1)
Rc_t		<- min(hh_exp[hh_exp$recession == 1, "biweek"])

if(make_plot){
	pdf(paste(plot.wd, "/graph_raw_share.pdf",sep=""), width=ww1, height = ww1*ar)
	print(
		ggplot(ggtmp, aes(biweek, value, fill = variable, order = as.numeric(variable))) + geom_bar(stat="identity") + 
				facet_wrap( ~ first_income) + 
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
ggtmp 		<- melt(hh_exp[,c("first_income","biweek",tmp_dv)], id.vars=c("first_income","biweek"))
ggtmp$variable <- factor(ggtmp$variable, levels = tmp_dv, labels = fmt_name)
if(make_plot){
	pdf(paste(plot.wd,"/graph_raw_share_gam.pdf",sep=""), width = ww, height = ww*.8)
	print(ggplot(ggtmp, aes(biweek, value, linetype = first_income)) + stat_smooth() + 
			facet_wrap(~ variable, ncol = 2, scales = "free_y") + 
			scale_y_continuous(labels=percent) + 
			labs(x = "Biweek", y = "Expenditure share") + 
			guides(linetype = guide_legend(title = "Segment")) + 
			theme_bw()
		)
	dev.off()
}

cat("------------------------------------------------------\n")
cat("Now fitting within-household trend nonparametrically.\n")
tmp_dv		<- paste("SHR_", gsub("\\s", "_", fmt_name), sep="")
plots		<- list(NULL)
reg_dat 	<- data.frame()
for(i in 1:length(fmt_name)){
	prc 	<- proc.time()
	ggtmp	<- within_loess(tmp_dv[i], "first_income", hh_exp)
	reg_dat <- rbind(reg_dat, ggtmp)
	plots[[i]] 	<- ggplot(ggtmp, aes(biweek, y, col=Segment)) + stat_smooth() + 
						geom_vline(xintercept = 100) + 
						labs(title = paste("Expenditure share at ", fmt_name[i], sep=""))
	use.time	<- proc.time() - prc
	cat("Finish", tmp_dv, "with", use.time[3]/60, "sec.\n")
}
reg_dat$retail <- gsub("SHR_","", reg_dat$DV)
reg_dat$retail <- gsub("_", " ", reg_dat$retail)

if(make_plot){
	pct		<- proc.time()
	pdf(paste(plot.wd,"/graph_within_gam.pdf",sep=""), width=ww, height = ww*.8)
	print( ggplot(reg_dat, aes(biweek, y, linetype=Segment)) + stat_smooth() + 
						geom_vline(xintercept = Rc_t, size = .5, linetype = 2) + 
						facet_wrap( ~ retail, ncol=2, scales = "free_y") + 
						scale_y_continuous(labels=percent) + 
						labs(x = "Biweek", y = "Expenditure share") + 
						theme_bw())
	dev.off()
	use.time <- proc.time() - prc
	cat("Ploting loess took", use.time[3]/60,"min.\n")	
}

# Trend of shopping strategies
tmp_dv1		<- c("dol_purchases", "num_day", "num_trip")
reg_dat1 	<- data.frame()
for(i in 1:length(tmp_dv1)){
	prc 		<- proc.time()
	ggtmp		<- within_loess(tmp_dv1[i], "first_income", hh_exp)
	reg_dat1 	<- rbind(reg_dat1, ggtmp)
	use.time	<- proc.time() - prc
	cat("Finish", tmp_dv1[i], "with", use.time[3]/60, "min.\n")
}
reg_dat1$dv		<- factor(reg_dat1$DV, levels=tmp_dv1, labels = c("Expenditure","Number of\n shopping days","Number\n trips"))

if(make_plot){
	pdf(paste(plot.wd,"/graph_within_strategy_gam.pdf",sep=""), width=ww, height = ww*ar)
	print( ggplot(subset(reg_dat1, DV!="num_trip"), aes(biweek, y, linetype=Segment)) + stat_smooth() + 
						geom_vline(xintercept = Rc_t, size = .5, linetype = 2) + 
						facet_wrap( ~ dv, ncol=2, scales = "free_y") + 
						labs(x = "Biweek", y = "Value")+ 
						theme_bw())
	dev.off()
}

######################
# Mediation analysis #
######################
# Plot conditional expenditure share. 
seq_y 	<- seq(10, 500, 10)
selcol	<- paste("SHR_", gsub("\\s","_", fmt_name),sep="")
tmp		<- subset(hh_exp[,c("first_income", "dol_purchases", "d", selcol)], dol_purchases > 0)

tmp$ybin<- cut(tmp$dol_purchases, c(0,seq_y), include.lowest = T, label=seq_y)
tmp[is.na(tmp$ybin), "ybin"] <- max(seq_y)
tmp		<- data.table(melt(tmp[,c("d","ybin",selcol)], id.vars=c("ybin","d")))
tmp		<- tmp[,list(value = mean(value)), by = list(ybin, d, variable)]
tmp$d	<- factor(tmp$d, levels=1:5, labels=c("1", "2", "3","4-5","6+"))
tmp1	<- tapply(tmp$value, tmp$variable, mean)
ord		<- order(tmp1)
tmp$variable <- factor(tmp$variable, levels = selcol[ord], labels=fmt_name[ord])

if(make_plot){
	tmp1	<- c(10, 250, 500)
	pdf(paste(plot.wd, "/graph_conditional_share.pdf",sep=""), width =ww1, height = ww1*ar)
	print(ggplot(tmp, aes(ybin, value, fill = variable, order = as.numeric(variable))) + geom_bar(stat ="identity") + 
				facet_wrap(~d, nrow = 1) + 
				scale_x_discrete(breaks = tmp1) + 
				scale_y_continuous(labels=percent) + 
				labs(x = "Expenditure", y = "Expenditure share") + 
				scale_fill_brewer() + 
				theme(legend.position = "bottom") + 
				guides(fill = guide_legend(title = "Retail format", reverse = TRUE, nrow=2)) 
				# theme_bw() 
		) 
	dev.off()
}

# Fixed effect regressions
cat("------------------------------------------------------\n")
cat("Now fitting within-household regressions for expenditure share.\n")
tmpdata 	<- hh_exp
for(i in tmp_dv){
	tmpdata[,i] <- tmpdata[,i]*100
}
tmpdata$ln_income	<- log(tmpdata$income_midvalue)
tmpdata$month		<- factor(tmpdata$month)
tmpdata$ln_dolpurchases <- log(tmpdata$dol_purchases)
tmpdata$recession 	<- factor(tmpdata$recession)

# Add income uncertainty
tmp		<- data.table(hh_exp)
tmp$income_nodes	<- factor(tmp$income_nodes)
tmp		<- tmp[, list(income_nodes = unique(income_nodes)), by=list(first_income, household_code, year)]
tmp$recession <- ifelse(tmp$year >= 2008, 1, 0)
tmp		<- tmp[, income_forward := my_forward(income_nodes), by=list(household_code)]
tmp		<- data.frame(tmp)

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
tmpdata	<- merge(tmpdata, trans.dat, by = c("first_income","year","income_nodes"), all.x=T)

reg_dat2	<- data.frame()
myfml		<- list(indirect 	= as.formula("y ~ ln_income + p_low + p_high + recession + month"), 
					direct		= as.formula("y ~ ln_dolpurchases + num_day + month"), 
					mediation	= as.formula("y ~ ln_dolpurchases + num_day + ln_income + p_low + p_high + recession + month"), 
					mediationiv	= as.formula("y ~ ln_dolpurchases + num_day + ln_income + p_low + p_high + recession + month|
											. - ln_dolpurchases - num_day + lag(dol_purchases,1) + lag(num_day, 1) + lag(dol_purchases,2) + lag(num_day, 2)"),
					indirect_int= as.formula("y ~ ln_income:first_income + recession:first_income + p_low:first_income + p_high:first_income + month"), 
					direct_int	= as.formula("y ~ ln_dolpurchases:first_income + num_day:first_income + month"), 
					mediation_int= as.formula("y ~ ln_dolpurchases:first_income + num_day:first_income + 
											ln_income:first_income + recession:first_income + p_low:first_income + p_high:first_income + month"), 
					mediationiv_int = as.formula("y ~ ln_dolpurchases:first_income + num_day:first_income + 
												ln_income:first_income + recession:first_income + p_low:first_income + p_high:first_income + month |
												. - ln_dolpurchases:first_income - num_day:first_income + 
												lag(dol_purchases,1):first_income+ lag(num_day,1):first_income + lag(dol_purchases,2):first_income+ lag(num_day,2):first_income")
					)
for(i in 1:length(fmt_name)){
	prc 		<- proc.time()
	for(j in 1:length(myfml)){
		tmp 	<- as.formula(substitute(y ~ x, list(y = as.name(tmp_dv[i]), x = terms(myfml[[j]])[[3]])) )
		myfit	<- plm(tmp, data = tmpdata, index = c("household_code","biweek"), model="within")
		tmp2	<- summary(myfit)
		tmp3	<- data.frame(retail = fmt_name[i], tmp2$coefficients, model = names(myfml)[j], adjrsq = tmp2$r.squared["adjrsq"])
		tmp3$Var<- rownames(tmp3)
		rownames(tmp3) <- NULL
		reg_dat2<- rbind(reg_dat2, tmp3)
	}
	use.time 	<- proc.time() - prc
	cat("Regressions for", fmt_name[i], "finishes, using", use.time[3]/60, "min.\n")
}

cat("------------------------------------------------------\n")
cat("Now fitting within-household regressions for shopping strategies.\n")
reg_dat3 	<- data.frame()
tmp			<- c("ln_dolpurchases", "num_day")
sel			<- !is.na(tmpdata$ln_dolpurchases)
for(i in 1:length(tmp)){
	prc			<- proc.time()
	for(j in c(1,5)){
		tmp1	<- update(myfml[[j]], as.formula(paste(tmp[i], "~ .")))
		myfit	<- plm(tmp1, data = tmpdata[sel,], index = c("household_code","biweek"), model="within")
		tmp2	<- summary(myfit)
		tmp3	<- data.frame(dv = tmp[i], tmp2$coefficients, model = names(myfml)[j], adjrsq = tmp2$r.squared["adjrsq"])
		tmp3$Var<- rownames(tmp3)
		rownames(tmp3) <- NULL
		reg_dat3<- rbind(reg_dat3, tmp3)		
	}
	use.time	<- proc.time() - prc
	cat("Regressions for", tmp[i], "finishes, using", use.time[3]/60, "min.\n")
}

# Organize regression results 
sel			<- substr(reg_dat2$Var, 1, 5) == "month"
tmp_df 		<- reg_dat2[!sel,]
tmp.ret		<- sort(unique(reg_dat2$retail))

# # Unify names
unique(tmp_df$Var)
tmp1	<- sapply(tmp_df$Var, function(x) strsplit(x, ":")[[1]][1])
tmp2	<- sapply(tmp_df$Var, function(x) strsplit(x, ":")[[1]][2])
sel		<- substr(tmp1, 1, 12) == "first_income" & !is.na(tmp2)
tmp_df[sel,"Var"] <- paste(tmp2[sel], tmp1[sel], sep=":")
myvarord	<- unique(tmp_df$Var)

sel			<- substr(reg_dat3$Var, 1, 5) == "month"
tmp			<- reg_dat3[!sel,]
names(tmp)[1] <- "retail"
tmp1	<- sapply(tmp$Var, function(x) strsplit(x, ":")[[1]][1])
tmp2	<- sapply(tmp$Var, function(x) strsplit(x, ":")[[1]][2])
sel		<- substr(tmp1, 1, 12) == "first_income" & !is.na(tmp2)
tmp[sel,"Var"] <- paste(tmp2[sel], tmp1[sel], sep=":")

tmpls		<- vector("list", length(tmp.ret))
names(tmpls) <- sort(unique(reg_dat2$retail))

for(i in 1:length(tmp.ret)){
	tmp1			<- rbind(subset(tmp_df, retail == tmp.ret[i]), tmp)
	tmp1$seg		<- with(tmp1, paste(retail, model, sep="_"))
	model_list		<- split(tmp1, tmp1$seg)
	tmp.tab			<- model_outreg(model_list, p.given = TRUE, head.name = c("Estimate","Std..Error","Pr...t..","Var"), digits = 4, varord = myvarord)
	tmpls[[i]]		<- tmp.tab
}

if(write2csv){
	f	<- file(tmpcsv, "at")
	for(i in 1:length(tmpls)){
		writeLines(paste("\n", names(tmpls)[i], "\n",sep=""), f)
		write.csv(tmpls[[i]], f)
	}
	close(f)
}

save.image(file = "4_sumstat_larger.rdata")

cat("This program is done.")


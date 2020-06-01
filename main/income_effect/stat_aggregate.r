# This file analyzes changes of the aggregated market shares

library(ggplot2)
library(reshape2)
library(data.table)
library(plm)
library(gridExtra)
library(scales)
library(systemfit)
library(stargazer)
library(gridExtra)
library(mvProbit)

setwd("~/Documents/Research/Store switching/processed data")
plot.wd	<- '/Users/chaoqunchen/Desktop'
source('~/Documents/Research/Store switching/Exercise/main/income_effect/plm_vcovHC.R')

week.price	<- FALSE
cpi.adj		<- TRUE
write2csv	<- FALSE
make_plot	<- TRUE
fname		<- "expenditure_reg_sub"
if(cpi.adj) { fname <- paste(fname, "_cpi", sep="")}
if(week.price){ fname <- paste(fname, "_wkprc", sep="")}

if(week.price){
	load("hh_biweek_exp_20150812.rdata")
}else{
	load("hh_biweek_exp.rdata")
}
codebook	<- read.csv("code_book.csv")
pan_org		<- read.csv("2_panelists_20160606.csv")
sel			<- grep("panel_year", names(pan_org))
names(pan_org)[sel]	<- "year"
age.code	<- setNames(c(0, 22, 27, 32, 37, 42, 47, 52, 60, 65), 
						0:9)
pan_org$female_head_age	<- age.code[as.character(pan_org$female_head_age)]						
pan_org$male_head_age	<- age.code[as.character(pan_org$male_head_age)]						

################# 
# Organize data # 
#################
# Retail formats 
fmt_name 	<- as.character(sort(unique(fmt_attr$channel_type)))
R			<- length(fmt_name)
fmt_order	<- c("Grocery", "Discount Store", "Warehouse Club", "Drug Store", "Dollar Store", "Convenience Store")

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
# Demographics 
pan_org$age <- rowMeans(cbind(ifelse(pan_org$female_head_age==0, NA, pan_org$female_head_age), 
						ifelse(pan_org$male_head_age==0, NA, pan_org$male_head_age)),na.rm=T)
pan_org$retire	<- ifelse((pan_org$female_head_occupation==12| pan_org$male_head_occupation ==12), 1, 0)

# Expenditure 
pan_yr		<- data.table(hh_exp)
pan_yr		<- pan_yr[,list(Income = unique(income_midvalue), dol = sum(dol), 
							DOL_Grocery = sum(DOL_Grocery), DOL_Discount_Store = sum(DOL_Discount_Store), 
							DOL_Warehouse_Club = sum(DOL_Warehouse_Club), DOL_Dollar_Store = sum(DOL_Dollar_Store), 
							DOL_Drug_Store = sum(DOL_Drug_Store), DOL_Convenience_Store = sum(DOL_Convenience_Store)), 
						by = list(first_incomeg, household_code, year, cpi)]
setkeyv(pan_yr, c("household_code", "year"))
pan_yr		<- pan_yr[, ':='(ln_income = log(Income), cpi_adj_income = Income/cpi, recession = 1*(year >= 2008))]
pan_yr    	<- pan_yr[,':='(first_income = Income[1], year_id= year - year[1] + 1, tenure = length(year), move = c(0, 1*(diff(Income)!=0)))
                    , by = list(household_code)]
pan_yr		<- pan_yr[,move_all := sum(move), by = list(household_code)]
pan_yr 	<- pan_yr[, last.income := Income[length(Income)], by = list(household_code)]
pan_yr$last_incomeg	<- cut(pan_yr$last.income, c(0, 37500, 65000, 200000), labels = c("Low", "Med", "High"))
pan_yr		<- pan_yr[, period:=ifelse(year_id==1, "first", ifelse(year_id==max(year_id), "last", "")), by = list(household_code)]
table(pan_yr$first_incomeg, pan_yr$last_incomeg)
table(pan_yr[tenure ==7,first_incomeg], pan_yr[tenure ==7, last_incomeg])

# Merge in demographics 
dim(pan_yr)
pan_yr		<- merge(data.frame(pan_yr), pan_org[,c("household_code", "year", "age", "household_size", "NumChild", "retire")],
					by = c("household_code","year"),all.x =T)					
dim(pan_yr)					
pan_yr		<- data.table(pan_yr)
setkeyv(pan_yr, c("first_incomeg", "household_code","year"))

panelist	<- pan_yr[, list(tenure = unique(tenure), p.change = unique(move_all/tenure), 
							large.init.change = max(abs(Income-Income[1])/Income[1] ), 
							delta.lastincome = Income[length(Income)] - Income[(length(Income)-1)], 
							growth = mean(diff(Income)/Income[-length(Income)], na.rm=TRUE), 
							age = age[1], retire = 1*any(retire==1) ), 
						by = list(first_incomeg, household_code)]
summary(panelist)

# Aggregate shares.
# sel			<- pan_yr$tenure == 7 
# sel			<- pan_yr$tenure == 7 & pan_yr$household_code %in% panelist[p.change>.3,household_code]
sel				<- with(pan_yr, tenure == 7 & ((first_incomeg=="Low"&last_incomeg=="High") | (first_incomeg=="High"&last_incomeg=="Low")))

length(unique(pan_yr[sel,household_code]))
agg.data	<- pan_yr[sel, list(DOL_Grocery = sum(DOL_Grocery), DOL_Discount_Store = sum(DOL_Discount_Store), 
							DOL_Warehouse_Club = sum(DOL_Warehouse_Club), DOL_Dollar_Store = sum(DOL_Dollar_Store), 
							DOL_Drug_Store = sum(DOL_Drug_Store), DOL_Convenience_Store = sum(DOL_Convenience_Store)), 
							by = list(year)]
tmp.tab	<- as.matrix(agg.data)
tmp.tab[,-1]	<- tmp.tab[,-1]/rowSums(tmp.tab[,-1])
colnames(tmp.tab)	<- gsub("DOL_", "", colnames(tmp.tab))
cat("Share over year:\n"); print(cbind(tmp.tab[,1], round(tmp.tab[,-1]*100, 2))); cat("\n")

# By income group and recession
agg.data	<- pan_yr[sel, list(DOL_Grocery = sum(DOL_Grocery), DOL_Discount_Store = sum(DOL_Discount_Store), 
							DOL_Warehouse_Club = sum(DOL_Warehouse_Club), DOL_Dollar_Store = sum(DOL_Dollar_Store), 
							DOL_Drug_Store = sum(DOL_Drug_Store), DOL_Convenience_Store = sum(DOL_Convenience_Store)), 
							by = list(first_incomeg, recession)]
tmp.tab1	<- as.matrix(agg.data[, paste("DOL_", gsub("\\s","_",fmt_name),sep=""), with=FALSE])
tmp.tab1	<- tmp.tab1/rowSums(tmp.tab1)
cat("Share before and after the recession:\n"); print(cbind(agg.data[,list(first_incomeg, recession)], round(tmp.tab1*100,2))); cat("\n")

# By income group and 07 vs. 08
agg.data	<- pan_yr[sel & year %in% c(2007, 2008), list(DOL_Grocery = sum(DOL_Grocery), DOL_Discount_Store = sum(DOL_Discount_Store), 
							DOL_Warehouse_Club = sum(DOL_Warehouse_Club), DOL_Dollar_Store = sum(DOL_Dollar_Store), 
							DOL_Drug_Store = sum(DOL_Drug_Store), DOL_Convenience_Store = sum(DOL_Convenience_Store)), 
							by = list(first_incomeg, year)]
tmp.tab1	<- as.matrix(agg.data[, paste("DOL_", gsub("\\s","_",fmt_name),sep=""), with=FALSE])
tmp.tab1	<- tmp.tab1/rowSums(tmp.tab1)
cat("Share before and after the 2008:\n"); print(cbind(agg.data[,list(first_incomeg, year)], round(tmp.tab1*100,2))); cat("\n")

# By income change direction and 07 vs. 08
agg.data	<- pan_yr[year %in% c(2007,2008),]
agg.data	<- agg.data[,sign:=sign(Income[year==2008]-Income[year==2007]), by = list(household_code)]
agg.data	<- agg.data[!is.na(sign),]
tmp1		<- agg.data[, list(DOL_Grocery = sum(DOL_Grocery), DOL_Discount_Store = sum(DOL_Discount_Store), 
							DOL_Warehouse_Club = sum(DOL_Warehouse_Club), DOL_Dollar_Store = sum(DOL_Dollar_Store), 
							DOL_Drug_Store = sum(DOL_Drug_Store), DOL_Convenience_Store = sum(DOL_Convenience_Store)), 
							by = list(sign, year)]
tmp2		<- agg.data[, list(DOL_Grocery = sum(DOL_Grocery), DOL_Discount_Store = sum(DOL_Discount_Store), 
							DOL_Warehouse_Club = sum(DOL_Warehouse_Club), DOL_Dollar_Store = sum(DOL_Dollar_Store), 
							DOL_Drug_Store = sum(DOL_Drug_Store), DOL_Convenience_Store = sum(DOL_Convenience_Store)), 
							by = list(first_incomeg, sign, year)]
tmp.tab1	<- as.matrix(tmp1[, paste("DOL_", gsub("\\s","_",fmt_name),sep=""), with=FALSE])
tmp.tab1	<- tmp.tab1/rowSums(tmp.tab1)
tmp.tab2	<- as.matrix(tmp2[, paste("DOL_", gsub("\\s","_",fmt_name),sep=""), with=FALSE])
tmp.tab2	<- tmp.tab2/rowSums(tmp.tab2)
ggtmp		<- rbind(data.frame(first_incomeg = "All", tmp1[,list(sign, year)], tmp.tab1), 
					data.frame(tmp2[, list(first_incomeg, sign, year)], tmp.tab2))
ggtmp		<- melt(ggtmp, id.var = c("first_incomeg", "sign","year"))
ggtmp$sign	<- factor(ggtmp$sign, levels = seq(-1, 1, 1), labels = c("Decrease", "No Change", "Increase"))
ggtmp$variable	<- factor(ggtmp$variable, levels = c("DOL_Grocery", "DOL_Discount_Store", "DOL_Warehouse_Club", 
													"DOL_Drug_Store", "DOL_Dollar_Store", "DOL_Convenience_Store"), 
						labels = c("Grocery", "Discount Store", "Warehouse Club", "Drug Store", "Dollar Store", "Convenience Store"))

ggplot(subset(ggtmp, first_incomeg=="All"), aes(year, value*100, fill = variable)) + 
			geom_bar(stat = "identity") + facet_wrap(~sign)
ggplot(subset(ggtmp, first_incomeg!="All"), aes(year, value*100, fill = variable)) + 
			geom_bar(stat = "identity") + facet_grid(first_incomeg~sign)
			
# Plot the difference
ggtmp1	<- data.frame(sign = factor(seq(-1, 1, 1), labels = c("Decrease", "No Change", "Increase")), tmp.tab1[tmp1$year==2008,] - tmp.tab1[tmp1$year==2007,])			
ggtmp1	<- melt(ggtmp1, id.var = "sign")
ggtmp2	<- data.frame(unique(tmp2[,list(first_incomeg,sign)]), tmp.tab2[tmp2$year==2008,]- tmp.tab2[tmp2$year==2007,])
ggtmp2	<- melt(ggtmp2, id.var = c("first_incomeg", "sign"))
ggtmp1$variable	<- factor(ggtmp1$variable, levels = paste("DOL_", gsub(" ", "_", fmt_order), sep=""), labels = fmt_order)
ggtmp2$variable	<- factor(ggtmp2$variable, levels = paste("DOL_", gsub(" ", "_", fmt_order), sep=""), labels = fmt_order)
ggtmp2$sign		<- factor(ggtmp2$sign, levels = seq(-1,1,1), labels = c("Decrease", "No Change", "Increase"))


ggplot(ggtmp1, aes(sign, value*100, fill = variable)) + geom_bar(position = "dodge", stat = "identity") + coord_flip()
ggplot(ggtmp2, aes(sign, value*100, fill = variable)) + geom_bar(position = "dodge", stat = "identity") + 
		coord_flip() + 
		facet_grid(.~first_incomeg)


# Check demographic changes
agg.data	<- pan_yr[sel, list(age = mean(age), famsize = mean(household_size), 
								NumChild = mean(NumChild), retire = mean(retire)), 
						by = list(first_incomeg, recession)]
agg.data	<- data.frame(agg.data)						
cat("Demograhpic shifts before and after the recession:\n"); print(cbind(agg.data[,1:2], round(agg.data[,-c(1:2)],2))); cat("\n")

agg.data	<- pan_yr[sel & year %in% c(2004, 2010), list(age = mean(age), famsize = mean(household_size), 
								NumChild = mean(NumChild), retire = mean(retire)), 
						by = list(first_incomeg, year)]
agg.data	<- data.frame(agg.data)						
cat("Demograhpic shifts between the first and last year:\n"); print(cbind(agg.data[,1:2], round(agg.data[,-c(1:2)],2))); cat("\n")


tmp 	<- dcast(data.frame(pan_yr[sel & year %in% c(2004, 2010), list(retire, household_code, first_incomeg, year)]), 
					first_incomeg+household_code~year, value.var = "retire")
by(tmp, tmp$first_incomeg, function(x) table(x[,3:4]))

# Delta changes
# Compute expenditure share
pan_yr		<- data.frame(pan_yr)
sel			<- paste("DOL_", gsub("\\s", "_", fmt_name), sep="")
sel1		<- gsub("DOL", "SHR", sel)
for(i in 1:length(fmt_name)){
	pan_yr[,sel1[i]] <- pan_yr[,sel[i]]/pan_yr$dol
}

# Demean share 
selcol		<- c("dol", sel1)
pan_yrd		<- data.table(pan_yr)
large.shock	<- 2
for(i in selcol){
	pan_yrd	<- pan_yrd[, eval(as.name(i)):=eval(as.name(i)) - eval(as.name(i))[1], by = list(household_code)]
}
pan_yrd		<- pan_yrd[, ':='(d_Income = Income - Income[1]), by = list(household_code)]
pan_yrd		<- pan_yrd[, d_Income := d_Income/10000]
pan_yrd		<- pan_yrd[, ':='(large = factor(1*(d_Income>large.shock)), retire_age = factor(1*(age>=60)),
 							 ADP_Convenience_Store 	= factor(1*(DOL_Convenience_Store[1]>0)), 
							ADP_Discount_Store 		= factor(1*(DOL_Discount_Store[1]>0)), 
							ADP_Dollar_Store 		= factor(1*(DOL_Dollar_Store[1]>0)),
							ADP_Drug_Store 			= factor(1*(DOL_Drug_Store[1]>0)), 
							ADP_Grocery 			= factor(1*(DOL_Grocery[1]>0)), 
							ADP_Warehouse_Club 		= factor(1*(DOL_Warehouse_Club[1]>0)) ), 
						by = list(household_code)]
pan_yrd$year <- factor(pan_yrd$year)

# sel		<- pan_yrd$tenure == 7 & pan_yrd$household_code %in% panelist[large.init.change<.8, household_code]
# sel		<- pan_yrd$tenure == 7 & pan_yrd$household_code %in% panelist[growth<.5, household_code]
sel 	<- pan_yrd$year_id > 1
fml.ls	<- list(IncGrp 		= y ~ d_Income*first_incomeg, 
				ShockSize	= , 
				Adoption	= , 
				RetireAge	= )
tmplist	<- setNames(vector("list", R+1), selcol)
for(i in 1:(1+R)){
	fml	<- as.formula(paste(selcol[i], "~d_Income", sep=""))
	tmplist[[i]] <- lm(fml, data = pan_yrd[sel,])
}
stargazer(tmplist, type = "text")

tmplist1	<- setNames(vector("list", R+1), selcol)
for(i in 1:(R+1)){
	fml	<- as.formula(paste(selcol[i], "~d_Income:first_incomeg", sep=""))
	tmplist1[[i]] <- lm(fml, data = pan_yrd[sel,])
}
stargazer(tmplist1, type = "text")

# Period-to-period difference 
selcol		<- c("dol", sel1)
pan_yrd		<- data.table(pan_yr)
setkeyv(pan_yrd, c("household_code", "year"))
for(i in selcol){
	# pan_yrd	<- pan_yrd[, eval(as.name(i)):=, by = list(household_code)]
	pan_yrd	<- pan_yrd[,eval(as.name(i)):=c(NA, diff(eval(as.name(i))))]
	pan_yrd	<- pan_yrd[year_id==1, eval(as.name(i)):=NA]
}
pan_yrd		<- pan_yrd[, ':='(d_Income1 = c(NA, diff(log(Income))), d_Income2 = c(NA, NA, diff(log(Income), 2))), by = list(household_code)]
pan_yrd		<- pan_yrd[, d_Income := (d_Income1+ d_Income2)/1000]

# State dependence of incidence 
tmp			<- as.matrix(pan_yr[,paste("DOL_", gsub("\\s","_", fmt_name), sep=""), with = FALSE])
colnames(tmp)	<- gsub("DOL", "ID", colnames(tmp))
tmp			<- 1*(tmp>0)
colMeans(tmp)
pan_yr		<- cbind(pan_yr, data.frame(tmp))
setkeyv(pan_yr, c("household_code", "year"))
pan_yr		<- pan_yr[, ':='(lag_ID_Convenience=c(NA, ID_Convenience_Store[-length(year)]), lag_ID_Discount=c(NA, ID_Discount_Store[-length(year)]), 
							 lag_ID_Dollar=c(NA, ID_Dollar_Store[-length(year)]), lag_ID_Drug=c(NA, ID_Drug_Store[-length(year)]), 
							 lag_ID_Grocery=c(NA, ID_Grocery[-length(year)]), lag_ID_Warehouse=c(NA, ID_Warehouse_Club[-length(year)])), 
							by = list(household_code)]

sel		<- pan_yr$tenure == 7
fit.ls	<- vector("list", 10)
tmpstart	<- rep(c("Low", "Med", "High"), each = 3)
tmpend		<- rep(c("Low", "Med", "High"), 3)
for(i in 1:length(tmpstart)){
	sel1	<- sel & pan_yr$first_incomeg==tmpstart[i] &pan_yr$last_incomeg==tmpend[i]
	fit.ls[[i]]	<- glm(ID_Warehouse_Club ~ lag_ID_Warehouse+Income, data = pan_yr[sel1,], family = binomial(link = "logit") )
}
fit.ls[[10]]	<- glm(ID_Warehouse_Club ~ lag_ID_Warehouse+Income, data = pan_yr[sel,], family = binomial(link = "logit") )
stargazer(fit.ls, type = "text", column.labels = c(paste(substr(tmpstart,1,1), "->", substr(tmpend,1,1), sep=""), "all"))

# Multivariate probit
sel		<- pan_yr$tenure == 7 & pan_yr$first_incomeg == "Low" & pan_yr$last_incomeg == "High"  & !is.na(pan_yr$lag_ID_Warehouse)
est0 	<- mvProbit( cbind(ID_Convenience_Store,ID_Dollar_Store,ID_Drug_Store,ID_Warehouse_Club) 
						~ lag_ID_Convenience + lag_ID_Dollar + lag_ID_Drug + lag_ID_Warehouse,
   					data = as.data.frame(pan_yr[sel,]), iterlim = 1, nGHK = 50 )
summary(est0)
matrix(coef(est0)[grep("b", names(coef(est0)))], nrow = 4, 
		dimnames = list(c("ID_Convenience_Store","ID_Dollar_Store","ID_Drug_Store","ID_Warehouse_Club"), 
						c("Intercept", "lag_ID_Convenience", "lag_ID_Dollar", "lag_ID_Drug", "lag_ID_Warehouse")  ))

est1 	<- mvProbit( cbind(ID_Convenience_Store,ID_Dollar_Store,ID_Drug_Store,ID_Warehouse_Club) 
						~ age + household_size + NumChild + factor(retire),
   					data = as.data.frame(pan_yr[sel,]), iterlim = 1, nGHK = 50 )
summary(est1)


sel		<- pan_yr$tenure == 7 & pan_yr$first_incomeg == "High" & pan_yr$last_incomeg == "Low"  & !is.na(pan_yr$lag_ID_Warehouse)
est2 	<- mvProbit( cbind(ID_Convenience_Store,ID_Dollar_Store,ID_Drug_Store,ID_Warehouse_Club) 
						~ lag_ID_Convenience + lag_ID_Dollar + lag_ID_Drug + lag_ID_Warehouse,
   					data = as.data.frame(pan_yr[sel,]), iterlim = 1, nGHK = 50 )
summary(est2)
est3 	<- mvProbit( cbind(ID_Convenience_Store,ID_Dollar_Store,ID_Drug_Store,ID_Warehouse_Club) 
						~ age + household_size + NumChild + factor(retire),
   					data = as.data.frame(pan_yr[sel,]), iterlim = 1, nGHK = 50 )
summary(est3)


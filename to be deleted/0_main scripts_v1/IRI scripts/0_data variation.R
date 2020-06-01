library(lme4)
library(ggplot2)
library(reshape2)
library(lfe)
library(data.table)

setwd("~/Documents/Research/Store switching/processed data")
my.wd <- "~/Desktop"
source('~/Documents/Research/Store switching/Exercise/outreg function.R')

setwd("E:/Users/ccv103/Desktop")
my.wd <- "E:/Users/ccv103/Desktop"

exp_month <- read.csv("exp_month.csv",header=T)
exp_week <- read.csv("exp_week.csv",header=T)
tmp.csv <- paste(my.wd,"/tmp.csv",sep="")

# Volatility data 
exp_month$recession <- factor(exp_month$recession,levels=c(0,1))
exp_month$PANID		<- factor(exp_month$PANID)
exp_month$FAMSIZE1	<- factor(exp_month$FAMSIZE1,levels=c("Single","Two","Three+"))
exp_month$INCOME1	<- factor(exp_month$INCOME1,levels=c("Low","Med","High"))

exp_week$recession 	<- factor(exp_week$recession,levels=c(0,1))
exp_week$PANID		<- factor(exp_week$PANID)
exp_week$FAMSIZE1	<- factor(exp_week$FAMSIZE1,levels=c("Single","Two","Three+"))
exp_week$INCOME1	<- factor(exp_week$INCOME1,levels=c("Low","Med","High"))

vol_month <- data.table(exp_month)
vol_month <- vol_month[,list(nobs=length(BSKTDOL),vol = sd(BSKTDOL)),by=list(PANID,INCOME1,FAMSIZE1,recession)]
vol_week <- data.table(exp_week)
vol_week <- vol_week[,list(nobs=length(BSKTDOL),vol = sd(BSKTDOL), vol_n0 = sd(BSKTDOL[BSKTDOL>0])),by=list(PANID,INCOME1,FAMSIZE1,recession)]

f <- file(tmp.csv,"w")
writeLines("Variation in the expenditure",f)
close(f)
f <- file(tmp.csv,"at")
pan_char <- data.frame(PANID=sort(unique(exp_month$PANID)))
Nh <- nrow(pan_char)

######################
# Expenditure change # 
######################
# Montly data
# Expenditure change for average households
fit0 <- felm(BSKTDOL ~ recession |PANID + MONTH, data=exp_month)
fit1 <- felm(BSKTDOL ~ recession*INCOME1 |PANID + MONTH, data=exp_month)
fit2 <- lmer(BSKTDOL ~ factor(MONTH) + recession + (1 + recession|PANID),exp_month)
pan_char$month_Dexp_RE <- ranef(fit2)[[1]][,"recession1"] + fixef(fit2)['recession1']
pan_char$month_Dexp_dr <- with(pan_char,ifelse(month_Dexp_RE>0,"Exp increase","Exp decrease"))
table(pan_char$month_Dexp_dr)
with(pan_char,by(month_Dexp_RE,month_Dexp_dr,mean))

# # Individual t-test 
# pan_char$month_Dexp		<- NA
# pan_char$month_Dexp_t 	<- NA
# pan_char$month_Dexp_TL 	<- NA
# pan_char$month_Dexp_TU 	<- NA
# for(i in 1:nrow(pan_char)){
# 	tmp1 <- subset(exp_month,PANID == pan_char[i,"PANID"] & recession == 0)[,"BSKTDOL"]
# 	tmp2 <- subset(exp_month,PANID == pan_char[i,"PANID"] & recession == 1)[,"BSKTDOL"]
# 	tmp  <- t.test(tmp2,tmp1)
# 	pan_char[i,"month_Dexp"] <- tmp$estimate[1] - tmp$estimate[2]
# 	pan_char[i,"month_Dexp_t"] <- tmp$statistic
# 	pan_char[i,c("month_Dexp_TL","month_Dexp_TU")] <- tmp$conf.int[1:2]
# 	if(i%%500==0){
# 		print(i)
# 	}
# }
# pan_char$month_Dexp_dr <- 0			# No change
# sel <- with(pan_char, month_Dexp_TL<0 & month_Dexp_TU<0 )
# pan_char[sel,"month_Dexp_dr"] <- -1
# sel <- with(pan_char, month_Dexp_TL>0 & month_Dexp_TU>0 )
# pan_char[sel,"month_Dexp_dr"] <- 1
# tmp.tab <- table(pan_char$month_Dexp_dr)/Nh
# tmp.tab <- rbind(tmp.tab,c(with(pan_char,by(month_Dexp,month_Dexp_dr,mean))))
# rownames(tmp.tab) <- c("Percentage","Change level")

writeLines("-----------------------------------------------",f)
writeLines("Change of expenditure with monthly data",f)
writeLines("Expenditure change for average households",f)
write.csv(summary(fit0)$coefficients,f)
writeLines("Expenditure change for income groups",f)
write.csv(outreg(summary(fit1)$coefficients[,"Estimate"],summary(fit0)$coefficients[,"Std. Error"]),f)
# writeLines("Individual T-test",f)
# write.csv(tmp.tab,f)
close(f)

#--------------------# 
# Weekly data
# Expenditure change for average households
fit0 <- felm(BSKTDOL ~ recession + factor(MONTH)|PANID, data=exp_week)
fit1 <- felm(BSKTDOL ~ recession*INCOME1 + factor(MONTH)|PANID, data=exp_week)
fit2 <- lmer(BSKTDOL ~ factor(MONTH) + recession + (1 + recession|PANID),exp_week)
pan_char$week_Dexp_RE <- ranef(fit2)[[1]][,"recession1"] + fixef(fit2)['recession1']
pan_char$week_Dexp_dr <- with(pan_char,ifelse(week_Dexp_RE>0,"Exp increase","Exp decrease"))

# # Individual t-test 
# pan_char$week_Dexp		<- NA
# pan_char$week_Dexp_t 	<- NA
# pan_char$week_Dexp_TL 	<- NA
# pan_char$week_Dexp_TU 	<- NA
# for(i in 1:nrow(pan_char)){
# 	tmp1 <- subset(exp_week,PANID == pan_char[i,"PANID"] & recession == 0)[,"BSKTDOL"]
# 	tmp2 <- subset(exp_week,PANID == pan_char[i,"PANID"] & recession == 1)[,"BSKTDOL"]
# 	tmp  <- t.test(tmp2,tmp1)
# 	pan_char[i,"week_Dexp"] <- tmp$estimate[1] - tmp$estimate[2]
# 	pan_char[i,"week_Dexp_t"] <- tmp$statistic
# 	pan_char[i,c("week_Dexp_TL","week_Dexp_TU")] <- tmp$conf.int[1:2]
# 	if(i%%500){
# 		print(i)
# 	}
# }
# pan_char$week_Dexp_dr <- 0			# No change
# sel <- with(pan_char, week_Dexp_TL<0 & week_Dexp_TU<0 )
# pan_char[sel,"week_Dexp_dr"] <- -1
# sel <- with(pan_char, week_Dexp_TL>0 & week_Dexp_TU>0 )
# pan_char[sel,"week_Dexp_dr"] <- 1
# tmp.tab <- table(pan_char$week_Dexp_dr)/Nh
# tmp.tab <- rbind(tmp.tab,c(with(pan_char,by(week_Dexp,week_Dexp_dr,mean))))
# rownames(tmp.tab) <- c("Percentage","Change level")

f <- file(tmp.csv,"at")
writeLines("-----------------------------------------------",f)
writeLines("Change of expenditure with weekly data",f)
writeLines("Expenditure change for average households",f)
write.csv(summary(fit0)$coefficients[1,],f)
writeLines("Expenditure change for income groups",f)
write.csv(outreg(summary(fit0)$coefficients[,"Estimate"],summary(fit0)$coefficients[,"Std. Error"]),f)
# writeLines("Individual T-test",f)
# write.csv(tmp.tab,f)
close(f)

##########################
# Expenditure volatility # 
##########################
# Monthly data 
vol_month$ln_vol <- log(vol_month$vol)
vol_month <- merge(vol_month,pan_char[,c("PANID","month_Dexp_dr")],by=c("PANID"),all.x=T)
fit0 <- felm(vol ~ recession|PANID, data=vol_month)
fit1 <- felm(vol ~ INCOME1*recession |PANID, data=vol_month)
fit2 <- lmer(ln_vol ~ recession + (1 + recession | PANID), data=vol_month)
pan_char$month_Dvol_RE <- ranef(fit2)[[1]][,"recession1"] + fixef(fit2)['recession1']
pan_char$month_Dvol_dr <- with(pan_char,ifelse(month_Dvol_RE>0,"Volatility increase","Volatility decrease"))
fit3 <- felm(vol ~ month_Dexp_dr*recession |PANID, data=vol_month)
table(pan_char$month_Dvol_dr)/Nh
with(pan_char,by(month_Dvol_RE,month_Dvol_dr,mean))

# # Individual F test
# tmp <- felm(BSKTDOL ~ factor(MONTH) + recession*INCOME1 |PANID, data=exp_month)
# month.ef <- c(0,coef(tmp)[1:11])
# names(month.ef) <- 1:12
# tmp.exp <- exp_month
# tmp.exp$month.ef <- month.ef[as.character(tmp.exp$MONTH)]
# tmp.exp$BSKTDOL <- tmp.exp$BSKTDOL - tmp.exp$month.ef
# pan_char$month_vol_F	<- NA
# par_char$month_vol_FL	<- NA
# par_char$month_vol_FU	<- NA
# for(i in 1:nrow(pan_char)){
# 	tmp1 <- subset(tmp.exp,PANID == pan_char[i,"PANID"] & recession == 0)[,"BSKTDOL"]
# 	tmp2 <- subset(tmp.exp,PANID == pan_char[i,"PANID"] & recession == 1)[,"BSKTDOL"]
# 	tmp  <- var.test(tmp2,tmp1)
# 	pan_char[i,"month_vol_F"] <- tmp$statistic
# 	pan_char[i,c("month_vol_FL","month_vol_FU")] <- tmp$conf.int[1:2]
# 	if(i%%500==0){
# 		print(i)
# 	}
# }

f <- file(tmp.csv,"at")
writeLines("-----------------------------------------------",f)
writeLines("Change of expenditure volatility with monthly data",f)
writeLines("Volatility change for average households",f)
write.csv(summary(fit0)$coefficients,f)
writeLines("Volatility change by income group",f)
write.csv(outreg(summary(fit1)$coefficients[,"Estimate"],summary(fit1)$coefficients[,"Std. Error"]),f)
writeLines("Volatility change by expenditure change group",f)
write.csv(outreg(summary(fit3)$coefficients[,"Estimate"],summary(fit3)$coefficients[,"Std. Error"]),f)
close(f)

###############
# Basket type #
###############
bskt 	<- subset(exp_week,DOL > 0 )
selcol 	<- c(grep("DOL_",colnames(bskt)),grep("QUANT_",colnames(bskt)))
bskt 	<- bskt[,setdiff(1:ncol(bskt),selcol)]
tmp 	<- exp_month[,c("PANID","YMONTH","BSKTDOL")]; 
names(tmp) <- c("PANID","YMONTH","MONTH_BSKTDOL") 
tmp1	<- data.frame(vol_month)[,c("PANID","recession","vol")]
names(tmp1) <- c("PANID","recession","MONTH_VOL")
bskt 	<- merge(bskt,tmp,by=c("PANID","YMONTH"),all.x=T)
bskt	<- merge(bskt,tmp1,by=c("PANID","recession"),all.x=T)
bskt	<- merge(bskt,pan_char[,c("PANID","month_Dexp_dr","month_Dvol_dr")],by="PANID",all.x=T)
bskt	<- data.table(bskt)
setkey(bskt,PANID,WEEK)
bskt	<- bskt[,':='(r_quant=quant_noned/quant_ed,r_num=num_nonedible/num_edible, ln_QUANT = log(QUANT), 
					ln_num_cat=log(num_cat),ln_MONTH_BSKTDOL=log(MONTH_BSKTDOL), ln_MONTH_VOL=log(MONTH_VOL))]
bskt 	<- bskt[,':='(QUANT_LAG=c(NA,QUANT[-length(WEEK)]),num_cat_lag=c(NA,num_cat[-length(WEEK)]),
					  ln_QUANT_lag = c(NA,ln_QUANT[-length(WEEK)]), ln_num_cat_lag = c(NA,ln_num_cat[-length(WEEK)]),
					  quant_ed_lag=c(NA,quant_ed[-length(WEEK)]),quant_noned_lag=c(NA,quant_noned[-length(WEEK)]),
					r_quant_lag=c(NA,r_quant[-length(WEEK)]),r_num_lag=c(NA,r_num[-length(WEEK)])),by=list(PANID)]

# Basket quantity
fit0 <- felm(ln_QUANT ~ ln_MONTH_BSKTDOL + ln_MONTH_VOL + ln_QUANT_lag | PANID+MONTH, data = bskt)
fit1 <- felm(ln_QUANT ~ ln_MONTH_BSKTDOL:INCOME1 + ln_MONTH_VOL:INCOME1 + ln_QUANT_lag | PANID+MONTH, data = bskt)
fit2 <- felm(ln_QUANT ~ ln_MONTH_BSKTDOL:month_Dexp_dr + ln_MONTH_VOL:month_Dexp_dr + ln_QUANT_lag | PANID+MONTH, data = bskt)
fit3 <- felm(ln_QUANT ~ ln_MONTH_BSKTDOL:month_Dvol_dr + ln_MONTH_VOL:month_Dvol_dr + ln_QUANT_lag | PANID+MONTH, data = bskt)

f <- file(tmp.csv,"at")
writeLines("-----------------------------------------------",f)
writeLines("Regression of log(basket quantity)",f)
write.csv(outreg(summary(fit0)$coefficients[,"Estimate"],summary(fit0)$coefficients[,"Std. Error"]),f)
writeLines("Regression of log(basket quantity) by income group",f)
write.csv(outreg(summary(fit1)$coefficients[,"Estimate"],summary(fit0)$coefficients[,"Std. Error"]),f)
writeLines("Regression of log(basket quantity) by expenditure change group",f)
write.csv(outreg(summary(fit2)$coefficients[,"Estimate"],summary(fit0)$coefficients[,"Std. Error"]),f)
writeLines("Regression of log(basket quantity) by volatility change group",f)
write.csv(outreg(summary(fit3)$coefficients[,"Estimate"],summary(fit0)$coefficients[,"Std. Error"]),f)

# Basket variety
fit0 <- felm(ln_num_cat ~ ln_MONTH_BSKTDOL + ln_MONTH_VOL + ln_num_cat_lag | PANID+MONTH, data = bskt)
fit1 <- felm(ln_num_cat ~ ln_MONTH_BSKTDOL:INCOME1 + ln_MONTH_VOL:INCOME1 + ln_num_cat_lag | PANID+MONTH, data = bskt)
fit2 <- felm(ln_num_cat ~ ln_MONTH_BSKTDOL:month_Dexp_dr + ln_MONTH_VOL:month_Dexp_dr + ln_num_cat_lag | PANID+MONTH, data = bskt)
fit3 <- felm(ln_num_cat ~ ln_MONTH_BSKTDOL:month_Dvol_dr + ln_MONTH_VOL:month_Dvol_dr + ln_num_cat_lag | PANID+MONTH, data = bskt)

writeLines("-----------------------------------------------",f)
writeLines("Regression of log(number of category)",f)
write.csv(outreg(summary(fit0)$coefficients[,"Estimate"],summary(fit0)$coefficients[,"Std. Error"]),f)
writeLines("Regression of log(number of category) by income group",f)
write.csv(outreg(summary(fit1)$coefficients[,"Estimate"],summary(fit0)$coefficients[,"Std. Error"]),f)
writeLines("Regression of log(number of category) by expenditure change group",f)
write.csv(outreg(summary(fit2)$coefficients[,"Estimate"],summary(fit0)$coefficients[,"Std. Error"]),f)
writeLines("Regression of log(number of category) by volatility change group",f)
write.csv(outreg(summary(fit3)$coefficients[,"Estimate"],summary(fit0)$coefficients[,"Std. Error"]),f)
close(f)

#############################
# Monthly expenditure share #
#############################
# I drop groumet stores
fm.name <- c("Convenience","Conventional_supermarkets","Dollar_variety","Drug","Limited_assortment_stores","Mass",
             "Other","Super_stores","Wholesale_clubs")
fm.lab <- c("Convenience","Conventional supermarkets","Dollar","Drug","Limited assortment stores","Mass",
             "Other","Super stores","Wholesale clubs")
shr.fm 	<- paste("SHARE_",fm.name,sep="")

incd.df <- data.frame(NULL)
coef.df	<- data.frame(NULL)
for(i in 1:length(fm.name)){
	fit0 	<- felm(substitute(y ~ YMONTH | PANID, list(y = as.name(shr.fm[i]))),data = exp_month)
	coef.df <- rbind(coef.df,data.frame(summary(fit0)$coefficients,Format = fm.lab[i], Segment = NA))
	# fit1 	<- glmer(1*(SHARE_Mass>0) ~ YMONTH + (1|PANID), data=exp_month, family = binomial(link="logit"))
	# incd.df <- rbind(incd.df,data.frame(summary(fit0)$coefficients,Format = fm.lab[i], Segment = NA))
}

# Segment by INCOME
tmp <- levels(exp_month$INCOME1)
for(i in 1:length(fm.name)){
	for(j in 1:length(levels(exp_month$INCOME1))){
		fit0 <- felm(substitute(y ~ YMONTH | PANID, list(y = as.name(shr.fm[i]))),data = subset(exp_month,INCOME1==tmp[j]))
		coef.df <- rbind(coef.df,data.frame(summary(fit0)$coefficients,Format = fm.lab[i], Segment = paste(tmp[j],"-income",sep="")))
	}	
}

# Segment by expenditure change direction
tmp <- unique(pan_char$month_Dexp_dr)
for(i in 1:length(fm.name)){
	for(j in 1:length(tmp)){
		fit0 <- felm(substitute(y ~ YMONTH | PANID, list(y = as.name(shr.fm[i]))),
					data = subset(exp_month,PANID %in% pan_char[pan_char$month_Dexp_dr==tmp[j],"PANID"]))
		coef.df <- rbind(coef.df,data.frame(summary(fit0)$coefficients,Format = fm.lab[i], Segment = tmp[j]))
	}	
}

# Segment by expenditure volatility change 
tmp <- unique(pan_char$month_Dvol_dr)
for(i in 1:length(fm.name)){
	for(j in 1:length(tmp)){
		fit0 <- felm(substitute(y ~ YMONTH | PANID, list(y = as.name(shr.fm[i]))),
					data = subset(exp_month,PANID %in% pan_char[pan_char$month_Dvol_dr==tmp[j],"PANID"]))
		coef.df <- rbind(coef.df,data.frame(summary(fit0)$coefficients,Format = fm.lab[i], Segment = tmp[j]))
	}	
}

#-----------------------------------------#
# Plot monthly trend 
ggtmp <- coef.df
ggtmp$Time <- as.Date(paste("01",gsub("YMONTH","",rownames(ggtmp)),sep=""),format="%d%b%Y")

tmp1 <- c("Low-income","Med-income","High-income")
tmp2 <- c("Exp increase","Exp decrease")
tmp3 <- c("Volatility increase","Volatility decrease")

pdf(paste(my.wd,"/graph_share_trend.pdf",sep=""),width=12,height=8)
print(ggplot(subset(ggtmp,is.na(Segment)),aes(Time,Estimate)) + geom_point() + geom_line() + 
		facet_wrap(~Format) + ylab("Change of expenditure share") +
		geom_errorbar(aes(ymin=Estimate-1.96*Std..Error,ymax=Estimate+1.96*Std..Error) )
print(ggplot(subset(ggtmp,Segment %in% tmp1),aes(Time,Estimate,col=Segment)) + 
		geom_point() + geom_line() + 
		facet_wrap(~Format) + ylab("Change of expenditure share") + 
		geom_errorbar(aes(ymin=Estimate-1.96*Std..Error,ymax=Estimate+1.96*Std..Error)) )
print(ggplot(subset(ggtmp,Segment %in% tmp2),aes(Time,Estimate,col=Segment)) + 
		geom_point() + geom_line() + 
		facet_wrap(~Format) + ylab("Change of expenditure share") + 
		geom_errorbar(aes(ymin=Estimate-1.96*Std..Error,ymax=Estimate+1.96*Std..Error)))
print(ggplot(subset(ggtmp,Segment %in% tmp3),aes(Time,Estimate,col=Segment)) + 
		geom_point() + geom_line() + 
		facet_wrap(~Format) + ylab("Change of expenditure share") )
dev.off()







# This script estimates model. The moduels follow as: 
# 1. Segment households based on initial income level and subset data; 
# 2. Prepare estimation data; 
# 3. Estimate the conditional allocation model at different initial values; 
# 4. Simulate inclusive values; 
# 5. Estimate upper level expenditure model. 

library(ggplot2)
library(reshape2)
library(Rcpp)
library(RcppArmadillo)
library(maxLik)
library(evd)
library(data.table)
options(error = quote({dump.frames(to.file = TRUE)}))

args <- commandArgs(trailingOnly = TRUE)
print(args)
if(length(args)>0){
    for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
    }
}

# setwd("~/Documents/Research/Store switching/processed data")
# plot.wd	<- '~/Desktop'
# sourceCpp(paste("../Exercise/Multiple_discrete_continuous_model/", model_name, ".cpp", sep=""))
# source("../Exercise/Multiple_discrete_continuous_model/0_Allocation_function.R")

setwd("/home/brgordon/ccv103/Exercise/run")
# setwd("/kellogg/users/marketing/2661703/Exercise/run")
# setwd("/sscc/home/c/ccv103/Exercise/run")
model_name 	<- "MDCEV_share"
run_id		<- 1
# seg_id		<- 1
make_plot	<- TRUE
plot.wd		<- paste(getwd(), "/estrun_",run_id, sep="")

sourceCpp(paste(model_name, ".cpp", sep=""))
source("0_Allocation_function.R")

#################################
# Read data and subsetting data #
#################################
load("hh_biweek_exp.rdata")

# Extract 5% random sample
# length(unique(hh_exp$household_code))
# sel			<- sample(unique(hh_exp$household_code), .01*length(unique(hh_exp$household_code)) )
# hh_exp_save	<- hh_exp
# hh_exp		<- subset(hh_exp, household_code %in% sel)

# Retail formats 
fmt_name 	<- as.character(sort(unique(fmt_attr$channel_type)))
R			<- length(fmt_name)

# Convert some date and factor variables
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
hh_exp$first_incomeg	<- factor(hh_exp$first_incomeg, levels = c("T1", "T2", "T3"))
hh_exp$ln_inc<- log(hh_exp$income_midvalue)
cat("Table of initial income distribution:\n"); print(table(panelist$first_incomeg)); cat("\n")
cat("Table of segments in the expenditure data:\n"); print(table(hh_exp$first_incomeg)); cat("\n")

# Subset data 
selcol		<- c("household_code", "biweek", "dol", "year", "month", "income_midvalue", "ln_inc","first_incomeg", "scantrack_market_descr",paste("SHR_", gsub("\\s", "_", fmt_name), sep=""))
mydata		<- subset(hh_exp[,selcol], as.numeric(first_incomeg) == seg_id & dol > .1 )

################################
# Organize data for estimation #
################################
# The data required for estimation: 
# outcome variables: y (expenditure), shr (expenditure share);
# explanatory variables: price, X_list

# Format attributes
fmt_attr 		<- subset(fmt_attr, year > 2003)
fmt_attr$name 	<- paste(fmt_attr$scantrack_market_descr,fmt_attr$year, sep="-")
fmt_attr$ln_num_module 	<- log(fmt_attr$num_module)
fmt_attr$ln_upc_per_mod <- log(fmt_attr$avg_upc_per_mod)

# Match biweekly prices
price_dat 		<- subset(price_dat, year > 2003)
tmp		<- price_dat
tmp$name<- gsub("\\s", "_", tmp$channel_type)
tmp$name<- paste("PRC_", tmp$name, sep="")
price 	<- dcast(tmp, scantrack_market_descr + biweek ~ name, value.var = "bsk_price_paid_2004")
tmpidx	<- 2

# Impute the missing price for some regioins;
sum(is.na(price))
tmp			<- data.table(price_dat)
tmp			<- tmp[,list(price = mean(bsk_price_paid_2004, na.rm=T)), by=list(biweek, channel_type)]
tmp			<- data.frame(tmp)
sel			<- which(is.na(as.matrix(price[,(tmpidx+1):ncol(price)])), arr.ind=T)
dim(sel)
fmt_name[unique(sel[,2])]
for(i in 1:length(unique(sel[,2]))){
	sel1		<- unique(sel[,2])[i]
	sel2		<- sel[sel[,2]==sel1,1]
	tmp1		<- merge(subset(tmp, channel_type==fmt_name[sel1]), price[sel2,1:3], by=c("biweek"))
	tmp2		<- as.vector(tmp1$price)
	names(tmp2)	<- with(tmp1, paste(scantrack_market_descr, biweek,sep="-"))
	price[sel2, (sel1+tmpidx)] <- tmp2[with(price[sel2,], paste(scantrack_market_descr, biweek,sep="-"))]
}

# Merge price data to the trip data
mydata <- merge(mydata, price, by=c("scantrack_market_descr", "biweek"), all.x=T)
ord		<- order(mydata$household_code, mydata$biweek)
mydata	<- mydata[ord,]

# Outcome variables as matrix
sel		<- paste("SHR_", gsub("\\s", "_", fmt_name), sep="")
shr		<- as.matrix(mydata[,sel])
y		<- as.vector(mydata$dol)
ln_inc	<- mydata$ln_inc
sel		<- paste("PRC_", gsub("\\s", "_", fmt_name), sep="")
price	<- as.matrix(mydata[,sel])

# Select the index that has the positive expenditure
tmp 	<- price * 1 *(shr> 0)
s1_index<- apply(tmp, 1, which.max)

# Match retailers' attributes
tmp1	<- unique(fmt_attr$year)
tmp2	<- unique(fmt_attr$scantrack_market_descr)
tmp		<- paste(rep(tmp2, each=length(tmp1)), rep(tmp1, length(tmp2)), sep="-")
tmpn	<- 1:length(tmp)		
names(tmpn) <- tmp
sel 	<- with(mydata, paste(scantrack_market_descr,year, sep="-"))
sel1	<- tmpn[sel]
selcol	<- c("size_index", "ln_upc_per_mod", "ln_num_module","overall_prvt")
nx 		<- length(selcol) * 2

X_list 	<- vector("list", length=length(fmt_name))
for(i in 1:length(fmt_name)){
	sel2		<- fmt_attr$channel_type == fmt_name[i]
	tmp			<- fmt_attr[sel2,selcol]
	tmp1 		<- as.matrix(tmp[sel1,])
	X_list[[i]]	<- cbind(tmp1, tmp1 * ln_inc)
}

##############
# Estimation #
##############
# Set estimation control: fixing parameters
fixgamma 	<- FALSE
fixsigma 	<- TRUE
if(fixgamma){
	myfix 	<- (nx+R):(nx+2*R-1)
}else{
	myfix 	<- NULL
}

if(fixsigma){
	myfix	<- c(myfix, nx+2*R)
}
beta0_base 	<- which(fmt_name == "Grocery")

MDCEV_wrapper <- function(param){
	MDCEV_ll_fnC(param, nx, shr, y, s1_index, price, X_list, beta0_base)$ll
}

# Estimation with multiple initial values. 
theta_init	<- list(c(-1, -.1, 1, -2, .1, .1, -.1, .1, -1, -1, -1, -.5, -2, 3, 15, 4, 3, 6, 33, 0), 
					c(-1, -.5, 1.5, -.1, .1, -.1, -1, .5, -3, -1, -2, -1, -1, 3, 15, 4, 4, 8, 40, 0  ) )
system.time(tmp <- MDCEV_wrapper(theta_init[[1]]) )

tmp_sol <- vector("list", length(theta_init))
pct <- proc.time()
for(i in 1:length(theta_init)){
	names(theta_init[[i]]) 	<- c(paste("beta_",1:nx, sep=""), paste("beta0_",setdiff(1:R, beta0_base),sep=""), 
								 paste("gamma_",1:R,sep=""), "ln_sigma")
	tmp_sol[[i]]	<- maxLik(MDCEV_wrapper, start=theta_init[[i]], method="BFGS", fixed=myfix)
}
tmp		<- sapply(tmp_sol, function(x) ifelse(is.null(x$maximum), NA, x$maximum))
cat("The ll from these estimates are:\n"); print(tmp); cat("\n")
tmp1	<- sapply(tmp_sol, coef)
cat("The coefficient estimates are:\n"); print(tmp1); cat("\n")

# Choose the estimation that max ll
sel 	<- which(abs(tmp - max(tmp, na.rm=T)) < 1e-4, arr.ind=T)
sel1	<- sapply(tmp_sol[sel], function(x) !any(is.na(summary(x)$estimate[,"Std. error"])) & !any(summary(x)$estimate[,"Std. error"]==Inf) )
sel		<- ifelse(sum(sel1)==1, sel[sel1], ifelse(sum(sel1)==0, sel[1], sel[sel1][1]))
# sel 	<- which.max(sapply(tmp_sol, function(x) ifelse(is.null(x$maximum), NA, x$maximum)))
sol		<- tmp_sol[[sel]]
if(!sol$code %in% c(0,1,2,8)){
	cat("MDCEV does NOT have normal convergence.\n")
}
use.time <- proc.time() - pct
cat("MDCEV estimation finishes with", use.time[3]/60, "min.\n")
print(summary(sol))
cat("--------------------------------------------------------\n")

############################
# Simulate inclusive value #
############################
set.seed(666)
numsim 		<- 100
eps_draw	<- matrix(rgev(numsim*R), numsim, R)

# Simulate inclusive value omega(Inc, y|bar{X})
# For simulation, we need elements: income level, expenditure, (average) price, (average) retail attributes
lnInc	<- sort(unique(hh_exp$ln_inc))
numnodes<- 30		# Number of interpolation nodes of expenditure
tmpy	<- quantile(mydata$dol, c(0:(numnodes-1)/numnodes) )
tmpy	<- c(tmpy, max(mydata$dol)*1.2)
numnodes<- numnodes + 1

# Average price
tmp 		<- dcast(price_dat,	 scantrack_market_descr + year + biweek ~ channel_type, value.var = "bsk_price_paid_2004") 
tmp_price	<- rep(1,numnodes) %*% t(colMeans(as.matrix(tmp[,4:(3+R)]), na.rm=T) )

# Average retail attributes
tmpX_list	<- lapply(fmt_name, function(x) colMeans(as.matrix(subset(fmt_attr, channel_type == x)[,selcol])))		

#-------------------------------------------- # 
# Compute the inclusive value with random draws
tmp_coef	<- coef(sol)
omega_draw	<- array(NA, c(numsim, length(lnInc), numnodes), dimnames = list(NULL, lnInc, tmpy))
pct			<- proc.time()
for(i in 1:numsim){
	for(j in 1:length(lnInc)){
		tmpX_list1	<- lapply(tmpX_list, function(x) rep(1, numnodes) %*% t(c(x, x*lnInc[j])))
		tmpsol 		<- incl_value_fn(param_est=tmp_coef, base= beta0_base, X_list=tmpX_list1, y=tmpy, Q=Inf, price=tmp_price, 
					R=R, Ra=R, qz_cons = 0, exp_outside = FALSE, quant_outside = FALSE, eps_draw = rep(1, numnodes) %*% t(eps_draw[i,]) )
		omega_draw[i,j,] <- tmpsol$omega
	}	
}
use.time	<- proc.time() - pct
omega1		<- apply(omega_draw, c(2,3), mean, na.rm=T)

cat("Simulation of omega with random draws finishes with", use.time[3]/60, "min.\n")
cat("--------------------------------------------------------\n")

#---------------------------------------- #
# Simulate inclusive value omega(Inc,y,X) #
mkt.yr 	<- sort(unique(with(mydata, paste(scantrack_market_descr,year, sep="-"))))
cat("Dimention of market-year combination is", length(mkt.yr), ".\n")

# Prices 
tmp 	<- dcast(price_dat,	 scantrack_market_descr + year ~ channel_type, value.var = "bsk_price_paid_2004", fun.aggregate = mean) 
tmpn	<- setNames(1:nrow(tmp), with(tmp, paste(scantrack_market_descr,year, sep="-")) )
sel		<- tmpn[mkt.yr]
tmp_price1	<- as.matrix(tmp[sel,-c(1:2)])

# Match retailers' attributes
X1_list 	<- setNames(vector("list", length=length(fmt_name)), fmt_name)
for(i in 1:length(fmt_name)){
	tmp			<- subset(fmt_attr, channel_type == fmt_name[i])
	tmpn		<- setNames(1:nrow(tmp), with(tmp, paste(scantrack_market_descr,year, sep="-")) )
	sel			<- tmpn[mkt.yr]
	tmp1 		<- as.matrix(tmp[sel,selcol])
	X1_list[[i]]	<- tmp1
}

# Simulate
omega_full	<- array(NA, c(length(lnInc), numnodes, length(mkt.yr)), dimnames = list(lnInc, tmpy, mkt.yr))
for(i in 1:length(lnInc)){
	pct			<- proc.time()
	for(j in 1:numnodes){
		tmpX		<- lapply(X1_list, function(x) cbind(x, x*lnInc[i]))
		tmp.draw	<- matrix(NA, numsim, length(mkt.yr))
		for(k in 1:numsim){
			tmpsol 		<- incl_value_fn(param_est=tmp_coef, base= beta0_base, X_list=tmpX, y= rep(tmpy[j], length(mkt.yr)), Q=Inf, price=tmp_price1, 
						R=R, Ra=R, qz_cons = 0, exp_outside = FALSE, quant_outside = FALSE, eps_draw = rep(1, length(mkt.yr)) %*% t(eps_draw[k,]) )
			tmp.draw[k,]	<- tmpsol$omega
		}
		omega_full[i,j,]	<- colMeans(tmp.draw, na.rm = T)
	}
	use.time	<- proc.time() - pct
	cat("Full simulation at income level", i, "finishes using", use.time[3]/60,"min.\n")
}

####################################################
# Compare omega(Inc,y|bar{X}) and omega(Inc, y, X) #
####################################################
# Plot the distribution of inclusive across market-year
if(make_plot){
	pdf(paste(plot.wd,"/graph_inclusive_box_seg", seg_id, ".pdf", sep=""), width = 8, height = 8)
	ggtmp	<- setNames(melt(omega_full), c("Inc", "y", "mkty", "value"))
	ggtmp1	<- setNames(melt(omega1), c("Inc", "y", "value") )
	ggtmp	<- rbind(ggtmp, data.frame(ggtmp1[,1:2], mkty = NA, value = ggtmp1[,3]))
	ggtmp$Inc	<- round(exp(ggtmp$Inc))
	tmp		<- unique(ggtmp$Inc)
	for(i in 1:length(lnInc)){
		print(
			ggplot(subset(ggtmp, Inc == tmp[i])) + 
				geom_point(data = subset(ggtmp, is.na(mkty) & Inc == tmp[i]), aes(x= y, y = value), col = 2) + 
				geom_boxplot(data = subset(ggtmp, !is.na(mkty) & Inc == tmp[i]), aes(x = y, y = value, group = y)) + 
				scale_x_log10() + 
				labs(title = paste("Income leve = ", tmp[i], sep=""))
			)
	}
	dev.off()
}

save.image(paste("estrun_",run_id,"/test_inclusive_seg",seg_id,".rdata",sep=""))

cat("This program is done. ")


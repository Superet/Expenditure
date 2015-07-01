library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(Rcpp)
library(RcppArmadillo)
library(maxLik)
library(evd)
library(data.table)

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
# source("../Exercise/Multiple_discrete_continuous_model/0_marginal_effect_function.R")

setwd("/home/brgordon/ccv103/Exercise/run")
# setwd("/kellogg/users/marketing/2661703/Exercise/run")
# setwd("/sscc/home/c/ccv103/Exercise/run")
model_name 	<- "MDCEV_a1b1"
# run_id		<- "sub"
# seg_id		<- 1
make_plot	<- FALSE

sourceCpp(paste(model_name, ".cpp", sep=""))
source("0_Allocation_function.R")
source("0_marginal_effect_function.R")

##########################
# Read data and checking # 
##########################
if(run_id=="sub"){
	hh_exp		<- read.csv("2_hh_biweek_exp_merge_sub.csv")
}else{
	hh_exp		<- read.csv("2_hh_biweek_exp_merge_larger.csv")
}
fmt_attr 	<- read.csv("1_format_year_attr.csv")
price_dat	<- read.csv("1_format_biweek_price.csv")

# Segment households 
hh_exp$income_group <- factor(hh_exp$income_group, levels=paste("Qt",1:4,sep=""))
hh_exp$famsize		<- factor(hh_exp$famsize, levels = c("Single","Two", "Three+"))
hh_exp$Q			<- hh_exp$food_quant + hh_exp$nonedible_quant
mycut 				<- c(0,2,3,4,5,14) + .5
hh_exp$d			<- as.numeric(cut(hh_exp$num_day, mycut, labels=1:5))

panelist	<- data.table(hh_exp)
setkeyv(panelist, c("household_code","year","biweek"))
panelist	<- panelist[,list(first_income = income_group[1], first_famsize = famsize[1]), by=list(household_code)]
table(panelist$first_income, panelist$first_famsize)
seg_index	<- 1:(length(levels(panelist$first_income))*length(levels(panelist$first_famsize)))
names(seg_index) <- paste(rep(levels(panelist$first_income), length(levels(panelist$first_famsize))), 
						  rep(levels(panelist$first_famsize), each=length(levels(panelist$first_income)) ), sep="-")
tmp			<- paste(panelist$first_income, panelist$first_famsize, sep="-")
panelist$segment <- seg_index[tmp]
summary(panelist)
# Drop the household with NA demo
panelist	<- subset(panelist, !is.na(first_income))

hh_exp		<- merge(hh_exp, data.frame(panelist)[,c("household_code","segment")], by = "household_code")
sel			<- panelist$segment == seg_id
mydata		<- subset(hh_exp, household_code %in% panelist[sel,household_code])

# Exclude obs with 0 DOL
fmt_name 	<- as.character(sort(unique(fmt_attr$channel_type)))
sum(mydata$net_dol==0)
sum(mydata$DOL<=0)
sum(is.na(mydata$dol_purchases))
mydata 		<- subset(mydata, DOL>0)
mydata 		<- subset(mydata, net_dol>0)
mydata		<- subset(mydata, !is.na(dol_purchases) & dol_purchases > 0)
selc 		<- c("food_quant","nonedible_quant", "num_food_module","num_noned_module",
				"num_storable_module","num_nonstr_module", "storable_quant", "nonstr_quant")
sel1 		<- is.na(mydata[,selc[1]])
mydata[sel1,selc] <- 0

quantile(mydata$DOL, c(.9, .99, .999, .9999))
sel <- mydata$DOL > 3000
length(unique(mydata[sel,"household_code"]))
selhh <- unique(mydata[sel,"household_code"])
mydata <- subset(mydata, !household_code %in% selhh)

# Format attributes
fmt_attr 		<- subset(fmt_attr, year > 2003)
fmt_attr$name 	<- paste(fmt_attr$scantrack_market_descr,fmt_attr$year, sep="-")
fmt_attr$ln_num_module <- log(fmt_attr$num_module)
fmt_attr$ln_upc_per_mod <- log(fmt_attr$avg_upc_per_mod)
price_dat 		<- subset(price_dat, year > 2003)

# Check within-household change of shopping strategries. 
tmp			<- data.table(mydata)
tmp			<- tmp[,list(N=length(biweek)), by=list(household_code, d)]
tmp1		<- tmp[,list(times = max(N)/sum(N), freq_HHI = (N/sum(N))^2), by=list(household_code)]
summary(tmp1)

# Match biweekly prices
tmp		<- price_dat
tmp$name<- gsub("\\s", "_", tmp$channel_type)
tmp$name<- paste("PRC_", tmp$name, sep="")
price 	<- dcast(tmp, scantrack_market_descr + year + biweek ~ name, value.var = "bsk_price_paid_2004")

# Impute the missing price for some regioins;
tmp			<- data.table(price_dat)
tmp			<- tmp[,list(price = mean(bsk_price_paid_2004, na.rm=T)), by=list(year, biweek, channel_type)]
sel			<- which(is.na(as.matrix(price[,4:ncol(price)])), arr.ind=T)
dim(sel); 
fmt_name[unique(sel[,2])]
for(i in 1:length(unique(sel[,2]))){
	sel1		<- unique(sel[,2])[i]
	sel2		<- sel[sel[,2]==sel1,1]
	tmp1		<- merge(subset(tmp, channel_type==fmt_name[i]), price[sel2,1:3], by=c("year","biweek"))
	tmp2		<- as.vector(tmp1$price)
	names(tmp2)	<- with(tmp1, paste(scantrack_market_descr, year, biweek,sep="-"))
	price[sel2, (sel1+3)] <- tmp2[with(price[sel2,], paste(scantrack_market_descr, year, biweek,sep="-"))]
}

mydata <- merge(mydata, price, by=c("scantrack_market_descr", "year", "biweek"), all.x=T)

##########################
# Organize modeling data #
##########################
D				<- length(unique(mydata$d))
R				<- length(fmt_name)
all_e1_index 	<- vector("list", D)
all_eR			<- vector("list", D)
all_price		<- vector("list", D)
all_X_list		<- vector("list", D)
all_Xint_list	<- vector("list", D)
selcol			<- c("size_index", "ln_upc_per_mod", "ln_num_module","overall_prvt")
nx				<- length(selcol)

for(j in 1:D){
	my_data1<- subset(mydata, d==j & !is.na(dol_purchases))
	ord 	<- with(my_data1, order(scantrack_market_descr,household_code, year, biweek))
	my_data1<- my_data1[ord,]
	sel		<- paste("DOLP_", gsub("\\s", "_", fmt_name), sep="")
	eR		<- as.matrix(my_data1[,sel])
	sel		<- paste("PRC_", gsub("\\s", "_", fmt_name), sep="")
	price	<- as.matrix(my_data1[,sel])
	Q		<- as.vector(my_data1$storable_quant + my_data1$nonstr_quant) 
	
	# Select the index that has the positive expenditure
	e1_index <- rep(NA, nrow(eR))
	for(i in 1:nrow(eR)){
		sel	<- eR[i,]>0
		m	<- max(price[i,sel])
		e1_index[i] <- which(price[i,] == m)
	}

	# Match retailers' attributes
	tmp1	<- unique(fmt_attr$year)
	tmp2	<- unique(fmt_attr$scantrack_market_descr)
	tmp		<- paste(rep(tmp2, each=length(tmp1)), rep(tmp1, length(tmp2)), sep="-")
	tmpn	<- 1:length(tmp)		
	names(tmpn) <- tmp
	
	# Principal component analysis
	sel 	<- with(my_data1, paste(scantrack_market_descr,year, sep="-"))
	sel1	<- tmpn[sel]

	X_list 	<- vector("list", length=length(fmt_name))
	X_list1 <- vector("list", length=length(fmt_name))
	for(i in 1:length(fmt_name)){
		sel2		<- fmt_attr$channel_type == fmt_name[i]
		tmp			<- fmt_attr[sel2,selcol]
		X_list[[i]] <- as.matrix(tmp[sel1,])
		X_list1[[i]]<- cbind(X_list[[i]], X_list[[i]] * my_data1$recession)
	}
	
	all_eR[[j]]		<- eR
	all_e1_index[[j]]<- e1_index
	all_price[[j]]	<- price
	all_X_list[[j]]	<- X_list
	all_Xint_list[[j]]	<- X_list1
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

sol_list	<- vector("list", D)
for(k in 1:D){
	MDCEV_wrapper <- function(param){
		MDCEV_LogLike_fnC(param, nx, all_e1_index[[k]], all_eR[[k]], all_price[[k]], all_X_list[[k]], beta0_base)
	}
	
	# Estimation with multiple initial values. 
	theta_init	<- list(c(.2, .5, .5,-.5, -5, -1, -3, -3, -2, .08*c(1:R), 0), 
						c(-.2, 0, 0,-.1, -3, -1, -2, -2, -1, 1/c(1:R), 0),
						c(.7, -1, -1,-2, -2, -.5, -1, -2, -2, .5*c(1:R), 0),
						c(-3, -2, -2,-2, -10, -2, -2, -3, -2, 1/c(1:R), 0) )
	system.time(tmp <- MDCEV_wrapper(theta_init[[1]]) )
	
	tmp_sol <- vector("list", length(theta_init))
	pct <- proc.time()
	for(i in 1:length(theta_init)){
		names(theta_init[[i]]) 	<- c(paste("beta_",1:nx, sep=""), paste("beta0_",setdiff(1:R, beta0_base),sep=""), 
									 paste("gamma_",1:R,sep=""), "ln_sigma")
		tmp_sol[[i]]	<- maxLik(MDCEV_wrapper, start=theta_init[[i]], method="BFGS", fixed=myfix)
	}
	tmp		<- sapply(tmp_sol, function(x) ifelse(is.null(x$maximum), NA, x$maximum))
	sel 	<- which(abs(tmp - max(tmp, na.rm=T)) < 1e-4, arr.ind=T)
	sel1	<- sapply(tmp_sol[sel], function(x) !any(is.na(summary(x)$estimate[,"Std. error"])) & !any(summary(x)$estimate[,"Std. error"]==Inf) )
	sel		<- ifelse(sum(sel1)==1, sel[sel1], ifelse(sum(sel1)==0, sel[1], sel[sel1][1]))
	# sel 	<- which.max(sapply(tmp_sol, function(x) ifelse(is.null(x$maximum), NA, x$maximum)))
	sol		<- tmp_sol[[sel]]
	if(!sol$code %in% c(0,1,2,8)){
		cat("Shopping strategy",k, "does NOT have normal convergence.\n")
	}
	use.time <- proc.time() - pct
	cat("Shopping strategy", k, "with time", use.time[3], "sec. \n")
	sol_list[[k]] <- sol
}

# Test if preference changes over time 
cat("------------------------------------------------------------\n")
cat("Now estimate the model with interaction terms.\n")
myfix 		<- myfix + nx
sol_list1	<- vector("list", D)
for(k in 1:D){
	MDCEV_wrapper <- function(param){
		MDCEV_LogLike_fnC(param, nx*2, all_e1_index[[k]], all_eR[[k]], all_price[[k]], all_Xint_list[[k]], beta0_base)
	}
	
	# Estimation with multiple initial values. 
	theta_init	<- list(c(.2, .5, .5,-.5, rep(0,4), -5, -1, -3, -3, -2, .8*c(1:R), 0), 
						c(-.2, 0, 0, -.1, rep(.1,4), -3, -1, -2, -2, -1, 10/c(1:R), 0),
						c(.7, -1, -1, -2, rep(-.1,4), -2, -.5, -1, -2, -2, 2*c(1:R), 0),
						c(-3, -2, -2, -2, rep(.1, 4), -10, -2, -2, -3, -2, 20/c(1:R), 0) )
	system.time(tmp <- MDCEV_wrapper(theta_init[[1]]) )
	
	tmp_sol <- vector("list", length(theta_init))
	pct <- proc.time()
	for(i in 1:length(theta_init)){
		names(theta_init[[i]]) 	<- c(paste("beta_",1:nx, sep=""), paste("beta_",1:nx,"*Rc",sep=""),paste("beta0_",setdiff(1:R, beta0_base),sep=""), 
									 paste("gamma_",1:R,sep=""), "ln_sigma")
		tmp_sol[[i]]	<- maxLik(MDCEV_wrapper, start=theta_init[[i]], method="BFGS", fixed=myfix)
	}
	tmp		<- sapply(tmp_sol, function(x) ifelse(is.null(x$maximum), NA, x$maximum))
	sel 	<- which(abs(tmp - max(tmp, na.rm=T)) < 1e-4, arr.ind=T)
	sel1	<- sapply(tmp_sol[sel], function(x) !any(is.na(summary(x)$estimate[,"Std. error"])) & !any(summary(x)$estimate[,"Std. error"]==Inf) )
	sel		<- ifelse(sum(sel1)==1, sel[sel1], ifelse(sum(sel1)==0, sel[1], sel[sel1][1]))
	# sel 	<- which.max(sapply(tmp_sol, function(x) ifelse(is.null(x$maximum), NA, x$maximum)))
	sol		<- tmp_sol[[sel]]
	if(!sol$code %in% c(0,1,2,8)){
		cat("Shopping strategy",k, "does NOT have normal convergence.\n")
	}
	use.time <- proc.time() - pct
	cat("Shopping strategy", k, "with time", use.time[3], "sec. \n")
	cat("The estimates are:\n"); print(summary(sol)); cat("\n")
	sol_list1[[k]] <- sol
}

################################
# Simulate the inclusive value # 
################################
#----------------------------------------------------------------- # 
# Compute the inclusive value with average price and no random draws
numnodes<- 30		# Number of interpolation nodes 
tmpy	<- quantile(mydata$dol_purchases, c(0:(numnodes-1)/numnodes) )
tmpy	<- c(tmpy, max(mydata$dol_purchases)*1.2)
numnodes<- numnodes + 1

tmpX_list	<- lapply(fmt_name, function(x) 
				rep(1,numnodes) %*% t(colMeans(as.matrix(subset(fmt_attr, channel_type == x)[,selcol]))) )		
tmp 		<- dcast(price_dat,	 scantrack_market_descr + year + biweek ~ channel_type, value.var = "bsk_price_paid_2004") 
tmp_price	<- rep(1,numnodes) %*% t(colMeans(as.matrix(tmp[,4:(3+R)]), na.rm=T) )
omega_nodraw <- matrix(NA, D, numnodes, dimnames = list(d = 1:D, y = tmpy))

for(i in 1:D){
	tmp_coef	<- as.vector(coef(sol_list[[i]]))
	sol 		<- incl_value_fn(param_est=tmp_coef, base= beta0_base, X_list=tmpX_list, y=tmpy, Q=Inf, price=tmp_price, 
					R=R, Ra=R, qz_cons = 0, exp_outside = FALSE, quant_outside = FALSE, eps_draw = matrix(0, numnodes,R))
	omega_nodraw[i,] <- sol$omega
}

#----------------------------------------------------------------- # 
# Compute the inclusive value with average price and random draws 
set.seed(666)
numsim 	<- 200		# Number of simulation draws
numnodes<- 30		# Number of interpolation nodes 
tmpy	<- quantile(mydata$dol_purchases, c(0:(numnodes-1)/numnodes) )

tmpX_list	<- lapply(fmt_name, function(x) 
				rep(1,numnodes) %*% t(colMeans(as.matrix(subset(fmt_attr, channel_type == x)[,selcol]))) )		
tmp 		<- dcast(price_dat,	 scantrack_market_descr + year + biweek ~ channel_type, value.var = "bsk_price_paid_2004") 
tmp_price	<- rep(1,numnodes) %*% t(colMeans(as.matrix(tmp[,4:(3+R)]), na.rm=T) )
omega_draw <- array(NA, c(numsim, D, numnodes))
for(i in 1:numsim){
	pct	<- proc.time()
	for(d in 1:D){
		tmp_coef	<- as.vector(coef(sol_list[[d]]))
		sol 		<- incl_value_fn(param_est=tmp_coef, base= beta0_base, X_list=tmpX_list, y=tmpy, Q=Inf, price=tmp_price, 
						R=R, Ra=R, qz_cons = 0, exp_outside = FALSE, quant_outside = FALSE)
		omega_draw[i,d,] <- sol$omega
	}
	use.time <- proc.time() - pct
	cat("Simulation", i, ", using", use.time[3]/60,"min.\n")
}

mytrim			<- .2
tmp 			<- apply(omega_draw, c(2,3), mean, na.rm=T, trim = mytrim)
dimnames(tmp)	<- list(d = 1:D, y = tmpy)
tmp				<- melt(tmp)
omega_draw_mean	<- data.frame(tmp[,1:2], omega = tmp$value)

# Plot omega function
# -------------------------------------------# 
# Inclusive value of purchase utility: omega(y,d)
omega_dat 	<- omega_draw_mean
omega_list	<- lapply(1:D, function(i) 
			splinefun(omega_dat[omega_dat$d==i,"y"], omega_dat[omega_dat$d==i,"omega"], method="natural"))
omega_fn 	<- function(y,d){
	if(d==0){
		out	<- rep(0, length(y)) 
	}else{
		out <- omega_list[[d]](y) 
	}
	return(as.vector(out))
}

tmp 		<- seq(10, 2000, by=10)
ggtmp		<- sapply(1:D, function(i) omega_fn(tmp,i))
dimnames(ggtmp) <- list(y=tmp, d=1:D)
ggtmp		<- melt(ggtmp)

pdf(paste("run_",run_id,"/graph_omega_seg",seg_id,".pdf",sep=""), width = 6.5, height = 6.5)
ggplot(ggtmp, aes(y, value, col=factor(d))) + geom_point() + geom_line() 
dev.off()

####################
# Save the results #
####################
ls()
rm(list=c("mydata", "my_data1", "Allocation_constr_fn","Allocation_fn","Allocation_nlop_fn","incl_value_fn","i","j","k","d","MDCEV_grad_fnC","MDCEV_ll_fnC",
		  "MDCEV_LogLike_fnC","MDCEV_wrapper","tmp","tmp1","tmp2","tmp_sol","sel","sel1","selcol","param_assign"))

save.image(paste("run_",run_id,"/5_est_MDCEV_seg",seg_id,".rdata",sep=""))

rm(list=c("run_id","seg_id","vid"))

cat("This program is done. ")

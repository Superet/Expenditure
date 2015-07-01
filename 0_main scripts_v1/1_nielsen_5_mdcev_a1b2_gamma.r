library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(Rcpp)
library(RcppArmadillo)
library(maxLik)
library(evd)
library(data.table)

# setwd("~/Documents/Research/Store switching/processed data")
# plot.wd	<- '~/Desktop'
# 
# setwd("//tsclient/Resear1/Store switching/processed data")
# plot.wd <- "//tsclient/Resear1/Store switching/processed data/temp_results"

setwd("/home/brgordon/ccv103/processed data")

sourceCpp("../Exercise/Multiple_discrete_continuous_model/MDCEV_a1b2_gamma.cpp")
source("../Exercise/Multiple_discrete_continuous_model/Allocation_function.R")
source("../Exercise/Multiple_discrete_continuous_model/marginal_effect_function.R")
model_name 	<- "MDCEV_a1b2"
run_id		<- 1

#############
# Functions #
#############
heat_map <- function(x, y, x.cut, y.cut, breaks=NULL, print=F){
	if(!is.null(breaks)){
		x.cut <- quantile(x, c(0:breaks)/breaks, na.rm=T)
		y.cut <- quantile(y, c(0:breaks)/breaks, na.rm=T)
	}
	x1 <- cut(x, x.cut, include.lowest=T)
	y1 <- cut(y, y.cut, include.lowest=T)
	tmp <- table(x1, y1)/length(x)*100
	ggtmp <- melt(tmp)
	if(print){
		p <- ggplot(ggtmp, aes(x1, y1, fill=value)) + geom_tile() + 
				scale_fill_gradientn(colours = myPalette(100)) + 
				theme_bw()
		print(p)
	}		
	return(ggtmp)
}
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")

##########################
# Read data and checking # 
##########################
hh_exp		<- read.csv("2_hh_month_exp_merge_sub.csv")
fmt_attr 	<- read.csv("1_format_year_attr.csv")
price_dat	<- read.csv("1_format_month_price.csv")

# Exclude obs with 0 DOL
fmt_name 	<- as.character(unique(fmt_attr$channel_type))
sum(hh_exp$net_dol==0)
sum(hh_exp$DOL<=0)
hh_exp 		<- subset(hh_exp, DOL>0)
hh_exp 		<- subset(hh_exp, net_dol>0)
selc 		<- c("food_quant","nonedible_quant", "num_food_module","num_noned_module",
				"num_storable_module","num_nonstr_module", "storable_quant", "nonstr_quant")
sel1 		<- is.na(hh_exp[,selc[1]])
hh_exp[sel1,selc] <- 0

quantile(hh_exp$DOL, c(.9, .99, .999, .9999))
sel <- hh_exp$DOL > 3000
length(unique(hh_exp[sel,"household_code"]))
selhh <- unique(hh_exp[sel,"household_code"])
hh_exp <- subset(hh_exp, !household_code %in% selhh)

# Format attributes
fmt_attr 	<- subset(fmt_attr, year > 2003)
fmt_attr$name <- paste(fmt_attr$scantrack_market_descr,fmt_attr$year, sep="-")
fmt_attr$ln_num_module <- log(fmt_attr$num_module)
fmt_attr$ln_upc_per_mod <- log(fmt_attr$avg_upc_per_mod)

# Basket type
selc <- c("food_quant","nonedible_quant", "num_food_module","num_noned_module",
			"num_storable_module","num_nonstr_module", "storable_quant", "nonstr_quant")
hh_exp$quant <- with(hh_exp, food_quant + nonedible_quant)

# Check distribution of basket quantiy
# ggtmp <- heat_map(hh_exp$storable_quant, hh_exp$nonstr_quant, x.cut, y.cut, breaks=3, T)

# Classify shopping strategies
my.break	<- 3
my_alpha	<- .05
my_scale	<- 100			# Scale up size index 
tmp 	 	<- cut(hh_exp$quant, quantile(hh_exp$quant, c(0:my.break)/my.break, na.rm=T), include.lowest=T)
tmp1 	 	<- cut(hh_exp$num_trip, quantile(hh_exp$num_trip, c(0:my.break)/my.break, na.rm=T), include.lowest=T)
table(tmp, tmp1)
bskt_type	<- 1:my.break^2
names(bskt_type) <- paste(rep(levels(tmp), my.break), rep(levels(tmp1), each=my.break), sep="-")
bskt_typeQ	<- my_scale * rep(quantile(hh_exp$quant, seq(1, 2*my.break - 1, by=2)/(my.break*2)), my.break)
bskt_typeTR	<- rep(quantile(hh_exp$num_trip, seq(1, 2*my.break - 1, by=2)/(my.break*2)), each=my.break)
names(bskt_typeQ) <- bskt_type
names(bskt_typeTR)<- bskt_type
tmp2 		<- paste(tmp, tmp1, sep="-")
hh_exp$basket_type	<- bskt_type[tmp2]

# Check within-household change of shopping strategries. 
tmp			<- data.table(hh_exp)
tmp			<- tmp[,list(N=length(month)), by=list(household_code, basket_type)]
tmp			<- tmp[,':='(Q = bskt_typeQ[tmp$basket_type], TR = bskt_typeTR[tmp$basket_type])]
tmp1		<- tmp[,list(times = max(N)/sum(N), freq_HHI = (N/sum(N))^2), by=list(household_code)]
summary(tmp1)
tmp1		<- tmp[,list(N = sum(N)), by=list(household_code, Q)]
tmp1		<- tmp1[,list(times = max(N)/sum(N), freq_HHI = (N/sum(N))^2), by=list(household_code)]
summary(tmp1)
tmp1		<- tmp[,list(N = sum(N)), by=list(household_code, TR)]
tmp1		<- tmp1[,list(times = max(N)/sum(N), freq_HHI = (N/sum(N))^2), by=list(household_code)]
summary(tmp1)

# Match monthly prices
tmp		<- price_dat
tmp$name<- gsub("\\s", "_", tmp$channel_type)
tmp$name<- paste("PRC_", tmp$name, sep="")
price 	<- dcast(tmp, scantrack_market_descr + year + month ~ name, value.var = "price_paid_norm_index")
my_data <- merge(hh_exp, price, by=c("scantrack_market_descr", "year", "month"), all.x=T)

# Check the the discrepency between sum of e/p and Q
q 		<- with(my_data, food_quant + nonedible_quant )
sel		<- paste("DOLP_", gsub("\\s", "_", fmt_name), sep="")
sel1	<- paste("PRC_", gsub("\\s", "_", fmt_name), sep="")
q1		<- rowSums( my_data[,sel]/my_data[,sel1] )
# hist(q - q1, breaks=100)
quantile(q-q1, c(.01, .99))

q		<- bskt_typeQ[my_data$basket_type]
# hist(q - q1, breaks=100)
quantile(q-q1, c(.01, .99))
delta_q <- min(q - q1) - .001

##########################
# Organize modeling data #
##########################
K				<- length(bskt_type)
all_e1_index 	<- vector("list", length(bskt_type))
all_eR			<- vector("list", length(bskt_type))
all_price		<- vector("list", length(bskt_type))
all_X_list		<- vector("list", length(bskt_type))
all_Q			<- vector("list", length(bskt_type))

for(j in 1:length(bskt_type)){
	bskt_idx<- j
	R		<- length(fmt_name)
	my_data1<- subset(my_data, basket_type == bskt_idx & !is.na(lag_dol))
	ord 	<- with(my_data1, order(scantrack_market_descr,household_code, year, month))
	my_data1 <- my_data1[ord,]
	sel		<- paste("DOLP_", gsub("\\s", "_", fmt_name), sep="")
	eR		<- as.matrix(my_data1[,sel])
	sel		<- paste("PRC_", gsub("\\s", "_", fmt_name), sep="")
	price	<- as.matrix(my_data1[,sel])
	# Q		<- as.vector(my_data1$storable_quant + my_data1$nonstr_quant) 
	Q		<- as.vector(bskt_typeQ[my_data1$basket_type])
	
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

	selcol	<- c("size_index", "ln_num_module", "ln_upc_per_mod", "overall_prvt")
	nx		<- length(selcol)
	sel 	<- with(my_data1, paste(scantrack_market_descr,year, sep="-"))
	X_list 	<- vector("list", length=length(fmt_name))
	for(i in 1:length(fmt_name)){
		tmp <- subset(fmt_attr, channel_type == fmt_name[i])
		sel1<- tmpn[sel]
		X_list[[i]] <- as.matrix(tmp[sel1, selcol])
	}

	
	all_eR[[j]]		<- eR
	all_e1_index[[j]]<- e1_index
	all_price[[j]]	<- price
	all_X_list[[j]]	<- X_list
	all_Q[[j]]		<- Q
}

##############
# Estimation #
##############
sol_list <- vector("list", length(unique(my_data$basket_type)))
for(k in 1:length(bskt_type)){
	MDCEV_wrapper <- function(param){
		MDCEV_LogLike_fnC(param, nx, c_q=-delta_q, Q=all_Q[[k]], e=all_eR[[k]], e1_index = all_e1_index[[k]], 
			p=all_price[[k]], X_list = all_X_list[[k]])
	}
	
	# Estimation with multiple initial values. 
	theta_init	<- list(c(.1, -.5, 1, -.5, 1, -1, 1, 1, -1, 1, 300, 20, 50, 10, 3, 10, 150), 
						c(-.2, 0, .5, -.1, 1, -1, 1, 1, -2, 1, 300, 10, 50, 20, 3, 5, 100),
						c(-3,  .2,  .5, -2, 1, -.5, 1, 2, -2, 1, 200, 20, 50, 50, 10, 10, 80),
						c(.5, .5, 1, -2, 5, -5, 5, 5, -5, 5, 200, 20, 5, 40, 5, 6, 100)  )
	system.time(tmp <- MDCEV_wrapper(theta_init[[1]]) )
	
	tmp_sol <- vector("list", length(theta_init))
	pct <- proc.time()
	for(i in 1:length(theta_init)){
		names(theta_init[[i]]) <- c(paste("beta_",1:nx, sep=""), paste("beta0_",2:R,sep=""),paste("gamma_", 1:R, sep=""))	
		tmp_sol[[i]]	<- maxLik(MDCEV_wrapper, start=theta_init[[i]], method="BFGS")
	}
	tmp		<- sapply(tmp_sol, function(x) ifelse(is.null(x$maximum), NA, x$maximum))
	sel 	<- which(abs(tmp - max(tmp, na.rm=T)) < 1e-6, arr.ind=T)
	sel1	<- sapply(tmp_sol[sel], function(x) !any(is.na(summary(x)$estimate[,"Std. error"])))
	sel		<- ifelse(sum(sel1)==1, sel[sel1], ifelse(sum(sel1)==0, sel[1], sel[sel1][1]))
	sel 	<- which.max(sapply(tmp_sol, function(x) ifelse(is.null(x$maximum), NA, x$maximum)))
	sol		<- tmp_sol[[sel]]
	if(!sol$code %in% c(0,1,2,8)){
		cat("Basket type",k, "does NOT have normal convergence.\n")
	}
	use.time <- proc.time() - pct
	cat("Basket type", k, "with time", use.time[3], "sec. \n")
	sol_list[[k]] <- sol
}

################################
# Simulate the inclusive value # 
################################
#----------------------------------------------------------------- # 
# Compute the inclusive value with average price and no random draws
numnodes<- 20		# Number of interpolation nodes 
K		<- length(bskt_type)
tmpy	<- quantile(hh_exp$dol_purchases, c(0:(numnodes-1)/numnodes) )

selcol	<- c("size_index", "ln_num_module", "ln_upc_per_mod", "overall_prvt")
tmpX_list	<- lapply(fmt_name, function(x) 
				rep(1,numnodes) %*% t(colMeans(as.matrix(subset(fmt_attr, channel_type == x)[,selcol]))) )		
tmp 		<- dcast(price_dat,	 scantrack_market_descr + year + month ~ channel_type, value.var = "price_paid_norm_index") 
tmp_price	<- rep(1,numnodes) %*% t(colMeans(as.matrix(tmp[,4:(3+R)]), na.rm=T) )
omega_nodraw <- matrix(NA, K, numnodes, dimnames = list(bskt_type = 1:K, y = tmpy))

for(k in 1:K){
	tmp_coef	<- as.vector(coef(sol_list[[k]]))
	sol 		<- incl_value_fn(param_est=tmp_coef, X_list=tmpX_list, y=tmpy, Q=bskt_typeQ[k], price=tmp_price, 
					R=R, Ra=R+1, qz_cons = -delta_q, exp_outside = FALSE, quant_outside = TRUE, eps_draw = 0)
	omega_nodraw[k,] <- sol$omega
}

#----------------------------------------------------------------- # 
# Compute the inclusive value with average price and random draws 
numsim 	<- 200		# Number of simulation draws
numnodes<- 20		# Number of interpolation nodes 
K		<- length(bskt_type)
tmpy	<- quantile(hh_exp$dol_purchases, c(0:(numnodes-1)/numnodes) )

selcol	<- c("size_index", "ln_num_module", "ln_upc_per_mod", "overall_prvt")
tmpX_list	<- lapply(fmt_name, function(x) 
				rep(1,numnodes) %*% t(colMeans(as.matrix(subset(fmt_attr, channel_type == x)[,selcol]))) )		
tmp 		<- dcast(price_dat,	 scantrack_market_descr + year + month ~ channel_type, value.var = "price_paid_norm_index") 
tmp_price	<- rep(1,numnodes) %*% t(colMeans(as.matrix(tmp[,4:(3+R)]), na.rm=T) )
omega_draw <- array(NA, c(numsim, K, numnodes))
for(i in 1:numsim){
	pct	<- proc.time()
	for(k in 1:length(bskt_type)){
		tmp_coef	<- as.vector(coef(sol_list[[k]]))
		sol 		<- incl_value_fn(param_est=tmp_coef, X_list=tmpX_list, y=tmpy, Q=bskt_typeQ[k], price=tmp_price, 
						R=R, Ra=R+1, qz_cons = -delta_q, exp_outside = FALSE, quant_outside = TRUE)
		omega_draw[i,k,] <- sol$omega
	}
	use.time <- proc.time() - pct
	cat("Simulation", i, ", using", use.time[3]/60,"min.\n")
}

tmp 			<- apply(omega_draw, c(2,3), median, na.rm=T)
dimnames(tmp)	<- list(bskt_type = 1:K, y = tmpy)
tmp				<- melt(tmp)
omega_draw_mean	<- data.frame(tmp[,1:2], omega = tmp$value)

# #--------------------------------------------------- # 
# # Compute the inclusve value within each scantrack_market_descr-year #
# numsim 	<- 100		# Number of simulation draws
# numnodes<- 20		# Number of interpolation nodes 
# K		<- length(bskt_type)
# mktyr_dat 	<- unique(fmt_attr[,c("scantrack_market_descr", "year", "name")] )
# mktyr_dat	<- mktyr_dat[order(mktyr_dat$name),]
# mktyr_dat$index <- 1:nrow(mktyr_dat)
# tmpy	<- quantile(hh_exp$net_dol, c(0:numnodes)/numnodes)
# 
# selcol	<- c("size_index", "ln_num_module", "ln_upc_per_mod", "overall_prvt")
# omega_draw_mkt <- array(NA, c(numsim, nrow(mktyr_dat),K, numnodes+1))
# for(i in 1:numsim){
# 	pct	<- proc.time()
# 	for(j in 1:nrow(mktyr_dat)){
# 		tmpX_list	<- lapply(fmt_name, function(x) 
# 						rep(1,numnodes+1) %*% as.matrix(subset(fmt_attr, name == mktyr_dat[j,"name"] & channel_type == x)[,selcol]) )		
# 		tmp 		<- dcast(subset(price_dat, scantrack_market_descr == mktyr_dat[j,"scantrack_market_descr"] & year == mktyr_dat[j,"year"]), 
# 							 scantrack_market_descr + year + month ~ channel_type, value.var = "price_paid_norm_index") 
# 		tmp_price	<- rep(1,numnodes+1) %*% t(colMeans(as.matrix(tmp[,4:(3+R)])) )
# 		for(k in 1:length(bskt_type)){
# 			tmp_coef	<- as.vector(coef(sol_list[[k]]))
# 			sol 		<- incl_value_fn(param_est=tmp_coef, X_list=tmpX_list, y=tmpy, Q=bskt_typeQ[k], price=tmp_price, 
# 							R=R, Ra=R+1, qz_cons = -delta_q, exp_outside = FALSE, quant_outside = TRUE)
# 			omega_draw_mkt[i,j,k,] <- sol$omega
# 		}
# 	}
# 	use.time <- proc.time() - pct
# 	cat("Simulation", i, ", using", use.time[3]/60,"min.\n")
# }
# 
# tmp 			<- apply(omega_draw_mkt, c(2,3,4), mean, na.rm=T)
# dimnames(tmp)	<- list(index = mktyr_dat$index, bskt_type = 1:K, y = tmpy)
# tmp				<- melt(tmp)
# names(tmp)[4] 	<- "omega"
# omega_dat		<- merge(mktyr_dat, tmp, by="index")
# omega_dat 		<- omega_dat[order(omega_dat$index, omega_dat$bskt_type, omega_dat$y),]


####################
# Save the results #
####################
ls()
rm(list=c("my_data", "my_data1", "Allocation_constr_fn","Allocation_fn","Allocation_nlop_fn","incl_value_fn","i","j","k","MDCEV_grad_fnC","MDCEV_ll_fnC",
		  "MDCEV_LogLike_fnC","MDCEV_wrapper","my_fn","src","tmp","tmp1","tmp2","tmp_sol","sel","sel1","selcol"))
save.image(paste("run_",run_id,"/result_5_",model_name,".rdata",sep=""))


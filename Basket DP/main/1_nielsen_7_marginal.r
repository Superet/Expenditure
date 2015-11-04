library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(Rcpp)
library(RcppArmadillo)
library(maxLik)
library(evd)
library(data.table)
library(scales)

options(error = quote({dump.frames(to.file = TRUE)}))
args <- commandArgs(trailingOnly = TRUE)
print(args)
if(length(args)>0){
    for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
    }
}
vid_save 	<- vid
run_id_save <- run_id
take_draw	<- TRUE

# setwd("~/Documents/Research/Store switching/processed data/Estimation/run_1")
# plot.wd 	<- "~/Desktop"
# sourceCpp("~/Documents/Research/Store switching/Exercise/Multiple_discrete_continuous_model/MDCEV_a1b1.cpp")
# source("~/Documents/Research/Store switching/Exercise/Multiple_discrete_continuous_model/0_Allocation_function.R")
# source("~/Documents/Research/Store switching/Exercise/Multiple_discrete_continuous_model/0_marginal_effect_function.R")

# setwd("/home/brgordon/ccv103/Exercise/run/run_5")
# setwd("/sscc/home/c/ccv103/Exercise/run/run_5")
setwd("/kellogg/users/marketing/2661703/Exercise/run/run_5")

#############
# Functions #
#############
convert_coef <- function(par, se, nx, R, beta0_base, pca){
	beta.old	<- par[1:nx]
	beta0.old	<- par[(nx+1):(nx+R-1)]
	M			<- pca$rotation
	beta.new	<- (M/pca$scale) %*% beta.old
	beta0.new	<- beta0.old - sum(pca$center * beta.new)
	return(c(beta.new, beta0.new, par[-(1:(nx+R-1))]))
}

##################
# Construct data #  
##################
load(paste("5_est_MDCEV_seg", seg_id, ".rdata",sep=""))

sourceCpp("../MDCEV_a1b1.cpp")
source("../0_Allocation_function.R")
source("../0_marginal_effect_function.R")

sel		<- paste("DOLP_", gsub("\\s", "_", fmt_name), sep="")
hh_exp$y<- hh_exp$dol_purchases
s		<- hh_exp[,sel]/hh_exp$y
summary(s)

# Explanatory variables
selcol		<- c("size_index", "ln_upc_per_mod", "ln_num_module","overall_prvt")
fmt_attr1	<- data.table(fmt_attr)
fmt_attr1	<- fmt_attr1[,list(size_index = median(size_index), num_module = median(num_module), 
								num_brands = median(num_brands), avg_brands_per_mod = median(avg_brands_per_mod), 
								num_upc	= median(num_upc), avg_upc_per_mod = median(avg_upc_per_mod), 
								avg_prvt_per_mod = median(avg_prvt_per_mod), overall_prvt = median(overall_prvt), 
								ln_num_module = log(median(num_module)), ln_upc_per_mod = log(median(avg_upc_per_mod))), 
						by = list(scantrack_market_descr, channel_type)]
fmt_attr1	<- data.frame(fmt_attr1)
fmt_attr1	<- fmt_attr1[order(fmt_attr1$scantrack_market_descr),]
summary(fmt_attr[,selcol])
summary(fmt_attr1[,selcol])
X0_list		<- lapply(fmt_name, function(x) as.matrix(subset(fmt_attr1, channel_type == x)[,selcol]) )

# Price 
tmp 		<- dcast(price_dat,	 scantrack_market_descr ~ channel_type, value.var = "bsk_price_paid_2004", mean) 				
tmp			<- tmp[order(tmp$scantrack_market_descr),]
mktyr		<- tmp[,1]
tmp_price	<- as.matrix(tmp[,-1])
summary(hh_exp$y)
numnodes	<- 30		# Number of interpolation nodes 
seq_y		<- quantile(hh_exp$dol_purchases, c(1:(numnodes-1)/numnodes), na.rm=T )

##############
# Simulation #
##############
set.seed(666)
mychange	<- .1
if(take_draw){
	numsim 		<- 20
	eps_draw	<- matrix(rgev(numsim*R), numsim, R)
}else{
	numsim		<- 1
	eps_draw	<- matrix(0, numsim, R)
}

#------------------------------------------------------#
# Simulate the effect of y
cat("Now simulate the effect of expenditure.\n")
marg_y		<- data.frame()
s0_base		<- array(0, c(D, length(seq_y), numsim, nrow(tmp_price), R))
s0_base_avg	<- array(0, c(D, length(seq_y), nrow(tmp_price), R))

for(d in 1:D){
	prc			<- proc.time()
	tmp_coef 	<- coef(sol_list[[d]])
	for(i in 1:length(seq_y)){
		tmp		<- rep(seq_y[i], nrow(tmp_price))
		s0_arr	<- array(0, c(numsim, length(tmp), R))
		for(k in 1:numsim){
			tmpsol	<- 	marginal_fn(tmp_coef, beta0_base, y = tmp, Inf, X0_list, X0_list, price0=tmp_price, price1=tmp_price, R=R, Ra=R, 
						qz_cons = 0, exp_outside = FALSE, quant_outside = FALSE, eps_draw = rep(1, length(tmp)) %*% t(eps_draw[k,]), single = TRUE)
			s0_arr[k,,]	<- tmpsol$s0
			s0_base[d,i,k,,] <- tmpsol$s0
		}
		
		# Trimming some draws
		if(take_draw){
			tmp_s0	<- apply(s0_arr, c(2,3), mean, na.rm=T, trim = .2)
		}else{
			tmp_s0	<- apply(s0_arr, c(2,3), mean, na.rm=T)
		}
		
		s0_base_avg[d,i,,] <- tmp_s0
		tmp		<- data.frame(seg_id = seg_id, d = d, Var = "y", y = seq_y[i], mktyr, s = tmp_s0)
		marg_y	<- rbind(marg_y, tmp)
	}
	use.time 	<- proc.time() - prc
	cat("Simulation for d =", d, "using", use.time[3]/60,"min.\n")
}

filename	<- ifelse(take_draw, paste("7_effect_seg",seg_id, ".rdata", sep=""), paste("7_effect_seg",seg_id, "_nodraw.rdata", sep=""))
save.image(file = filename)

#------------------------------------------------------#
cat("-------------------------------------------------\n")
cat("Now simulate the marginal effect of retail attributes.\n")
cond_marg	<- data.frame()
for(d in 1:D){
	tmp_coef 	<- coef(sol_list[[d]])
	for(i in 1:length(seq_y)){
		tmpy	<- rep(seq_y[i], nrow(tmp_price))
		prc			<- proc.time()
		for(j in 1:length(selcol)){
			for(k in 1:length(fmt_name)){
				# Changes in covariates
				tmp <- fmt_attr1
				sel	<- tmp$channel_type == fmt_name[k]
				if(selcol[j] == "ln_upc_per_mod"){
					tmp[sel, selcol[j]]	<- tmp[sel, selcol[j]] + log(1 + mychange)
				}else{
					tmp[sel, selcol[j]]	<- tmp[sel, selcol[j]] * (1 + mychange)
				}
				X1_list	<- lapply(fmt_name, function(x) as.matrix(subset(tmp, channel_type == x)[,selcol]) )
				
				# Simulate allocation for each epsilon draw
				s1_arr	<- array(0, c(numsim, nrow(tmp_price), R))				
				for(l in 1:numsim){
					tmpsol	<- 	marginal_fn(tmp_coef, beta0_base, y = tmpy, Inf, X0_list, X1_list, price0=tmp_price, price1=tmp_price, R=R, Ra=R, 
								qz_cons = 0, exp_outside = FALSE, quant_outside = FALSE, 
								eps_draw = rep(1, length(tmp)) %*% t(eps_draw[l,]), e0 = s0_base[d,i,l,,]*seq_y[i])
					s1_arr[l,,] <- tmpsol$s1
				}	
				if(take_draw){
					tmp_s1	<- apply(s1_arr, c(2,3), mean, na.rm=T, trim=.2)
				}else{
					tmp_s1	<- apply(s1_arr, c(2,3), mean, na.rm=T)
				}							
				
				tmp		<- data.frame(seg_id = seg_id, d = d, y=seq_y[i],Var = selcol[j], X = fmt_name[k], mktyr, 
										s0 = s0_base_avg[d,i,,], s1 = tmp_s1)
				cond_marg <- rbind(cond_marg, tmp)							
			}
		}
		use.time 	<- proc.time() - prc
		cat("Simulation for d =", d, "and y =", seq_y[i], "with", use.time[3]/60,"min.\n")
	}
}

# Price effect 
cat("-------------------------------------------------\n")
cat("Now simulate the price effect.\n")
for(d in 1:D){
	tmp_coef 	<- coef(sol_list[[d]])
	for(i in 1:length(seq_y)){
		tmpy	<- rep(seq_y[i], nrow(tmp_price))
		prc		<- proc.time()
		for(j in 1:length(fmt_name)){
			tmp		<- tmp_price
			tmp[,j]	<- tmp[,j] * (1 + mychange)
			
			# Simulate allocation for each epsilon draw
			s1_arr	<- array(0, c(numsim, nrow(tmp_price), R))				
			for(l in 1:numsim){
				tmpsol	<- 	marginal_fn(tmp_coef, beta0_base, y = tmpy, Inf, X0_list, X0_list, price0=tmp_price, price1=tmp, R=R, Ra=R, 
							qz_cons = 0, exp_outside = FALSE, quant_outside = FALSE, 
							eps_draw = rep(1, length(tmp)) %*% t(eps_draw[l,]), e0 = s0_base[d,i,l,,]*seq_y[i])
				s1_arr[l,,] <- tmpsol$s1
			}
			if(take_draw){
				tmp_s1	<- apply(s1_arr, c(2,3), mean, na.rm=T, trim=.2)
			}else{
				tmp_s1	<- apply(s1_arr, c(2,3), mean, na.rm=T)
			}								
			
			tmp		<- data.frame(seg_id = seg_id, d = d, y = seq_y[i], Var = "price", X = fmt_name[j], mktyr, s0 = s0_base_avg[d,i,,], s1 = tmp_s1)
			cond_marg <- rbind(cond_marg, tmp)							
		}
		use.time 	<- proc.time() - prc
		cat("Simulation for d =", d, "and y =", seq_y[i], "with", use.time[3]/60,"min.\n")
	}
}

filename	<- ifelse(take_draw, paste("7_effect_seg",seg_id, ".rdata", sep=""), paste("7_effect_seg",seg_id, "_nodraw.rdata", sep=""))
save.image(file = filename)


cat("File has been saved.\n")
cat("This program is done. \n")

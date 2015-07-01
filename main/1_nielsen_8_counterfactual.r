library(ggplot2)
library(reshape2)
library(scatterplot3d)
library(Rcpp)
library(RcppGSL)
library(glmmML)
library(plm)
library(zoo)
library(data.table)
library(VGAM)

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

# setwd("~/Documents/Research/Store switching/processed data/Estimation/run_2")
sourceCpp("~/Documents/Research/Store switching/Exercise/main/1_CDAMT_functions.cpp")
source('~/Documents/Research/Store switching/Exercise/Multiple_discrete_continuous_model/0_Allocation_function.R')
source('~/Documents/Research/Store switching/Exercise/Multiple_discrete_continuous_model/0_marginal_effect_function.R')

# setwd("/home/brgordon/ccv103/Exercise/run")
# setwd("/kellogg/users/marketing/2661703/Exercise/run")
setwd("/sscc/home/c/ccv103/Exercise/run")

load("7_effect_aggr.rdata")
load("6_est_aggr.rdata")
sourceCpp("1_CDAMT_functions.cpp")
source("0_Allocation_function.R")
source("0_marginal_effect_function.R")

make_plot	<- FALSE

#############
# Functions #
#############
param_assign <- function(param, nx, R, base = 0){
	beta	<- param[1:nx]
	beta0	<- param[(nx+1):(nx+R-1)]
	gamma	<- param[(nx+R):(nx+2*R-1)]
	sigma	<- exp(param[(nx+2*R)])
	if(base==0){
		beta0 <- c(0, beta0)
	}else{
		beta0	<- c(beta0[1:(base-1)], 0, beta0[base:length(beta0)])
	}
	return(list(beta = beta, beta0=beta0, gamma = gamma, sigma = sigma))
}

sim_oneseg	<- function(hh_index, TT_vec, init_state, DataState, DP_list, control_list, logit_draw, margin_y_arr, y_grid){
	choice_arr		<- matrix(NA, nrow(DataState), 3)
	share_arr		<- matrix(NA, nrow(DataState), R)
	seqout 		<- sim_hhseqC(hh_index, TT_vec, init_state, DataState, choice_arr, DP_list, control_list, rep(NA, nrow(DataState)), 
						simidx = 3, logit_draw_mat = logit_draw)
	choice_out	<- seqout$choice_seq
	choice_arr <- choice_out
	
	for(d in 1:D){
		sel 	<- choice_out[,1] == d
		if(length(sel) > 0){
			idx				<- ceiling(choice_out[sel,2]/10)
			idx				<- cut(choice_out[sel,2], c(0,y_grid), include.lowest = TRUE, labels=1:length(y_grid))
			idx[is.na(idx)]	<- length(y_grid)
			share_arr[sel,]	<- margin_y_arr[d,idx,]
		}
	}
	share_arr[is.na(share_arr)] <- 0
	
	return(list(choice = choice_out, share = share_arr))
}

sim_allseg	<- function(TT_vec_all, init_state_all, DataState_all, DP_list_all, control_list, t_all, logit_draw_all, marg_y){
	state_out	<- NULL	
	choice_out	<- NULL
	share_out	<- NULL
	y_grid		<- sort(unique(marg_y$y))
	
	for(i in 1:length(TT_vec_all)){
		# Organze the array of marginal effect of y 
		tmpy		<- unique(marg_y$y)
		tmp			<- subset(marg_y, seg_id == sel_seg[[i]])
		ord 		<- with(tmp, order(d, y))
		selcol		<- which(names(marg_y) %in% paste("s.",1:R,sep=""))
		marg_y_arr	<- array(0, c(D, length(tmpy), R))
		for(d in 1:D){
			sel				<- seq_len(length(tmpy)) + length(tmpy)*(d-1)
			marg_y_arr[d,,]	<- as.matrix(tmp[sel,selcol])
		}
		
		tmpsol		<- sim_oneseg(1:length(TT_vec_all[[i]]), TT_vec_all[[i]], init_state_all[[i]], DataState_all[[i]], 
								DP_list_all[[i]], control_list, logit_draw_all[[i]], marg_y_arr, y_grid)
		choice_out	<- rbind(choice_out, tmpsol$choice)
		share_out	<- rbind(share_out, tmpsol$share)						
		state_out	<- rbind(state_out, DataState_all[[i]])
	}
	
	# Aggregation of simulation
	seg			<- do.call("c", sapply(1:length(DataState_all), function(i) rep(i, nrow(DataState_all[[i]]))))
	tmpdata		<- data.frame(seg = seg, t = t_all, choice_out, share_out)
	tmpdata$Rc	<- 1*(tmpdata$t >= Rc_t)
	names(tmpdata) <- c("seg","t","d","y","c",paste("s.",1:R,sep=""),"Rc")
	tmpdt		<- data.table(tmpdata)
	
	# 1. Aggregate over all the observations: before and after 
	bf_af		<- tmpdt[,list(d = mean(d), y=mean(y), c = mean(c),
							   s.1 = mean(s.1), s.2 = mean(s.2), s.3 = mean(s.3), 
							   s.4 = mean(s.4), s.5 = mean(s.5), s.6 = mean(s.6)), by=list(Rc)]
	
	# 2. Aggregate over segments: before and after
	bf_af_seg	<- tmpdt[,list(d = mean(d), y=mean(y), c = mean(c),
							   s.1 = mean(s.1), s.2 = mean(s.2), s.3 = mean(s.3), 
							   s.4 = mean(s.4), s.5 = mean(s.5), s.6 = mean(s.6)), by=list(seg, Rc)]
	
	# 3. Aggregate over each period
	sim_t		<- tmpdt[,list(d = mean(d), y=mean(y), c = mean(c),
							   s.1 = mean(s.1), s.2 = mean(s.2), s.3 = mean(s.3), 
							   s.4 = mean(s.4), s.5 = mean(s.5), s.6 = mean(s.6)), by=list(t)]
							
	# 4. Aggregate over each perod by segment
	sim_t_seg	<- tmpdt[,list(d = mean(d), y=mean(y), c = mean(c),
							   s.1 = mean(s.1), s.2 = mean(s.2), s.3 = mean(s.3), 
							   s.4 = mean(s.4), s.5 = mean(s.5), s.6 = mean(s.6)), by=list(seg, t)]
	
	return(list(bf_af = bf_af, bf_af_seg = bf_af_seg, sim_t = sim_t, sim_t_seg = sim_t_seg))
}

# Plot policy function
for(i in 1:length(DP_list_all)){
	DP_list_all[[i]]	<- policy_iterC(DP_list_all[[i]], control_list, print = 1)
}

# -------------------------------------------# 
# Fill in the no-trip weeks
mydata		<- subset(hh_exp, !is.na(segment))
tmp			<- data.table(mydata)
tmp			<- tmp[,list(minw=min(biweek), maxw=max(biweek)), by=list(segment, household_code)]
tmp			<- tmp[,list(biweek=minw:maxw), by=list(segment, household_code)]
cat("Before filling in non-shoppiing periods, dim(mydata) =", dim(mydata),"\n") 
tmpb		<- with(mydata, paste(household_code, biweek, sep="*"))
mydata		<- merge(data.frame(tmp), mydata, by=c("segment","household_code","biweek"), all.x=T)
cat("After filling in non-shopping periods, dim(mydata)=", dim(mydata),"\n")
tmpa		<- with(mydata, paste(household_code, biweek, sep="*"))
sel			<- !tmpa %in% tmpb
cat("Sum of non-shopping biweeks =", sum(sel), "\n")

# Expenditure duing noshoping trips is 0
mydata[sel,"dol_purchases"] <- 0
mydata[sel,"d"] 			<- 0
sum(is.na(mydata$d))
mydata[is.na(mydata$dol_purchases),"dol_purchases"] <- 0
mydata[is.na(mydata$d),"d"] 						<- 0

# Fill in the missing income and recession 
tmp			<- min(hh_exp$biweek):max(hh_exp$biweek)
bdate		<- as.Date("2004-01-01", format="%Y-%m-%d") + 14*(tmp-1)		
byear		<- year(bdate)
tmp			<- byear[mydata[sel,"biweek"]]
mydata[sel,"recession"]	<- ifelse(tmp>=2008, 1, 0)
mydata[sel,"year"]		<- tmp
ord			<- order(mydata$household_code, mydata$year, mydata$biweek)
mydata		<- mydata[ord,]

sel1 		<- is.na(mydata$income_nodes)
tmp1 		<- data.table(mydata[!sel1,])
tmp1		<- tmp1[,list(income_nodes=unique(income_nodes), income_real=unique(income_real)), by=list(household_code, year)]
tmp2		<- merge(mydata[sel1,c("household_code","year")], tmp1, by=c("household_code","year"),all.x=T)
ord			<- order(tmp2$household_code, tmp2$year)
tmp2		<- tmp2[ord,]
mydata[sel1,"income_real"]	<- tmp2$income_real
mydata[sel1,"income_nodes"]	<- tmp2$income_nodes
mydata$income_nodes 		<- factor(mydata$income_nodes)

# Sort data 
ord			<- order(mydata$segment, mydata$household_code, mydata$year, mydata$biweek)
mydata		<- mydata[ord,]
dim(mydata)
t_all		<- mydata$biweek

choice_seq	<- data.frame(d=mydata$d, y=mydata$dol_purchases, c=NA)
DataState 	<- data.frame(I = rep(NA, nrow(mydata)), Inc = as.numeric(mydata$income_nodes), 
							   Rc=as.numeric(mydata$recession), lnZ = NA)
cat("Summary of choice variables:\n"); print(summary(choice_seq)); cat("\n")
cat("Summary of state variables:\n"); print(summary(DataState)); cat("\n")

# -------------------------------------------# 
# Household index and initial stata #
myidx		<- data.table(mydata)
myidx		<- myidx[,list(minw=min(biweek), maxw=max(biweek), stay=max(biweek)-min(biweek)+1, avg_Q = mean(Q,na.rm=T), 
					first_Inc = as.numeric(income_nodes[1]), first_Rc = recession[1]), 
					by=list(segment, household_code)]
myidx		<- myidx[,id:=1:nrow(myidx)]
cat("Check cross-household distribution of stay patten:\n"); summary(myidx); cat("\n")

hh_index 	<- as.vector(myidx$id)
TT_vec		<- as.vector(myidx$stay)
init_state	<- data.frame(myidx[,list(avg_Q, first_Inc, first_Rc)],lnZ=rep(NA,nrow(myidx)))
cat("Initial state:\n"); summary(init_state); cat("\n")

# Check data consistency
cat("dim(mydata) =", dim(mydata),";")
cat("sum(TT_vec) =", sum(TT_vec),";")
cat("dim(DataState) =", dim(DataState),";")
cat("dim(choice_seq) =",dim(choice_seq),"\n")

# -------------------------------------------# 
# Split data states by segment 
TT_vec_all		<- split(TT_vec, myidx$segment)
init_state_all	<- split(init_state, myidx$segment)
init_state_all	<- lapply(init_state_all, as.matrix)
DataState_all	<- split(DataState, mydata$segment)
DataState_all	<- lapply(DataState_all, as.matrix)

set.seed(666)
sel_seg		<- 1:12
logit_draw_all <- lapply(1:length(DataState_all), function(i) matrix(rgev(nrow(DataState_all[[i]])*(D+1)), nrow(DataState_all[[i]]), D+1) )
Rc_t		<- min(hh_exp[hh_exp$recession == 1, "biweek"])

###################################################################################
# Counterfactual 1: decomposing the effect of income level and income uncertainty #
###################################################################################
# Scenario: baseline, no uncertainty, no incone change
# Baseline simulation 
sim_base	<-sim_allseg(TT_vec_all, init_state_all, DataState_all, DP_list_all, control_list, t_all, logit_draw_all, marg_y)

# No income uncertainty
DataState_new1	<- DataState_all
for(i in 1:length(DataState_new1)){
	DataState_new1[[i]][,3] <- 0
}

sim1	<-sim_allseg(TT_vec_all, init_state_all, DataState_new1, DP_list_all, control_list, t_all, logit_draw_all, marg_y)

# No income changes
DataState_new2	<- DataState_all
for(i in 1:length(DataState_new2)){
	cnt	<- 0
	for(h in 1:length(TT_vec_all[[i]])){
		sel		<- cnt + seq_len(TT_vec_all[[i]][h])
		sel1	<- sel[DataState_new2[[i]][sel,3] == 1]
		sel0	<- setdiff(sel,sel1)
		if(length(sel1)>0 & length(sel0)>0){
			tmp	<- DataState_new2[[i]][sel0[length(sel0)],2]
			DataState_new2[[i]][sel1,2]	<- tmp
		}
		cnt		<- cnt + TT_vec_all[[i]][h]
	}
}

sim2	<-sim_allseg(TT_vec_all, init_state_all, DataState_new2, DP_list_all, control_list, t_all, logit_draw_all, marg_y)

# Combine the data together
bf_af		<- NULL
bf_af_seg	<- NULL
sim_t		<- NULL
bf_af		<- rbind(bf_af, data.frame(Scenario = "Base", sim_base$bf_af))
bf_af_seg	<- rbind(bf_af_seg, data.frame(Scenario = "Base", sim_base$bf_af_seg))
sim_t		<- rbind(sim_t, data.frame(Scenario = "Base", sim_base$sim_t))

bf_af		<- rbind(bf_af, data.frame(Scenario = "NoUncertainty", sim1$bf_af))
bf_af_seg	<- rbind(bf_af_seg, data.frame(Scenario = "NoUncertainty", sim1$bf_af_seg))
sim_t		<- rbind(sim_t, data.frame(Scenario = "NoUncertainty", sim1$sim_t))

bf_af		<- rbind(bf_af, data.frame(Scenario = "IncomeSame", sim2$bf_af))
bf_af_seg	<- rbind(bf_af_seg, data.frame(Scenario = "IncomeSame", sim2$bf_af_seg))
sim_t		<- rbind(sim_t, data.frame(Scenario = "IncomeSame", sim2$sim_t))

bf_af		<- subset(bf_af, !(Scenario!="Base" & Rc ==0))
tmp.tab		<- t(as.matrix(bf_af[,-1]))
colnames(tmp.tab) <- c("Pre-recession","Base post", "NoUncertainty post", "IncomeSame post")
cat("The effect of income reduction and income uncertainty:\n"); print(tmp.tab); cat("\n")

#####################################
# Counterfactual 2: Grocery changes # 
#####################################
# Grocery increase 10% size index during the recession
#----------------------------------------------------------------- # 
# Simulate omega for each segment
tmpy		<- unique(marg_y$y)
numnodes	<- length(tmpy)
fmt_name 	<- as.character(sort(unique(fmt_attr$channel_type)))

beta0_base	<- 5
sel_fmt		<- "Grocery"
sel_col		<- "overall_prvt"
mychange 	<- .1

tmp			<- fmt_attr
sel			<- tmp$channel_type == sel_fmt
tmp[sel,sel_col] <- tmp[sel,sel_col] * (1 + mychange)
selcol		<- c("size_index", "ln_upc_per_mod", "overall_prvt")
nx			<- length(selcol)
tmpX_list	<- lapply(fmt_name, function(x) 
				rep(1,numnodes) %*% t(colMeans(as.matrix(subset(tmp, channel_type == x)[,selcol]))) )		
tmp 		<- dcast(price_dat,	 scantrack_market_descr + year + biweek ~ channel_type, value.var = "price_paid_norm_index") 
tmp_price	<- rep(1,numnodes) %*% t(colMeans(as.matrix(tmp[,4:(3+R)]), na.rm=T) )

# For each segment: solve for DP using the new values of omega
DP_list_new	<- DP_list_all
for(i in seg_index){
	omega_nodraw <- matrix(NA, D, numnodes, dimnames = list(d = 1:D, y = tmpy))
	for(d in 1:D){
		tmp_coef	<- alpha_list[[i]][,d]
		sol 		<- incl_value_fn(param_est=tmp_coef, base= beta0_base, X_list=tmpX_list, y=tmpy, Q=Inf, price=tmp_price, 
						R=R, Ra=R, qz_cons = 0, exp_outside = FALSE, quant_outside = FALSE, eps_draw = matrix(0, numnodes,R))
		omega_nodraw[d,] <- sol$omega
	}
	DP_list_new[[i]]$omega 	<- list(omega_y = tmpy, omega = omega_nodraw)
	DP_list_new[[i]]		<- policy_iterC(DP_list_new[[i]], control_list, print_level = 1)
	cat("Having update DP_list for segment",i,"\n")
}

# Simulate the share accounting for shopping strategies
sim_size	<-sim_allseg(TT_vec_all, init_state_all, DataState, DP_list_new, control_list, t_all, logit_draw_all, marg_y)

bf_af		<- rbind(bf_af, data.frame(Scenario = "Size", sim_size$bf_af))
bf_af_seg	<- rbind(bf_af_seg, data.frame(Scenario = "Size", sim_size$bf_af_seg))
sim_t		<- rbind(sim_t, data.frame(Scenario = "Size", sim_size$sim_t))





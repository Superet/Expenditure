library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(Rcpp)
library(RcppArmadillo)
library(maxLik)
library(evd)
library(data.table)
library(scales)
library(gridExtra)

options(error = quote({dump.frames(to.file = TRUE)}))
args <- commandArgs(trailingOnly = TRUE)
print(args)
if(length(args)>0){
    for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
    }
}

# setwd("~/Documents/Research/Store switching/processed data/Estimation/run_2")
# plot.wd 	<- "~/Desktop"
# sourceCpp("~/Documents/Research/Store switching/Exercise/Multiple_discrete_continuous_model/MDCEV_a1b1.cpp")
# source("~/Documents/Research/Store switching/Exercise/Multiple_discrete_continuous_model/0_Allocation_function.R")
# source("~/Documents/Research/Store switching/Exercise/Multiple_discrete_continuous_model/0_marginal_effect_function.R")

# setwd("/home/brgordon/ccv103/Exercise/run/run_5")
# setwd("/sscc/home/c/ccv103/Exercise/run/run_5")
setwd("/kellogg/users/marketing/2661703/Exercise/run/run_5")
plot.wd <- getwd()

tmpcsv		<- paste(plot.wd,"/7_marginal_effect.csv", sep="")
write2csv	<- TRUE
make_plot	<- FALSE
ww			<- 8
ar 			<- .6

if(write2csv){
	f 		<- file(tmpcsv, "w")
	writeLines("Marginal effects of MDCEV.\n", f)
	close(f)
}

# ############################################
# # Contrasting marginal effects by segments # 
# ############################################
# # Read in data 
# sel_seg 	<- c(9, 12)
# cond_marg2	<- data.frame()
# marg_price2	<- data.frame()
# marg_y2		<- data.frame()
# 
# for(ii in 1:length(sel_seg)){
# 	load(paste("7_effect_seg",sel_seg[ii], ".rdata", sep=""))
# 	cond_marg2	<- rbind(cond_marg2, cond_marg)
# 	marg_price2	<- rbind(marg_price2, marg_price)
# 	marg_y2		<- rbind(marg_y2, marg_y)
# 	rm(list = c("mydata", "cond_marg", "marg_price","marg_y"))
# }
# 
# # Compute market weight
# mkt_wt		<- data.table(hh_exp)
# mkt_wt		<- mkt_wt[,list(nhh = length(unique(household_code))), by = list(scantrack_market_descr)]
# mkt_wt		<- mkt_wt[,weight:=nhh/sum(nhh)]
# ny			<- length(unique(cond_marg2$year))
# mkt_wt		<- data.frame(mkt_wt)
# mkt_wt$mkt_wt	<- mkt_wt$weight / ny
# 
# # Conditional marginal effect
# tmp 	<- unique(cond_marg2$Var)
# for(i in 1:length(sel_seg)){
# 	for(j in 1:length(tmp)){
# 		ggtmp		<- subset(cond_marg2, Var == tmp[j] & seg_id == sel_seg[i])
# 		selcol0		<- paste("s0.", 1:R, sep="")
# 		selcol1		<- paste("s1.", 1:R, sep="")
# 		ggtmp		<- cbind(ggtmp[,1:6], ggtmp[,selcol1] - ggtmp[,selcol0])
# 		ggtmp		<- melt(ggtmp, id.vars = c("seg_id","d","Var","X","scantrack_market_descr","year"))
# 		ggtmp$variable <- factor(ggtmp$variable, levels=paste("s1.", 1:R,sep=""), label = gsub("\\s", "\n",fmt_name))
# 		ggtmp$value <- ggtmp$value * 100
# 		
# 		if(make_plot){
# 			pdf(paste(plot.wd,"/graph_condm_seg",sel_seg[i],"_",tmp[j], ".pdf", sep=""), width = ww, height = ww*ar)
# 			print(ggplot(ggtmp, aes(X, value, color = factor(d))) + 
# 					stat_summary(fun.y=mean, fun.ymin=function(x) quantile(x, .25), fun.ymax = function(x) quantile(x, .75), 
# 								geom ="pointrange", size=.5, position=position_dodge(width=0.5)) + 
# 					coord_flip() + 
# 					facet_grid(.~variable) + 
# 					scale_colour_grey() + theme_bw() + 
# 					# scale_y_continuous(labels=percent) + 
# 					labs(x = "Retail format", y = "Expenditure share difference(%)")
# 			)
# 			dev.off()
# 		}
# 	}
# }
# 
# # Price effect 
# for(i in 1:length(sel_seg)){
# 	ggtmp		<- subset(marg_price2, seg_id == sel_seg[i])
# 	selcol0		<- paste("s0.", 1:R, sep="")
# 	selcol1		<- paste("s1.", 1:R, sep="")
# 	ggtmp		<- cbind(ggtmp[,1:6], ggtmp[,selcol1] - ggtmp[,selcol0])
# 	ggtmp		<- melt(ggtmp, id.vars = c("seg_id","d","Var","X","scantrack_market_descr","year"))
# 	ggtmp$variable <- factor(ggtmp$variable, levels=paste("s1.", 1:R,sep=""), label = gsub("\\s", "\n",fmt_name))
# 	ggtmp$value <- ggtmp$value * 100
# 	
# 	if(make_plot){
# 		pdf(paste(plot.wd,"/graph_condm_seg",sel_seg[i],"_price.pdf", sep=""), width = ww, height = ww*ar)
# 		print(ggplot(ggtmp, aes(X, value, color = factor(d))) + 
# 				stat_summary(fun.y=mean, fun.ymin=function(x) quantile(x, .25), fun.ymax = function(x) quantile(x, .75), 
# 							geom ="pointrange", size=.5, position=position_dodge(width=0.5)) + 
# 				coord_flip() + 
# 				facet_grid(.~variable) + 
# 				scale_colour_grey() + theme_bw() + 
# 				# scale_y_continuous(labels=percent) + 
# 				labs(x = "Retail format", y = "Expenditure share difference(%)")
# 		)
# 		dev.off()
# 	}
# }
# 
# # Expenditure effect 
# ggtmp 		<- melt(marg_y2, id.vars=c("seg_id","d","Var","X","scantrack_market_descr","year"))
# ggtmp		<- merge(ggtmp, mkt_wt[,c("scantrack_market_descr","mkt_wt")], all.x=T)
# ggtmp		<- data.table(ggtmp)
# ggtmp		<- ggtmp[,list(value = sum(value*mkt_wt)), by=list(seg_id, d, X, variable)]
# ord			<- with(ggtmp, order(seg_id, d, X, variable))
# ggtmp		<- ggtmp[ord, ]
# ggtmp$variable <- factor(ggtmp$variable, levels=paste("s.",1:R,sep=""), labels=fmt_name)
# ggtmp$seg	<- names(seg_index)[ggtmp$seg_id]
# 
# if(make_plot){
# 	pdf(paste(plot.wd,"/graph_effecty_seg", paste(sel_seg,collapse="_"),".pdf", sep=""), width = ww, height = ww*ar)
# 	print(ggplot(ggtmp, aes(X, value, fill = variable)) + geom_bar(stat="identity") + 
# 			facet_grid(seg ~ d) + 
# 			labs(x = "Expenditure y", y = "Expenditure share")
# 		)
# 	dev.off()	
# }

############################################################
# Aggregating the marginal effects across all demographics # 
############################################################
sel_seg		<- 7:12
take_draw	<- TRUE

# Read in data 
cond_marg_all	<- data.frame()
marg_y_all		<- data.frame()
alpha_list		<- vector("list", length(sel_seg))

for(ii in 1:length(sel_seg)){
	f			<- ifelse(take_draw, paste("7_effect_seg",sel_seg[ii], ".rdata", sep=""), 
									paste("7_effect_seg",sel_seg[ii], "_nodraw.rdata", sep=""))
	load(f)
	cond_marg_all	<- rbind(cond_marg_all, subset(cond_marg, y>0))
	marg_y_all		<- rbind(marg_y_all, subset(marg_y, y>0))
	alpha_list[[ii]]<- sapply(sol_list, coef)
	rm(list = c("marg_y", "cond_marg", "sol_list"))
	cat("Having load segment",ii,".\n")
}

#----------------------------------------------------------# 
# Compute weight 
# Compute market weight
mkt_wt		<- data.table(hh_exp)
mkt_wt		<- mkt_wt[,list(nhh = length(unique(household_code))), by = list(scantrack_market_descr)]
mkt_wt		<- mkt_wt[,mkt_wt:=nhh/sum(nhh)]
mkt_wt		<- data.frame(mkt_wt)

# F(y,d | s)
seq_y		<- unique(marg_y_all$y)
tmp			<- hh_exp[,c("segment", "dol_purchases", "d")]
tmp$ybin	<- cut(tmp$dol_purchases, c(0,seq_y), include.lowest = T, label=seq_y)
tmp[is.na(tmp$ybin), "ybin"] <- max(seq_y)
yd_wt		<- array(0, c(length(seg_index), length(seq_y), D))
for(i in seg_index){
	sel		<-  tmp$segment == i
	tmp.tab	<- table(tmp[sel,"ybin"], tmp[sel,"d"])
	yd_wt[i,,] <- tmp.tab/sum(tmp.tab)
}

# Segment weight
seg_wt		<- table(panelist$first_income,panelist$first_famsize)
seg_wt		<- seg_wt/sum(seg_wt)

#----------------------------------------------------------# 
# Elasticty of retail attributes
tmpdata		<- merge(cond_marg_all, mkt_wt[,c("scantrack_market_descr","mkt_wt")], by.x="mktyr", by.y="scantrack_market_descr", all.x=T)
tmp			<- melt(yd_wt)
names(tmp)	<- c("seg_id","y","d","yd_wt")
tmp[,"y"]	<- seq_y[tmp$y]
tmpdata		<- merge(tmpdata, tmp, by=c("seg_id","y","d"), all.x=T)
tmp			<- c(seg_wt)
tmpdata$seg_wt	<- tmp[tmpdata$seg_id]

# Compute elasticity conditional on positive s0
mys_min		<- 1e-4
sel0		<- paste("s0.", 1:R, sep="")
sel1		<- paste("s1.", 1:R, sep="")
tmp0		<- tmpdata[,sel0]
summary(tmp0)
tmp0[tmp0<mys_min]	<- NA

tmp.e		<- as.matrix(((tmpdata[,sel1] - tmp0)/tmp0)/(mychange/.01) * 100)
summary(tmp.e)
myalpha		<- .2
tmp			<- apply(tmp.e, 2, quantile, c(.5*myalpha, 1 - .5*myalpha), na.rm=T)
for(i in 1:ncol(tmp.e)){
	sel		<- tmp.e[,i] < tmp[1,i] | tmp.e[,i] > tmp[2,i]
	tmp.e[sel,i] <- NA
}
summary(tmp.e)
tmpdata		<- cbind(tmpdata, e = tmp.e)

# Elasticty by segment, d, y
cond_marg	<- data.table(tmpdata)
cond_marg	<- cond_marg[,list(e1 = sum(e.s1.1*mkt_wt, na.rm=T), e2 = sum(e.s1.2*mkt_wt, na.rm=T),
							   e3 = sum(e.s1.3*mkt_wt, na.rm=T), e4 = sum(e.s1.4*mkt_wt, na.rm=T),
							   e5 = sum(e.s1.5*mkt_wt, na.rm=T), e6 = sum(e.s1.6*mkt_wt, na.rm=T)), 
						by = list(seg_id, seg_wt, y, d, yd_wt, Var, X)]
setkeyv(cond_marg, c("seg_id","d","y","Var","X"))	
cond_marg	<- data.frame(cond_marg)					

# Overall elasticity
tmp.dt		<- data.table(tmpdata)
tmp.dt		<- tmp.dt[,list(e1 = sum(e.s1.1*mkt_wt*yd_wt*seg_wt, na.rm=T)/sum(1*(!is.na(e.s1.1))*mkt_wt*yd_wt*seg_wt), 
							e2 = sum(e.s1.2*mkt_wt*yd_wt*seg_wt, na.rm=T)/sum(1*(!is.na(e.s1.2))*mkt_wt*yd_wt*seg_wt),
							e3 = sum(e.s1.3*mkt_wt*yd_wt*seg_wt, na.rm=T)/sum(1*(!is.na(e.s1.3))*mkt_wt*yd_wt*seg_wt), 
							e4 = sum(e.s1.4*mkt_wt*yd_wt*seg_wt, na.rm=T)/sum(1*(!is.na(e.s1.4))*mkt_wt*yd_wt*seg_wt),
							e5 = sum(e.s1.5*mkt_wt*yd_wt*seg_wt, na.rm=T)/sum(1*(!is.na(e.s1.5))*mkt_wt*yd_wt*seg_wt), 
							e6 = sum(e.s1.6*mkt_wt*yd_wt*seg_wt, na.rm=T)/sum(1*(!is.na(e.s1.6))*mkt_wt*yd_wt*seg_wt)), 
						by = list(Var, X)]
tmp.dt		<- tmp.dt[order(tmp.dt$Var, tmp.dt$X)]
cat("Aggregated elasticity matrix:\n"); print(tmp.dt); cat("\n")

# Elasticity by income group
tmp			<- c(seg_wt/rowSums(seg_wt))
tmpdata$Qt_wt	<- tmp[tmpdata$seg_id]
tmp			<- substr(names(seg_index), 1, 3)
tmpdata$Qt	<- tmp[tmpdata$seg_id]
tmp.dt1		<- data.table(tmpdata)
tmp.dt1		<- tmp.dt1[,list(e1 = sum(e.s1.1*mkt_wt*yd_wt*Qt_wt, na.rm=T)/sum(1*(!is.na(e.s1.1))*mkt_wt*yd_wt*Qt_wt), 
							 e2 = sum(e.s1.2*mkt_wt*yd_wt*Qt_wt, na.rm=T)/sum(1*(!is.na(e.s1.2))*mkt_wt*yd_wt*Qt_wt),
							 e3 = sum(e.s1.3*mkt_wt*yd_wt*Qt_wt, na.rm=T)/sum(1*(!is.na(e.s1.3))*mkt_wt*yd_wt*Qt_wt), 
							 e4 = sum(e.s1.4*mkt_wt*yd_wt*Qt_wt, na.rm=T)/sum(1*(!is.na(e.s1.4))*mkt_wt*yd_wt*Qt_wt),
							 e5 = sum(e.s1.5*mkt_wt*yd_wt*Qt_wt, na.rm=T)/sum(1*(!is.na(e.s1.5))*mkt_wt*yd_wt*Qt_wt), 
							 e6 = sum(e.s1.6*mkt_wt*yd_wt*Qt_wt, na.rm=T)/sum(1*(!is.na(e.s1.6))*mkt_wt*yd_wt*Qt_wt)), 
						by = list(Var, Qt, X)]
setkeyv(tmp.dt1, c("Var","Qt","X"))
cat("Aggregated elasticity matrix by segment:\n"); print(tmp.dt1); cat("\n")

if(write2csv){
	f		<- file(tmpcsv, "at")
	writeLines("Expected elasticity:\n", f)
	write.csv(tmp.dt, f)
	writeLines("Elasticty by segment:\n", f)
	write.csv(tmp.dt1, f)
	close(f)
}

# Plot the elasticity difference by segment

#----------------------------------------------------------# 
# Conditional allocation rule 
marg_y	<- merge(marg_y_all, mkt_wt[,c("scantrack_market_descr","mkt_wt")], by.x="mktyr", by.y="scantrack_market_descr", all.x=T)
tmp		<- c(seg_wt)
marg_y$seg_wt	<- tmp[marg_y$seg_id]
marg_y	<- data.table(marg_y)
marg_y	<- marg_y[,list(	s.1 = sum(s.1*mkt_wt), s.2 = sum(s.2*mkt_wt), 
							s.3 = sum(s.3*mkt_wt), s.4 = sum(s.4*mkt_wt),
							s.5 = sum(s.5*mkt_wt), s.6 = sum(s.6*mkt_wt)), 
					by=list(seg_id, d, y, seg_wt)]
setkeyv(marg_y, c("seg_id","d","y"))
					
marg_y1	<- marg_y[,list(s.1 = sum(s.1*seg_wt), s.2 = sum(s.2*seg_wt), s.3 = sum(s.3*seg_wt), 
						s.4 = sum(s.4*seg_wt), s.5 = sum(s.5*seg_wt), s.6 = sum(s.6*seg_wt)), 
					by = list(d,y)]
marg_y	<- data.frame(marg_y)	
marg_y1	<- data.frame(marg_y1)				

plots	<- list(NULL)
ggtmp	<- melt(marg_y1, id.vars=c("d","y"))
ggtmp$variable <- factor(ggtmp$variable, levels = paste("s.", 1:R, sep=""), labels=fmt_name)
plots[[1]] <- ggplot(ggtmp, aes(y, value, fill = variable)) + geom_bar(stat ="identity") + 
					facet_wrap(~d)

# Verify the share allocation in the data
tmp		<- hh_exp[,c("segment", "dol_purchases", "d")]
selcol	<- paste("DOLP_", gsub("\\s","_", fmt_name), sep="")
tmp		<- cbind(tmp, hh_exp[,selcol]/hh_exp[,"dol_purchases"])
tmp$ybin<- cut(tmp$dol_purchases, c(0,seq_y), include.lowest = T, label=seq_y)
tmp[is.na(tmp$ybin), "ybin"] <- max(seq_y)
tmp		<- data.table(melt(tmp[,c("d","ybin",selcol)], id.vars=c("ybin","d")))
tmp		<- tmp[,list(value = mean(value)), by = list(ybin, d, variable)]
tmp$variable <- factor(tmp$variable, levels = selcol, labels=fmt_name)

plots[[2]] <- ggplot(tmp, aes(ybin, value, fill = variable)) + geom_bar(stat ="identity") + 
					facet_wrap(~d)

save(cond_marg_all, cond_marg, marg_y, price_dat, fmt_attr, mkt_wt, yd_wt, seg_wt, alpha_list, panelist, R, fmt_name, seg_index, 
	file = "7_effect_aggr.rdata")

cat("This program is done.\n")


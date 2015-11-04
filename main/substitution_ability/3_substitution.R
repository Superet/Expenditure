library(r2excel)
library(data.table)
library(ggplot2)
library(Hmisc)
library(reshape2)
library(gridExtra)
library(plm)
options(error = quote({dump.frames(to.file = TRUE)}))

# setwd("~/Documents/Research/Store switching/processed data")
# source('~/Documents/Research/Store switching/Exercise/main/outreg function.R')
setwd("/kellogg/users/marketing/2661703/Exercise")

# Plot setting
# plot.wd	<- "~/Desktop"
plot.wd	<- getwd()
ww		<- 6.5
ar 		<- .6
incg.col<- c("red","grey30", "grey50")			# Color of low, median, high income bar
make_plot	<- TRUE

# Function to change factor levels and label
changelab	<- function(x, ord, lab){
	y	<- rep(NA, length(x))
	for(i in 1:length(ord)){
		sel		<- grep(ord[i], x)
		y[sel]	<- lab[i]
	}
	y	<- factor(y, levels = unique(lab))
	return(y)
}

load("substitution_hh_exp.rdata")
hh_exp$first_income	<- factor(hh_exp$first_income, levels = c("T1", "T2", "T3"), labels = c("Low", "Med", "High"))
hh_exp$dol_rt		<- hh_exp$dol/hh_exp$income_cpi
hh_exp$ln_dol		<- log(hh_exp$dol)
hh_exp$ln_income	<- log(hh_exp$income_cpi)

###########################################################
# Distribution cross income groups prior to the recession # 
###########################################################
tmp.data	<- data.table(subset(hh_exp, Rc == 0))
tmp.data	<- tmp.data[,list(dol = mean(dol), 
						PCTL_UPRICE = mean(PCTL_UPRICE, na.rm =T), PCTL_SIZE = mean(PCTL_SIZE, na.rm =T), 
						SHR_Convenience_Store = mean(SHR_Convenience_Store, na.rm =T) * 100, 
						SHR_Discount_Store = mean(SHR_Discount_Store, na.rm =T) * 100, 
						SHR_Dollar_Store = mean(SHR_Dollar_Store, na.rm =T) * 100, 
						SHR_Drug_Store = mean(SHR_Drug_Store, na.rm =T) * 100, 
						SHR_Grocery = mean(SHR_Grocery, na.rm =T) * 100, 
						SHR_Warehouse_Club = mean(SHR_Warehouse_Club, na.rm =T) * 100), 
						by = list(first_income, household_code)]

ggtmp		<- melt(tmp.data, id.vars = c("first_income","household_code"))
ggtmp$variable	<- factor(ggtmp$variable, levels = levels(ggtmp$variable), 
						labels = c("Biweek Exp ($)", "Price tier (Rank)", "Size tier (Rank)", 
								"Conveience Store (%)", "Discount Store (%)", "Dollar Store (%)", 
								"Drug Store (%)", "Grocery Store (%)", "Warehouse Club (%)"))
tmp.tab		<- ggtmp[, list(value = mean(value)), by = list(first_income, variable)]
tmp.tab		<- dcast(tmp.tab, variable ~ first_income, value.var = "value")
cat("The average level prior to the recession:\n"); print(tmp.tab); cat("\n")

# Plot the distribution by income group
if(make_plot){
	pdf(paste(plot.wd, "/graph_sbst_prior_dist.pdf", sep=""), width = ww, height = ww*ar)
	# print(ggplot(ggtmp, aes(first_income, value)) + 
	# 		stat_summary(fun.data = "mean_cl_boot", geom = "errorbar") + 
	# 		stat_summary(fun.y = "mean", geom = "point") + 
	# 		facet_wrap(~ variable, scales = "free_y") + 
	# 		labs(x = "Income Group", title = "Distribution by income group prior to the recession")
	# 		)
	
	print(ggplot(subset(ggtmp,variable %in% levels(ggtmp$variable)[1:3]), aes(first_income, value)) + 
			stat_summary(fun.data = "mean_cl_boot", geom = "errorbar") + 
			stat_summary(fun.y = "mean", geom = "point") + 
			facet_wrap(~ variable, scales = "free_y") + 
			labs(x = "Income Group")
			)
	print(ggplot(subset(ggtmp,variable %in% levels(ggtmp$variable)[4:9]), aes(first_income, value)) + 
			stat_summary(fun.data = "mean_cl_boot", geom = "errorbar") + 
			stat_summary(fun.y = "mean", geom = "point") + 
			facet_wrap(~ variable, scales = "free_y") + 
			labs(x = "Income Group")
			)		
	dev.off()
}

#####################################
# Benefits and Cost of substitution #
#####################################
my.fml	<- y ~ x:first_income + famsize + factor(year) + factor(month)
x.vec	<- c("ln_dol", "PCTL_UPRICE", "PCTL_SIZE", "SHR_Discount_Store", "SHR_Dollar_Store", "SHR_Grocery", "SHR_Warehouse_Club")
names(x.vec)	<- c("Exp", rep("Price/size", 2), rep("Share", 4))
y.vec	<- c("num_mod", "quant", "dol_rt")
names(y.vec)	<- c("Num. categories", "Basket quantity", "Exp/Income")

coef.df	<- data.frame()
for(i in 1:length(x.vec)){
	for(j in 1:length(y.vec)){
		tmp.fml	<- do.call("substitute", list(my.fml, setNames(list(as.name(y.vec[j]), as.name(x.vec[i])), c("y", "x")) ))
		tmp.fit	<- lm(tmp.fml, data = subset(hh_exp, dol > 0))
		b		<- coef(summary(tmp.fit))
		sel		<- grep(x.vec[i], rownames(b))
		tmp		<- data.frame(b[sel,], IV = x.vec[i], DV = y.vec[j], Var = rownames(b)[sel])
		rownames(tmp)	<- NULL
		coef.df	<- rbind(coef.df, tmp)
	}
	cat("Last formular: "); print(tmp.fml); cat("\n")
}

# Plot the parameters 
ggtmp			<- coef.df
ggtmp$IncGrp	<- changelab(ggtmp$Var, c("Low", "Med", "High"), c("Low", "Med", "High"))
ggtmp$Seg		<- changelab(ggtmp$IV, x.vec, names(x.vec))
ggtmp$Retail	<- ifelse(ggtmp$Seg == "Share", substr(ggtmp$IV, 5, nchar(as.character(ggtmp$IV))), as.character(ggtmp$IV))
ggtmp$DV		<- factor(ggtmp$DV, levels = y.vec, labels = names(y.vec))

plots 	<- list(NULL)
plots[[1]]	<- ggplot(subset(ggtmp, Seg == levels(ggtmp$Seg)[1]), aes(IncGrp, Estimate, fill = IncGrp)) + 
				geom_bar(stat = "identity", position = position_dodge(width=0.9)) + 
				geom_errorbar(aes(ymin = Estimate-1.96*Std..Error, ymax = Estimate + 1.96 * Std..Error), position=position_dodge(width=0.9), width=0.25) + 
				scale_fill_manual(values = incg.col) + 
				xlab("") + 
				facet_wrap(~ DV, scales = "free")
for(i in 2:3){
	plots[[i]]	<- ggplot(subset(ggtmp, Seg == levels(ggtmp$Seg)[i]), aes(Retail, Estimate, fill = IncGrp)) + 
					geom_bar(stat = "identity", position = position_dodge(width=0.9)) + 
					geom_errorbar(aes(ymin = Estimate-1.96*Std..Error, ymax = Estimate + 1.96 * Std..Error), position=position_dodge(width=0.9), width=0.25) + 
					facet_wrap(~ DV, scales = "free") + 
					theme(axis.text.x = element_text(angle = 45, hjust = .8)) + xlab("") + 
					scale_fill_manual(values = incg.col)
}


if(make_plot){
	pdf(paste(plot.wd, "/graph_sbst_cost.pdf",sep=""), width = ww, height = ww*ar)
	print(plots[[1]])
	print(plots[[2]])
	print(plots[[3]])
	dev.off()
}

######################################################
# Within-household substitution during the recession # 
######################################################
my.fml	<- y ~ first_income:ln_income + factor(month)

coef.within	<- data.frame()
for(i in 1:length(x.vec)){
	pct		<- proc.time()
	tmp.fml	<- as.formula(do.call("substitute", list(my.fml, setNames(list(as.name(x.vec[i])), c("y")) )) )
	tmp.fit	<- plm(tmp.fml, data = subset(hh_exp, dol>0 & !is.na(quant)), index = c("household_code", "biweek"), model = "within")
	b		<- coef(summary(tmp.fit))
	sel		<- grep("ln_income", rownames(b))
	tmp		<- data.frame(b[sel,], DV = x.vec[i], Var = rownames(b)[sel])
	rownames(tmp)	<- NULL
	coef.within	<- rbind(coef.within, tmp)
	use.time	<- proc.time() - pct
	cat("FE regressions with", x.vec[i], "finishes using", use.time[3]/60, "min.\n")
}

# Plot the parameters 
ggtmp			<- coef.within
ggtmp$IncGrp	<- changelab(ggtmp$Var, c("Low", "Med", "High"), c("Low", "Med", "High"))
ggtmp$Seg		<- changelab(ggtmp$DV, x.vec, names(x.vec))
ggtmp$Retail	<- substr(ggtmp$DV, 5, nchar(as.character(ggtmp$DV)))
tmplab			<- c("log(Exp)", "Price tier", "Size", "Discount_Store", "Dollar_Store", "Grocery", "Warehouse_Club")
ggtmp$DV		<- factor(ggtmp$DV, levels = x.vec, labels = tmplab)

plots 	<- list(NULL)
for(i in 1:2){
	plots[[i]]	<- ggplot(subset(ggtmp, Seg == levels(ggtmp$Seg)[i]), aes(DV, Estimate, fill = IncGrp)) + 
					geom_bar(stat = "identity", position = position_dodge(width=0.9)) + 
					geom_errorbar(aes(ymin = Estimate-1.96*Std..Error, ymax = Estimate + 1.96 * Std..Error), position=position_dodge(width=0.9), width=0.25) + 
					scale_fill_manual(values = incg.col) + 
					xlab("") 
}
plots[[3]]	<- ggplot(subset(ggtmp, Seg == "Share"), aes(Retail, Estimate, fill = IncGrp)) + 
				geom_bar(stat = "identity", position = position_dodge(width=0.9)) + 
				geom_errorbar(aes(ymin = Estimate-1.96*Std..Error, ymax = Estimate + 1.96 * Std..Error), position=position_dodge(width=0.9), width=0.25) + 
				theme(axis.text.x = element_text(angle = 45, hjust = .8)) + xlab("") + 
				scale_fill_manual(values = incg.col)

if(make_plot){
	pdf(paste(plot.wd, "/graph_sbst_income_effect.pdf",sep=""), width = ww, height = ww*ar)
	print(plots[[1]])
	print(plots[[2]])
	print(plots[[3]])
	dev.off()
}

save.image(paste(plot.wd, "/substitution_cost.rdata", sep=""))

cat("This program is done.\n")

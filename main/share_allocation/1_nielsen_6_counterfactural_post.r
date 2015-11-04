library(ggplot2)
library(reshape2)
library(Rcpp)
library(RcppArmadillo)
library(maxLik)
library(evd)
library(data.table)
library(scales)


setwd("~/Documents/Research/Store switching/processed data")
source("../Exercise/main/outreg function.R")
plot.wd		<- "/Users/chaoqunchen/Desktop"
ww			<- 6.5
ar			<- .8

outfile 	<- paste(plot.wd, "/counter1.xlsx", sep="")
mywb		<- createWorkbook()
sht1		<- createSheet(mywb, "Counterfactual1")
make_plot	<- TRUE

# Read estimation from all segments #
nseg	<- 3
sim.ls	<- vector("list", length = nseg)
my.simdata	<- data.frame()
names(sim.ls)	<- c("Low", "Med", "High")
for(ii in 1:nseg){
	load(paste("Estimation/counter_v2_seg",ii,".rdata",sep=""))
	sel			<- tmp[as.character(sim.data.all$income)]
	sim.ls[[ii]] <- list(	Base07 	= sim.base07[sel,], 
							Base08 	= sim.base08[sel,], 
							FixExp	= sim.fixy08[sel,] )
	my.simdata	<- rbind(my.simdata, data.frame(sim.data.all, IncGrp = names(sim.ls)[ii]))
	rm(list = c("sim.base07", "sim.base08", "sim.fixy08", "sim.fixp08"))						
}

# Decompose the effects 
tmp.tab	<- data.frame(NULL)
tmpsim	<- list(NULL)
for(i in 1:nseg){
	tmp		<- rbind(Baseline = colMeans(sim.ls[[i]]$Base07), Base08 = colMeans(sim.ls[[i]]$Base08), 
					FixExp = colMeans(sim.ls[[i]]$FixExp))
	tmp1	<- rbind(Baseline = tmp["Baseline",], PreferenceEff = tmp["FixExp",] - tmp["Baseline",], 
					 ExpEff	  = tmp["Base08",] - tmp["FixExp",], IncomeDrop = tmp["Base08",])

	tmp.tab	<- rbind(tmp.tab, tmp1)
	if(i == 1){
		tmpsim	<- sim.ls[[i]]
	}else{
		tmpsim	<- lapply(1:length(tmpsim), function(x) rbind(tmpsim[[x]], sim.ls[[i]][[x]] ))
	}
}
tmp.tab$Var		<- gsub("[0-9]", "", rownames(tmp.tab))
names(tmpsim)	<- names(sim.ls[[1]])
tmp		<- rbind(Baseline = colMeans(tmpsim$Base07), Base08 = colMeans(tmpsim$Base08), 
				FixExp = colMeans(tmpsim$FixExp))
tmp1	<- rbind(Baseline = tmp["Baseline",], PreferenceEff = tmp["FixExp",] - tmp["Baseline",], 
				 ExpEff	  = tmp["Base08",] - tmp["FixExp",], IncomeDrop = tmp["Base08",])
tmp1	<- data.frame(as.matrix(tmp1), Var = rownames(tmp1))
colnames(tmp1)	<- colnames(tmp.tab)

sel		<- sapply(1:ncol(tmp.tab), function(x) is.numeric(tmp.tab[,x]))
xlsx.addHeader(mywb, sht1, value = "Difference of avearge by income group ", level = 2)						
xlsx.addTable(mywb, sht1, cbind(Var = tmp.tab[,!sel], round(tmp.tab[,sel], 4)), row.names = F)
xlsx.addHeader(mywb, sht1, value = "Difference of avearge of overall effect across all households", level = 2)						
xlsx.addTable(mywb, sht1, cbind(Var = tmp1[,!sel], round(tmp1[,sel], 4)), row.names = F)

saveWorkbook(mywb, outfile)

# Plot the decomposed effeccts 
ggtmp	<- rbind(tmp.tab, tmp1)
ggtmp$IncGrp	<- rep(c("Low", "Med", "High", "Overall"), each = nrow(tmp1))
ggtmp$IncGrp	<- factor(ggtmp$IncGrp, levels = c("Low", "Med", "High", "Overall"))
ggtmp	<- subset(ggtmp, Var %in% c("PreferenceEff", "ExpEff"))
ggtmp	<- melt(ggtmp[,-(1:3)], id.vars = c("IncGrp","Var"))

ggplot(ggtmp, aes(variable, value, fill = Var)) + geom_bar(stat = "identity", position = "dodge") + 
		facet_grid(IncGrp ~ .)

# Plot the distribution of difference
tmp1	<- data.frame(tmpsim$FixExp - tmpsim$Base07, IncGrp = my.simdata$IncGrp)	# Preference effect 
tmp2	<- data.frame(tmpsim$Base08 - tmpsim$FixExp, IncGrp = my.simdata$IncGrp)	# Expenditure effect
ggtmp	<- rbind( data.frame( melt(tmp1[,-(1:3)], id.vars = "IncGrp"), Effect = "Preference"), 
				data.frame( melt(tmp2[,-(1:3)], id.vars = "IncGrp"), Effect = "Expenditure"))

qt50 <- function(x){
	out <- quantile(x, probs = c(0.25,0.75))
	names(out) <- c("ymin","ymax")
	return(out) 
}

pdf(paste(plot.wd, "/graph_counter1.pdf", sep=""), width = ww, height = ww*ar)
ggplot(ggtmp, aes(variable, value, col = Effect, fill = Effect, alpha = .5)) + 
		stat_summary(fun.data = qt50, geom = "errorbar", position = position_dodge(width=0.8), width = .2) +
		stat_summary(fun.y = mean, geom = "bar", position = position_dodge(width=0.8), width=0.2) +
		facet_grid(IncGrp ~ ., margins = TRUE) + 
		coord_cartesian(ylim = c(-.01, .01)) + 
		guides(alpha = FALSE) + 
		labs(x = "Retail", y = "Share difference") + 
		scale_y_continuous(labels = percent) +
		theme(axis.text.x = element_text(angle = 60, vjust = .8))
dev.off()		
		
		
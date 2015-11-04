library(ggplot2)
library(reshape2)
library(Rcpp)
library(RcppArmadillo)
library(maxLik)
library(evd)
library(data.table)
library(scales)
library(r2excel)

setwd("~/Documents/Research/Store switching/processed data")
source("../Exercise/main/outreg function.R")
plot.wd		<- "/Users/chaoqunchen/Desktop"
ww			<- 6.5
ar			<- .8

outfile 	<- paste(plot.wd, "/counterfact.xlsx", sep="")
mywb		<- createWorkbook()
sht1		<- createSheet(mywb, "Counterfactual1")
make_plot	<- TRUE

qt50 <- function(x){
	out <- quantile(x, probs = c(0.25,0.75))
	names(out) <- c("ymin","ymax")
	return(out) 
}

# Read estimation from all segments #
nseg	<- 3
sim.ls	<- vector("list", length = nseg)
my.simdata	<- data.frame()
names(sim.ls)	<- c("Low", "Med", "High")
sim.ls1	<- setNames(vector("list", length = nseg), c("Low", "Med", "High"))
for(ii in 1:nseg){
	load(paste("Estimation/ctrfact_seg",ii,".rdata",sep=""))
	tmp			<- 1:nrow(sim.unq)
	names(tmp)	<- sim.unq$income2007
	sel			<- tmp[as.character(sim.data$income)]
	sim.ls[[ii]] <- list(Base07 = sim.base07[sel,], 
						Base08 	= sim.base08[sel,], Price = sim.prc[sel,])
	names(sim.X.ls)	<- selcol
	sim.X.ls		<- lapply(sim.X.ls, function(x) x[sel,])
	sim.ls[[ii]]	<- c(sim.ls[[ii]], sim.X.ls)						
	my.simdata	<- rbind(my.simdata, data.frame(sim.data, IncGrp = names(sim.ls)[ii]))
	sim.ls1[[ii]]	<- tmp.tab
	rm(list = c("tmp.tab", "sim.prc", "sim.X.ls","sim.base07", "sim.base08"))						
}

# Income effect 
grp.wt	<- table(my.simdata$IncGrp)
tmp.tab	<- do.call(cbind, lapply(sim.ls1, function(x) x[,c("Base07","Base08")]))
tmp1	<- apply(tmp.tab[,seq(1, 6, 2)], 1, function(x) sum(x*grp.wt)/sum(grp.wt) )		# Overall effect
tmp2	<- apply(tmp.tab[,seq(2, 6, 2)], 1, function(x) sum(x*grp.wt)/sum(grp.wt) )
tmp.tab	<- cbind(tmp.tab, Base07 = tmp1, Base08 = tmp2)
sel		<- setdiff(1:nrow(tmp.tab), grep("omega",rownames(tmp.tab)))
tmp.tab	<- tmp.tab[sel,]
tmp.tab[-1,]	<- tmp.tab[-1,] * 100				# Percentage as unit

xlsx.addHeader(mywb, sht1, value = "Income effect by income group", level = 2)	
xlsx.addTable(mywb, sht1, round(tmp.tab, 2))

# The effectiveness of each strategy. 
selcol1	<- c(selcol, "Price")
ggtmp	<- data.frame()
for(i in 1:nseg){
	for(j in 1:length(selcol1)){
		tmp	<- sim.ls[[i]][[selcol1[j]]] - sim.ls[[i]][["Base08"]]
		tmp	<- data.frame(tmp, Var = selcol1[j], IncGrp = names(sim.ls)[i])
		ggtmp	<- rbind(ggtmp, tmp)
	}
}
ggtmp	<- ggtmp[,setdiff(colnames(ggtmp), c("omega","spl.omega"))]
ggtmp1	<- melt(ggtmp, id.vars = c("IncGrp", "Var"))
ggtmp1$Var	<- factor(ggtmp1$Var, levels = selcol1, 
						labels = c("Size", "Num. UPC/module", "Num. module", "Private label", "Price"))
sel		<- ggtmp1$variable != "Exp"
ggtmp1[sel,"value"] <- ggtmp1[sel,"value"] * 100		# Change expenditure share to percentage 

var.sec	<- list(Expenditure = "Exp", FocalRetail = "Grocery", OtherRetail = gsub("\\s", ".",setdiff(fmt_name, "Grocery")))						
plots	<- setNames(vector("list", length = length(var.sec)), names(var.sec))
for(i in 1:length(var.sec)){
	tmplab		<- ifelse(i==1, "Expenditure change ($)", "Share change (%)")
	plots[[i]]	<- ggplot(subset(ggtmp1, variable %in% var.sec[[i]]), aes(IncGrp, value)) + 
					stat_summary(fun.data = qt50, geom = "errorbar", width = .3 ) + 
					stat_summary(fun.y = mean, geom = "point") + 
					geom_hline(yintercept = 0, size = .25, linetype = 2) + 
					facet_grid(variable ~ Var, scales = "free" ) + 
					labs(x = "Income group", y = tmplab)
}

saveWorkbook(mywb, outfile)
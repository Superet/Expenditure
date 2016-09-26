library(ggplot2)
library(reshape2)
library(maxLik)
library(evd)
library(data.table)
library(doParallel)
library(foreach)
# library(chebpol)
library(nloptr)
library(mgcv)
options(error = quote({dump.frames(to.file = TRUE)}))

# setwd("~/Documents/Research/Store switching/Processed_data/processed data_20160809")
# plot.wd	<- '~/Desktop'
# sourceCpp(paste("../../Exercise/Multiple_discrete_continuous_model/", model_name, ".cpp", sep=""))
# source("../../Exercise/Multiple_discrete_continuous_model/0_Efficient_Allocation_function.R")
# source("../../Exercise/main/ctrfact_sim_functions.r")

# setwd("/home/brgordon/ccv103/Exercise/run")
# setwd("/kellogg/users/marketing/2661703/Exercise/run")
# setwd("/sscc/home/c/ccv103/Exercise/run")
setwd("U:/Users/ccv103/Documents/Research/Store switching/run")

run_id		<- 2
plot.wd		<- getwd()
make_plot	<- TRUE
ww			<- 10
ar			<- .6

# Set plotting parameters
annual.month	<- 12
annual		<- FALSE
draw.par		<- FALSE
sim.omega		<- FALSE

# Load estimation data 
cpi.adj		<- TRUE
setwd(paste("estrun_", run_id, sep=""))

###########################
# Read counterfactal data # 
###########################
# Read estimation from all segments #
nseg		<- 3
segname		<- c("Low", "Med", "High")
my.simdata	<- data.frame()
numsim		<- 500
ver.date1	<- "2016-09-20"	# "2016-08-29"	#"2016-08-24"  
base.08		<- sim08	<- data.frame()
sim.prc		<- data.frame()
attr.data	<- data.frame()
tmpls		<- ls()

for(ii in 1:nseg){
	# Load data 
	fname			<- paste("newfmt_seg",ii, sep="")
	fname			<- paste(fname, "_sim", numsim, "_", ver.date1, sep="")
	load(paste(fname, ".rdata",sep=""))
	
	# Compare the difference 
	tmp1	<- colSums(cbind(sim.base08$Average[,fmt_name],0), na.rm=T)
	tmp2	<- colSums(newf.sim[[1]][[1]][,fmt_name1], na.rm=T)
	mkt.dif	<- setNames(tmp2/sum(tmp2)-tmp1/sum(tmp1), fmt_name1)
	cat("Difference of market share:\n"); print(round(mkt.dif*100,2)); cat("\n")

	tmp2	<- colSums(newf.sim[[2]][[1]][,fmt_name1], na.rm=T)
	mkt.dif	<- setNames(tmp2/sum(tmp2)-tmp1/sum(tmp1), fmt_name1)
	cat("Difference of market share in scenario", sim.scn[2], "\n"); print(round(mkt.dif*100,2)); cat("\n")
	
	# 
	tmp1	<- sim.base08$Average
	tmp2	<- data.frame(do.call(rbind, newf.sim[[1]]), scn = sim.scn[1])
	tmp3	<- data.frame(do.call(rbind, newf.sim[[2]]), scn = sim.scn[2])
	
	# If we want to convert the biweekly simulation to annual level
	if(annual){
		sel	<- c("Exp", fmt_name)
		tmp1[,sel]	<- tmp1[,sel] * annual.month
		tmp2[,sel]	<- tmp2[,c(sel,"New")] * annual.month
		tmp3[,sel]	<- tmp3[,c(sel,"New")] * annual.month
	}
	base.08	<- rbind(base.08, tmp1)	
	sim08	<- rbind(sim08, tmp2, tmp3)		
		
	my.simdata	<- rbind(my.simdata, data.frame(sim.unq, IncGrp = segname[ii]))
	attr.data	<- rbind(attr.data, attr_mat)
		
	# Delete unused objects
	rm(list = setdiff(ls(), c(tmpls, "tmpls", "ii", "fmt_name", "R", "iv.change", "sim.scn")))
	print(ii)						
}

# Overall market share
bar.col	<- "#4F81BD"

fmt_name1	<- c(fmt_name, "New")
fmt_lab1	<- c(fmt_name, "Discount Store Express")
fmt_lab1[fmt_lab1=="Grocery"]	<- "Grocery Store"
tmp0	<- colSums(cbind(base.08[,fmt_name],0), na.rm=T)
tmp0	<- tmp0/sum(tmp0)
tmpls	<- split(sim08, sim08$scn)
tmp1	<- lapply(tmpls, function(x) colSums(x[x$d==0,gsub("\\s", ".",fmt_name1)], na.rm=T))
tmp1	<- lapply(tmp1, function(x) x/sum(x))
ggtmp	<- sapply(tmp1, function(x) x-tmp0)
cat("Difference of market share:\n"); print(round(ggtmp*100,2)); cat("\n")
ord 	<- order(ggtmp[,1])
ggtmp	<- melt(ggtmp)
ggtmp$retailer <- factor(ggtmp$Var1, levels = gsub("\\s", ".",fmt_name1)[ord], labels = fmt_lab1[ord])
text.col	<- rep("black", R+1)
text.col[c(2,7)]	<- bar.col
ggtmp$colind	<- ifelse(ggtmp$retailer %in% c("Discount Store Express", "Discount Store"), 1, 0) 

plots	<- list(NULL)
plots[[1]]	<- ggplot(subset(ggtmp, Var2 == "Convenience Store"), aes(retailer, value)) + 
				geom_bar(stat = "identity") + 
				scale_y_continuous(labels = scales::percent) +
				xlab("") + ylab("Percentage point change")+ coord_flip() + 
				theme_bw() + theme(legend.position = "none")

plots[[2]]	<- ggplot(subset(ggtmp, Var2 == "Convenience Store"), aes(retailer, value, fill = factor(colind))) + 
				geom_bar(stat = "identity") + 
				scale_y_continuous(labels = scales::percent) +
				scale_fill_manual(values = c("grey70", bar.col)) + 
				xlab("") + ylab("Percentage point change")+ coord_flip() + 
				theme_bw() + theme(legend.position = "none", axis.text.y = element_text(colour=text.col))

plots[[3]]	<- ggplot(subset(ggtmp, Var2 == "Convenience Store"), aes(retailer, value, fill = factor(colind), alpha = factor(colind))) + 
				geom_bar(stat = "identity", fill = bar.col) + 
				scale_y_continuous(labels = scales::percent) +
				scale_fill_manual(values = c("grey70", bar.col)) +
				scale_alpha_manual(values = c(0, 1), guide='none') + 
				xlab("") + ylab("Percentage point change")+ coord_flip() + 
				theme_bw() + theme(legend.position = "none", axis.text.y = element_text(colour=text.col))

pdf(paste(plot.wd, "/graph_newformat.pdf",sep=""), width =4.5, height = 4.5*ar)
for(i in 1:length(plots)){ print(plots[[i]])}
dev.off()

# Price range
sel		<- sim08$scn == "Convenience Store"
tmp		<- split(sim08[sel,], sim08[sel, "d"])
tmp		<- lapply(tmp, function(x) colSums(x[,gsub("\\s", ".",fmt_name1)]))
tmp		<- lapply(tmp, function(x) x/sum(x))
ggtmp	<- sapply(tmp, function(x) x - tmp0)
ggtmp1	<- ggtmp["Discount.Store",] + ggtmp["New",]
ggtmp1	<- data.frame(d = as.numeric(names(ggtmp1)), Total_share = ggtmp1)

pdf(paste(plot.wd, "/graph_newformat_pricing.pdf",sep=""), width =4.5, height = 4.5*ar)
ggplot(ggtmp1, aes(d, Total_share)) + geom_line() + 
		scale_y_continuous(labels = scales::percent) + 
		labs(x = expression(delta), y = "Net share gain") + 
		theme_bw()
dev.off()

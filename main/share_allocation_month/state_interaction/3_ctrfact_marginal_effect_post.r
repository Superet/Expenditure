library(ggplot2)
library(reshape2)
library(Rcpp)
library(RcppArmadillo)
library(maxLik)
library(evd)
library(data.table)
library(scales)
library(gridExtra)
library(stargazer)

plot.wd		<- "/Users/chaoqunchen/Desktop"

# setwd("/home/brgordon/ccv103/Exercise/run")
run_id		<- 10
ww			<- 6.5
ar			<- .6
ww1			<- 8.5
ar1			<- .5
# setwd(paste("~/Documents/Research/Store switching/processed data/Estimation/estrun_", run_id, sep=""))
setwd(paste("U:/Users/ccv103/Documents/Research/Store switching/run/estrun_", run_id, sep=""))

# Set plotting parameters
annual.month	<- 12
annual		<- TRUE
draw.par		<- FALSE
sim.omega		<- FALSE

###########################
# Read counterfactal data # 
###########################
# Read estimation from all segments #
nseg	<- 3
segname	<- c("Low", "Med", "High")
sim.ls	<- setNames(vector("list", length = nseg), c("Low", "Med", "High"))
my.simdata	<- data.frame()
sim.df		<- data.frame()
numsim		<- 100
ver.date1	<- "2016-07-17"  
tmpls	<- ls()

for(ii in 1:nseg){
	# Load data 
	fname			<- paste("marginal_seg",ii, sep="")
	if(draw.par){ fname <- paste(fname, "_pardraw", sep="") }
	if(sim.omega){fname	<- paste(fname, "_simomega", sep="") }
	fname			<- paste(fname, "_sim", numsim, "_", ver.date1, sep="")
	load(paste(fname, ".rdata",sep=""))
	
	# Compute the shopping incidence for baseline simulation
	# tmp1	<- apply(sim.base07$Allocation, c(2,3), function(x) mean(x>.01))
	# colnames(tmp1)	<- paste("p_", fmt_name, sep="")
	# tmp2	<- data.frame(lnInc = as.numeric(rownames(tmp1)), sim.base07$Average, tmp1)
	# names(tmp2)	<- gsub("\\.", " ", names(tmp2))
	# tmp1	<- apply(sim.base08$Allocation, c(2,3), function(x) mean(x>.01))
	# colnames(tmp1)	<- paste("p_", fmt_name, sep="")
	# tmp3	<- data.frame(lnInc = as.numeric(rownames(tmp1)), sim.base08$Average, tmp1)
	# names(tmp3)	<- gsub("\\.", " ", names(tmp3))
	tmp2	<- data.frame(sim.unq[,c("income2007","Inc07","zip3")], sim.base07)
	tmp3	<- data.frame(sim.unq[,c("income2007","Inc07","zip3")], sim.base08)
	names(tmp2)	<- gsub("\\.", " ", names(tmp2))
	names(tmp3)	<- gsub("\\.", " ", names(tmp3))
		
	# Compute the shopping incidence for price simualtion in 2007
	sim08			<- do.call(rbind, lapply(sim.prc.ls08, function(x) cbind(sim.unq[,c("income2007","zip3")],x$Average)))
	names(sim08)	<- gsub("\\.", " ", names(sim08))
	# tmp1	<- melt(lapply(sim.prc.ls08, function(x) apply(x$Allocation, c(2,3), function(xx) mean(xx>0.01))) )
	# tmp1	<- dcast(tmp1, L1+lnInc ~ retailer, value.var = "value")
	# names(tmp1)	<- c("retailer", "change", "lnInc", paste("p_", names(tmp1)[-(1:3)], sep=""))
	# tmp1$retailer	<- fmt_name[tmp1$retailer]
	# sim08		<- merge(sim08, tmp1, by = c("retailer", "lnInc"))
				
	# If we want to convert the biweekly simulation to annual level
	if(annual){
		sel	<- c("Exp", fmt_name)
		tmp2[,sel]	<- tmp2[,sel] * annual.month
		tmp3[,sel]	<- tmp3[,sel] * annual.month
		sim08[,sel]	<- sim08[,sel] * annual.month
	}
		
	sim.ls[[ii]]<- list(Base07 = tmp2, Base08 = tmp3, price08 = sim08)
						
	my.simdata	<- rbind(my.simdata, data.frame(sim.data, IncGrp = names(sim.ls)[ii]))
	
	# Delete unused objects
	rm(list = setdiff(ls(), c(tmpls, "tmpls", "ii", "fmt_name", "R", "iv.change.vec")))
	print(ii)						
}


#############
# Functions #
#############
weight_summary	<- function(x, freq, fun, ...){
	x.full	<- unlist(lapply(1:length(x), function(i) rep(x[i], freq[i])))
	out		<- fun(x.full, ...)
	return(out)
}

mean_qt <- function(x, alpha = .025){
	out	<- c(mean(x), quantile(x, c(alpha, 1-alpha)))
	names(out)	<- c("y", "ymin", "ymax")
	return(out)
}

weight_se	<- function(se, freq){
# This function computes the standard error of x.bar = sum_i (freq_i*x_i)/sum_i freq_i
	wt	<- freq/sum(freq)
	v	<- sum(wt^2 * se^2, na.rm =T)
	return(sqrt(v))
}

mytrimfun	<- function (x, alpha = 0.05){
    if (alpha == 0) {
        return(x)
    }else {
        n <- length(x)
        lo <- floor(n * alpha) + 1
        hi <- n + 1 - lo
		if(sum(is.na(x)) >= (lo-1)){
			out	<- x
		}else{
			out <- sort.int(x, partial = unique(c(lo, hi)), na.last = FALSE)[lo:hi]
		}      
        return(out)
    }
}

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

compare_fn <- function(dat1, dat2, byvar, weight.dat, weight.by, fml, f = "mean", 
				val1 = "value", val2 = "value", agg.first = TRUE, y.name = "y"){
	if(nrow(dat1)!=nrow(dat2)){
		stop("Error: nrow(dat1) not equal to nrow(dat2).")
	}
	if(!identical(dat1[,byvar], dat2[, byvar])){
		stop("Error: two data are not sorted equally.")
	}
	if(agg.first){
		tmp	<- paste("list(", paste(byvar, collapse=","), ")")
		x1	<- data.table(merge(dat1, weight.dat, by.var = weight.by, all.x = T))
		x2	<- data.table(merge(dat2, weight.dat, by.var = weight.by, all.x = T))
		if(f=="mean"){
			x1	<- x1[,list(x1 = weight_summary(get(val1), nhh, mean)), by = eval(parse(text = tmp))]
			x2	<- x2[,list(x2 = weight_summary(get(val2), nhh, mean)), by = eval(parse(text = tmp))]
		}else if(f=="share"){
			x1	<- x1[,list(x1 = weight_summary(get(val1), nhh, sum)), by = eval(parse(text = tmp))]
			x1	<- x1[,x1:=x1/sum(x1)]
			x2	<- x2[,list(x2 = weight_summary(get(val2), nhh, sum)), by = eval(parse(text = tmp))]
			x2	<- x2[,x2:=x2/sum(x2)]
		}else if(f == "sum"){
			x1	<- x1[,list(x1 = weight_summary(get(val1), nhh, sum)), by = eval(parse(text = tmp))]
			x2	<- x2[,list(x2 = weight_summary(get(val2), nhh, sum)), by = eval(parse(text = tmp))]
		}
		out	<- merge(x1, x2, by = byvar)
		out	<- out[, y:=eval(parse(text = fml))]
		out	<- data.frame(out)
		out	<- out[,c(byvar, "y")]
		if(y.name != "y"){
			names(out)	<- gsub("y", y.name, names(out))
		}
	}else{
		x1	<- dat1[,val1]
		x2	<- dat2[,val2]
		y	<- eval(parse(text = fml))
		dat.sum 	<- data.frame(dat1[,unique(c(byvar, weight.by))], y = y)
		# names(dat.sum)	<- c(byvar, "y")
		dat.sum		<- merge(dat.sum, weight.dat, by.var = weight.by, all.x=T)
		dat.sum		<- data.table(dat.sum)
		tmp			<- paste("list(", paste(byvar, collapse=","), ")")
		out			<- dat.sum[,list(y = weight_summary(y, nhh, mean), y.sd = weight_summary(y, nhh, sd)), 
								by = eval(parse(text = tmp))]
		out			<- data.frame(out)
		if(y.name != "y"){
			names(out)	<- gsub("y", y.name, names(out))
		}
	}
	return(out)
}

#################
# Income effect #
#################
iv.change 	<- -.1
p.idx		<- 0
plots		<- list(NULL)
test.alpha	<- .95						# Significance level
crtv		<- qnorm(test.alpha)		# Critical value
ret.ord		<- c("Dollar Store", "Grocery", "Discount Store", "Drug Store", "Convenience Store", "Warehouse Club")
dodge.width	<- .5						# The width parameter in position_dodge

# Income frequency
tab.inc	<- data.table(my.simdata)
tab.inc	<- tab.inc[, list(nhh = length(unique(household_code))), by = list(IncGrp, income, zip3)]
setnames(tab.inc, "income", "income2007")
setkeyv(tab.inc, c("IncGrp", "income2007","zip3"))
# tab.inc$income2007	<- round(tab.inc$income2007, 2)

# Merge the weights into the baseline simulation
base07	<- do.call(rbind, lapply(1:nseg, function(i) cbind(IncGrp=names(sim.ls)[i], sim.ls[[i]]$Base07)))
base07	<- merge(base07, tab.inc, by = c("IncGrp", "income2007","zip3"), all.x = T)
base08	<- do.call(rbind, lapply(1:nseg, function(i) cbind(IncGrp=names(sim.ls)[i], sim.ls[[i]]$Base08)))
base08	<- merge(base08, tab.inc, by = c("IncGrp", "income2007","zip3"), all.x = T)

# Check if data is ordered 
base07	<- base07[order(base07$IncGrp, base07$income2007,base07$zip3),]
base08	<- base08[order(base08$IncGrp, base08$income2007,base08$zip3),]
identical(base07$income2007, base08$income2007)

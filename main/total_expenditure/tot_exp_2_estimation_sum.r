cat("This program begins to run at", as.character(Sys.time()), ".\n")

library(ggplot2)
library(reshape2)
library(maxLik)
library(evd)
library(data.table)
library(mgcv)
library(stargazer)
library(lmtest)
library(scatterplot3d)
options(error = quote({dump.frames(to.file = TRUE)}))

run_id			<- 1
make_plot		<- TRUE
interp.method	<- "spline"			# "cheb"
trim.alpha		<- 0 #0.05
cpi.adj			<- TRUE
week.price		<- FALSE
ver.date		<- "2016-09-11"
nseg			<- 9


# setwd("~/Documents/Research/Store switching/Processed_data")
# plot.wd	<- '~/Desktop'
# source("../Exercise/Multiple_discrete_continuous_model/0_Efficient_Allocation_function.R")
# source("../Exercise/main/ctrfact_sim_functions.r")

# setwd("/home/brgordon/ccv103/Exercise/run")
# setwd("/kellogg/users/marketing/2661703/Exercise/run")
setwd("/sscc/home/c/ccv103/Exercise/run")
# setwd("U:/Users/ccv103/Documents/Research/Store switching/run")

# Read estimation from all segments # 
mydata.full		<- data.frame()
shr.est.ls		<- vector("list", length = nseg)
exp.est.ls		<- vector("list", length = nseg)
omega.ls		<- vector("list", length = nseg)
tmp.ls			<- ls()

for(ii in 1:nseg){
	fname	<- paste("estrun_",run_id,"/MDCEV_est_seg", ii, "_", ver.date,".rdata",sep="")
	cat(fname, "\n")
	load(fname)
	
	mydata.full	<- rbind(mydata.full, cbind(mydata, segment = ii))
	shr.est.ls[[ii]]	<- sol
	exp.est.ls[[ii]]	<- lsol	
	omega.ls[[ii]]		<- omega_deriv
	
	rm(list = setdiff(ls(), c(tmp.ls, "tmp.ls", "ii", "run_id", "hh_dist", "hh_month", "pan", "fmt_name", "R")))
}

# Summarize the parameters
tmp		<- vector("list", length = nseg)
for(i in 1:length(tmp)){
	tmp[[i]]	<- summary(shr.est.ls[[i]])$estimate
	class(tmp[[i]])	<- "coeftest"
}
stargazer(tmp, type = "text")
stargazer(tmp, type = "html", out = paste("estrun_", run_id,"/est.html", sep=""))

# Plot omega function
tmpdat	<- data.frame(y = seq(50, 500, 20))
selcol	<- c("distance_haver","size_index", "ln_upc_per_mod", "ln_num_module", "prvt_overall")
tmp		<- dcast(hh_dist[,c("household_code", "year", "channel_type","unitprice_paid")], 
				household_code + year ~ channel_type, value.var = "unitprice_paid")
tmp1	<- as.vector(colMeans(tmp[,fmt_name], na.rm=T))
tmp2	<- unlist(lapply(fmt_name, function(x) colMeans(hh_dist[hh_dist$channel_type==x, selcol], na.rm=T)))
tmpdat$price	<- rep(1, nrow(tmpdat)) %*% t(tmp1)
tmpdat$X<- rep(1, nrow(tmpdat)) %*% t(tmp2)

ggtmp	<- data.frame()
for(i in 1:nseg){
	ggtmp <- rbind(ggtmp, data.frame(u = omega.ls[[i]](tmpdat), y = tmpdat[,"y"], seg = i))
}
ggplot(ggtmp, aes(y, u)) + geom_line() + 
	facet_wrap(~seg, scales = "free_y")

# Save the estimation 
save(exp.est.ls, shr.est.ls, omega.ls, fmt_name, hh_dist, mydata.full, nseg, pan, R, run_id, selcol, 
	file = paste("estrun_",run_id,"/MDCEV_est_", ver.date,".rdata",sep=""))



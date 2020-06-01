cat("This program begins to run at", as.character(Sys.time()), ".\n")

library(ggplot2)
library(reshape2)
library(Rcpp)
library(RcppArmadillo)
library(maxLik)
library(evd)
library(data.table)
library(doParallel)
library(foreach)
library(nloptr)
library(mgcv)
options(error = quote({dump.frames(to.file = TRUE)}))

args <- commandArgs(trailingOnly = TRUE)
print(args)
if(length(args)>0){
    for(i in 1:length(args)){
      eval(parse(text=args[[i]]))
    }
}
fmt_name<- c("Convenience Store","Discount Store","Dollar Store","Drug Store", "Grocery", "Warehouse Club")
arrid	<- as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))
selret	<- fmt_name[ceiling(arrid/3)]
seg_id	<- ifelse(ceiling(arrid%%3)==0, 3, ceiling(arrid%%3))
cat("seg_id =",seg_id, ".\n")
cat("selret =", selret, "\n")

model_name 		<- "MDCEV_share_gamma_outside"
run_id			<- 6
make_plot		<- TRUE
trim.alpha		<- 0.05
# selret			<- "Grocery"

# setwd("~/Documents/Research/Store switching/Processed_data")
# plot.wd	<- '~/Desktop'
# sourceCpp(paste("../Exercise/Multiple_discrete_continuous_model/", model_name, ".cpp", sep=""))
# source("../Exercise/Multiple_discrete_continuous_model/0_Efficient_Allocation_function.R")
# source("../Exercise/main/ctrfact_sim_functions.r")

# setwd("/home/brgordon/ccv103/Exercise/run")
# setwd("/kellogg/users/marketing/2661703/Exercise/run")
setwd("/sscc/home/c/ccv103/Exercise/run")
# setwd("U:/Users/ccv103/Documents/Research/Store switching/run")

sourceCpp(paste(model_name, ".cpp", sep=""))
source("0_Efficient_Allocation_function.R")
source("ctrfact_sim_functions.r")

(fname	<- paste("estrun_",run_id,"/retcat_est_", selret, "_seg", seg_id, "_", Sys.Date(),".rdata",sep=""))

#################################
# Read data and subsetting data #
#################################
load("hh_month_exp_2016-09-26.rdata")

# Extract 5% random sample
# length(unique(hh_month$household_code))
# sel			<- sample(unique(hh_month$household_code), .05*length(unique(hh_month$household_code)) )
# hh_month_save	<- hh_month
# hh_month		<- subset(hh_month, household_code %in% sel)

# Retail formats 
(fmt_name 	<- as.character(sort(unique(hh_dist$channel_type))))
R			<- length(fmt_name)
hh_month$rec		<- 1*(hh_month$year >= 2008)
cat("Table of segments in the expenditure data:\n"); print(table(pan$segment)); cat("\n")
dpt.name	<- c("DG", "GM","NFG","RF", "HBC", "OTHER")
dpt.lab		<- c("DRY GROCERY","GENERAL MERCHANDISE", "NON-FOOD GROCERY", "REFRIGERATED FROZEN", "HEALTH BEAUTY CARE","PRODUCE OTHER")
nc			<- length(dpt.name)						

# Fill in the missing distance for the channels that are not available for households; 
tmp			<- data.table(hh_dist)
tmp			<- tmp[, list(nozip = all(is.na(panelist_zip_code))), by = list(household_code, year)]
tmp			<- tmp[nozip==TRUE,]
unique(tmp$household) %in% hh_month$household_code

max.dist	<- 100
sum(hh_dist$distance_haver >max.dist|is.na(hh_dist$distance_haver), na.rm=T)/nrow(hh_dist)
sel			<- hh_dist$distance_haver >max.dist|is.na(hh_dist$distance_haver)
hh_dist[sel,"distance_haver"]	<- max.dist

# ------------#
# Subset data #
# Subset the segment 
selcol		<- c("household_code", "year", "month", "rec", paste(dpt.name, ".", selret, sep=""))
tmp1		<- subset(pan, ceiling(as.numeric(segment)/3) == seg_id)
mydata		<- subset(hh_month[,selcol], household_code %in% tmp1$household_code)		

# Drop the years for which income is not matched
sel			<- with(mydata, paste(household_code, year, sep="*")) %in%
				with(pan, paste(household_code, panel_year, sep="*"))				
sum(sel)/length(sel)
mydata		<- mydata[sel,]			
mydata		<- merge(mydata, pan[,c("household_code", "panel_year", "income_mid","scantrack_market_descr")], 
				by.x = c("household_code", "year"), by.y = c("household_code", "panel_year"), all.x=T)
				
# Filter the households who did not have all the channel
tmp 		<- data.table(fmt_dpt[fmt_dpt$channel_type == selret,])
tmp 		<- tmp[,list(n = length(department)), by = list(scantrack_market_descr, year)]			
tmp			<- subset(tmp, year > 2003 & year < 2011)
sel			<- tmp$n == nc
sel			<- paste(mydata$scantrack_market_descr, mydata$year, sep="*") %in% 
			   paste(tmp[sel,scantrack_market_descr], tmp[sel,year], sep="*")	
mydata		<- mydata[sel,]			

# Only keep the positive expenditure 
mydata$dol	<- rowSums(as.matrix(mydata[,grep(selret, names(mydata))]))								
sum(mydata$dol < .1)
mydata.ful	<- mydata
sum(mydata$dol > .1)/nrow(mydata)
mydata		<- subset(mydata, dol > .1)							# Subset positive spending
cat("dim(mydata) =", dim(mydata), "\n")

################################
# Organize data for estimation #
################################
# The data required for estimation: 
# outcome variables: y (expenditure), shr (expenditure share);
# explanatory variables: price, X_list

# Match price
sel		<- fmt_dpt$channel_type == selret
tmp		<- fmt_dpt[sel,c("scantrack_market_descr", "year", "department","unitprice_paid")]
tmp$department	<- factor(tmp$department, levels = dpt.lab, labels = dpt.name)
tmp		<- dcast(tmp, scantrack_market_descr + year ~ department, value.var = "unitprice_paid")				
names(tmp)	<- c("scantrack_market_descr", "year", paste("PRC_", dpt.name, sep=""))
mydata	<- merge(mydata, tmp, by = c("scantrack_market_descr", "year"), all.x=T)
ord		<- order(mydata$household_code, mydata$year, mydata$month)
mydata	<- mydata[ord,]

# Outcome variables as matrix
beta0.base	<- which(dpt.name == "OTHER")
y		<- as.vector(mydata$dol)
shr		<- as.matrix(mydata[,paste(dpt.name,".", selret,sep="")])/y
price	<- cbind(as.matrix(mydata[,paste("PRC_", dpt.name[-beta0.base], sep="")]), PRC_OTHER = 1)
ln_inc	<- log(mydata$income_mid)
shr_sgn	<- 1 * (shr > 0)
M		<- rowSums(shr_sgn)
cat("Summary of expenditure share:\n"); print(summary(shr)); cat("\n")

# Match retailers' attributes
selcol	<- c("size_index", "ln_upc_per_mod", "ln_num_module", "prvt_overall")
fmt_dpt$ln_num_module	<- log(fmt_dpt$num_module)
fmt_dpt$ln_upc_per_mod 	<- log(fmt_dpt$num_upc_per_mod)
sel		<- fmt_dpt$channel_type == selret
attr_mat	<- merge(mydata[,c("household_code","scantrack_market_descr", "year","month")], fmt_dpt[sel,c("scantrack_market_descr", "year", "department", selcol)], 
				by = c("scantrack_market_descr", "year"), all.x=T)
attr_mat	<- attr_mat[order(attr_mat$household_code, attr_mat$year, attr_mat$month),]	
cat("dim(mydata) =", dim(mydata), ".\n")
cat("dim(attr_mat) =", dim(attr_mat), ", dim(mydata)*nc =",nrow(mydata)*nc, ".\n")
cat("corr(attr_mat):\n"); print(cor(attr_mat[,-(1:5)])); cat("\n")
summary(attr_mat[,-c(1:5)])
summary(price)

X1	<- vector("list", length=nc)
names(X1)	<- dpt.name
for(i in 1:(nc-1)){
	tmp0<- rep(0, nc-1)
	tmp0[i]	<- 1
	tmp1<- as.matrix(cbind(1, rec=mydata$rec))	
	tmp2<- kronecker(tmp1, t(tmp0))
	colnames(tmp2)	<- paste(rep(c("Intercept", "Rec"), each = nc-1), dpt.name[-beta0.base], sep=":")
	sel	<- attr_mat$department == dpt.lab[i]
	tmp	<- as.matrix(attr_mat[sel,selcol])
	X1[[i]]	<- cbind(tmp2, tmp, tmp*mydata$rec)
	colnames(X1[[i]])	<- c(colnames(tmp2), colnames(tmp), paste(colnames(tmp), ":Rec", sep=""))
}
X1[[nc]]	<- matrix(0, nrow(X1[[1]]), ncol(X1[[1]]))
colnames(X1[[nc]])	<- colnames(X1[[1]])
nx 		<- ncol(X1[[1]])
cat("nx =", nx,"\n")
rm(list = intersect(ls(), c("tmp1", "tmp2", "tmp3", "tmp", "tmp0","ord")))

##############
# Estimation #
##############
MDCEV_wrapper <- function(param){
	MDCEV_ll_fnC(param, nx=nx, shr = shr, y = y, p = price, X1=X1, X2=X1, shr_sgn=shr_sgn, M = M)
}

# Guess initial values
inits	<- NULL
inits2	<- NULL
pct		<- proc.time()
for(i in 1:(nc-1)){
	tmpy	<- 1 * (shr[,i]>0)
	tmpx	<- cbind(Rec = mydata$rec, X1[[i]][,((nx-length(selcol)*2+1):nx)])
	tmpfit	<- glm(y ~ ., data = data.frame(y = tmpy, tmpx), family=binomial(link='logit') )
	inits	<- cbind(inits, coef(tmpfit))
	
	sel		<- shr[,i] > 0
	tmpy	<- log(shr[,i]/price[,i] + 1) - log(shr[beta0.base] + 1)
	tmpfit	<- lm(tmpy[sel] ~ as.matrix(tmpx)[sel,])
	inits2	<- cbind(inits2, coef(tmpfit))
	print(i)
}
use.time	<- proc.time() - pct
cat("Initial logit finishes with", use.time[3]/60, "min.\n")
rm(list = intersect(ls(), c("tmpy", "tmpx", "tmpfit")))

# Normalize the coefficients for demographics relative to those with grocery 
beta.init	<- c(t(inits[c(1:2),]))
beta.init	<- c(beta.init, rowMeans(inits[-c(1:2),]))
gamma.init	<- c(t(inits2[c(1:2),]))
gamma.init	<- c(gamma.init, rowMeans(inits2[-c(1:2),]))
theta.init	<- list(c(beta.init/10, gamma.init/10, -.5))
for(i in 1:length(theta.init)){
	names(theta.init[[i]])	<- c(colnames(X1[[1]]), paste("gamma0_", colnames(X1[[1]]), sep=""), "ln_sigma")
}					
cat("Initial values are:\n"); print(theta.init); cat("\n")
system.time(tmp <- MDCEV_wrapper(theta.init[[1]]) )

max.method	<- "BFGS"	#"NM"
tmpsol	<- vector("list", length = length(theta.init))
for(i in 1:length(theta.init)){
	pct <- proc.time()
	tmpsol[[i]]		<- maxLik(MDCEV_wrapper, start=theta.init[[i]], method=max.method)
	use.time <- proc.time() - pct
	cat("MDCEV estimation finishes with", use.time[3]/60, "min.\n")
}
(tmp		<- sapply(tmpsol, function(x) x$maximum))
sel		<- which.max(tmp)
sol 	<- tmpsol[[sel]]
print(summary(sol))
cat("--------------------------------------------------------\n")

################################
# Simulate the inclusive value # 
################################
#-------------------------------------------- # 
# Compute the inclusive value with random draws
shr.par	<- coef(sol)
set.seed(687)
numsim 		<- 500
eps_draw	<- matrix(rgev(numsim*nc, scale = exp(shr.par["ln_sigma"])), numsim, nc)

# Subset simulation data
sel			<- duplicated(mydata[,c("household_code", "year")]) & mydata$dol <= quantile(mydata$dol, .95)
sht.X1		<- lapply(X1, function(x) x[!sel,])
sht.price	<- price[!sel,]
sht.y		<- y[!sel]
delta		<- .1

pct			<- proc.time()
omega_draw	<- foreach(i = 1:numsim, .packages= c("nloptr"), .combine = 'rbind', .errorhandling='remove') %do% {
	incl_value_fn(param_est = shr.par, nx, base=beta0_base, X1 = sht.X1, X2= sht.X1, y = sht.y, price = sht.price, nc, 
					eps_draw = rep(1,length(sht.y)) %*% t(eps_draw[i,]), returnmax = TRUE)$omega
}

omega_draw1	<- foreach(i = 1:numsim, .packages= c("nloptr"), .combine = 'rbind', .errorhandling='remove') %do% {
	incl_value_fn(param_est = shr.par, nx, base=beta0_base, X1 = sht.X1, X2= sht.X1, y = sht.y+delta, price = sht.price, R, 
					eps_draw = rep(1,length(sht.y)) %*% t(eps_draw[i,]), returnmax = TRUE)$omega
}
use.time	<- proc.time() - pct
cat("Simulation of omega with random draws finishes with", use.time[3]/60, "min.\n")
cat("--------------------------------------------------------\n")
print(summary(c(omega_draw)))

if(trim.alpha == 0){
	tmpdat	<- data.frame(omega = (colMeans(omega_draw1) - colMeans(omega_draw))/delta, y = sht.y, ln_inc = ln_inc[!sel])
}else{
	tmp1	<- apply(omega_draw, 2, mean, trim = trim.alpha, na.rm=T)
	tmp2	<- apply(omega_draw1, 2, mean, trim = trim.alpha, na.rm=T)
	tmpdat	<- data.frame(omega = (tmp2-tmp1)/delta, y = sht.y, ln_inc = ln_inc[!sel])
}
summary(tmpdat$omega)
tmpdat$price	<- sht.price[,-beta0.base]
tmpdat$X		<- as.matrix(do.call(cbind, lapply(sht.X1[1:(nc-1)], function(x) x[,selcol])))	
omega_deriv		<- spl2dfun(tmpdat, fml = omega ~ y + I(y^2) + price + X)

# Plot omega function
tmp		<- dcast(fmt_dpt[fmt_dpt$channel_type == selret, c("scantrack_market_descr", "year", "department","unitprice_paid")], 
				scantrack_market_descr + year ~ department, value.var = "unitprice_paid")
tmp1	<- as.vector(colMeans(tmp[,dpt.lab[-beta0.base]], na.rm=T))
tmp		<- attr_mat[,c("household_code", "year", "department", selcol)]
tmp2	<- unlist(lapply(dpt.lab[-beta0.base], function(x) colMeans(tmp[tmp$department==x, selcol], na.rm=T)))
(tmp		<- quantile(y, c(.1, .99)))
tmpdat	<- data.frame(y = seq(tmp[1], tmp[2], length = 100))
tmpdat$price	<- rep(1, nrow(tmpdat)) %*% t(tmp1)
tmpdat$X<- rep(1, nrow(tmpdat)) %*% t(tmp2)
ggtmp	<- data.frame(y = tmpdat$y, d1 = omega_deriv(tmpdat), d2 = omega_deriv(tmpdat, deriv = 1))
ggtmp	<- melt(ggtmp, id.vars = "y")

pdf(gsub(".rdata", ".pdf",fname), width = 6, height = 4.5)
ggplot(ggtmp, aes(y, value)) + geom_line() + facet_wrap(~variable, scales = "free")
dev.off()

rm(list = intersect(ls(), c("hh_month","args","tmpdat", "shr.x", "sht.y", "sht.X1", "sht.X2", "sht.price", "sel", "pct", "use.time", "i", 
			"Allocation_constr_fn","Allocation_fn","Allocation_nlop_fn","incl_value_fn", "omega.lm", "mytrimfun", "ggtmp","simExp_fn",
			"exp_fn", "expFOC_fn", "f", "SimWrapper_fn","solveExp_fn", "spl2dfun", "tmp", "tmp1", "tmp2", "uP_fn","uPGrad_fn", "mysplfun", 
			"inits", "inits2","MDCEV_ll_fnC","MDCEV_wrapper","model_name","tmpsol"))
	)

save.image(file = fname)

cat("This program is done. \n")


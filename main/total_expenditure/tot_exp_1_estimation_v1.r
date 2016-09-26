# This script estimates model. The moduels follow as: 
# 1. Segment households based on initial income level and subset data; 
# 2. Prepare estimation data; 
# 3. Estimate the conditional allocation model at different initial values; 
# 4. Simulate inclusive values; 
# 5. Estimate upper level expenditure model. 

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

seg_id	<- as.numeric(Sys.getenv("PBS_ARRAY_INDEX"))

model_name 		<- "MDCEV_share_gamma"
run_id			<- 1
make_plot		<- TRUE
interp.method	<- "spline"			# "cheb"
trim.alpha		<- 0.05
cpi.adj			<- TRUE
week.price		<- FALSE

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

(fname	<- paste("estrun_",run_id,"/MDCEV_est_seg", seg_id, "_", Sys.Date(),".rdata",sep=""))

#################################
# Read data and subsetting data #
#################################
load("hh_month_exp_2016-09-10.rdata")

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

# Fill in the missing distance for the channels that are not available for households; 
tmp			<- data.table(hh_dist)
tmp			<- tmp[, list(nozip = all(is.na(panelist_zip_code))), by = list(household_code, year)]
tmp			<- tmp[nozip==TRUE,]
unique(tmp$household) %in% hh_month$household_code

max.dist	<- 100
sum(hh_dist$distance_haver >max.dist|is.na(hh_dist$distance_haver), na.rm=T)/nrow(hh_dist)
sel			<- hh_dist$distance_haver >max.dist|is.na(hh_dist$distance_haver)
hh_dist[sel,"distance_haver"]	<- max.dist

# Subset data 
# Subset the segment 
selcol		<- c("household_code", "year", "month", "rec", paste("DOL_", gsub("\\s", "_", fmt_name), sep=""))
tmp1		<- subset(pan, as.numeric(segment) == seg_id)
mydata		<- subset(hh_month[,selcol], household_code %in% tmp1$household_code)					
# Filter the households who did not have all the channel
tmp			<- data.table(hh_dist)
tmp			<- tmp[,list(nchannel = length(channel_type)), by = list(household_code, year)]
sel			<- tmp$nchannel == R
sel			<- paste(mydata$household_code, mydata$year, sep="*") %in% paste(tmp[sel,household_code], tmp[sel,year], sep="*")	
mydata		<- mydata[sel,]
# Drop the years for which income is not matched
sel			<- with(mydata, paste(household_code, year, sep="*")) %in%
				with(pan, paste(household_code, panel_year, sep="*"))
sum(sel)
mydata		<- mydata[sel,]
mydata$dol	<- rowSums(as.matrix(mydata[,grep("DOL", names(mydata))]))								
sum(mydata$dol < .1)
mydata		<- subset(mydata, dol > .1)							# Subset positive spending

################################
# Organize data for estimation #
################################
# The data required for estimation: 
# outcome variables: y (expenditure), shr (expenditure share);
# explanatory variables: price, X_list

# Match price
tmp		<- dcast(hh_dist[,c("household_code", "year", "channel_type","unitprice_paid")], 
				household_code + year ~ channel_type, value.var = "unitprice_paid")
names(tmp)	<- c("household_code", "year", paste("PRC_", gsub("\\s", "_", fmt_name), sep=""))
mydata	<- merge(mydata, tmp, by = c("household_code", "year"), all.x=T)
mydata	<- merge(mydata, pan[,c("household_code", "panel_year", "income_mid")], 
				by.x = c("household_code", "year"), by.y = c("household_code", "panel_year"), all.x=T)
ord		<- order(mydata$household_code, mydata$year, mydata$month)
mydata	<- mydata[ord,]

# Outcome variables as matrix
y		<- as.vector(mydata$dol)
shr		<- as.matrix(mydata[,paste("DOL_", gsub("\\s", "_", fmt_name), sep="")])/y
price	<- as.matrix(mydata[,paste("PRC_", gsub("\\s", "_", fmt_name), sep="")])
ln_inc	<- log(mydata$income_mid)

# Select the index that has the positive expenditure
tmp 	<- price * 1 *(shr> 0); colnames(tmp)	<- NULL
s1_index<- c(apply(tmp, 1, which.max))

# Match retailers' attributes
selcol	<- c("distance_haver","size_index", "ln_upc_per_mod", "ln_num_module", "prvt_overall")
hh_dist$ln_num_module	<- log(hh_dist$num_module)
hh_dist$ln_upc_per_mod 	<- log(hh_dist$num_upc_per_mod)
beta0_base 	<- which(fmt_name == "Grocery")

# Normalize distance by zipcode
tmp		<- data.table(hh_dist[,c("household_code", "panelist_zip_code","year", "channel_type", selcol)])
tmp		<- tmp[, maxd:=max(distance_haver, na.rm=T), by = list(panelist_zip_code)]
tmp		<- tmp[, distance_haver := distance_haver/maxd]
tmp		<- data.frame(tmp)

attr_mat	<- merge(mydata[,c("household_code", "year", "month")], tmp[,c("household_code", "year", "channel_type", selcol)], 
				by = c("household_code", "year"), all.x=T)
attr_mat	<- attr_mat[order(attr_mat$household_code, attr_mat$year, attr_mat$month),]	
cat("dim(mydata) =", dim(mydata), ".\n")
cat("dim(attr_mat) =", dim(attr_mat), ", dim(mydata)*R =",nrow(mydata)*R, ".\n")
cat("corr(attr_mat):\n"); print(cor(attr_mat[,-(1:4)])); cat("\n")
summary(attr_mat[,-c(1:4)])
summary(price)

X2	<- X1 	<- vector("list", length=length(fmt_name))
for(i in 1:length(fmt_name)){
	tmp0<- rep(0, R-1)
	if(i<beta0_base){
		tmp0[i]	<- 1
	}else if(i > beta0_base){
		tmp0[(i-1)]	<- 1
	}
	tmp1<- as.matrix(cbind(1, rec=mydata$rec))	
	tmp2<- kronecker(tmp1, t(tmp0))
	colnames(tmp2)	<- paste(rep(c("Intercept", "Rec"), each = R-1), fmt_name[-beta0_base], sep=":")
	tmp3<- rep(0, R)
	tmp3[i]	<- 1
	tmp3<- kronecker(tmp1, t(tmp3))
	colnames(tmp3)	<- paste(rep(c("Intercept", "Rec"), each = R), fmt_name, sep=":")
	sel	<- attr_mat$channel_type == fmt_name[i]	
	tmp	<- as.matrix(attr_mat[sel,selcol])
	X1[[i]]	<- cbind(tmp2, tmp, tmp*mydata$rec)
	X2[[i]]	<- cbind(tmp3, tmp, tmp*mydata$rec)
	colnames(X1[[i]])	<- c(colnames(tmp2), colnames(tmp), paste(colnames(tmp), ":Rec", sep=""))
	colnames(X2[[i]])	<- c(colnames(tmp3), colnames(tmp), paste(colnames(tmp), ":Rec", sep=""))
}
nx 		<- ncol(X1[[1]])
cat("nx =", nx,"\n")

##############
# Estimation #
##############
MDCEV_wrapper <- function(param){
	MDCEV_ll_fnC(param, nx, shr, y, s1_index, price, X1, X2, beta0_base)$ll
}

# Guess initial values
inits	<- NULL
inits2	<- NULL
pct		<- proc.time()
for(i in 1:R){
	tmpy	<- 1 * (shr[,i]>0)
	tmpx	<- cbind(Rec = mydata$rec, X1[[i]][,((nx-length(selcol)*2+1):nx)])
	tmpfit	<- glm(y ~ ., data = data.frame(y = tmpy, tmpx), family=binomial(link='logit') )
	inits	<- cbind(inits, coef(tmpfit))
	
	tmpfit	<- lm(shr[,i]*y ~ as.matrix(tmpx))
	inits2	<- cbind(inits2, coef(tmpfit))
	print(i)
}
use.time	<- proc.time() - pct
cat("Initial logit finishes with", use.time[3]/60, "min.\n")

# Normalize the coefficients for demographics relative to those with grocery 
beta.init	<- inits[c(1:2),] - inits[c(1:2),beta0_base]
beta.init	<- c(t(beta.init[,-beta0_base]))
beta.init	<- c(beta.init, rowMeans(inits[-c(1:2),]))
gamma.init	<- c(t(inits2[c(1:2),]))
gamma.init	<- c(gamma.init, rowMeans(inits2[-c(1:2),]))
# theta.init	<- list(c(beta.init, gamma.init/100, -.5), 
# 					c(-1, -.8, -1, -.5, - 1, 	.3, .05,.4, .2, .4,		-.4, 0, .5, .1, 2,		.01, -.015, .1, -.2, -.6, 
# 					  2, 3, 2.5, 2.5, 3, 4.3, 	-.2, -.5, -.5, -.5, -.5, -.5, 	.1, -.003, .1, -.07, -2.5, 	-.02, .002, .07, .02, 1.2, -.5))
theta.init	<- list(c(beta.init/10, gamma.init/100, -.5), 
					c(-.1, -.08, -.1, -.05, -.1, 	.03, .005,.04, .02, .04,		-.04, 0, .05, .01, .2,		.001, -.0015, .01, -.02, -.06, 
					  .2, .3, .25, .25, .3, .4, 	-.02, -.05, -.05, -.05, -.05, -.05, 	.01, -.003, .01, -.007, -.25, 	-.002, .0002, .007, .002, .12, -.6))
for(i in 1:length(theta.init)){
	names(theta.init[[i]])	<- c(colnames(X1[[1]]), paste("gamma0_", colnames(X2[[1]]), sep=""), "ln_sigma")
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
# Register parallel computing
# mycore 	<- 2
# cl		<- makeCluster(mycore, type = "FORK")
# # cl		<- makeCluster(mycore, type = "PSOCK")			# Register cores in Windows
# registerDoParallel(cl)
# cat("Register", mycore, "core parallel computing. \n")

#-------------------------------------------- # 
# Compute the inclusive value with random draws
shr.par	<- coef(sol)
set.seed(687)
numsim 		<- 500
eps_draw	<- matrix(rgev(numsim*R, scale = exp(shr.par["ln_sigma"])), numsim, R)

# Subset simulation data
sel			<- duplicated(mydata[,c("household_code", "year")]) & mydata$dol <= quantile(mydata$dol, .95)
sht.X1		<- lapply(X1, function(x) x[!sel,])
sht.X2		<- lapply(X2, function(x) x[!sel,])
sht.price	<- price[!sel,]
sht.y		<- y[!sel]

pct			<- proc.time()
omega_draw	<- foreach(i = 1:numsim, .packages= c("nloptr"), .combine = 'rbind', .errorhandling='remove') %do% {
	incl_value_fn(param_est = shr.par, nx, base=beta0_base, X1 = sht.X1, X2= sht.X2, y = sht.y, price = sht.price, R, 
					eps_draw = rep(1,length(sht.y)) %*% t(eps_draw[i,]), returnmax = TRUE)$omega
}
use.time	<- proc.time() - pct
cat("Simulation of omega with random draws finishes with", use.time[3]/60, "min.\n")
cat("--------------------------------------------------------\n")
# stopCluster(cl)
# cat("Stop clustering. \n")

# Fit spline funciton
# cat("Summary of omega draws:\n"); print(summary(c(omega_draw))); cat("\n")
# tmp0	<- omega_draw
# colnames(tmp0)		<- 1:ncol(omega_draw)
# tmp0	<- melt(tmp0)				# The column of Var2 is the index of data
# tmpdat	<- data.frame(omega = tmp0$value, y = sht.y[tmp0$Var2], ln_inc = ln_inc[!sel][tmp0$Var2])
# tmpdat$price	<- sht.price[tmp0$Var2,]
# tmp		<- as.matrix(do.call(cbind, lapply(sht.X1, function(x) x[tmp0$Var2,selcol])))	
# tmpdat$X<- tmp			# Retail attributes

# # Check convexity of omega function 
# tmp1	<- seq(50, 1000, 50)
# tmp2	<- sort(unique(mydata$income_mid))
# tmpm	<- mesh(tmp1, tmp2)
# newx	<- data.frame(y = c(tmpm[[1]]), ln_inc = c(tmpm[[2]]))
# newx$price	<- rep(1, nrow(newx)) %*% t(colMeans(sht.price))
# newx$X<- rep(1, nrow(newx)) %*% t(colMeans(tmpdat$X))
# u 		<- matrix(omega_deriv(newx), dim(tmpm[[1]]))
# 
# pdf(paste(fname, "_omega.pdf", sep=""), width =6, height = 6)
# surf3D(tmpm[[1]], tmpm[[2]], u, bty = "b2", xlab = "y", ylab = "ln_inc", zlab = "omega", phi = 0, theta = 45)
# dev.off()

if(trim.alpha == 0){
	tmpdat	<- data.frame(omega = colMeans(omega_draw), y = sht.y, ln_inc = ln_inc[!sel])
}else{
	tmpdat	<- data.frame(omega = apply(omega_draw, 2, mean, trim = trim.alpha, na.rm=T), y = sht.y, ln_inc = ln_inc[!sel])
}
tmpdat$price	<- sht.price
tmpdat$X		<- as.matrix(do.call(cbind, lapply(sht.X1, function(x) x[,selcol])))	
omega_deriv	<- spl2dfun(tmpdat, fml = omega ~ s(y) + price + price*y)

# Check convexity of omega function 
newx	<- data.frame(y = seq(50, 1000, 20))
newx$price	<- rep(1, nrow(newx)) %*% t(colMeans(sht.price))
newx$X<- rep(1, nrow(newx)) %*% t(colMeans(tmpdat$X))
u 		<- omega_deriv(newx)

pdf(paste(gsub(".rdata", "", fname), "_omega.pdf", sep=""), width =6, height = 6)
plot(newx$y, u, type = "l")
dev.off()

# save.image(file = fname)
rm(list = c("sht.X1", "sht.X2", "sht.price", "sht.y", "tmpdat", "newx", "tmpm"))

###################################
# Upper level expenditue decision # 
###################################
rec		<- mydata$rec
tmpdat	<- data.frame(y = y, ln_inc = ln_inc)
tmpdat$price	<- price
tmp		<- do.call(cbind, lapply(X1, function(x) x[,selcol]))
tmpdat$X<- as.matrix(tmp)
o.deriv	<- omega_deriv(tmpdat, deriv = 1)
cat("sum(is.na(o.deriv)) =", sum(is.na(o.deriv)), ".\n")
rm(list = "tmpdat")

M_fn	<- function(lambda, omega_deriv, y, ln_inc, dT = NULL){
# Moment function: lambda * omega'(y,Inc) - 1 = 0
	m1	<- (lambda[1] + lambda[2] * rec)* o.deriv - 1
	m	<- cbind(m1, m1*ln_inc)
	mbar<- colMeans(m)
	return(list(moment = mbar, omega.derive = o.deriv, m = m))
}

GMM_fn	<- function(lambda, omega_deriv, y, ln_inc, dT, W = NULL){
# We use unit diagnial matrix as the weighting matric in GMM
	mm 	<- M_fn(lambda, omega_deriv, y, ln_inc, dT = dT)
	m	<- mm$moment
	if(is.null(W)){
		W	<- diag(length(m))
	}
	obj	<- - t(m) %*% W %*% m				# negative moment function for maxLik
	return(obj)
}

lsol	<- lm(1/o.deriv ~ factor(rec))
pct		<- proc.time()
init	<- matrix(c(coef(lsol), 
					.1, -.2, 
					.05, -.1), 
					3, 2, byrow = T)
colnames(init)	<- c("lambda1", "lambda2")
cat("Initial values are\n"); print(init); cat("\n")
dT		<- NULL
system.time(tmp <- GMM_fn(init[1,], omega_deriv, y, ln_inc))

tmp_sol	<- tmp_sol1 <- list(NULL)
for(i in 1:nrow(init)){
	tmp_sol[[i]]	<- maxLik(GMM_fn, start=init[i,], method="BFGS", omega_deriv = omega_deriv, y = y, ln_inc = ln_inc, dT = dT)
	tmp_sol1[[i]]	<- maxLik(GMM_fn, start=init[i,], method="NM", omega_deriv = omega_deriv, y = y, ln_inc = ln_inc, dT = dT)
}

# Choose the solution with max(-MM)
tmp		<- sapply(tmp_sol, function(x) ifelse(is.null(x$maximum), NA, x$maximum))
tmp1	<- sapply(tmp_sol1, function(x) ifelse(is.null(x$maximum), NA, x$maximum))
if(max(tmp1, na.rm =T) > max(tmp, na.rm=T)){
	tmp_sol	<- tmp_sol1
	cat("NM optimization has better results.\n")
}
sel 	<- which(abs(tmp - max(tmp, na.rm=T)) < 1e-4, arr.ind=T)
sel1	<- sapply(tmp_sol[sel], function(x) !any(is.na(summary(x)$estimate[,"Std. error"])) & !any(summary(x)$estimate[,"Std. error"]==Inf) )
sel		<- ifelse(sum(sel1)==1, sel[sel1], ifelse(sum(sel1)==0, sel[1], sel[sel1][1]))
sol.top	<- tmp_sol[[sel]]
if(!sol.top$code %in% c(0,1,2,8)){
	cat("Top level estimation does NOT have normal convergence.\n")
}
use.time <- proc.time() - pct
cat("--------------------------------------------------------\n")
summary(sol.top)
cat("Top level estimation finishes with", use.time[3]/60, "min.\n")

# 2nd step with updated W 
m			<- M_fn(lambda = coef(sol.top), omega_deriv = omega_deriv, y = y , ln_inc = ln_inc, dT = dT)
W			<- solve(var(m$m))
cat("Optimal weighting matrix is:\n"); print(W); cat("\n")
sol.top2	<- maxLik(GMM_fn, start=coef(sol.top), method="BFGS", omega_deriv = omega_deriv, y = y, ln_inc = ln_inc, dT = dT, W = W)
summary(sol.top2)
cat("Top level estimation finishes.\n")

####################
# Save the results #
####################
rm(list= intersect(ls(), c("hh_month", "Allocation_constr_fn","Allocation_fn","Allocation_nlop_fn","incl_value_fn","i","MDCEV_ll_fnC",
		  "MDCEV_LogLike_fnC","MDCEV_wrapper","tmp","tmp1","tmp2","tmp_sol","sel","sel1","sel2","param_assign","tmpX_list", 
		  "use.time", "pct", "uP_fn","uPGrad_fn", "theta_init", "make_plot", "ord","panelist","tmpidx","tmpn","cl", 
		  "tmp_sol1", "GMM_fn", "M_fn", "init", "m", "mycore", "param_assignR", "ggtmp", "ggtmp1", "tmpdat",
		  "interp.method", "trim.alpha", "mysplfun", "mytrimfun", "expFOC_fn", "exp_fn", "solveExp_fn","dT","i", "j", "k",  
		  "simExp_fn", "SimWrapper_fn", "SimOmega_fn", "cheb.1d.basis", "cheb.basis", "chebfun", "omega_parallel", 
		  "lastFuncGrad", "lastFuncParam", "args", "plot.wd", "tmpv1", "tmpv2")))

save.image(file = fname)

cat("This program is done. ")
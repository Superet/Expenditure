library(MCMCglmm)
library(parallelMCMCcombine)
library(r2excel)

# Save an aggregated data file
Nsec	<- 30
sim.ls	<- vector("list", length = Nsec)
simhet.ls	<- vector("list", length = Nsec)
for(i in 1:Nsec){
	load(paste("glmm_sub", i, ".rdata", sep=""))
	sim.ls[[i]]	<- sim[c("Sol","VCV")]
	simhet.ls[[i]]	<- sim_het[c("Sol","VCV")]
	rm("sim")
}
save(sim.ls, simhet.ls, file = paste("2_income_effect_glmm_", gsub("-", "", as.character(Sys.Date())), ".rdata", sep=""))

# Load parameters 
load("2_income_effect_glmm_20150813.rdata")
source("outreg function.R")
make_plot	<- TRUE
plot.wd		<- getwd()

# Make the table of coefficient estimates
outxls		<- paste(plot.wd, "/2_het_income_effect_incidence", gsub("-", "", as.character(Sys.Date())), ".csv", sep="")
# mywb		<- createWorkbook()
# sht1		<- createSheet(mywb, "Incidence")

############################
# Homogenous income effect # 
############################
Nsec	<- 30
numsim	<- 1000
nbeta	<- 18
nvcv	<- 42
fmt_name <- c("Convenience_Store", "Discount_Store", "Dollar_Store", "Drug_Store", "Grocery", "Warehouse_Club")
R		<- length(fmt_name)
beta	<- array(NA, c(nbeta, numsim, Nsec), 
				dimnames = list(c(paste(fmt_name,":Intercept", sep=""), paste(fmt_name, ":", "ln_income", sep=""), paste(fmt_name,":", "lag_dol", sep = "")), NULL, NULL))
vcv		<- array(NA, c(nvcv, numsim, Nsec),
				dimnames = list(c(paste(fmt_name,":hetvar", sep=""), paste(rep(fmt_name, R), ":", rep(fmt_name, each= R), sep="") ), NULL, NULL))
for(i in 1:Nsec){
	beta[,,i]	<- t(sim.ls[[i]]$Sol[,1:nbeta])
	vcv[,,i]	<- t(sim.ls[[i]]$VCV[,1:nvcv])
}

# Drop the fixed variance terms from vcv
sel		<- setdiff(dimnames(vcv)[[1]], paste(fmt_name, ":", fmt_name, sep=""))
vcv		<- vcv[sel,,]

# Make the trace plots
b.mcmc	<- do.call("mcmc.list", lapply(1:Nsec, function(i) as.mcmc(t(beta[,,i]))))
v.mcmc	<- do.call("mcmc.list", lapply(1:Nsec, function(i) as.mcmc(t(vcv[,,i]))))
if(make_plot){
	pdf(paste(plot.wd, "/graph_trace_multisec.pdf", sep=""), width = 8, height = 8)
	plot(b.mcmc)
	plot(v.mcmc)
	dev.off()
}

# Consensus method 
full.beta	<- consensusMCcov(subchain=beta, shuff=FALSE)
full.vcv	<- consensusMCindep(subchain=vcv, shuff=FALSE)

# Semiparametric method
norm.beta.est 	<- rep(0, nbeta)
norm.vcv.est	<- rep(0, nvcv - R)
for(i in 1:nbeta){
	for(j in 1:Nsec){
		norm.beta.est[i] 	<- norm.beta.est[i] + var(beta[i,,j])
	}
}
for(i in 1:(nvcv-R)){
	for(j in 1:Nsec){
		norm.vcv.est[i] 	<- norm.vcv.est[i] + var(vcv[i,,j])
	}
}

norm.beta.sd	<- sqrt(norm.beta.est/Nsec)
norm.vcv.sd		<- sqrt(norm.vcv.est/Nsec)
h_opt_b			<- (4/(nbeta+2))^(1/(4+nbeta)) * (numsim^(-1/(4+nbeta))) * norm.beta.sd
h_opt_v			<- (4/(nvcv+2))^(1/(4+nvcv)) * (numsim^(-1/(4+nvcv))) * norm.vcv.sd

full.beta1 		<- semiparamDPE(subchain = beta, bandw = h_opt_b * 2, anneal = FALSE)
full.vcv 		<- semiparamDPE(subchain = vcv, bandw = h_opt_v * 2, anneal = FALSE)

# Plot the two distribution
if(make_plot){
	ggtmp	<- rbind(data.frame(melt(full.beta), method = "census"), data.frame(melt(full.beta1), method = "semi"))
	ggplot(ggtmp, aes(value, color = method)) + geom_density(position = "identity") + 
			facet_wrap(~ Var1, scales = "free", ncol = 6)
}

# Coefficients 
full.beta1	<- as.mcmc(t(full.beta1))
tmp.tab		<- cbind(HPDinterval(full.beta1), postior.mode = posterior.mode(full.beta1), sd = apply(full.beta1, 2, sd))
rownames(tmp.tab)	<- dimnames(beta)[[1]]

# Covariance
full.vcv1	<- as.mcmc(t(full.vcv))
tmp.tab1		<- cbind(HPDinterval(full.vcv1), postior.mode = posterior.mode(full.vcv1), sd = apply(full.vcv1, 2, sd))
rownames(tmp.tab1)	<- dimnames(vcv)[[1]]

# Reshape the regression results 
tmp.tab2	<- data.frame(rbind(tmp.tab, tmp.tab1[1:R,]))
tmp			<- sapply(rownames(tmp.tab2), function(x) strsplit(x, ":")[[1]])
tmp.tab2$retailer <- tmp[1,]
tmp.tab2$Variable	<- tmp[2,]
tmp			<- split(tmp.tab2, tmp.tab2$retailer)
tmp.tab2	<- model_outreg(tmp, p.given = FALSE, head.name = c("postior.mode","sd","","Variable"), digits = 4)

# Covariance matrics
tmp.tab1[-c(1:R), "postior.mode"]
tmp			<- rep(1, R)
names(tmp)	<- paste(fmt_name,":",fmt_name, sep="")
tmp.tab3 	<- c(tmp.tab1[-c(1:R), "postior.mode"], tmp)
ord <- sort(names(tmp.tab3))
tmp.tab3	<- matrix(tmp.tab3[ord], R, R, byrow = T, dimnames = list(fmt_name, fmt_name))

# # Export to excel
# xlsx.addHeader(mywb, sht1, value = "Coefficients from Multivariate-logit regressions of shopping incidence", level = 2)
# xlsx.addTable(mywb, sht1, tmp.tab, row.names=FALSE)
# xlsx.addHeader(mywb, sht1, value = "Variance estimates from Multivariate-logit regressions of shopping incidence", level = 2)
# xlsx.addTable(mywb, sht1, tmp.tab1, row.names=FALSE)
# saveWorkbook(mywb, outxls)

# Export to csv
f	<- file(outxls, "w")
writeLines("Coefficients from Multivariate-logit regressions of shopping incidence", f)
write.csv(tmp.tab, f)
writeLines("\nVariance estimates from Multivariate-logit regressions of shopping incidence", f)
write.csv(tmp.tab1, f)

writeLines("\nRegresion coefficients:", f)
write.csv(tmp.tab2, f)
writeLines("\nCovariance matrix:", f)
write.csv(round(tmp.tab3, 4), f)
close(f)

#############################
# Heteogenous income effect # 
#############################
Nsec	<- 30
numsim	<- 1000
nbeta	<- 32
nvcv	<- 42
fmt_name <- c("Convenience_Store", "Discount_Store", "Dollar_Store", "Drug_Store", "Grocery", "Warehouse_Club")
R		<- length(fmt_name)
vcv		<- array(NA, c(nvcv, numsim, Nsec),
				dimnames = list(c(paste(fmt_name,":hetvar", sep=""), paste(rep(fmt_name, R), ":", rep(fmt_name, each= R), sep="") ), NULL, NULL))
sel		<- c( grep("month", colnames(simhet.ls[[1]]$Sol)), 
			  grep("household_code", colnames(simhet.ls[[1]]$Sol)) )
sel		<- setdiff(c(1:ncol(simhet.ls[[1]]$Sol)), sel)	
colnames(simhet.ls[[1]]$Sol)[sel]
beta	<- array(NA, c(nbeta, numsim, Nsec), 
				dimnames = list(c(paste("All:",c("Med", "Low","High"), sep=""), paste(fmt_name[-1],":Intercept", sep=""), 
								paste(fmt_name, ":", "ln_income", sep=""), 
								paste(fmt_name, ":", "lag_dol", sep=""), 
								paste(rep(fmt_name, 2),":", "ln_income-", rep(c("Low", "High"), each = R), sep = "")), NULL, NULL))
			
for(i in 1:Nsec){
	beta[,,i]	<- t(simhet.ls[[i]]$Sol[,sel])
	vcv[,,i]	<- t(simhet.ls[[i]]$VCV[,1:nvcv])
}

# Drop the fixed variance terms from vcv
sel		<- setdiff(dimnames(vcv)[[1]], paste(fmt_name, ":", fmt_name, sep=""))
vcv		<- vcv[sel,,]

# Make the trace plots
b.mcmc	<- do.call("mcmc.list", lapply(1:Nsec, function(i) as.mcmc(t(beta[,,i]))))
v.mcmc	<- do.call("mcmc.list", lapply(1:Nsec, function(i) as.mcmc(t(vcv[,,i]))))
if(make_plot){
	pdf(paste(plot.wd, "/graph_trace_multisec_het.pdf", sep=""), width = 8, height = 8)
	plot(b.mcmc)
	plot(v.mcmc)
	dev.off()
}

# Consensus method 
full.beta	<- consensusMCcov(subchain=beta, shuff=FALSE)
full.vcv	<- consensusMCindep(subchain=vcv, shuff=FALSE)

# Semiparametric method
norm.beta.est 	<- rep(0, nbeta)
norm.vcv.est	<- rep(0, nvcv - R)
for(i in 1:nbeta){
	for(j in 1:Nsec){
		norm.beta.est[i] 	<- norm.beta.est[i] + var(beta[i,,j])
	}
}
for(i in 1:(nvcv-R)){
	for(j in 1:Nsec){
		norm.vcv.est[i] 	<- norm.vcv.est[i] + var(vcv[i,,j])
	}
}

norm.beta.sd	<- sqrt(norm.beta.est/Nsec)
norm.vcv.sd		<- sqrt(norm.vcv.est/Nsec)
h_opt_b			<- (4/(nbeta+2))^(1/(4+nbeta)) * (numsim^(-1/(4+nbeta))) * norm.beta.sd
h_opt_v			<- (4/(nvcv+2))^(1/(4+nvcv)) * (numsim^(-1/(4+nvcv))) * norm.vcv.sd

full.beta1 		<- semiparamDPE(subchain = beta, bandw = h_opt_b * 2, anneal = FALSE)
full.vcv1 		<- semiparamDPE(subchain = vcv, bandw = h_opt_v * 2, anneal = FALSE)

# Plot the two distribution
if(make_plot){
	ggtmp	<- rbind(data.frame(melt(full.beta), method = "census"), data.frame(melt(full.beta1), method = "semi"))
	ggplot(ggtmp, aes(value, color = method)) + geom_density(position = "identity") + 
			facet_wrap(~ Var1, scales = "free", ncol = 6)	
}

# Make the table of coefficient estimates
# Coefficients 
full.beta1	<- as.mcmc(t(full.beta1))
tmp.tab		<- cbind(HPDinterval(full.beta1), postior.mode = posterior.mode(full.beta1), sd = apply(full.beta1, 2, sd))
rownames(tmp.tab)	<- dimnames(beta)[[1]]

# Covariance
full.vcv1	<- as.mcmc(t(full.vcv))
tmp.tab1		<- cbind(HPDinterval(full.vcv1), postior.mode = posterior.mode(full.vcv1), sd = apply(full.vcv1, 2, sd))
rownames(tmp.tab1)	<- dimnames(vcv)[[1]]

# Reshape the regression results 
tmp.tab2	<- data.frame(rbind(tmp.tab, tmp.tab1[1:R,]))
tmp			<- sapply(rownames(tmp.tab2), function(x) strsplit(x, ":")[[1]])
tmp.tab2$retailer <- tmp[1,]
tmp.tab2$Variable	<- tmp[2,]
tmp			<- split(tmp.tab2, tmp.tab2$retailer)
tmp.tab2	<- model_outreg(tmp, p.given = FALSE, head.name = c("postior.mode","sd","","Variable"), digits = 4)

# Covariance matrics
tmp.tab1[-c(1:R), "postior.mode"]
tmp			<- rep(1, R)
names(tmp)	<- paste(fmt_name,":",fmt_name, sep="")
tmp.tab3 	<- c(tmp.tab1[-c(1:R), "postior.mode"], tmp)
ord <- sort(names(tmp.tab3))
tmp.tab3	<- matrix(tmp.tab3[ord], R, R, byrow = T, dimnames = list(fmt_name, fmt_name))

# # Export to excel
# xlsx.addHeader(mywb, sht1, value = "Coefficients from Multivariate-logit regressions of shopping incidence", level = 2)
# xlsx.addTable(mywb, sht1, tmp.tab, row.names=FALSE)
# xlsx.addHeader(mywb, sht1, value = "Variance estimates from Multivariate-logit regressions of shopping incidence", level = 2)
# xlsx.addTable(mywb, sht1, tmp.tab1, row.names=FALSE)
# saveWorkbook(mywb, outxls)

# Export to csv
f	<- file(outxls, "at")
writeLines("Coefficients from Multivariate-logit regressions of shopping incidence with heterogeneous income effects", f)
write.csv(tmp.tab, f)
writeLines("\nVariance estimates from Multivariate-logit regressions of shopping incidence with heterogeneous income effects", f)
write.csv(tmp.tab1, f)

writeLines("\nRegresion coefficients:", f)
write.csv(tmp.tab2, f)
writeLines("\nCovariance matrix:", f)
write.csv(round(tmp.tab3, 4), f)
close(f)


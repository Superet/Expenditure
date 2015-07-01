# ###################
# # Marginal effect #
# ###################
# marginal_dat<- data.frame(NULL)
# 
# # Marginal effect of price 
# my_change <- .1
# selr	<- c(2,5,7); fmt_name[selr]
# for(k in 1:length(bskt_type)){
# 	tmp_coef	<- as.vector(coef(sol_list[[k]]))
# 	for(r in selr){
# 		pct 	<- proc.time()
# 		tmp_X0	<- all_X_list[[k]]
# 		tmpPrice0<- all_price[[k]]
# 		tmpPrice1<- tmpPrice0
# 		tmpPrice1[,r] <- tmpPrice1[,r] + my_change
# 		tmpy	<- rowSums(all_eR[[k]])
# 		tmpQ	<- bskt_typeQ[k]
# 		tmp		<- marginal_fn(tmp_coef, tmpy, tmpQ, tmp_X0, tmp_X0, tmpPrice0, tmpPrice1, random = FALSE,
# 								R=R, Ra=R, qz_cons = 0, exp_outside = FALSE, quant_outside = FALSE)
# 		marginal_dat <- rbind(marginal_dat, data.frame(bskt_type = rep(k,R), variable = "price", var_channel = fmt_name[r], 
# 														delta = my_change, channel_type = fmt_name, tmp))
# 		use.time<- proc.time() - pct
# 		cat("Maringal effect conditional on basket",k, ": changeing price for", fmt_name[r], 
# 			", time =", use.time[3]/60,"min.\n")
# 	}
# }
# 
# # Marginal effect of covariates
# selcol	<- c("size_index", "ln_num_module", "ln_upc_per_mod", "overall_prvt")
# my_change	<- c(1, log(1.1), log(1.1), .1)
# selr	<- c(2,5,7); fmt_name[selr]
# for(k in 1:length(bskt_type)){
# 	tmp_coef	<- as.vector(coef(sol_list[[k]]))
# 	for(j in 1:nx){		
# 		for(r in selr){
# 			pct 	<- proc.time()
# 			tmp_X0	<- all_X_list[[k]]
# 			tmp_X1	<- all_X_list[[k]]
# 			tmp_X1[[r]][,j]	<- tmp_X1[[r]][,j] + my_change[j]
# 			tmpy	<- rowSums(all_eR[[k]])
# 			tmpQ	<- bskt_typeQ[k]
# 			tmp		<- marginal_fn(tmp_coef, tmpy, tmpQ, tmp_X0, tmp_X1, all_price[[k]], all_price[[k]], random = FALSE,
# 									R=R, Ra=R, qz_cons = 0, exp_outside = FALSE, quant_outside = FALSE)
# 			marginal_dat <- rbind(marginal_dat, data.frame(bskt_type = rep(k,R), variable = selcol[j], var_channel = fmt_name[r], 
# 															delta = my_change[j], channel_type = fmt_name, tmp))
# 			use.time<- proc.time() - pct
# 			cat("Maringal effect conditional on basket",k, ": changeing variable", selcol[j], "for", fmt_name[r], 
# 				", time =", use.time[3]/60,"min.\n")
# 		}
# 	}
# }
# 
# # Plot the matrix of marginal effect. 
selcol	<- c("price", selcol)
pdf(paste(plot.wd, "/graph_marginal_effect_a1b1.pdf",sep=""), width = 10, height = 7)
for(i in 1:length(selcol)){
	ggtmp	<- subset(marginal_dat, variable == selcol[i])
	plots <- ggplot(ggtmp, aes(channel_type, dif_mean)) + geom_bar(stat="identity", position = "dodge") + 
				facet_grid( var_channel ~ bskt_type) + 
				theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
				labs(x = "Affected channels", y = "Change in expenditure share", 
					title = paste("Change of expenditure share from increasing \n",selcol[i], " by ",round(unique(ggtmp$delta),2),sep=""))
	print(plots)
}
dev.off()

# Compute elasticity
tmp		<- sapply(all_eR, function(x) colMeans(x/rowSums(x)))
tmp		<- melt(tmp)
tmp$channel_type <- gsub("_"," ", substr(tmp$Var1, 5, 25))
tmp		<- tmp[,-1]
names(tmp)	<- c("bskt_type","share","channel_type")
ggtmp	<- merge(marginal_dat, tmp, by=c("bskt_type","channel_type"), all.x=T)

# Append the mean of retailers' attribute mean
# NOTE: the merging variable is the var_channel
tmp		<- tapply(price_dat$price_paid_norm_index, price_dat$channel_type, mean)
sel		<- ggtmp$variable == "price"
ggtmp[sel,"var_mean"]	<- tmp[ggtmp[sel,"var_channel"]]
tmp		<- tapply(fmt_attr$size_index, fmt_attr$channel_type, mean)
sel		<- ggtmp$variable == "size_index"
ggtmp[sel,"var_mean"]	<- tmp[ggtmp[sel,"var_channel"]]
tmp		<- tapply(fmt_attr$overall_prvt, fmt_attr$channel_type, mean)
sel		<- ggtmp$variable == "overall_prvt"
ggtmp[sel,"var_mean"]	<- tmp[ggtmp[sel,"var_channel"]]
sel		<- is.na(ggtmp$var_mean)
ggtmp[sel,"var_mean"] 	<- 1
ggtmp$elasticity 	<- with(ggtmp, (dif_mean/share)/(delta/var_mean))

pdf(paste(plot.wd, "/graph_elasticity_a1b1.pdf",sep=""), width = 10, height = 7)
for(i in 1:length(selcol)){
	plots <- ggplot(subset(ggtmp, variable == selcol[i]), aes(channel_type, elasticity)) + geom_bar(stat="identity", position = "dodge") + 
				facet_grid( var_channel ~ bskt_type) + 
				theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
				labs(x = "Affected channels", y = "Change in expenditure share", 
					title = paste("Elasticity of ", selcol[i], sep=""))
	print(plots)
}
dev.off()

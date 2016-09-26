model_outreg <- function(model.list,p.given=FALSE,varname=NULL,varlabs=NULL,digits=2,modelname=NULL,latex.superscript=F,
						head.name=c("estimate","std. error")){
# Return result table wit p-value level indication star
# =======================================================
# model.list 		...	A list of models, each element has "Estimate" and "Std.Error"
# p.given			...	A logical indicator, TRUE if the model list includes p-values, 
#						 if FALSE then the function computes z-value and two-sided p-value
# varname			...	A vector of variable characters
# varlabs			...	A vector of variable 
# digits			...	rounding number
# modelname			...	A vector of model names
# latex.superscript	... Logical values, if T, it returns Latex style. 
	nm <- length(model.list)
	if(is.null(modelname)){
		if(is.null(names(model.list))){
			modelname <- paste("Model",1:nm,sep="")
		}else{
			modelname <- names(model.list)
		}
	}
	for(i in 1:nm){
		colnames(model.list[[i]]) <- tolower(colnames(model.list[[i]]))
	}
	if(p.given){
		out.list <- lapply(1:nm,function(i) 
					outreg(model.list[[i]][,head.name[1]],model.list[[i]][,head.name[2]],pval=model.list[[i]][,"pr(> t)"],
						   varname,varlabs,digits,modelname=modelname[i],latex.superscript) )
	}else{
		out.list <- lapply(1:nm,function(i) 
					outreg(model.list[[i]][,head.name[1]],model.list[[i]][,head.name[2]],pval=NULL,
						   varname,varlabs,digits,modelname=modelname[i],latex.superscript) )
	}
	
	out <- do.call(cbind,out.list)
	if(nm>1){
		sel <- names(out)=="Variable" & 1:ncol(out)>1
		out <- out[,!sel]
	}
	out
}


outreg <- function(estimate,se,pval=NULL,varname=NULL,varlabs=NULL,digits=2,modelname=NULL,latex.superscript=F){
# Return result table with p-value level indication star
# NOTE: the p-value is computed by assuming normal distribution, not t distribution	 
#==================================================================================
# estimate			...	A vector of model estimate
# se				... A vector of standard error
# varname			...	A vector of variable characters
# varlabs			...	A vector of variable 
# digits			...	rounding number
# modelname			...	A character of model name
# latex.superscript	... Logical values, if T, it returns Latex style.
	
	
	n <- length(estimate)
	if(is.null(varname)){
		if(is.null(names(estimate))){
			varname <- paste("X",1:n,sep="")
		}else{
			varname <- names(estimate)
		}
	}else{
		estimate <- estimate[varname]
		se <- se[varname]
		n <- length(estimate)
	}
	if(is.null(varlabs)){
		varlabs <- varname
	}
	if(is.null(pval)){
		zval <- estimate/se
		pval <- (1-pnorm(abs(zval)))*2
	}else{
		if(any(!is.null(varname))){
			pval <- pval[varname]
		}
	}
		
	sig.lev <- c(.01,.05,.1)
	sig.lab <- c("***","**","*","")
	names(sig.lab) <- c(1:3,0)
	sig <- apply(sapply(1:3,function(x) x*(pval<sig.lev[x])),1,min)
	
	out <- rep(NA,2*n)
	out[seq(1,2*n,2)] <- paste(formatC(estimate,digits=digits,format="f"),sig.lab[as.character(sig)],sep="")
	out[seq(2,2*n,2)] <- paste("(",formatC(se,digits=digits,format="f"),")",sep="")
	if(latex.superscript){
		sel <- which(sig.lab[as.character(sig)]!="")
		out[(sel*2-1)] <- paste("$",formatC(estimate[sel],digits=digits,format="f"),"^{",sig.lab[as.character(sig)][sel],"}$",sep="")
	}
	out <- data.frame(Variable=rep(NA,2*n),out)
	colnames(out) <- c("Variable",modelname)
	out[seq(1,2*n,2),"Variable"] <- as.character(varname)
	out$Variable <- factor(out$Variable,levels=varname,labels=varlabs)
	out
}
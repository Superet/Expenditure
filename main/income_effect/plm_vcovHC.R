# The code is taken from CRAN git website. 
# The original code of vcovHC in the package plm only allows for clustering over group or time, 
# where group and time are the index in the panel.data. 

# I modified the code such that it can allow for any clustering. 


vcovHC.plm.new <-function(x,method=c("arellano","white1","white2"),
                    type=c("HC0", "HC1", "HC2", "HC3", "HC4"),
                    cluster=c("group","time","define"), clustervec = NULL, ...) {
  ## Robust vcov for panel models (pooling, random, within or
  ## fd -type plm obj.)
  ##
  ## This version: October 20th 2009; allows choosing the clustering dimension
  ## so as to have serial- or x-sectional-correlation robustness;
  ## indices redone to cope with unbalanced data.
  ## 'pcse' SEs à la Beck and Katz (1995) moved to another function.
  ##
  ## This function takes the demeaned data from the
  ## plm object, then vcovHC as in Greene, Ec. An. (2003), pag. 315
  ## (Greene mentions the RE case on page 316, next par.)
  ## The Arellano estimator is referenced in Wooldridge as well:
  ## Wooldridge 2002, Econometrics of cross-section and
  ## panel data:
  ## 10.5.4 (and formula (10.59)) for FE/within,
  ## 10.4.2 (referring to formula 7.49) for RE.
  ##
  ## This version: compliant with plm.1.2-0; lmtest.
  ##
  ## Usage:
  ## myplm <- plm(<model>,<data>,...)
  ## # default (Arellano (1987)):
  ## coeftest(myplm, vcov=vcovHC)
  ## # White1:
  ## coeftest(myplm, vcov=function(x) pvcovHC(x,method="white1"))
  ## # idem, HC3 weighting:
  ## coeftest(myplm, vcov=function(x) pvcovHC(x,method="arellano",type="HC3"))
  ## waldtest(myplm,update(myplm,<new formula>),vcov=pvcovHC)
  ##
  ## This weighted version implements a system of weights as
  ## in vcovHC/meatHC. Sure this makes sense for white1, but it
  ## is open to question for white2 and arellano. We'll see.
  ##
  ## This version allows for choosing the grouping of covariance
    ## blocks in Arellano and White2. 'cluster="group"' means that vcov(e) is
    ## grouped according to individuals and thus the estimator is robust vs.
    ## heteroskedasticity *and serial correlation* (requires large N).
    ## Else it is robust vs. hetero *and xsectional correlation* (req. large T).
    ##
    ## Notice that this is independent of the type of effects, as everything is
    ## done on demeaned data.
  ##
  ## Results OK vs. Eviews, vcov=White cross-section/period, as long as
  ## the data are simple, else Eviews screws up... Check unbal.

    method <- match.arg(method)
    type <- match.arg(type)
    model <- x$args$model
    # if (!model %in% c("random", "within", "pooling", "fd")) {
    if (!model %in% c("random", "within", "fd", "pooling")) {
        stop("Model has to be either random, within or pooling model")
    }
	if(cluster == "define" & is.null(clustervec)){
		stop("User needs to provide clustervec when cluster is set 'define'.")
	}

  ## extract demeaned data
    demX <- model.matrix(x, model = model)
    demy <- pmodel.response(x, model = model)
    dimnames(demX)[[2]][1] <- attr(vcov(x), "dimnames")[[1]][1]

  ## extract dimensions not dependent on clustering
    pdim <- pdim(x)
    nT <- pdim$nT$N
    Ti <- pdim$Tint$Ti
    k <- dim(demX)[[2]]

  ## extract residuals
    uhat <- x$residuals

  ## robustifying against either serial or xs intragroup dependence:
  ## if 'group' then keep current indexing, if 'time' then swap i<->t
  ## so that residuals get 'clustered' by time period instead of by
  ## group (i.e. the vcov estimator is robust vs. xsectional dependence)

  ## extract indices
    groupind<-as.numeric(attr(x$model, "index")[,1])
    timeind<-as.numeric(attr(x$model, "index")[,2])
	if(!is.null(clustervec)){
		clusterind	<- as.numeric(factor(clustervec))
	}

    ## Achim's fix
     if(model == "fd") {
       groupind <- groupind[timeind > 1]
       timeind <- timeind[timeind > 1]
	   if(!is.null(clustervec)){clusterind	<- clusterind[timeind > 1]}
     }

    
  ## set grouping indexes
    switch(match.arg(cluster), group = {
           relevant.ind <- groupind
         }, time = {
           relevant.ind <- timeind
         }, define = {
			relevant.ind <- clusterind
		 })
    n <- length(unique(relevant.ind))
    tind <- vector("list", n)
    for (i in 1:n) {
        tind[[i]]<-which(relevant.ind==i)
    }

  ## extract residuals
  uhat<-x$residuals

  ## define residuals weighting function omega(res)
  ## (code taken from meatHC and modified)
  ##
  ## theor. comment:
  ## here we are decomposing the corrected (uhat_i)^2 in White and
  ## MacKinnon in uhat_i*uhat_i, in order to construct the submatrix
  ## Omega_i as:
  ## white1: crossprod(uhat_i)
  ## white2: crossprod(rep(mean(uhat_i),length(uhat_i))
  ## arellano: outer(uhat_i)
  ##
  ## Let cf be the correction factor (HC0: 1, HC1:(n/n-k), HC2: 1/(1-h_ii),
  ## HC3: 1/(1-h_ii)^2 etc.); it is easily seen that cf>0, cf->1 as nT->Inf.
  ## Looks like the diagonal values of the
  ## hat matrix X(X'X)^(-1)X' are always <1, thus (1-h_ii) in HC2 is
  ## always positive as well. The same for HC3 where it is squared.
  ## Note that for HC>2, cf=cf(i).
  ## We transform the residuals by rcf=sqrt(cf), so that the diag elements
  ## of Omegai for the White est.s are:
  ## [Omegai]_ii = (rcf_i*uhat_i)^2 = cf_i * uhat_i^2 (White 1)
  ## [Omegai]_ii = 1/n*sum((rcf_i*uhat_i)^2) = 1/n*sum(rcf_i * uhat_i^2) (White 2)
  ## while the generic element of Omegai for Arellano is:
  ## [Omegai]_ij = rcf_i*uhat_i * rcf_j*uhat_j
  ## and the diagonal is as in White 1
  ## [Omegai]_ii = cf_i * uhat_i^2
  ##
  ## in HC1-2 we are correcting by non-negative quantities and taking
  ## only squares of errors.

    ## diaghat function for matrices
    dhat <- function(x) {tx<-t(x)
                         diag(crossprod(tx,solve(crossprod(x),tx)))}

    ## this is computationally heavy, do only if needed
    switch(match.arg(type), HC0 = {diaghat<-NULL},
                            HC1 = {diaghat<-NULL},
                            HC2 = {diaghat<-try(dhat(demX), silent = TRUE)},
                            HC3 = {diaghat<-try(dhat(demX), silent = TRUE)},
                            HC4 = {diaghat<-try(dhat(demX), silent = TRUE)})
    df <- nT - k
    switch(type, HC0 = {
            omega <- function(residuals, diaghat, df) residuals
        }, HC1 = {
            omega <- function(residuals, diaghat, df) residuals *
                sqrt(length(residuals)/df)
        }, HC2 = {
            omega <- function(residuals, diaghat, df) residuals /
                sqrt(1 - diaghat)
        }, HC3 = {
            omega <- function(residuals, diaghat, df) residuals /
                (1 - diaghat)
        }, HC4 = {
            omega <- function(residuals, diaghat, df) residuals/sqrt(1 -
                diaghat)^pmin(4, length(residuals) * diaghat/as.integer(round(sum(diaghat),
                digits = 0)))
        })

  ## transform residuals by weights
  uhat<-omega(uhat,diaghat,df)


  ## define Omegai(e_i) function for Omega_i diag. blocks in E^2
  ## in Greene's formula (top of page 315)
  switch(match.arg(method),
               white1 = {Omegai<-function(x) diag(x^2)},
               white2 = {Omegai<-function(x) {
                         n<-length(x)
                         diag(sum(x^2)/n,n)
                         }},
               arellano = {Omegai<-function(x) outer(x,x)})

  salame<-array(dim=c(k,k,n))
  for(i in 1:n) {
      groupinds<-tind[[i]]
      xi<-demX[groupinds,,drop=FALSE]
      ui<-uhat[groupinds]
      salame[,,i]<-crossprod(xi,Omegai(ui))%*%xi
      }

  ## meat
  salame<-apply(salame,1:2,sum)

  ## bread
  pane<-solve(crossprod(demX))

  ## sandwich
  mycov <- pane %*% salame %*% pane
  return(mycov)
}


vcovHC.sur	<- function(x,method=c("arellano","white1","white2"), type=c("HC0", "HC1", "HC2", "HC3", "HC4"), clustervec, ...){
	method <- match.arg(method)
	type <- match.arg(type)

	## extract demeaned data
	demX <- model.matrix(x)
	demy <- model.response(x$model)
	dimnames(demX)[[2]][1] <- attr(x$coefCov, "dimnames")[[1]][1]

	## extract dimensions not dependent on clustering
	nT	<- nrow(demX)
	# nT	<- 
	# pdim <- pdim(x)
	# nT <- pdim$nT$N
	# Ti <- pdim$Tint$Ti
	k <- dim(demX)[[2]]

	## extract residuals
	uhat <- x$residuals
	
	## extract indices
	clusterind	<- as.numeric(factor(clustervec))

  	## set grouping indexes
	relevant.ind <- clusterind	
    n 	<- length(unique(relevant.ind))
    tind <- vector("list", n)
    for (i in 1:n) {
        tind[[i]]<-which(relevant.ind==i)
    }

  	## diaghat function for matrices
    dhat <- function(x) {tx<-t(x)
                         diag(crossprod(tx,solve(crossprod(x),tx)))}

    ## this is computationally heavy, do only if needed
    switch(match.arg(type), HC0 = {diaghat<-NULL},
                            HC1 = {diaghat<-NULL},
                            HC2 = {diaghat<-try(dhat(demX), silent = TRUE)},
                            HC3 = {diaghat<-try(dhat(demX), silent = TRUE)},
                            HC4 = {diaghat<-try(dhat(demX), silent = TRUE)})
    df <- nT - k
    switch(type, HC0 = {
            omega <- function(residuals, diaghat, df) residuals
        }, HC1 = {
            omega <- function(residuals, diaghat, df) residuals *
                sqrt(length(residuals)/df)
        }, HC2 = {
            omega <- function(residuals, diaghat, df) residuals /
                sqrt(1 - diaghat)
        }, HC3 = {
            omega <- function(residuals, diaghat, df) residuals /
                (1 - diaghat)
        }, HC4 = {
            omega <- function(residuals, diaghat, df) residuals/sqrt(1 -
                diaghat)^pmin(4, length(residuals) * diaghat/as.integer(round(sum(diaghat),
                digits = 0)))
        })

	## transform residuals by weights
	uhat<-omega(uhat,diaghat,df)


	## define Omegai(e_i) function for Omega_i diag. blocks in E^2
	## in Greene's formula (top of page 315)
	switch(match.arg(method),
	             white1 = {Omegai<-function(x) diag(x^2)},
	             white2 = {Omegai<-function(x) {
	                       n<-length(x)
	                       diag(sum(x^2)/n,n)
	                       }},
	             arellano = {Omegai<-function(x) outer(x,x)})

	salame<-array(dim=c(k,k,n))
	for(i in 1:n) {
	    groupinds<-tind[[i]]
	    xi<-demX[groupinds,,drop=FALSE]
	    ui<-uhat[groupinds]
	    salame[,,i]<-crossprod(xi,Omegai(ui))%*%xi
	    }

	## meat
	salame<-apply(salame,1:2,sum)

	## bread
	pane<-solve(crossprod(demX))

	## sandwich
	mycov <- pane %*% salame %*% pane
	return(mycov)
}
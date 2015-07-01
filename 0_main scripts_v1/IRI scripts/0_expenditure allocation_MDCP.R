library(ggplot2)
library(reshape2)
library(plyr)
library(xlsx)
library(data.table)
library(maxLik)

iri.wd <- "~/Documents/Research/Data/IRI/kelloggirir"
setwd("~/Documents/Research/Store switching/processed data")
source('~/Documents/Research/Store switching/Exercise/MDCP_function.R')

iri.wd <- "E:/Users/Projects/Project_Marketing_IRI/kelloggirir"
setwd("//tsclient/Resear1/Store switching/processed data")
source('//tsclient/Resear1/Store switching/Exercise/MDCP_function.R')

### Read demographics ###
demo <- read.delim(file=paste(iri.wd,"/DEMO.txt",sep=""), header=F,blank.lines.skip = T, sep="|", row.names=NULL,comment.char = '', quote = '')
sel <- colMeans(is.na(demo))!=1
demo <- demo[,sel]
colnames(demo) <- c("PANID", "FAMSIZE", "INCOME", "RACE", "CHILDREN", "FMLE_AGE", "FMLE_EDU", "FMLE_HRS", "FMLE_OCC", 
					"MALE_AGE", "MALE_EDU", "MALE_HRS", "MALE_OCC", "M_STATUS", "RENTOWN", "NUM_CATS", "NUM_DOGS",
					"REGION", "MKT", "PROJ09", "PROJ08", "PROJ07", "PROJ06", "ZIPCODE")
					
# Read trip and demo code
demo.code <- read.xlsx(paste(iri.wd,"/Academic Household File.xls",sep=""),sheetName="Demo Code")
demo.code <- demo.code[-c(1,2),]
colnames(demo.code) <- c("VARIABLE","CODE","DEMO_DSC")
sel <- apply(demo.code,1,function(x) all(is.na(x)|x==""))
demo.code <- demo.code[!sel,]

# Convert demographics codes
tmp <- unique(demo.code$VARIABLE)
for(i in 1:length(tmp)){
	sel <- demo.code$VARIABLE == tmp[i]
	tmp1 <- gsub("^\\s+|\\s+$","",as.character(demo.code[sel,"DEMO_DSC"]))
	demo[,as.character(tmp[i])] <- factor(demo[,as.character(tmp[i])],levels=demo.code[sel,"CODE"],labels=tmp1)
}

# Regroup some demographic variables
cumsum(table(demo$INCOME))/nrow(demo)
demo$INCOME1 <- cut(as.numeric(demo$INCOME),c(0,6,9,12),label=c("Low","Mid","High"))
demo$FAMSIZE1 <- ifelse(as.numeric(demo$FAMSIZE)>=3,3,as.numeric(demo$FAMSIZE))
demo$FAMSIZE1 <- factor(demo$FAMSIZE1,levels=1:3,labels=c("One","Two","Three and more"))

### Read the processed household trip data
trip <- read.csv("0_subset.csv",header=T)
fm.mkt <- read.csv("0_format_mkt_static.csv",header=T)
fm.price <- read.csv("0_format_price.csv",header=T)

sel <- demo.code$VARIABLE == "MKT"
tmp1 <- gsub("^\\s+|\\s+$","",as.character(demo.code[sel,"DEMO_DSC"]))
fm.mkt[,"MKT"] <- factor(fm.mkt[,"MKT"],levels=demo.code[sel,"CODE"],labels=tmp1)
fm.price[,"MKT"] <- factor(fm.price[,"MKT"],levels=demo.code[sel,"CODE"],labels=tmp1)
sel.mkt <- "Chicago"
fm.mkt <- subset(fm.mkt, MKT==sel.mkt)
fm.price <- subset(fm.price, MKT==sel.mkt)

fm.name <- c("Convenience","Natural_gourmet_food_stor","Drug","Dollar_variety","Limited_assortment_stores","Conventional_supermarkets",
			"Super_stores", "Mass", "Wholesale_clubs","Other")
dol.col <- paste("DOL_",fm.name,sep="")
qnt.col <- paste("QUANT_",fm.name,sep="")
prc.col <- paste("PRICE_",fm.name,sep="")
ast.col <- paste("AST_",fm.name,sep="")
demo.col <- c("INCOME1","FAMSIZE1")

#######################################
# Construct estimation data structure #
#######################################
LagSum <- function(x,t,lag.window){
	n <- length(t)
	R <- matrix(0,n,n)
	i <- 1
	while((t[i]+lag.window)<=t[n]){
		sel <- which(t>=t[i] & t<=t[i] + lag.window)
		R[i,sel] <- 1
		i <- i + 1
	}
	R[((n-lag.window+1):n),] <- NA
	
	# Currently drop the last few observation for each hh
	return(as.numeric(R %*% x))
}

# Classify basket 
cat.cut <- median(trip$num_cat)
qnt.cut <- quantile(trip$QUANT,.5)
trip$k1 <- with(trip,ifelse(num_cat<=cat.cut,"few cat","many cat"))
tmp <- 1*(trip[,qnt.col]>0)
by(tmp,trip$k1,colMeans)

trip$k2 <- with(trip,ifelse(QUANT<=qnt.cut,"small","large"))
trip$k <- 0
for(i in 1:length(unique(trip$k2))){
	for(j in 1:length(unique(trip$k1))){
		sel <- trip$k1 == unique(trip$k1)[j] & trip$k2 == unique(trip$k2)[i]
		trip[sel,"k"] <- (i-1)*length(unique(trip$k1))+j
	}
}

# Compute the budget and outside spending
trip.full <- merge(trip,demo[,c("PANID",demo.col)],by="PANID",all.x=T)
trip.full <- data.table(trip.full,key=c("PANID","WEEK"))
trip.full <- trip.full[,':='(y_free=LagSum(BSKTDOL,WEEK,lag.window=3),y_cstr=LagSum(BSKTDOL,WEEK,lag.window=1)),by=list(PANID)]
trip.full <- subset(trip.full,!is.na(y_free))

trip.full <- trip.full[,':='(y = ifelse(INCOME1=="Low",y_cstr,y_free))]
trip.full <- data.frame(trip.full)
trip.full$m0 <- trip.full$y - trip.full$DOL

N 	<- nrow(trip.full)
S 	<- length(fm.name)
eS	<- as.matrix(trip.full[,dol.col])
e0 	<- as.vector(trip.full[,"m0"])
e0[e0==0] <- 1
y 	<- as.vector(trip.full[,"y"])
n1S	<- eS>0
k.idx <- as.vector(trip.full[,"k"])
K 	<- 4

# price.mat <- as.matrix(trip.full[,prc.col])
tmp <- as.vector(fm.mkt[,"price"])
names(tmp) <- fm.mkt$FORMAT
tmp2 <- gsub("\\s","_",names(tmp))
tmp2 <- gsub("([\\/])","_",tmp2)
names(tmp) <- paste("PRICE_",tmp2,sep="")
price.mat <- rep(1,N) %*% t(tmp[prc.col])

### Format attributes ###
fmt.attr.list <- vector("list",3)
names(fmt.attr.list) <- c("ASSORT","ASSORT:TWO","ASSORT:Three") 
nx <- length(fmt.attr.list) 

# Format attribute: assortment
tmp <- data.frame(t(log(fm.mkt[,"allnum_upc"])))
names(tmp) <- fm.mkt$FORMAT
tmp2 <- gsub("\\s","_",colnames(tmp))
tmp2 <- gsub("([\\/])","_",tmp2)
names(tmp) <- paste("AST_",tmp2,sep="")
fmt.attr.list[["ASSORT"]] <- rep(1,nrow(trip.full)) %*% as.matrix(tmp[ast.col])

# Assortment * FamSize
fmt.attr.list[["ASSORT:TWO"]] <- as.matrix(1*(trip.full$FAMSIZE1=="TWO") %*% as.matrix(tmp[ast.col]))
fmt.attr.list[["ASSORT:Three"]] <- as.matrix(1*(trip.full$FAMSIZE1==3) %*% as.matrix(tmp[ast.col]))

### Choice set ###
CS.list <- list(fm.name,
				c("Convenience","Natural_gourmet_food_stor","Drug","Dollar_variety","Limited_assortment_stores"),
				c("Conventional_supermarkets","Super_stores", "Mass", "Wholesale_clubs","Other"),
				fm.name
				)
R.base <- diag(length(fm.name))
rownames(R.base) <- colnames(R.base) <- fm.name
R.list <- vector("list",length(unique(trip$k)))
for(i in 1:4){
	R.list[[i]] <- R.base[,CS.list[[i]]]
}

##############
# Estimation #
##############
my.iterlim <- 1000
my.printl <- 1

set.seed(999)
param.init <- c(rep(10,S),rep(-1,S),runif(nx))
names(param.init) <- c(paste("gamma:",fm.name,sep=""),paste("beta0:",fm.name),names(fmt.attr.list))
system.time(tmp <- MDCP_LogLik(param.init))

pct <- proc.time()
sol.list <- maxLik(MDCP_LogLik,start=param.init,method="BHHH",iterlim=my.iterlim,print.level=my.printl)
timeuse <- proc.time() - pct



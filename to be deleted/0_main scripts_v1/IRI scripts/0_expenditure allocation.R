library(ggplot2)
library(reshape2)
library(plyr)
library(xlsx)
library(data.table)
library(gridExtra)
library(scales)

iri.wd <- "~/Documents/Research/Data/IRI/kelloggirir"
setwd("~/Documents/Research/Store switching/processed data")
plot.wd <- "~/Desktop"
ww <- 6.5
ar <- .6

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

fm.cat <- read.csv("0_format_category.csv",header=T)
edible.dpt <- 5
fm.cat$dpt <- with(fm.cat,ifelse(DEPARTMENT<=edible.dpt,'Edible','Non-edible'))

## Subsetting a single market 
# sel <- demo.code$VARIABLE == "MKT"
# tmp1 <- gsub("^\\s+|\\s+$","",as.character(demo.code[sel,"DEMO_DSC"]))
# fm.mkt[,"MKT"] <- factor(fm.mkt[,"MKT"],levels=demo.code[sel,"CODE"],labels=tmp1)
# fm.price[,"MKT"] <- factor(fm.price[,"MKT"],levels=demo.code[sel,"CODE"],labels=tmp1)
# sel.mkt <- "Chicago"
# fm.mkt <- subset(fm.mkt, MKT==sel.mkt)
# fm.price <- subset(fm.price, MKT==sel.mkt)

fm.name1 <- c("Convenience","Drug","Dollar_variety","Limited_assortment_stores","Conventional_supermarkets",
			"Super_stores", "Mass", "Wholesale_clubs")
fm.name <- c("Convenience","Drug","Dollar/variety","Limited assortment stores","Conventional supermarkets",
			"Super stores", "Mass", "Wholesale clubs")
dol.col	<- paste("DOL_",fm.name1,sep="")
names(fm.name) <- dol.col

#################################			
# Retail format characteristics #
#################################
# Assortment 
tmp <- data.table(fm.cat)
tmp <- subset(tmp,FORMAT %in% fm.name)
tmp <- tmp[,list(num_cat=length(CATEGORY),num_brand=sum(num_brand),num_upc=sum(num_upc)),by=list(MKT,FORMAT,dpt)]
tmp <- tmp[,list(num_brand=mean(num_brand),num_upc=mean(num_upc)),by=list(FORMAT,dpt)]
ggtmp <- melt(tmp,id=c("FORMAT","dpt"))

plots <- list(NULL)
tmp1 <- unique(ggtmp$variable)
tmp.lab <- c("Number of brands","Number of UPC")
for(i in 1:length(unique(ggtmp$variable))){
	sel <- ggtmp$variable == tmp1[i]
	tmp2 <- tapply(ggtmp[sel,"value"],as.character(ggtmp[sel,"FORMAT"]),sum)
	ggtmp$FORMAT <- factor(ggtmp$FORMAT,levels=names(tmp2)[order(tmp2)])
	plots[[i]] <- ggplot(subset(ggtmp,variable==tmp1[i]),aes(FORMAT,value,fill=dpt)) + geom_bar(stat="identity") + 
					coord_flip() + ylab(tmp.lab[i]) +
					scale_fill_grey(name="Category group",start=.7,end=.2)
}

pdf(paste(plot.wd,"/graph_assortment.pdf",sep=""),width=ww,height=ww*.7)
do.call(grid.arrange,plots)
dev.off()

# Quantity
tmp <- data.table(fm.cat)
tmp <- subset(tmp,FORMAT %in% fm.name)
tmp <- tmp[,list(size_index=median(cat_size_index),unit=median(cat_unit)),by=list(FORMAT,dpt,CATEGORY)]
ggplot(tmp,aes(size_index,col=FORMAT)) + geom_density() + 
		facet_wrap(~dpt) + xlim(c(0,500))
ggplot(tmp,aes(unit,col=FORMAT)) + geom_density() + 
		facet_wrap(~dpt) + xlim(c(0,5))
ggplot(tmp,aes(FORMAT,size_index)) + stat_summary(fun.data=median_hilow, conf.int = .1, geom = "pointrange") + 
		facet_wrap(~dpt) + coord_flip()
		
cut.pt <- c(0,50,100,150,Inf)
tmp.lab <- c("<50",'50-99','100-149','150+')
tmp$size_index_grp <- cut(tmp$size_index,cut.pt,labels=tmp.lab)

ggtmp <- tmp[,list(count=length(CATEGORY)),by=list(FORMAT,dpt,size_index_grp)]
ggtmp <- ggtmp[,':='(percent=count/sum(count)),by=list(FORMAT,dpt)]
sel <- ggtmp$dpt == "Edible" & ggtmp$size_index_grp == tmp.lab[1]
ord <- order(ggtmp[sel,percent])
ggtmp$FORMAT <- factor(ggtmp$FORMAT,levels=as.character(ggtmp[sel,FORMAT])[ord])
my.color <- c("#cccccc","#969696","#636363","#252525")

pdf(paste(plot.wd,"/graph_packagesize.pdf",sep=""),width=ww*1.2,height=ww*ar)
ggplot(ggtmp,aes(FORMAT,percent,fill=size_index_grp)) + geom_bar() + 
		facet_wrap(~dpt) + coord_flip() + 
		scale_y_continuous(labels = percent) + ylab("Percent of categories") + 
		scale_fill_grey(name="Quantity index",start=.8,end=.2)
dev.off()
		
################################
# Bakset type and store choice #
################################
trip$edible_pct <- with(trip,quant_ed/(quant_ed+quant_noned))

par(mfrow=c(2,1))
hist(trip$QUANT)
hist(trip$edible_pct)
quantile(trip$QUANT,c(.33,.5,.67),na.rm=T)
quantile(trip$edible_pct,c(.33,.5,.67),na.rm=T)

variety.cut <- c(-Inf,.9,1)
quant.cut <- c(-Inf,31,Inf)
lab1 <- c("Mixed","Mainly edible")
lab2 <- c("Small","Large")
trip$variety <- cut(trip$edible_pct,variety.cut,labels=lab1)
trip$quant <- cut(trip$QUANT,quant.cut,labels=lab2)
trip$bskt <- with(trip,paste(quant,variety,sep="-"))
table(trip$bskt)

### Plot the trip incidence and dollar share contingent on basket type ###
# Compute average dollar share
tmp <- subset(trip,!is.na(DOL))
tmp1 <- tmp[,dol.col]/tmp$BSKT
ggtmp <- cbind(variety=tmp$variety,quant=tmp$quant,tmp1)
ggtmp <- data.table(melt(ggtmp,id=c("variety","quant")))
setnames(ggtmp,c("Variety","Quantity","variable","dol.share"))
ggtmp$Format <- fm.name[ggtmp$variable]
ggtmp$trip <- 1*(ggtmp$dol.share>0)
ggtmp <- ggtmp[,list(dol.share=mean(dol.share),dol.share.n0=mean(dol.share[dol.share>0])),by=list(Variety,Quantity,Format)]

# Compute trip incidence frequencey
tmp2 <- data.frame(variety=tmp$variety,quant=tmp$quant,1*(tmp1>0))
ggtmp1 <- data.table(melt(tmp2,id=c("variety","quant")))
setnames(ggtmp1,c("Variety","Quantity","variable","Incidence"))
ggtmp1 <- subset(ggtmp1,Incidence > 0)
ggtmp1$Format <- fm.name[ggtmp1$variable]
ggtmp1 <- ggtmp1[,list(count=length(Incidence)),by=list(Variety,Quantity,Format)]
ggtmp1 <- ggtmp1[,':='(Percent=count/sum(count)),by=list(Variety,Quantity)]
ggtmp <- merge(ggtmp,ggtmp1,by=c("Variety","Quantity","Format"))

# Plot 
sel <- ggtmp$Variety == lab1[1] & ggtmp$Quantity == lab2[1]
ord <- order(ggtmp[sel,dol.share])
ggtmp$Format <- factor(ggtmp$Format, levels=ggtmp[sel,Format][ord])
ggtmp1 <- melt(ggtmp,id=c("Variety","Quantity","Format"))

sel.var <- c("dol.share","Percent","dol.share.n0")
tmp.lab <- c("Average dollar share","Frequency of trip incidence","Average dollar share conditional on trip visit")
plots <- list(NULL)
for(i in 1:length(sel.var)){
	sel <- ggtmp1$variable == sel.var[i]
	plots[[i]] <- ggplot(ggtmp1[sel,],aes(Format,value)) + geom_bar(stat="identity") + 
					facet_grid(Variety ~ Quantity,labeller = "label_both") + coord_flip() +
					# ylim(c(0,.45)) + 
					scale_y_continuous(labels = percent) + ylab(tmp.lab[i]) 
}

pdf(paste(plot.wd,"/graph_contingent_share.pdf",sep=""),width=6.5, height=9.5)
do.call(grid.arrange,plots)
dev.off()









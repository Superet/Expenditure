library(ggplot2)
library(data.table)
library(RDSTK)
library(stringr)

setwd("~/Documents/Research/Store switching/Geographic data/store_location")
license	<- read.csv("Business_Licenses.csv", header = T)
head(license)
license$DOING.BUSINESS.AS.NAME	<- as.character(license$DOING.BUSINESS.AS.NAME)
license$DOING.BUSINESS.AS.NAME	<- toupper(gsub("\\s+", " ", str_trim(license$DOING.BUSINESS.AS.NAME)))		# Compress spaces
license$ADDRESS                 <- as.character(license$ADDRESS)
license$DATE.ISSUED				      	<- as.Date(as.character(license$DATE.ISSUED), format = "%m/%d/%Y")
license$LICENSE.TERM.START.DATE			<- as.Date(as.character(license$LICENSE.TERM.START.DATE), format = "%m/%d/%Y")
license$LICENSE.TERM.EXPIRATION.DATE	<- as.Date(as.character(license$LICENSE.TERM.EXPIRATION.DATE), format = "%m/%d/%Y")
range(license$DATE.ISSUED)

# Key words to identify retailers
mykey	<- list(CVS 			= c("CVS"),
 				Walgreens		= c("WALGREEN"),										# "WALGREENS", "WALGREEN'S"
				Dollar.General 	= "DOLLAR GENERAL",
				Family.Dollar	= "FAMILY DOLLAR",
				Dollar.Tree		= "DOLLAR TREE",
				WalMart			= c("WALMART #", "WAL-MART"),			# There is neighborhood market and Walmart express
				Target			= c("TARGET #", "TARGET STORE", "TARGET (STORE", "TARGET - STORE", "TARGET GROUP"),
				Jewel			= c("JEWEL FOOD STORE", "JEWEL-OSCO"),
				Osco			= "OSCO DRUG",
				Ultra.Foods		= "ULTRA FOODS", 
				Dominicks		= c("DOMINICK'S", "DOMINICKS"),
				Treasure		= c("TREASURE ISLAND FOODS"), 
				Marianos		= c("MARIANO'S"),
				Costco			= c("COSTCO WHOLESALE", "COSTCO BUSINESS CENTER"),
				Sams			= "SAM'S CLUB"
				)
exclude.keys	<- c("MARIANO'S GARAGE DOORS & CONSTRUCTION INCORPORATED", "MARIANO'S BURGER & GRILL", "CVS INDIANA, L.L.C.", 
					"CVS FARMS.INC", "FAMILY DOLLAR PLUS STORE", "FAMILY DOLLAR UP", "FAMILY DOLLAR PLUS #1", 
					"WALGREEN HOME MEDICAL CENTER")
					#, "TARGET & RESPONSE INCORPORATED", "TARGETCOM,INC", "TARGETCOM LLC.", "TARGET ZERO HVACR INC.", "TARGET MOBILE", "TARGET AIRE INC", "TARGET DATA")

# Extract inforamtion for each keyword
mysub	<- data.frame(NULL)
nm		<- data.frame(NULL)
for(i in 1:length(mykey)){
	for(j in 1:length(mykey[[i]])){
		sel	<- substr(license$DOING.BUSINESS.AS.NAME, 1, nchar(mykey[[i]][j])) == mykey[[i]][j]
		tmp	<- sort(table(license[sel,"DOING.BUSINESS.AS.NAME"]), decreasing = T)
		nm	<- rbind(nm, data.frame(retailer = names(mykey)[i], business = names(tmp), freq = tmp))
		mysub	<- rbind(mysub, license[sel,])
	}
}
nm			<- subset(nm, !business %in% exclude.keys)
nm$business	<- as.character(nm$business)
mysub		<- subset(mysub, !DOING.BUSINESS.AS.NAME %in% exclude.keys)
mysub 		<- merge(mysub, nm[,c("retailer", "business")], by.x = "DOING.BUSINESS.AS.NAME", by.y = "business", all.x = T)
mysub 		<- mysub[order(mysub$DOING.BUSINESS.AS.NAME),]

# Check the business names
tmp 	<- regexpr('#', nm$business)
tmp1	<- ifelse(tmp < 0, nm$business, substr(nm$business, 1, tmp-2))
sort(unique(tmp1))
dim(nm); length(unique(nm$business)); length(unique(mysub$DOING.BUSINESS.AS.NAME))
mysub 	<- mysub[order(mysub$DOING.BUSINESS.AS.NAME, mysub$LICENSE.TERM.START.DATE),]

# Construct store profiles for each unique stores
# ‘ISSUE’ is the record associated with the initial license application. 
# ‘RENEW’ is a subsequent renewal record. All renewal records are created with a term start date and term expiration date. 
# ‘C_LOC’ is a change of location record. It means the business moved. 
# ‘C_CAPA’ is a change of capacity record. Only a few license types may file this type of application. 
# ‘C_EXPA’ only applies to businesses that have liquor licenses. It means the business location expanded. 
# 'C_SBA' is a change of business activity record. It means that a new business activity was added or an existing business activity was marked as expired.
table(mysub$APPLICATION.TYPE)
tmp		<- data.table(mysub)
tmp		<- tmp[,list(naddress = length(unique(ADDRESS)), nzip = length(unique(ZIP.CODE)), 
					issue = sum(APPLICATION.TYPE == "ISSUE"), renew = sum(APPLICATION.TYPE == "RENEW"), 
					app.other = sum(! APPLICATION.TYPE %in% c("ISSUE", "RENEW")), 
					date.min = min(LICENSE.TERM.START.DATE), date.max = max(LICENSE.TERM.EXPIRATION.DATE)), 
				by = list(DOING.BUSINESS.AS.NAME)]
summary(tmp)
tmp[tmp$naddress > 1 |tmp$nzip > 1,]
tmp1 <- tmp[naddress > 1 & nzip>1,]
mysub$storeid <- mysub$DOING.BUSINESS.AS.NAME

# For stores with multiple address, asign different store numbers to them
for(i in 1:nrow(tmp1)){
  sel 	<- mysub$DOING.BUSINESS.AS.NAME == tmp1[i,DOING.BUSINESS.AS.NAME]
  tmp 	<- mysub[sel,"ADDRESS"]
  tmp2 	<- setNames(1:length(unique(tmp)), unique(tmp))
  mysub[sel,"storeid"] <- paste(mysub[sel,"DOING.BUSINESS.AS.NAME"], " No.", tmp2[as.character(tmp)], sep="")
}

# Check again if the store id is unique and if each store has unique location
tmp		<- data.table(mysub)
tmp		<- tmp[,list(naddress = length(unique(ADDRESS)), nzip = length(unique(ZIP.CODE)), 
                  issue = sum(APPLICATION.TYPE == "ISSUE"), renew = sum(APPLICATION.TYPE == "RENEW"), 
                  app.other = sum(! APPLICATION.TYPE %in% c("ISSUE", "RENEW")), 
                  nlat = length(unique(LATITUDE)), nlong = length(unique(LONGITUDE)),
                  lat.na = 1*all(is.na(LATITUDE)), 
                  date.min = min(LICENSE.TERM.START.DATE), date.max = max(LICENSE.TERM.EXPIRATION.DATE)), 
            by = list(storeid)]
summary(tmp)
tmp[naddress > 1,]			  # 4 observations that just have minor typo			
summary(tmp[naddress==1,])
sum(tmp$lat.na == 1)          # 14 stores have missing latitudes and longitudes.
mysub <- mysub[order(mysub$storeid, mysub$LICENSE.TERM.EXPIRATION.DATE),]

# Construct store profile
dim(mysub); length(unique(mysub$storeid))
store 	<- data.table(mysub)
store 	<- store[,':='(date.min = min(LICENSE.TERM.START.DATE), date.max = max(LICENSE.TERM.EXPIRATION.DATE), 
                     loc.change = sum(any(APPLICATION.TYPE == "C_LOC"))), 
                by = list(storeid)]
store 	<- subset(store, LICENSE.TERM.EXPIRATION.DATE == date.max)				# Keep the observation with the latest termination date
dim(store)
length(unique(mysub$storeid))
store 	<- unique(store, by = c("retailer", "storeid", "ADDRESS", "ZIP.CODE", "date.min", "date.max"))
store 	<- store[,list(retailer, storeid, LEGAL.NAME, DOING.BUSINESS.AS.NAME, ADDRESS, ZIP.CODE, CITY, STATE, LATITUDE, LONGITUDE,
                     loc.change,date.min, date.max)]
sum(is.na(store$LATITUDE))
table(as.character(store$LEGAL.NAME))

# Find the latitude and longitudes for the missing ones from address
sel   <- is.na(store$LATITUDE)
tmp   <- paste(store[sel,ADDRESS], store[sel,CITY], store[sel,STATE], store[sel,ZIP.CODE], sep=",")
tmp   <- gsub("  ", " ", tmp)
tmp1  <- do.call(rbind, lapply(tmp, function(x) street2coordinates(x)))
identical(tmp, tmp1$full.address)
store <- data.frame(store)
store[sel,"LATITUDE"]   <- tmp1$latitude
store[sel,"LONGITUDE"]  <- tmp1$longitude

# Format the data 
names(store)  <- tolower(names(store))
store$retailer <- as.character(store$retailer)
store$legal.name <- as.character(store$legal.name)
store$city    <- as.character(store$city)
store$state <- as.character(store$state)

# Change some retailer names
table(store$retailer)
(tmp   <- setNames(c("CVS", "Walgreen's", "Dollar General", "Family Dollar", "Dollar Tree", "WalMart", "Target", 
                     "Jewel", "Jewel", "Ultra Foods", "Dominick's", "Treasure Island", "Mariano's", "Costco","Sam's"), 
                    names(mykey)) )
store$retailer  <- tmp[as.character(store$retailer)]
table(store$retailer)  

# Number of locations for each retailer in 2016
tmp		<- data.table(subset(store, date.max >= as.Date("2016-01-01", format = "%Y-%m-%d")))
tmp		<- tmp[,list(n_loc = length(unique(storeid)), n_add = length(unique(address))), by = list(retailer)]
tmp

# There are some Osco drug stores that are the same with Jewel stores and Target
sel		<- which(store$retailer == "Jewel")
sel		<- sel[duplicated(store[sel,"address"])]
length(sel)
store	<- store[setdiff(1:nrow(store), sel),]

############################################################
# Combine the license data with location data from AggData #
############################################################
(fname	<- setdiff(list.files(), "Business_Licenses.csv"))
aggloc	<- data.frame(NULL)
for(i in 1:length(fname)){
	tmp	<- read.csv(fname[i], header = T, stringsAsFactors = F)
	names(tmp) 	<- tolower(names(tmp))
	names(tmp)	<- gsub("_", ".", names(tmp))
	if(fname[i] %in% c("target.csv", "walmart_08-01-16.csv")){
		tmp		<- subset(tmp, !store.type %in% c("Walmart Neighborhood Market", "Walmart Express", "City", "TargetExpress"))
	}
	aggloc	<- rbind(aggloc, 
					cbind(retailer = gsub(".csv","", fname[i]), tmp[,c("address", "city", "state", "zip.code", "latitude", "longitude")]))
}
aggloc	<- subset(aggloc, state == "IL" & toupper(city) %in% unique(as.character(license[license$STATE == "IL","CITY"])))
(tmp		<- setNames(c("7-eleven", "Costco", "CVS", "Dollar General", "Dollar Tree", "Family Dollar", "Kmart", "Sam's", "Target", "Walgreen's", "WalMart"), 
					unique(aggloc$retailer)))
aggloc$retailer	<- tmp[as.character(aggloc$retailer)]
table(aggloc$retailer)

# Compare the number of locations
# The numbers are similar, but the main difference is that AggData has more stores in suburbs 
tmp1	<- data.table(aggloc)
tmp1	<- tmp1[,list(agg = length(unique(address))), by = list(retailer)]
tmp2	<- data.table(store)
tmp2	<- tmp2[,list(license = length(unique(address))), by = list(retailer)]
tmp		<- merge(tmp1, tmp2, by = "retailer", all = T)
tmp
table(store$city); table(aggloc$city)

# Append the stores in surburbs to the main data 
tmp		<- subset(aggloc, city != "Chicago" & !retailer %in% c("7-eleven", "Kmart"))
tmp$city<- toupper(tmp$city)
store 	<- merge(store, tmp, by = names(tmp), all = T)

# Save data
# save(store, file = "~/Documents/Research/Store switching/Processed_data/store_location_license.rdata")
write.csv(store, file = "~/Documents/Research/Store switching/Processed_data/store_location_license.csv")

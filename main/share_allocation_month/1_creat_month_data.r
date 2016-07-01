# Create monthly expenditure data from biweekly data
library(data.table)

setwd("~/Documents/Research/Store switching/processed data")
load("hh_biweek_exp.rdata")

hh_exp	<- data.table(hh_exp)
hh_exp$month <- month(as.Date("2004-1-1", format="%Y-%m-%d") + 14*(hh_exp$biweek-1))
hh_exp	<- hh_exp[, list(dol = sum(dol), num_trip = sum(num_trip), num_day = sum(num_day), 
						 DOL_Convenience_Store = sum(DOL_Convenience_Store), DOL_Discount_Store = sum(DOL_Discount_Store), 
						 DOL_Dollar_Store = sum(DOL_Dollar_Store), DOL_Drug_Store = sum(DOL_Drug_Store), 
						 DOL_Grocery = sum(DOL_Grocery), DOL_Warehouse_Club = sum(DOL_Warehouse_Club), 
						 income_midvalue = unique(income_midvalue), first_income = unique(first_income), 
						 income_real = unique(income_real), first_incomeg = unique(first_incomeg), 
						 scantrack_market_descr = unique(scantrack_market_descr), household_size = unique(household_size), 
						 condo = unique(condo), employed = unique(employed), NumChild = unique(NumChild), 
						 stone_price = mean(stone_price), cpi = mean(cpi)), 
					by = list(household_code, year, month)]

# check if the observations are unique
dim(hh_exp)
dim(unique(hh_exp, by = c("household_code", "year", "month")))
hh_exp	<- data.frame(hh_exp)

ls()
save(hh_exp, price_dat, fmt_attr, file = "hh_month_exp.rdata")

setwd("~/Documents/Research/Store switching/processed data")

# Read data 
hh_exp		<- read.csv("2_hh_biweek_exp_merge.csv")
fmt_attr 	<- read.csv("1_format_year_attr.csv")
price_dat	<- read.csv("1_format_biweek_price.csv")

# Shuffle the order of household_code
set.seed(666)
panid		<- unique(hh_exp$household_code)
cat("Number of households =", length(panid), "\n")
ord			<- sample(1:length(panid))
names(ord)	<- panid
hh_exp$household_sample	<- ord[as.character(hh_exp$household_code)]

save(hh_exp, fmt_attr, price_dat, file = "hh_biweek_exp.rdata")
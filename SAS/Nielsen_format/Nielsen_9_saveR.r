library(data.table)

setwd("U:/Users/ccv103/Desktop")
setwd("U:/Users/ccv103/Documents/Research/Store switching/run")

hh_dist		<- read.csv("hh_format.csv", header = T)
hh_dpt		<- read.csv("hh_format_dpt.csv", header = T)
hh_month	<- read.csv("hh_month_dpt.csv", header = T)
pan			<- read.csv("panelists.csv", header = T)
fmt_dpt		<- read.csv("format_department.csv", header = T)
codebook	<- read.csv("U:/Users/ccv103/Documents/Research/Store switching/run/code_book.csv", header = T)
cpi			<- read.csv("CPI.csv")

names(hh_dist)
names(hh_dpt)
names(hh_month)
names(pan)
names(fmt_dpt)

# Rename the expenditure data
fmt_name	<- levels(hh_dist$channel_type)
dpt_name	<- c("DG", "GM", "HBC", "NFG", "OTHER", "RF")

tmp			<- names(hh_month)
ord			<- matrix(NA, length(fmt_name), length(dpt_name))
tmpn		<- setNames(fmt_name, substr(gsub("\\s", "_", fmt_name), 1, 8))
for(i in 1:length(dpt_name)){
	tmp1	<- paste("EXP_", gsub("\\s", "_", dpt_name[i]), sep="")
	sel		<- grep(tmp1, tmp)
	ord[,i]	<- sel
	tmp2	<- gsub(tmp1, "", tmp[sel])
	tmp2	<- tmpn[substr(tmp2, 1, 8)]
	tmp[sel]<- paste(dpt_name[i], ".", tmp2, sep="")
}
tmp		<- gsub("REPORT", "DOL", tmp)
names(hh_month)	<- tmp
hh_month	<- hh_month[,c(setdiff(1:ncol(hh_month), ord), ord)]

# Adjust expenditure with CPI
tmp			<- setNames(cpi$Annual, cpi$Year)
tmp			<- tmp/tmp[1]
tmp1		<- tmp[as.character(hh_month$year)]
for(i in 4:ncol(hh_month)){
	hh_month[,i]	<- hh_month[,i]/tmp1
}

###########################################
# Convert household demographic variables #
###########################################
# Income: use middle value of the bin
tmp				<- setNames(codebook$mid_range, codebook$code_value)
pan$income_mid	<- tmp[as.character(pan$income)]

# Number of kids
pan$num_kids	<- with(pan, ifelse(age_and_presence_of_children==9, 0, ifelse(age_and_presence_of_children<=3, 1,
								ifelse(age_and_presence_of_children<=6, 2, 3))))

# Use average age as the age for the household								
tmp				<- setNames(c(NA, 25, 27, 32, 37, 42, 47, 52, 59, 65), 0:9)
pan$age			<- rowMeans(cbind(tmp[as.character(pan$male_head_age)], tmp[as.character(pan$female_head_age)]), na.rm=T)

# A household is unemployed if nobody is employed for pay
# pan$unemploy	<- ifelse(pan$female_head_employment!=9 | pan$male_head_employment!=9, 0, 1)
summary(pan[,c("income_mid", "num_kids", "household_size", "age")])

# Segment household based on their first-year income and family size
tmp	<- data.table(pan)
setkeyv(tmp, c("household_code", "panel_year"))
tmp	<- tmp[,list(income = income_mid[1], household_size = household_size[1]), by = list(household_code)]
tmp1	<- cut(tmp$income, quantile(tmp$income, c(0:3)/3), include.lowest = TRUE, labels = c("Low", "Medium", "High"))
tmp2	<- cut(tmp$household_size, c(0, 1, 2, 9), include.lowest = TRUE, labels = c("Single", "Two", "Three+"))
table(tmp1, tmp2)
tmp3	<- paste(tmp1, tmp2, sep="*")
tmp$segment	<- factor(tmp3, levels = paste(rep(c("Low", "Medium", "High"), each = 3), c("Single", "Two", "Three+"), sep="*"))
tmp		<- setNames(tmp$segment, tmp$household_code)
pan$segment 	<- tmp[as.character(pan$household_code)]

# save(pan, hh_dist, file = "hh_fmt.rdata")
save(pan, hh_dist, hh_dpt, hh_month, fmt_dpt,file = paste("hh_month_exp_", as.character(Sys.Date()),".rdata", sep=""))

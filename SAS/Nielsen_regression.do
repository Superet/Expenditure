cd "E:\Users\ccv103\Desktop"
use "E:\Users\ccv103\Documents\Research\Store switching\SAS_temp\Nielsen\all_hh_month_exp.dta", clear
set more off

sort household_code year month
gen ymonth = ym(year, month)
format ymonth %tm
tabulate month, generate (dum_m)
xi i.ymonth

* Bad city; 
areg dol dum_m1-dum_m12 _Iymonth_540-_Iymonth_635 if mkt_type=="Badcity", absorb(household_code) cluster(market)
outreg2 _Iymonth_540-_Iymonth_635 using stata_badcity, excel replace noaster noparen

foreach x of varlist net_dol dol_purchase dolp_convenience_store dolp_dollar_store dolp_drug_store dolp_discount_store dolp_grocery dolp_health_food_store dolp_warehouse_club   q n_mod storable_quant nonstr_quant num_day num_trip pct_coupon unit_price2004{
	areg `x' dum_m1-dum_m12 _Iymonth_540-_Iymonth_635 if mkt_type=="Badcity", absorb(household_code) cluster(market)
	outreg2 _Iymonth_540-_Iymonth_635 using stata_badcity, excel append noaster noparen
}


* Good city; 
areg dol dum_m1-dum_m12 _Iymonth_540-_Iymonth_635 if mkt_type=="Goodcity", absorb(household_code) cluster(market)
outreg2 _Iymonth_540-_Iymonth_635 using stata_goodcity, excel replace noaster noparen

foreach x of varlist net_dol dol_purchase dolp_convenience_store dolp_dollar_store dolp_drug_store dolp_discount_store dolp_grocery dolp_health_food_store dolp_warehouse_club   q n_mod storable_quant nonstr_quant num_day num_trip pct_coupon unit_price2004{
	areg `x' dum_m1-dum_m12 _Iymonth_540-_Iymonth_635 if mkt_type=="Goodcity", absorb(household_code) cluster(market)
	outreg2 _Iymonth_540-_Iymonth_635 using stata_goodcity, excel append noaster noparen
}



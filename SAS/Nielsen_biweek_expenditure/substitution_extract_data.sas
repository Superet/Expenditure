libname mylib "E:\Users\ccv103\Documents\Research\Store switching\SAS_temp\Nielsen";
libname mysel "E:\Users\ccv103\Documents\Research\Store switching\SAS_temp\Nielsen\Work";
libname mainex "E:\Users\ccv103\Documents\Research\Store switching\SAS_temp\Nielsen\Work\main";

options compress=yes reuse=yes; 
proc datasets library=WORK kill noprint; run;

%let dirname= E:\Users\ccv103\Documents\Research\Nielsen\nielsen_extracts\HMS;

*****************; 
* Organize data *; 
*****************; 
/*
Four main datasets: 
1. demo_init: the initial states of households. KEY (household_code)
2. demo_annual: the annual profiles of households. KEY (household_code, year)
3. hh_exp: the biweekly expenditure data of households. KEY (household_code, biweek)
4. purchase: individual purchases of products. KEY (household_code, date, UPC)
*/

* Read in CPI data;
PROC IMPORT OUT= CPI
            DATAFILE= "\\tsclient\Resear1\Store switching\processed data\CPI.csv" 
            DBMS=csv REPLACE;
     		GETNAMES=YES; 
RUN;

proc sql noprint;
	select annual into:cpi_base from cpi where year = 2004; 
quit; 
%put &cpi_base; 

data cpi(keep = year annual cpi); 
	set cpi; 
	cpi = annual/&cpi_base; 
run;

* Household demographics; 
data demo_annual; 
	set mylib.panelists(keep = household_code panel_year income_real income_midvalue famsize scantrack_market_descr); 
run;

data tmp1; set demo_annual; run; 

* Adjust income with CPI; 
proc sql noprint; 
	create table demo_annual as 
	select A.*, income_midvalue/cpi as income_cpi
	from tmp1 as A left join cpi as B
	on A.panel_year = B.year; 
quit; 

* Segment households using intial states;
proc sort data = demo_annual; by household_code panel_year; run;

data demo_init; 
	set demo_annual; 
	by household_code; 
	if not first.household_code then delete; 
run;

* The distribution of intial income; 
proc freq data = demo_init; table income_real; run; 
data demo_init;
	set demo_init; 
	first_income = 'T1'; 
	if income_real >17 & income_real <=23 then first_income = 'T2'; 
	if income_real >23 then first_income = 'T3'; 
run; 

*--------------------------------------------------------*; 
* Compute percentile of outcome variables to hh_exp -- expenditure and share; 
proc sql noprint;
	create table tmp_rank as 
	select household_code, biweek, year, dol, 
		DOL_Convenience_Store, DOL_Discount_Store, DOL_Dollar_Store, 
		DOL_Drug_Store, DOL_Grocery, DOL_Warehouse_Club,
		DOL_Convenience_Store/dol as SHR_Convenience_Store, DOL_Discount_Store/dol as SHR_Discount_Store, 
		DOL_Dollar_Store/dol as SHR_Dollar_Store, DOL_Drug_Store/dol as SHR_Drug_Store, 
		DOL_Grocery/dol as SHR_Grocery, DOL_Warehouse_Club/dol as SHR_Warehouse_Club,
		num_day, num_trip, num_food_module, num_noned_module, food_quant, nonedible_quant, 
		scantrack_market_descr, income_real, income_midvalue
	from mylib.Hh_biweek_exp_merge
	where household_code in (select household_code from demo_init); 
quit; 

* Append initial states to hh_exp and adjust income with cpi;
proc sql noprint;
	create table tmp as 
	select A.*, B.first_income, B.famsize
	from tmp_rank as A inner join demo_init as B
	on A.household_code = B.household_code;
	
	create table hh_exp as 
	select A.*, income_midvalue/B.cpi as income_cpi
	from tmp as A left join cpi as B
	on A.year = B.year
	order by household_code, biweek; 
	
	drop table tmp_rank; 
quit; 

* NOTE: we use CPI adjusted income; 
data hh_exp;
	set hh_exp; 
	week_first 	= mdy(1,1,2004) + (biweek-1)*14; 
	month 		= month(week_first);
	quant 		= food_quant + nonedible_quant; 
	num_mod		= num_food_module + num_noned_module; 
	Rc 			= 1 * (year >= 2008); 
run; 

*--------------------------------------------------------*; 
* Compute percentile of outcome variables in purchase data -- price tier and size tier; 
data dict; 
	set mylib.products(keep = upcv department_descr product_group_descr product_module_descr
		 	brand_descr multi size size1_amount ); 
run; 

* Price tier of UPC within category; 
proc sql noprint;
	create table tmp1 as 
	select upcv, mean(price) as price
	from mainex.purchases
	group by upcv;
	
	create table tmp as 
	select A.*, B.product_module_descr, B.size, price/size as unit_price
	from tmp1 as A left join dict as B
	on A.upcv = B.upcv
	order by product_module_descr; 
quit; 

* Rank price within category; 
proc rank data=tmp out=myprod ties=low PERCENT;
	by product_module_descr;
	var unit_price size;
	ranks PCTL_UPRICE PCTL_SIZE; 
run;

* Append price tier and household segment to purchase data; 
proc sql noprint;
	create table tmp_purch as 
	select household_code, A.year, purchase_date, A.biweek, upcv, quantity, price
	from mainex.purchases as A 
	inner join (select * from mainex.trips where household_code in (select household_code from demo_init)) as B
	on A.trip_code_uc = B.trip_code_uc;
	
	create table tmp as 
	select A.*, month(purchase_date) as month, 1 * (year >= 2008) as Rc, 
		B.product_module_descr, B.size, B.unit_price, B.PCTL_UPRICE, B.PCTL_SIZE
	from tmp_purch as A left join myprod as B
	on A.upcv = B.upcv
	order by household_code; 
	
	drop table tmp_purch; 

/*NOTE: SAS may have insufficient space for this step. In that case, I use the following the data step.*/
	create table purchases as 
	select A.*, B.first_income, B.famsize
	from tmp as A left join demo_init as B
	on A.household_code = B.household_code
	order by household_code, purchase_date; 
quit; 

/*data purchases;
	merge tmp demo_init; 
	by household_code; 
run;*/

proc datasets noprint; delete tmp tmp_purch tmp1; run; 

* Aggregate size tier and price tier from purchase data to biweek level; 
proc sql noprint;
	* Combine multiple purchases from the same category; 
	create table tmp as 
	select household_code, biweek, product_module_descr, 
		sum(quantity*PCTL_UPRICE)/sum(quantity) as PCTL_UPRICE, 
		sum(quantity*PCTL_SIZE)/sum(quantity) as PCTL_SIZE
	from purchases
	group by household_code, biweek, product_module_descr; 
	
	* Append category weights;
	create table tmp1 as 
	select A.*, B.weight
	from tmp as A left join mylib.basket as B
	on A.product_module_descr = B.product_module_descr;
	
	* Compute weighted average of size/price tier during biweek; 
	create table tmp2 as 
	select household_code, biweek, sum(PCTL_UPRICE*weight)/sum(weight) as PCTL_UPRICE, 
			sum(PCTL_SIZE*weight)/sum(weight) as PCTL_SIZE
	from tmp1 
	group by household_code, biweek; 
	
	* Add price/size tier to hh_exp; 
	create table tmp3 as 
	select A.*, B.PCTL_UPRICE, B.PCTL_SIZE
	from hh_exp as A left join tmp2 as B
	on A.household_code = B.household_code and A.biweek = B.biweek; 
quit;
data hh_exp; set tmp3; run; 

proc datasets noprint; delete tmp tmp1 tmp2 tmp3 tmp_rank; run; 

*****************************;
* Export data to csv format *;
*****************************;
PROC EXPORT DATA= hh_exp
            OUTFILE= "E:\Users\ccv103\Desktop\substitution_hh_exp.csv"
            DBMS=csv REPLACE;
     PUTNAMES=YES;
RUN;

/*
# R code to import as .rdata
setwd("E:/Users/ccv103/Desktop")
hh_exp	<- read.csv("substitution_hh_exp.csv", header = T)
head(hh_exp)
save(hh_exp, file = "substitution_hh_exp.rdata")

*/
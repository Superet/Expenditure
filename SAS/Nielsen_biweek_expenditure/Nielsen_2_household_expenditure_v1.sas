libname mylib "E:\Users\ccv103\Documents\Research\Store switching\SAS_temp\Nielsen";
libname mysel "E:\Users\ccv103\Documents\Research\Store switching\SAS_temp\Nielsen\Work";
libname mainex "E:\Users\ccv103\Documents\Research\Store switching\SAS_temp\Nielsen\Work\main";

%let dirname= E:\Users\ccv103\Documents\Research\Nielsen\nielsen_extracts\HMS;
options compress=yes reuse=yes;


*************;
* Read data *;
*************;
* Read in codebook;
PROC IMPORT OUT= mysel.my_codebook
            DATAFILE= "\\tsclient\Resear1\Data\Nielsen\reference_documentation\code_book.csv" 
            DBMS=csv REPLACE;
     		GETNAMES=YES; 
RUN;


* Read in the data of deleted UPC;
PROC IMPORT OUT= my_delete_upc
            DATAFILE= "\\tsclient\Resear1\Store switching\processed data\deleted_upc.csv" 
            DBMS=csv REPLACE;
     		GETNAMES=YES; 
RUN;

* Read in CPI data; 
PROC IMPORT OUT= cpi
            DATAFILE= "\\tsclient\Resear1\Store switching\processed data\CPI.csv" 
            DBMS=csv REPLACE;
     		GETNAMES=YES; 
RUN;

* Classify product department by food and non-edible;
data products;
	set mylib.products;
	* Classify food and non-edible departments; 
	Dpmt_type = "Non-edible";
	if department_descr in ("ALCOHOLIC BEVERAGES","DAIRY","DELI","DRY GROCERY","FRESH PRODUCE","FROZEN FOODS","PACKAGED MEAT","MAGNET DATA") 
	then do;
		Dpmt_type = "Food";
		Food_module = product_module_descr;
	end;
	else Nonedible_module = product_module_descr;
	* Classify storable and perishable departments; 
	Storable = 1;
	if department_descr in ("DAIRY","DELI","FRESH PRODUCE","MAGNET DATA") then Storable = 0;
	if department_descr = "DRY GROCERY" and substr(product_module_descr,1,6) = "BAKERY" then Storable = 0; 
	if Storable = 1 then storable_module = product_module_descr;
	if Storable = 0 then nonstorable_module = product_module_descr;
	if my_magnet = 1 then magnet_module = product_module_descr;
run;

* NOTE: deleted UPCs are already dropped in mylib.products; 
* Read master product data;
%let fpath = &dirname\Master_Files\Latest\products.tsv;
%put &fpath;
PROC IMPORT OUT= WORK.orig_products
            DATAFILE= "&fpath" 
            DBMS=dlm REPLACE;
			DELIMITER ='09'x;
     		DATAROW=2; 
			GUESSINGROWS=5000;
RUN;

* Find the UPC for deleted product modules;
proc sql noprint;
	create table tmp1 as 
	select  A.product_module_descr, B.upc, B.upc_ver_uc, cats(B.upc, "-", B.upc_ver_uc) as upcv, B.upc_descr
	from (select * from my_delete_upc where upc_descr^="") as A inner join orig_products as B
	on A.upc_descr=B.upc_descr;
	
	create table tmp2 as 
	select A.product_module_descr, B.upc, B.upc_ver_uc, cats(B.upc, "-", B.upc_ver_uc) as upcv, B.upc_descr
	from (select * from my_delete_upc where upc_descr="") as A inner join orig_products as B
	on A.product_module_descr=B.product_module_descr;
quit;

data mysel.my_delete_upc; set tmp1 tmp2; run;
proc datasets noprint; delete orig_products; run;

* Define the retailer scope; 
data mysel.my_channel;
	input channel_type $40.	channel_group $10.;
	datalines;
	Grocery								Grocery
	Drug Store							Grocery
	Dollar Store						Grocery
	Convenience Store					Grocery
	Discount Store						Grocery
	Warehouse Club						Grocery
	;
run;

* Read master retailer data;
%let fpath = &dirname\Master_Files\Latest\retailers.tsv;
%put &fpath;
PROC IMPORT OUT= WORK.orig_retailers
            DATAFILE= "&fpath" 
            DBMS=dlm REPLACE;
			DELIMITER ='09'x;
     		DATAROW=2; 
			GUESSINGROWS=5000;
RUN;

* NOTE: we extract purchase data from only grocery stores;  
proc sql noprint;
	create table mysel.retailers as
	select A.*, B.channel_group
	from orig_retailers as A inner join mysel.my_channel as B
	on A.channel_type = B.channel_type;
quit;

* Append products attributes to purchase data; 
proc sql noprint; 	
	create table tmp as 
	select A.*, B.size_index, B.size, B.Dpmt_type, B.Food_module, B.Nonedible_module, 
			B.storable, B.storable_module, B.nonstorable_module, B.product_module_descr, B.my_magnet, B.weight
	from mainex.purchases as A inner join 
	(select A.*, B.weight from products as A inner join mylib.basket as B on A.product_module_descr = B.product_module_descr) as B
	on A.upcv = B.upcv; 
quit; 

* Summarize a shopping basket in a trip; 
proc sql noprint; 
	create table tmp1 as 
	select trip_code_uc, sum(total_price_paid - coupon_value) as dol_purchases, sum(coupon_value) as coupon_dol, 
			sum((total_price_paid - coupon_value)*my_magnet) as dol_magnet, sum(ifn(my_magnet=1, size_index*weight*quantity, 0)) as magnet_quant,
			count(distinct Food_module) as num_food_module, count(distinct Nonedible_module) as num_noned_module,
			sum(ifn(Dpmt_type="Food", size_index*weight*quantity, 0)) as food_quant, sum(ifn(Dpmt_type^="Food", size_index*weight*quantity, 0)) as nonedible_quant, 
			count(distinct storable_module) as num_storable_module, count(distinct nonstorable_module) as num_nonstr_module,
			sum(storable*size_index*weight*quantity) as storable_quant, sum((1-storable)*size_index*weight*quantity) as nonstr_quant
	from tmp
	group by trip_code_uc; 
			
	create table trips(drop = old_trip) as
	select *, ceil((purchase_date - MDY(12,31,2003))/14) as biweek
	from mainex.trips as A left join tmp1(rename = (trip_code_uc = old_trip)) as B
	on A.trip_code_uc = B.old_trip; 
quit; 

* Summarize a shopping basket in bi-weekly basis; 
proc sql noprint; 
	create table biweek_basket as 
	select household_code, min(year) as year, biweek, sum(total_price_paid - coupon_value) as dol_purchases,
			count(unique(purchase_date)) as num_day, count(unique(trip_code_uc)) as num_trip, sum(coupon_value) as coupon_dol, 
			sum((total_price_paid - coupon_value)*my_magnet) as dol_magnet, sum(ifn(my_magnet=1, size_index*weight*quantity, 0)) as magnet_quant,
			count(distinct Food_module) as num_food_module, count(distinct Nonedible_module) as num_noned_module,
			sum(ifn(Dpmt_type="Food", size_index*quantity*weight, 0)) as food_quant, sum(ifn(Dpmt_type^="Food", size_index*quantity*weight, 0)) as nonedible_quant,
			count(distinct storable_module) as num_storable_module, count(distinct nonstorable_module) as num_nonstr_module,
			sum(storable*size_index*quantity*weight) as storable_quant, sum((1-storable)*size_index*quantity*weight) as nonstr_quant
	from (select A.*, B.household_code, B.purchase_date from tmp as A inner join trips as B on A.trip_code_uc = B.trip_code_uc) 
	group by household_code, biweek
	order by household_code, biweek; 
quit;

proc datasets noprint; delete tmp tmp1 tmp2; run;
	
************************;
* Household selection * ;
************************;
/*Selection critetion:
1. The maximum inter-purchase gap is less than 30 days; 
2. Households have at least 6 observations before and after the recession; */

%let gap_thresh = 30;
%let min_n		= 6;

* Compute the inter-purchase date;
proc sort data = trips; by household_code purchase_date; run;
data tmp;
	set trips;
	inter_date = dif(purchase_date);
	by household_code; 
	if first.household_code then inter_date = .;
	recession = 0; 
	recdate = MDY(1,1,2008);
	if purchase_date > recdate then recession = 1; 
	biweek = ceil((purchase_date - MDY(12,31,2003))/14);
run;

* Find the maximum interpurchase date for each household and plot the histogram;
* NOTE: households may move to a different geo market; 
proc sql noprint;
	create table tmp1 as 
	select scantrack_market_descr, household_code, max(inter_date) as gap, sum(biweek*(1-recession) > 0) as n1, sum(biweek*recession > 0) as n2
	from tmp 
	group by scantrack_market_descr, household_code;
	
	create table tmp2 as 
	select household_code, max(gap) as gap, count(distinct scantrack_market_descr) as n_mkt, sum(n1) as n1, sum(n2) as n2
	from tmp1
	group by household_code;
quit;

* Check how many households will be dropped; 
proc univariate data=tmp2 noprint; histogram gap; run;
proc sql; 
	select count(gap) from tmp2 where gap <= &gap_thresh; 
	select sum(n_mkt>1) from tmp2;
	select sum(n1 > &min_n and n2 > &min_n) from tmp2;
quit;

* Drop the households who do not satisfy our selection criterion; 
proc sql noprint;
	create table mysel.keeppan as 
	select household_code 
	from tmp2 
	where gap <= &gap_thresh and n_mkt=1 and (n1 > &min_n and n2 > &min_n)
	order by household_code;
quit;

*******************************************;
* Construct stone price index every biweek *;
*******************************************;
* COMPUTE STONE PRICE INDEX FROM RETAIL FORMAT-LEVEL PRICE INDEX;
proc sql noprint;
	* Expenditure share of retail formats from all households in a scantrack_market_descr each biweek;
	create table tmp as 
	select *, revenue/sum(revenue) as share
	from (select  scantrack_market_descr, min(year) as year, biweek, channel_type,
			sum(total_spent) as revenue from trips
			group by  scantrack_market_descr, biweek, channel_type)
	group by  scantrack_market_descr, biweek;
	
	* Compute stone price index;
	create table mainex.biweek_stone_price as 
	select  scantrack_market_descr, min(year) as year, biweek, 
			sum(share*log(bsk_price_paid)) as stone_price
	from (select A.*, B.share from mylib.format_biweek_price as A inner join tmp as B 
		  on A.scantrack_market_descr=B.scantrack_market_descr 
			 and A.biweek=B.biweek and A.channel_type=B.channel_type)
	group by  scantrack_market_descr, biweek
	order by  scantrack_market_descr, biweek;
quit;

proc datasets noprint; delete tmp tmp1 tmp2; run;

**************************************************************;
* Construct the data for the process of biweekly expenditure *;
**************************************************************;
*---------------------------*;
* Fill in nonshopping weeks *;
proc sql noprint;
	* Aggregate daily trip data to biweekly expenditure data; 
	create table tmp as 
	select household_code, min(year) as year, biweek, sum(total_spent) as DOL, sum(net_spent) as net_dol, sum(dol_purchases) as dol_purchases
	from (select * from trips where household_code in (select household_code from mysel.keeppan))
	group by household_code, biweek;
	
	create table tmp1 as 
	select household_code, min(biweek) as min_n, max(biweek) as max_n
	from tmp
	group by household_code
	order by household_code; 
	
	create table biweek_year as 
	select biweek, min(year) as year
	from tmp
	group by biweek
	order by biweek; 
quit;

data tmp1(drop=min_n max_n);
	set tmp1;
	by household_code;
	biweek = min_n - 1;
	if first.household_code then do; 	
		do while(biweek<max_n);
			biweek + 1; output;
		end;
	end;
run;

* Merge together expenditure with complete biweek index; 	
proc sql noprint;
	create table tmp_exp as 
	select B.*, A.DOL, A.net_dol, A.dol_purchases
	from tmp as A right join tmp1 as B
	on A.household_code=B.household_code and A.biweek=B.biweek
	order by biweek;
quit;

* Add year vairable; 
data tmp_exp; 
	merge tmp_exp biweek_year;
	by biweek;
proc sort; by household_code biweek; 
run; 

* Check the non-shopping weeks; 
proc sql; select nmiss(dol)/count(biweek) from tmp2; quit;
proc sql noprint;
	create table tmp as 
	select household_code, nmiss(dol_purchases)/count(biweek) as n_missbiweek
	from tmp_exp
	group by household_code
	order by household_code;
quit;
proc univariate data = tmp; histogram n_missbiweek; run;

proc datasets noprint; delete tmp tmp1 tmp2; run;

*---------------------------------------------------------------------*;
* Compute household expenditure share across retail format each biweek ;
proc sql noprint;
	create table tmp as
	select scantrack_market_descr, household_code, min(year) as year, biweek, channel_type, sum(net_spent) as DOL, sum(dol_purchases) as dolpurch
	from (select * from trips where household_code in (select household_code from mysel.keeppan))
	group by scantrack_market_descr, household_code, biweek, channel_type;
quit;

proc transpose data = tmp out = tmp_share prefix=DOLP_;
	by scantrack_market_descr household_code biweek;
	id channel_type;
	var dolpurch;
run;

data tmp_share; 
	merge tmp_share(drop=_NAME_);
	by scantrack_market_descr household_code biweek;
	array dolvar{*} DOL:;
	do i = 1 to dim(dolvar);
		if dolvar[i]=. then dolvar[i] = 0;
	end;
	drop i;
run;

*---------------------------------------*;
* Constrct data of household demographics;
data panelists; set mainex.panelists; run;
proc sort data=panelists; by household_code panel_year; run;

* Adjust the income reporting; 
data tmp; 
	set panelists(keep = household_code panel_year household_income);
	income_year = panel_year - 2;
	drop panel_year; 
run;

* NOTE: we lose data from the last two years because of the lag income reporting;
proc sql noprint;
	create table tmp1 as 
	select A.*, B.household_income as income_real
	from panelists as A inner join tmp as B
	on A.household_code = B.household_code and A.panel_year = B.income_year
	order by household_code, panel_year;
quit;

* The initial income status of households; 
data tmp2(keep = household_code first_income); 
	set tmp1;
	by household_code panel_year; 
	first_income = income_real; 
	if not first.panel_year then delete;
run; 

data tmp;
	merge tmp1(keep = household_code panel_year scantrack_market_descr 
			income_real household_size type_of_residence male_head_employment female_head_employment age_and_presence_of_children)
		  tmp2; 
	by household_code; 		
	famsize = "Three+";
	if household_size = 2 then famsize = "Two";
	if household_size = 1 then famsize = "Single";
	condo = 0; 
	if mod(type_of_residence, 2)=0 then condo = 1;
	employed = 0;
	if (male_head_employment <= 3 | female_head_employment <=3) then employed = 1;
	NumChild = 0;
	if age_and_presence_of_children <= 3 then NumChild = 1;
	if age_and_presence_of_children >3 and age_and_presence_of_children <= 7 then NumChild = 2;
run;

* Add the middle value of income range; 
proc sql noprint;
	create table panelists as
	select A.*, B.mid_range as income_midvalue
	from tmp as A left join (select code_value, description, mid_range from mysel.my_codebook where variable='household_income') as B
	on A.income_real = B.code_value
	order by household_code;	
quit;

proc datasets noprint; delete tmp tmp1 tmp2; run; 

*-------------------------------------------------------------------------*;
* Merge household demongraphics and price index to biweekly expenditure data; 
proc sql noprint;
	* Combine expenditure share, and basket data together; 
	create table tmp1(drop = old_household_code old_biweek dol_purchases year) as 
	select * 
	from tmp_share as A left join biweek_basket(rename = (household_code = old_household_code biweek = old_biweek)) as B 
	on A.household_code = B.old_household_code and A.biweek = B.old_biweek; 

	* Combine complete shopping expenditure and expenditure share; 
	create table tmp2(drop = old_household_code old_biweek)  as 
	select *
	from tmp_exp as A left join tmp1(rename = (household_code = old_household_code biweek = old_biweek)) as B
	on A.household_code = B.old_household_code and A.biweek = B.old_biweek;
	
	* Combine price index; 
	create table tmp3 as 
	select A.*, B.stone_price
	from tmp2 as A left join mainex.biweek_stone_price as B
	on A.scantrack_market_descr = B.scantrack_market_descr and A.biweek = B.biweek;
	
	* Merge in household demographics; 
	create table hh_biweek_exp_merge(drop = household_code_old panel_year_old) as 
	select *
	from tmp3 as A inner join 
	panelists(rename=(household_code = household_code_old panel_year=panel_year_old) drop= scantrack_market_descr) as B
	on A.household_code = B.household_code_old and A.year = B.panel_year_old
	order by household_code, biweek;
quit; 

proc datasets noprint; delete tmp1 tmp2 tmp3; run;

*---------------------------------*;
* Adjust all dollar measure for CPI; 
proc sql noprint;
	select Annual into: price_base from cpi where year = 2004; 
	
	create table tmp as 
	select A.*, B.Annual/&price_base as cpi
	from hh_biweek_exp_merge as A left join cpi as B
	on A.year = B.year
	order by household_code, year, biweek;
quit;

* Adjust CPI;
data hh_biweek_exp_merge(drop=i);
	set tmp;
	array dolvar{*} dol: coupon_dol stone_price net_dol;
	do i=1 to dim(dolvar);
		dolvar[i] = dolvar[i]/cpi;
	end;
	recession = 0;
/*	rec_week = ceil((MDY(11,30,2007)- MDY(12,31,2003))/14); 
	if biweek > rec_week then recession = 1;*/
	if year >= 2008 then recession = 1; 
run;	

* Check if the data has duplicated by key;
* NOTE: the data is not unique in keys because some households move; 
proc sort data=hh_biweek_exp_merge nodupkeys out=tmp; 
by household_code biweek;
run;

proc contents data = hh_biweek_exp_merge; run;

* Save data; 
data mylib.panelists; set panelists; run;
data mylib.hh_biweek_exp_merge; set hh_biweek_exp_merge; run;

proc datasets noprint; delete tmp tmp_share tmp_exp; run;

****************************;
* Export the data to excel *;
****************************;
PROC EXPORT DATA=mylib.panelists
			OUTFILE= "\\tsclient\Resear1\Store switching\processed data\2_panelists.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;

PROC EXPORT DATA=mylib.hh_biweek_exp_merge
			OUTFILE= "\\tsclient\Resear1\Store switching\processed data\2_hh_biweek_exp_merge.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;

/*
proc datasets library=mysel kill noprint;
proc datasets library=mainex kill noprint;
run;
*/



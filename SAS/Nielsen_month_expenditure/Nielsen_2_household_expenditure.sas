libname mylib "E:\Users\ccv103\Documents\Research\Store switching\SAS_temp\Nielsen";
libname mysel "E:\Users\ccv103\Documents\Research\Store switching\SAS_temp\Nielsen\Work";
libname mainex "E:\Users\ccv103\Documents\Research\Store switching\SAS_temp\Nielsen\Work\main";

%let dirname= E:\Users\ccv103\Documents\Research\Nielsen\nielsen_extracts\HMS;

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
	Health Food Store					Grocery
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

* Read trip data;
%macro read_groc_fun;
	%do i=2004 %to 2012;
		* Read panelists data for geographic scantrack_market_descr information; 
		%let fpath = &dirname\&i\Annual_Files\panelists_&i..tsv;
		%put &fpath;
		PROC IMPORT OUT= WORK.panelists&i
		            DATAFILE= "&fpath" 
		            DBMS=dlm REPLACE;
					DELIMITER ='09'x;
		     		DATAROW=2; 
					GUESSINGROWS=5000;
		RUN;
		
		* Import trip data;
		%let fpath = &dirname\&i\Annual_Files\trips_&i..tsv;
		%put &fpath;
		PROC IMPORT OUT= WORK.trips&i
		            DATAFILE= "&fpath" 
		            DBMS=dlm REPLACE;
					DELIMITER ='09'x;
		     		DATAROW=2; 
		RUN;
		
		* Import purchase data;
		%let fpath = &dirname\&i\Annual_Files\purchases_&i..tsv;
		%put &fpath;
		PROC IMPORT OUT= WORK.purchases&i
		            DATAFILE= "&fpath" 
		            DBMS=dlm REPLACE;
					DELIMITER ='09'x;
		     		DATAROW=2; 
		RUN;
		
		data panelists&i(drop = projection_factor_magnet_char);
			length scantrack_market_descr $50.;
			set panelists&i(rename = (projection_factor_magnet = projection_factor_magnet_char)); 
			projection_factor_magnet = input(projection_factor_magnet_char, best12.);
			magnet_flag = 0;
			if projection_factor_magnet^=. then magnet_flag = 1;
			drop wic_indicator_current wic_indicator_ever_notcurrent Member:;
		run;
		
		proc sql noprint; 
		 	* Add households market to trip data;
			create table tmp1 as 
			select A.*, year(A.purchase_date) as year, B.scantrack_market_descr, 
					B.projection_factor_magnet, B.magnet_flag 
			from trips&i as A inner join panelists&i as B
			on A.household_code = B.household_code;	
				
			* Restrict the trip data to the selected retailers;
			create table tmp2 as
			select A.*, B.channel_group, B.channel_type
			from tmp1 as A inner join mysel.retailers as B
			on A.retailer_code = B.retailer_code
			order by trip_code_uc;	
			
			* Subset purchases data within the restrictions above;
			create table tmp3 as 
			select A.*, cats(upc, "-", upc_ver_uc) as upcv, (total_price_paid-coupon_value) as dol_paid, total_price_paid/quantity as price,
					B.household_code, B.scantrack_market_descr, B.channel_type, B.retailer_code,
					year(B.purchase_date) as year, month(B.purchase_date) as month, B.purchase_date
			from purchases&i as A inner join tmp2 as B
			on A.trip_code_uc = B.trip_code_uc;	

			* A subset of purchases of deleted products, then aggregated by trip; 
			create table tmp4 as 
			select trip_code_uc, sum(total_price_paid-coupon_value) as delete_dol
			from (select * from tmp3 where upcv in(select upcv from mysel.my_delete_upc))
			group by trip_code_uc;
			
			* Merge in the delete dollar for each trip;
			create table trips&i as
			select A.*, B.delete_dol
			from tmp2 as A left join tmp4 as B
			on A.trip_code_uc = B.trip_code_uc
			order by trip_code_uc;
		quit;
		
		* Summarize the basket each month;
		* NOTE: products data set already drops the deleted UPC;
		proc sql noprint;
			create table tmp as
			select A.*, B.size_index, B.size, B.Dpmt_type, B.Food_module, B.Nonedible_module, 
					B.storable, B.storable_module, B.nonstorable_module, B.product_module_descr, B.my_magnet
			from tmp3 as A inner join products as B on A.upcv = B.upcv;
			
			* Add the module share (varying by market );
			create table tmp1 as
			select A.*, B.share_paid
			from tmp as A left join mylib.module_wallet as B
			on A.product_module_descr = B.product_module_descr and A.scantrack_market_descr = B.scantrack_market_descr;
			
			create table tmp_trip as 
			select trip_code_uc, sum(total_price_paid - coupon_value) as dol_purchases, sum(coupon_value) as coupon_dol, 
					sum((total_price_paid - coupon_value)*my_magnet) as dol_magnet, sum(ifn(my_magnet=1, size_index*share_paid*quantity, 0)) as magnet_quant,
					count(distinct Food_module) as num_food_module, count(distinct Nonedible_module) as num_noned_module,
					sum(ifn(Dpmt_type="Food", size_index*share_paid*quantity, 0)) as food_quant, sum(ifn(Dpmt_type^="Food", size_index*share_paid*quantity, 0)) as nonedible_quant, 
					count(distinct storable_module) as num_storable_module, count(distinct nonstorable_module) as num_nonstr_module,
					sum(storable*size_index*share_paid*quantity) as storable_quant, sum((1-storable)*size_index*share_paid*quantity) as nonstr_quant
			from tmp1
			group by trip_code_uc
			order by trip_code_uc;

			create table month_basket&i as 
			select household_code, year, month, sum(total_price_paid - coupon_value) as dol_purchases,
					count(unique(purchase_date)) as num_day, count(unique(trip_code_uc)) as num_trip, sum(coupon_value) as coupon_dol, 
					sum((total_price_paid - coupon_value)*my_magnet) as dol_magnet, sum(ifn(my_magnet=1, size_index*share_paid*quantity, 0)) as magnet_quant,
					count(distinct Food_module) as num_food_module, count(distinct Nonedible_module) as num_noned_module,
					sum(ifn(Dpmt_type="Food", size_index*quantity*share_paid, 0)) as food_quant, sum(ifn(Dpmt_type^="Food", size_index*quantity*share_paid, 0)) as nonedible_quant,
					count(distinct storable_module) as num_storable_module, count(distinct nonstorable_module) as num_nonstr_module,
					sum(storable*size_index*quantity*share_paid) as storable_quant, sum((1-storable)*size_index*quantity*share_paid) as nonstr_quant
			from tmp1
			group by household_code, year, month
			order by household_code, year, month;
			
			* NOTE: we would have duplicated december basket due to the inconsistence between calendar year and nielsen year; 
			* 		Therefore, we keep December and January data separately for later aggregation; 
			create table dec_jan&i as
			select *
			from tmp1
			where month in (1, 12);
		quit;
		
		data purchases&i; set tmp1(drop=Dpmt_type Food_module Nonedible_module storable storable_module nonstorable_module); run;
		
		data trips&i; 
			merge trips&i tmp_trip; 
			by trip_code_uc; 
			if delete_dol=. then delete_dol=0;
			net_spent = total_spent - delete_dol;
			month	= month(purchase_date);
			week 	= week(purchase_date, 'u');
			year 	= year(purchase_date);
			ymonth	= put(purchase_date, monyy7.);
			yweek 	= put(purchase_date, weeku5.);
		run;
			
		proc datasets noprint; delete tmp tmp1 tmp2 tmp3 tmp_trip; run;
	%end;
	
	* Combein the panelists data together; 
	data mainex.panelists;
		set panelists2004 - panelists2012;
	run;
	
	* Combine trips data together;
	data mainex.trips;
		set trips2004 - trips2012;
	run;
	
	* Combine monthly basket for each household;
	data month_basket_tmp;
		set month_basket2004 - month_basket2012;
		where month not in (1, 12);
	run;
	
	data dec_jan;
		set dec_jan2004 - dec_jan2012;
	run;
	
	proc sql noprint;
		create table tmp as 
		select household_code, year, month, sum(total_price_paid - coupon_value) as dol_purchases,
				count(unique(purchase_date)) as num_day, count(unique(trip_code_uc)) as num_trip, sum(coupon_value) as coupon_dol, 
				sum((total_price_paid - coupon_value)*my_magnet) as dol_magnet, sum(ifn(my_magnet=1, size_index*share_paid*quantity, 0)) as magnet_quant,
				count(distinct Food_module) as num_food_module, count(distinct Nonedible_module) as num_noned_module,
				sum(ifn(Dpmt_type="Food", size_index*share_paid*quantity, 0)) as food_quant, sum(ifn(Dpmt_type^="Food", size_index*share_paid*quantity, 0)) as nonedible_quant,
				count(distinct storable_module) as num_storable_module, count(distinct nonstorable_module) as num_nonstr_module,
				sum(storable*size_index*quantity*share_paid) as storable_quant, sum((1-storable)*size_index*quantity*share_paid) as nonstr_quant
		from dec_jan
		group by household_code, year, month
		order by household_code, year, month;
	quit;
	
	data mainex.month_basket;
		set month_basket_tmp tmp;
	proc sort; by household_code year month;
	run;
	
	data mainex.purchases;
		set purchases2004 - purchases2012;
	run;
		
	proc datasets noprint; 
		delete panelists2004 - panelists2012 trips2004 - trips2012 
		month_basket2004 - month_basket2012 purchases2004 - purchases2012
		dec_jan2004 - dec_jan2012 tmp month_basket_tmp dec_jan; 
	run;
%mend read_groc_fun;

%read_groc_fun;

proc contents data=mainex.trips; run;
proc contents data=mainex.month_basket; run;

* Drop the data from 2003 due to the discrepency btw panel_year and calendar year; 
data mainex.trips;
	set mainex.trips;
	where year>2003;
run;

proc datasets library=WORK kill noprint; run;


* Keep the data from all households; 
data mainex.all_trips; set mainex.trips; run;
data mainex.all_panelists; set mainex.panelists; run;
data mainex.all_month_basket; set mainex.month_basket; run;

************************;
* Household selection * ;
************************;
* Compute the inter-purchase date;
proc sort data = mainex.all_trips; by household_code purchase_date; run;
data tmp;
	set mainex.all_trips;
	inter_date = dif(purchase_date);
	by household_code; 
	if first.household_code then inter_date = .;
run;

* Find the maximum interpurchase date for each household and plot the histogram;
proc sql noprint;
	create table tmp1 as 
	select scantrack_market_descr, household_code, max(inter_date) as gap
	from tmp 
	group by scantrack_market_descr, household_code;
	
	create table tmp2 as 
	select household_code, max(gap) as gap, count(distinct scantrack_market_descr) as n_mkt
	from tmp1
	group by household_code;
quit;

%let gap_thresh = 30;
proc univariate data=tmp2 noprint; histogram gap; run;
proc sql; 
	select count(gap) from tmp2 where gap <= &gap_thresh; 
	select sum(n_mkt>1) from tmp2;
quit;

* Drop the households who have interpurchase date greater than the threshold;
proc sql noprint;
	create table mysel.keeppan as 
	select household_code 
	from tmp2 
	where gap <= &gap_thresh and n_mkt=1
	order by household_code;
quit;

* Restrict all the data to the selected households;
proc sql noprint;
	create table tmp1 as
	select *
	from mainex.panelists 
	where household_code in (select household_code from mysel.keeppan);
	
	create table tmp2 as
	select *
	from mainex.trips 
	where household_code in (select household_code from mysel.keeppan);

	create table tmp3 as
	select *
	from mainex.month_basket 
	where household_code in (select household_code from mysel.keeppan);
quit;

data mainex.panelists; 	set tmp1 ; run;
data mainex.trips; 		set tmp2; run;
data mainex.month_basket; set tmp3; run;

proc datasets noprint; delete tmp tmp1 tmp2 tmp3; run;

*******************;
* Variable format *;
*******************;
* The proportion of magnet households;
proc freq data=mainex.panelists;
	table magnet_flag;
	table scantrack_market_descr*magnet_flag;
	table panel_year*magnet_flag;
run;

* Check if a household being consistenly in the magnet system;
proc sql noprint;
	create table tmp as 
	select household_code, mean(magnet_flag) as magnet_flag
	from mainex.panelists
	group by household_code
	order by household_code;
quit;

proc univariate data=tmp; histogram magnet_flag; run;
data tmp; 
	set tmp; 
	magnet_const_flag = 0;
	if magnet_flag = 1 then magnet_const_flag = 1;
proc freq; table magnet_const_flag;
run;

* If the household is consistently in the magnet system, magnet_const_flag = 1;
proc sort data=mainex.panelists; by household_code; run;
data mainex.panelists;
	merge mainex.panelists tmp(drop = magnet_flag);
	by household_code;
run;

*******************************************;
* Construct stone price index every month *;
*******************************************;
* COMPUTE STONE PRICE INDEX FROM RETAIL FORMAT-LEVEL PRICE INDEX;
proc sql noprint;
	* Expenditure share of retail formats from all households in a scantrack_market_descr each month;
	create table tmp as 
	select *, revenue/sum(revenue) as share
	from (select  scantrack_market_descr, year, month, channel_type,
			sum(net_spent) as revenue from mainex.all_trips
			group by  scantrack_market_descr, year, month, channel_type)
	group by  scantrack_market_descr, year, month;
	
	* Compute stone price index;
	create table mainex.month_stone_price as 
	select  scantrack_market_descr, year, month, 
			sum(share*log(price_tag_norm_index)) as stone_price_tag, 
			sum(share*log(price_paid_norm_index)) as stone_price
	from (select A.*, B.share from mylib.format_month_price as A inner join tmp as B 
		  on A.scantrack_market_descr=B.scantrack_market_descr and A.year=B.year 
			 and A.month=B.month and A.channel_type=B.channel_type)
	group by  scantrack_market_descr, year, month
	order by  scantrack_market_descr, year, month;
quit;

/*
* RUN THE CODES BELOW TO COMPUTE STONE PRICE INDEX FROM AGGREGATING UPC *;
* Wallet share of upc;
proc sql noprint;
	create table tmp as 
	select upcv, sum(dol_paid) as dol_paid
	from mainex.month_upc_price 
	group by upcv
	order by upcv;
	
	create table mainex.wallet_share as
	select upcv, dol_paid, dol_paid/sum(dol_paid) as share
	from tmp
	order by upcv;
quit;

* Construct scantrack_market_descr leve stone price index;
proc sql noprint;
	create table tmp as 
	select A.*, B.share
	from mainex.month_upc_price as A left join mainex.wallet_share as B
	on A.upcv = B.upcv;

	create table mainex.month_stone_price as
	select  scantrack_market_descr, year, month, sum(share*log(mean_price)) as stone_price
	from tmp
	group by  scantrack_market_descr, year, month
	order by  scantrack_market_descr, year, month;
quit;
*/

*************************************************************;
* Construct the data for the process of monthly expenditure *;
*************************************************************;
* Check whether to have outside options and determine the threshold for household selection ;
proc sql noprint;
	* Aggregate daily trip data to monthly expenditure data; 
	create table tmp as 
	select household_code, year, month, sum(total_spent) as DOL, sum(net_spent) as net_dol, sum(dol_purchases) as dol_purchases
	from mainex.trips 
	group by household_code, year, month;

	create table ymonth as 
	select ymonth, mean(year) as year, mean(month) as month
	from mainex.trips
	group by ymonth
	order by year, month;
quit;

data ymonth;
	set ymonth;
	n = _N_; 
run;

proc sql noprint;
	create table tmp1 as 
	select household_code, min(n) as min_n, max(n) as max_n
	from (select A.*, B.n from tmp as A left join ymonth as B on A.year=B.year and A.month=B.month)
	group by household_code
	order by household_code;
quit;

data tmp1;
	set tmp1;
	by household_code;
	n = min_n - 1;
	if first.household_code then do; 	
		do while(n<max_n);
			n + 1; output;
		end;
	end;
run;
proc sort data=tmp1; by n; run;
data tmp1(drop=min_n max_n n); 
	merge tmp1 ymonth;
	by n;
run;
	
proc sql noprint;
	create table tmp2 as 
	select B.*, A.DOL, A.net_dol, A.dol_purchases
	from tmp as A full join tmp1 as B
	on A.household_code=B.household_code and A.year=B.year and A.month=B.month
	order by household_code, year, month;
quit;

proc sql; select nmiss(dol)/count(ymonth) from tmp2; quit;
proc sql noprint;
	create table tmp as 
	select household_code, nmiss(dol_purchases)/count(ymonth) as n_missmonth
	from tmp2
	group by household_code
	order by household_code;
quit;

proc univariate data = tmp; histogram n_missmonth; run;

*--------------------------------------------------------------------* 
* Delete the households who dont have continuous spending every month;
* This happens because Feburary has less than 30 days; 
proc sql noprint;
	create table mysel.keeppan as 
	select household_code 
	from tmp 
	where n_missmonth>0;
	
	create table tmp1 as
	select *
	from mainex.panelists 
	where household_code not in (select household_code from mysel.keeppan);
	
	create table tmp2 as
	select *
	from mainex.trips 
	where household_code not in (select household_code from mysel.keeppan);
	
	create table tmp3 as
	select *
	from mainex.month_basket
	where household_code not in (select household_code from mysel.keeppan);
quit;


data mainex.panelists; 	set tmp1 ; run;
data mainex.trips; 		set tmp2; run;
data mainex.month_basket; set tmp3; run;

*---------------------------------------*;
* Constrct data of household demographics;
proc sort data=mainex.panelists; by household_code panel_year; run;

* Adjust the income reporting; 
data tmp; 
	set mainex.panelists(keep = household_code panel_year household_income);
	income_year = panel_year - 2;
	drop panel_year; 
run;

* NOTE: we lose data from the last two years because of the lag income reporting;
proc freq data = mainex.all_panelists; table household_income; run;

proc sql noprint;
	create table tmp1 as 
	select A.*, B.household_income as income_real
	from mainex.panelists as A inner join tmp as B
	on A.household_code = B.household_code and A.panel_year = B.income_year;
quit;

data tmp;
	set tmp1(keep = household_code panel_year scantrack_market_descr magnet_flag magnet_const_flag 
			income_real household_size type_of_residence male_head_employment female_head_employment age_and_presence_of_children); 
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
	length income_group $5.;
	income_group = '';
	if income_real <= 15 and income_real ^='' then income_group = 'Qt1';
	if income_real>15 and income_real<=19 then income_group = 'Qt2';
	if income_real>19 and income_real<=23 then income_group = 'Qt3';
	if income_real>23 then income_group = 'Qt4';	
run;

* Add the middle value of income range; 
proc sql noprint;
	create table panelists as
	select A.*, B.mid_range as income_midvalue
	from tmp as A left join (select code_value, description, mid_range from mysel.my_codebook where variable='household_income') as B
	on A.income_real = B.code_value
	order by household_code;	
quit;

*---------------------------------------------------------------------*;
* Compute household expenditure share across retail format each month ;
proc sql noprint;
	create table tmp as
	select scantrack_market_descr, household_code, year, month, channel_type, sum(net_spent) as DOL, sum(dol_purchases) as dolpurch
	from mainex.trips
	group by scantrack_market_descr, household_code, year, month, channel_type;

	create table tmp1 as 
	select *, DOL/sum(DOL) as share
	from tmp 
	group by scantrack_market_descr, household_code,year, month
	order by scantrack_market_descr, household_code,year, month;
quit;

proc transpose data = tmp1 out = tmp_share1 prefix=DOL_;
	by scantrack_market_descr household_code year month;
	id channel_type;
	var DOL;
	* var share;
run;

proc transpose data = tmp1 out = tmp_share2 prefix=DOLP_;
	by scantrack_market_descr household_code year month;
	id channel_type;
	var dolpurch;
run;

data tmp_share; 
	merge tmp_share1(drop=_NAME_) tmp_share2(drop=_NAME_);
	by scantrack_market_descr household_code year month;
	array dolvar{*} DOL:;
	do i = 1 to dim(dolvar);
		if dolvar[i]=. then dolvar[i] = 0;
	end;
	drop i;
run;

*-------------------------------------------------------------------------*;
* Merge household demongraphics and price index to monthly expenditure data; 	
proc sql noprint;
	* Aggregate daily trip data to monthly expenditure data; 
	create table tmp as 
	select household_code, year, month, sum(total_spent) as DOL, sum(net_spent) as net_dol, sum(dol_purchases) as dol_purchases
	from mainex.trips 
	group by household_code, year, month;
	
	* Merge in expenditure share; 
	create table tmp1 as 
	select B.*, A.DOL, A.net_dol, A.dol_purchases
	from tmp as A inner join tmp_share as B
	on A.household_code=B.household_code and A.year=B.year and A.month=B.month;
	
	* Merge in price index; 
	create table tmp2 as
	select A.*, B.stone_price
	from tmp1 as A left join mainex.month_stone_price as B
	on A.scantrack_market_descr = B.scantrack_market_descr and A.year=B.year and A.month = B.month
	order by household_code, year, month; 
	
	* Merge in basket data;
	create table hh_month_exp as 
	select A.*, B.num_day, B.num_trip, B.coupon_dol,
				B.dol_magnet, B.magnet_quant, 
				B.num_food_module, B.num_noned_module, B.food_quant, B.nonedible_quant, 
				B.num_storable_module, B.num_nonstr_module, B.storable_quant, B.nonstr_quant
	from tmp2 as A left join mainex.month_basket as B
	on A.household_code=B.household_code and A.year=B.year and A.month=B.month;
	
	* Merge in demographics; 
	* NOTE: WE LOSE ABOUT HALF DATA HERE WHEN MERGING WITH DEMOGRAPHICS; 
	create table hh_month_exp_merge(drop = household_code_old panel_year_old) as 
	select A.*, B.*
	from hh_month_exp as A inner join 
	panelists(rename=(household_code = household_code_old panel_year=panel_year_old) drop= scantrack_market_descr) as B
	on A.household_code = B.household_code_old and A.year = B.panel_year_old
	order by household_code;
quit;

* Adjust all dollar measure for CPI; 
proc sql noprint;
	select Annual into: price_base from cpi where year = 2004; 
	
	create table tmp as 
	select A.*, B.Annual/&price_base as cpi
	from hh_month_exp_merge as A left join cpi as B
	on A.year = B.year
	order by household_code, year, month;
quit;
proc contents data=hh_month_exp_merge; run;

* Adjust CPI;
data hh_month_exp_merge;
	set tmp;
	array dolvar{*} dol: coupon_dol stone_price net_dol;
	do i=1 to dim(dolvar);
		dolvar[i] = dolvar[i]/cpi;
	end;
	recession = 0;
	if year >= 2008 then recession = 1;
run;	

* Update the data in the library of mainex;
data mylib.panelists; set panelists; run;
data mylib.hh_month_exp_merge;	set hh_month_exp_merge; run;

***********************************************************************************;
* Construct the expenditure data from all the samples (with no household selection);
***********************************************************************************;
*---------------------------------------*;
* Constrct data of household demographics;
proc sort data=mainex.all_panelists; by household_code panel_year; run;

* Adjust the income reporting; 
data tmp; 
	set mainex.all_panelists(keep = household_code panel_year household_income);
	income_year = panel_year - 2;
	drop panel_year; 
run;

* NOTE: we keep the data even without income reports; 
proc sql noprint;
	create table tmp1 as 
	select A.*, B.household_income as income_real
	from mainex.all_panelists as A left join tmp as B
	on A.household_code = B.household_code and A.panel_year = B.income_year;
quit;

data tmp;
	set tmp1(keep = household_code panel_year scantrack_market_descr magnet_flag 
			income_real household_size type_of_residence male_head_employment female_head_employment age_and_presence_of_children); 
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
	length income_group $5.;
	income_group = '';
	if income_real <= 15 and income_real ^='' then income_group = 'Qt1';
	if income_real>15 and income_real<=19 then income_group = 'Qt2';
	if income_real>19 and income_real<=23 then income_group = 'Qt3';
	if income_real>23 then income_group = 'Qt4';
run;

* Add the middle value of income range; 
proc sql noprint;
	create table panelists as
	select A.*, B.mid_range as income_midvalue
	from tmp as A left join (select code_value, description, mid_range from mysel.my_codebook where variable='household_income') as B
	on A.income_real = B.code_value
	order by household_code;	
quit;

*---------------------------------------------------------------------*;
* Compute household expenditure share across retail format each month ;
proc sql noprint;
	create table tmp as
	select scantrack_market_descr, household_code, year, month, channel_type, sum(net_spent) as DOL, sum(dol_purchases) as dolpurch
	from mainex.all_trips
	group by scantrack_market_descr, household_code, year, month, channel_type;

	create table tmp1 as 
	select *, DOL/sum(DOL) as share
	from tmp 
	group by scantrack_market_descr, household_code,year, month
	order by scantrack_market_descr, household_code,year, month;
quit;

proc transpose data = tmp1 out = tmp_share1 prefix=DOL_;
	by scantrack_market_descr household_code year month;
	id channel_type;
	var DOL;
	* var share;
run;

proc transpose data = tmp1 out = tmp_share2 prefix=DOLP_;
	by scantrack_market_descr household_code year month;
	id channel_type;
	var dolpurch;
run;

data tmp_share; 
	merge tmp_share1(drop=_NAME_) tmp_share2(drop=_NAME_);
	by  scantrack_market_descr household_code year month;
	array dolvar{*} DOL:;
	do i = 1 to dim(dolvar);
		if dolvar[i]=. then dolvar[i] = 0;
	end;
	drop i;
run;

*-------------------------------------------------------------------------*;
* Merge household demongraphics and price index to monthly expenditure data; 	
proc sql noprint;
	* Aggregate daily trip data to monthly expenditure data; 
	create table tmp as 
	select household_code, year, month, sum(total_spent) as DOL, sum(net_spent) as net_dol, sum(dol_purchases) as dol_purchases
	from mainex.all_trips 
	group by household_code, year, month;
	
	* Merge in expenditure share; 
	create table tmp1 as 
	select B.*, A.DOL, A.net_dol, A.dol_purchases
	from tmp as A inner join tmp_share as B
	on A.household_code=B.household_code and A.year=B.year and A.month=B.month;
	
	* Merge in price index; 
	create table tmp2 as
	select A.*, B.stone_price
	from tmp1 as A left join mainex.month_stone_price as B
	on A.scantrack_market_descr = B.scantrack_market_descr and A.year=B.year and A.month = B.month
	order by household_code, year, month; 
	
	* Merge in basket data;
	create table all_hh_month_exp as 
	select A.*, B.num_day, B.num_trip, B.coupon_dol, B.dol_magnet, B.magnet_quant,
				B.num_food_module, B.num_noned_module, B.food_quant, B.nonedible_quant, 
				B.num_storable_module, B.num_nonstr_module, B.storable_quant, B.nonstr_quant
	from tmp2 as A left join mainex.all_month_basket as B
	on A.household_code=B.household_code and A.year=B.year and A.month=B.month;
	
	* Merge in demographics, we still have all the expenditure data; 
	create table all_hh_month_exp_merge(drop = household_code_old panel_year_old) as 
	select A.*, B.*
	from all_hh_month_exp as A inner join 
	panelists(rename=(household_code = household_code_old panel_year=panel_year_old) drop= scantrack_market_descr) as B
	on A.household_code = B.household_code_old and A.year = B.panel_year_old
	order by household_code, year, month;
quit;

* Adjust all dollar measure for CPI; 
proc sql noprint;
	select Annual into: price_base from cpi where year = 2004; 
	
	create table tmp as 
	select A.*, B.Annual/&price_base as cpi
	from all_hh_month_exp_merge as A left join cpi as B
	on A.year = B.year
	order by household_code, year, month;
quit;
proc contents data=all_hh_month_exp_merge; run;

* Adjust CPI;
data all_hh_month_exp_merge;
	set tmp;
	array dolvar{*} dol: coupon_dol stone_price net_dol;
	do i=1 to dim(dolvar);
		dolvar[i] = dolvar[i]/cpi;
	end;
	recession = 0;
	if year >= 2008 then recession = 1;
run;	

* Update the data in the library of mainex;
data mylib.all_panelists; set panelists; run;
data mylib.all_hh_month_exp_merge; set all_hh_month_exp_merge;	
proc sort; by scantrack_market_descr household_code year month;
run;

proc datasets noprint; delete tmp1 tmp2 tmp3 tmp_share1 tmp_share2 tmp_share; run;

*****************************************;
* Export the data to csv and			*;
* Kill the tables in physical libraries *;
*****************************************;
* A subsample of households for testing; 
data tmp; set mylib.panelists; where scantrack_market_descr="Chicago"; run;
PROC EXPORT DATA=tmp
			OUTFILE= "\\tsclient\Resear1\Store switching\processed data\2_panelists_sub.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;

data tmp; set mylib.hh_month_exp_merge; where scantrack_market_descr="Chicago"; run;
PROC EXPORT DATA=tmp
			OUTFILE= "\\tsclient\Resear1\Store switching\processed data\2_hh_month_exp_merge_sub.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;

*---------------------------------------*;
PROC EXPORT DATA=mylib.panelists
			OUTFILE= "\\tsclient\Resear1\Store switching\processed data\2_panelists.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;

PROC EXPORT DATA=mylib.hh_month_exp_merge
			OUTFILE= "\\tsclient\Resear1\Store switching\processed data\2_hh_month_exp_merge.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;

/*PROC EXPORT DATA=mylib.all_panelists
			OUTFILE= "\\tsclient\Resear1\Store switching\processed data\2_all_panelists.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;

PROC EXPORT DATA=mylib.all_hh_month_exp_merge
			OUTFILE= "\\tsclient\Resear1\Store switching\processed data\2_all_hh_month_exp_merge.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;
*/

/*
proc datasets library=mysel kill noprint;
proc datasets library=mainex kill noprint;
run;
*/



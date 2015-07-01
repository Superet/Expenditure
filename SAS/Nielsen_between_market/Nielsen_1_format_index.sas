libname mylib "E:\Users\ccv103\Documents\Research\Store switching\SAS_temp\Nielsen";
%let dirname= E:\Users\ccv103\Documents\Research\Nielsen\nielsen_extracts\HMS;

*************;
* Read data *;
*************;
* Read in selected markets;
PROC IMPORT OUT= my_market
            DATAFILE= "\\tsclient\Resear1\Store switching\processed data\selected market.csv" 
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
PROC IMPORT OUT= CPI
            DATAFILE= "\\tsclient\Resear1\Store switching\processed data\CPI.csv" 
            DBMS=csv REPLACE;
     		GETNAMES=YES; 
RUN;

* Define the retailer scope; 
data my_channel;
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

proc sql noprint;
	create table retailers as
	select *
	from orig_retailers
	where channel_type in (select channel_type from my_channel where channel_group = "Grocery");
quit;

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

data my_delete_upc; set tmp1 tmp2; run;	

* Read UPC-level purchase data;
%macro read_groc_fun;
	%do i=2004 %to 2012;
		* Read panelists data for geographic market information; 
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
			drop wic_indicator_current wic_indicator_ever_notcurrent Member:;
		run;
		
		* Restirct to the focal markets;
		proc sql noprint;
			create table tmp as 
			select A.*, B.market, B.mkt_type
			from panelists&i as A inner join my_market as B
			on A.scantrack_market_descr = B.scantrack_market_descr;
		quit;
		data panelists&i; set tmp; run;
		
		proc sql noprint; 
			* Restrict households to be in the focal market;		
			* Merge in panelists geographic data to trip data;
			create table tmp as 
			select A.*, year(A.purchase_date) as year, B.scantrack_market_descr, B.market, B.mkt_type
			from trips&i as A inner join panelists&i as B
			on A.household_code = B.household_code;	
					
			* Restrict the trip data to the selected retailers;
			create table tmp1 as
			select A.*, B.channel_type
			from tmp as A inner join retailers as B
			on A.retailer_code = B.retailer_code;			
			
			* Subset purchases data within restrictions;
			create table tmp2 as 
			select A.*, cats(upc, "-", upc_ver_uc) as upcv, (total_price_paid-coupon_value) as dol_paid, total_price_paid/quantity as price,
					B.market, B.mkt_type, B.scantrack_market_descr, B.channel_type, B.retailer_code, B.panel_year, 
					year(B.purchase_date) as year, month(B.purchase_date) as month
			from purchases&i as A inner join tmp1 as B
			on A.trip_code_uc = B.trip_code_uc;
			
			* Delete the undesirable UPC;
			create table tmp3 as 
			select *
			from tmp2 
			where upcv not in (select upcv from my_delete_upc);
			
			* A subset of purchases of deleted products, then aggregated by trip; 
			create table tmp4 as 
			select trip_code_uc, sum(total_price_paid-coupon_value) as delete_dol
			from (select * from tmp2 where upcv in(select upcv from my_delete_upc))
			group by trip_code_uc;
			
			* Merge in the delete dollar for each trip;
			create table trips&i as 
			select A.*, B.delete_dol
			from tmp1 as A left join tmp4 as B
			on A.trip_code_uc = B.trip_code_uc;
			
			* Aggregate individual transaction to annual data at upc level;
			create table annual_purchase&i as
			select market, year, channel_type, retailer_code, upcv, 
					mean(upc) as upc, mean(upc_ver_uc) as upc_ver_uc, 
					sum(total_price_paid) as dol, sum(quantity) as quantity,
					sum(dol_paid) as dol_paid, sum(total_price_paid)/sum(quantity) as price_tag, sum(dol_paid)/sum(quantity) as price_paid
			from (select * from tmp3 where mkt_type ^= "price_mkt")
			group by market, channel_type, retailer_code, upcv, year;
			
			* Aggretate individual transaction to monthly data at upc level;
			create table month_purchase&i as
			select market, channel_type, retailer_code, upcv, year, month,
					mean(upc) as upc, mean(upc_ver_uc) as upc_ver_uc, 
					sum(total_price_paid) as dol, sum(quantity) as quantity,
					sum(dol_paid) as dol_paid, sum(total_price_paid)/sum(quantity) as price_tag, sum(dol_paid)/sum(quantity) as price_paid
			from (select * from tmp3 where mkt_type ^= "price_mkt" OR 
					(mkt_type="price_mkt" and channel_type in ("Convenience Store","Health Food Store","Dollar Store" )))
/*			from tmp3*/
			group by market, channel_type, retailer_code, upcv, year, month;
		quit;
		
		data trips&i; 
			set trips&i; 
			if delete_dol=. then delete_dol=0;
			net_spent = total_spent - delete_dol;
			if net_spent<0 then net_spent = 0;
		run;		
		proc datasets noprint; delete purchases&i tmp tmp1 tmp2 tmp3 tmp4; run;
	%end;
	
	* Combein the panelists data together; 
	data panelists;
		set panelists2004 - panelists2012;
	run;
	
	* Combine the annual purchase data together;
	data annual_purchase;
		set annual_purchase2004 - annual_purchase2012;
	run;
	
	* Combine trips data together;
	data trips;
		set trips2004 - trips2012;
	run;
	
	* Combine monthly transaction;
	data month_transaction;
		set month_purchase2004 - month_purchase2012;
	run;
	
	proc datasets noprint; 
		delete panelists2004 - panelists2012 annual_purchase2004 - annual_purchase2012 trips2004 - trips2012 month_purchase2004 - month_purchase2012; 
	run;
%mend read_groc_fun;

%read_groc_fun;

************************************************;
* Compute product size index within categories *; 
************************************************;
* ------------------------------------------------------------------------------------------ *;
* In the original data, the unit measurements (size1_units) are not unique for a single module;
* Assign a unit measurement for each product module; 
proc sort data = orig_products; by department_descr product_group_descr product_module_descr;
proc freq data = orig_products noprint; 
	by department_descr product_group_descr product_module_descr;
	tables size1_units/ out=tmp;
run;

proc sort data=tmp; by department_descr product_group_descr product_module_descr COUNT;
data module_measure;
	set tmp;
	by department_descr product_group_descr product_module_descr; 
	if last.product_module_descr then output;
run;

* ------------------------------------------------------------------------------------------ *;
* Imputing size for the UPCs that have different unit measurement from the majority UPCs in the category; 
* Use medium unit measurement to replace 1 CT;
* First merge in the main unit measure for each module; 
proc sql noprint;
	create table tmp as 
	select A.*, A.size1_amount*A.multi as size, B.size1_units as module_units
	from orig_products as A left join module_measure as B
	on A.product_module_descr = B.product_module_descr
	order by department_descr, product_group_descr, product_module_descr;
quit;

* Check what percent of products have different unit measure;
proc sql; select sum(size1_units ^= module_units)/count(size1_units) from tmp; quit;

* Compute the median size of each category;
proc univariate data=tmp(where = (size1_units=module_units)) noprint;
	by department_descr product_group_descr product_module_descr;
	var size;
	output out = tmp1 MEDIAN = main_med_size;
run;

* Merge the product data with median size for other unit measurements;
proc sql noprint;
	create table tmp2 as 
	select A.*, B.main_med_size
	from tmp as A left join tmp1 as B
	on A.product_module_descr = B.product_module_descr
	order by department_descr, product_group_descr, product_module_descr, upc;
quit;

* If the unit measurement of a product is not the main measurement of the module, then we assign the medium size with known measurement; 
data tmp2;
	set tmp2;
	if size1_units ^= module_units then size = main_med_size;
	upcv = cats(upc, "-", upc_ver_uc);
run;

proc rank data=tmp2 out=products ties=mean PERCENT;
	by department_descr product_group_descr product_module_descr;
	var size;
	ranks size_percentile;
run;

* Drop some UPC from magnet data;
proc sql noprint;
	create table tmp_products as 
	select *
	from products 
	where upcv not in (select upcv from my_delete_upc)
	order by department_descr, product_group_descr, product_module_descr;
quit;

*---------------------------------------------------------------------------------*;
* Compute product size index by normalizing the size with median size in a category;
data products;
	merge tmp_products tmp1;
	by department_descr product_group_descr product_module_descr;
	size_index = size/main_med_size; 
	PRVT_flag = 0; 
	if brand_descr="CTL BR" then PRVT_flag = 1;
	label	module_units = "Unit measure of category"
			main_med_size = "The median size of category measured at unit measure"
			upcv	= "UPC-version"
			size_percentile = "Size percentile within category"
			size_index = "Size normalized by median size within category";
run;

data mylib.products; set products; run;
proc datasets noprint; delete tmp tmp1 tmp2 module_measure tmp_products; run;

/*
*********************************************************************;
* Run the following codes if we want to delete some product modules *;
*********************************************************************;

* Compute the size index of each upc within a category;
* NOTE: MULTI-PACK IS INCORPORATED HERE;
proc sql noprint; 
	create table products as 
	select cats(A.upc, "-", A.upc_ver_uc) as upcv, A.*, B.med_size, 
			100 + 100*(size - B.med_size)/B.med_size as size_index
	from tmp2 as A left join tmp as B
	on A.product_module_descr = B.product_module_descr
	order by department_descr, product_group_descr, product_module_descr;
quit;

* Drop the modules that have missing values, and drop magnet data;
data products; 
	set products;
	WHERE product_module_descr not in (' ', 'MAGNET DATA') and product_group_descr ^= ' ' and department_descr^=' ';
run;

* Check if any extreme value of size_index;
proc univariate data=products;
	var size_index;
run;

* Find the modules that have the most extreme values;
proc freq data=products(where = (size_index>5000));
tables product_module_descr / out = tmp;
run;
proc sort data=tmp; by descending count; run;

* Drop the three modules that have too many extreme size index;
proc sql noprint;
	create table tmp1 as 
	select *
	from products
	where product_module_descr not in (select product_module_descr from tmp(obs=3));
quit;
data products; set tmp1; run;

proc datasets noprint; delete tmp tmp1 tmp2; run;
*/

***********************************;
* Wallet share of product modules *;
***********************************;
proc sql noprint;
	* Apend product module to purchase data;
	create table tmp as 
	select A.*, B.product_module_descr
	from annual_purchase as A left join products as B
	on A.upcv = B.upcv;
	
	* Sum up dollars spent on each module in every market;
	create table tmp1 as 
	select market, year, product_module_descr, sum(dol) as dol, sum(dol_paid) as dol_paid
	from tmp
	group by market, year, product_module_descr;
	
	* Compute wallet share of each module in every market;
	create table module_wallet as
	select market, year, product_module_descr, dol/sum(dol) as share, dol_paid/sum(dol_paid) as share_paid
	from tmp1 
	group by market, year
	order by market, year, product_module_descr;
quit;

data mylib.module_wallet; set module_wallet; run;
proc datasets noprint; delete tmp tmp1; run;

********************************************************************;
* Channel-leve size index averaged from retailer-level calculation *;
********************************************************************;
* Compute the product size at each retailer in each market ;
proc sql noprint;
	* Use average transaction volume to compute weights associated with each upc;
	create table tmp as
	select  market, year, channel_type, retailer_code, upcv, mean(upc) as upc, 
			mean(upc_ver_uc) as upc_ver_uc, mean(quantity) as quantity, mean(dol_paid) as revenue, 
			mean(price_paid) as price_paid, min(year) as start_year, max(year) as end_year
	from annual_purchase
	group by market, year, channel_type, retailer_code, upcv;

	* Append the upc-related information from products data;
	create table retailer_upc as 
	select A.*, B.department_descr, B.product_group_descr, B.product_module_descr, B.brand_descr, B.size_index,
			B.size1_amount, B.size1_units, B.multi, B.PRVT_flag, quantity*size1_amount*multi as vol
	from tmp as A inner join products as B
	on A.upcv = B.upcv;
	
	* Aggregate to retailer-module level weighted by UPC volumn at each retailer;
	create table retailer_module as 
	select market, year, channel_type, retailer_code, department_descr, product_group_descr, product_module_descr,
			sum(size_index*vol)/sum(vol) as size_index_volwt, 
			count(distinct brand_descr) as num_brands, count(distinct upcv) as num_upc, 
			mean(PRVT_flag) as PRVT, sum(PRVT_flag) as num_prvt,
			sum(quantity) as module_units, sum(vol) as module_vol, sum(revenue) as revenue
	from retailer_upc
	group by market, year, channel_type, retailer_code, department_descr, product_group_descr, product_module_descr;
	
	* Aggregate retailer level index using wallet share as module weight;
	create table retailer_size as 
	select market, year, channel_type, retailer_code,
		sum(size_index_volwt*share_paid)/sum(share_paid) as size_index,
		sum(num_brands) as num_brands, sum(num_brands*share_paid)/sum(share_paid) as avg_brands_unitwt,
		sum(num_upc) as num_upc, sum(num_upc*share_paid)/sum(share_paid) as avg_upc_unitwt,
		sum(PRVT*share_paid)/sum(share_paid) as avg_prvt,
		sum(num_prvt)/sum(num_brands) as overall_prvt,
		count(distinct product_module_descr) as num_module
	from (select A.*, B.share_paid from retailer_module as A left join module_wallet as B 
		  on A.market=B.market and A.year = B.year 
			and A.product_module_descr = B.product_module_descr)
	group by market, year, channel_type, retailer_code
	order by market, year, channel_type, retailer_code;
	
/*	* Aggregate to retailer level weighted by product module scan units;
	create table retailer_size as
	select market, channel_type, retailer_code,
			sum(size_index_volwt*module_units)/sum(module_units) as size_index,
			sum(num_brands) as num_brands, sum(num_brands*module_units)/sum(module_units) as avg_brands_unitwt,
			sum(num_upc) as num_upc, sum(num_upc*module_units)/sum(module_units) as avg_upc_unitwt,
			sum(PRVT*module_units)/sum(module_units) as avg_prvt,
			count(distinct product_module_descr) as num_module
	from retailer_module 
	group by market, channel_type, retailer_code;*/
quit;

* Compute the format level index from retailer-level index weighted by transaction number;
proc sql noprint;	
	* Annual number of trips and revenue for each retailer;
	create table retailer_transaction as
	select market, year, channel_type, retailer_code, 
			count(trip_code_uc) as num_trips, sum(net_spent) as revenue
	from trips
	group by market, year, channel_type, retailer_code;

	* Append retailer transaction information to the size index at retailer level;
	create table tmp as
	select A.*, B.num_trips, B.revenue
	from retailer_size as A left join retailer_transaction as B
	on A.market=B.market and A.year=B.year and A.retailer_code = B.retailer_code;
	
	* Channel aggregation 2: Aggregate to market-channel level from retailer data weighted by store revenue;
	create table channel_size as 
	select market, year, channel_type, sum(size_index*revenue)/sum(revenue) as size_index, 
			sum(num_module*revenue)/sum(revenue) as num_module,
			sum(num_brands*revenue)/sum(revenue) as num_brands,
			sum(avg_brands_unitwt*revenue)/sum(revenue) as avg_brands_per_mod,
			sum(num_upc*revenue)/sum(revenue) as num_upc,
			sum(avg_upc_unitwt*revenue)/sum(revenue) as avg_upc_per_mod,
			sum(avg_prvt*revenue)/sum(revenue) as avg_prvt_per_mod, 
			sum(overall_prvt*revenue)/sum(revenue) as overall_prvt
	from tmp
	group by market, year, channel_type
	order by market, year, channel_type;
quit;

data mylib.format_year_attr; 
	set channel_size; 
	label	size_index		= "Size index"
			num_module		= "Number of modules"
			num_brands		= "Total number of brands"
			avg_brands_per_mod = "Number of brands per module"
			num_upc			= "Total number of UPC"
			avg_upc_per_mod	= "Number of UPC per module"
			avg_prvt_per_mod= "Porportion of private labels per module"
			overall_prvt		= "Proportion of private labels among all brands";
run;
proc datasets noprint; delete tmp retailer_upc retailer_module retailer_size retailer_transaction; run;

/*
*************************************************************************;	
* Directly compute channel leve size index without aggregating retailers*;
*************************************************************************;
proc sql noprint;
	* Average annual quantity of UPC in a channel type in a market;
	create table tmp as 
	select market, channel_type, upcv, mean(upc) as upc, mean(upc_ver_uc) as upc_ver_uc,
			mean(quantity) as quantity, mean(dol) as revenue
	from annual_purchase
	group by market, channel_type, upcv;
	
	* Append with UPC-related information from products;
	create table tmp1 as 
	select A.*, B.department_descr, B.product_group_descr, B.product_module_descr, B.brand_descr, B.size_index,
			B.size1_amount, B.size1_units, B.multi, quantity*size1_amount*multi as vol
	from tmp as A inner join products as B
	on A.upcv = B.upcv;
	
	* Compute channel-module index weighted by module units;
	create table channel_mod1 as 
	select market, channel_type, department_descr, product_group_descr, product_module_descr, 
			sum(size_index*vol)/sum(vol) as size_index_volwt, count(distinct brand_descr) as num_brands, 
			count(distinct upcv) as num_upc, sum(quantity) as module_units, sum(vol) as module_vol, sum(revenue) as module_revenue
	from tmp1 
	group by market, channel_type, department_descr, product_group_descr, product_module_descr;
	
	* Aggregate channel-level index using wallet share as module weight;
	create table channel_size1 as 
	select market, channel_type, sum(size_index_volwt*share_paid)/sum(share_paid) as size_index,
		sum(num_brands) as num_brands, sum(num_brands*share_paid)/sum(share_paid) as avg_brands_unitwt,
		sum(num_upc) as num_upc, sum(num_upc*share_paid)/sum(share_paid) as avg_upc_unitwt,
		count(distinct product_module_descr) as num_module
	from (select A.*, B.share_paid from channel_mod1 as A left join module_wallet as B 
		  on A.market=B.market and A.product_module_descr = B.product_module_descr)
	group by market, channel_type
	order by market, channel_type;
quit; 

* Plot distribution of product sizes within a module across channel types;
%let sel_mkt = "Chicago";
%let sel_mod = "BEER";

proc sql noprint;
	create table tmp_plot as
	select *
	from tmp1 
	where market=&sel_mkt and product_module_descr = &sel_mod
	order by channel_type;
quit;

proc sgpanel data=tmp_plot;
	panelby channel_type;
	histogram size_index;
	density size_index / type = kernel;
run;

*/

*********************************;
* Construct monthly price index *;
*********************************;
* Compute the monthly price index for each retailer at each market;
proc sql noprint;
	* Append the upc-related information from products data;
	create table retailer_upc as 
	select A.*, B.department_descr, B.product_group_descr, B.product_module_descr, B.brand_descr, B.size_index,
			B.size1_amount, B.size1_units, B.multi, 
			quantity*size_index as vol, price_tag/size_index as price_tag_norm, price_paid/size_index as price_paid_norm
	from month_transaction as A inner join products as B
	on A.upcv = B.upcv;
	
	* Aggregate to retailer-module level weighted by UPC volumn at each retailer;
	create table retailer_module as 
	select market, year, month, channel_type, retailer_code, department_descr, product_group_descr, product_module_descr,
			sum(price_tag*vol)/sum(vol) as price_tag_volwt, sum(price_paid*vol)/sum(vol) as price_paid_volwt,
			sum(price_tag_norm*vol)/sum(vol) as price_tag_volnorm, sum(price_paid_norm*vol)/sum(vol) as price_paid_volnorm,			
			sum(quantity) as module_units, sum(vol) as module_vol, sum(dol_paid) as revenue, sum(dol) as dol
	from retailer_upc
	group by market, year, month, channel_type, retailer_code, department_descr, product_group_descr, product_module_descr;
	
	* Aggregate retailer level index using wallet share as module weight;
	* NOTE: WE USE ANNUAL WALLET SHARE TO WEIGHT MONTHLY PRICES ACROSS CATEGORIES; 
	create table retailer_index as 
	select market, year, month, channel_type, retailer_code,
		sum(price_tag_volwt*share_paid)/sum(share_paid) as price_tag_index, 
		sum(price_paid_volwt*share_paid)/sum(share_paid) as price_paid_index, 
		sum(price_tag_volnorm*share_paid)/sum(share_paid) as price_tag_normindex, 
		sum(price_paid_volnorm*share_paid)/sum(share_paid) as price_paid_normindex
	from (select A.*, B.share_paid from retailer_module as A left join module_wallet as B 
		  on A.market=B.market and A.year = B.year and 
		  A.product_module_descr = B.product_module_descr)
	group by market, year, month, channel_type, retailer_code
	order by market, year, month, channel_type, retailer_code;
quit;

* Compute the format level index from retailer-level index weighted by transaction number;
proc sql noprint;	
	* Annual number of trips and revenue for each retailer;
	create table retailer_transaction as
	select market, year, month, channel_type, retailer_code, 
			count(trip_code_uc) as num_trips, sum(net_spent) as revenue
	from (select *, month(purchase_date) as month from trips )
	group by market, year, month, channel_type, retailer_code;

	* Append retailer transaction information to the price index at retailer level;
	create table tmp as
	select A.*, B.num_trips, B.revenue
	from retailer_index as A left join retailer_transaction as B
	on A.market=B.market and A.year=B.year 
		and A.month=B.month and A.retailer_code = B.retailer_code;
	
	* Channel aggregation: Aggregate to market-channel level from retailer data weighted by store revenue;
	create table channel_index as 
	select market, year, month, channel_type, 
			sum(price_tag_index*revenue)/sum(revenue) as price_tag_index,
			sum(price_paid_index*revenue)/sum(revenue) as price_paid_index,
			sum(price_tag_normindex*revenue)/sum(revenue) as price_tag_norm_index,
			sum(price_paid_normindex*revenue)/sum(revenue) as price_paid_norm_index
	from tmp
	group by market, year, month, channel_type
	order by market, year, month, channel_type;
quit;

* Merge in the CPI to have price at common dollars; 
proc sql noprint;
	select annual into: cpi2004 from cpi where year=2004 ;
	
	create table tmp as 
	select A.*, (B.annual/&cpi2004) as cpi
	from channel_index as A left join CPI as B
	on A.year = B.year
	order by market, year, month, channel_type;
quit;

data mylib.format_month_price; 
	set tmp; 
	price_tag_2004		= price_tag_index/cpi;
	price_paid_2004		= price_paid_index/cpi;
	price_tag_norm_2004	= price_tag_norm_index/cpi;
	price_paid_norm_2004= price_paid_norm_index/cpi;
	label	price_tag_index			= "Average price tag per item"
			price_paid_index		= "Average price paid per item"
			price_tag_norm_index	= "Average price tag per size index"
			price_paid_norm_index	= "Average price paid per size index"
			price_tag_2004			= "Average price tag per item at 2004 dollar"
			price_paid_2004			= "Average price paid per item at 2004 dollar"
			price_tag_norm_2004		= "Average price tag per size index at 2004 dollar"
			price_paid_norm_2004	= "Average price paid per size index at 2004 dollar";
run;
	
proc datasets noprint; delete tmp retailer_upc retailer_module retailer_index retailer_transaction channel_index; run;

**************************;
* Export the data to csv *;
**************************;
PROC EXPORT DATA=mylib.format_year_attr
			OUTFILE= "\\tsclient\Resear1\Store switching\processed data\1_format_year_attr.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;

PROC EXPORT DATA=mylib.format_month_price
			OUTFILE= "\\tsclient\Resear1\Store switching\processed data\1_format_month_price.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;




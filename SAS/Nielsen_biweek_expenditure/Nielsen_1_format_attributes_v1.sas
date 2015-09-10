LIBNAME mysqllib OLEDB
OLEDB_SERVICES=NO
Datasource="kdc01\kdwh02"
PROVIDER=SQLOLEDB.1
Properties=('Initial Catalog'=USRDB_ccv103
                  'Integrated Security'=SSPI)
SCHEMA=DBO
BULKLOAD=YES
bl_options='ROWS_PER_BATCH = 500000';

libname mylib "E:\Users\ccv103\Documents\Research\Store switching\SAS_temp\Nielsen";
libname mysel "E:\Users\ccv103\Documents\Research\Store switching\SAS_temp\Nielsen\Work";
libname mainex "E:\Users\ccv103\Documents\Research\Store switching\SAS_temp\Nielsen\Work\main";

%let dirname= E:\Users\ccv103\Documents\Research\Nielsen\nielsen_extracts\HMS;
options compress=yes reuse=yes; 

*************;
* Read data *;
*************;
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

* Read the conversion of MAGNET DATA;
PROC IMPORT OUT= magnet_convt
            DATAFILE= "\\tsclient\Resear1\Store switching\processed data\convert_magnet_product.csv" 
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
proc datasets noprint; delete tmp1 tmp2; run;

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
			drop dma_descr fips:; 
			drop kitchen_appliances tv_items household_internet_connection wic_indicator_current wic_indicator_ever_notcurrent Member:;
		run;
		
		proc sql noprint; 
			* Merge in panelists geographic data to trip data;
			create table tmp as 
			select A.*, year(A.purchase_date) as year, B.scantrack_market_descr
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
					B.scantrack_market_descr, B.channel_type, B.retailer_code, 
					year(B.purchase_date) as year, ceil((purchase_date - MDY(12,31,2003))/14) as biweek
			from purchases&i as A inner join tmp1 as B
			on A.trip_code_uc = B.trip_code_uc;
			
			* Delete the undesirable UPC;
			create table purchases&i as 
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
		quit;
		
		data trips&i; 
			set trips&i; 
			if delete_dol=. then delete_dol=0;
			net_spent = total_spent - delete_dol;
			if net_spent<0 then net_spent = 0;
		run;		
		
		%if &i=2004 %then %do; 
			data panelists;	set panelists&i; run; 
			data purchases; 	set purchases&i; run;
			data trips; 	set trips&i; run;
		%end;
		%else %do;
			data panelists;	set panelists panelists&i; run;
			data purchases; 	set purchases purchases&i; run; 
			data trips; 	set trips trips&i; run;
		%end; 
		proc datasets noprint; delete purchases&i tmp tmp1 tmp2 tmp3 tmp4 trips&i panelists&i; run;
	%end;
%mend read_groc_fun;
%read_groc_fun;

* Save data; 
data mainex.panelists;	set panelists; run;
data mainex.trips; 		set trips; run;
data mainex.purchases; 	set purchases; run;

data mysqllib.trips (BULKLOAD=yes); 	
set mainex.trips; 
run;

data mysqllib.purchases (BULKLOAD=yes); 
set mainex.purchases; 
run; 

* Aggregate annual purchases; 
proc sql noprint;
	create table annual_purchase as 
	select scantrack_market_descr, year, channel_type, retailer_code, upcv, 
			sum(total_price_paid) as dol, sum(quantity) as quantity,
			sum(dol_paid) as dol_paid, sum(total_price_paid)/sum(quantity) as price_tag, sum(dol_paid)/sum(quantity) as price_paid
	from purchases
	group by scantrack_market_descr, channel_type, retailer_code, upcv, year;
quit; 

* Aggregate biweek purchases; 
proc sql noprint;	
	create table biweek_purchase as 
	select scantrack_market_descr, channel_type, retailer_code, upcv, biweek,
			min(year) as year,
			mean(upc) as upc, mean(upc_ver_uc) as upc_ver_uc, 
			sum(total_price_paid) as dol, sum(quantity) as quantity,
			sum(dol_paid) as dol_paid, sum(total_price_paid)/sum(quantity) as price_tag, sum(dol_paid)/sum(quantity) as price_paid
	from purchases
	group by scantrack_market_descr, channel_type, retailer_code, upcv, biweek;
quit; 

***********************;
* Set basket baseline *; 
***********************;
* Examine what belongs to magnet modules; 
* Change product modules in MAGNET DATA to other DEPARMENTS;
proc sql noprint;
	create table tmp1 as 
	select *
	from orig_products
	where department_descr = "MAGNET DATA" and upc not in (select upc from my_delete_upc);
quit;

* Note that starting in 2011, Magnet data is now also found in Product Modules 0445-0468;
* but these Product Modules are assigned to many different Product; 
* Just see what those modules are; 
data tmp2; 
	set orig_products;
	where product_module_code between 445 and 468;
run;

* Change some reference card products from magnet data to specific departments; 
proc sql noprint;
	create table my_products as 
	select A.*, B.module_new, B.group_new, B.department_new
	from orig_products as A left join 
		magnet_convt(rename=(product_module_descr=module_new product_group_descr=group_new department_descr=department_new)) as B
	on A.upc = B.upc
	order by department_new;
quit;

data my_products;
	set my_products;
	my_magnet = 0; 
	if department_descr = "MAGNET DATA" then my_magnet = 1;
	if product_module_code >= 445 & product_module_code <= 468 then my_magnet = 1; 
/*	if module_new ^= '' then do;
		product_module_descr=module_new;
		product_group_descr=group_new;
		department_descr=department_new;
	end;*/
	drop module_new group_new department_new;
run;
proc datasets noprint; delete tmp1 tmp2; run;

* ------------------------------------------------------------------------------------------ *;
* In the original data, the unit measurements (size1_units) are not unique for a single module;
* Assign a unit measurement for each product module; 
proc sort data = my_products; by department_descr product_group_descr product_module_descr;
proc freq data = my_products noprint; 
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
	from my_products as A left join module_measure as B
	on A.product_module_descr = B.product_module_descr
	order by department_descr, product_group_descr, product_module_descr;
quit;

* Check what percent of products have different unit measure;
proc sql; select sum(size1_units ^= module_units)/count(size1_units) from tmp; quit;

* Compute the median size of each category;
proc univariate data=tmp(where = (size1_units=module_units)) noprint;
	by department_descr product_group_descr product_module_descr;
	var size;
	output out = basket MEDIAN = module_med_size;
run;

* Merge the product data with median size for other unit measurements;
proc sql noprint;
	create table tmp2 as 
	select A.*, B.module_med_size
	from tmp as A left join basket as B
	on A.product_module_descr = B.product_module_descr
	order by department_descr, product_group_descr, product_module_descr, upc;
quit;

* If the unit measurement of a product is not the main measurement of the module, then we assign the medium size with known measurement; 
data tmp2;
	set tmp2;
	if size1_units ^= module_units then size = module_med_size;
	upcv = cats(upc, "-", upc_ver_uc);
run;

* Drop some UPC from magnet data;
proc sql noprint;
	create table tmp_products as 
	select *
	from tmp2 
	where upcv not in (select upcv from my_delete_upc)
	order by department_descr, product_group_descr, product_module_descr;
quit;

* Compute relative size index; 
data products;
	set tmp_products; 
	size_index = size/module_med_size; 
	PRVT_flag = 0; 
	if brand_descr="CTL BR" then PRVT_flag = 1;
	label	module_units 	= "Unit measure of category"
			module_med_size = "The median size of category measured at unit measure"
			upcv			= "UPC-version"
			size_index 		= "Size normalized by median size within category";
run;	
data mylib.products; set products; run; 

*-------------------------------------; 
* Compute wallet share of all modules ; 
*------------------------------------------------------*;
* Compute the dollar spent on non-scannable at each trip; 
* Save temporary datasets products to sql server; 
data mysqllib.products; 
	set products; 
run; 

* Compute dollars spent on magnet and nonscannables; 
proc sql; 
	connect to oledb as mydb
	(OLEDB_SERVICES=NO Datasource="kdc01\kdwh02" PROVIDER=SQLOLEDB.1 
	Properties=('Initial Catalog'=USRDB_ccv103 'Integrated Security'=SSPI)
	SCHEMA=DBO);
	
	create table nonscan_all as 	
	select * from connection to mydb (
		SELECT A.trip_code_uc, year, channel_type, retailer_code, total_spent, net_spent, dol_magnet, dol_purchase,
				(net_spent - dol_purchase + dol_magnet) as dol_nonscan
		FROM dbo.trips as A LEFT JOIN 
		(SELECT trip_code_uc, sum(dol_paid) as dol_purchase, sum(dol_paid * my_magnet) as dol_magnet
		FROM (select T.*, P.my_magnet from dbo.purchases as T inner join dbo.products as P on T.upcv = P.upcv) as tmp
		GROUP by trip_code_uc ) as B
		on A.trip_code_uc = B.trip_code_uc
	);
	
	* Apend product module to purchase data;
	create table tmp as
	select * from connection to mydb (
		select A.*, B.product_module_descr
		from dbo.annual_purchase as A inner join dbo.products as B
		on A.upcv = B.upcv;
	); 
	disconnect from mydb; 
quit; 

* For unmatched trips, set dollars spent to 0; 
data nonscan_all; 
	set nonscan_all;
	if dol_magnet = . then dol_magnet = 0; 
	if dol_purchase = . then dol_purchase = 0; 
	if dol_nonscan = . then dol_nonscan = 0; 
run;

proc sql noprint;
	* Sum up the dollars on non-scannables; 
	select sum(dol_nonscan) into: nonscan from nonscan_all;
	
	* Insert the dollars on non-scannables; 
	insert into tmp 
	set product_module_descr = "Nonscan", dol_paid = &nonscan; 
	
	* Compute wallet share of each module;
	create table module_wallet as
	select product_module_descr, dol_paid/sum(dol_paid) as weight
	from (select product_module_descr, sum(dol_paid) as dol_paid from tmp group by product_module_descr) 
	order by product_module_descr;
quit;

* Find the median value of dollar spent on non-scannables over all trips; 
proc univariate data = nonscan_all;
	var dol_nonscan; 
	output out = tmp1 MEDIAN = module_med_size;
run;
data tmp1; set tmp1; product_module_descr = "Nonscan"; run;

* Merge module weight to basket; 
proc sort data = basket; by product_module_descr; 
data basket; 
	merge basket module_wallet tmp1;
	by product_module_descr;
run; 
proc contents data = basket; run;

data mylib.basket; set basket; run; 
data mysqllib.basket; set basket; run;

proc datasets noprint; delete tmp tmp1 tmp2 tmp_products module_measure module_wallet my_products; run;

*********************************************************************************;
* Compute retail attributes at market-retail format level using annual purchase *; 
*********************************************************************************;
proc sql noprint;
	create table retailer_transaction as
	select scantrack_market_descr, year, channel_type, retailer_code, 
			sum(total_spent) as revenue
	from trips
	group by scantrack_market_descr, year, channel_type, retailer_code;
quit; 

* Compute format-level attrubutes; 
proc sql noprint;
	* Append the upc-related information from products data;
	create table retailer_upc as 
	select A.*, B.department_descr, B.product_group_descr, B.product_module_descr, B.brand_descr, B.size_index,
			B.size1_amount, B.size1_units, B.multi, B.PRVT_flag, quantity*size1_amount*multi as vol
	from annual_purchase as A inner join products as B
	on A.upcv = B.upcv;
	
	* Aggregate to retailer-module level weighted by UPC volumn at each retailer;
	create table retailer_module as 
	select scantrack_market_descr, year, channel_type, retailer_code, department_descr, product_group_descr, product_module_descr,
			sum(size_index*vol)/sum(vol) as size_index, 
			count(distinct brand_descr) as num_brands, count(distinct upcv) as num_upc, 
			mean(PRVT_flag) as PRVT, sum(PRVT_flag) as num_prvt,
			sum(quantity) as module_units, sum(vol) as module_vol, sum(dol_paid) as revenue
	from retailer_upc
	group by scantrack_market_descr, year, channel_type, retailer_code, department_descr, product_group_descr, product_module_descr;
	
	* Aggregate retailer level index using wallet share as module weight;
	create table retailer_size as 
	select scantrack_market_descr, year, channel_type, retailer_code,
		sum(size_index*weight)/sum(weight) as size_index,
		sum(num_brands) as num_brands, sum(num_brands*weight)/sum(weight) as avgn_brands,
		sum(num_upc) as num_upc, sum(num_upc*weight)/sum(weight) as avgn_upc,
		sum(PRVT*weight)/sum(weight) as avg_prvt,
		sum(num_prvt)/sum(num_brands) as overall_prvt,
		count(distinct product_module_descr) as num_module
	from (select A.*, B.weight from retailer_module as A left join basket as B 
		  on A.product_module_descr = B.product_module_descr)
	group by scantrack_market_descr, year, channel_type, retailer_code
	order by scantrack_market_descr, year, channel_type, retailer_code;
quit; 

* Compute the format level index from retailer-level index weighted by transaction number;
proc sql noprint;		
	* Channel aggregation 2: Aggregate to market-channel level from retailer data weighted by store revenue;
	create table channel_size as 
	select scantrack_market_descr, year, channel_type, sum(size_index*revenue)/sum(revenue) as size_index, 
			sum(num_module*revenue)/sum(revenue) as num_module,
			sum(num_brands*revenue)/sum(revenue) as num_brands,
			sum(avgn_brands*revenue)/sum(revenue) as avg_brands_per_mod,
			sum(num_upc*revenue)/sum(revenue) as num_upc,
			sum(avgn_upc*revenue)/sum(revenue) as avg_upc_per_mod,
			sum(avg_prvt*revenue)/sum(revenue) as avg_prvt_per_mod, 
			sum(overall_prvt*revenue)/sum(revenue) as overall_prvt
	from (select A.*, B.revenue
			from retailer_size as A left join retailer_transaction as B
			on A.scantrack_market_descr=B.scantrack_market_descr and A.year=B.year and A.retailer_code = B.retailer_code)
	group by scantrack_market_descr, year, channel_type
	order by scantrack_market_descr, year, channel_type;
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
			overall_prvt	= "Proportion of private labels among all brands";
run;
proc datasets noprint; delete tmp retailer_upc retailer_module retailer_size; run;

********************************************************;
* Compute biweekly price at market-retail format level *;
********************************************************;
* Compute average category price for each market-chain; 
proc sql noprint;
	* Append the upc-related information from products data;
	create table retailer_upc as 
	select A.*, B.department_descr, B.product_group_descr, B.product_module_descr, B.brand_descr, B.size_index,
			B.size1_amount, B.size1_units, B.multi, 
			quantity*size_index as vol, price_tag/size_index as unitprice_tag, price_paid/size_index as unitprice_paid
	from biweek_purchase as A inner join products as B
	on A.upcv = B.upcv;
	
	* Aggregate to retailer-module level weighted by UPC volumn at each retailer;
	create table retailer_module as 
	select scantrack_market_descr, biweek, channel_type, retailer_code, department_descr, product_group_descr, product_module_descr,
			min(year) as year, 
			sum(price_tag*vol)/sum(vol) as price_tag_volwt, sum(price_paid*vol)/sum(vol) as price_paid_volwt,
			sum(unitprice_tag*vol)/sum(vol) as unitprice_tag, sum(unitprice_paid*vol)/sum(vol) as unitprice_paid,			
			sum(quantity) as module_units, sum(vol) as module_vol, sum(dol_paid) as revenue, sum(dol) as dol
	from retailer_upc
	group by scantrack_market_descr, biweek, channel_type, retailer_code, department_descr, product_group_descr, product_module_descr;
quit; 

* Compute chain level price index; 
proc sql noprint; 
	* Aggregate retailer level index using wallet share as module weight;
	create table retailer_index as 
	select scantrack_market_descr, biweek, channel_type, retailer_code, min(year) as year, 
		sum(price_tag_volwt*weight)/sum(weight) as price_tag_index, 
		sum(price_paid_volwt*weight)/sum(weight) as price_paid_index, 
		sum(unitprice_tag*weight)/sum(weight) as bsk_price_tag, 
		sum(unitprice_paid*weight)/sum(weight) as bsk_price_paid 
	from (select A.*, B.weight from retailer_module as A left join basket as B 
		  on A.product_module_descr = B.product_module_descr)
	group by scantrack_market_descr, biweek, channel_type, retailer_code
	order by scantrack_market_descr, year, biweek, channel_type, retailer_code;
quit;

* Aggregate from chain-level to format-level price index; 
proc sql noprint;	
	* Channel aggregation: Aggregate to market-channel level from retailer data weighted by store revenue;
	create table channel_index as 
	select scantrack_market_descr, min(year) as year, biweek, channel_type, 
			sum(price_tag_index*revenue)/sum(revenue) as price_tag_index,
			sum(price_paid_index*revenue)/sum(revenue) as price_paid_index,
			sum(bsk_price_tag*revenue)/sum(revenue) as bsk_price_tag, 
			sum(bsk_price_paid*revenue)/sum(revenue) as bsk_price_paid
	from (select A.*, B.revenue from retailer_index as A left join retailer_transaction as B
			on A.scantrack_market_descr=B.scantrack_market_descr and A.year=B.year 
			and A.retailer_code = B.retailer_code)
	group by scantrack_market_descr, biweek, channel_type
	order by scantrack_market_descr, year, biweek, channel_type;
quit;

* Merge in the CPI to have price at common dollars; 
proc sql noprint;
	select annual into: cpi2004 from cpi where year=2004 ;
	
	create table tmp as 
	select A.*, (B.annual/&cpi2004) as cpi
	from channel_index as A left join CPI as B
	on A.year = B.year
	order by scantrack_market_descr, year, biweek, channel_type;
quit;

data mylib.format_biweek_price; 
	set tmp; 
	price_tag_2004		= price_tag_index/cpi;
	price_paid_2004		= price_paid_index/cpi;
	bsk_price_tag_2004	= bsk_price_tag/cpi; 
	bsk_price_paid_2004	= bsk_price_paid/cpi; 
	label	price_tag_index			= "Average price tag per item"
			price_paid_index		= "Average price paid per item"
			bsk_price_tag			= "Price tag for a whole basket"
			bsk_price_paid			= "Price paid for a whole basket"
			price_tag_2004			= "Average price tag per item at 2004 dollar"
			price_paid_2004			= "Average price paid per item at 2004 dollar"
			bsk_price_tag_2004		= "Price tag for a whole basket at 2004 dollar"
			bsk_price_paid_2004		= "Price paid for a whole basket at 2004 dollar";
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

PROC EXPORT DATA=mylib.format_biweek_price
			OUTFILE= "\\tsclient\Resear1\Store switching\processed data\1_format_biweek_price.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;
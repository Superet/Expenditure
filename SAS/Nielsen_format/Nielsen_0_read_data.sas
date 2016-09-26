/*
This script reads data from Nielsen Homescan Panel. 
Outupt library: mainex
Output data: panelists, trips, purchases, products
*/

LIBNAME mysqllib OLEDB
OLEDB_SERVICES=NO
Datasource="kdc01\kdwh02"
PROVIDER=SQLOLEDB.1
Properties=('Initial Catalog'=USRDB_ccv103
                  'Integrated Security'=SSPI)
SCHEMA=DBO
BULKLOAD=YES
bl_options='ROWS_PER_BATCH = 500000';

libname mainex "U:\Users\ccv103\Documents\Research\Store switching\SAS_temp\Nielsen\Work\main";

%let dirname= U:\Users\ccv103\Documents\Research\Nielsen\nielsen_extracts\HMS;
options compress=yes reuse=yes; 

*************;
* Read data *;
*************;
* Read in the data of deleted UPC;
PROC IMPORT OUT= my_delete_upc
            DATAFILE= "U:\Users\ccv103\Documents\Research\Store switching\SAS_temp\deleted_upc.csv" 
            DBMS=csv REPLACE;
     		GETNAMES=YES; 
RUN;

* Read in CPI data;
PROC IMPORT OUT= CPI
            DATAFILE= "U:\Users\ccv103\Documents\Research\Store switching\SAS_temp\CPI.csv" 
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
	biweek = ceil((purchase_date - MDY(12,31,2003))/14); 
run;
data mysqllib.purchases (BULKLOAD=yes); set mainex.purchases; run;
data mysqllib.panelists (BULKLOAD=yes); set mainex.panelists; run;

***********************;
* Organize product data; 
***********************;
* See what belongs to magnet modules; 
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

* Flag magnet data; 
data my_products;
	set orig_products;
	if product_module_code >= 445 & product_module_code <= 468 then department_descr = "MAGNET DATA"; 
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
	prvt_flag = 0; 
	if brand_descr="CTL BR" then prvt_flag = 1;
	label	module_units 	= "Unit measure of category"
			module_med_size = "The median size of category measured at unit measure"
			upcv			= "UPC-version"
			size_index 		= "Size normalized by median size within category";
run;	

* Save product data; 
data mainex.products; set products; run; 
data mysqllib.products; set products; run;

* Adjust unit difference (Revised -- 20160910): LI - ML, QT - OZ;
data products;
	set mysqllib.products;
	if size1_units = 'LI' and module_units = 'ML' then size = size1_amount*multi*1000; 
	if size1_units = 'ML' and module_units = 'LI' then size = size1_amount*multi/1000;
	if size1_units = 'QT' and module_units = 'OZ' then size = size1_amount*multi*32;
	if size1_units = 'OZ' and module_units = 'QT' then size = size1_amount*multi/32;
	size_index = size/module_med_size; 
run;

proc datasets library=mysqllib; delete products; run;
proc datasets library=mainex; delete products; run;
data mainex.products; set products; run; 
data mysqllib.products; set products; run;


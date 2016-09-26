/*
This script generates the household-format-department data; 
Output library: mylib
Output data: hh_format_dpt

Input data: mylib.hh_basket, mylib.channel_module. 

NOTE: when aggregating from product moduels to departments, module weights are household-specific. 

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
libname mylib "U:\Users\ccv103\Documents\Research\Store switching\SAS_temp\Nielsen\Work\mylib";

options compress=yes reuse=yes; 

data mysqllib.hh_basket (BULKLOAD=yes); set mylib.hh_basket; run;
/*data mysqllib.retailer_location (BULKLOAD=yes); set mylib.retailer_location; run;
*/
* Merge household-specific basket to format-category data; 
proc sql noprint; 
	connect to oledb as mydb
	(OLEDB_SERVICES=NO Datasource="kdc01\kdwh02" PROVIDER=SQLOLEDB.1 
	Properties=('Initial Catalog'=USRDB_ccv103 'Integrated Security'=SSPI)
	SCHEMA=DBO);
	
	EXECUTE(
		select A.household_code, A.wallet_share, B.*
		into dbo.tmp
		from (select C.*, D.scantrack_market_descr, D.panel_year 
			 from dbo.hh_basket as C left join dbo.panelists as D on C.household_code = D.household_code
			 where wallet_share > 0) as A 
		full join (select *, 
				CASE WHEN department_descr in ('DAIRY', 'DELI', 'FROZEN FOODS', 'PACKAGED MEAT') THEN 'REFRIGERATED FROZEN'
					 WHEN department_descr = 'ALCOHOLIC BEVERAGES' THEN 'DRY GROCERY'
					 WHEN department_descr in ('FRESH PRODUCE', 'MAGNET DATA') THEN 'PRODUCE OTHER'
					 WHEN department_descr = 'HEALTH & BEAUTY CARE' THEN 'HEALTH BEAUTY CARE'
					 WHEN department_descr = '' THEN 'UNCLASSIFIED'
				ELSE department_descr
				END as department
			 from dbo.channel_module where department_descr != '') as B 
		on A.scantrack_market_descr = B.scantrack_market_descr and A.panel_year = B.year and 
		A.product_module_descr = B.product_module_descr;

		/** Aggregate across categories within a department; 	*/
		select household_code, scantrack_market_descr, year, channel_type, department, 
				sum(unitprice_tag*module_med_size*wallet_share)/sum(wallet_share) as unitprice_tag, 
				sum(unitprice_paid*module_med_size*wallet_share)/sum(wallet_share) as unitprice_paid, 
				count(distinct product_module_descr) as num_module, 
				sum(num_brands) as num_brands, sum(num_upc) as num_upc, 
				sum(size_index*wallet_share)/sum(wallet_share) as size_index, 
				sum(num_brands*wallet_share)/sum(wallet_share) as num_brands_per_mod, 
				sum(num_upc*wallet_share)/sum(wallet_share) as num_upc_per_mod, 
				sum(prvt_pen*wallet_share)/sum(wallet_share) as prvt_pen_mod,  
				sum(num_prvt_upc)/sum(num_upc) as prvt_overall
		into dbo.hh_format_dpt		
		from (select A.*,module_med_size from dbo.tmp as A left join 
				(select distinct product_module_descr, module_med_size from dbo.products) as B 
				on A.product_module_descr = B.product_module_descr) as C
		group by household_code, scantrack_market_descr, year, channel_type, department
		order by household_code, scantrack_market_descr, year, channel_type, department; 

		drop table dbo.tmp, dbo.hh_basket; 
	) by mydb; 
	disconnect from mydb;	
quit; 

proc sort data = mysqllib.hh_format_dpt out = mylib.hh_format_dpt; 
	by household_code year channel_type department; 
run;

proc datasets library = mysqllib noprint; delete hh_format_dpt; run;


/*
This script generates the depement-level metrics for each format in a market in a year; 
Output library: mylib
Output data: format_department

Input data: mysqllibr.annual_purchase, mysqllib.products, mysqllib.channel_module
Intermediary data: mysqllib.basket

NOTE: Module weight is the wallet share of revenue sales from all the data, and it is fixed across stores,formats,markets and time. 

Original department: ALCOHOLIC BEVERAGES, DAIRY, DELI, DRY GROCERY, FRESH PRODUCE, FROZEN FOODS, GENERAL MERCHANDISE, 
					HEALTH & BEAUTY CARE, MAGNET DATA, NON-FOOD GROCERY, PACKAGED MEAT. 
Re-classified department: 	REFRIGERATED & FROZEN (DAIRY, DELI, FROZEN FOODS, PACKAGED MEAT), 
							DRY GROCERY (ALCOHOLIC BEVERAGES), 
							NON-FOOD GROCERY, 
							GENERAL MERCHANDISE, 
							HEALTH & BEAUTY CARE, 
							PRODUCE & OTHER (FRESH PRODUCE, MAGNET DATA)
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

%let dirname= U:\Users\ccv103\Documents\Research\Nielsen\nielsen_extracts\HMS;
options compress=yes reuse=yes;

* Fixed basket -- common wallet share *; 
* Compute wallet share of each module;
proc sql noprint; 
	connect to oledb as mydb
	(OLEDB_SERVICES=NO Datasource="kdc01\kdwh02" PROVIDER=SQLOLEDB.1 
	Properties=('Initial Catalog'=USRDB_ccv103 'Integrated Security'=SSPI)
	SCHEMA=DBO);
	
	EXECUTE(
		select department, product_module_descr, module_med_size, dol_paid/sum(dol_paid) as weight
		into dbo.basket
		from (select department_descr, product_module_descr, module_med_size, sum(dol_paid) as dol_paid 
			from (select A.*, B.product_module_descr, B.module_med_size, B.department_descr
					from dbo.annual_purchase as A inner join dbo.products as B on A.upcv = B.upcv)
		group by product_module_descr)
		order by department_descr, product_module_descr; 
	) by mydb; 	
quit; 

/*Within a department, summarize channel characteristics by aggregating across all modules; */
proc sql noprint;
	create table format_department as 
	select scantrack_market_descr, year, channel_type, department, 
		sum(unitprice_tag*module_med_size*weight)/sum(weight) as unitprice_tag, 
		sum(unitprice_paid*module_med_size*weight)/sum(weight) as unitprice_paid, 
		count(unique(product_module_descr)) as num_module, 
		sum(num_brands) as num_brands, sum(num_upc) as num_upc, 
		sum(size_index*weight)/sum(weight) as size_index, 
		sum(num_brands*weight)/sum(weight) as num_brands_per_mod, 
		sum(num_upc*weight)/sum(weight) as num_upc_per_mod, 
		sum(prvt_pen*weight)/sum(weight) as prvt_pen_mod,  
		sum(num_prvt_upc)/sum(num_upc) as prvt_overall
	from (select A.*, B.weight, B.module_med_size
		from mysqllib.channel_module as A inner join 
		(select *, CASE WHEN department_descr in ('DAIRY', 'DELI', 'FROZEN FOODS', 'PACKAGED MEAT') THEN 'REFRIGERATED FROZEN'
			 WHEN department_descr = 'ALCOHOLIC BEVERAGES' THEN 'DRY GROCERY'
			 WHEN department_descr in ('FRESH PRODUCE', 'MAGNET DATA') THEN 'PRODUCE OTHER'
			 WHEN department_descr = 'HEALTH & BEAUTY CARE' THEN 'HEALTH BEAUTY CARE'
			 WHEN department_descr = '' THEN 'UNCLASSIFIED'
		ELSE department_descr
		END as department from mysqllib.basket) as B 
		on A.product_module_descr = B.product_module_descr
		where A.department_descr ^= "")
	group by scantrack_market_descr, year, channel_type, department
	order by scantrack_market_descr, year, channel_type, department; 
	
/*Add number of stores in each market; */
	create table mylib.format_department as
	select A.*, num_store
	from format_department as A left join
	(select scantrack_market_descr, year, channel_type, count(unique(store_code_uc)) as num_store 
	from mylib.store_location group by scantrack_market_descr, year, channel_type ) as B
	on A.scantrack_market_descr = B.scantrack_market_descr and A.year = B.year and A.channel_type = B.channel_type
	order by scantrack_market_descr, year, channel_type, department;
quit; 

proc univariate data = format_department;
	var unitprice_paid;
run;

proc corr data = mylib.format_department(drop = scantrack_market_descr year channel_type department); run; 


/*
This script calculates the category-level metrics in each format;
Output library: mylib
Output data: channel_module

Intermediate data: mysqllib.annual_purchase, mysqllib.retailer_transaction
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

*********************************************;
* Summarize annual sales from purchase data *; 
*********************************************;
* Compute the dollar spent on non-scannalbes for each trip;
proc sql; 
	connect to oledb as mydb
	(OLEDB_SERVICES=NO Datasource="kdc01\kdwh02" PROVIDER=SQLOLEDB.1 
	Properties=('Initial Catalog'=USRDB_ccv103 'Integrated Security'=SSPI)
	SCHEMA=DBO);

	EXECUTE(
		/*Aggregate annual sales for each upc at each store; */
		select scantrack_market_descr, year, channel_type, retailer_code, upcv, 
				sum(total_price_paid) as dol, sum(quantity) as quantity,
				sum(dol_paid) as dol_paid, sum(total_price_paid)/sum(quantity) as price_tag, sum(dol_paid)/sum(quantity) as price_paid
		into dbo.annual_purchase
		from dbo.purchases
		group by scantrack_market_descr, channel_type, retailer_code, upcv, year;
	) by mydb;
	
	%put &SQLXRC. &SQLXMSG. ; 
	disconnect from mydb;
quit; 

* Check the number of households in each scantrack market; 
proc sql noprint;
	create table tmp as 
	select scantrack_market_descr, count(unique(household_code)) as freq
	from mainex.panelists
	group by scantrack_market_descr
	order by freq; 
quit; 
proc print data = tmp; run;
proc datasets noprint; delete tmp; run;

******************************;
* Summarize category metrics *; 
******************************;
proc sql;
	connect to oledb as mydb
	(OLEDB_SERVICES=NO Datasource="kdc01\kdwh02" PROVIDER=SQLOLEDB.1 
	Properties=('Initial Catalog'=USRDB_ccv103 'Integrated Security'=SSPI)
	SCHEMA=DBO);
	 
	EXECUTE(
		/*Annual chain revenue as weight weithin market format weight;*/
		select scantrack_market_descr, year, channel_type, retailer_code, 
				sum(total_spent) as revenue
		into dbo.retailer_transaction
		from dbo.trips
		group by scantrack_market_descr, year, channel_type, retailer_code;
		
		/*Aggregate to retailer-module level weighted by UPC volumn at each retailer;*/
		select scantrack_market_descr, year, channel_type, retailer_code, department_descr, product_module_descr,
				sum(dol)/sum(vol) as unitprice_tag, sum(dol_paid)/sum(vol) as unitprice_paid,
				sum(size_index*vol)/sum(vol) as size_index, 
				count(distinct brand_descr) as num_brands, count(distinct upcv) as num_upc, 
				sum(PRVT_flag) as num_prvt_upc, avg(prvt_flag) as prvt_pen
		into dbo.retailer_module
		from (select A.*, B.department_descr, B.product_module_descr, B.brand_descr, B.size_index,
				B.size1_amount, B.size1_units, B.prvt_flag, quantity*size as vol
				from dbo.annual_purchase as A 
				inner join 
				dbo.products as B on A.upcv = B.upcv) as retailer_upc
		group by scantrack_market_descr, year, channel_type, retailer_code, department_descr, product_module_descr;
		
		select scantrack_market_descr, year, channel_type, department_descr, product_module_descr, 
				sum(unitprice_tag*revenue)/sum(revenue) as unitprice_tag, 
				sum(unitprice_paid*revenue)/sum(revenue) as unitprice_paid, 
				sum(size_index*revenue)/sum(revenue) as size_index, 
				sum(num_brands*revenue)/sum(revenue) as num_brands,
				sum(num_upc*revenue)/sum(revenue) as num_upc, 
				sum(num_prvt_upc*revenue)/sum(revenue) as num_prvt_upc, 
				sum(prvt_pen*revenue)/sum(revenue) as prvt_pen
		into dbo.channel_module
		from (select A.*, B.revenue
				from dbo.retailer_module as A left join dbo.retailer_transaction as B
				on A.scantrack_market_descr=B.scantrack_market_descr and A.year=B.year and A.retailer_code = B.retailer_code) as tmp1
		group by scantrack_market_descr, year, channel_type, department_descr, product_module_descr
		order by scantrack_market_descr, year, channel_type, department_descr, product_module_descr;

		drop table dbo.retailer_module; 
	) by mydb; 
	
	%put &SQLXRC. &SQLXMSG. ; 
	disconnect from mydb;
quit; 	

proc univairate data = mysqllib.channel_module; var unitprice_paid; run;

data mylib.channel_module; set mysqllib.channel_module; run;


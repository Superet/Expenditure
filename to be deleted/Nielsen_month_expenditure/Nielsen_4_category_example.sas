libname mylib "E:\Users\ccv103\Documents\Research\Store switching\SAS_temp\Nielsen";
libname mainex "E:\Users\ccv103\Documents\Research\Store switching\SAS_temp\Nielsen\Work\main";

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

* NOTE: we extract purchase data from only grocery stores;  
proc sql noprint;
	create table retailers as
	select A.*, B.channel_group
	from orig_retailers as A inner join my_channel as B
	on A.channel_type = B.channel_type;
quit;

* UPC from the selected product module; 
data my_module;
	input product_module_descr $30.;
	datalines;
	TOILET TISSUE
	CEREAL - READY TO EAT
	CHEESE - SHREDDED
	SOFT DRINKS - CARBONATED
	;
run;

proc sql noprint;
	create table my_products as 
	select *
	from mylib.products 
	where product_module_descr in (select product_module_descr from my_module);
quit;

* Read trip data;
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
		 	* Subset the panelists to the selected areas;
			create table tmp as 
			select A.*, B.market, B.mkt_type
			from panelists&i as A inner join my_market as B
			on A.scantrack_market_descr = B.scantrack_market_descr;
						
			create table tmp1 as 
			select A.*, year(A.purchase_date) as year, B.market, B.scantrack_market_descr, B.mkt_type, 
					B.projection_factor_magnet, B.magnet_flag 
			from trips&i as A inner join tmp as B
			on A.household_code = B.household_code;	
				
			* Restrict the trip data to the selected retailers;
			create table tmp2 as
			select A.*, B.channel_group, B.channel_type
			from tmp1 as A inner join retailers as B
			on A.retailer_code = B.retailer_code
			order by trip_code_uc;	
			
			* Subset purchases data within the restrictions above;
			create table tmp3 as 
			select A.*, cats(upc, "-", upc_ver_uc) as upcv, (total_price_paid-coupon_value) as dol_paid, 
					total_price_paid/quantity as price_tag, (total_price_paid-coupon_value)/quantity as price,
					B.household_code, B.mkt_type, B.market, B.scantrack_market_descr, B.channel_type, B.retailer_code,
					year(B.purchase_date) as year, month(B.purchase_date) as month
			from purchases&i as A inner join tmp2 as B
			on A.trip_code_uc = B.trip_code_uc;	
			
			* Subset the selected product module;
			create table purchases&i as 
			select A.*, B.product_module_descr, B.brand_descr, B.upc_descr, B.multi, B.size1_amount, B.size_index, B.size
			from tmp3 as A inner join my_products as B
			on A.upcv = B.upcv;
		quit;

		proc datasets noprint; delete tmp tmp1 tmp2 tmp3 panelists2004 - panelists2012 trips2004 - trips2012 ; run;
	%end;
		
	data purchases;
		set purchases2004 - purchases2012;
	run;
	
	proc datasets noprint; 
		delete purchases2004 - purchases2012;
	run;
%mend read_groc_fun;

%read_groc_fun;

data purchases;
	set purchases;
	date = mdy(month, 1, year);
	unit_price = price/size;
run;

**********************************;
* Price series from a single UPC *;
**********************************;
%let sel_upcv = 1200000129-1;
%let sel_mkt  = Boston;

proc sql noprint;
	create table tmp as 
	select market, year, month, date, channel_type, mean(price) as price
	from (select * from purchases where upcv = "&sel_upcv" and market = "&sel_mkt")
	group by market, year, month, date, channel_type
	order by market, year, month, channel_type; 
	
	select price into: price_base from tmp where year=2004 and month = 1 and channel_type="Grocery";
quit;

data tmp;
	set tmp;
	price_base = &price_base;
	price_norm = price/price_base;
run;

proc sgplot data=tmp;
	series x=date y=price_norm / group = channel_type;
	title "Price series for UPC=&sel_upcv in &sel_mkt";
run;

proc datasets noprint; delete tmp; run;

***********************************************************;
* Price series of average unit price of a single category *;
***********************************************************;
%let sel_module = CEREAL - READY TO EAT;
%let sel_mkt	= Denver;

* Subset the data;
proc sql noprint;
	* Compute the monthly price and quantity from the subset;
	create table tmp as 
	select market, year, month, date, channel_type, upcv, mean(price/size) as unit_price, sum(dol_paid) as dol_paid, sum(quantity*size) as vol
	from (select * from purchases where product_module_descr = "&sel_module" and market = "&sel_mkt")
	group by market, year, month, date, channel_type, upcv;
	
	* Compute the market-level volume for each UPC;
	create table tmp1 as 
	select market, date, upcv, sum(vol) as mkt_vol
	from tmp
	group by market, date, upcv;
	
	* Mergin in volume data; 
	create table tmp2 as 
	select A.*, B.mkt_vol
	from tmp as A left join tmp1 as B
	on A.upcv=B.upcv and A.date=B.date;
	
	create table tmp3 as 
	select market, year, month, date, channel_type,
	 		mean(unit_price) as unit_price, 
			sum(unit_price*vol)/sum(vol) as unit_price_channelvol, 
			sum(unit_price*mkt_vol)/sum(mkt_vol) as unit_price_mktvol
	from tmp2
	group by market, year, month, date, channel_type;
quit;

proc sgplot data=tmp3;
	series x = date y = unit_price / group = channel_type; 
	title "Average unit price for &sel_module in &sel_mkt";
run;

proc sgplot data=tmp3;
	series x = date y = unit_price_channelvol / group = channel_type; 
	title "Average unit price weighted by channel volume for &sel_module in &sel_mkt";
run;

proc sgplot data=tmp3;
	series x = date y = unit_price_mktvol / group = channel_type; 
	title "Average unit price weighted by market volume for &sel_module in &sel_mkt";
run;

proc datasets noprint; delete tmp tmp1 tmp2 tmp3; run;

****************************************;
* Product size histogram in a category *;
****************************************;
%let sel_module = CEREAL - READY TO EAT;
%let sel_mkt	= Denver;

* Subset the data;
proc sql noprint;
	* Compute the monthly price and quantity from the subset;
	create table tmp as 
	select market, year, channel_type, upcv, mean(price) as price, mean(size) as size
	from (select * from purchases where product_module_descr = "&sel_module" and market = "&sel_mkt")
	group by market, year, channel_type, upcv;
quit;

data tmp;
	set tmp;
	module = "&sel_module";
run;

PROC EXPORT DATA=tmp
			OUTFILE= "\\tsclient\Resear1\Store switching\processed data\sas_category_example.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;


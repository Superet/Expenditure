/*
This script generates the household-format data (Distance between households and the nearest stores within a format)
Output library: mylib
Output data: hh_format

Input data: mainex.panelists, mylib.keephh, mylib.store_location, mylib.store_location_app_excl

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

********************************;
* Organize store location data *; 
********************************;
* Append year for the static store location; 
data tmp1; 
	set mylib.store_location_app_excl(where = (year=.)); 
	year = 2003; 
	do while (year <= 2012); 
		year+1; output; 
	end; 
run;
data store_location_app; 
	set mylib.store_location_app_excl(where = (year^=.)) tmp1; 
run;
proc sql; 
	select nmiss(zipcode)/count(*) as nzip, nmiss(latitude)/count(*) as ncoord, sum(zipcode=. and latitude=.)/count(*) as nboth
	from store_location_app;
quit; 

* Use the coordinates of the zipcode to impute the missing latitude; 
proc sql noprint;
	create table tmp as 
	select store_name, zipcode,  B.Y as latitude, B.X as longitude, year, store_id, store_zip3, channel_type
	from (select * from store_location_app where latitude = . and zipcode^=.) as A 
	left join SASHELP.zipcode as B
	on A.zipcode = B.zip;
quit; 

data tmp; 
	set store_location_app(where = (latitude^=.)) tmp; 
run;
		
* Add fips code and state to each zipcode; 
proc sql noprint;
	create table store_location_app as 
	select A.*, state, cats(put(state, z2.), put(county, z3.)) as fips, city
	from (select * from tmp where not (zipcode=. and latitude=.)) as A left join SASHELP.zipcode as B
	on A.zipcode = B.zip;
quit;
	
* Cross join the full set of fips-channel; 		
proc sql noprint;	
	create table tmp1 as 
	select * from 
	(select distinct state, fips from store_location_app where state^=.) as A 
	cross join
	(select distinct channel_type from store_location_app) as B
	order by state, fips, channel_type;  
	
	* Match fips with location; 
	create table tmp2 as 
	select A.*, store_id, latitude
	from tmp1 as A left join store_location_app as B
	on A.fips = B.fips and A.channel_type = B.channel_type; 
	
	* For the fips that do not have all the channels, use the stores in the same state; 
	create table tmp3 as 
	select store_name, zipcode, B.latitude, longitude, year, store_id, store_zip3, A.channel_type, A.state, A.fips, B.city
	from (select * from tmp2(keep = state fips channel_type latitude) where latitude =. ) as A
	left join store_location_app as B
	on A.state = B.state and A.channel_type = B.channel_type; 
quit; 

* Fill the full set of channels;
data store_location_app; 
	set store_location_app(drop=store_zip3) tmp3(drop=store_zip3); 
	store_zip3 = substr(put(zipcode, z5.), 1, 3);
run;
proc datasets noprint; delete tmp tmp1 tmp2 tmp3; run;

* ----------------------------------------*
* Merge panelist data with store location *; 
proc sql noprint; 
	* Merge households zip code with coordinates; 
	create table hh as 
	select household_code, scantrack_market_descr, panel_year, panelist_zip_code, hh_zip3, 
	D.X as hh_longitude, D.Y as hh_latitude, city as hh_city, cats(put(state, z2.), put(county, z3.)) as hh_fips
	from (	select *, substr(put(panelist_zip_code, z5.), 1, 3) as hh_zip3 
			from mainex.panelists where household_code in (select household_code from mylib.keephh)) as C 
	left join SASHELP.ZIPCODE as D on C.panelist_zip_code = D.zip; 
	
	* Subset the unique zipcode and coordinates for a channel in a market in a year from the store location data; 
	* NOTE: SAS zipde can not match all the zipcodes in household panel, so I drop those that could not get matched; 
	create table tmp1 as
	select distinct year, channel_type, store_zip3,	zipcode, longitude, latitude, 
			fips, city
	from store_location_app
	where fips^="";  

	* Merge household data and subsetted store location data; 	
	* NOTE: there are a few fips that do not have any stores and those fips are in the rare areas;	
	create table tmp2 as 
	select household_code, panelist_zip_code, hh_zip3, hh_longitude, hh_latitude, hh_city, B.*
	from hh as A inner join tmp1 as B 
	on A.hh_fips = B.fips and A.panel_year = B.year
	order by household_code, year, channel_type; 
	
	* For a household, flag if there is any store in the same zipcode or in the same areacode; 
	create table tmp3 as 
	select *, sum(hh_zip3 = store_zip3) as n_anyzip3, sum(panelist_zip_code = zipcode) as n_anyzip, sum(hh_city = city) as n_anycity
	from tmp2
	group by household_code, fips, year, channel_type
	order by household_code, fips, year, channel_type;
	
/*	Subset the merged data: 
	If there is at least one store in the same zipcode, then remove all the other zipcodes; 
	If there are no stores in the same zipcode but there is at least one in the same area code, then remove all the other areacodes; */
	create table hh_allzip as 
	select *
	from tmp3
/*	where (n_anycity > 0 and hh_city = city) or 
		  (n_anycity = 0 and n_anyzip3 > 0 and hh_zip3 = store_zip3) or 
		  (n_anycity = 0 and n_anyzip3 = 0)*/
	where (n_anyzip > 0 and panelist_zip_code = zipcode) or 
		  (n_anyzip = 0 and n_anycity > 0 and hh_city = city) or 
		  (n_anyzip = 0 and n_anycity = 0 and n_anyzip3 > 0 and hh_zip3 = store_zip3) or 
		  (n_anyzip = 0 and n_anycity = 0 and n_anyzip3 = 0)
	order by household_code, year, channel_type; 
quit; 
proc datasets noprint; delete tmp tmp1 tmp2 tmp3; run;

* Calculate the distance between household zipcode and the store location; 
* distance_haver is the Haversine distance between two coordinates;
* Need to convert from radian to degree, multiply by pi/180; 
%let pi180=0.0174532925199433;

proc sql noprint;
	create table channel_dist as 
	select distinct panelist_zip_code, year, channel_type, zipcode, latitude, longitude
	from hh_allzip(where = (hh_latitude^=.));
	
	create table dist as 
	select distinct panelist_zip_code, hh_latitude, hh_longitude, year, zipcode, latitude, longitude
	from hh_allzip; 
quit; 

data dist; 
	set dist; 
	if panelist_zip_code = zipcode then distance_zipcode = 0; 
	else do;
		if panelist_zip_code^=. and zipcode^=. then distance_zipcode = zipcitydistance(panelist_zip_code, zipcode); 	
	end; 	
	if latitude ^=. and hh_latitude ^=. then do; 
		a = sin((hh_latitude*&pi180-latitude*&pi180)/2)**2 + 
			cos(latitude*&pi180)*cos(hh_latitude*&pi180)*sin((longitude*&pi180-hh_longitude*&pi180)/2)**2;
		b =  2*atan2(sqrt(a), sqrt(1-a));
		distance_haver = 3959*b; 
	end; 
	drop a b; 
run;

proc corr data = dist; var distance_zipcode distance_haver; run;

* Find the shortest distance between a household and a retail format; 
proc sql noprint; 
	* Find the shortest distance between each zipcode and retailers for a format in a given year; 
	create table short_dist as 
	select panelist_zip_code, year, channel_type, min(distance_zipcode) as distance_zipcode, min(distance_haver) as distance_haver
	from (	select A.*, B.channel_type 
			from dist as A right join channel_dist as B
			on A.panelist_zip_code = B.panelist_zip_code and A.year = B.year and A.zipcode = B.zipcode and 
				A.latitude = B.latitude and A.longitude = B.longitude )
	group by panelist_zip_code, year, channel_type; 
	
	* Merge back zipcode distance and households; 
	* NOTE: 76 households can not matched zipcodes; 
	create table hh_format as 
	select A.household_code, B.*
	from hh as A inner join short_dist as B
	on A.panelist_zip_code = B.panelist_zip_code and A.panel_year = B.year
	order by household_code, year, channel_type; 
quit; 

* ------------------------------------------------------- *; 
* Check the correlation between two distance measurements *; 
proc corr data = hh_format; 
	var distance_zipcode distance_haver; 
run; 

data mylib.hh_dist; set hh_format; run; 

*********************************************************************************;
* Calculate the retail attributes for each format for household-specific bakset *; 
*********************************************************************************;
data mysqllib.hh_basket (BULKLOAD=yes); set mylib.hh_basket; run;

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
		full join (select * from dbo.channel_module where department_descr != '') as B 
		on A.scantrack_market_descr = B.scantrack_market_descr and A.panel_year = B.year and 
		A.product_module_descr = B.product_module_descr;

		/** Aggregate across categories within a department; 	*/
		select household_code, scantrack_market_descr, year, channel_type, 
				sum(unitprice_tag*module_med_size*wallet_share)/sum(wallet_share) as unitprice_tag, 
				sum(unitprice_paid*module_med_size*wallet_share)/sum(wallet_share) as unitprice_paid, 
				count(distinct product_module_descr) as num_module, 
				sum(num_brands) as num_brands, sum(num_upc) as num_upc, 
				sum(size_index*wallet_share)/sum(wallet_share) as size_index, 
				sum(num_brands*wallet_share)/sum(wallet_share) as num_brands_per_mod, 
				sum(num_upc*wallet_share)/sum(wallet_share) as num_upc_per_mod, 
				sum(prvt_pen*wallet_share)/sum(wallet_share) as prvt_pen_mod,  
				sum(num_prvt_upc)/sum(num_upc) as prvt_overall
		into dbo.hh_format_retattr		
		from (select A.*, module_med_size from dbo.tmp as A left join 
				(select distinct product_module_descr, module_med_size from dbo.products) as B 
				on A.product_module_descr = B.product_module_descr) as C
		group by household_code, scantrack_market_descr, year, channel_type
		order by household_code, scantrack_market_descr, year, channel_type; 
		
		drop table dbo.tmp, dbo.hh_basket, dbo.tmp_channel_module; 
	) by mydb; 
	disconnect from mydb;	
quit; 

proc sql; 
	select distinct channel_type from mylib.hh_dist; 
	select distinct channel_type from mysqllib.hh_format_retattr; 
quit; 

proc univariate data = mysqllib.hh_format_retattr; 
	var unitprice_paid; 
run;


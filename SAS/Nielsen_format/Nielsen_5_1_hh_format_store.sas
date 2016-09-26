/*
This script generates the household-format data (Distance between households and the nearest stores within a format)
Output library: mylib
Output data: hh_format

Input data: mainex.panelists, mylib.keephh, mylib.store_location

*/
libname mainex "U:\Users\ccv103\Documents\Research\Store switching\SAS_temp\Nielsen\Work\main";
libname mylib "U:\Users\ccv103\Documents\Research\Store switching\SAS_temp\Nielsen\Work\mylib";

options compress=yes reuse=yes; 

proc contents data = mainex.panelists; run;

* Merge panelist data with store location; 
proc sql noprint; 
	* Merge households zip code with coordinates; 
	create table tmp as 
	select household_code, scantrack_market_descr, panel_year, panelist_zip_code, 
	D.X as hh_longitude, D.Y as hh_latitude, D.AREACODE as hh_areacode
	from (select * from mainex.panelists where household_code in (select household_code from mylib.keephh)) as C 
	left join SASHELP.ZIPCODE as D on C.panelist_zip_code = D.zip; 
	
	* Subset the unique zipcode and coordinates for a channel in a market in a year from the store location data; 
	create table tmp1 as
	select A.*, B.areacode
	from (select distinct scantrack_market_descr, year, channel_type, 
		zipcode, round(longitude, 0.1) as longitude, round(latitude, 0.1) as latitude
		from mylib.store_location ) as A 
	left join SASHELP.ZIPCODE as B 
	on A.zipcode = B.zip; 
			
	* Merge household data and subsetted store location data; 		
	create table tmp2 as 
	select household_code, panelist_zip_code, hh_longitude, hh_latitude, hh_areacode, B.*
	from tmp as A left join tmp1 as B 
	on A.scantrack_market_descr = B.scantrack_market_descr and A.panel_year = B.year
	order by household_code, year; 
	
	* For a household, flag if there is any store in the same zipcode or in the same areacode; 
	create table tmp3 as 
	select *, sum(hh_areacode = areacode) as n_anyarea, sum(panelist_zip_code = zipcode) as n_anyzip
	from tmp2
	group by household_code, scantrack_market_descr, year, channel_type
	order by household_code, scantrack_market_descr, year, channel_type;
	
/*	Subset the merged data: 
	If there is at least one store in the same zipcode, then remove all the other zipcodes; 
	If there are no stores in the same zipcode but there is at least one in the same area code, then remove all the other areacodes; */
	create table hh_allzip as 
	select *
	from tmp3
	where (n_anyzip > 0 and panelist_zip_code = zipcode) or 
		  (n_anyzip = 0 and n_anyarea > 0 and hh_areacode = areacode) or 
		  (n_anyarea = 0)
	order by household_code, year, channel_type; 
	
	* Check how many households-year-channel has at least one store in the same zipcode;
	create table tmp as 
	select distinct household_code, scantrack_market_descr, year, channel_type, n_anyarea, n_anyzip
	from tmp3; 
quit; 

proc sql; select sum(n)
proc datasets noprint; delete tmp tmp1 tmp2 tmp3; run;

* Calculate the distance between household zipcode and the store location; 
* distance_haver is the Haversine distance between two coordinates;
* Need to convert from radian to degree, multiply by pi/180; 
%let pi180=0.0174532925199433;

data hh_allstore; 
	set hh_allstore; 
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
run;

* Find the shortest distance between a household and a retail format; 
proc sql noprint; 
	create table hh_format as 
	select household_code, year, channel_type, min(distance_zipcode) as distance_zipcode, min(distance_haver) as distance_haver
	from hh_allstore
	group by household_code, year, channel_type
	order by channel_type, household_code, year; 
quit; 

* ------------------------------------------------------- *; 
* Check the correlation between two distance measurements *; 
proc corr data = hh_format; 
	var distance_zipcode distance_haver; 
run; 

proc univariate data = hh_format; 
	group by channel_type; 
	var distance_zipcode distance_haver; 
	histogram; 
run; 

proc sort data = hh_format; by household_code year channel_type; run;
data mylib.hh_format; set hh_format; run; 


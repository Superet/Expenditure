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

* Check the number of households in each zip code; 
proc sql noprint; 
	create table tmp as 
	select panelist_zip_code, count(unique(household_code)) as num_hh
	from mainex.panelists
	group by panelist_zip_code; 
quit; 
proc univariate data = tmp; var num_hh; histogram; run;

* There are only 7 households in a zip code, not enough to impute penetration; 
* Consider the 3-digit zipcode; 
data tmp; 
	set tmp; 
	zip3 = substr(put(panelist_zip_code,z5.), 1, 3); 
run;
proc sql noprint;
	create table tmp1 as 
	select zip3, sum(num_hh) as num_hh
	from tmp group by zip3; 
quit; 
proc univariate data = tmp1; var num_hh; histogram; run;
proc datasets noprint; delete tmp tmp1; run;

* Summarize household-year trips to each store together with store zip code; 
proc sql noprint; 
	create table tmp as
	select household_code, year, channel_type, count(unique(purchase_date)) as num_trip
	from mainex.trips
	group by household_code, year, channel_type; 
	
	create table tmp1 as 
	select A.*, B.panelist_zip_code
	from tmp as A left join (select household_code, panelist_zip_code from mainex.panelists) as B
	on A.household_code = B.household_code;  
quit; 

data tmp1; 
	set tmp1; 
	zip3 = substr(put(panelist_zip_code,z5.), 1, 3);
run;

* Compute the channel penetration for each zip3; 
proc sql noprint;
	create table tmp as 
	select zip3, channel_type, sum(num_trip) as num_trip
	from tmp1
	group by zip3, channel_type;
	
	create table pen as 
	select zip3, channel_type, num_trip, num_trip/sum(num_trip) as penetration
	from tmp 
	group by zip3; 
quit; 

proc univariate data = pen(where=(channel_type = "Warehouse Club")); var penetration; histogram; run;

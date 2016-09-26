/*

This script imputes the store location from household trip data
Output library: mylib
Output data: store_location

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
proc printto log="U:\Users\ccv103\Desktop\templog.log"; run;

/*Zipcode data in 2010; */
/*proc cimport infile='U:\Users\ccv103\Downloadszipcode_oct10_v9.cpt'
lib=work; run;

proc append base=Zipcode_10q4_unique data=zipmil_10q4;
run;

proc sort data=zipcode_10q4_unique out=sasuser.zipcode; by zip zip_class; run;

proc datasets lib=sasuser;
modify zipcode;
index create zip;
run;*/

**********************************;
* List zipcode around each store *; 
**********************************;
* Check how many trips do not a store_code; 
* NOTE: store_code_uc = 0 when a store is unknown; 
proc sql; 
	select sum(store_code_uc=0)/count(*) from mainex.trips; 
	select sum(store_code_uc=0)/count(*) from (select * from mainex.trips where channel_type = "Discount Store"); 
	select sum(store_code_uc=0)/count(*) from (select * from mainex.trips where channel_type = "Drug Store"); 
	select sum(store_code_uc=0)/count(*) from (select * from mainex.trips where channel_type = "Grocery"); 
	select sum(store_code_uc=0)/count(*) from (select * from mainex.trips where channel_type = "Warehouse Club"); 
	select sum(store_code_uc=0)/count(*) from (select * from mainex.trips where channel_type = "Convenience Store"); 
	select sum(store_code_uc=0)/count(*) from (select * from mainex.trips where channel_type = "Dollar Store"); 
quit; 

* Find the customer size from each zip code area that is around a given store ;
proc sql noprint;
	/*First merge trip data with household zipcode, then find all the zipcodes around a store*/
	create table tmp as 	
	select scantrack_market_descr, store_code_uc, year, panelist_zip_code, count(distinct household_code) as num_hh
	from (select A.*, B.panelist_zip_code from mainex.trips(where=(year>2003 and store_code_uc > 0)) as A 
		left join mainex.panelists(where=(panelist_zip_code^=.)) as B
		 on A.household_code = B.household_code and A.year=B.panel_year) 
	group by scantrack_market_descr, store_code_uc, year, panelist_zip_code
	order by scantrack_market_descr, store_code_uc, year, num_hh;
	
	* count the number of households in each zipcode area; 
	create table store_allzip as 
	select *, sum(num_hh) as nhh_allyear
	from tmp(where = (panelist_zip_code^=.)) 
	group by scantrack_market_descr, store_code_uc, panelist_zip_code
	order by scantrack_market_descr, store_code_uc, nhh_allyear, panelist_zip_code; 
	
	* Extract 3-digit zip code in the data; 
	create table store_zip3 as 
	select distinct scantrack_market_descr, store_code_uc, year, store_zip3, retailer_code, channel_type
	from mainex.trips(where=(year>2003 and store_code_uc > 0)); 
quit;	

* Check if a store has multiple zip3; 
* There are 883 duplicated observations; 
proc sort data store_zip3 nodupkeys out = tmp; by scantrack_market_descr store_code_uc year; run;	

* Check the number of households from each zipcode; 
proc univariate data = store_allzip; var num_hh; histogram; run;

* Check the number of households and number of zipcodes for each store; 
proc sql noprint; 
	create table tmp as 
	select store_code_uc, year, count(unique(panelist_zip_code)) as n_zip, sum(num_hh) as n_hh
	from store_allzip
	group by store_code_uc, year; 
quit; 
proc univariate data = tmp; var n_zip n_hh; histogram; run;

* Find the zipcode with the largest consumer size; 
proc sql noprint; 
	create table store_zipcode as 
	select scantrack_market_descr, store_code_uc, panelist_zip_code as zipcode
	from store_allzip
	group by scantrack_market_descr, store_code_uc
	having nhh_allyear = max(nhh_allyear);
quit; 

************************************************;
* Calculate the central longitude and latitude *; 
************************************************;
* Need to convert from radian to degree, multiply by pi/180; 
%let pi180=0.0174532925199433;

* A function to find the centroid of stores around neighboring zipcodes; 
%macro zipcenter; 
	proc sql noprint; select max(id) into: n from store_zipmulti; quit; 
	%put &n; 
	
	* Loop over unique store id; 	
	%do i=1 %to &n; 
		* Subset the zipcodes around a store; 
		proc sql noprint; 
			create table tmpz(keep = x y) as 
			select * from store_zipmulti where id = &i; 
			
			select x, y into :initx, :inity from tmpz where monotonic() = 1; 
		quit; 	
		
		* Set initial values for the longitude(x0) and latitude(y0); 
		* NOTE: initial value of coordinates can not be exactly the same with surrounding coordinates, so I use floor; 
		data par(type=est);
		   keep _type_ x0 y0;
		   _type_='parms'; x0 = floor(&initx); y0=floor(&inity); output; 
		run;
		
		* Solve the coordinates that minimize the sum squared Haversine distance to all the surounding zip codes; 
		/** Reference: http://www.movable-type.co.uk/scripts/latlong.html; */
		proc nlp data=tmpz inest = par outest=opar noprint;
		   lsq d;
		   parms x0 y0;
		   a = sin((y0*&pi180-y*&pi180)/2)**2 + cos(y*&pi180)*cos(y0*&pi180)*sin((x*&pi180-x0*&pi180)/2)**2;
		   b =  2*atan2(sqrt(a), sqrt(1-a));
		   d =  6371*b; 
		run;
		
		* Stack the solution; 
		%if &i=1 %then %do; 
			data outd(keep = id x0 y0); set opar(where = (_type_="PARMS")); id = &i; run; 
			proc datasets noprint; delete opar; run;
		%end; 
		%else %do; 
			data tmp1(keep = id x0 y0); set opar(where = (_type_="PARMS")); id = &i; run; 
			data outd; set outd tmp1; run; 
			proc datasets noprint; delete opar tmp1; run;
		%end; 
	%end; 
	proc datasets noprint; delete tmpz par; run;
%mend zipcenter; 

* Merge zipcode with latitude and logitude; 
proc sort data = store_allzip nodupkeys out = tmp(keep = scantrack_market_descr store_code_uc panelist_zip_code); 
	by scantrack_market_descr store_code_uc panelist_zip_code; 
run;

proc sql noprint; 
	create table store_zipll as
	select *, count(unique(panelist_zip_code)) as nzip
	from (select A.*, B.X, B.Y, B.AREACODE from tmp as A left join SASHELP.ZIPCODE as B
		on A.panelist_zip_code = B.zip
		)
	group by scantrack_market_descr, store_code_uc
	order by scantrack_market_descr, store_code_uc; 
quit; 
proc sql; select sum(X=.) from store_zipll; quit; 

* Solve the minimization algorithm only for the stores with multiple surrounding zipcodes; 	
data store_zipmulti; 
	set store_zipll(where = (nzip>1)); 
	by scantrack_market_descr store_code_uc;
	retain id 0;
	if first.store_code_uc then id = id+1; 
run; 

* Run the macro; 
%zipcenter; 

* Merge the coordinate data; 
proc sql noprint; 
	create table tmp1(drop = id) as 
	select A.*, x0 as longitude, y0 as latitude 
	from (select distinct scantrack_market_descr, store_code_uc, id from store_zipmulti) as A 
	inner join outd as B on A.id = B.id; 
quit; 

data tmp2(keep = scantrack_market_descr store_code_uc longitude latitude); 
	set store_zipll(where = (nzip=1)); 
	longitude = x;
	latitude = y; 
run;

data store_coord; 
	set tmp1 tmp2; 
run;
proc datasets noprint; delete tmp tmp1 tmp2 store_zipmulti store_zipll; run;

* Merge zipcodes and coordinates; 
proc sql noprint; 
	create table tmp1 as 
	select A.*, B.longitude, B.latitude
	from store_zipcode as A left join store_coord as B
	on A.scantrack_market_descr = B.scantrack_market_descr and A.store_code_uc = B.store_code_uc
	order scantrack_market_descr, store_code_uc; 
	
	create table store_location as 
	select A.*, B.zipcode, B.longitude, B.latitude
	from store_zip3 as A full join tmp1 as B
	on A.scantrack_market_descr = B.scantrack_market_descr and A.store_code_uc = B.store_code_uc 
	order scantrack_market_descr, store_code_uc; 
quit; 

* Check any number of missing coordinates; 
proc sql; select nmiss(latitude)/count(*) from store_location; run;

data mylib.store_location; set store_location; run;




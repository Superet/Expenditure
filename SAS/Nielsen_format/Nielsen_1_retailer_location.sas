/*

This script imputes the store location from household trip data
Output library: mylib
Output data: retailer_location

Assumption: there is a "representative" store in a 3-digit zipcode area for a retailer.
Retailer identifier: zip3-retailer_code

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

**********************************;
* List zipcode around each retailer *; 
**********************************;
* Find the customer size from each zip code area that is around a given retailer ;
proc sql noprint;
	/*First merge trip data with household zipcode, then find all the zipcodes around a retailer*/
	create table tmp as 	
	select scantrack_market_descr, zip3, retailer_code, channel_type, year, panelist_zip_code, 
		count(distinct household_code) as num_hh, sum(projection_factor) as projection
	from (select A.*, B.panelist_zip_code, substr(put(panelist_zip_code, z5.), 1, 3) as zip3, projection_factor
		  from mainex.trips(where=(year>2003)) as A 
		  left join mainex.panelists(where=(panelist_zip_code^=.)) as B
		  on A.household_code = B.household_code and A.year=B.panel_year) 
	group by scantrack_market_descr, zip3, retailer_code, channel_type, year, panelist_zip_code
	order by scantrack_market_descr, zip3, retailer_code, channel_type, year, num_hh;
	
	* count the number of households in each zipcode area; 
	create table retailer_allzip as 
	select *, sum(num_hh) as nhh_allobs, sum(projection) as tot_projection
	from tmp(where = (panelist_zip_code^=.)) 
	group by scantrack_market_descr, zip3, retailer_code, panelist_zip_code
	order by scantrack_market_descr, zip3, retailer_code, nhh_allobs, panelist_zip_code; 
	
	* Extract 3-digit zip code - retailer code in the data; 
	create table retailer_zip3 as 
	select distinct scantrack_market_descr, zip3, retailer_code, year, channel_type
	from tmp
	order by scantrack_market_descr, zip3, retailer_code, year, channel_type; 
quit;	

* Check the number of households from each zipcode; 
proc univariate data = retailer_allzip; var num_hh; histogram; run;

* Check the number of households and number of zipcodes for each retailer; 
proc sql noprint; 
	create table tmp as 
	select zip3, retailer_code, year, count(unique(panelist_zip_code)) as n_zip, sum(num_hh) as n_hh
	from retailer_allzip
	group by zip3, retailer_code, year; 
quit; 
proc univariate data = tmp; var n_zip n_hh; histogram; run;

* Find the zipcode with the largest consumer size; 
proc sql noprint; 
	create table retailer_zipcode as 
	select scantrack_market_descr, zip3, retailer_code, panelist_zip_code as zipcode
	from (select distinct scantrack_market_descr, zip3, retailer_code, panelist_zip_code, nhh_allobs, tot_projection from retailer_allzip)
	group by scantrack_market_descr, zip3, retailer_code
	/*having nhh_allobs = max(nhh_allobs);*/
	having tot_projection = max(tot_projection);
quit; 

************************************************;
* Calculate the central longitude and latitude *; 
************************************************;
* Need to convert from radian to degree, multiply by pi/180; 
%let pi180=0.0174532925199433;

* A function to find the centroid of retailers around neighboring zipcodes; 
%macro zipcenter; 
	proc sql noprint; select max(id) into: n from retailer_zipmulti; quit; 
	%put &n; 
	
	* Loop over unique retailer id; 	
	%do i=1 %to &n; 
		* Subset the zipcodes around a retailer; 
		proc sql noprint; 
			create table tmpz(keep = x y w) as 
			select * from retailer_zipmulti where id = &i; 
			
			select x, y into :initx, :inity from tmpz where monotonic() = 1; 
		quit; 	
		
		* Set initial values for the longitude(x0) and latitude(y0); 
		* NOTE: initial value of coordinates can not be exactly the same with surrounding coordinates, so I use floor; 
		data par(type=est);
		   keep _type_ x0 y0;
		   _type_='parms'; x0 = floor(&initx)-0.1; y0=floor(&inity)-0.1; output; 
		run;
		
		* Solve the coordinates that minimize the weighted sum squared Haversine distance to all the surounding zip codes; 
		/** Reference: http://www.movable-type.co.uk/scripts/latlong.html; */
		proc nlp data=tmpz inest = par outest=opar noprint;
		   lsq d;
		   parms x0 y0;
		   a = sin((y0*&pi180-y*&pi180)/2)**2 + cos(y*&pi180)*cos(y0*&pi180)*sin((x*&pi180-x0*&pi180)/2)**2;
		   b =  2*atan2(sqrt(a), sqrt(1-a));
		   d =  6371*b*sqrt(w); 
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
/*proc sort data = retailer_allzip nodupkeys out = tmp(keep = scantrack_market_descr zip3 retailer_code panelist_zip_code nhh_allobs); 
	by scantrack_market_descr zip3 retailer_code panelist_zip_code nhh_allobs; 
run;*/
proc sort data = retailer_allzip nodupkeys out = tmp(keep = scantrack_market_descr zip3 retailer_code panelist_zip_code tot_projection); 
	by scantrack_market_descr zip3 retailer_code panelist_zip_code tot_projection; 
run;

proc sql noprint; 
	create table retailer_zipll as
	select *, count(unique(panelist_zip_code)) as nzip, 
			/*nhh_allobs/sum(nhh_allobs) as w*/
			tot_projection/sum(tot_projection) as w
	from (select A.*, B.X, B.Y from tmp as A left join SASHELP.ZIPCODE as B
		on A.panelist_zip_code = B.zip
		where x^=. 
		)
	group by scantrack_market_descr, zip3, retailer_code
	order by scantrack_market_descr, zip3, retailer_code; 
quit; 
proc sql; select sum(X=.) from retailer_zipll; quit; 

* Solve the minimization algorithm only for the retailers with multiple surrounding zipcodes; 	
data retailer_zipmulti; 
	set retailer_zipll(where = (nzip>1 and X^=.)); 
	by scantrack_market_descr zip3 retailer_code;
	retain id 0;
	if first.retailer_code then id = id+1; 
run; 
proc sql; select Nmiss(w) from retailer_zipmulti; quit; 

proc printto log="U:\Users\ccv103\Desktop\templog.log"; run;
* Run the macro; 
%zipcenter; 

* Merge the coordinate data; 
proc sql noprint; 
	create table tmp1(drop = id) as 
	select A.*, x0 as longitude, y0 as latitude 
	from (select distinct scantrack_market_descr, zip3, retailer_code, id from retailer_zipmulti) as A 
	inner join outd as B on A.id = B.id; 
quit; 

data tmp2(keep = scantrack_market_descr zip3 retailer_code longitude latitude); 
	set retailer_zipll(where = (nzip=1)); 
	longitude = x;
	latitude = y; 
run;

data retailer_coord; 
	set tmp1 tmp2; 
run;
proc datasets noprint; delete tmp tmp1 tmp2 retailer_zipmulti retailer_zipll; run;

* Merge zipcodes and coordinates; 
proc sql noprint; 
	create table tmp1 as 
	select A.*, B.longitude, B.latitude
	from retailer_zipcode as A left join retailer_coord as B
	on A.scantrack_market_descr = B.scantrack_market_descr and A.zip3 = B.zip3 and A.retailer_code = B.retailer_code
	order scantrack_market_descr, retailer_code; 
	
	create table retailer_location as 
	select A.*, B.zipcode, B.longitude, B.latitude
	from retailer_zip3 as A left join tmp1 as B
	on A.scantrack_market_descr = B.scantrack_market_descr and A.zip3 = B.zip3
	and A.retailer_code = B.retailer_code 
	order scantrack_market_descr, zip3, retailer_code, year; 
quit; 

* Check any number of missing coordinates; 
proc sql; select nmiss(latitude)/count(*) from retailer_location; run;

data mylib.retailer_location; set retailer_location; run;


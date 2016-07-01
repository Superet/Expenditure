libname orig "U:\Users\ccv103\Documents\Research\Store switching\SAS_temp\Orig";

* Check the number of households in each zip code; 
proc sql noprint; 
	create table tmp as 
	select zipcode, count(unique(panid)) as num_hh
	from orig.demo
	group by zipcode; 
quit; 
proc univariate data = tmp; var num_hh; histogram; run;

* There are only 6 households in a zip code, not enough to impute penetration; 
* Consider the 3-digit zipcode; 
data tmp; 
	set tmp; 
	zip3 = substr(put(zipcode,z5.), 1, 3); 
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
	select A.*, b.format
	from (select panid, year, chain, count(unique(date)) as num_trip from orig.trip_sum group by panid, year, chain) as A
	left join orig.storefm as B
	on A.chain = B.chain; 
	
	create table tmp1 as 
	select A.*, B.zipcode
	from (select panid, format, sum(num_trip) as num_trip from tmp group by panid, format) as A 
	left join orig.demo as B
	on A.panid = B.panid
quit; 

data tmp1; 
	set tmp1; 
	zip3 = substr(put(zipcode,z5.), 1, 3);
run;

* Compute the channel penetration for each zip3; 
proc sql noprint;
	create table tmp as 
	select zip3, format, sum(num_trip) as num_trip
	from tmp1
	group by zip3, format;
	
	create table pen as 
	select zip3, format, num_trip, num_trip/sum(num_trip) as penetration
	from tmp 
	group by zip3; 
quit; 

* Fill 0 to the channels without trip observations; 
proc transpose data = pen out = pen_wide; 
	by zip3; id format; var penetration; 
run;
data pen_wide; 
	set pen_wide; 
	array alln {*} _numeric_; 
	do i = 1 to dim(alln); 
		if alln[i] = . then alln[i] = 0; 
	end;
	drop i; 
run;

proc transpose data = pen_wide out = pen_long(rename=(_name_ = channel_type)); 
	by zip3; 
run;

* Histogram of penetration of wholesale clubs; 
proc univariate data = pen_long(where=(channel_type = "Wholesale_clubs")); var penetration; histogram; run;

* Export the data; 
PROC EXPORT DATA= pen_long
            OUTFILE= "U:\Users\ccv103\Documents\Research\Store switching\iri_channel_penetration.csv" 
            DBMS=csv REPLACE;
     PUTNAMES=YES;
RUN;



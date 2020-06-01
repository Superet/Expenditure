/*
This script identifies the retialer_code for major retailer in Chicago area.
The retailers under focus are: 
	Drug stores: CVS, Walgreen's, 
	Discount stores: WalMart, Target
	Dollar stores: Dollar General, Dollar tree, Family Dollar
	Grocery stores: Jewel, Dominick's, Mariano's, Treasure Island
	Wholesale clubs: Costco, Sam's
	
From other analysis, we have known that the retailer for CVS and Dominick's are	. 
We determine the retailer codes for other as follows: 
1. Within drug stores, the other top drug store is Walgreen's. 
2. Within discount stores, the top two retailers are Walmart and Target, respectively. 
		(We also verify this with national retailer revenue.)
3. Within dollar stores: 
4. Within grocery stores, the top grocery retailer with the largest revenue is Jewel, 
		the one that enters in 2014 is Mariano's. 
		The third one to which households in Chicago often go to is Treasure Island. 
5. Within wholesale clubs, the top 2 retailers are Costco and Sam's, respectively. 
		(We also verify this with national retailer revenue.)
*/

libname mylib "U:\Users\ccv103\Documents\Research\Store switching\SAS_temp\Nielsen";
libname iri "U:\Users\ccv103\Documents\Research\Store switching\SAS_temp\Orig"; 
%let dirname= U:\Users\ccv103\Documents\Research\Nielsen\nielsen_extracts\HMS;
options compress=yes reuse=yes; 

* Define the retailer scope; 
data my_channel;
	input channel_type $40.	channel_group $10.;
	datalines;
	Grocery								Grocery
	Drug Store							Grocery
	Dollar Store						Grocery
	Convenience Store					Grocery
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

* Import trip data;
%let fpath = &dirname\2013\Annual_Files\trips_2013.tsv;
%put &fpath;
PROC IMPORT OUT= WORK.trips2013
            DATAFILE= "&fpath" 
            DBMS=dlm REPLACE;
			DELIMITER ='09'x;
     		DATAROW=2; 
RUN;

%let fpath = &dirname\2014\Annual_Files\trips_2014.tsv;
%put &fpath;
PROC IMPORT OUT= WORK.trips2014
            DATAFILE= "&fpath" 
            DBMS=dlm REPLACE;
			DELIMITER ='09'x;
     		DATAROW=2; 
RUN;

data trips; 
	set trips2013 trips2014; 
	year = year(purchase_date); 
	if year = 2012 then delete; 
run;
proc datasets noprint; delete trips2013 trips2014; run;

%let fpath = &dirname\2013\Annual_Files\panelists_2013.tsv;
%put &fpath;
PROC IMPORT OUT= WORK.panelists
            DATAFILE= "&fpath" 
            DBMS=dlm REPLACE;
			DELIMITER ='09'x;
     		DATAROW=2; 
			GUESSINGROWS=5000;
RUN;

proc sql noprint; 
	* Subset trip data to Chicago; 
	create table trips_chicago  as 
	select *
	from trips 
	where household_code in 
	(select household_code from panelists where scantrack_market_descr = "Chicago");
	
	* Subset households to Chicago; 
	create table pan_chicago as 
	select A.*, B.city, B.city2, countynm
	from (	select * 
			from panelists where scantrack_market_descr = "Chicago") as A
		left join sashelp.zipcode as B
		on A.panelist_zip_code = B.zip;
	
	* Sum up retailer revenue in 2013 and 2014 for Chicago market; 
	create table retailer_rev as 
	select * 
	from (	select A.*, B.channel_type
			from (select year, retailer_code, sum(total_spent) as revenue
			from trips_chicago
			group by retailer_code, year) as A 
			left join 
			orig_retailers as B
			on A.retailer_code = B.retailer_code
			)
	where channel_type in (select channel_type from my_channel)
	order by year, channel_type, revenue descending	;	
	
	* Sum up retailer revenue in 2013 and 2014 for national market; 
	create table retailer_rev_nation as 
	select * 
	from (	select A.*, B.channel_type
			from (select year, retailer_code, sum(total_spent) as revenue
			from trips
			group by retailer_code, year) as A 
			left join 
			orig_retailers as B
			on A.retailer_code = B.retailer_code)
	where channel_type in (select channel_type from my_channel) 	
	order by year, channel_type, revenue descending;
quit;

* Subset top retailer for each channel; 
proc rank data = retailer_rev_nation descending out = retailer_rev_rk;
	by year channel_type; 
	var revenue;
	ranks rev_rank; 
run;

data retailer_rev_nation;
	set retailer_rev_rk;
	if rev_rank > 20 then delete; 
run;

proc rank data = retailer_rev descending out = retailer_rev_rk;
	by year channel_type; 
	var revenue;
	ranks rev_rank; 
run;

data retailer_rev;
	set retailer_rev_rk;
	if rev_rank > 20 then delete; 
run;

* --------------------------------*; 
* Summary of chain revenue in IRI *; 
* Mkt = 10 for Chicago; 
proc sql noprint;
	create table tmp as 
	select A.*, B.chain_name, B.format
	from (	select chain, sum(BSKTDOL) as revenue
			from iri.trip_sum
			where year = 2009
			group by chain) as A 
	left join iri.storefm as B
	on A.chain = B.chain
	order by format, revenue descending; 
	
	create table tmp1 as 
	select A.*, B.chain_name, B.format
	from (	select chain, sum(BSKTDOL) as revenue
			from iri.trip_sum
			where year = 2009 and panid in (select panid from iri.demo where mkt = 10)
			group by chain) as A 
	left join iri.storefm as B
	on A.chain = B.chain
	order by format, revenue descending;	
quit;

proc rank data = tmp descending out = iri_rev_nation; 
	by format; 
	var revenue;
	ranks rev_rank;
run;

proc rank data = tmp1 descending out = iri_rev; 
	by format; 
	var revenue;
	ranks rev_rank;
run;
proc datasets noprint; delete tmp tmp1; run;

*---------------*; 
* Stores in RMS *; 
data orig_stores; 
	infile datalines; 
	length file2read $300; 
	input file2read $;
	infile dummy filevar=file2read delimiter='09'x MISSOVER DSD LRECL=32767 FIRSTOBS=2 end = done;
	do while(not done);
		informat channel_code $1. fips_state_descr $2. fips_county_descr $50. dma_descr $50.;
		input 	store_code_uc year parent_code retailer_code channel_code $ store_zip3 fips_state_code fips_state_descr $
				fips_county_code fips_county_descr $ dma_code dma_descr $;
		output; 
	end;
	datalines;
	U:\Users\ccv103\Documents\Research\Nielsen\nielsen_extracts\RMS\2013\Annual_Files\stores_2013.tsv
	U:\Users\ccv103\Documents\Research\Nielsen\nielsen_extracts\RMS\2014\Annual_Files\stores_2014.tsv
	;
run;

***************************;
* Graph retailer revenue * ;
***************************;
*-------------------------------------------*; 
* Drug stores: Walgreens = 4904, CVS = 4914 *; 
* Osco drug = 4908, Rite Aid = 4901; 
* RMS tracks 4904; 
proc gchart data = retailer_rev(where = (channel_type = "Drug Store")); 
	vbar year/discrete sumvar = revenue group = retailer_code; 
	title "Annual revenue for drug stores in Chicago"; 
run;

proc gchart data = retailer_rev_nation(where = (channel_type = "Drug Store")); 
	vbar year/discrete sumvar = revenue group = retailer_code; 
	title "Annual revenue for drug stores in the U.S."; 
run;

data tmp; 
	set orig_stores; 
	where retailer_code in (4904, 4914); 
run;
proc freq data = tmp; table retailer_code; run;

*-----------------------------------------------------------------*; 
* Discount stores: WalMart = 6920 (super stores) and 6905 (regular), ; 
*		Target = 6901 (regular) and 6921 (super store)*;	
* RMS tracks 6901 and 6921; 
proc gchart data = retailer_rev(where = (channel_type = "Discount Store")); 
	vbar year/discrete sumvar = revenue group = retailer_code; 
	title "Annual revenue for discount stores in Chicago"; 
run;

proc gchart data = iri_rev(where = (format = "Mass")); 
	vbar chain_name/discrete sumvar = revenue ; 
	title "Annual revenue for discount stores from IRI in Chicago "; 
run;

proc gchart data = retailer_rev_nation(where = (channel_type = "Discount Store")); 
	vbar year/discrete sumvar = revenue group = retailer_code; 
	title "Annual revenue for discount stores in the U.S."; 
run;

proc gchart data = iri_rev_nation(where = (format = "Mass")); 
	vbar chain_name/discrete sumvar = revenue ; 
	title "Annual revenue for discount stores from IRI in U.S."; 
run;

*------------------------------------------------------------------------------------*; 
* Dollar stores: Dollar General = 5850, Dollar tree = 5853, Family dollar = 5851	*;	
* RMS tracks 5851; 
proc gchart data = retailer_rev(where = (channel_type = "Dollar Store")); 
	vbar year/discrete sumvar = revenue group = retailer_code; 
	title "Annual revenue for dollar stores in Chicago"; 
run;

proc gchart data = iri_rev(where = (format = "Dollar/variety")); 
	vbar chain_name/discrete descending sumvar = revenue; 
	title "Annual revenue for dollar stores from IRI in Chicago"; 
run;

proc gchart data = retailer_rev_nation(where = (channel_type = "Dollar Store")); 
	vbar year/discrete sumvar = revenue group = retailer_code; 
	title "Annual revenue for dollar stores in the U.S."; 
run;

proc gchart data = iri_rev_nation(where = (format = "Dollar/variety")); 
	vbar chain_name/discrete descending sumvar = revenue; 
	title "Annual revenue for dollar stores from IRI in Chicago"; 
run;

data tmp; 
	set orig_stores; 
	where retailer_code in (5850, 5851, 5852, 5853); 
run;
proc freq data = tmp; table retailer_code; run;

*-----------------------------------------------*; 
* Grocery stores: Jewel = 6, Dominicks = 61, *;	
* Others, Aldi = 24, Food 4 less = 858, Ultra foods = 770, Meijer = 151, A/O 2MM+ GROCER = 3997; 
* RMS tracks 6, 61, 770, 858, 889, 3997; 
proc gchart data = retailer_rev(where = (channel_type = "Grocery")); 
	vbar year/discrete sumvar = revenue group = retailer_code; 
	title "Annual revenue for grocery stores in Chicago"; 
run;

proc gchart data = iri_rev(where = (format in ("Conventional supermarkets", "Super stores"))); 
	vbar chain_name/discrete descending sumvar = revenue; 
	title "Annual revenue for grocery stores from IRI in Chicago"; 
run;

proc gchart data = retailer_rev_nation(where = (channel_type = "Grocery")); 
	vbar year/discrete sumvar = revenue group = retailer_code; 
	title "Annual revenue for grocery stores in the U.S."; 
run;

proc gchart data = iri_rev_nation(where = (format in ("Conventional supermarkets", "Super stores"))); 
	vbar chain_name/discrete descending sumvar = revenue; 
	title "Annual revenue for grocery stores from IRI in U.S."; 
run;

data tmp; 
	set orig_stores; 
	where retailer_code in (6, 61, 151, 24, 296, 455, 770, 858, 889, 916, 3997, 3999); 
run;
proc freq data = tmp; table retailer_code; run;

*--------------------------------------------------------*; 
* Warehouse clubs: Costco = 9103, Sam's = 9101, BJ = 9104;	
* RMS tracks 9104;
proc gchart data = retailer_rev(where = (channel_type = "Warehouse Club")); 
	vbar year/discrete sumvar = revenue group = retailer_code; 
	title "Annual revenue for warehouse clubs in Chicago"; 
run;

proc gchart data = retailer_rev_nation(where = (channel_type = "Warehouse Club")); 
	vbar year/discrete sumvar = revenue group = retailer_code; 
	title "Annual revenue for warehouse clubs in the U.S."; 
run;

proc gchart data = iri_rev_nation(where = (format = "Wholesale clubs")); 
	vbar chain_name/discrete descending sumvar = revenue; 
	title "Annual revenue for grocery stores from IRI in U.S."; 
run;

data tmp; 
	set orig_stores; 
	where retailer_code in (9103, 9101, 9104); 
run;
proc freq data = tmp; table retailer_code; run;

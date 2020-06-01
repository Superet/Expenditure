/*
This file generates the geographic market for each store;

Store data in RMS is at store-year level. It contains zip3, fips code, and dma, but not scantrack market code. 
To make sure it matches with the scantrack in HMS, we merge together RMS and HMS store data;

Channel type has different classification in HMS and RMS.
*/

libname mainex "E:\Users\ccv103\Documents\Research\Store switching\SAS_temp\Nielsen\Work\main";

*--------------------------*; 
* Read master retailer data ;
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
	E:\Users\ccv103\Documents\Research\Nielsen\nielsen_extracts\RMS\2006\Annual_Files\stores_2006.tsv
	E:\Users\ccv103\Documents\Research\Nielsen\nielsen_extracts\RMS\2007\Annual_Files\stores_2007.tsv
	E:\Users\ccv103\Documents\Research\Nielsen\nielsen_extracts\RMS\2008\Annual_Files\stores_2008.tsv
	E:\Users\ccv103\Documents\Research\Nielsen\nielsen_extracts\RMS\2009\Annual_Files\stores_2009.tsv
	E:\Users\ccv103\Documents\Research\Nielsen\nielsen_extracts\RMS\2010\Annual_Files\stores_2010.tsv
	E:\Users\ccv103\Documents\Research\Nielsen\nielsen_extracts\RMS\2011\Annual_Files\stores_2011.tsv
	E:\Users\ccv103\Documents\Research\Nielsen\nielsen_extracts\RMS\2012\Annual_Files\stores_2012.tsv
	E:\Users\ccv103\Documents\Research\Nielsen\nielsen_extracts\RMS\2013\Annual_Files\stores_2013.tsv
	;
run;

*-----------------------------*;
* Obstain store code from HMS *;
* NOTE: store_code_uc show up only in trip data, we need to match the fips code using household as a link; 

* Read panelist data;
%let dirname= E:\Users\ccv103\Documents\Research\Nielsen\nielsen_extracts\HMS;
%macro read_panelists;
	%do i=2004 %to 2012;
		* Read panelists data for geographic market information; 
		%let fpath = &dirname\&i\Annual_Files\panelists_&i..tsv;
		%put &fpath;
		PROC IMPORT OUT= WORK.panelists&i
		            DATAFILE= "&fpath" 
		            DBMS=dlm REPLACE;
					DELIMITER ='09'x;
		     		DATAROW=2; 
					GUESSINGROWS=5000;
		RUN;
		
		data panelists&i;
			length scantrack_market_descr $50.;
			set panelists&i(keep = household_code panel_year dma_code dma_descr fips_state_code fips_county_code scantrack_market_descr); 
		run;
	%end; 
	data panelists;
		set panelists2004 - panelists2012;
	run;
	proc datasets noprint; 
		delete panelists2004 - panelists2012;
	run;
%mend read_panelists;

%read_panelists;

* Examine if scantrack_market_code corresponds to one unique fips code -- Yes; 
proc sql noprint;
	create table tmp as 
	select fips_state_code, fips_county_code, count(unique(scantrack_market_descr)) as n_market
	from panelists
	group by fips_state_code, fips_county_code;
quit;
proc freq data = tmp; table n_market;run;

* Merge fips code into the trip data; 
proc sql noprint;
	create table joint_trip as
	select A.*, fips_state_code, fips_county_code
	from mainex.trips(keep = household_code store_code_uc panel_year channel_type scantrack_market_descr)
	as A left join panelists as B
	on A.household_code = B.household_code;
	
	* Determine the fips of each store-year; 
	create table tmp as 
	select store_code_uc, panel_year, fips_state_code, fips_county_code, count(unique(household_code)) as nhh
	from joint_trip(keep = store_code_uc panel_year fips_state_code fips_county_code household_code)
	group by store_code_uc, panel_year, fips_state_code, fips_county_code
	order by store_code_uc, panel_year, nhh; 
quit;

* Determin the fips code for each store in HMS, by matching fips code with RMS; 
proc sql noprint;
	create table tmp1 as 
	select A.*
	from tmp as A, orig_stores as B
	where A.store_code_uc = B.store_code_uc and A.panel_year = B.year 
		and A.fips_state_code = B.fips_state_code and A.fips_county_code = B.fips_county_code;
quit; 
* Check if store has a unique fips code -- Yes; 
proc sort data = tmp1 nodupkeys out = tmp2; by store_code_uc panel_year; run;

* Subset the unique pair of store_code_uc and scantrack_market_uc that can match fips code; 
proc sql noprint; 
	create table store_hms as 
	select distinct *
	from (select A.*
		from joint_trip(keep = store_code_uc panel_year channel_type scantrack_market_descr fips_state_code fips_county_code) as A, 
		tmp1 as B
		where A.store_code_uc = B.store_code_uc and A.panel_year = B.panel_year 
			and A.fips_state_code = B.fips_state_code and A.fips_county_code = B.fips_county_code)
	order by store_code_uc, panel_year;
	
	drop table joint_trip; 
quit; 

* Check if channel_type, and market is unique for a store-year -- No; 	
proc sql noprint;
	create table tmp as 
	select store_code_uc, panel_year, 
			count(unique(channel_type)) as n_channel, 
			count(unique(scantrack_market_descr)) as n_market
	from store_hms
	group by store_code_uc, panel_year;
quit; 
proc means data = tmp; var n_channel n_market; run;

* For the stores that have multiple scantrack market, we set the scantrack_market as the one that show up the most; 
proc sql noprint; 
	create table tmp1 as 
	select *
	from store_hms 
	where store_code_uc in (select store_code_uc from tmp where n_market > 1);
	
	create table tmp2 as
	select store_code_uc, scantrack_market_descr, count(panel_year) as n
	from tmp1 
	group by store_code_uc, scantrack_market_descr
	order by store_code_uc, n; 
quit;

* Select the scantrack_market_uc that show up most year; 
data tmp2; 	
	set tmp2;
	by store_code_uc n; 
	if not last.store_code_uc then delete; 
run;

data tmp3(drop = scantrack); 
	merge tmp1 tmp2(drop = n rename=(scantrack_market_descr=scantrack)); 
	by store_code_uc; 
	if scantrack_market_descr^=scantrack then delete; 
proc sort; by store_code_uc panel_year; 
run;

* Re-stack the stores that have unique scantrack_market_descr;
proc sql noprint; 
	create table tmp2 as 
	select * from store_hms where store_code_uc not in (select store_code_uc from tmp where n_market>1);
quit; 
data store_hms;
	set tmp2 tmp3; 
proc sort; by store_code_uc panel_year; 
run;

* Re-check the store number; 
* NOTE: we lose 34 store-year pair; 
proc contents data = tmp; run;
proc sql noprint;
	create table tmp as 
	select store_code_uc, panel_year, 
			count(unique(channel_type)) as n_channel, 
			count(unique(scantrack_market_descr)) as n_market
	from store_hms
	group by store_code_uc, panel_year;
quit; 
proc means data = tmp; var n_channel n_market; run;

* Match together the data; 
proc sql noprint;
	create table my_stores as 
	select A.*, B.scantrack_market_descr, B.channel_type
	from orig_stores as A
	left join store_hms as B
	on A.store_code_uc = B.store_code_uc and A.year = B.panel_year
	order by store_code_uc, year;
quit; 

* NOTE: there are many stores that do not show up in HMS, we fill in the missing channel_type using channel_code;
proc freq data = my_stores; table channel_code*channel_type; run;
data my_stores; 
	set my_stores; 
	if channel_type = "" then do; 
		if channel_code = "C" then channel_type = "Convenience Store"; 
		if channel_code = "D" then channel_type = "Drug Store";
		if channel_code = "F" then channel_type = "Grocery"; 
		if channel_code = "M" then channel_type = "Discount Store"; 
	end; 
run;

* Export matched store data; 
PROC EXPORT DATA= my_stores
            OUTFILE= "E:\Users\ccv103\Documents\Research\Nielsen\RMS_attributes\Master\stores_both.csv"
            DBMS=csv REPLACE;
     PUTNAMES=YES;
RUN;


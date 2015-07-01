%let dirname= E:\Users\ccv103\Documents\Research\Nielsen\nielsen_extracts\RMS;
libname mylib "E:\Users\ccv103\Documents\Research\Store switching\SAS_temp\Nielsen";

*************;
* Read data *;
*************;
* Read master retailer data;
data orig_stores; 
	infile datalines; 
	length file2read $300; 
	input file2read $;
	infile dummy filevar=file2read delimiter='09'x MISSOVER DSD LRECL=32767 FIRSTOBS=2 end = done;
	do while(not done);
		informat channel_code $1. fips_state_descr $2. fips_county_descr $12. dma_descr $40.;
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
	;
run;

* The movement name; 
/*
%let fpath = &dirname\2006\Movement_Files\2506_2006\3625_2006.tsv;
%put &fpath;
options obs=100;
PROC IMPORT OUT= WORK.move06
            DATAFILE= "&fpath" 
            DBMS=dlm REPLACE;
			DELIMITER ='09'x;
     		DATAROW=2; 
			GUESSINGROWS=5000;
RUN;
options obs=max; 
proc contents data=move06; run;
*/

data movename; 
	input department module;
	datalines; 
	2506	3625
	;
run;

/* Yogurt
	2510	3603
	
*/

data tmp_file; 
	length filep $300. folder $300.; 
	set movename; 
	do year=2006 to 2012;
		folder = cats("E:\Users\ccv103\Documents\Research\Nielsen\nielsen_extracts\RMS\",year,"\Movement_Files\",department,"_",year);
		filen = cats(module,'_',year,'.tsv');
		filep = cats(folder, "\",filen);
		output;
	end;
run;

data movement;
	set tmp_file(keep=filep);
	infile dummy filevar = filep delimiter='09'x MISSOVER DSD LRECL=32767 FIRSTOBS=2 end = done;
	do while(not done);
		input store_code_uc upc week_end units prmult price feature display; 
		output; 
	end;
	format week_end date9.;
	drop filep;
run;
proc contents data=movement; run;

proc datasets noprint; delete movename tmp_file; run;

************************;
* Store entry and exit *;
************************;
%let mywindow = 90; 		
proc sql noprint;
	create table tmp as 
	select store_code_uc, min(week_end) as first_week, max(week_end) as last_week
	from movement 
	group by store_code_uc
	order by store_code_uc;
	
	create table my_store as
	select A.*,B.first_week, B.last_week
	from orig_stores as A left join tmp as B
	on A.store_code_uc = B.store_code_uc
	order by store_code_uc, year; 
quit;

data my_store; 
	set my_store(rename = (first_week = first_week1 last_week=last_week1));
	first_week = input(put(first_week1, 8.), YYMMDD10.);
	last_week = input(put(last_week1, 8.), YYMMDD10.);
	format first_week YYMMDD10. last_week YYMMDD10.; 
	drop first_week1 last_week1;
run;
proc contents data=my_store; run;

* Check if store records exist when no sales ;
proc sql; 
	select count(year) from my_store where year < year(first_week); 
	select count(year) from my_store where year > year(last_week); 
quit;

data tmp; 
	set my_store;
	where year > year(last_week) and first_week ^=. ; 
run;

data tmp1; 
	set my_store;
	where year < year(first_week) and first_week ^=. ; 
run;

* Unique stores with locations; 
proc sort data=my_store nodupkeys out=prim_store;
	by store_code_uc ;
run;

* Subset the stores that did not have sales until the last week; 
proc sql noprint; 
	select max(last_week) into: week_end from prim_store; 
	select min(first_week) into: week_first from prim_store; 
quit;
%put &week_first; 
%put &week_end;


data prim_store_exit;
	set prim_store(drop = year);
	where last_week < &week_end - &mywindow and first_week^=.; 
proc sort; by fips_state_descr dma_descr store_zip3 store_code_uc; 
run;

data prim_store_entry;
	set prim_store(drop = year);
	where first_week > &week_first + &mywindow and last_week^=.; 
proc sort; by fips_state_descr dma_descr store_zip3 store_code_uc; 
run;

* Summarize the stores over market-year; 
proc sql noprint;
	create table tmp as 
	select dma_descr, year, count(unique(store_code_uc)) as num_store
	from my_store
	group by dma_descr, year
	order by dma_descr, year;
quit;

proc transpose data=tmp out = overall_sum(drop=_NAME_) prefix=num;
	by dma_descr;
	id year;
	var num_store; 
run;

* Summary by chain-market: number of stores, number of entry, number of exists;
proc sql noprint; 
	create table tmp as 
	select parent_code, dma_descr, year, count(distinct(store_code_uc)) as num_store
	from my_store 
	group by parent_code, dma_descr, year
	order by parent_code, dma_descr, year;
quit;

proc transpose data=tmp out=tmp1(drop=_NAME_) prefix = num;
	by parent_code dma_descr;
	id year;
	var num_store;
run;

proc transpose data=prim_store

************************;
* Save and export data *;
************************;
data mylib.rms_stores; set my_store; run;

PROC EXPORT DATA=prim_store_exit
			OUTFILE= "\\tsclient\Resear1\Store switching\processed data\99_exist_stores.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;

PROC EXPORT DATA=prim_store_entry
			OUTFILE= "\\tsclient\Resear1\Store switching\processed data\99_entry_stores.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;




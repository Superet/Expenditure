/*
This script exports datasets to csv files. 
Files to be exported: format_department, panelists, hh_format, hh_format_dpt, hh_month_dpt

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

proc contents data = mylib.format_department;
proc contents data = mylib.hh_format_dpt;
proc contents data = mylib.hh_month_dpt;
proc contents data = mainex.panelists; 
run;

%let minhh = 500; 
%let minprice = 0.5;
%let maxprice = 20;

* Restrict to the markets where there is a large number of households;
proc sql noprint;
	create table tmp as
	select scantrack_market_descr, count(unique(household_code)) as nhh
	from mainex.panelists
	where panel_year = 2008
	group by scantrack_market_descr; 
quit; 
proc univariate data = tmp; var nhh; histogram; run;

proc univariate data = mylib.hh_format;
	var unitprice_paid;
run;

* Subset the households in the large markets; 
proc sql; 
select count(*) from tmp where nhh >= &minhh and scantrack_market_descr ^="" ;
select count(unique(household_code)) from mylib.hh_format;
select count(unique(household_code)) from mylib.hh_month_dpt;
select count(unique(household_code)) from mysqllib.hh_format_retattr where unitprice_paid < &minprice or unitprice_paid > &maxprice;
quit; 

proc sql noprint;
	create table mkt as 
	select scantrack_market_descr from tmp where nhh >= &minhh and scantrack_market_descr ^="" ; 
	
	create table hh1 as 
	select distinct A.household_code
	from  
	(select * from mylib.keephh where
	household_code not in (select household_code from mysqllib.hh_format_retattr where unitprice_paid < &minprice or unitprice_paid > &maxprice)) as A
	inner join
 	(select household_code, scantrack_market_descr 
	from mainex.panelists 
	where scantrack_market_descr in (select scantrack_market_descr from mkt)) as B
	on A.household_code= B.household_code;
quit; 

* Adjust household income; 
proc sql noprint; 
	create table panelists as 
	select household_code, panel_year, scantrack_market_descr, panelist_zip_code, 
			household_income, household_size, age_and_presence_of_children, type_of_residence, 
			female_head_age, female_head_employment, male_head_age, male_head_employment, projection_factor
	from mainex.panelists where scantrack_market_descr in (select scantrack_market_descr from mkt) and 
		household_code in (select household_code from hh1)
	order by household_code, panel_year descending;
quit; 
proc sql; 
	select count(unique(scantrack_market_descr)) from panelists; 
	select count(unique(household_code)) from panelists; 
quit; 

* Adjust income for panelists data: income is reported two years prior to the concurrent date; 
* The number of households decreases from 49865 to 34567 after matching income; 
data panelists(drop=n household_income); 
	set panelists; 
	by household_code; 
	income = lag2(household_income); 
	if first.household_code then n=0;
	n+1; 
	if n<=2 then income = .;
	if income=. then delete; 
run;
proc sql; 
	select count(unique(household_code)) from panelists; 
quit;

proc sql noprint;
	create table hh as 
	select distinct A.household_code
	from hh1 as A inner join panelists as B
	on A.household_code = B.household_code;
	
	drop table hh1; 
	
	* Add MSA and FIPS code to the panel data; 
	create table tmp as 
	select A.*, cats(put(state, z2.), put(county,z3.)) format = $5. as fips, msa
	from panelists as A left join SASHELP.zipcode as B
	on A.panelist_zip_code = B.zip
	order by household_code, panel_year;
quit; 
data panelists; set tmp; run; 

* Subset the files to in the large markets; 
proc sql noprint;
	create table format_department as
	select *
	from mylib.format_department 
	where scantrack_market_descr in (select scantrack_market_descr from mkt); 
	
	create table hh_format as
	select *
	from (select B.*, panelist_zip_code, distance_zipcode, distance_haver
		from 
		(select *, case when channel_type = 'Dollar store' then 'Dollar Store' else channel_type end as new_channel
		from mylib.hh_dist) as A 
		full join 
		(select * from mysqllib.hh_format_retattr where channel_type^="") as B
		on A.household_code = B.household_code and A.year = B.year and A.new_channel = B.channel_type
	) 
	where household_code in (select household_code from hh) and scantrack_market_descr in (select scantrack_market_descr from mkt)
	order by household_code, year, channel_type;
	
	create table hh_format_dpt as
	select *
	from mylib.hh_format_dpt 
	where household_code in (select household_code from hh) and 
		scantrack_market_descr in (select scantrack_market_descr from mkt); 
	
	create table hh_month_dpt as 
	select *
	from mylib.hh_month_dpt 
	where household_code in (select household_code from hh);

	* Check if undesirable scantrack_market_descr is dropped out; 
	create table tmp2 as 
	select A.household_code, A.year, scantrack_market_descr
	from hh_month_dpt as A left join panelists as B
	on A.household_code = B.household_code and A.year = B.panel_year;	
quit; 
proc sql;
	select count(unique(scantrack_market_descr)) from panelists;
	select count(unique(scantrack_market_descr)) from format_department;
	select count(unique(scantrack_market_descr)) from hh_format;
	select count(unique(scantrack_market_descr)) from hh_format_dpt;
	select count(unique(scantrack_market_descr)) from tmp2;

	select count(unique(household_code)) from panelists;
	select count(unique(household_code)) from hh_format;
	select count(unique(household_code)) from hh_format_dpt;
	select count(unique(household_code)) from hh_month_dpt;
	
	drop table tmp2; 
quit;

* Check the distribution of retail attributes and see if there are still big outliers; 
proc contents data = hh_format; run;
proc univariate data = hh_format; 
	var unitprice_paid distance_haver num_module; 
run;


* For monthly purchase data, exclude the first 6 months that are used to calculate basket share; 
data hh_month_dpt(drop = n); 
	set hh_month_dpt; 
	by household_code;
	if first.household_code then n = 0; 
	n+1; 
	if n <=6 then delete;
run;

PROC EXPORT DATA=format_department
			OUTFILE= "U:\Users\ccv103\Desktop\format_department.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;

PROC EXPORT DATA=hh_format
			OUTFILE= "U:\Users\ccv103\Desktop\hh_format.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;

PROC EXPORT DATA=hh_format_dpt
			OUTFILE= "U:\Users\ccv103\Desktop\hh_format_dpt.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;

PROC EXPORT DATA=hh_month_dpt
			OUTFILE= "U:\Users\ccv103\Desktop\hh_month_dpt.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;

PROC EXPORT DATA=panelists
			OUTFILE= "U:\Users\ccv103\Desktop\panelists.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;


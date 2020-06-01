libname mylib "U:\Users\ccv103\Documents\Research\Store switching\SAS_temp\Nielsen";
%let dirname= U:\Users\ccv103\Documents\Research\Nielsen\nielsen_extracts\HMS;
options compress=yes reuse=yes; 

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
	create table trips_chicago  as 
	select *
	from trips 
	where household_code in 
	(select household_code from panelists where scantrack_market_descr = "Chicago");

	create table retailer_rev as 
	select A.*, B.channel_type
	from (select year, retailer_code, sum(total_spent) as revenue
	from trips_chicago
	group by retailer_code, year) as A 
	left join 
	orig_retailers as B
	on A.retailer_code = B.retailer_code
	order by year, channel_type, revenue descending;
quit;

* Compare the changes in revenue for retailers in two year; 
* retailer_code = 61 is Dominicks; 
* Top selling Grocery in 2013 include: 6, 3997, 151, 61, 24, 770, 858, 3999, 889, 916; 
data tmp;
	set retailer_rev; 
	where channel_type = "Grocery"; 
run; 

* Select households; 
proc sql noprint;
	create table tmp as 
	select household_code, year, count(unique(purchase_date)) as purchase_n, sum(retailer_code = 61) as domn_n
	from trips_chicago 
	group by household_code, year; 	
	
	create table tmp1 as 
	select household_code, count(year) as n_yr, sum(purchase_n > 5) as n1, 1*(sum(domn_n)>0) as domnick_ever
	from tmp
	group by household_code; 
quit; 
proc freq data = tmp1; 
table n_yr n1 domnick_ever; 
run;

proc sql noprint;
	create table panelist_chi as 
	select A.*, B.domnick_ever
	from (select * from panelists where scantrack_market_descr = "Chicago") as A
	inner join 
	(select * from tmp1 where n_yr = 2 and n1 = 2) as B
	on A.household_code = B.household_code
	order by household_code; 
	
	create table tmp as 
	select A.*, B.domnick_ever
	from trips_chicago as A
	inner join 
	panelist_chi as B
	on A.household_code = B.household_code;
	
	create table trips_chi_sub as 
	select A.*, B.channel_type
	from tmp as A 
	inner join 
	(select * from orig_retailers where channel_type in 
		('Grocery' 'Discount Store' 'Warehouse Club' 'Drug Store' 'Dollar Store' 'Convenience Store')) as B
	on A.retailer_code = B.retailer_code; 
quit; 

***************;
* Export data *; 
***************;
PROC EXPORT DATA = trips_chi_sub
			OUTFILE= "U:\Users\ccv103\Desktop\domnick_trips.csv"
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;

PROC EXPORT DATA = panelist_chi
			OUTFILE= "U:\Users\ccv103\Desktop\domnick_panelists.csv"
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;

*************************;
* Plot share difference *; 
*************************;
* Aggregate expenditure share at formats; 
proc sql noprint;
	create table tmp as 
	select *, sum(dol) as exp, dol/sum(dol) as share
	from (select year, channel_type, sum(total_spent) as dol 
		from trips_chi_sub
		group by year, channel_type)
	group by year
	order by year, share descending; 
	
	* For the households who have shopped at Domnicks; 
	create table tmp1 as 
	select *, sum(dol) as exp, dol/sum(dol) as share
	from (select year, channel_type, sum(total_spent) as dol 
		from trips_chi_sub
		where domnick_ever = 1
		group by year, channel_type)
	group by year
	order by year, share descending;
quit; 

proc gchart data = tmp; 
	vbar year/discrete sumvar = share group = channel_type; 
	title "Share difference between two years"; 
run;

proc gchart data = tmp; 
	vbar year/discrete sumvar = share group = channel_type; 
	title "Share difference between two years for the Domnicks customers"; 
run;

* Aggregate expenditure within grocery store; 
proc sql noprint;
	create table tmp as 
	select *, sum(dol) as exp, dol/sum(dol) as share
	from (select year, channel_type, retailer_code, sum(total_spent) as dol 
		from trips_chi_sub
		group by year, channel_type, retailer_code)
	group by year
	order by year, channel_type, share descending; 
quit; 

data tmp1; 
	set tmp;
	where channel_type = "Grocery"; 
	by year; 
	if first.year then rank = 1; 
	else rank + 1; 
	if rank > 10 then delete; 
run;

proc gchart data = tmp1; 
	vbar year/discrete sumvar = share group = retailer_code; 
	title "Share difference within top grocery stores"; 
run;	

* The annual spending of average consumer; 
proc sql noprint;
	select count(unique(household_code)) into: n1 from trips_chi_sub; 
	select count(unique(household_code)) into: n2 from trips_chi_sub where domnick_ever = 1; 
	
	create table tmp as 
	select household_code, domnick_ever, year, channel_type, retailer_code, sum(total_spent) as dol 
	from trips_chi_sub
	group by household_code, domnick_ever, year, channel_type, retailer_code;
	
	* Expenditure share at formats in different years; 
	create table tmp1 as 
	select *, sum(dol) as exp, dol/sum(dol) as share
	from (select year, channel_type, sum(dol)/&n1 as dol 
		from tmp
		group by year, channel_type)
	group by year
	order by year, share descending; 
	
	* Expenditure shares at top grocery stores; 
	create table tmp2 as 
	select *, sum(dol) as exp, dol/sum(dol) as share
	from (select year, channel_type, retailer_code,  sum(dol)/&n1 as dol 
		from tmp
		group by year, channel_type, retailer_code)
	group by year
	order by year, channel_type, share descending;
	
	* For the households who have shopped at Domnicks; 
	create table tmp3 as 
	select *, sum(dol) as exp, dol/sum(dol) as share
	from (select year, channel_type, sum(dol)/&n2 as dol
		from tmp
		where domnick_ever = 1
		group by year, channel_type)
	group by year
	order by year, channel_type, share descending;
	
	* For the households who have shopped at Domnicks, share at top grocery stores; 
	create table tmp4 as 
	select *, sum(dol) as exp, dol/sum(dol) as share
	from (select year, channel_type, retailer_code, sum(dol)/&n2 as dol
		from tmp
		where domnick_ever = 1
		group by year, channel_type, retailer_code)
	group by year
	order by year, channel_type, share descending;
quit;

data tmp2; 
	set tmp2;
	where channel_type = "Grocery"; 
	by year; 
	if first.year then rank = 1; 
	else rank + 1; 
	if rank > 10 then delete; 
run;

data tmp4; 
	set tmp4;
	where channel_type = "Grocery"; 
	by year; 
	if first.year then rank = 1; 
	else rank + 1; 
	if rank > 10 then delete; 
run;

* Expenditure shift for an average consumer; 
proc gchart data = tmp2; 
	vbar year/discrete sumvar = dol group = retailer_code; 
	title "Dollar difference within top grocery stores"; 
run;
proc gchart data = tmp1; 
	vbar year/discrete sumvar = dol group = channel_type; 
	title "Dollar difference across retail formats"; 
run;
proc gchart data = tmp2; 
	vbar year/discrete sumvar = share group = retailer_code; 
	title "Share difference within top grocery stores"; 
run;
proc gchart data = tmp1; 
	vbar year/discrete sumvar = share group = channel_type; 
	title "Share difference across retail formats"; 
run;

* Expenditure shift for a Domnimck consumer; 
proc gchart data = tmp4; 
	vbar year/discrete sumvar = dol group = retailer_code; 
	title "Dollar difference within top grocery stores for an average Domnicks consumer"; 
run;
proc gchart data = tmp3; 
	vbar year/discrete sumvar = dol group = channel_type; 
	title "Dollar difference across retail formats for an average Domnicks consumer"; 
run;
proc gchart data = tmp4; 
	vbar year/discrete sumvar = share group = retailer_code; 
	title "Share difference within top grocery stores for an average Domnicks consumer"; 
run;
proc gchart data = tmp3; 
	vbar year/discrete sumvar = share group = channel_type; 
	title "Share difference across retail formats for an average Domnicks consumer"; 
run;


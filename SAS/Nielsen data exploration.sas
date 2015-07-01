libname mylib "E:\Users\ccv103\Documents\Research\Store switching\SAS_temp\Nielsen";
%let dirname= E:\Users\ccv103\Documents\Research\Nielsen\nielsen_extracts\HMS;

*************;
* Read data *;
*************;
* Read master retailer data;
%let fpath = &dirname\Master_Files\Latest\retailers.tsv;
%put &fpath;
PROC IMPORT OUT= WORK.retailers
            DATAFILE= "&fpath" 
            DBMS=dlm REPLACE;
			DELIMITER ='09'x;
     		DATAROW=2; 
			GUESSINGROWS=5000;
RUN;

%let fpath = &dirname\Master_Files\Latest\products.tsv;
%put &fpath;
PROC IMPORT OUT= WORK.products
            DATAFILE= "&fpath" 
            DBMS=dlm REPLACE;
			DELIMITER ='09'x;
     		DATAROW=2; 
			GUESSINGROWS=5000;
RUN;

* Read in selected markets;
PROC IMPORT OUT= my_market
            DATAFILE= "\\tsclient\Resear1\Store switching\processed data\selected market.csv" 
            DBMS=csv REPLACE;
     		GETNAMES=YES; 
RUN;

data my_channel;
	input channel_type $40.	channel_group $10.;
	datalines;
	Grocery								Grocery
	Drug Store							Grocery
	Dollar Store						Grocery
	Convenience Store					Grocery
	Health Food Store					Grocery
	Discount Store						Grocery
	Warehouse Club						Grocery
	Pizzeria							Restaurant
	Restaurant							Restaurant
	Quick Serve Restaurants				Restaurant
	Bakery								Restaurant
	;
run;	

proc sql noprint;
	create table my_retailers as
	select A.*, B. channel_group
	from retailers as A inner join my_channel as B
	on A.channel_type = B.channel_type;
quit;

* Read in 2004-2012 panelists data and trip data; 
%macro read_all_year;
	%do i=2004 %to 2012;
		%let fpath = &dirname\&i\Annual_Files\trips_&i..tsv;
		%put &fpath;
		PROC IMPORT OUT= WORK.trip&i
		            DATAFILE= "&fpath" 
		            DBMS=dlm REPLACE;
					DELIMITER ='09'x;
		     		DATAROW=2; 
		RUN;
		
		%let fpath = &dirname\&i\Annual_Files\panelists_&i..tsv;
		%put &fpath;
		PROC IMPORT OUT= WORK.panelists&i
		            DATAFILE= "&fpath" 
		            DBMS=dlm REPLACE;
					DELIMITER ='09'x;
		     		DATAROW=2; 
		RUN;	
		
		data panelists&i(drop = projection_factor_magnet_char);
			length scantrack_market_descr $50.;
			set panelists&i(rename = (projection_factor_magnet = projection_factor_magnet_char)); 
			projection_factor_magnet = input(projection_factor_magnet_char, best12.);
			drop wic_indicator_current wic_indicator_ever_notcurrent Member:;
		run;
		
		* Restirct to the focal markets;
		proc sql noprint;
			create table tmp as 
			select * from panelists&i 
			where scantrack_market_descr in (select scantrack_market_descr from my_market);
		quit;
		data panelists&i; set tmp; run;
		
		
		proc sql noprint; 
			* Restrict households to be in the focal market;		
			* Merge in panelists geographic data to trip data;
			create table tmp as 
			select A.*, year(A.purchase_date) as year, B.scantrack_market_descr
			from trips&i as A inner join panelists&i as B
			on A.household_code = B.household_code;	
					
			* Restrict the trip data to the selected retailers;
			create table trips&i as
			select A.*, B.channel_type
			from tmp as A inner join my_retailers as B
			on A.retailer_code = B.retailer_code;			
		quit;
	%end;
	
	* Combine panelists data together;
	data panelists;
		set panelists2004 - panelists2012;
	run;
	proc datasets noprint; delete panelists2004-panelists2012;
	
	* Combine trip data together;
	data trip;
		set trip2004 - trip2012;
	run;
	proc datasets noprint; delete trip2004-trip2012; run;
%mend read_all_year;

%read_all_year;

******************************************************************;
* Check the restaurant and gorcery spending from the full sample *;
******************************************************************;
data hh_rest_groc;
	input household_code panel_year channel_group $10. obs DOL;
run;

%macro rest_groc_fun;
	%do i=2004 %to 2012;
		%let fpath = &dirname\&i\Annual_Files\trips_&i..tsv;
		%put &fpath;
		PROC IMPORT OUT= WORK.trip&i
		            DATAFILE= "&fpath" 
		            DBMS=dlm REPLACE;
					DELIMITER ='09'x;
		     		DATAROW=2; 
		RUN;
		
		proc sql noprint; 
			create table tmp as 
			select A.*, B.channel_type, B.channel_group
			from trip&i as A inner join my_retailers as B
			on A.retailer_code = B.retailer_code;
			
			create table tmp1 as 
			select household_code, panel_year, channel_group, count(distinct purchase_date) as obs, sum(total_spent) as DOL
			from tmp
			group by household_code, panel_year, channel_group;
		quit;
		
		data hh_rest_groc;
			set hh_rest_groc tmp1;
		run;
		
		proc datasets noprint; delete trip&i tmp tmp1; run;
	%end;
%mend rest_groc_fun;

%rest_groc_fun;

* Transform the long data to wide form;
proc sort data = hh_rest_groc; by household_code panel_year channel_group; run;
proc transpose data = hh_rest_groc out = hh_wide;
	by household_code panel_year;
	id channel_group;
	var DOL;
run;

data hh_wide;
	set hh_wide;
	rat_rg = restaurant/Grocery;
run;

proc means data = hh_wide;
	class panel_year;
	var rat_rg;
run;

* Histogram of dollars spent in grocery and restaurants by year;
proc sort data = hh_wide; by panel_year; run;

goptions reset=all;
ods pdf file = "\\tsclient\Resear1\Store switching\processed data\graph_nielsen_restaurant_coverage.pdf";
proc sgpanel data=hh_wide;
  panelby panel_year / columns = 3;
  colaxis label = "Annual spending" MIN = 0 MAX = 10000;
  histogram Restaurant /transparency = .5 legendlabel = "Restaurant";
  histogram Grocery / transparency = .5 legendlabel = "Grocery";
run;

proc sgpanel data=hh_wide;
	panelby panel_year/ columns = 3;
	colaxis label = "Ratio of annaul restaurant spending to annual grocery expenditure" MIN = 0 MAX = 1;
	histogram rat_rg ;
	density rat_rg / type = Kernel;
run;
ods pdf close;

**************************;
* Magnet households only *;
**************************;
* The distribution of projection_factor_magnet;
proc univariate data=panelists;
	class panel_year;
	var projection_factor_magnet;
	output out=tmp NOBS=Nobs NMISS=Nmiss;
run;

data tmp;
	set tmp;
	magnet_prop = 1- Nmiss/Nobs;
proc print data=tmp; var panel_year magnet_prop; title 'The proportion of magnet households over year'; run;

* Distribution of projection_factor_magnet;
goptions reset=all;
proc sgpanel data=panelists; 
	panelby panel_year/ columns = 3;
	histogram projection_factor_magnet;
run;	

*-------------------------------------------------------------------;
* Aggregate the restaurant and gorcery spending by household by year; 
proc sql noprint; 
	* Restrict to the trips from households who are in magnet projection;
	create table tmp_trip as 
	select A.*
	from trip as A inner join (select household_code, panel_year from panelists where projection_factor_magnet^=.) as B
	on A.household_code = B.household_code and A.panel_year = B.panel_year;
	
	create table tmp as 
	select A.*, B.channel_type, B.channel_group
	from tmp_trip as A inner join my_retailers as B
	on A.retailer_code = B.retailer_code;
	
	create table hh_rest_groc as 
	select household_code, panel_year, channel_group, count(distinct purchase_date) as obs, sum(total_spent) as DOL
	from tmp
	group by household_code, panel_year, channel_group;
quit;

*-------------------------------------------------------------;
* The distribution of restaurant spending vs. grocery spending;
* Transform the long data to wide form;
proc sort data = hh_rest_groc; by household_code panel_year channel_group; run;
proc transpose data = hh_rest_groc out = hh_wide;
	by household_code panel_year;
	id channel_group;
	var DOL;
run;

data hh_wide;
	set hh_wide;
	rat_rg = restaurant/Grocery;
run;

proc means data = hh_wide;
	class panel_year;
	var rat_rg;
run;

* Histogram of dollars spent in grocery and restaurants by year;
proc sort data = hh_wide; by panel_year; run;

proc sgpanel data=hh_wide;
  panelby panel_year / columns = 3;
  colaxis label = "Annual spending" MIN = 0 MAX = 10000;
  histogram Restaurant /transparency = .5 legendlabel = "Restaurant";
  histogram Grocery / transparency = .5 legendlabel = "Grocery";
run;

proc sgpanel data=hh_wide;
	panelby panel_year/ columns = 3;
	colaxis label = "Ratio of annaul restaurant spending to annual grocery expenditure" MIN = 0 MAX = 1;
	histogram rat_rg ;
	density rat_rg / type = Kernel;
run;

*---------------------------------------------------------------;
* Restrict to households reporing sufficent restaurant incidence;
%let my_rest_trip_prop = 0.25;
proc transpose data = hh_rest_groc out = tmp1;
	by household_code panel_year;
	id channel_group;
	var obs;
run;

data tmp1;
	set tmp1;
	if restaurant=. then delete;
	if restaurant < &my_rest_trip_prop*(grocery + restaurant) then delete;
	magnet_flag = 1;
run;

proc sql noprint;
	create table hh_wide1 as
	select A.*
	from hh_wide as A inner join tmp1 as B
	on A.household_code = B.household_code and A.panel_year = B.panel_year
	order by panel_year;
	
	create table sub_panelists as
	select A.*, B.magnet_flag
	from panelists as A left join tmp1 as B
	on A.household_code = B.household_code and A.panel_year = B.panel_year
	order by household_code, panel_year;
quit;

* The sample proportion after slicing;
data sub_panelists;
	set sub_panelists;
	if magnet_flag =. then magnet_flag = 0;
run;
proc sort data=sub_panelists; by panel_year; run;
proc freq data=sub_panelists; 
	tables panel_year*magnet_flag;
run;

* The stay lenth of remaining households;
proc sql noprint;
	create table tmp as
	select household_code, min(panel_year) as start, max(panel_year) as endyear, (max(panel_year) - min(panel_year)) as stay_length
	from sub_panelists
	where magnet_flag = 1
	group by household_code
	order by household_code;
quit;

proc freq data=tmp; tables stay_length; run;

* Plot the histogram of restaurant vs. grocery spending;
proc sgpanel data=hh_wide1;
  panelby panel_year / columns = 3;
  colaxis label = "Annual spending" MIN = 0 MAX = 10000;
  histogram Restaurant /transparency = .5 legendlabel = "Restaurant";
  histogram Grocery / transparency = .5 legendlabel = "Grocery";
run;

proc sgpanel data=hh_wide1;
	panelby panel_year/ columns = 3;
	colaxis label = "Ratio of annaul restaurant spending to annual grocery expenditure" MIN = 0 MAX = 1;
	histogram rat_rg ;
	density rat_rg / type = Kernel;
run;

************************************************************************;
* Substituition between restaurant and grocery for the panelist subset *;
************************************************************************;
* The stay lenth of remaining households;
proc sql noprint;
	create table tmp as
	select household_code, min(panel_year) as start, max(panel_year) as endyear, (max(panel_year) - min(panel_year)) as stay_length
	from sub_panelists
	where magnet_flag = 1
	group by household_code
/*	having stay_length>1*/
	order by household_code;
quit;
proc freq data=tmp; tables stay_length; run;

* Select the detailed trip data for the qualified households;
proc sql noprint;
	create table tmp1 as 
	select A.*
	from trip as A inner join 
	(select household_code, panel_year from sub_panelists where magnet_flag=1 and household_code in (select household_code from tmp)) as B
	on A.household_code = B.household_code and A.panel_year = B.panel_year;
	
	create table tmp_trip as
	select A.*, B.channel_type, B.channel_group
	from tmp1 as A inner join my_retailers as B
	on A.retailer_code = B.retailer_code
	order by household_code, purchase_date;
quit;
proc contents data=tmp_trip; run;

* Generate date variables;
data tmp_trip;
	set tmp_trip;
	month	= month(purchase_date);
	week 	= week(purchase_date, 'u');
	year 	= year(purchase_date);
	ymonth	= put(purchase_date, monyy7.);
	yweek 	= put(purchase_date, weeku5.);
run;

* Aggregate the detailed trip data to weekly data;
proc sql noprint;
	create table tmp_rest_groc as
	select household_code, ymonth, yweek, channel_group, mean(year) as year, mean(month) as month, 
			sum(total_spent) as DOL
	from tmp_trip
	group by household_code, ymonth, yweek, channel_group
	order by household_code, ymonth, yweek;
quit;

* Transform the long data to wide data;
proc transpose data=tmp_rest_groc out=reg_data;
	by household_code ymonth yweek;
	id channel_group;
var DOL;

* Linear regression with household fixed effects;
data reg_data;
	set reg_data;
	if grocery=. then grocery= 0;
	if restaurant =. then restaurant = 0;
	month = substr(ymonth, 1, 3);
	year = substr(ymonth, 4, 4);
run;

proc glm data=reg_data;
	absorb household_code  ymonth;
	model restaurant = grocery/solution;
run;

* Go back and check demographic distribution of the qualified households;
ods graphics on;
proc freq data=sub_panelists;
	tables panel_year*magnet_flag*household_income /plots=freqplot; 
run;

*---------------------------------------------;
* Mosiac plot ;
filename mosaic  '\\tsclient\Resear1\Store switching\Exercise\SAS\mosaics';
* storage location of compiled macros;
libname  mosaic   '\\tsclient\Resear1\Store switching\Exercise\SAS\mosaics';

* Code to read in, compile and store the macros;
proc iml ;
   reset storage=mosaic.mosaic;
   %include mosaic(mosaics) ;
   store module=_all_;
   show storage;
quit;

* Read in the wrapper macro;
%include "\\tsclient\Resear1\Store switching\Exercise\SAS\mosaics\mosaic.sas";

* Prep: create the table, save the cell counts;
proc freq data = "c:\book\help.sas7bdat";
tables homeless * female / out=outhelp;
run;
* Make the plot;
%mosaic(data=outhelp,var = panel_year magnet_flag, 
        sort=panel_year, space = 1 1);


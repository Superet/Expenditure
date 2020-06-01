/*
This scripts selects households for analysis, and construct household-specific basket.
Output library: mylib, mysqllib
Output data: hh_bskt, keephh

Input data: mysqllib.trips, mysqllib.purchases, mysqllib.products. 

Household selection criterion: 
1. Stay in the panel for at least one year; 
2. The maximum interpurchase gap is less than 2 months/one month; 
3. Have at least 10 shopping trips. 

Household basket: 
The basket is the wallet share for all the categories, computed from the first 6 months from each household
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

*********************;
* Select households *; 
*********************;
proc sql noprint; 
	connect to oledb as mydb
	(OLEDB_SERVICES=NO Datasource="kdc01\kdwh02" PROVIDER=SQLOLEDB.1 
	Properties=('Initial Catalog'=USRDB_ccv103 'Integrated Security'=SSPI)
	SCHEMA=DBO);
	
	* Calculate the interpurchase days; 
	create table tmp as 
	select * from connection to mydb (
		select household_code, year, purchase_date, 
		DATEDIFF(day, LAG(purchase_date, 1) OVER (PARTITION BY household_code ORDER BY purchase_date), purchase_date) as inter_date
		from dbo.trips;
	);
	disconnect from mydb;

	* For each household, summarize stay lengths, number of shopping trips, and max interpurchase days; 
	create table tmp1 as 
	select household_code, intck('day', DATEPART(min(purchase_date)), DATEPART(max(purchase_date))) as days, count(unique(purchase_date)) as ntrips, max(inter_date) as gap
	from tmp
	group by household_code; 
quit; 

* Calcualte the second largest interpurchase days; 
proc sort data = tmp out = tmp2; by household_code descending inter_date; run;
data tmp3; 
	set tmp2; 
	by household_code descending inter_date; 
	if first.household_code then n = 0; 
		n+1;
	if n=2 then output;
run;

data tmp1; 
	merge tmp1 tmp3(keep = household_code inter_date); 	
	by household_code; 
run;
proc datasets noprint; delete tmp2 tmp3; run;

* Check how many households need to be dropped; 
* NOTE: some households leaf the panel for a year or two and then came back, ; 
* I still keep them if the second largest interpurchase date is small then the cutoff; 
%let min_stay	= 365; 
%let min_trip 	= 10; 
%let max_gap 	= 30; 

proc sql; 
	select sum(days<&min_stay)/count(*) as drop_length, sum(ntrips<&min_trip)/count(*) as drop_trips, 
		sum(gap>&max_gap)/count(*) as drop_gap, 
		sum(gap>&max_gap and inter_date >&max_gap)/count(*) as drop_2gap, 
		sum(days<&min_stay or ntrips<&min_trip or gap>&max_gap)/count(*) as drop_overall, 
		sum(days<&min_stay or ntrips<&min_trip or (gap>&max_gap and inter_date>&max_gap))/count(*) as drop_overall2
	from tmp1;
quit; 

* Filter out households who do not satisify the selection criterion; 
proc sql noprint; 
	create table keephh(rename = (inter_date = gap2)) as 
	select *
	from tmp1 
	where days >= &min_stay and ntrips >= &min_trip and ^(gap>&max_gap and inter_date>&max_gap)
	order by household_code; 
quit; 

data mylib.keephh; set keephh; run; 
data mysqllib.keephh; set keephh; run;

************************************************************; 
* Construct household basket using the first 6-month trips *; 
************************************************************; 
proc sql noprint; 
	* Subset trips within the frist 6 months; 
	create table tmp as 
	select household_code, trip_code_uc
	from (select * from mainex.trips where household_code in (select household_code from keephh)) 
	group by household_code
	having purchase_date <= min(purchase_date) + 180;
	
	create table mysqllib.tmp as select * from tmp; 
	
	connect to oledb as mydb
	(OLEDB_SERVICES=NO Datasource="kdc01\kdwh02" PROVIDER=SQLOLEDB.1 
	Properties=('Initial Catalog'=USRDB_ccv103 'Integrated Security'=SSPI)
	SCHEMA=DBO);	

	EXECUTE(
		/*Subset the purchases from the subset of trips; */
		select A.*, B.product_module_descr, B.department_descr
		into dbo.tmp1
		from (select C.*, D.household_code from dbo.purchases as C inner join dbo.tmp as D on C.trip_code_uc = D.trip_code_uc) as A
		left join dbo.products as B
		on A.upcv = B.upcv; 
	
		/*Sum up the expenditure for each department; */
		select household_code, department_descr, product_module_descr, sum(dol_paid) as dol
		into dbo.hh_dpt
		from dbo.tmp1
		group by household_code, department_descr, product_module_descr; 
	
		drop table dbo.tmp, dbo.tmp1;  
	
	) by mydb; 
	disconnect from mydb;	
quit; 

* Calculate wallet share by each household; 
proc sql noprint; 
	create table hh_basket as
	select household_code, product_module_descr, dol/sum(dol) as wallet_share
	from (select * from mysqllib.hh_dpt where department_descr^="")
	group by household_code
	order by household_code; 
quit; 

data mylib.hh_basket; set hh_basket; run;
proc datasets library = mysqllib noprint; delete hh_dpt; run;


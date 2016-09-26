/*
This script generates the monthly expenditure spent on 6 big departments
Output library: mylib
Output data: hh_month_dpt

Original department: ALCOHOLIC BEVERAGES, DAIRY, DELI, DRY GROCERY, FRESH PRODUCE, FROZEN FOODS, GENERAL MERCHANDISE, 
					HEALTH & BEAUTY CARE, MAGNET DATA, NON-FOOD GROCERY, PACKAGED MEAT. 
Re-classified department: 	REFRIGERATED & FROZEN (DAIRY, DELI, FROZEN FOODS, PACKAGED MEAT), 
							DRY GROCERY (ALCOHOLIC BEVERAGES), 
							NON-FOOD GROCERY, 
							GENERAL MERCHANDISE, 
							HEALTH & BEAUTY CARE, 
							PRODUCE & OTHER (FRESH PRODUCE, MAGNET DATA)

NOTE: 
1. Expenditure should delete UPCs that are not commonly carried by supermarkets, mysqllib.products has delete them.  							
2. Nielsen adds new products on product & other since 2007, so tha expenditure on produce and other increased after 2007. 
	I use total expenditure reported from trip minus the undesired categories as product & other. 
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

proc contents data = mainex.trips; run;
proc contents data = mainex.purchases; run;

* ---------------------------------------------------------------------------------*; 
* Aggregate purchase data to household monthly purchases: by department by channel *; 
proc sql noprint; 
	connect to oledb as mydb
	(OLEDB_SERVICES=NO Datasource="kdc01\kdwh02" PROVIDER=SQLOLEDB.1 
	Properties=('Initial Catalog'=USRDB_ccv103 'Integrated Security'=SSPI)
	SCHEMA=DBO);

	EXECUTE(		
		/*Append household and purchase date to the data and subset households */
		select household_code, purchase_date, DATEPART(yyyy,purchase_date) as year,DATEPART(month,purchase_date) as month, 
			   B.trip_code_uc, B.channel_type, upcv, A.dol_paid
		into dbo.tmp
		from dbo.purchases as A 
		right join 
		(select * from dbo.trips where household_code in (select household_code from dbo.keephh) ) as B
		on A.trip_code_uc = B.trip_code_uc;
	
		/*Combine purchase data and product category data by UPC, note right join products trop drop undesireable products */
		/*Aggregate across household month; */
		select household_code, year, month, department, channel_type, sum(dol_paid) as dol_paid
		into dbo.tmp1
		from (	select A.*, 
				CASE WHEN department_descr in ('DAIRY', 'DELI', 'FROZEN FOODS', 'PACKAGED MEAT') THEN 'RF'
					 WHEN department_descr in ('ALCOHOLIC BEVERAGES', 'DRY GROCERY') THEN 'DG'
					 WHEN department_descr in ('FRESH PRODUCE', 'MAGNET DATA') THEN 'PRODUCE OTHER'
					 WHEN department_descr = 'HEALTH & BEAUTY CARE' THEN 'HBC'
					 WHEN department_descr = 'NON-FOOD GROCERY' THEN 'NFG'
					 WHEN department_descr = 'GENERAL MERCHANDISE' THEN 'GM'
					 WHEN department_descr = '' THEN 'UNCLASSIFIED'
				ELSE department_descr
				END as department
			  	from dbo.tmp as A right join dbo.products as B on A.upcv = B.upcv) as C
		group by household_code, year, month, department, channel_type
		order by household_code, year, month, department, channel_type; 
		
		/*Calculate the total expenditure and expenditure on focal categories at each channe; */
		select household_code, year, month, channel_type, sum(total_spent) as total_spent, 
			sum(net_spent) as focal_spent
		into dbo.tmp2
		from (	select household_code, DATEPART(yyyy,purchase_date) as year,DATEPART(month,purchase_date) as month, 
				channel_type, total_spent, net_spent 
				from dbo.trips where household_code in (select household_code from dbo.keephh) ) as C
		group by household_code, year, month, channel_type; 
		
	) by mydb; 
	%put &SQLXRC. &SQLXMSG. ; 
	disconnect from mydb;	
quit;		
proc freq data = mysqllib.tmp1; table department; run;

* ------------------------------------------- *;
* Check the expenditure on 'other' department *; 
proc sql noprint;
	create table tmp as 
	select A.*, ifn(dol_produce=., 0, dol_produce) as dol_produce, ifn(dol_dpt=., focal_spent, focal_spent - dol_dpt) as dol_other
	from mysqllib.tmp2 as A left join 
	(select household_code, year, month, channel_type, 
		sum(ifn(department='PRODUCE OTHER', 1, 0)*dol_paid) as dol_produce, 
		sum(ifn(department='PRODUCE OTHER', 0, 1)*dol_paid) as dol_dpt
	from mysqllib.tmp1 group by household_code, year, month, channel_type) as B
	on A.household_code=B.household_code and A.year = B.year and A.month = B.month and A.channel_type = B.channel_type; 
	
	create table tmp1 as 
	select year, month, mdy(month, 1, year) format = date9. as t, mean(dol_other) as dol_other, mean(dol_produce) as dol_produce
	from tmp
	group by year, month;
quit; 

proc univariate data = tmp; var dol_other dol_produce; histogram; run;
proc sgplot data = tmp1; 
	series x = t y = dol_produce;
	series x = t y = dol_other; 
	xaxis values=('01JAN2004'd to '01DEC2012'd by year);
run;

* Adjust other expenditure; 
data tmp; 
	set tmp; 
	if dol_other <0 then dol_other=0;
run;

* -------------------------------------- *; 
* Transpose the long data to wide format *; 
proc transpose data=mysqllib.tmp1(where = (department^='PRODUCE OTHER')) out=tmp1(drop = _NAME_ _LABEL_) prefix = EXP_; 
	by household_code year month; 
	id department channel_type; 
	var dol_paid; 
run; 

proc transpose data = tmp out=tmp2(drop = _NAME_ _LABEL_) prefix = REPORT_;
	by household_code year month; 
	id channel_type; 
	var focal_spent; 
run;

proc transpose data = tmp out = tmp3(drop = _NAME_ _LABEL_) prefix = EXP_OTHER; 
	by household_code year month; 
	id channel_type; 
	var dol_other; 
run;

proc sql noprint;
	create table tmp(drop = old_hh old_year old_month) as 
	select *
	from (	select A.*, B.EXP_OTHERConvenience_Store, B.EXP_OTHERDiscount_Store, B.EXP_OTHERDollar_Store, 
			B.EXP_OTHERDrug_Store, B.EXP_OTHERGrocery, B.EXP_OTHERWarehouse_Club
			from tmp2 as A left join tmp3 as B
		 	on A.household_code=B.household_code and A.year = B.year and A.month = B.month ) as C
	left join tmp1(rename = (household_code = old_hh year=old_year month=old_month )) as D
	on C.household_code=D.old_hh and C.year = D.old_year and C.month = D.old_month 
	order by household_code, year, month;
quit; 
* Check if data has the full household-month: if the final observation is equivalent to the number of observation of tmp1; 

* Set missing values to zero; 
data hh_month_dpt; 
	set tmp; 
	array dol{*} REPORT: EXP_:;
	do i = 1 to dim(dol); 
		if dol{i} = . then dol{i} = 0; 
	end; 
	drop i; 
run; 

proc contents data = hh_month_dpt; run;
proc datasets library=mysqllib noprint; delete tmp tmp1 tmp2; run;

data mylib.hh_month_dpt; 
	set hh_month_dpt; 
	where household_code ^= . and year > 2003; 
run;


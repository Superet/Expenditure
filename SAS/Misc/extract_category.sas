%let dirname= E:\Users\ccv103\Documents\Research\Nielsen\nielsen_extracts\HMS;

* Read master product data;
%let fpath = &dirname\Master_Files\Latest\products.tsv;
%put &fpath;
PROC IMPORT OUT= WORK.orig_products
            DATAFILE= "&fpath" 
            DBMS=dlm REPLACE;
			DELIMITER ='09'x;
     		DATAROW=2; 
			GUESSINGROWS=5000;
RUN;

data myprod; 
	set orig_products;
	where product_module_descr in ("YOGURT-REFRIGERATED");
	upcv = cats(upc, "-", upc_ver_uc);
run;

* Read UPC-level purchase data;
%macro read_groc_fun;
	%do i=2004 %to 2012;
		* Import trip data;
		%let fpath = &dirname\&i\Annual_Files\trips_&i..tsv;
		%put &fpath;
		PROC IMPORT OUT= WORK.trips&i
		            DATAFILE= "&fpath" 
		            DBMS=dlm REPLACE;
					DELIMITER ='09'x;
		     		DATAROW=2; 
		RUN;
	
		* Import purchase data;
		%let fpath = &dirname\&i\Annual_Files\purchases_&i..tsv;
		%put &fpath;
		PROC IMPORT OUT= WORK.purchases&i
		            DATAFILE= "&fpath" 
		            DBMS=dlm REPLACE;
					DELIMITER ='09'x;
		     		DATAROW=2; 
		RUN;
		
		proc sql noprint;
			create table tmp as 
			select A.*, B.upc_descr, B.brand_descr, B.multi, B.size1_amount, B.size1_units
			from (select *, cats(upc, "-", upc_ver_uc) as upcv from purchases&i) as A inner join myprod as B
			on A.upcv = B.upcv;
			
			create table purchases&i as 
			select B.household_code, B.purchase_date, A.* 
			from tmp as A left join trips&i as B
			on A.trip_code_uc = B.trip_code_uc;
		quit;
	%end;
	data purchases;
		set purchases2004 - purchases2012;
	run;
	proc sort data=purchases; by hosuehold_code purchase_date; run;
%mend read_groc_fun;

%read_groc_fun;

proc contents data = purchases; run;
proc datasets noprint; delete purchases2004 - purchases2012 trips2004-trips2012; run;

* Export the data;
PROC EXPORT DATA = purchases
			OUTFILE= "E:\Users\ccv103\Desktop\yogurt_purchases.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;



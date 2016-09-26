

libname mainex "U:\Users\ccv103\Documents\Research\Store switching\SAS_temp\Nielsen\Work\main";
libname mylib "U:\Users\ccv103\Documents\Research\Store switching\SAS_temp\Nielsen\Work\mylib";
%let dirname= U:\Users\ccv103\Documents\Research\Store switching\store_data;
options compress=yes reuse=yes;

proc contents data = mylib.store_location; run;

data my_channel;
	input retailer_name $30.	channel_type $30.;
	datalines;
	7_eleven						Convenience Store
	walmart							Discount Store
	target							Discount Store
	kmart							Discount Store
	bj_s							Warehouse Club
	costco							Warehouse Club
	sam_s_club						Warehouse Club
	family_dollar_stores			Dollar Store
	dollar_tree						Dollar Store
	dollar_general					Dollar Store
	cvs								Drug Store	
	rite_aid						Drug Store
	walgreens						Drug Store					
	;
run;

* A function of reading data;
%macro readfile(dirlist,outname);
	proc sql noprint;
		select count(fname) into :num_files
		from &dirlist;
	quit;

	%do j=1 %to &num_files;
		proc sql noprint;
		select fname into :fname
		from &dirlist
		where n=&j;
		quit;
	
		%let fpath=&dirname\&fname;
		%put &fpath;
		PROC IMPORT OUT= tmp&j
		            DATAFILE= "&fpath" 
		            DBMS=csv REPLACE;
		     		GETNAMES=YES; 
					GUESSINGROWS=5000;
		RUN;
		
		%if &j=&num_files %then %do; 
			data tmp&j; 
				format store_name $60. zipcode 5. opening_date MMDDYY10. ; 
				set tmp&j(keep = zip_code latitude longitude opening_date); 
				store_name = "walmart"; 
				store_number = _n_; 
				zipcode = zip_code; 
				drop zip_code; 
			run;
		%end; 
		%else %do; 
			data tmp&j; 
				format store_name $60. zipcode 5. opening_date MMDDYY10. ; 				
				set tmp&j(keep = zip_code latitude longitude); 
				opening_date = .;
				store_name = tranwrd("&fname",".csv",""); 
				store_number = _n_; 
				if vtype(zip_code)='C' then do;
					zipcode = input(scan(zip_code, 1, '-'), 5.); 
				end;
				else zipcode = zip_code;
				drop zip_code;   
			run;
		%end; 
	%end;

	%let tmpend=%sysfunc(catx(%str(),tmp,&num_files));
	%put &tmpend;
	data &outname;
		set tmp1-&tmpend;
	run;
	
	* Delete separate files;
	proc datasets noprint; delete tmp1-&tmpend;run;
%mend readfile;

filename dirlist pipe 'dir /b "U:\Users\ccv103\Documents\Research\Store switching\store_data\*.csv" ';
data dirlist;
	length fname $30;
	infile dirlist length=reclen;
	input fname $varying30. reclen;
	n=_n_;
	*call symput ('num_files',_n_);
run;
proc print data=dirlist;run;

%readfile(dirlist, orig_store);
proc contents data = orig_store; run;
proc sql; select nmiss(zipcode) from orig_store; quit;

* Extract only location and opening dates from data; 
* Expand the data to panel data from 2004-2012; 
data tmp; 
	set orig_store(where = (store_name='walmart')); 
	year = 2003; 
	do while (year <= 2012); 
		year + 1; output; 
	end; 
run;

* Drop the stores that were not open in the concurrent year; 
data tmp; 
	set tmp;
	oyear = year(opening_date); 
	if oyear > year then delete; 
	drop oyear; 
run;

data store; 
	merge orig_store(where = (store_name^='walmart')) tmp; 
	by store_name store_number; 
	store_id = cats(store_name, '-', store_number); 
	store_zip3 = input(substr(put(zipcode, z5.),1, 3), 3.);
run;

* Merge with the store location data imputed from data; 
proc sql noprint;
	create table tmp1(drop=opening_date store_number) as 
	select A.*, B.channel_type
	from store as A left join my_channel as B
	on A.store_name = B.retailer_name; 
	
/*	create table tmp2 as 
	select put(retailer_code, 10.) as store_name, zipcode, latitude, longitude, year, 
		put(store_code_uc, 10.) as store_id, store_zip3, channel_type
	from mylib.store_location; */
	
	* Exclude WalMart Target and KMart stores from the store_location data imputed from the data; 
	* The top 3 retailers in revenue: 6920, 6901, 6905; 
	create table tmp2 as 
	select put(retailer_code, 10.) as store_name, zipcode, latitude, longitude, year, 
		put(store_code_uc, 10.) as store_id, store_zip3, channel_type
	from mylib.store_location
/*	where retailer_code not in (6920, 6901, 6905);*/
	where channel_type ^= 'Discount Store'; 
quit; 

data store_location_app; set tmp1 tmp2; run;
proc contents data = store_location_app; run;

data mylib.store_location_app; set store_location_app; run;

* If exluding discount store store_location data, save a new file; 
data mylib.store_location_app_excl; set store_location_app; run;



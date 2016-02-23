* This file is to generate group_year data for the product groups which do not have sales in 2006; 
* NOTE: for product group 5520, I change the name from 5520_2007.csv to 5520_2006.csv.; 

libname mymaster "/sscc/home/c/ccv103/rms_attributes/master";
libname	grpyear "/sscc/home/c/ccv103/rms_attributes/group_year";
libname	newprod "/sscc/home/c/ccv103/rms_attributes/newprod";
%let dirname= /sscc/home/c/ccv103/rms_attributes;
options compress=yes reuse=yes; 

%macro assign_value(single); 
	%global sel_group sel_year; 
	%if (&single = 1) %then %do; 
		%let sel_group 	= %sysget(sel_group); 
		%let sel_year	= %sysget(sel_year); 
	%end; 
	%else %do; 
		%let grp_idx =  %sysget(PBS_ARRAY_INDEX); 
		%let sel_year	= %sysget(sel_year); 
		proc sql noprint;
			select put(product_group_code, z4.) into: sel_group from mymaster.prod_group_index
			where n = &grp_idx; 
		quit; 				
	%end; 
%mend assing_value; 

%assign_value(%sysget(single));
%put &sel_group; 
%put &sel_year;

%let jmdata=joint_move&sel_group;
%let upcwt	= upcwt&sel_group; 
%put &jmdata; 
%put &upcwt;

* Read in Movement data for a single product_group and a single year *; 
%let subdir	= &dirname/nielsen_extracts/RMS/&sel_year/Movement_Files/&sel_group._&sel_year;
%put &subdir; 
filename dirlist pipe "ls &subdir";
data dirlist;
	length fname $30 filep $300.;
	infile dirlist length=reclen;
	input fname $varying30. reclen;
	filep = cats("&subdir", "/",fname); 
run;	

data movement;
	set dirlist(keep = filep);
	infile dummy filevar = filep delimiter='09'x MISSOVER DSD LRECL=32767 FIRSTOBS=2 end = done;
	do while(not done);
		informat week_end YYMMDD10.;
		input store_code_uc upc week_end units prmult price feature display; 
		output; 
	end;
	drop filep;
	format week_end YYMMDD10.;
run;

%let fpath = %sysfunc(catx(%str(), &dirname, /nielsen_extracts/RMS/&sel_year/Annual_Files/rms_versions_, &sel_year, .tsv));
%put &fpath;
PROC IMPORT OUT= WORK.upc_ver
            DATAFILE= "&fpath" 
            DBMS=dlm REPLACE;
			DELIMITER ='09'x;
     		DATAROW=2; 
			GUESSINGROWS=5000;
RUN;

* Merge in UPC version; 
proc sql noprint;
	create table &jmdata as 
	select C.*, D.product_module_code, D.brand_descr, D.size, D.size_index, D.prvt_flag
	from (select A.*, B.upc_ver_uc
		from movement as A left join upc_ver as B on A.upc = B.upc) as C
	left join mymaster.products as D
	on C.upc = D.upc and C.upc_ver_uc = D.upc_ver_uc; 
quit; 
proc datasets noprint; delete dirlist movement stores upc_ver; run;

*----------------------*;
* Generate UPC weights ; 
proc sql noprint;
	create table &upcwt as 
	select store_code_uc, upc, upc_ver_uc, product_module_code, sum(size*units) as vol_weight
	from &jmdata 
	group by store_code_uc, product_module_code, upc, upc_ver_uc
	order by store_code_uc, product_module_code, upc, upc_ver_uc; 
quit; 

* Export data; 
%let outfile = %sysfunc(catx(%str(), &dirname,/master/upc_weights/,&sel_group,_, &sel_year, .csv));
%put &outfile;
PROC EXPORT DATA= &upcwt
            OUTFILE= "&outfile"
            DBMS=csv REPLACE;
     PUTNAMES=YES;
RUN;

*---------------------------------------------------------------*;
* Generate retail attributes for a product module within a stre *;
%let outf = grpyear.swm_&sel_group._&sel_year; 
%put &outf;
proc sql noprint;
* Merge the scan data with the fixed UPC weights; 
	create table new&jmdata as 
	select A.*, B.vol_weight 
	from &jmdata as A left join &upcwt as B
	on A.store_code_uc = B.store_code_uc and A.upc = B.upc and A.upc_ver_uc = B.upc_ver_uc; 
quit; 

proc datasets noprint; 
	delete &jmdata; 
	change new&jmdata = &jmdata; 
run;

proc sql noprint;
* Compute store-week-module attributes; 
	create table &outf as 
	select store_code_uc, week_end, product_module_code, 
			sum(price/prmult*units) as revenue, sum(units) as units, sum(units*size) as quantity,
			sum(price/prmult/size*vol_weight)/sum(vol_weight) as unit_price_fixwt, 
			sum(size_index*units*size)/sum(units*size) as size_index,
			count(unique(brand_descr)) as num_brand, 
			count(unique(upc)) as num_upc, 
			sum(prvt_flag) as num_prvtlb
	from &jmdata
	group by store_code_uc, week_end, product_module_code
	order by store_code_uc, week_end, product_module_code;
quit; 

*-----------------------------*;
* Generate UPC national sales *; 
* Aggregate UPC sales over stores for a given week;
proc sql noprint;
	create table upc_sale as 
	select *, revenue/sum(revenue) as mkt_share
	from (select upc, upc_ver_uc, week_end, product_module_code, 
			sum(price/prmult*units) as revenue, sum(units) as units,
			sum(price/prmult*units)/sum(units) as price, 
			count(unique(store_code_uc)) as num_store
			from &jmdata
			group by upc, upc_ver_uc, week_end, product_module_code) 
	group by product_module_code
	order by upc, upc_ver_uc, week_end;
quit;

* Export data; 
%let outfile = %sysfunc(catx(%str(), &dirname,/newprod/upc_sales/&sel_year/,&sel_group,_, &sel_year, .csv));
%put &outfile;
PROC EXPORT DATA= upc_sale
            OUTFILE= "&outfile"
            DBMS=csv REPLACE;
     PUTNAMES=YES;
RUN;

proc datasets lib = work kill memtype=data noprint; run; 

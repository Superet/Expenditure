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
%let outf = grpyear.swm_&sel_group._&sel_year; 
%put &jmdata; 
%put &upcwt;
%put &outf;

* Read in UPC weight data; 
%let fpath = %sysfunc(catx(%str(), &dirname,/master/upc_weights/,&sel_group,_2006.csv));
%put &fpath;
PROC IMPORT OUT= &upcwt
            DATAFILE= "&fpath" 
            DBMS=csv REPLACE;
     		GETNAMES=YES;
			GUESSINGROWS=5000; 
RUN;
proc sort data = &upcwt; by product_module_code store_code_uc upc upc_ver_uc; run;
proc datasets library = work noprint;
	modify &upcwt;
	index create product_module_code; 
run;

* Read in UPC version; 
%let fpath = %sysfunc(catx(%str(), &dirname, /nielsen_extracts/RMS/&sel_year/Annual_Files/rms_versions_, &sel_year, .tsv));
%put &fpath;
PROC IMPORT OUT= WORK.upc_ver
            DATAFILE= "&fpath" 
            DBMS=dlm REPLACE;
			DELIMITER ='09'x;
     		DATAROW=2; 
			GUESSINGROWS=5000;
RUN;
proc datasets library = work noprint;
	modify upc_ver;
	index create upc;
run;

* Get the file names under the folder of group year; 
%let subdir	= &dirname/nielsen_extracts/RMS/&sel_year/Movement_Files/&sel_group._&sel_year;
%put &subdir; 
filename dirlist pipe "ls &subdir";
data dirlist;
	length fname $30 filep $300.;
	infile dirlist length=reclen;
	input fname $varying30. reclen;
	filep = cats("&subdir", "/",fname); 
	product_module_code = substr(fname, 1, 4); 
	n = _N_;
run;

* Get the loop number; 
proc sql noprint; 
	select max(n) into: num_mod from dirlist; 
quit; 
%put &num_mod; 

* Loop over product modules; 
%macro readmerg;
	%do i=1 %to &num_mod; 
		* Read movement file for one module; 
		data movement(drop=feature display); 
			set dirlist(where = (n = &i) keep = filep n);
			infile dummy filevar = filep delimiter='09'x MISSOVER DSD LRECL=32767 FIRSTOBS=2 end = done;
			do while(not done);
				informat week_end YYMMDD10.;
				input store_code_uc upc week_end units prmult price feature display; 
				output; 
			end;
			drop filep n;
			format week_end YYMMDD10.;
		run;

		proc sql noprint;
			* Select the focal product_module_code; 
			select compress(product_module_code) into: sel_mod from dirlist where n = &i; 
			%put &sel_mod;

			* Create subset of products for the focal product module; 
			create table subprod as 
			select *
			from mymaster.products where product_module_code = &sel_mod
			order by upc, upc_ver_uc; 

			* Merge movement data and product characteristics; 
			create table tmp as 
			select C.*, D.product_module_code, D.brand_descr, D.size, D.size_index, D.prvt_flag
			from (select A.*, B.upc_ver_uc
				from movement as A left join (select * from upc_ver where upc in (select upc from subprod)) as B on A.upc = B.upc) as C
			left join subprod as D
			on C.upc = D.upc and C.upc_ver_uc = D.upc_ver_uc;

			drop table movement; 

			* Merge the joint data with baseline upc weight in each store; 
			create table &jmdata as 
			select A.*, B.vol_weight 
			from tmp as A left join &upcwt(where=(product_module_code = &sel_mod)) as B
			on A.store_code_uc = B.store_code_uc and A.upc = B.upc and A.upc_ver_uc = B.upc_ver_uc;

			drop table tmp; 

			* Compute store-week-module attributes; 
			create table mod&i as 
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

			* Aggregate UPC national sales for a given week; 
			create table upcsale&i as 
			select *, revenue/sum(revenue) as mkt_share
			from (select upc, upc_ver_uc, week_end, product_module_code, 
					sum(price/prmult*units) as revenue, sum(units) as units,
					sum(price/prmult*units)/sum(units) as price, 
					count(unique(store_code_uc)) as num_store
					from &jmdata
					group by upc, upc_ver_uc, week_end, product_module_code) 
			group by product_module_code
			order by upc, upc_ver_uc, week_end;

			drop table &jmdata, subprod; 
		quit; 

		* Stack the processed data; 
		%if &i=1 %then %do; 
			data store_upc; set mod1; run;
			data upc_sale; set upcsale1; run;
			proc datasets noprint; delete mod1 upcsale1; run;
		%end; 
		%else %do; 
			data store_upc; set store_upc mod&i; run;
			data upc_sale; set upc_sale upcsale&i; run;
			proc datasets noprint; delete mod&i upcsale&i; run;
		%end;
	%end; 
%mend readmerg; 

%readmerg; 

* Export data; 
data &outf; set store_upc; run;
%let outfile = %sysfunc(catx(%str(), &dirname,/newprod/upc_sales/&sel_year/,&sel_group,_, &sel_year, .csv));
%put &outfile;
PROC EXPORT DATA= upc_sale
            OUTFILE= "&outfile"
            DBMS=csv REPLACE;
     PUTNAMES=YES;
RUN;

proc datasets lib = work kill memtype=data noprint; run;

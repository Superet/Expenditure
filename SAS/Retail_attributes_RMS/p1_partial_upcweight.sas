libname mymaster "/sscc/home/c/ccv103/rms_attributes/master";
libname	grpyear "/sscc/home/c/ccv103/rms_attributes/group_year";
libname	newprod "/sscc/home/c/ccv103/rms_attributes/newprod";
%let dirname= /sscc/home/c/ccv103/rms_attributes;

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



*******************;
* Macro Functions *; 
*******************;
* Read in Movement data for a single product_group and a single year *; 
%macro read_merge; 
/*	%let subdir	= %sysfunc(catx(%str(), &dirname,/nielsen_extracts/RMS/, &sel_year,/Movement_Files/, &sel_group,_, &sel_year));*/
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
	
	%let fpath = %sysfunc(catx(%str(), &dirname, /nielsen_extracts/RMS/&sel_year/Annual_Files/stores_, &sel_year, .tsv));
	%put &fpath;
	PROC IMPORT OUT= WORK.stores
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
%mend read_merge; 

*-----------------------------------------------------------*
* Generate UPC volume weight within a module within a store *; 
%macro gen_upcwt;  
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
%mend gen_upcwt;


******************************************************************;
* Run the functions for a single product group and a single year *;
******************************************************************;
%read_merge;
%gen_upcwt; 

proc datasets lib = work kill memtype=data noprint; run; 

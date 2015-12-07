libname mymaster "E:\Users\ccv103\Documents\Research\Nielsen\RMS_attributes\Master";
libname	grpyear "E:\Users\ccv103\Documents\Research\Nielsen\RMS_attributes\Group_year";
libname	newprod "E:\Users\ccv103\Documents\Research\Nielsen\RMS_attributes\Newprod";
%let dirname= E:\Users\ccv103\Documents\Research\Nielsen\RMS_attributes;

%let sel_group 	= 4001; 
%let sel_year	= 2006; 

*******************;
* Macro Functions *; 
*******************;
* Read in Movement data for a single product_group and a single year *; 
%macro read_merge(add_fixwt = 1); 
	/*%let subdir	= %sysfunc(catx(%str(), &dirname,\nielsen_extracts\RMS\, &sel_year,\Movement_Files\, &sel_group,_, &sel_year));*/
	%let subdir	= &dirname/nielsen_extracts/RMS/&sel_year/Movement_Files/&sel_group._&sel_year;
	%put &subdir; 
	filename dirlist pipe "dir /B &subdir\*.tsv";
	data dirlist;
		length fname $30 filep $300.;
		infile dirlist length=reclen;
		input fname $varying30. reclen;
		filep = cats("&subdir", "\",fname); 
	run;
	proc print data = dirlist; run;

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
	
	%let fpath = %sysfunc(catx(%str(), &dirname, \nielsen_extracts\RMS\&sel_year\Annual_Files\rms_versions_, &sel_year, .tsv));
	%put &fpath;
	PROC IMPORT OUT= WORK.upc_ver
	            DATAFILE= "&fpath" 
	            DBMS=dlm REPLACE;
				DELIMITER ='09'x;
	     		DATAROW=2; 
				GUESSINGROWS=5000;
	RUN;
	
	%let fpath = %sysfunc(catx(%str(), &dirname, \nielsen_extracts\RMS\&sel_year\Annual_Files\stores_, &sel_year, .tsv));
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
		create table joint_move as 
		select C.*, D.product_module_code, D.brand_descr, D.size, D.size_index, D.prvt_flag
		from (select A.*, B.upc_ver_uc
			from movement as A left join upc_ver as B on A.upc = B.upc) as C
		left join mymaster.products as D
		on C.upc = D.upc and C.upc_ver_uc = D.upc_ver_uc; 
	quit; 
	
	%if(add_fixwt = 1) %then %do; 
		proc sql noprint;
			create table joint_move1 as 
			select A.*, B.vol_weight 
			from joint_move as A left join mymaster.upc_weight as B
			on A.store_code_uc = B.store_code_uc and A.upc = B.upc and A.upc_ver_uc = B.upc_ver_uc; 
		quit; 
		
		proc datasets noprint; 
			delete joint_move; 
			change joint_move1 = joint_move; 
		run; 
	%end; 
	
	proc datasets noprint; delete dirlist movement; run;
%mend read_merge; 

*-----------------------------------------------------------*
* Generate UPC volume weight within a module within a store *; 
%macro gen_upcwt(append_data = 1); 
	%read_merge(add_fixwt = 0); 
	proc sql noprint;
		create table upcwt as 
		select store_code_uc, upc, upc_ver_uc, product_module_code, sum(size*units) as vol_weight
		from joint_move 
		group by store_code_uc, product_module_code, upc, upc_ver_uc
		order by store_code_uc, product_module_code, upc, upc_ver_uc; 
	quit; 
	
	%if (append_data = 1) %then %do; 
		data mymaster.upc_weight; 
			set mymaster.upc_weight upcwt;  
		run;
	%end;
	%else %do; 
		data mymaster.upc_weight; 
			set upcwt; 
		run;
	%end; 
%mend gen_upcwt; 

*---------------------------------------------------------------*;
* Generate retail attributes for a product module within a stre *;
%macro gen_store_module_week(joint_move); 
	%let outf = %sysfunc(catx(%str(), grpyear., swm_, &sel_group, _, &sel_year)); 
	%put &outf;
	proc sql noprint;
	* Merge the scan data with the fixed UPC weights; 
		create table joint_move1 as 
		select A.*, B.vol_weight 
		from &joint_move as A left join mymaster.upc_weight as B
		on A.store_code_uc = B.store_code_uc and A.upc = B.upc and A.upc_ver_uc = B.upc_ver_uc; 
	
	* Compute store-week-module attributes; 
		create table &outf as 
		select store_code_uc, week_end, product_module_code, 
				sum(price/prmult*units) as revenue, sum(units) as units, sum(units*size) as quantity,
				sum(price/prmult*units)/sum(units*size) as unit_price, 
				sum(price/prmult/size*vol_weight)/sum(vol_weight) as unit_price_fixwt, 
				sum(size_index*units*size)/sum(units*size) as size_index,
				count(unique(brand_descr)) as num_brand, 
				count(unique(upc)) as num_upc, 
				sum(prvt_flag) as num_prvtlb
		from joint_move1 
		group by store_code_uc, week_end, product_module_code
		order by store_code_uc, week_end, product_module_code;
	quit; 
	
	proc datasets noprint; delete joint_move1; run;
	
/*	* Export data; 
	%let outfile = %sysfunc(catx(%str(), &dirname,\Group_year\&sel_year\,&sel_group,_, &sel_year, .csv));
	%put &outfile;
	PROC EXPORT DATA= store_week_mod
	            OUTFILE= "&outfile"
	            DBMS=csv REPLACE;
	     PUTNAMES=YES;
	RUN;*/
%mend gen_store_module_week; 

*-----------------------------*;
* Generate UPC national sales *; 
* Aggregate UPC sales over stores for a given week;
%macro gen_upcweek(joint_move); 
	proc sql noprint;
		create table upc_sale as 
		select *, revenue/sum(revenue) as mkt_share
		from (select upc, upc_ver_uc, week_end, product_module_code, 
				sum(price/prmult*units) as revenue, sum(units) as units,
				sum(price/prmult*units)/sum(units) as price, 
				count(unique(store_code_uc)) as num_store
				from &joint_move
				group by upc, upc_ver_uc, week_end, product_module_code) 
		group by product_module_code
		order by upc, upc_ver_uc, week_end;
	quit;
	
	* Export data; 
	%let outfile = %sysfunc(catx(%str(), &dirname,\Newprod\UPC_sales\&sel_year\,&sel_group,_, &sel_year, .csv));
	%put &outfile;
	PROC EXPORT DATA= upc_sale
	            OUTFILE= "&outfile"
	            DBMS=csv REPLACE;
	     PUTNAMES=YES;
	RUN;
%mend gen_upcweek; 

******************************************************************;
* Run the functions for a single product group and a single year *;
******************************************************************;
%gen_upcwt(append_data = 0); 
%gen_store_module_week(joint_move);
%gen_upcweek(joint_move); 



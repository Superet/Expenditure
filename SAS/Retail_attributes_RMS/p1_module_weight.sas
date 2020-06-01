/*This file compute the module wallet share.
Wallet share is computed using the dollar revenue from all store in 2006. 

NOTE: 
1. There are 3 modules that have no sales in 2006 (5520, 5521, 9599). 
2. MAGNET DATA has product_group_code as 0 in 2006.

*/

libname mymaster "/sscc/home/c/ccv103/rms_attributes/master";
libname	grpyear "/sscc/home/c/ccv103/rms_attributes/group_year";
libname	newprod "/sscc/home/c/ccv103/rms_attributes/newprod";
%let dirname= /sscc/home/c/ccv103/rms_attributes;

* Product group index data; 
data prod_group_index;
	set mymaster.prod_group_index(drop = n);
	where product_group_code not in (5520, 5521, 9599); 
	if product_group_code=99 then product_group_code=9999; 
	n = _N_; 
run;
proc sql noprint;
	select max(n) into: num_grp from prod_group_index; 
quit;
%put &num_grp; 

* Compute the annual dollar sales of each module; 
* Loop over modules; 
%macro agrg_mod; 
	%do i = 1 %to &num_grp; 
		proc sql noprint;
			select put(product_group_code, z4.) into: sel_group from prod_group_index
			where n = &i; 
		quit;
		%let outf = grpyear.swm_&sel_group._2006; 
		%put &outf;
		
		proc sql noprint; 
			create table grp&i as
			select product_module_code, sum(revenue) as revenue
			from &outf
			group by product_module_code
			order by product_module_code;
		quit; 
	%end; 
%mend agrg_mod; 

%agrg_mod; 

* Stack together to compute moduel weight; 
%let outf = %sysfunc(catx(%str(),grp,&num_grp));
%put &outf;
data module_wt;
	set grp1-&outf; 
run;
proc sql noprint;
	create table mymaster.module_weight as
	select *, revenue/sum(revenue) as wallet_share
	from module_wt
	order by product_module_code; 
quit;
proc contents data = mymaster.module_weight; run;

/*
This file aggregates channel-level retail attributes. 
Input data: 
	master/module_weight
	master/products	
	Group_year/swm_XXXX_YYYY.csv
	master/store
Output data: 
	store_rev_year.csv
	channel_year.csv
	
*/

options compress=yes reuse=yes; 
libname mymaster "/sscc/home/c/ccv103/rms_attributes/master";
libname	grpyear "/sscc/home/c/ccv103/rms_attributes/group_year";
libname	newprod "/sscc/home/c/ccv103/rms_attributes/newprod";
%let dirname= /sscc/home/c/ccv103/rms_attributes;
%let sel_year	= %sysget(sel_year); 
%put &sel_year; 

/*%let sel_year = 2006; */

* Prepare moduel data: wallet share and module size; 
proc sql noprint;
	create table module_dat as 
	select B.product_module_code, wallet_share, module_size
	from mymaster.module_weight as A
	full join (select product_module_code, mean(module_size) as module_size
				from mymaster.products group by product_module_code) as B
	on A.product_module_code = B.product_module_code;
quit; 

* Read store data; 
%let fpath = &dirname/master/stores_both.csv;
%put &fpath;
PROC IMPORT OUT= stores_both
            DATAFILE= "&fpath" 
            DBMS=csv REPLACE;
     		GETNAMES=YES;
			GUESSINGROWS=5000; 
RUN;
* Drop all the liqure stoes;
data stores_both;
	set stores_both; 
	where channel_type ^= "";
run;

* Obstain file names from Group_year folder; 
%let subdir	= &dirname/group_year;
%put &subdir; 
filename dirlist pipe "ls &subdir";
data dirlist;
	length fname $30;
	infile dirlist length=reclen;
	input fname $varying30. reclen;
	product_group_code = substr(fname, 5, 4); 
	year = substr(fname, 10, 4); 
	if year ^= &sel_year then delete; 
	n = _N_;
run;	

proc sql noprint; 
	select max(n) into: num_group from dirlist; 
quit; 
%put &num_group;

****************************;
* Loop over product groups *; 
****************************;
%macro agrg_store_group; 
	%do i=1 %to &num_group;
		* Select product_group_code; 
		proc sql noprint;	
			select product_group_code into: sel_group from dirlist
			where n = &i; 
		quit;
		%put &sel_group; 
		%let indat = grpyear.swm_%trim(&sel_group)_&sel_year; 
		%put &indat;
	
		proc sql noprint; 
		* Merge store-moduel data with module summary data; 		
			create table tmp as 
			select A.*, B.wallet_share, B.module_size
			from &indat as A left join module_dat as B
			on A.product_module_code = B.product_module_code; 
			
		* Sum up partial basket price and other retail attributes from the product_group;
			create table grp&i as 
			select store_code_uc, week_end,
				sum(wallet_share) as sum_wallet_share,
				sum(revenue/quantity*module_size*wallet_share) as pt_bskt_price, 
				sum(unit_price_fixwt*module_size*wallet_share) as pt_bskt_price_fixwt,
				sum(size_index*wallet_share) as pt_size_index,
				count(unique(product_module_code)) as pt_num_module,
				sum(num_brand) as pt_num_brand, 
				sum(num_upc) as pt_num_upc,
				sum(wallet_share*num_brand) as pt_num_brand_wt, 
				sum(wallet_share*num_upc) as pt_num_upc_wt,
				sum(num_prvtlb) as pt_num_prvtlb,
				sum(num_prvtlb/num_upc*wallet_share) as pt_prvt_per_mod
			from tmp 
			group by store_code_uc, week_end;

		* Compute annual revenue of the product_group; 
			create table store&i as 
			select store_code_uc, &sel_year as year, sum(revenue) as pt_revenue
			from tmp 
			group by store_code_uc; 
			
		* Delet joint table; 
			drop table tmp; 
		quit; 
	%end; 
%mend agrg_store_group;

%agrg_store_group;

***************;
* Aggregation *; 
***************;
*-------------------------*; 
* Summarize store revenue *; 
* Stack all the revenue data; 
%let lastd = %sysfunc(catx(%str(), store, &num_group)); 
%put &lastd;
data store_tmp; 
	set store1-&lastd;
run;
proc datasets noprint; delete store1-&lastd; run;

* Compute store revenue; 
proc sql noprint;
	create table store_rev as 
	select store_code_uc, year, sum(pt_revenue) as revenue
	from store_tmp
	group by store_code_uc, year
	order by store_code_uc; 
	
	create table store_rev_joint as 
	select A.*, channel_type, dma_descr, scantrack_market_descr
	from store_rev as A inner join stores_both as B
	on A.store_code_uc = B.store_code_uc and A.year = B.year;
quit; 

proc datasets noprint; delete store_tmp; run;

*-----------------------------*;
* Summarize retail attributes *;
* Stack all the data of retail attributes from product_groups; 
%let lastd = %sysfunc(catx(%str(), grp, &num_group)); 
%put &lastd;
data grp_tmp; 
	set grp1-&lastd;
run;
proc datasets noprint; delete grp1-&lastd; run;

* Compute store-level retail attributes; 
proc sql noprint;
	create table store_attr as 
	select store_code_uc, week_end, 
		sum(pt_bskt_price)/sum(sum_wallet_share) as bskt_price, 
		sum(pt_bskt_price_fixwt)/sum(sum_wallet_share) as bskt_price_fixwt, 
		sum(pt_size_index)/sum(sum_wallet_share) as size_index, 
		sum(pt_num_module) as num_module, 
		sum(pt_num_brand) as num_brand, 
		sum(pt_num_upc) as num_upc, 
		sum(pt_num_brand_wt)/sum(sum_wallet_share) as num_brand_per_mod, 
		sum(pt_num_upc_wt)/sum(sum_wallet_share) as num_upc_per_mod, 
		sum(pt_num_prvtlb)/sum(pt_num_upc) as overall_prvt,
		sum(pt_prvt_per_mod)/sum(sum_wallet_share) as prvt_per_mod
	from grp_tmp
	group by store_code_uc, week_end
	order by store_code_uc, week_end; 
	
	drop table grp_tmp; 
	
* Merge store-level retail attributes and store revenue; 
	create table store_attr_joint as 
	select A.*, revenue, channel_type, scantrack_market_descr, dma_descr
	from store_attr as A inner join store_rev_joint as B
	on A.store_code_uc = B.store_code_uc;
	
	drop table store_attr; 

* Aggregate channel-level retail attributes: scantrack_market_code-channel_type-week; 
	create table channel_attr_scantrack as 
	select scantrack_market_descr, channel_type, week_end, 
		sum(bskt_price*revenue)/sum(revenue) as bskt_price, 
		sum(bskt_price_fixwt*revenue)/sum(revenue) as bskt_price_fixwt, 
		sum(size_index*revenue)/sum(revenue) as size_index, 
		sum(num_module*revenue)/sum(revenue) as num_module, 
		sum(num_brand*revenue)/sum(revenue) as num_brand, 
		sum(num_upc*revenue)/sum(revenue) as num_upc, 
		sum(num_brand_per_mod*revenue)/sum(revenue) as num_brand_per_mod, 
		sum(num_upc_per_mod*revenue)/sum(revenue) as num_upc_per_mod, 
		sum(overall_prvt*revenue)/sum(revenue) as overall_prvt, 
		sum(prvt_per_mod*revenue)/sum(revenue) as prvt_per_mod
	from store_attr_joint
	group by scantrack_market_descr, channel_type, week_end
	order by scantrack_market_descr, channel_type, week_end;

* Aggregate channel-level retail attributes: dma-channel_type-week; 
	create table channel_attr_dma as 
	select dma_descr, channel_type, week_end, 
		sum(bskt_price*revenue)/sum(revenue) as bskt_price, 
		sum(bskt_price_fixwt*revenue)/sum(revenue) as bskt_price_fixwt, 
		sum(size_index*revenue)/sum(revenue) as size_index, 
		sum(num_module*revenue)/sum(revenue) as num_module, 
		sum(num_brand*revenue)/sum(revenue) as num_brand, 
		sum(num_upc*revenue)/sum(revenue) as num_upc, 
		sum(num_brand_per_mod*revenue)/sum(revenue) as num_brand_per_mod, 
		sum(num_upc_per_mod*revenue)/sum(revenue) as num_upc_per_mod, 
		sum(overall_prvt*revenue)/sum(revenue) as overall_prvt, 
		sum(prvt_per_mod*revenue)/sum(revenue) as prvt_per_mod
	from store_attr_joint
	group by dma_descr, channel_type, week_end
	order by dma_descr, channel_type, week_end;
	
	drop table store_attr_joint, store_rev_joint; 
quit;		

**********************;		
* Export the results *;
**********************;
%let outfile = %sysfunc(catx(%str(), &dirname,/annual_files/channel_dma/,channel_attr_,&sel_year,.csv));
%put &outfile;
PROC EXPORT DATA= channel_attr_dma
            OUTFILE= "&outfile"
            DBMS=csv REPLACE;
     PUTNAMES=YES;
RUN;

%let outfile = %sysfunc(catx(%str(), &dirname,/annual_files/channel_scantrack/,channel_attr_,&sel_year,.csv));
%put &outfile;
PROC EXPORT DATA= channel_attr_scantrack
            OUTFILE= "&outfile"
            DBMS=csv REPLACE;
     PUTNAMES=YES;
RUN;

%let outfile = %sysfunc(catx(%str(), &dirname,/annual_files/store_revenue/,store_revenue_,&sel_year,.csv));
%put &outfile;
PROC EXPORT DATA= store_rev
            OUTFILE= "&outfile"
            DBMS=csv REPLACE;
     PUTNAMES=YES;
RUN;
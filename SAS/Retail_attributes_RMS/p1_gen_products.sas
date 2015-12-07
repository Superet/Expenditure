libname mymaster "E:\Users\ccv103\Documents\Research\Nielsen\RMS_attributes\Master";

%let dirname= E:\Users\ccv103\Documents\Research\Nielsen\nielsen_extracts\RMS;
options compress=yes reuse=yes;

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

* Add magnet indicator; 
data my_products;
	set orig_products;
	magnet = 0; 
	if department_descr = "MAGNET DATA" then magnet = 1;
	if product_module_code >= 445 & product_module_code <= 468 then magnet = 1; 
run;

* ------------------------------------------------------------------------------------------ *;
* In the original data, the unit measurements (size1_units) are not unique for a single module;
* Assign a unit measurement for each product module; 
proc sort data = my_products; by department_descr product_group_descr product_module_descr;
proc freq data = my_products noprint; 
	by department_descr product_group_descr product_module_descr;
	tables size1_units/ out=tmp;
run;

proc sort data=tmp; by department_descr product_group_descr product_module_descr COUNT;
data module_measure;
	set tmp;
	by department_descr product_group_descr product_module_descr; 
	if last.product_module_descr then output;
run;

* ------------------------------------------------------------------------------------------ *;
* Imputing size for the UPCs that have different unit measurement from the majority UPCs in the category; 
* Use medium unit measurement to replace 1 CT;
* First merge in the main unit measure for each module; 
proc sql noprint;
	create table tmp as 
	select A.*, A.size1_amount*A.multi as size, B.size1_units as module_units
	from my_products as A left join module_measure as B
	on A.product_module_descr = B.product_module_descr
	order by department_descr, product_group_descr, product_module_descr;
quit;

* Check what percent of products have different unit measure;
proc sql; select sum(size1_units ^= module_units)/count(size1_units) from tmp; quit;

* Compute the median size of each category;
proc univariate data=tmp(where = (size1_units=module_units)) noprint;
	by department_descr product_group_descr product_module_descr;
	var size;
	output out = basket MEDIAN = module_size;
run;

* Merge the product data with median size for other unit measurements;
proc sql noprint;
	create table tmp2 as 
	select A.*, B.module_size
	from tmp as A left join basket as B
	on A.product_module_descr = B.product_module_descr
	order by upc;
quit;

* If the unit measurement of a product is not the main measurement of the module, then we assign the medium size with known measurement; 
data tmp2;
	set tmp2;
	if size1_units ^= module_units then size = module_size;
run;

* Compute relative size index; 
data products;
	retain upc upc_ver_uc 
		department_code department_descr product_group_code product_group_descr product_module_code product_module_descr magnet
		brand_descr module_units module_size size size_index prvt_flag; 
	set tmp2 (keep = 	upc upc_ver_uc 
			department_code department_descr product_group_code product_group_descr product_module_code product_module_descr magnet
			brand_descr module_units module_size size); 
	size_index = size/module_size; 
	PRVT_flag = 0; 
	if brand_descr="CTL BR" then PRVT_flag = 1;
	label	module_units 	= "Unit measure of product module"
			module_size 	= "The median size of module measured at unit measure"
			size_index 		= "Size normalized by median size within module";
run;
proc contents data = products; run;	

data mymaster.products; set products; run; 

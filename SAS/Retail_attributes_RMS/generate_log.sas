%let dirname= E:\Users\ccv103\Documents\Research\Nielsen\nielsen_extracts\HMS;
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

* Table of UPC by product group;
proc sql noprint;
	create table prod_group_index as 
	select department_code, product_group_code, product_group_descr, 
		count(unique(product_module_descr)) as num_module,
		count(unique(upc)) as num_upc
	from mymaster.products
	group by department_code, product_group_code, product_group_descr
	order by department_code, product_group_code, product_group_descr;
quit; 

data mymaster.prod_group_index; 
	set prod_group_index; 
	where department_code is not missing; 
	n = _N_; 
run; 
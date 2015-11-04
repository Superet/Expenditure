drop table dbo.tmp, dbo.mydata, dbo.mydata1; 

select trip_code_uc, household_code, purchase_date, year, C.upcv, quantity, price, product_module_descr, size_index, size_index*quantity as size
into dbo.tmp
from (
	select A.trip_code_uc, household_code, B.year, purchase_date, upcv, quantity, price
	from (select trip_code_uc, upcv, quantity,price from dbo.purchases) as A inner join 
	dbo.trips as B on A.trip_code_uc = B.trip_code_uc
)as C inner join dbo.products as D on C.upcv = D.upcv; 

/*Run panelist ; */
select A.*, B.famsize, B.income_real, B.first_income, log(B.income_midvalue) as ln_income, 
	first_incomeg = CASE WHEN first_income <= 17 THEN 'T1'
						 WHEN first_income > 17 and first_income <= 23 THEN 'T1'
						 WHEN first_income > 23 THEN 'T3'
					END
	into dbo.mydata
	from dbo.tmp as A inner join dbo.panelists as B
	on A.household_code = B.household_code and A.year = B.panel_year
	order by household_code, product_module_descr;

/*Restrict to storable goods; */
select *
into dbo.mydata1
from dbo.mydata
where product_module_descr in (select product_module_descr from dbo.tmp_products where storable = 1)
order by household_code; 

* In SAS; 
proc sql noprint;
	create table mydata1 as
	select *
	from mysqllib.mydata1
	order by household_code; 
quit;

proc glm data = mydata1;
	absorb household_code; 
	class product_module_descr famsize first_incomeg; 
	model size = ln_income*first_incomeg price product_module_descr/ solution; 
	ods output ParameterEstimates=est_out;
run; 


proc glm data = mydata1;
	absorb household_code; 
	class product_module_descr year first_incomeg; 
	model size = year*first_incomeg price product_module_descr/ solution; 
	ods output ParameterEstimates=est_out;
run; 



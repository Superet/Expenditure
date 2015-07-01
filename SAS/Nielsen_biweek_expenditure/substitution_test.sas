libname mylib "E:\Users\ccv103\Documents\Research\Store switching\SAS_temp\Nielsen";
libname mysel "E:\Users\ccv103\Documents\Research\Store switching\SAS_temp\Nielsen\Work";
libname mainex "E:\Users\ccv103\Documents\Research\Store switching\SAS_temp\Nielsen\Work\main";

options compress=yes reuse=yes; 
proc datasets library=WORK kill noprint; run;

%let dirname= E:\Users\ccv103\Documents\Research\Nielsen\nielsen_extracts\HMS;

* The output files; 

*************;
* Functions *; 
*************;
* A function that runs fixed-effect regression; 
%macro run_FE(data, dv, iv, classvar, FEvar, outname);
	proc glm data = &data; 
		absorb &FEvar; 
		class &classvar; 
		model &dv = &iv/solution; 
		ods output ParameterEstimates=&outname; 
	run;
%mend run_FE; 

* A function that merge and delete additional data; 
%macro my_merge(base_data, add_data, byvar, outname, delete_add = 1);
	proc sort data = &add_data; by &byvar; run; 
	data &outname; 
		merge &base_data &add_data; 
		by &byvar; 
		proc sort; by &byvar; 
	run;
	
	%if &delete_add = 1 %then %do; 
		proc datasetes noprint; delete &add_data; run; 
	%end; 
%mend my_merge; 

*****************; 
* Organize data *; 
*****************; 
/*
Four main datasets: 
1. demo_init: the initial states of households. KEY (household_code)
2. demo_annual: the annual profiles of households. KEY (household_code, year)
3. hh_exp: the biweekly expenditure data of households. KEY (household_code, biweek)
4. purchase: individual purchases of products. KEY (household_code, date, UPC)
*/
data demo_annual; 
	set mylib.panelists(keep = household_code year income fam); 
run;

* Segment households using intial states;
proc sort data = demo_annual; by household_code year; run;
data demo_init(drop = income_group); 
	set demo_annual; 
	by household_code; 
	if not first.household_code then delete; 
	first_income = income_group; 
run;

*--------------------------------------------------------*; 
* Compute percentile of outcome variables to hh_exp -- expenditure and share; 


* Append initial states to hh_exp; 


*--------------------------------------------------------*; 
* Compute percentile of outcome variables in purchase data -- price tier and size tier; 
* Price tier of UPC within category; 
proc sql noprint;
	* Average unit price of UPC; 
	
quit; 

* Rank price within category; 

* Append price tier to purchase data; 



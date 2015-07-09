libname mylib "E:\Users\ccv103\Documents\Research\Store switching\SAS_temp\Nielsen";
libname mysel "E:\Users\ccv103\Documents\Research\Store switching\SAS_temp\Nielsen\Work";
libname mainex "E:\Users\ccv103\Documents\Research\Store switching\SAS_temp\Nielsen\Work\main";

options compress=yes reuse=yes; 
proc datasets library=WORK kill noprint; run;

%let dirname= E:\Users\ccv103\Documents\Research\Nielsen\nielsen_extracts\HMS;

* The output files; 
%let reg_outfile = %sysfunc(catx(%str(), E:\Users\ccv103\Desktop\substitution_test_,&sysdate9.,.csv));
%put &reg_outfile; 


*************;
* Functions *; 
*************;
* A function that runs fixed-effect regression; 
%macro run_FE(data, dv, iv, classvar, FEvar, outname);
	proc glm data = &data order = data; 
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
		proc datasets noprint; delete &add_data; run; 
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

* Read in CPI data;
PROC IMPORT OUT= CPI
            DATAFILE= "\\tsclient\Resear1\Store switching\processed data\CPI.csv" 
            DBMS=csv REPLACE;
     		GETNAMES=YES; 
RUN;

proc sql noprint;
	select annual into:cpi_base from cpi where year = 2004; 
quit; 
%put &cpi_base; 

data cpi(keep = year annual cpi); 
	set cpi; 
	cpi = annual/&cpi_base; 
run;

* Household demographics; 
data demo_annual; 
	set mylib.panelists(keep = household_code panel_year income_group income_real income_midvalue famsize scantrack_market_descr); 
run;

*Select only 5% households to test the codes; 
/*proc sql noprint;
	create table tmp as 
	select *
	from (select unique(household_code) from demo_annual)
	where RANUNI(4321) < .05; 
	
	create table tmp1 as 
	select *
	from demo_annual 
	where household_code in (select * from tmp); 
quit; */
data tmp1; set demo_annual; run; 

* Adjust income with CPI; 
proc sql noprint; 
	create table demo_annual as 
	select A.*, income_midvalue/cpi as income_cpi
	from tmp1 as A left join cpi as B
	on A.panel_year = B.year; 
quit; 

* Segment households using intial states;
proc sort data = demo_annual; by household_code panel_year; run;

data demo_init; 
	set demo_annual; 
	by household_code; 
	if not first.household_code then delete; 
run;

* The distribution of intial income; 
proc freq data = demo_init; table income_real; run; 
data demo_init(drop = income_group);
	set demo_init; 
	first_income = 'T1'; 
	if income_real >16 & income_real <=21 then first_income = 'T2'; 
	if income_real >21 then first_income = 'T3'; 
run; 

*--------------------------------------------------------*; 
* Compute percentile of outcome variables to hh_exp -- expenditure and share; 
proc sql noprint;
	create table tmp as 
	select household_code, biweek, year, dol_purchases, 
		DOLP_Convenience_Store, DOLP_Discount_Store, DOLP_Dollar_Store, 
		DOLP_Drug_Store, DOLP_Grocery, DOLP_Warehouse_Club,
		DOLP_Convenience_Store/dol_purchases as SHR_Convenience_Store, DOLP_Discount_Store/dol_purchases as SHR_Discount_Store, 
		DOLP_Dollar_Store/dol_purchases as SHR_Dollar_Store, DOLP_Drug_Store/dol_purchases as SHR_Drug_Store, 
		DOLP_Grocery/dol_purchases as SHR_Grocery, DOLP_Warehouse_Club/dol_purchases as SHR_Warehouse_Club,
		num_day, num_trip, num_food_module, num_noned_module, food_quant, nonedible_quant, 
		scantrack_market_descr, income_real, income_midvalue
	from mylib.Hh_biweek_exp_merge
	where household_code in (select household_code from demo_init); 
quit; 

proc rank data=tmp out=tmp_rank ties=mean PERCENT;
	var dol_purchases SHR_Convenience_Store SHR_Discount_Store SHR_Dollar_Store
		SHR_Drug_Store SHR_Grocery SHR_Warehouse_Club;
	ranks PCTL_DOL PCTL_Convenience_Store PCTL_Discount_Store PCTL_Dollar_Store
		PCTL_Drug_Store PCTL_Grocery PCTL_Warehouse_Club;
run;

* Append initial states to hh_exp and adjust income with cpi;
proc sql noprint;
	create table tmp as 
	select A.*, B.first_income, B.famsize
	from tmp_rank as A inner join demo_init as B
	on A.household_code = B.household_code;
	
	create table hh_exp as 
	select A.*, income_midvalue/B.cpi as income_cpi
	from tmp as A left join cpi as B
	on A.year = B.year
	order by household_code, biweek; 
quit; 

data hh_exp;
	set hh_exp; 
	first_date 	= mdy(1,1,2004) + (biweek-1)*14; 
	month 		= month(first_date);
	ymonth 		= cats(year, '-', month); 
	dol_rt 		= dol_purchases/(income_cpi/26);
	quant 		= food_quant + nonedible_quant; 
	num_mod		= num_food_module + num_noned_module; 
	ln_y		= .; 
	if dol_purchases > 0 then ln_y		= log(dol_purchases); 
	drop first_date;
run; 

*--------------------------------------------------------*; 
* Compute percentile of outcome variables in purchase data -- price tier and size tier; 
data dict; 
	set mylib.products(keep = upcv department_descr product_group_descr product_module_descr
		 	brand_descr multi size size1_amount size_percentile ); 
run; 

* Price tier of UPC within category; 
proc sql noprint;
	* Subset purchase data; 
	create table tmp_purch as 
	select household_code, year, purchase_date, month(purchase_date) as month, biweek, upcv, product_module_descr, 
			size, quantity, price, price/size as unit_price
	from (select * from mainex.purchases where household_code in (select household_code from demo_init))
	order by product_module_descr; 

	* Average unit price of UPC; 	
	create table tmp as 
	select product_module_descr, upcv, mean(unit_price) as unit_price, mean(size) as size 
	from tmp_purch
	group by product_module_descr, upcv;
quit; 

* Rank price within category; 
proc rank data=tmp out=myprod ties=mean PERCENT;
	by product_module_descr;
	var unit_price size;
	ranks PCTL_UPRICE PCTL_SIZE; 
run;

* Append price tier and household segment to purchase data; 
proc sql noprint;
	create table tmp as 
	select A.*, B.PCTL_UPRICE, B.PCTL_SIZE
	from tmp_purch as A left join myprod as B
	on A.upcv = B.upcv; 

	create table purchases as 
	select A.*, B.first_income, B.famsize
	from tmp as A left join demo_init as B
	on A.household_code = B.household_code
	order by household_code, purchase_date; 
quit; 

* Aggregate size tier and price tier from purchase data to biweek level; 
proc sql noprint;
	* Combine multiple purchases from the same category; 
	create table tmp as 
	select household_code, biweek, product_module_descr, 
		sum(quantity*PCTL_UPRICE)/sum(quantity) as PCTL_UPRICE, 
		sum(quantity*PCTL_SIZE)/sum(quantity) as PCTL_SIZE
	from purchases
	group by household_code, biweek, product_module_descr; 
	
	* Append category weights;
	create table tmp1 as 
	select A.*, B.share_paid
	from tmp as A left join (select product_module_descr, mean(share_paid) as share_paid from mylib.module_wallet group by product_module_descr) as B
	on A.product_module_descr = B.product_module_descr;
	
	* Compute weighted average of size/price tier during biweek; 
	create table tmp2 as 
	select household_code, biweek, sum(PCTL_UPRICE*share_paid)/sum(share_paid) as PCTL_UPRICE, 
			sum(PCTL_SIZE*share_paid)/sum(share_paid) as PCTL_SIZE
	from tmp1 
	group by household_code, biweek; 
	
	* Add price/size tier to hh_exp; 
	create table tmp3 as 
	select A.*, B.PCTL_UPRICE, B.PCTL_SIZE
	from hh_exp as A left join tmp2 as B
	on A.household_code = B.household_code and A.biweek = B.biweek; 
quit;
data hh_exp; set tmp3; run; 

proc datasets noprint; delete tmp tmp1 tmp2 tmp3 tmp_rank; run; 

*****************************;  
* Table of raw distribution *; 
*****************************; 
%let Rc_year = 2008; 
proc means data = hh_exp(where = (year < &Rc_year));
	class famsize first_income;
	var PCTL_DOL PCTL_Discount_Store PCTL_Grocery PCTL_Warehouse_Club;
	output out = tmp1; 
run;

proc means data=purchases(where = (year < &Rc_year)); 
	class famsize first_income; 
	var PCTL_UPRICE PCTL_SIZE; 
	output out = tmp2; 
run; 

proc sort data = tmp1; 	by famsize first_income _STAT_; 
proc sort data = tmp2; 	by famsize first_income _STAT_; 
run; 
data raw_mean; 
	merge tmp1 tmp2; 
	by famsize first_income _STAT_; 
run;

proc datasets noprint; delete tmp1 tmp2; run;

**********************************************************************;
* The position of outcome variable by segment prior to the recession *; 
**********************************************************************;
proc sort data = hh_exp; by month descending first_income descending famsize descending year; run; 

* The position within expenditure distribution by segment; 
%run_FE(data=hh_exp(where = (year < &Rc_year)), dv=PCTL_DOL, iv= first_income famsize year, 
		classvar = first_income famsize year, FEvar = month, outname = tmp);
data tmp; length test $30 Dependent $50 Parameter $50 Parameter $50 ; set tmp; test = 'Pre-recession'; run; 
data myreg; set tmp; run;

* Distribution of expenditure share; 
%run_FE(data=hh_exp(where = (year < &Rc_year)), dv=PCTL_Discount_Store, iv= first_income famsize year, 
		classvar = first_income famsize year, FEvar = month, outname = tmp);
data tmp; length test $30 Dependent $50 Parameter $50 ; set tmp; test = 'Pre-recession'; run; 
%my_merge(base_data=myreg, add_data=tmp, byvar=test Dependent, outname = myreg);

%run_FE(data=hh_exp(where = (year < &Rc_year)), dv=PCTL_Grocery, iv= first_income famsize year, 
		classvar = first_income famsize year, FEvar = month, outname = tmp);
data tmp; length test $30 Dependent $50 Parameter $50 ; set tmp; test = 'Pre-recession'; run; 
%my_merge(base_data=myreg, add_data=tmp, byvar=test Dependent, outname = myreg);

%run_FE(data=hh_exp(where = (year < &Rc_year)), dv=PCTL_Warehouse_Club, iv= first_income famsize year, 
		classvar = first_income famsize year, FEvar = month, outname = tmp);
data tmp; length test $30 Dependent $50 Parameter $50 ; set tmp; test = 'Pre-recession'; run; 
%my_merge(base_data=myreg, add_data=tmp, byvar=test Dependent, outname = myreg);

* Price tier; 
proc sort data = purchases; by product_module_descr descending first_income descending year descending month; 

%run_FE(data=purchases(where = (year < &Rc_year)), dv=PCTL_UPRICE, iv= first_income famsize year month, 
		classvar = first_income famsize year month, FEvar = product_module_descr, outname = tmp);
data tmp; length test $30 Dependent $50 Parameter $50 ; set tmp; test = 'Pre-recession'; run; 
%my_merge(base_data=myreg, add_data=tmp, byvar=test Dependent, outname = myreg);

* Size tier; 
%run_FE(data=purchases(where = (year < &Rc_year)), dv=PCTL_SIZE, iv= first_income famsize year month, 
		classvar = first_income famsize year month, FEvar = product_module_descr, outname = tmp);
data tmp; length test $30 Dependent $50 Parameter $50 ; set tmp; test = 'Pre-recession'; run; 
%my_merge(base_data=myreg, add_data=tmp, byvar=test Dependent, outname = myreg);

***************************************;
* The cost of substitution by segment *; 
***************************************;
* The cost of moving expenditure percentile up -- expenditure to income ratio; 
%run_FE(data=hh_exp, dv=dol_purchases, iv= first_income*PCTL_DOL famsize year, 
		classvar = first_income famsize year, FEvar = month, outname = tmp);
data tmp; length test $30 Dependent $50 Parameter $50 ; set tmp; test = 'Cost-expenditurePCT'; run; 
%my_merge(base_data=myreg, add_data=tmp, byvar=test Dependent, outname = myreg);

%run_FE(data=hh_exp, dv=dol_rt, iv= first_income*PCTL_DOL famsize year, 
		classvar = first_income famsize year, FEvar = month, outname = tmp);
data tmp; length test $30 Dependent $50 Parameter $50 ; set tmp; test = 'Cost-expenditurePCT'; run; 
%my_merge(base_data=myreg, add_data=tmp, byvar=test Dependent, outname = myreg);

* The cost of moving expenditure down -- quantity and variety loss; 
%run_FE(data=hh_exp, dv=quant, iv= first_income*PCTL_DOL famsize year, 
		classvar = first_income famsize year, FEvar = month, outname = tmp);
data tmp; length test $30 Dependent $50 Parameter $50 ; set tmp; test = 'Cost-expenditurePCT'; run; 
%my_merge(base_data=myreg, add_data=tmp, byvar=test Dependent, outname = myreg);

%run_FE(data=hh_exp, dv=num_mod, iv= first_income*PCTL_DOL famsize year, 
		classvar = first_income famsize year, FEvar = month, outname = tmp);
data tmp; length test $30 Dependent $50 Parameter $50 ; set tmp; test = 'Cost-expenditurePCT'; run; 
%my_merge(base_data=myreg, add_data=tmp, byvar=test Dependent, outname = myreg);

*----------------------------------------------------------------*; 
* The cost of moving expenditure up -- expenditure to income ratio; 
%run_FE(data=hh_exp, dv=dol_rt, iv= first_income*ln_y famsize year, 
		classvar = first_income famsize year, FEvar = month, outname = tmp);
data tmp; length test $30 Dependent $50 Parameter $50 ; set tmp; test = 'Cost-expenditure'; run; 
%my_merge(base_data=myreg, add_data=tmp, byvar=test Dependent, outname = myreg);

%run_FE(data=hh_exp, dv=dol_purchases, iv= first_income*ln_y famsize year, 
		classvar = first_income famsize year, FEvar = month, outname = tmp);
data tmp; length test $30 Dependent $50 Parameter $50 ; set tmp; test = 'Cost-expenditure'; run; 
%my_merge(base_data=myreg, add_data=tmp, byvar=test Dependent, outname = myreg);

* The cost of moving expenditure down -- quantity and variety loss; 
%run_FE(data=hh_exp, dv=quant, iv= first_income*ln_y famsize year, 
		classvar = first_income famsize year, FEvar = month, outname = tmp);
data tmp; length test $30 Dependent $50 Parameter $50 ; set tmp; test = 'Cost-expenditure'; run; 
%my_merge(base_data=myreg, add_data=tmp, byvar=test Dependent, outname = myreg);

%run_FE(data=hh_exp, dv=num_mod, iv= first_income*ln_y famsize year, 
		classvar = first_income famsize year, FEvar = month, outname = tmp);
data tmp; length test $30 Dependent $50 Parameter $50 ; set tmp; test = 'Cost-expenditure'; run; 
%my_merge(base_data=myreg, add_data=tmp, byvar=test Dependent, outname = myreg);

*----------------------------------------------------------------*; 
* The cost of moving expenditure share up -- expenditure to income ratio; 
%run_FE(data=hh_exp, dv=dol_rt, iv= first_income*SHR_Discount_Store famsize year, 
		classvar = first_income famsize year, FEvar = month, outname = tmp);
data tmp; length test $30 Dependent $50 Parameter $50 ; set tmp; test = 'Cost-Discount_Store'; run; 
%my_merge(base_data=myreg, add_data=tmp, byvar=test Dependent, outname = myreg);

%run_FE(data=hh_exp, dv=dol_purchases, iv= first_income*SHR_Discount_Store famsize year, 
		classvar = first_income famsize year, FEvar = month, outname = tmp);
data tmp; length test $30 Dependent $50 Parameter $50 ; set tmp; test = 'Cost-Discount_Store'; run; 
%my_merge(base_data=myreg, add_data=tmp, byvar=test Dependent, outname = myreg);

* The cost of moving expenditure down -- quantity and variety loss; 
%run_FE(data=hh_exp, dv=quant, iv= first_income*SHR_Discount_Store famsize year, 
		classvar = first_income famsize year, FEvar = month, outname = tmp);
data tmp; length test $30 Dependent $50 Parameter $50 ; set tmp; test = 'Cost-Discount_Store'; run; 
%my_merge(base_data=myreg, add_data=tmp, byvar=test Dependent, outname = myreg);

%run_FE(data=hh_exp, dv=num_mod, iv= first_income*SHR_Discount_Store famsize year, 
		classvar = first_income famsize year, FEvar = month, outname = tmp);
data tmp; length test $30 Dependent $50 Parameter $50 ; set tmp; test = 'Cost-Discount_Store'; run; 
%my_merge(base_data=myreg, add_data=tmp, byvar=test Dependent, outname = myreg);

*----------------------------------------------------------------*; 
* The cost of moving expenditure up -- expenditure to income ratio; 
%run_FE(data=hh_exp, dv=dol_rt, iv= first_income*SHR_Grocery famsize year, 
		classvar = first_income famsize year, FEvar = month, outname = tmp);
data tmp; length test $30 Dependent $50 Parameter $50 ; set tmp; test = 'Cost-Grocery'; run; 
%my_merge(base_data=myreg, add_data=tmp, byvar=test Dependent, outname = myreg);

%run_FE(data=hh_exp, dv=dol_purchases, iv= first_income*SHR_Grocery famsize year, 
		classvar = first_income famsize year, FEvar = month, outname = tmp);
data tmp; length test $30 Dependent $50 Parameter $50 ; set tmp; test = 'Cost-Grocery'; run; 
%my_merge(base_data=myreg, add_data=tmp, byvar=test Dependent, outname = myreg);

* The cost of moving expenditure down -- quantity and variety loss; 
%run_FE(data=hh_exp, dv=quant, iv= first_income*SHR_Grocery famsize year, 
		classvar = first_income famsize year, FEvar = month, outname = tmp);
data tmp; length test $30 Dependent $50 Parameter $50 ; set tmp; test = 'Cost-Grocery'; run; 
%my_merge(base_data=myreg, add_data=tmp, byvar=test Dependent, outname = myreg);

%run_FE(data=hh_exp, dv=num_mod, iv= first_income*SHR_Grocery famsize year, 
		classvar = first_income famsize year, FEvar = month, outname = tmp);
data tmp; length test $30 Dependent $50 Parameter $50 ; set tmp; test = 'Cost-Grocery'; run; 
%my_merge(base_data=myreg, add_data=tmp, byvar=test Dependent, outname = myreg);

*----------------------------------------------------------------*; 
* The cost of moving expenditure up -- expenditure to income ratio; 
%run_FE(data=hh_exp, dv=dol_rt, iv= first_income*SHR_Warehouse_Club famsize year, 
		classvar = first_income famsize year, FEvar = month, outname = tmp);
data tmp; length test $30 Dependent $50 Parameter $50 ; set tmp; test = 'Cost-Warehouse_Club'; run; 
%my_merge(base_data=myreg, add_data=tmp, byvar=test Dependent, outname = myreg);

%run_FE(data=hh_exp, dv=dol_purchases, iv= first_income*SHR_Warehouse_Club famsize year, 
		classvar = first_income famsize year, FEvar = month, outname = tmp);
data tmp; length test $30 Dependent $50 Parameter $50 ; set tmp; test = 'Cost-Warehouse_Club'; run; 
%my_merge(base_data=myreg, add_data=tmp, byvar=test Dependent, outname = myreg);

* The cost of moving expenditure down -- quantity and variety loss; 
%run_FE(data=hh_exp, dv=quant, iv= first_income*SHR_Warehouse_Club famsize year, 
		classvar = first_income famsize year, FEvar = month, outname = tmp);
data tmp; length test $30 Dependent $50 Parameter $50 ; set tmp; test = 'Cost-Warehouse_Club'; run; 
%my_merge(base_data=myreg, add_data=tmp, byvar=test Dependent, outname = myreg);

%run_FE(data=hh_exp, dv=num_mod, iv= first_income*SHR_Warehouse_Club famsize year, 
		classvar = first_income famsize year, FEvar = month, outname = tmp);
data tmp; length test $30 Dependent $50 Parameter $50 ; set tmp; test = 'Cost-Warehouse_Club'; run; 
%my_merge(base_data=myreg, add_data=tmp, byvar=test Dependent, outname = myreg);

*----------------------------------------------------------------*; 
* The cost of change size tier;
%run_FE(data=purchases, dv=unit_price, iv= first_income*PCTL_SIZE famsize year month, 
		classvar = first_income famsize year month, FEvar = product_module_descr, outname = tmp);
data tmp; length test $30 Dependent $50 Parameter $50 ; set tmp; test = 'Cost-SIZE'; run; 
%my_merge(base_data=myreg, add_data=tmp, byvar=test Dependent, outname = myreg);

%run_FE(data=hh_exp, dv=dol_rt, iv= first_income*PCTL_SIZE famsize year, 
		classvar = first_income famsize year, FEvar = month, outname = tmp);
data tmp; length test $30 Dependent $50 Parameter $50 ; set tmp; test = 'Cost-SIZE'; run; 
%my_merge(base_data=myreg, add_data=tmp, byvar=test Dependent, outname = myreg);

%run_FE(data=hh_exp, dv=dol_purchases, iv= first_income*PCTL_SIZE famsize year, 
		classvar = first_income famsize year, FEvar = month, outname = tmp);
data tmp; length test $30 Dependent $50 Parameter $50 ; set tmp; test = 'Cost-SIZE'; run; 
%my_merge(base_data=myreg, add_data=tmp, byvar=test Dependent, outname = myreg);

%run_FE(data=hh_exp, dv=num_mod, iv= first_income*PCTL_SIZE famsize year, 
		classvar = first_income famsize year, FEvar = month, outname = tmp);
data tmp; length test $30 Dependent $50 Parameter $50 ; set tmp; test = 'Cost-SIZE'; run; 
%my_merge(base_data=myreg, add_data=tmp, byvar=test Dependent, outname = myreg);

%run_FE(data=hh_exp, dv=quant, iv= first_income*PCTL_SIZE famsize year, 
		classvar = first_income famsize year, FEvar = month, outname = tmp);
data tmp; length test $30 Dependent $50 Parameter $50 ; set tmp; test = 'Cost-SIZE'; run; 
%my_merge(base_data=myreg, add_data=tmp, byvar=test Dependent, outname = myreg);

* Export regression results to excel; 
PROC EXPORT DATA= myreg
            OUTFILE= "&reg_outfile" 
            DBMS=csv REPLACE;
     PUTNAMES=YES;
RUN;

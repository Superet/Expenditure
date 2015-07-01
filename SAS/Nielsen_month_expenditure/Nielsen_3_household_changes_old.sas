**************************************************************;
* Check the time variation of demographics within households *;
**************************************************************;
* Merge in the codebook to the panelist data; 
proc sql noprint;
	create table my_data as 
	select A.mkt_type, A.scantrack_market_descr, A.household_code, A.panel_year, A.magnet_flag, A.magnet_const_flag, A.household_income, B.description as INCOME_descr, 
			B.mid_range as INCOME_midvalue
	from mainex.panelists as A left join (select code_value, description, mid_range from mysel.my_codebook where variable='household_income') as B
	on A.household_income = B.code_value
	order by household_code, panel_year descending;
quit;

* Note that income reported was the income 2 years prior to the panel_year;
data my_data;
	set my_data;
	income_year = panel_year - 2;
run;

*------------------------------------------------------------------------------*;
* By market Regression to test any changes in income within households over time;
proc sort data=my_data; by household_code descending income_year; run;
proc glm data=my_data(where = (mkt_type="&mkt1")) ;
	absorb household_code;
	class income_year(ref=first);
	model income_midvalue = income_year / solution;
	ods output ParameterEstimates=tmp1;
run;

proc glm data=my_data(where = (mkt_type="&mkt2"));
	absorb household_code;
	class income_year (ref=first);
	model income_midvalue = income_year / solution;
	ods output ParameterEstimates=tmp2;
run;

data tmp1; set tmp1(rename= (Estimate=&mkt1)); proc sort; by Parameter; run;
data tmp2; set tmp2(rename= (Estimate=&mkt2)); proc sort; by Parameter; run;

data tmp;
	merge tmp1(keep=Parameter &mkt1) tmp2(keep=Parameter &mkt2); 
	by Parameter;
	year = input(substr(Parameter, 13),4.);
run;

* Plot the year trend of income;
proc sgplot data=tmp;
	series x=year y=Goodcity;
	series x=year y=Badcity;
	title 'The year trend of average income by market';
run;

*------------------------------------------------------------------------------*;
* By magnet flag Regression to test any changes in income within households over time;
proc sort data=my_data; by household_code descending income_year; run;
proc glm data=my_data(where = (magnet_const_flag =0));
	absorb household_code;
	class income_year (ref=first);
	model income_midvalue = income_year / solution;
	ods output ParameterEstimates=tmp1;
run;

proc glm data=my_data(where = (magnet_const_flag =1));
	absorb household_code;
	class income_year (ref=first);
	model income_midvalue = income_year / solution;
	ods output ParameterEstimates=tmp2;
run;

data tmp1; set tmp1(rename= (Estimate=Other)); proc sort; by Parameter; run;
data tmp2; set tmp2(rename= (Estimate=Magnet)); proc sort; by Parameter; run;

data tmp;
	merge tmp1(keep=Parameter Other) tmp2(keep=Parameter Magnet); 
	by Parameter;
	year = input(substr(Parameter, 13),4.);
run;

* Plot the year trend of income;
proc sgplot data=tmp;
	series x=year y=Other;
	series x=year y=Magnet;
	title 'The year trend of average income by magnet households';
run;

**************************************;
* Trend of monthly total expenditure *;
**************************************;
* Group household based on first-year income;
proc sort data=mainex.panelists; by household_code panel_year; run;
data panelists;
	set mainex.panelists;
	by household_code panel_year;
	first_year = 0;
	income_group = 'Low';
	if first.panel_year then do 
		first_year=1;
		if household_income>=17 and household_income<=21 then income_group = 'Med';
		if household_income>21 then income_group = 'Hig';
	end;
	income_year = panel_year - 2;
run;

* Aggregate daily trip data to monthly expenditure data; 	
proc sql noprint;
	create table tmp as 
	select  household_code, year, month, ymonth,channel_group, mean(ymonth_num) as ymonth_num, sum(total_spent) as DOL
	from mainex.trips 
	group by household_code, year, month, ymonth, channel_group;
	
	create table hh_month_exp as 
	select A.*, B.mkt_type, B.scantrack_market_descr, B.magnet_const_flag, B.magnet_flag, B.household_income, B.income_group
	from tmp as A inner join panelists as B
	on A.household_code = B.household_code and A.year = B.panel_year
	order by scantrack_market_descr, household_code, year, month;
quit;

*------------------------------------------------------------------------;
* Run a regression of monthly expenditure against time with all panelists;
proc sort data=hh_month_exp; by household_code month descending ymonth; run;
proc glm data=hh_month_exp(where=(channel_group="Grocery")) ;
	absorb household_code month;
	class ymonth_num; 
	model DOL = ymonth_num/ solution;
	ods output ParameterEstimates=tmp;
run;

data tmp;
	set tmp;
	start = input(substr(compress(Parameter, '','A'),2), 8.);
proc sort; by start;
run;

data tmp1; 
	merge tmp ymonth_num_fmt(keep=ymonth_date start);
	by start;
proc sort; by ymonth_date;
run;
data tmp1; set tmp1; where Estimate^=0; run;

* Plot the trend;
proc sgplot data=tmp1;
	series x = ymonth_date y = Estimate;
	title 'The trend of grocery expenditure with monthly data';
run;

*-------------------------------------;
* Trend of total expenditure by market; 
proc sort data = hh_month_exp; by household_code month; run;
proc glm data=hh_month_exp(where=(channel_group="Grocery" & mkt_type="&mkt1"));
	absorb household_code month;
	class ymonth_num; 
	model DOL = ymonth_num / solution;
	ods output ParameterEstimates=tmp1;
run;

data tmp1;
	set tmp1(rename=(Estimate=&mkt1));
	start = input(substr(compress(Parameter, '','A'),2), 8.);
proc sort; by start;
run;

data tmp1; 
	merge tmp1 ymonth_num_fmt(keep=ymonth_date start);
	by start;
proc sort; by ymonth_date;
run;
data tmp1; set tmp1; where &mkt1^=0; run;

* Regression for the second market;
proc glm data=hh_month_exp(where=(channel_group="Grocery" & mkt_type="&mkt2")) ;
	absorb household_code month;
	class ymonth_num; 
	model DOL = ymonth_num / solution;
	ods output ParameterEstimates=tmp2;
run;

data tmp2;
	set tmp2(rename=(Estimate=&mkt2));
	start = input(substr(compress(Parameter, '','A'),2), 8.);
proc sort; by start;
run;

data tmp2; 
	merge tmp2 ymonth_num_fmt(keep=ymonth_date start);
	by start;
proc sort; by ymonth_date;
run;
data tmp2; set tmp2; where &mkt2^=0; run;

* Put together the estimates from two regressions; 
data tmp; 
	merge tmp1(keep=ymonth_date &mkt1) tmp2(keep=ymonth_date &mkt2);
	by ymonth_date;
run;

* Plot the trend;
proc sgplot data=tmp;
	series x = ymonth_date y = Goodcity;
	series x = ymonth_date y = Badcity;
	title 'The trend of monthly grocery expenditure by market';
run;

*-----------------------------------------------;
* Trend of total exnepdutre of magnet households;
* Sum up grocery expenditure and restaurant expenditure for magnet households;
proc sql noprint;
	create table tmp_data as
	select mkt_type, scantrack_market_descr,household_code, year, month, ymonth, mean(ymonth_num) as ymonth_num, sum(DOL) as DOL
	from hh_month_exp
	where magnet_const_flag = 1
	group by scantrack_market_descr, household_code, year, month, ymonth
	order by scantrack_market_descr, household_code, year, month, ymonth;
quit;

* Regression of grocery expenditure;
proc glm data=hh_month_exp(where=(channel_group="Grocery" & magnet_const_flag = 1));
	absorb household_code month;
	class ymonth_num ; 
	model DOL = ymonth_num / solution;
	ods output ParameterEstimates=tmp1(rename=(Estimate=Grocery));
run;

* Regression of the sum of grocern and restaurant expenditure;
proc sort data=tmp_data; by household_code month; run;
proc glm data=tmp_data;
	absorb household_code month;
	class ymonth_num; 
	model DOL = ymonth_num / solution;
	ods output ParameterEstimates=tmp2(rename=(Estimate=Restaurant_grocery));
run;

* Merge the estiamtes from the two regressions;
proc sort data=tmp1; by parameter; 
proc sort data=tmp2; by parameter; run;
data tmp; 
	merge tmp1(keep=Parameter Grocery) tmp2(keep=Parameter Restaurant_grocery);
	by Parameter;
	start = input(substr(compress(Parameter, '','A'),2), 8.);
proc sort; by start;
run;

data tmp; 
	merge tmp ymonth_num_fmt(keep=ymonth_date start);
	by start;
proc sort; by ymonth_date;
run;
data tmp; set tmp; where Grocery^=0; run;

* Plot the trend;
proc sgplot data=tmp;
	series x = ymonth_date y = Grocery;
	series x = ymonth_date y = Restaurant_grocery;
	title 'The trend of monthly grocery expenditure for magnet households';
run;

******************************************************************;
* Change of monthly expenditure share in channel types over time *;
******************************************************************;
* Aggregate the daily purchase to monthly expenditure; 
proc sql noprint;
	create table tmp as
	select  household_code, year, ymonth, month, channel_type, mean(ymonth_num) as ymonth_num, sum(total_spent) as DOL
	from mainex.trips
	group by household_code, year, ymonth, month, channel_type;

	create table tmp1 as 
	select *, DOL/sum(DOL) as share
	from tmp 
	group by household_code,year, ymonth, month
	order by household_code,year, ymonth, month;
	
	create table month_share as 
	select A.*, B.mkt_type, B.scantrack_market_descr, B.magnet_flag, B.magnet_const_flag, B.income_group
	from tmp1 as A left join panelists as B
	on A.household_code = B.household_code and A.year = B.panel_year;
quit;	

*---------------------------------------------------------------------*;
* Regression to see the trend of expenditure share from the entire data;
proc sort data=month_share; by household_code month descending ymonth; run;
proc glm data=month_share(where = (channel_type = 'Discount Store'));
	absorb household_code month;
	class  ymonth_num;
	model share =  ymonth_num /solution;
	ods output ParameterEstimates=tmp;
run;

data tmp1; 
	set tmp;
	start = input(substr(compress(Parameter, '','A'),2), 8.);
proc sort; by start;
run;

data tmp1; 
	merge tmp1 ymonth_num_fmt(keep=ymonth_date start);
	by start;
proc sort; by ymonth_date;
run;
data tmp1; set tmp1; where Estimate^=0; run;
run;

* Plot the trend;
proc sgplot data=tmp1;
	series x = ymonth_date y = Estimate;
	title 'The trend of expenditure share at discount stores with monthly data';
run;

*--------------------------------------------------------------------;
* Trend of expenditure share of discount stores for magnet households;
proc sort data=month_share; by household_code month descending ymonth; run;
proc glm data=month_share(where = (channel_type = 'Discount Store' & magnet_const_flag=1));
	absorb household_code month;
	class  ymonth_num;
	model share =  ymonth_num /solution;
	ods output ParameterEstimates=tmp1(rename=(Estimate=Magnet));
run;

proc glm data=month_share(where = (channel_type = 'Discount Store' & magnet_const_flag=0 & magnet_flag=0) );
	absorb household_code month;
	class  ymonth_num;
	model share =  ymonth_num /solution;
	ods output ParameterEstimates=tmp2(rename=(Estimate=NonMagnet));
run;

* Put together the estimates from two regressions;
proc sort data=tmp1; by parameter;
proc sort data=tmp2; by parameter; 
run;
data tmp;
	merge tmp1(keep=Parameter Magnet) tmp2(keep=Parameter NonMagnet);
	by Parameter;
	start = input(substr(compress(Parameter, '','A'),2), 8.);
proc sort; by start;
run;

data tmp; 
	merge tmp ymonth_num_fmt(keep=ymonth_date start);
	by start;
proc sort; by ymonth_date;
run;
data tmp; set tmp; where NonMagnet^=0; run;
run;

* Plot the trend;
proc sgplot data=tmp;
	series x = ymonth_date y = Magnet;
	series x = ymonth_date y = NonMagnet;
	title 'The trend of expenditure share at discount stores by Magnet households';
run;

*------------------------------------------------------------------------------------------;
* Trend of expenditure share of discount stores by market (Only within grcoery expenditure);
* Ignore magnet households;
proc sort data=month_share; by household_code month descending ymonth; run;
proc glm data=month_share(where = (channel_type = 'Discount Store' & magnet_const_flag=0 & mkt_type="&mkt1"));
	absorb household_code month;
	class  ymonth_num;
	model share =  ymonth_num /solution;
	ods output ParameterEstimates=tmp1;
run;

proc glm data=month_share(where = (channel_type = 'Discount Store' & magnet_const_flag=0 & mkt_type="&mkt2") );
	absorb household_code month;
	class  ymonth_num;
	model share =  ymonth_num /solution;
	ods output ParameterEstimates=tmp2;
run;

proc sort data=tmp1(rename=(Estimate=&mkt1)); by parameter;
proc sort data=tmp2(rename=(Estimate=&mkt2)); by parameter;
run;
data tmp;
	merge tmp1(keep=Parameter &mkt1) tmp2(keep=Parameter &mkt2);
	by Parameter;
	start = input(substr(compress(Parameter, '','A'),2), 8.);
proc sort; by start;
run;

data tmp; 
	merge tmp ymonth_num_fmt(keep=ymonth_date start);
	by start;
proc sort; by ymonth_date;
run;
data tmp; set tmp; where &mkt1^=0; run;
run;

* Plot the trend;
proc sgplot data=tmp;
	series x = ymonth_date y = Goodcity;
	series x = ymonth_date y = Badcity;
	title 'The trend of expenditure share at discount stores by market among NonMagnet households';
run;

*******************************;
* Run expenditure regressions *;
*******************************;
data mydata_rec0; set mainex.hh_month_exp; where recession = 0; 
data mydata_rec1; set mainex.hh_month_exp; where recession = 1; 
run;

*----------------------------------------------------------------------------------*;
* A function that run GLM regressions by settting the data with a grouping variable ;
%macro exp_reg(dataset, dv, iv, classvar, groupvar, resultname);
	%local i numvar tmplast resultname1 tmpalast;
	proc sql noprint; 
		select count(distinct &groupvar) into: numvar from &dataset;
		create table varlevel as select distinct &groupvar as level from &dataset;
	quit;
	data varlevel; set varlevel; n = _N_; run;
	
	%do i=1 %to &numvar;
		* Subset the data;
		proc sql noprint; 
			create table mydata_run as
			select *
			from &dataset
			where &groupvar in (select level from varlevel where n= &i);
			
			select level into: varvalue from varlevel where n= &i;
		quit;
		
		* Run regression with the subset;
		proc glm data=mydata_run;
			absorb household_code;
			class &classvar;
			model &dv = &iv /solution; 
			ods output ParameterEstimates=tmp&i OverallANOVA=tmpa&i;
			title "GLM results for &groupvar is &varvalue with data &dataset";
		run;
		
		data tmp&i; set tmp&i; &groupvar = "&varvalue";
		data tmpa&i; set tmpa&i; &groupvar = "&varvalue";
	%end;
	
	* Put together the results;
	%let tmplast = tmp%trim(&numvar);
	data &resultname;
		set tmp1 - &tmplast;
		dataset = "&dataset";
	run;
	
	%let resultname1 = %trim(&resultname)_anova;
	%let tmpalast = tmpa%trim(&numvar);
	data &resultname1;
		set tmpa1 - &tmpalast;
		dataset = "&dataset";
	run;
	
	* Delete working data tables;
	proc datasets noprint; delete tmp1-&tmplast tmpa1 - &tmpalast varlevel mydata_run; run;
%mend exp_reg;

*----------------------------------------------------------*;
* A function that run the macro exp_reg with different IVs *; 
%macro exp_reg_outer(vardata, result_name);
	%local ivlist classlist varlab;
	proc sql noprint;
		select count(distinct n) into: numvar from &vardata;
	quit;
	
	%do i=1 %to &numvar;
		* Select the independent variables; 
		proc sql noprint;
			select iv into :ivlist separated by ' ' from &vardata where n=&i;
			select iv into :classlist separated by ' ' from &vardata where n=&i and classvar=1;
			select distinct income_use into :varlab from &vardata where n=&i;
		quit;
		%put &ivlist;
		%put &classlist;
		%put &varlab;
		
		* Run regressions by groupvar;
		%exp_reg(dataset=mydata_rec0, dv = ln_dol, iv = &ivlist, classvar = &classlist, groupvar = mkt_type, resultname = result_rec0);
		%exp_reg(dataset=mydata_rec1, dv = ln_dol, iv = &ivlist, classvar = &classlist, groupvar = mkt_type, resultname = result_rec1);
		
		* Combine results together;
		data param&i;
			set result_rec0 result_rec1;
			Income_value = "&varlab";
		run;
		data anova&i;
			set result_rec0_anova result_rec1_anova;
			Income_value = "&varlab";
		run;
		
		* Delete working tables;
		proc datasets noprint; delete result_rec0 result_rec1 result_rec0_anova result_rec1_anova; run;
	%end;
	
	%let tmplast = param%trim(&numvar);
	data &result_name;
		set param1 - &tmplast;
	run;
	
	%let result1 = %trim(&result_name)_anova;
	%let tmp1last = anova%trim(&numvar);
	data &result1;
		set anova1 - &tmp1last;
	run;
	
	proc datasets noprint; delete param1 - &tmplast anova1 - &tmp1last; run;
%mend exp_reg_outer;

*-------------------------------------------------------*;
* Try regression with different values of income variable; 
data myvar; 
	input iv $20. classvar 2.;
	datalines;
	ln_lag_dol		0
	stone_price		0
	income_real		1
	famsize			1
	condo			1
	employed		1
	NumChild		1
	month			1
	;
run;

data myvar;
	set myvar;
	do n = 1 to 3;
		output;
	end;
proc sort; by n;
run;

data myvar;
	length income_use $ 20;
	set myvar;
	income_use = "income_real";
	if n=2 then	income_use = "income_group";
	else if n=3 then income_use = "income_mid";
	if n=2 and iv="income_real" then iv = "income_group";
	if n=3 and iv="income_real" then do;
		iv = "ln_income_midvalue";
		classvar = 0;
	end;
run;

%exp_reg_outer(vardata=myvar, result_name=resultExp);

libname mylib "E:\Users\ccv103\Documents\Research\Store switching\SAS_temp\Nielsen";
libname mysel "E:\Users\ccv103\Documents\Research\Store switching\SAS_temp\Nielsen\Work";
libname mainex "E:\Users\ccv103\Documents\Research\Store switching\SAS_temp\Nielsen\Work\main";

%let dirname= E:\Users\ccv103\Documents\Research\Nielsen\nielsen_extracts\HMS;

************************************************;
* Fixed effect regression comparision function *;
************************************************;
%macro trend_reg(data, dv, groupvar, base, outname);
	proc sql noprint; select unique(&groupvar) into :comp from &data where &groupvar ^= "&base"; quit; 
	
	proc glm data=&data(where = (&groupvar="&base"));
		absorb household_code; 
		class month ymonth_num; 
		model &dv = month ymonth_num / solution noint;
		ods output ParameterEstimates=tmp1;
		title "Fixed effect regression for &dv in base";
	run;
	
	proc glm data=&data(where = (&groupvar^="&base"));
		absorb household_code; 
		class month ymonth_num; 
		model &dv = month ymonth_num / solution noint;
		ods output ParameterEstimates=tmp2;
		title "Fixed effect regression for &dv in comparison";
	run;
	
	data tmp1; length Dependent $30. &groupvar $15.; set tmp1; &groupvar = "&base";
	data tmp2; length Dependent $30. &groupvar $15.; set tmp2; &groupvar = "&comp";
	data &outname; set tmp1 tmp2; run;
	
	proc datasets noprint; delete tmp1 tmp2; run;
%mend trend_reg; 

***************;
* Create data *;
***************;
* Read in selected markets;
PROC IMPORT OUT= cpi
            DATAFILE= "\\tsclient\Resear1\Store switching\processed data\CPI.csv" 
            DBMS=csv REPLACE;
     		GETNAMES=YES; 
RUN;

*----------------------------------------*;
* Format year-month for later regressions;
proc sql noprint;
	create table ymonth as
	select distinct ymonth
	from mainex.trips
	order by ymonth;
quit;

* Define year-month format; 
data ymonth_fmt;
	set ymonth(rename = (ymonth=start));
	retain fmtname '$ymonth' type 'C';
	ymonth_date = input(strip('01'||start), date9.);
	format ymonth_date date9.;
proc sort; by ymonth_date;
run;

data ymonth_fmt;
	set ymonth_fmt;
	label = start; 
run;

data ymonth_num_fmt;
	set ymonth_fmt(drop=start fmtname type);
	year = year(ymonth_date);
	month = month(ymonth_date);
proc sort; by descending year descending month; 
run;

data ymonth_num_fmt;
	set ymonth_num_fmt;
	retain fmtname 'ymonth_num' type 'N';
	start = _N_;
run;

proc format library = WORK cntlin = ymonth_num_fmt; run;
proc datasets noprint; delete ymonth ymonth_fmt; run;

*--------------------------*;
* Create the regression data;
data tmp; 
	set mylib.all_hh_month_exp_merge;
	where market in ("Buffalo - Rochester","Denver","Sacramento", "Detroit");
	Q = storable_quant + nonstr_quant; 
	n_mod = num_food_module + num_noned_module;
	unit_price = dol_purchases/Q; 
	exp_income = dol_purchases/(income_midvalue/12); 
	array dolvar{*} dolp_:;
	do i=1 to dim(dolvar);
		if dol_purchases ^=. then dolvar[i] = dolvar[i]/dol_purchases; 
	end;
	drop i;
run;

proc sql noprint;
	create table tmp1 as 
	select A.*, B.start as ymonth_num, B.ymonth_date
	from tmp as A left join ymonth_num_fmt as B
	on A.year = B.year and A.month=B.month
	order by household_code; 
	
	select Annual into: price_base from cpi where year = 2004; 
	
	create table mydata as 
	select A.*, B.Annual/&price_base as cpi, A.unit_price/(B.Annual/&price_base) as unit_price2004
	from tmp1 as A left join cpi as B
	on A.year = B.year
	order by household_code; 
quit;
proc contents data=mydata; run;
proc datasets noprint; delete tmp tmp1; run;


*******************;
* Run regressions *;
*******************;
* Monthly expenditure ;
%trend_reg(mydata, DOL, mkt_type, Badcity, myfit);
%trend_reg(mydata, dol_purchases, mkt_type, Badcity, tmp);
data myfit; set myfit tmp; run;

%trend_reg(mydata, exp_income, mkt_type, Badcity, tmp);
data myfit; set myfit tmp; run;

* Expenditure share; 
%trend_reg(mydata, dolp_Convenience_Store, mkt_type, Badcity, tmp);
data myfit; set myfit tmp; run;

%trend_reg(mydata, DOLP_Discount_Store, mkt_type, Badcity, tmp);
data myfit; set myfit tmp; run;

%trend_reg(mydata, DOLP_Dollar_Store, mkt_type, Badcity, tmp);
data myfit; set myfit tmp; run;

%trend_reg(mydata, DOLP_Drug_Store, mkt_type, Badcity, tmp);
data myfit; set myfit tmp; run;

%trend_reg(mydata, DOLP_Grocery, mkt_type, Badcity, tmp);
data myfit; set myfit tmp; run;

%trend_reg(mydata, DOLP_Health_Food_Store, mkt_type, Badcity, tmp);
data myfit; set myfit tmp; run;

%trend_reg(mydata, DOLP_Warehouse_Club, mkt_type, Badcity, tmp);
data myfit; set myfit tmp; run;

* Basket quantity; 
%trend_reg(mydata, Q, mkt_type, Badcity, tmp);
data myfit; set myfit tmp; run;

%trend_reg(mydata, storable_quant, mkt_type, Badcity, tmp);
data myfit; set myfit tmp; run;

%trend_reg(mydata, nonstr_quant, mkt_type, Badcity, tmp);
data myfit; set myfit tmp; run;

%trend_reg(mydata, n_mod, mkt_type, Badcity, tmp);
data myfit; set myfit tmp; run;

* Trip frequency; 
%trend_reg(mydata, num_day, mkt_type, Badcity, tmp);
data myfit; set myfit tmp; run;

%trend_reg(mydata, num_trip, mkt_type, Badcity, tmp);
data myfit; set myfit tmp; run;

* Coupon usage; 
%trend_reg(mydata, pct_coupon, mkt_type, Badcity, tmp);
data myfit; set myfit tmp; run;

%trend_reg(mydata, unit_price2004, mkt_type, Badcity, tmp);
data myfit; set myfit tmp; run;

* Overall sample; 
proc glm data=mydata;
	absorb household_code; 
	class month ymonth_num; 
	model  DOLP_Discount_Store= month ymonth_num / solution noint;
	ods output ParameterEstimates=tmp;
	title "Fixed effect regression for DOLP_Discount_Store";
run;
data tmp; length Dependent $30. mkt_type $15.; set tmp; mkt_type = ""; run;

data myfit; set myfit tmp; run;

*----------------------*;
* Panlist income trend *;
data tmp; 
	set mylib.all_panelists; 
	where market in ("Buffalo - Rochester","Denver","Sacramento", "Detroit");
proc sort; by household_code  descending panel_year; 
run;

proc sql noprint;
	create table tmp_year as
	select unique(panel_year) as panel_year
	from mydata
	order by panel_year descending;
quit;

data tmp_year; set tmp_year; year_num = _N_; run;

proc sql noprint;
	create table mydata as 
	select A.*, B.year_num
	from tmp as A left join tmp_year as B
	on A.panel_year = B.panel_year
	order by household_code; 
quit;

* Run separately; 
proc glm data=mydata(where=(mkt_type ="Badcity"));
	absorb household_code; 
	class year_num; 
	model income_midvalue = year_num / solution noint;
	ods output ParameterEstimates=tmp;
	title "Fixed effect regression for income";
run;
data tmp; length Dependent $30. mkt_type $15.; set tmp; mkt_type = "Badcity"; run;
data myfit; set myfit tmp; run;

proc glm data=mydata(where=(mkt_type ="Goodcity")) ;
	absorb household_code; 
	class year_num; 
	model income_midvalue = year_num / solution noint;
	ods output ParameterEstimates=tmp;
	title "Fixed effect regression for income";
run;
data tmp; length Dependent $30. mkt_type $15.; set tmp; mkt_type = "Goodcity"; run;
data myfit; set myfit tmp; run;

*-----------------------------------------------*;
* Organize regression results and export the data; 
data tmp1; 
	set myfit; 
	where substr(Parameter,1,6)="ymonth";
	start = input(substr(compress(Parameter, '','A'),2), 8.);
proc sort; by start;
run;

data tmp1; 
	merge tmp1 ymonth_num_fmt(keep=ymonth_date start);
	by start;
	proc sort; by Dependent mkt_type ymonth_date;
run;

data tmp2;
	set myfit;
	where substr(Parameter,1,4)="year";
	year_num = input(substr(compress(Parameter, '','A'),2), 8.);
proc sort; by year_num;
run;

data tmp2;
	merge tmp2 tmp_year;
	by year_num;
proc sort; by Dependent mkt_type panel_year; 
run;

data myfit_export; 
	merge tmp1 tmp2(drop=year_num);
	by Dependent;
run;

PROC EXPORT DATA=myfit_export
			OUTFILE= "\\tsclient\Resear1\Store switching\processed data\trend_fit.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;

************************************************************************;
* Segment households according to the direction of their income change *;
************************************************************************;
*--------------------------*;
* Create the regression data;
data tmp; 
	set mylib.all_hh_month_exp_merge;
	where market in ("Buffalo - Rochester","Denver","Sacramento", "Detroit");
	Q = storable_quant + nonstr_quant; 
	n_mod = num_food_module + num_noned_module;
	unit_price = dol_purchases/Q; 
	array dolvar{*} dolp_:;
	do i=1 to dim(dolvar);
		if dol_purchases ^=. then dolvar[i] = dolvar[i]/dol_purchases; 
	end;
	drop i;
run;

proc sql noprint;
	create table tmp1 as 
	select A.*, B.start as ymonth_num, B.ymonth_date
	from tmp as A left join ymonth_num_fmt as B
	on A.year = B.year and A.month=B.month
	order by household_code; 
	
	select Annual into: price_base from cpi where year = 2004; 
	
	create table mydata as 
	select A.*, B.Annual/&price_base as cpi, A.unit_price/(B.Annual/&price_base) as unit_price2004
	from tmp1 as A left join cpi as B
	on A.year = B.year
	order by household_code; 
quit;
proc datasets noprint; delete tmp tmp1; run;

* Flag household income change; 
data tmp;
	set panelists;
	where market in ("Buffalo - Rochester","Denver","Sacramento", "Detroit");
proc sort; by household_code panel_year;
run;

data tmp;
	set tmp;
	income_dif = dif(income_midvalue);
run;

proc sql noprint;
	create table tmp1 as 
	select household_code, sum(income_dif) as income_dif, 1*(sum(income_dif)<0) as income_down
	from tmp
	group by household_code;
quit;

proc freq data=tmp1;
table income_down;
run;

proc sql noprint;
	create table tmp as 
	select A.*, B.income_down, ifc(B.income_down=1, "Decreasing","Non-decreasing") as income_direction
	from mydata as A left join tmp1 as B
	on A.household_code = B.household_code;
quit;
data mydata; set tmp; run;

*******************;
* Run regressions *;
*******************;
* Monthly expenditure ;
%trend_reg(mydata, DOL, income_direction, Decreasing, myfit);
%trend_reg(mydata, dol_purchases, income_direction, Decreasing, tmp);
data myfit; set myfit tmp; run;

%trend_reg(mydata, exp_income, mkt_type, Badcity, tmp);
data myfit; set myfit tmp; run;

* Expenditure share; 
%trend_reg(mydata, dolp_Convenience_Store, income_direction, Decreasing, tmp);
data myfit; set myfit tmp; run;

%trend_reg(mydata, DOLP_Discount_Store, income_direction, Decreasing, tmp);
data myfit; set myfit tmp; run;

%trend_reg(mydata, DOLP_Dollar_Store, income_direction, Decreasing, tmp);
data myfit; set myfit tmp; run;

%trend_reg(mydata, DOLP_Drug_Store, income_direction, Decreasing, tmp);
data myfit; set myfit tmp; run;

%trend_reg(mydata, DOLP_Grocery, income_direction, Decreasing, tmp);
data myfit; set myfit tmp; run;

%trend_reg(mydata, DOLP_Health_Food_Store, income_direction, Decreasing, tmp);
data myfit; set myfit tmp; run;

%trend_reg(mydata, DOLP_Warehouse_Club, income_direction, Decreasing, tmp);
data myfit; set myfit tmp; run;

* Basket quantity; 
%trend_reg(mydata, Q, income_direction, Decreasing, tmp);
data myfit; set myfit tmp; run;

%trend_reg(mydata, storable_quant, income_direction, Decreasing, tmp);
data myfit; set myfit tmp; run;

%trend_reg(mydata, nonstr_quant, income_direction, Decreasing, tmp);
data myfit; set myfit tmp; run;

%trend_reg(mydata, n_mod, income_direction, Decreasing, tmp);
data myfit; set myfit tmp; run;

* Trip frequency; 
%trend_reg(mydata, num_day, income_direction, Decreasing, tmp);
data myfit; set myfit tmp; run;

%trend_reg(mydata, num_trip, income_direction, Decreasing, tmp);
data myfit; set myfit tmp; run;

* Coupon usage; 
%trend_reg(mydata, pct_coupon, income_direction, Decreasing, tmp);
data myfit; set myfit tmp; run;

%trend_reg(mydata, unit_price2004, income_direction, Decreasing, tmp);
data myfit; set myfit tmp; run;

*-----------------------------------------------*;
* Organize regression results and export the data; 
data myfit_export; 
	set myfit; 
	where substr(Parameter,1,5)^="month";
	start = input(substr(compress(Parameter, '','A'),2), 8.);
proc sort; by start;
run;

data myfit_export; 
	merge myfit_export ymonth_num_fmt(keep=ymonth_date start);
	by start;
	proc sort; by Dependent income_direction ymonth_date;
run;

PROC EXPORT DATA=myfit_export
			OUTFILE= "\\tsclient\Resear1\Store switching\processed data\trend_fit.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;






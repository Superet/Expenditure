libname mylib "E:\Users\ccv103\Documents\Research\Store switching\SAS_temp\Nielsen";
libname mysel "E:\Users\ccv103\Documents\Research\Store switching\SAS_temp\Nielsen\Work";
libname mainex "E:\Users\ccv103\Documents\Research\Store switching\SAS_temp\Nielsen\Work\main";

proc datasets library=WORK kill noprint; run;

%let dirname= E:\Users\ccv103\Documents\Research\Nielsen\nielsen_extracts\HMS;
%let plotfile1 = E:\Users\ccv103\Desktop\graph_raw_trend.pdf;
%let plotfile2 = E:\Users\ccv103\Desktop\graph_within_trend.pdf;

*************;
* Functions *;
*************;
*-------------------------*;
* Trend plot from raw data ;
%macro raw_mean_trend(data, byv, classv, v);
	proc univariate data = &data noprint;
		by &byv; 
		class &classv; 
		var &v; 
		output out = tmp mean=mean STD=std pctlpts = 2.5 25 75 97.5 pctlpre = P;
	run;

	proc sgplot data=tmp;
		*band 	x = &classv lower=P25 upper=P75 / group = &byv transparency=0.7;
		series 	x = &classv y=mean / group = &byv transparency=0.7;
		title "&v over year month by &byv";
	run;
%mend raw_mean_trend;

*-------------------------*;
* Run fixed effect model; 
%macro trend_reg(data, dv, iv, classv, groupvar, fixed, outname);
	proc sql noprint; 
		select unique(&groupvar) into :grouplevel separated by ' ' from &data; 
		select count(unique(&groupvar)) into :n from &data; 
	quit; 
	
	%do i=1 %to &n;
		proc glm data=&data(where = (&groupvar="%scan(&grouplevel, &i, ' ')" ));
			absorb &fixed; 
			class &classv; 
			model &dv = &iv / solution noint;
			ods output ParameterEstimates=tmp&i;
			title "Fixed effect regression of &dv for &groupvar = %scan(&grouplevel, &i, ' ')";
		run;
		data tmp&i; length Dependent $30. &groupvar $15.; set tmp&i; &groupvar = "%scan(&grouplevel, &i, ' ')";
	%end;
	
	%let tmplast = tmp%trim(&n);
	data &outname; set tmp1-&tmplast; run;
	proc univariate data=&data noprint;
		class &groupvar;
		var &dv;
		output out = tmp_mean MEAN = mean P50=median;
	run;
	
	proc sort data=tmp_mean; by &groupvar; run;
	data &outname;
		merge &outname tmp_mean;
		by &groupvar;
	run;
	
	proc datasets noprint; delete tmp1-&tmplast tmp_mean; run;
%mend trend_reg; 

*-------------------------*;
* Plot the trend; 
%macro trendplot(data);
	data tmp;
		set &data;
		where substr(Parameter,1,6)="ymonth";
		start = input(substr(compress(Parameter, '','A'),2), 8.);
	proc sort; by start;
	run;

	data tmp; 
		merge tmp ymonth_num_fmt(keep=ymonth_date start);
		by start;
		y = Estimate + mean; 
		low = y - 1.96*StdErr;
		up = y + 1.96*StdErr;
		proc sort; by Dependent;
	run;
	
	proc sql noprint;
		select unique(Dependent) into : deplevel separated by ' ' from tmp;
		select count(unique(Dependent)) into :n from tmp;
	quit;
	
	%do i=1 %to &n;
		proc sgplot data=tmp(where=(Estimate^=0 & Estimate^=. & Dependent = "%scan(&deplevel, &i, ' ')"));
			band x = ymonth_date lower=low upper = up/group = seg_income transparency=0.5;
			series x = ymonth_date y = y/group = seg_income transparency=0.5; 
			title "Trend plot for %scan(&deplevel, &i, ' ')";
		run;
	%end;
%mend trenplot;
 
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

data year_fmt;
	do i=2004 to 2012;
	Year = i; 
	output;
	end;
	drop i;
proc sort; by descending year;run;
data year_fmt;
	set year_fmt;
	year_num = _N_;
run;

*------------------------------------------*;
* Segment households based on initial income; 
proc sort data=mylib.all_panelists; by household_code panel_year; run;

data tmp_panelists;
	set mylib.all_panelists; 
	where income_real ^= .;
	by household_code; 
	if ^first.household_code then delete; 
	if income_real <= 15 and income_real ^='' then seg_income = 'Qt1';
	if income_real>15 and income_real<=19 then seg_income = 'Qt2';
	if income_real>19 and income_real<=23 then seg_income = 'Qt3';
	if income_real>23 then seg_income = 'Qt4';
run;

proc sql noprint;
	create table tmp as
	select A.*, B.seg_income
	from mylib.all_panelists as A left join tmp_panelists as B
	on A.household_code = B.household_code;
	
	create table panelists as 
	select A.*, B.year_num
	from tmp as A left join year_fmt as B
	on A.panel_year = B.year; 
quit;

*--------------------------*;
* Create the regression data;
data tmp; 
	set mylib.all_hh_month_exp_merge;
	Q = storable_quant + nonstr_quant; 
	n_mod = num_food_module + num_noned_module;
	if Q^=. then unit_price = dol_purchases/Q; 
	if income_midvalue ^= . then exp_income = dol_purchases/(income_midvalue/12); 
	array dolvar{*} DOLP_Convenience_Store DOLP_Discount_Store DOLP_Dollar_Store DOLP_Drug_Store DOLP_Grocery DOLP_Health_Food_Store DOLP_Warehouse_Club;
	array sharevar{*} share_Convenience_Store share_Discount_Store share_Dollar_Store share_Drug_Store share_Grocery share_Health_Food_Store share_Warehouse_Club;
	do i=1 to dim(dolvar);
		if (dol_purchases ^=. & dol_purchases^ = 0 )then sharevar[i] = dolvar[i]/dol_purchases; 
		else sharevar[i] = 0; 
	end;
	drop i;
run;

* Reshape price data;
proc transpose data=mylib.Format_month_price out=tmp_price prefix = PRICE_; 
	by scantrack_market_descr year month;
	id channel_type;
	var price_paid_norm_index;
run;

proc sql noprint;
	* Add year-month number; 
	create table tmp1 as 
	select A.*, B.start as ymonth_num, B.ymonth_date
	from tmp as A left join ymonth_num_fmt as B
	on A.year = B.year and A.month=B.month
	order by household_code; 
	
	* Add CPI to normalize the prices; 
	select Annual into: price_base from cpi where year = 2004; 
	
	create table tmp2 as 
	select A.*, B.Annual/&price_base as cpi, A.unit_price/(B.Annual/&price_base) as unit_price2004
	from tmp1 as A left join cpi as B
	on A.year = B.year
	order by household_code; 
	
	* Add price data;
	create table tmp3(drop = old_mkt old_year old_month) as 
	select A.*, B.*
	from tmp2 as A left join tmp_price(drop=_NAME_ _LABEL_ rename=(scantrack_market_descr=old_mkt year=old_year month=old_month)) as B
	on A.scantrack_market_descr=B.old_mkt and A.year =B.old_year and A.month=B.old_month;
	
	* Add the segment variable; 
	create table mydata as
	select A.*, B.seg_income
	from tmp3 as A left join tmp_panelists as B
	on A.household_code = B.household_code
	order by household_code, year, month; 
quit;

* Compute quantity via expenditure and aggregate price; 
data mydata;
	set mydata;
	Q_comp = 0;
	array dolvar{*} dolp_:;
	array prcvar{*} price_:;
	do i=1 to dim(dolvar);
		Q_comp = Q_comp + dolvar[i]/prcvar[i];
	end;
	Q_dry = Q - magnet_quant;
run;

proc contents data=mydata; run;

proc datasets noprint; delete tmp tmp1 tmp2 tmp3; run;

********************;
* Raw data pattern *;
********************;
proc sort data=panelists nodupkeys out = tmp_year; by panel_year; run;

* --------------------------------------------------------*;
* Income shocks *;
ods graphics on;
ods select SmoothingComponentPlot;
proc gam data=mylib.all_panelists plots=components(clm);
	model income_midvalue= spline(panel_year);
run;
ods graphics off;
ods pdf close;

* Income trend by segment;
proc sort data=panelists; by seg_income; run;
ods graphics on;
ods select SmoothingComponentPlot;
proc gam data=panelists plots=components(clm commonaxes);
	by seg_income; 
	model income_midvalue= spline(panel_year);
	title "GAM smoothing of income over year";
run;

* --------------------------------------------------------*;
* Expenditure trend *;
proc sort data=mydata; by seg_income ymonth_date; run;

ods listing close;
ods pdf file = "&plotfile1";

%raw_mean_trend(mydata, seg_income, ymonth_date, dol_purchases);
	
* --------------------------------------------------------*;
* Outcome: trend of expenditure share ;
proc sql noprint;
	create table tmp as 
	select seg_income, ymonth_date, 
		sum(DOLP_Convenience_Store) as DOLP_Convenience_Store, sum(DOLP_Discount_Store) as DOLP_Discount_Store, 
		sum(DOLP_Dollar_Store) as DOLP_Dollar_Store, sum(DOLP_Drug_Store) as DOLP_Drug_Store, sum(DOLP_Grocery) as DOLP_Grocery, 
		sum(DOLP_Health_Food_Store) as DOLP_Health_Food_Store, sum(DOLP_Warehouse_Club) as DOL_Warehouse_Club
	from mydata
	group by seg_income, ymonth_date
	order by seg_income, ymonth_date;
quit;

proc transpose data=tmp out=tmp1(rename=(_NAME_=channel_type)) prefix=dol;
	by seg_income ymonth_date;
	var DOLP_Convenience_Store DOLP_Discount_Store DOLP_Dollar_Store DOLP_Drug_Store DOLP_Grocery DOLP_Health_Food_Store DOL_Warehouse_Club;
run;

proc sql noprint; 
	create table tmp2 as 
	select seg_income, ymonth_date, channel_type, dol1/sum(dol1) as share
	from tmp1
	group by seg_income,  ymonth_date
	order by seg_income;  
quit;

proc sgpanel data = tmp2;
	panelby seg_income;
	vbar ymonth_date / group = channel_type response = share; 
	*colaxis values=('01JAN2004'd to '01DEC2012'd by year);
	colaxis FITPOLICY=THIN interval=Year;
	title "Aggregate channel share over time";
run;
	
* --------------------------------------------------------*;
* Trend of shopping strategies;
* Number of shopping days; 
%raw_mean_trend(mydata(where=(seg_income^='')), seg_income, ymonth_date, num_day);
%raw_mean_trend(mydata(where=(seg_income^='')), seg_income, ymonth_date, num_trip);
%raw_mean_trend(mydata(where=(seg_income^='')), seg_income, ymonth_date, Q);
%raw_mean_trend(mydata(where=(seg_income^='')), seg_income, ymonth_date, n_mod);

* Expenditure per trips;
proc sql noprint; 
	create table tmp_trips as 
	select A.*, B.seg_income
	from mainex.all_trips as A inner join tmp_panelists as B 
	on A.household_code = B.household_code
	order by seg_income, ymonth;
quit;

%raw_mean_trend(tmp_trips, seg_income, ymonth, dol_purchases);

ods pdf close; 
ods listing; 

proc datasets noprint; delete tmp tmp1 tmp2; run;

**************************;
* Within household trend *;
**************************;
* Within household income change; 
proc sort data=panelists; by seg_income household_code; run;
%trend_reg(panelists, income_midvalue, year_num, year_num, seg_income, household_code, tmp);
data tmp1;
	set tmp;
	where Estimate^=0;
	year_num = input(substr(compress(Parameter, '','A'),2), 8.);
	y = Estimate + mean; 
	up 	= y + 1.96*StdErr;
	low = y - 1.96*StdErr;
proc sort; by year_num;
run;

data tmp1;
	merge tmp1 year_fmt;
	by year_num;
run;

proc sgplot data=tmp1(where=(Estimate^=.));
	band x = Year upper = up lower = low/group = seg_income transparency = 0.6;
	series x = Year y = y / group=seg_income;
	title "Within-household trend of income";
run;

* Monthly expenditure changes; 
proc sort data=mydata; by household_code; run;

%trend_reg(mydata, dol_purchases, month ymonth_num, month ymonth_num, seg_income, household_code, myfit);

* Expenditure share; 
%trend_reg(mydata, share_Convenience_Store, month ymonth_num, month ymonth_num, seg_income, household_code, tmp);	
data myfit; set myfit tmp; run;

%trend_reg(mydata, share_Discount_Store, month ymonth_num, month ymonth_num, seg_income, household_code, tmp);	
data myfit; set myfit tmp; run;

%trend_reg(mydata, share_Dollar_Store, month ymonth_num, month ymonth_num, seg_income, household_code, tmp);	
data myfit; set myfit tmp; run;

%trend_reg(mydata, share_Drug_Store, month ymonth_num, month ymonth_num, seg_income, household_code, tmp);	
data myfit; set myfit tmp; run;

%trend_reg(mydata, share_Grocery, month ymonth_num, month ymonth_num, seg_income, household_code, tmp);	
data myfit; set myfit tmp; run;

%trend_reg(mydata, share_Warehouse_Club, month ymonth_num, month ymonth_num, seg_income, household_code, tmp);	
data myfit; set myfit tmp; run;
 
* Shopping strategies; 
%trend_reg(mydata, num_day, month ymonth_num, month ymonth_num, seg_income, household_code, tmp);	
data myfit; set myfit tmp; run;

%trend_reg(mydata, num_trip, month ymonth_num, month ymonth_num, seg_income, household_code, tmp);	
data myfit; set myfit tmp; run;

%trend_reg(mydata, Q, month ymonth_num, month ymonth_num, seg_income, household_code, tmp);	
data myfit; set myfit tmp; run;

%trend_reg(mydata, Q_comp, month ymonth_num, month ymonth_num, seg_income, household_code, tmp);	
data myfit; set myfit tmp; run;

%trend_reg(mydata, Q_dry, month ymonth_num, month ymonth_num, seg_income, household_code, tmp);	
data myfit; set myfit tmp; run;

%trend_reg(mydata, n_mod, month ymonth_num, month ymonth_num, seg_income, household_code, tmp);	
data myfit; 
	set myfit tmp; 
	freq = "month";
run;

* --------------------------------------------------------*;
* Expenditure and quantity at trip level;
proc sql noprint; 
	create table tmp as 
	select A.*, B.start as ymonth_num, (food_quant + nonedible_quant) as Q
	from tmp_trips as A left join ymonth_num_fmt as B
	on A.year=B.year and A.month=B.month
	order by household_code;
quit;
data tmp_trips; set tmp; run;

%trend_reg(tmp_trips, dol_purchases, month ymonth_num, month ymonth_num, seg_income, household_code, tmp);	
data tmp; set tmp; Dependent = "trip_dol"; freq = "trip"; run;
data myfit; set myfit tmp; run;

%trend_reg(tmp_trips, Q, month ymonth_num, month ymonth_num, seg_income, household_code, tmp);	
data tmp; set tmp; Dependent = "trip_Q"; freq = "trip"; run;
data myfit; set myfit tmp; run;

***********************************;
* Size at household-product level *;
***********************************;
/* I run 6 regressions to verify my operation of size index, at each purchase: 
	UNITS(no vs. yes) * product size operation (absolute size, within-category size index, across-category size index)
*/

* Pick up 20 categories that take up wallet share; 
data tmp; set mylib.module_wallet; run;
proc sort data=tmp; by product_module_descr; run;
proc univariate data=tmp noprint; 
	by product_module_descr;
	var share_paid; 
	output out = tmp1 mean=share_mean P50=share_med; 
run;

* Order the wallet share and compute the cumulative share; 
proc sort data=tmp1; by descending share_med; run;
* Select the top 100 product modules but drop magnet and reference card; 
data tmp2;
	set tmp1; 
	where (product_module_descr ^= "MAGNET DATA") & (substr(product_module_descr, 1, 14) ^= "REFERENCE CARD");
	if _N_ > 100 then delete; 
	cum_share + share_med;
run;

* Pick one category; 
data tmp1; 
	set mylib.module_wallet;
	where product_module_descr = "DETERGENTS - HEAVY DUTY - LIQUID";
	*where product_module_descr = "DAIRY-MILK-REFRIGERATED";
proc sort nodupkeys out=tmp2; by product_module_descr;
run;

*-------------------------------------*

* Use purchases data, add household-category index, year-month number variable and household segment variable; 
proc sql noprint;
	* NOTE: comment out the where statement if we include all the products; 
	create table tmp3 as 
	select *
	from mainex.purchases(drop=retailer_code) 
	where product_module_descr in (select product_module_descr from tmp2);

	create table tmp as 
	select A.*, cats(household_code, "-", product_module_descr) as HH_module, size_index*share_paid as across_size, 
			quantity*size as size_quant, quantity*size_index as size_index_quant, quantity*size_index*share_paid as across_size_quant,
			B.start as ymonth_num
	from tmp3 as A left join ymonth_num_fmt as B
	on A.year = B.year and A.month = B.month;
	
	create table mydata_purchases as 
	select A.*, B.seg_income
	from tmp as A left join tmp_panelists as B
	on A.household_code = B.household_code
	order by HH_module, year, month;
quit;

* Save the data; 
* data mainex.mydata_purchases; set mydata_purchases; run;

* --------------------------------------------------------*;
* Pure product size selection, i.e., ignoring quantity; 
%trend_reg(mydata_purchases, size, month ymonth_num, month ymonth_num, seg_income, HH_module, tmp);	
data tmp; set tmp; Dependent = "purchase_size"; freq = "purchases"; run;
data myfit; set myfit tmp; run;

%trend_reg(mydata_purchases, size_index, month ymonth_num, month ymonth_num, seg_income, HH_module, tmp);	
data tmp; set tmp; Dependent = "purchase_size_index"; freq = "purchases"; run;
data myfit; set myfit tmp; run;

%trend_reg(mydata_purchases, across_size, month ymonth_num, month ymonth_num, seg_income, HH_module, tmp);	
data tmp; set tmp; Dependent = "purchase_across_size"; freq = "purchases"; run;
data myfit; set myfit tmp; run;

* --------------------------------------------------------*;
* Account for units, size* quantity; 
%trend_reg(mydata_purchases, size_quant, month ymonth_num, month ymonth_num, seg_income, HH_module, tmp);	
data tmp; set tmp; Dependent = "purchase_size_quant"; freq = "purchases"; run;
data myfit; set myfit tmp; run;

%trend_reg(mydata_purchases, size_index_quant, month ymonth_num, month ymonth_num, seg_income, HH_module, tmp);	
data tmp; set tmp; Dependent = "purchase_size_index_quant"; freq = "purchases"; run;
data myfit; set myfit tmp; run;

%trend_reg(mydata_purchases, across_size_quant, month ymonth_num, month ymonth_num, seg_income, HH_module, tmp);	
data tmp; set tmp; Dependent = "purchase_across_size_quant"; freq = "purchases"; run;
data myfit; set myfit tmp; run;

* --------------------------------------------------------*;
* Plot the within-household changes; 
ods listing close;
ods pdf file = "&plotfile2";
%trendplot(myfit);
ods pdf close; 
ods listing;

* Save the regression results; 
data mylib.output_3_estimation; set myfit; run;

*****************************;
* Direct regression testing *;
*****************************;
data mydata;
	set mydata;
	ln_income_midvalue = log(income_midvalue);
run;

* Baseline: share ~ income + price + month; 
proc glm data=mydata; 
	absorb household_code; 
	class month; 
	model share_Discount_Store = PRICE_Discount_Store ln_income_midvalue month/ solution noint;
	title "Baselline fixed effect regression of expenditure share at discount stores";
	ods output ParameterEstimates=tmp1;
run;
data tmp; length model $10. parameter $50.; set tmp1; model = "Base"; run;

* Year trend; 
proc glm data=mydata; 
	absorb household_code; 
	class month year; 
	model share_Discount_Store = PRICE_Discount_Store ln_income_midvalue year month/ solution noint;
	title "Model 2: fixed effect regression of expenditure share at discount stores";
	ods output ParameterEstimates=tmp2;
run;
data tmp2;  length model $10. parameter $50.; set tmp2; model = "year"; run;
data tmp;  set tmp tmp2; run;

* Add shopping strategies;
proc glm data=mydata; 
	absorb household_code; 
	class month year; 
	model share_Discount_Store = PRICE_Discount_Store ln_income_midvalue Q num_trip year month/ solution noint;
	title "Model 3: fixed effect regression of expenditure share at discount stores";
	ods output ParameterEstimates=tmp3;
run;
data tmp3; length model $10. parameter $50.; set tmp3; model = "shopping"; run;
data tmp; set tmp tmp3; run;

* Report the parameters together; 
proc format;           
	picture stderrf (round)       
     	low-high=' 9.9999)' (prefix='(')                                
            .=' ';                                                  
run;

proc tabulate data=tmp out = share_reg noseps;      
  class model Parameter;                      
  var estimate stderr;     
  table Parameter=''*(estimate =' '*sum=' ' stderr=' '*sum=' '*F=stderrf.),                
         model=' '/ box=[label="Parameter"] rts=15 row=float misstext=' '; 
run;



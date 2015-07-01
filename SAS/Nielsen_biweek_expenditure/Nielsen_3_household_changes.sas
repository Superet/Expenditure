libname mylib "E:\Users\ccv103\Documents\Research\Store switching\SAS_temp\Nielsen";
libname mysel "E:\Users\ccv103\Documents\Research\Store switching\SAS_temp\Nielsen\Work";
libname mainex "E:\Users\ccv103\Documents\Research\Store switching\SAS_temp\Nielsen\Work\main";

options compress=yes reuse=yes; 
proc datasets library=WORK kill noprint; run;

%let dirname= E:\Users\ccv103\Documents\Research\Nielsen\nielsen_extracts\HMS;
%let plotfile1 = E:\Users\ccv103\Desktop\graph_raw_trend.pdf;
%let plotfile2 = E:\Users\ccv103\Desktop\graph_within_trend.pdf;
%let plotfile2l = E:\Users\ccv103\Desktop\graph_within_trend_loess.pdf;

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
%macro trend_reg(data, dv, iv, classv, groupvar, fixed, levelvar, base, newdata, outname);
	proc sql noprint; 
		select unique(&groupvar) into :grouplevel separated by ' ' from &data; 
		select count(unique(&groupvar)) into :n from &data; 
	quit; 
	
	%do i=1 %to &n;
		proc glm data=&data(where = (&groupvar="%scan(&grouplevel, &i, ' ')" ));
			absorb &fixed; 
			class &classv; 
			model &dv = &iv / solution noint;
			store est_save;
/*			ods output ParameterEstimates=tmp&i;*/
			title "Fixed effect regression of &dv for &groupvar = %scan(&grouplevel, &i, ' ')";
		run;
		
		proc plm restore=est_save; 
			score data=&newdata out = tmp&i predicted stderr LCLM UCLM;; 
		run;
		data tmp&i; length Dependent $30. &groupvar $15.; set tmp&i; Dependent = "&dv"; &groupvar = "%scan(&grouplevel, &i, ' ')";
	%end;
	
	%let tmplast = tmp%trim(&n);
	data &outname; set tmp1-&tmplast; run;
	proc univariate data=&data(where=(&levelvar=&base)) noprint;
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
* Lowess plots of within household trend; 
%macro trend_loess(data, dv, iv1, iv2, classv, groupvar, fixed, levelvar, base, newdata, outname); 
	proc sql noprint; 
		select unique(&groupvar) into :grouplevel separated by ' ' from &data; 
		select count(unique(&groupvar)) into :n from &data; 
	quit; 
	
	proc sql noprint; 
		create table tmpdata as
		select *, (&dv - mean(&dv)) as y
		from &data
		group by &fixed;
	quit; 
	
	%do i=1 %to &n;
		proc glm data=tmpdata(where = (&groupvar="%scan(&grouplevel, &i, ' ')" )) noprint;
			class &classv; 
			model y = &iv1 / solution;
			output out = res RESIDUAL=res;
		run;
		
		proc loess data=res plots=none; 
			model res = &iv2 /degree=2 smooth=&mysmooth; 
			score data=&newdata / CLM; 
			ods output ScoreResults = tmp&i; 
		run;
		data tmp&i; length Dependent $30. &groupvar $15.; set tmp&i(drop=ScoreData obs); Dependent = "&dv"; &groupvar = "%scan(&grouplevel, &i, ' ')";
		proc datasets noprint; delete res; run; 
	%end;
	
	%let tmplast = tmp%trim(&n);
	data &outname; set tmp1-&tmplast; run;
	proc univariate data=&data(where=(&levelvar=&base)) noprint;
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
%mend trend_loess; 

*-------------------------*;
* Plot the trend; 
%macro trendplot(data, yvar, levelvar, lowvar, upvar, plotcl);
	proc sql noprint;
		select unique(Dependent) into : deplevel separated by ' ' from &data;
		select count(unique(Dependent)) into :n from &data;
	quit;
	
	%if &plotcl=1 %then %do; 
		data tmp; 
			set &data; 
			y 	= &yvar + &levelvar; 
			low = &lowvar + &levelvar;
			up  = &upvar + &levelvar; 
		run;
		
		%do i=1 %to &n;
			proc sgplot data=tmp(where=( Dependent = "%scan(&deplevel, &i, ' ')"));
				band x = biweek lower=low upper = up/group = seg_income transparency=0.5;
				series x = biweek y = y/group = seg_income transparency=0.5; 
				title "Trend plot for %scan(&deplevel, &i, ' ')";
			run;
		%end;
	%end;	
	%else %do; 
		data tmp; 
			set &data; 
			y 	= &yvar + &levelvar; 
		run;
		
		%do i=1 %to &n;
			proc sgplot data=tmp(where=( Dependent = "%scan(&deplevel, &i, ' ')"));
				series x = biweek y = y/group = seg_income transparency=0.5; 
				title "Trend plot for %scan(&deplevel, &i, ' ')";
			run;
		%end;
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
	if income_real <= 15 then seg_income = 'Qt1';
	if income_real>15 and income_real<=19 then seg_income = 'Qt2';
	if income_real>19 and income_real<=23 then seg_income = 'Qt3';
	if income_real>23 then seg_income = 'Qt4';
run;

proc sql noprint;
	create table tmp as
	select A.*, B.seg_income
	from (select * from mylib.all_panelists where income_real ^=. )as A left join tmp_panelists as B
	on A.household_code = B.household_code;
	
	create table panelists as 
	select A.*, B.year_num
	from tmp as A left join year_fmt as B
	on A.panel_year = B.year; 
quit;

*--------------------------*;
* Create the regression data;
data tmp; 
	set mylib.all_hh_biweek_exp_merge;
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
proc transpose data=mylib.Format_biweek_price out=tmp_price prefix = PRICE_; 
	by scantrack_market_descr year biweek;
	id channel_type;
	var price_paid_norm_index;
run;

proc sql noprint;
	* Add CPI to normalize the prices; 
	select Annual into: price_base from cpi where year = 2004; 
	
	create table tmp2 as 
	select A.*, B.Annual/&price_base as cpi, A.unit_price/(B.Annual/&price_base) as unit_price2004
	from tmp as A left join cpi as B
	on A.year = B.year
	order by household_code; 
	
	* Add price data;
	create table tmp3(drop = old_mkt old_year old_biweek) as 
	select A.*, B.*
	from tmp2 as A left join tmp_price(drop=_NAME_ _LABEL_ rename=(scantrack_market_descr=old_mkt year=old_year biweek=old_biweek)) as B
	on A.scantrack_market_descr=B.old_mkt and A.year =B.old_year and A.biweek=B.old_biweek;
	
	* Add the segment variable; 
	create table mydata as
	select A.*, B.seg_income
	from tmp3 as A left join tmp_panelists as B
	on A.household_code = B.household_code
	order by household_code, year, biweek; 
quit;

* Compute quantity via expenditure and aggregate price; 
data mydata;
	set mydata;
/*	Q_comp = 0;
	array dolvar{*} dolp_:;
	array prcvar{*} price_:;
	do i=1 to dim(dolvar);
		Q_comp = Q_comp + dolvar[i]/prcvar[i];
	end;
	Q_dry = Q - magnet_quant;*/
	first_date = mdy(1,1,2004) + (biweek-1)*14; 
	month = month(first_date); 
run;

proc contents data=mydata; run;

* Create biweek data for prediction; 
proc sql noprint;
	create table biweek_dat as 
	select distinct biweek
	from mydata;
quit;

data biweek_dat; 
	set biweek_dat; 
	month = 1; 
run;
proc sql noprint; select max(biweek) into: base from biweek_dat ; quit; 
%put &base; 

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
proc sort data=mydata; by seg_income biweek; run;

ods listing close;
ods pdf file = "&plotfile1";

%raw_mean_trend(mydata, seg_income, biweek, dol_purchases);
	
* --------------------------------------------------------*;
* Outcome: trend of expenditure share ;
proc sql noprint;
	create table tmp as 
	select seg_income, biweek, 
		sum(DOLP_Convenience_Store) as DOLP_Convenience_Store, sum(DOLP_Discount_Store) as DOLP_Discount_Store, 
		sum(DOLP_Dollar_Store) as DOLP_Dollar_Store, sum(DOLP_Drug_Store) as DOLP_Drug_Store, sum(DOLP_Grocery) as DOLP_Grocery, 
		sum(DOLP_Warehouse_Club) as DOL_Warehouse_Club
	from mydata
	group by seg_income, biweek
	order by seg_income, biweek;
quit;

proc transpose data=tmp out=tmp1(rename=(_NAME_=channel_type)) prefix=dol;
	by seg_income biweek;
	var DOLP_Convenience_Store DOLP_Discount_Store DOLP_Dollar_Store DOLP_Drug_Store DOLP_Grocery DOL_Warehouse_Club;
run;

proc sql noprint; 
	create table tmp2 as 
	select seg_income, biweek, channel_type, dol1/sum(dol1) as share
	from tmp1
	group by seg_income,  biweek
	order by seg_income;  
quit;

proc sgpanel data = tmp2;
	panelby seg_income;
	vbar biweek / group = channel_type response = share; 
	*colaxis values=('01JAN2004'd to '01DEC2012'd by year);
	colaxis FITPOLICY=THIN interval=Year;
	title "Aggregate channel share over time";
run;
	
* --------------------------------------------------------*;
* Trend of shopping strategies;
* Number of shopping days; 
%raw_mean_trend(mydata(where=(seg_income^='')), seg_income, biweek, num_day);
%raw_mean_trend(mydata(where=(seg_income^='')), seg_income, biweek, num_trip);
%raw_mean_trend(mydata(where=(seg_income^='')), seg_income, biweek, Q);
%raw_mean_trend(mydata(where=(seg_income^='')), seg_income, biweek, n_mod);

* Expenditure per trips;
proc sql noprint; 
	create table tmp_trips as 
	select A.*, B.seg_income
	from mainex.all_trips as A inner join tmp_panelists as B 
	on A.household_code = B.household_code
	order by seg_income, biweek;
quit;

%raw_mean_trend(tmp_trips, seg_income, biweek, dol_purchases);

ods pdf close; 
ods listing; 

proc datasets noprint; delete tmp tmp1 tmp2; run;

**************************;
* Within household trend *;
**************************;
* Within household income change; 
data tmp(drop=i); 
	do i=2004 to 2010; 
		panel_year = i; output; 
	end; 
run;

proc sort data=panelists; by seg_income household_code; run;
%trend_reg(panelists, income_midvalue, panel_year, panel_year, seg_income, household_code, panel_year, 2004, tmp, tmpout);

data tmpout; 
	set tmpout; 
	y 	= Predicted + mean; 
	low	= LCLM + mean; 
	up 	= UCLM + mean; 
run;

proc sgplot data=tmp1(where=(Estimate^=.));
	band x = panel_year upper = up lower = low/group = seg_income transparency = 0.6;
	series x = panel_year y = y / group=seg_income;
	title "Within-household trend of income";
run;

* --------------------------------------------------------*;
* Fixed effects of biweek; 
* --------------------------------------------------------*;
* Bi-weekly expenditure changes; 
proc sort data=mydata; by household_code; run;

%trend_reg(mydata, dol_purchases, month biweek, month biweek, seg_income, household_code, biweek, &base, biweek_dat, myfit);

* Expenditure share; 
%trend_reg(mydata, share_Convenience_Store, month biweek, month biweek, seg_income, household_code, biweek, &base, biweek_dat,tmp);	
data myfit; set myfit tmp; run;

%trend_reg(mydata, share_Discount_Store, month biweek, month biweek, seg_income, household_code, biweek, &base, biweek_dat,tmp);	
data myfit; set myfit tmp; run;

%trend_reg(mydata, share_Dollar_Store, month biweek, month biweek, seg_income, household_code, biweek, &base, biweek_dat,tmp);	
data myfit; set myfit tmp; run;

%trend_reg(mydata, share_Drug_Store, month biweek, month biweek, seg_income, household_code, biweek, &base, biweek_dat,tmp);	
data myfit; set myfit tmp; run;

%trend_reg(mydata, share_Grocery, month biweek, month biweek, seg_income, household_code, biweek, &base, biweek_dat,tmp);	
data myfit; set myfit tmp; run;

%trend_reg(mydata, share_Warehouse_Club, month biweek, month biweek, seg_income, household_code, biweek, &base, biweek_dat,tmp);	
data myfit; set myfit tmp; run;
 
* Shopping strategies; 
%trend_reg(mydata, num_day, month biweek, month biweek, seg_income, household_code, biweek, &base, biweek_dat,tmp);	
data myfit; set myfit tmp; run;

%trend_reg(mydata, num_trip, month biweek, month biweek, seg_income, household_code, biweek, &base, biweek_dat,tmp);	
data myfit; set myfit tmp; run;

%trend_reg(mydata, Q, month biweek, month biweek, seg_income, household_code, biweek, &base, biweek_dat,tmp);	
data myfit; set myfit tmp; run;

%trend_reg(mydata, Q_comp, month biweek, month biweek, seg_income, household_code, biweek, &base, biweek_dat,tmp);	
data myfit; set myfit tmp; run;

%trend_reg(mydata, Q_dry, month biweek, month biweek, seg_income, household_code, biweek, &base, biweek_dat,tmp);	
data myfit; set myfit tmp; run;

%trend_reg(mydata, n_mod, month biweek, month biweek, seg_income, household_code, biweek, &base, biweek_dat,tmp);	
data myfit; 
	set myfit tmp; 
	freq = "biweek";
run;

* --------------------------------------------------------*;
* Expenditure and quantity at trip level;
data tmp_trips; 
	set tmp_trips; 
	Q = food_quant + nonedible_quant;
run; 

%trend_reg(tmp_trips, dol_purchases, month biweek, month biweek, seg_income, household_code, biweek, &base, biweek_dat, tmp);	
data tmp; set tmp; Dependent = "trip_dol"; freq = "trip"; run;
data myfit; set myfit tmp; run;

%trend_reg(tmp_trips, Q, month biweek, month biweek, seg_income, household_code, biweek, &base, biweek_dat,tmp);	
data tmp; set tmp; Dependent = "trip_Q"; freq = "trip"; run;
data myfit; set myfit tmp; run;

* --------------------------------------------------------*;
* Cubic effects of biweek; 
* --------------------------------------------------------*;
%trend_reg(mydata, dol_purchases, month biweek|biweek|biweek, month , seg_income, household_code, biweek, &base, biweek_dat, myfit_cubic);

* Expenditure share; 
%trend_reg(mydata, share_Convenience_Store, month biweek|biweek|biweek, month, seg_income, household_code, biweek, &base, biweek_dat,tmp);	
data myfit_cubic; set myfit_cubic tmp; run;

%trend_reg(mydata, share_Discount_Store, month biweek|biweek|biweek, month, seg_income, household_code, biweek, &base, biweek_dat,tmp);	
data myfit_cubic; set myfit_cubic tmp; run;

%trend_reg(mydata, share_Dollar_Store, month biweek|biweek|biweek, month, seg_income, household_code, biweek, &base, biweek_dat,tmp);	
data myfit_cubic; set myfit_cubic tmp; run;

%trend_reg(mydata, share_Drug_Store, month biweek|biweek|biweek, month, seg_income, household_code, biweek, &base, biweek_dat,tmp);	
data myfit_cubic; set myfit_cubic tmp; run;

%trend_reg(mydata, share_Grocery, month biweek|biweek|biweek, month, seg_income, household_code, biweek, &base, biweek_dat,tmp);	
data myfit_cubic; set myfit_cubic tmp; run;

%trend_reg(mydata, share_Warehouse_Club, month biweek|biweek|biweek, month, seg_income, household_code, biweek, &base, biweek_dat,tmp);	
data myfit_cubic; set myfit_cubic tmp; run;
 
* Shopping strategies; 
%trend_reg(mydata, num_day, month biweek|biweek|biweek, month, seg_income, household_code, biweek, &base, biweek_dat,tmp);	
data myfit_cubic; set myfit_cubic tmp; run;

%trend_reg(mydata, num_trip, month biweek|biweek|biweek, month, seg_income, household_code, biweek, &base, biweek_dat,tmp);	
data myfit_cubic; set myfit_cubic tmp; run;

%trend_reg(mydata, Q, month biweek|biweek|biweek, month, seg_income, household_code, biweek, &base, biweek_dat,tmp);	
data myfit_cubic; set myfit_cubic tmp; run;

%trend_reg(mydata, Q_comp, month biweek|biweek|biweek, month, seg_income, household_code, biweek, &base, biweek_dat,tmp);	
data myfit_cubic; set myfit_cubic tmp; run;

%trend_reg(mydata, Q_dry, month biweek|biweek|biweek, month, seg_income, household_code, biweek, &base, biweek_dat,tmp);	
data myfit_cubic; set myfit_cubic tmp; run;

%trend_reg(mydata, n_mod, month biweek|biweek|biweek, month, seg_income, household_code, biweek, &base, biweek_dat,tmp);	
data myfit_cubic; 
	set myfit_cubic tmp; 
	freq = "biweek";
run;

* --------------------------------------------------------*;
* Expenditure and quantity at trip level;
%trend_reg(tmp_trips, dol_purchases, month biweek|biweek|biweek, month, seg_income, household_code, biweek, &base, biweek_dat, tmp);	
data tmp; set tmp; Dependent = "trip_dol"; freq = "trip"; run;
data myfit_cubic; set myfit_cubic tmp; run;

%trend_reg(tmp_trips, Q, month biweek|biweek|biweek, month, seg_income, household_code, biweek, &base, biweek_dat,tmp);	
data tmp; set tmp; Dependent = "trip_Q"; freq = "trip"; run;
data myfit_cubic; set myfit_cubic tmp; run;

* --------------------------------------------------------*;
* Loess plot of the within-household trend; 
* --------------------------------------------------------*;
%let mysmooth = 0.5; 
%let base= 1; 
%trend_loess(mydata, dol_purchases, month, biweek, month, seg_income, household_code, biweek, &base, biweek_dat, myfit_loess);

* Expenditure share; 
%trend_loess(mydata, share_Convenience_Store, month, biweek, month, seg_income, household_code, biweek, &base, biweek_dat,tmp);	
data myfit_loess; set myfit_loess tmp; run;

%trend_loess(mydata, share_Discount_Store, month, biweek, month, seg_income, household_code, biweek, &base, biweek_dat,tmp);	
data myfit_loess; set myfit_loess tmp; run;

%trend_loess(mydata, share_Dollar_Store, month, biweek, month, seg_income, household_code, biweek, &base, biweek_dat,tmp);	
data myfit_loess; set myfit_loess tmp; run;

%trend_loess(mydata, share_Drug_Store, month, biweek, month, seg_income, household_code, biweek, &base, biweek_dat,tmp);	
data myfit_loess; set myfit_loess tmp; run;

%trend_loess(mydata, share_Grocery, month, biweek, month, seg_income, household_code, biweek, &base, biweek_dat,tmp);	
data myfit_loess; set myfit_loess tmp; run;

%trend_loess(mydata, share_Warehouse_Club, month, biweek, month, seg_income, household_code, biweek, &base, biweek_dat,tmp);	
data myfit_loess; set myfit_loess tmp; run;
 
* Shopping strategies; 
%trend_loess(mydata, num_day, month, biweek, month, seg_income, household_code, biweek, &base, biweek_dat,tmp);	
data myfit_loess; set myfit_loess tmp; run;

%trend_loess(mydata, num_trip, month, biweek, month, seg_income, household_code, biweek, &base, biweek_dat,tmp);	
data myfit_loess; set myfit_loess tmp; run;

%trend_loess(mydata, Q, month, biweek, month, seg_income, household_code, biweek, &base, biweek_dat,tmp);	
data myfit_loess; set myfit_loess tmp; run;

%trend_loess(mydata, Q_comp, month, biweek, month, seg_income, household_code, biweek, &base, biweek_dat,tmp);	
data myfit_loess; set myfit_loess tmp; run;

%trend_loess(mydata, Q_dry, month, biweek, month, seg_income, household_code, biweek, &base, biweek_dat,tmp);	
data myfit_loess; set myfit_loess tmp; run;

%trend_loess(mydata, n_mod, month, biweek, month, seg_income, household_code, biweek, &base, biweek_dat,tmp);	
data myfit_loess; 
	set myfit_loess tmp; 
	freq = "biweek";
run;

* --------------------------------------------------------*;
* Expenditure and quantity at trip level;
%trend_loess(tmp_trips, dol_purchases, month, biweek, month, seg_income, household_code, biweek, &base, biweek_dat, tmp);	
data tmp; set tmp; Dependent = "trip_dol"; freq = "trip"; run;
data myfit_loess; set myfit_loess tmp; run;

%trend_loess(tmp_trips, Q, month, biweek, month, seg_income, household_code, biweek, &base, biweek_dat,tmp);	
data tmp; set tmp; Dependent = "trip_Q"; freq = "trip"; run;
data myfit_loess; set myfit_loess tmp; run;


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
	if _N_ > 10 then delete; 
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
	create table tmp as 
	select *
	from mainex.purchases(drop=retailer_code) 
	where product_module_descr in (select product_module_descr from tmp2);
	
	create table mydata_purchases as 
	select A.*, cats(household_code, "-", product_module_descr) as HH_module, size_index*share_paid as across_size, 
			quantity*size as size_quant, quantity*size_index as size_index_quant, quantity*size_index*share_paid as across_size_quant,
			B.seg_income
	from tmp as A left join tmp_panelists as B
	on A.household_code = B.household_code
	order by HH_module, year, month;
quit;

* Save the data; 
* data mainex.mydata_purchases; set mydata_purchases; run;

* --------------------------------------------------------*;
* Pure product size selection, i.e., ignoring quantity; 
%trend_reg(mydata_purchases, size, month biweek, month biweek, seg_income, HH_module, biweek, &base, biweek_dat,tmp);	
data myfit1; set tmp; Dependent = "purchase_size"; freq = "purchases"; run;

%trend_reg(mydata_purchases, size_index, month biweek, month biweek, seg_income, HH_module, biweek, &base, biweek_dat,tmp);	
data tmp; set tmp; Dependent = "purchase_size_index"; freq = "purchases"; run;
data myfit1; set myfit1 tmp; run;

%trend_reg(mydata_purchases, across_size, month biweek, month biweek, seg_income, HH_module, biweek, &base, biweek_dat, tmp);	
data tmp; set tmp; Dependent = "purchase_across_size"; freq = "purchases"; run;
data myfit1; set myfit1 tmp; run;

* --------------------------------------------------------*;
* Account for units, size* quantity; 
%trend_reg(mydata_purchases, size_quant, month biweek, month biweek, seg_income, HH_module, biweek, &base, biweek_dat,tmp);	
data tmp; set tmp; Dependent = "purchase_size_quant"; freq = "purchases"; run;
data myfit1; set myfit1 tmp; run;

%trend_reg(mydata_purchases, size_index_quant, month biweek, month biweek, seg_income, HH_module, biweek, &base, biweek_dat,tmp);	
data tmp; set tmp; Dependent = "purchase_size_index_quant"; freq = "purchases"; run;
data myfit1; set myfit1 tmp; run;

%trend_reg(mydata_purchases, across_size_quant, month biweek, month biweek, seg_income, HH_module, biweek, &base, biweek_dat,tmp);	
data tmp; set tmp; Dependent = "purchase_across_size_quant"; freq = "purchases"; run;
data myfit1; set myfit1 tmp; run;

* --------------------------------------------------------*;
* Plot the within-household changes; 
ods listing close;
ods pdf file = "&plotfile2";
%trendplot(myfit, Estimate, mean, LCLM, UCLM, 1);
%trendplot(myfit_cubic,Estimate, mean, LCLM, UCLM, 1); 
ods pdf close; 
ods listing;

* --------------------------------------------------------*;
* Plot loess of the within-household changes; 
ods listing close; 
ods pdf file = "&plotfile2l";
%trendplot(myfit_loess, p_res, mean, LCLM, UCLM, 0);
ods pdf close; 
ods listing;

* Save the regression results; 
data mylib.output_3_estimation_biweek; set myfit; run;
data mylib.output_3_estimation_biweek_loess; set myfit_loess; run;

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

data mylib.share_reg; set share_reg; run;

roc sort data=reg_share; by parameter; 
proc transpose data=reg_share out=est;
by parameter ;
id model; 
var Estimate_Sum;
run;

proc transpose data=reg_share out=std;
by parameter ;
id model; 
var StdErr_Sum;
run;

data tmp; 
	merge est std; 
	by parameter _NAME_; 
run;


* Mediation analysis; 
* Baseline: share ~ income + price + month; 
proc contents data=mydata; run;
data mydata; 
	set mydata;
	ln_dol = log(dol_purchases); 
	ln_income_midvalue = log(income_midvalue);
run;

proc glm data=mydata; 
	absorb household_code; 
	class month ; 
	model share_Discount_Store = ln_income_midvalue recession month/ solution noint;
	title "Baselline fixed effect regression of expenditure share at discount stores";
	ods output ParameterEstimates=tmp1;
run;
data tmp; length model $10. parameter $50.; set tmp1; model = "indirect"; run;

* Direct steps; 
proc glm data=mydata; 
	absorb household_code; 
	class month; 
	model share_Discount_Store = ln_dol num_day month/ solution noint;
	ods output ParameterEstimates=tmp1;
run;
data tmp1; length model $10. parameter $50.; set tmp1; model = "direct"; run; 
data tmp; set tmp tmp1; run;

* Mediation step; 
proc glm data=mydata; 
	absorb household_code; 
	class month; 
	model share_Discount_Store = ln_income_midvalue recession ln_dol num_day month/ solution noint;
	ods output ParameterEstimates=tmp1;
run;
data tmp1; length model $10. parameter $50.; set tmp1; model = "mediation"; run; 
data tmp; set tmp tmp1; run;
data tmp; set tmp; interaction = 0; run;

* Interacte with income groups; 
proc glm data=mydata; 
	absorb household_code; 
	class month income_group; 
	model share_Discount_Store = ln_income_midvalue*income_group recession*income_group month/ solution noint;
	ods output ParameterEstimates=tmp1;
run;
data tmp1; length model $10. parameter $50.; set tmp1; model = "indirect"; interaction = 1; run;
data tmp; set tmp tmp1; run;

* Direct steps; 
proc glm data=mydata; 
	absorb household_code; 
	class month income_group; 
	model share_Discount_Store = ln_dol*income_group num_day*income_group month/ solution noint;
	ods output ParameterEstimates=tmp1;
run;
data tmp1; length model $10. parameter $50.; set tmp1; model = "direct"; interaction = 1; run; 
data tmp; set tmp tmp1; run;

* Mediation step; 
proc glm data=mydata; 
	absorb household_code; 
	class month income_group; 
	model share_Discount_Store = ln_income_midvalue*income_group recession*income_group ln_dol*income_group num_day*income_group month/ solution noint;
	ods output ParameterEstimates=tmp1;
run;
data tmp1; length model $10. parameter $50.; set tmp1; model = "mediation"; interaction = 1; run; 
data tmp; set tmp tmp1; run;

* Export results; 
PROC EXPORT DATA=tmp
			OUTFILE= "\\tsclient\Resear1\Store switching\processed data\tmpresults.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;






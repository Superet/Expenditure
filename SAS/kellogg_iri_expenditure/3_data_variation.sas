libname mylib "E:\Users\ccv103\Documents\Research\Store switching\SAS_temp\RetailFormat";
libname ORIG "E:\Users\ccv103\Documents\Research\Store switching\SAS_temp\Orig";

%macro multiple_glm(set_data,outname);
	proc sql noprint; select count(n) into: num_model from &set_data;quit;
	%do i=1 %to &num_model;
		proc sql noprint;
			select data, dv, absorbvar, classvar, X, model_label
			into :data_m,:dv_m, :absorbvar_m,:classvar_m, :X_m, :model_labelm
			from &set_data where n=&i;
		quit;
		%put &data_m, &dv_m, &absorbvar_m, &classvar_m, &model_labelm;
		
		proc glm data=&data_m order=data;
			absorb &absorbvar_m;
			class &classvar_m;
			model &dv_m = &X_m / solution;
			ods output ParameterEstimates=tmp&i ;
		run;
		
		data tmp&i;
			set tmp&i;
			modelnum = &i;
			modelname = "&model_labelm";
		run;
	%end;
	
	%let tmplast = tmp%trim(&num_model);
	%put &tmplast;
	data &outname;
		length Dependent $30.;
		length Parameter $30.;
		set tmp1-&tmplast;
	run;
	
	proc datasets noprint; delete tmp1-&tmplast; run;
%mend multiple_glm;

* Format standard error;
proc format;           
	picture stderrf (round)       
     	low-high=' 9.99)' (prefix='(') .=' ';    
	picture stderrfh (round)       
     	low-high=' 9.9999)' (prefix='(') .=' ';                                              
run;

*****************;
* Organize data *;
*****************;
*** Monthly expenditure data ***;
proc sql noprint;
	create table exp_month as 
	select A.*, B.INCOME, B.FAMSIZE
	from mylib.pan_month as A left join ORIG.demo as B
	on A.PANID = B.PANID
	order by PANID, YMONTH;
quit;

data exp_month;
	set exp_month;
	if INCOME^=13;
	INCOME1 = 'High';
	if INCOME>7 & INCOME<=10 then INCOME1 = "Med";
	if INCOME<=7 then INCOME1 = "Low";
	FAMSIZE1 = 'Single';
	if FAMSIZE = 2 then FAMSIZE1 = 'Two';
	if FAMSIZE > 2 then FAMSIZE1 = 'Three+';
	recession = 0;
	if YMONTH >= mdy(12,1,2007) then recession = 1;
	ln_BSKTDOL = log(BSKTDOL);
	array dolvar{*} DOL_:;
	do i = 1 to dim(dolvar);
		dolvar[i] = dolvar[i]/DOL;
	end;
	drop i; 
	rename	DOL_Convenience 				= SHARE_Convenience
		   	DOL_Conventional_supermarkets 	= SHARE_Conventional_supermarkets
		  	DOL_Dollar_variety				= SHARE_Dollar_variety 
			DOL_Drug						= SHARE_Drug 		
			DOL_Limited_assortment_stores 	= SHARE_Limited_assortment_stores
			DOL_Mass						= SHARE_Mass 		
			DOL_Super_stores 				= SHARE_Super_stores
			DOL_Wholesale_clubs				= SHARE_Wholesale_clubs
			DOL_Other						= SHARE_Other
			;
run;

proc sql noprint;
	create table tmp1 as
	select MKT,PANID,INCOME1,FAMSIZE1,recession, count(YMONTH) as N_obs, std(BSKTDOL) as std_bsktdol, 
			log(std(BSKTDOL)) as ln_std_bsktdol
	from exp_month
	group by MKT,PANID,INCOME1,FAMSIZE1,recession;
	
	create table vol_month as
	select *
	from tmp1 
	where PANID in (select PANID from tmp1 
					group by PANID
					having sum(N_obs>=6) = 2);
quit;

proc means data=vol_month N MIN MEAN STD P25 P50 P75 P90 MAX; var std_bsktdol ln_std_bsktdol; run;	
proc univariate data=vol_month noprint; histogram std_bsktdol; run;

*** Weekly data ***;
proc sql noprint;
	create table tmp as 
	select MKT,PANID, min(WEEK) as min_WEEK, max(WEEK) as max_WEEK
	from mylib.pan_week 
	group by MKT,PANID
	order by PANID;
	
	create table tmp1 as
	select A.*,B.FAMSIZE, B.INCOME
	from tmp as A left join ORIG.demo as B
	on A.PANID = B.PANID
	order by PANID;
quit;

data tmp1;
	set tmp1;
	by PANID;
	WEEK = min_WEEK - 1;
	if first.PANID then do; 	
		do while(WEEK<max_WEEK);
			WEEK + 1; output;
		end;
	end;
	drop min_week max_week;	
run;
data tmp1;
	set tmp1;
	DATE = (WEEK - 400)*7 + 31900 - 21916;
	FORMAT DATE YYMMDD10.;
	MONTH 	= month(DATE);
	YMONTH	= mdy(month(DATE),1,year(DATE));
	FORMAT YMONTH monyy7.;
run;

proc sql noprint;
	create table exp_week(drop=old_panid old_week) as  
	select *
	from tmp1 as A left join mylib.pan_week(rename=(PANID=old_panid WEEK=old_week) drop=DATE MKT) as B
	on A.PANID = B.old_panid and A.WEEK = B.old_week
	order by PANID,WEEK;
quit;

data exp_week;
	set exp_week;
	ln_bsktdol = log(BSKTDOL);
	array allnum{*} BSKTDOL--QUANT;
	do i = 1 to dim(allnum);
		if allnum[i] = . then allnum[i] = 0;
	end;
	drop i;
	if INCOME^=13;
	INCOME1 = 'High';
	if INCOME>7 & INCOME<=10 then INCOME1 = "Med";
	if INCOME<=7 then INCOME1 = "Low";
	FAMSIZE1 = 'Single';
	if FAMSIZE = 2 then FAMSIZE1 = 'Two';
	if FAMSIZE > 2 then FAMSIZE1 = 'Three+';
	recession = 0;
	if YMONTH >= mdy(12,1,2007) then recession = 1;
run;

* Volitility of weekly expenditure;
proc sql noprint;
	create table tmp as
	select PANID,INCOME1,FAMSIZE1,recession, count(WEEK) as N_obs, std(BSKTDOL) as std_bsktdol
	from exp_week
	group by PANID,INCOME1,FAMSIZE1,recession
	order by PANID;
	
	create table vol_week as 
	select *
	from tmp
	where PANID in (select PANID from tmp
					group by PANID
					having sum(N_obs>=10) = 2);
	
	create table tmp1 as
	select PANID,INCOME1,FAMSIZE1,recession, count(WEEK) as N_obs, std(BSKTDOL) as std_bsktdol
	from (select * from exp_week where BSKTDOL>0)
	group by PANID,INCOME1,FAMSIZE1,recession
	order by PANID;	
	
	create table vol_week_excl0 as 
	select *
	from tmp1
	where PANID in (select PANID from tmp1 
					group by PANID
					having sum(N_obs>=10) = 2);
quit;

proc means data=vol_week N MIN MEAN STD P25 P50 P75 P90 MAX; var std_bsktdol ; run;	
proc means data=vol_week_excl0 N MIN MEAN STD P25 P50 P75 P90 MAX; var std_bsktdol ; run;	
proc univariate data=vol_week noprint; histogram std_bsktdol;
proc univariate data=vol_week_excl0 noprint; histogram std_bsktdol;
run;

proc datasets noprint; delete tmp tmp1 ; run;

***********************;
* Export data ;
PROC EXPORT DATA= exp_month
            OUTFILE= "\\tsclient\Resear1\Store switching\processed data\exp_month.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;

PROC EXPORT DATA= exp_week
            OUTFILE= "\\tsclient\Resear1\Store switching\processed data\exp_week.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;

PROC EXPORT DATA= exp_month
            OUTFILE= "E:\Users\ccv103\Desktop\exp_month.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;

PROC EXPORT DATA= exp_week
            OUTFILE= "E:\Users\ccv103\Desktop\exp_week.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;

***********************;
/*Expenditure change */
***********************;
/*
Change in expenditure 
Y (Dependent variables): monthly expenditure, monthly expenditure volitility, weekly expenditure, weekly expenditure volitility
X (Independent variables): only time dummies, time dummies + demographics 
T (Trend): recession dummy, year-month dummies
*/
proc sort data=exp_month;	by PANID MONTH descending YMONTH; 
proc sort data=exp_week;	by PANID MONTH descending YMONTH; 
proc sort data=vol_month;	by PANID descending receesion; 
proc sort data=vol_week;	by PANID descending receesion; 
run;

data model_run;
	input n data $6-22 dv $23-35 absorbvar $36-50 classvar $51-65 X $66-100 model_label $101-110;
	datalines;
	1	exp_month		BSKTDOL		PANID MONTH		recession			recession 						Y1X0T0
	2	exp_month		ln_BSKTDOL	PANID MONTH		recession			recession 						Y1lX0T0
	3	vol_month		std_BSKTDOL	PANID 			recession			recession 						Y2X0T0
	4	exp_week		BSKTDOL		PANID MONTH		recession			recession 						Y3X0T0
	5	exp_week		ln_BSKTDOL	PANID MONTH		recession			recession 						Y3lX0T0
	6	vol_week		std_BSKTDOL PANID 			recession			recession 						Y4X0T0
	6	vol_week_excl0	std_BSKTDOL PANID 			recession			recession 						Y4e0X0T0
	;
run;

%multiple_glm(model_run,est1);

data est1;
	set est1;
	where StdErr ^= .;
	level = 'Monthly';
	if substr(modelname,1,2) = "Y3" | substr(modelname,1,2) = "Y4" then level='Weekly';
	if substr(modelname,1,3) = "Y4e" then level='WeeklyE0';
run;

* Use year-month trend;
data model_run;
	input n data $6-20 dv $21-32 absorbvar $33-45 classvar $46-65 X $66-100 model_label $101-110;
	datalines;
	1	exp_month	BSKTDOL		PANID 			YMONTH				YMONTH	 						Y1X0T1
	2	exp_month	ln_BSKTDOL	PANID 			YMONTH				YMONTH	 						Y1lX0T1
	3	exp_week	BSKTDOL		PANID 			YMONTH				YMONTH	 						Y3X0T1
	4	exp_week	ln_BSKTDOL	PANID 			YMONTH				YMONTH	 						Y3lX0T1
	;
run;

%multiple_glm(model_run,est2);

data est2;
	set est2;
	where StdErr ^= .;
	level = 'Monthly';
	if substr(modelname,1,2) = "Y3" | substr(modelname,1,2) = "Y4" then level='Weekly';
run;

* Use recession*demographics trend;
data model_run;
	input n data $6-20 dv $21-32 absorbvar $33-45 classvar $46-65 X $66-100 model_label $101-110;
	datalines;
	1	exp_month	BSKTDOL		PANID MONTH		recession INCOME1	recession recession*INCOME1		Y1X1T0
	2	exp_month	ln_BSKTDOL	PANID MONTH		recession INCOME1	recession recession*INCOME1		Y1lX1T0
	3	vol_month	std_BSKTDOL	PANID 			recession INCOME1	recession recession*INCOME1		Y2X1T0
	4	exp_week	BSKTDOL		PANID MONTH		recession INCOME1	recession recession*INCOME1		Y3X1T0
	5	exp_week	ln_BSKTDOL	PANID MONTH		recession INCOME1	recession recession*INCOME1		Y3lX1T0
	6	vol_week	std_BSKTDOL PANID 			recession INCOME1	recession recession*INCOME1		Y4X1T0
	;
run;

%multiple_glm(model_run,est3);

data est3;
	set est3;
	where StdErr ^= .;
	level = 'Monthly';
	if substr(modelname,1,2) = "Y3" | substr(modelname,1,2) = "Y4" then level='Weekly';
run;

* Use year-month*demographics trend;
data model_run;
	input n data $6-20 dv $21-32 absorbvar $33-45 classvar $46-65 X $66-100 model_label $101-110;
	datalines;
	1	exp_month	BSKTDOL		PANID 			YMONTH INCOME1		YMONTH YMONTH*INCOME1			Y1X1T1
	2	exp_month	ln_BSKTDOL	PANID 			YMONTH INCOME1		YMONTH YMONTH*INCOME1			Y1lX1T1
	3	exp_week	BSKTDOL		PANID 			YMONTH INCOME1		YMONTH YMONTH*INCOME1			Y3X1T1
	4	exp_week	ln_BSKTDOL	PANID 			YMONTH INCOME1		YMONTH YMONTH*INCOME1			Y3lX1T1
	;
run;

%multiple_glm(model_run,est4);

data est4;
	set est4;
	where StdErr ^= .;
	level = 'Monthly';
	if substr(modelname,1,2) = "Y3" | substr(modelname,1,2) = "Y4" then level='Weekly';
run;

*****************************;
* Expenditure share changes *;
*****************************;
/*
Change in expenditure share 
X (Independent variables): only time dummies, time dummies + demographics 
T (Trend): recession dummy, year-month dummies
*/

* Only recession dummies;
data model_run;
	input n data $6-22 dv $23-55 absorbvar $56-65 classvar $66-85 X $86-105 model_label $106-135;
	datalines;
	1	exp_month		SHARE_Convenience					PANID MONTH		recession			recession 						X0T0
	2	exp_month		SHARE_Conventional_supermarkets		PANID MONTH		recession			recession 						X0T0
	3	exp_month		SHARE_Dollar_variety				PANID MONTH		recession			recession 						X0T0
	4	exp_month		SHARE_Drug							PANID MONTH		recession			recession 						X0T0
	5	exp_month		SHARE_Limited_assortment_stores		PANID MONTH		recession			recession 						X0T0
	6	exp_month		SHARE_Mass							PANID MONTH		recession			recession 						X0T0
	7	exp_month		SHARE_Super_stores					PANID MONTH		recession			recession 						X0T0
	8	exp_month		SHARE_Wholesale_clubs				PANID MONTH		recession			recession 						X0T0	
	9	exp_month		SHARE_Other							PANID MONTH		recession			recession 						X0T0	
	;
run;

%multiple_glm(model_run,share1);

data share1;
	set share1;
	where StdErr ^= .;
	level = 'Monthly';
	trend = 'recession';
	Dependent = substr(Dependent,7);
run;

* Only year-month dummies;
data model_run;
	input n data $6-22 dv $23-55 absorbvar $56-65 classvar $66-85 X $86-105 model_label $106-135;
	datalines;
	1	exp_month		SHARE_Convenience					PANID 			YMONTH				YMONTH 							X0T1
	2	exp_month		SHARE_Conventional_supermarkets		PANID 			YMONTH				YMONTH 							X0T1
	3	exp_month		SHARE_Dollar_variety				PANID 			YMONTH				YMONTH 							X0T1
	4	exp_month		SHARE_Drug							PANID 			YMONTH				YMONTH 							X0T1
	5	exp_month		SHARE_Limited_assortment_stores		PANID 			YMONTH				YMONTH 							X0T1
	6	exp_month		SHARE_Mass							PANID 			YMONTH				YMONTH 							X0T1
	7	exp_month		SHARE_Super_stores					PANID 			YMONTH				YMONTH 							X0T1
	8	exp_month		SHARE_Wholesale_clubs				PANID 			YMONTH				YMONTH 							X0T1
	9	exp_month		SHARE_Other							PANID 			YMONTH				YMONTH 							X0T1
	;
run;

%multiple_glm(model_run,share2);

data share2;
	set share2;
	where StdErr ^= .;
	level = 'Monthly';
	trend = 'Year-month';
	Dependent = substr(Dependent,7);
run;

* recession dummies*demographics ;
data model_run;
	input n data $6-22 dv $23-55 absorbvar $56-75 classvar $76-95 X $96-125 model_label $126-135;
	datalines;
	1	exp_month		SHARE_Convenience					PANID MONTH		recession INCOME1	recession recession*INCOME1		X1T0
	2	exp_month		SHARE_Conventional_supermarkets		PANID MONTH		recession INCOME1	recession recession*INCOME1		X1T0
	3	exp_month		SHARE_Dollar_variety				PANID MONTH		recession INCOME1	recession recession*INCOME1		X1T0
	4	exp_month		SHARE_Drug							PANID MONTH		recession INCOME1	recession recession*INCOME1		X1T0
	5	exp_month		SHARE_Limited_assortment_stores		PANID MONTH		recession INCOME1	recession recession*INCOME1		X1T0
	6	exp_month		SHARE_Mass							PANID MONTH		recession INCOME1	recession recession*INCOME1		X1T0
	7	exp_month		SHARE_Super_stores					PANID MONTH		recession INCOME1	recession recession*INCOME1		X1T0
	8	exp_month		SHARE_Wholesale_clubs				PANID MONTH		recession INCOME1	recession recession*INCOME1		X1T0
	9	exp_month		SHARE_Other							PANID MONTH		recession INCOME1	recession recession*INCOME1		X1T0
	;
run;

%multiple_glm(model_run,share3);

data share3;
	set share3;
	where StdErr ^= .;
	level = 'Monthly';
	trend = 'recession*demo';
	Dependent = substr(Dependent,7);
run;

* year-month dummies * demo;
data model_run;
	input n data $6-22 dv $23-55 absorbvar $56-65 classvar $66-90 X $91-120 model_label $121-135;
	datalines;
	1	exp_month		SHARE_Convenience					PANID 			YMONTH INCOME1		YMONTH YMONTH*INCOME1			X1T1
	2	exp_month		SHARE_Conventional_supermarkets		PANID 			YMONTH INCOME1		YMONTH YMONTH*INCOME1			X1T1
	3	exp_month		SHARE_Dollar_variety				PANID 			YMONTH INCOME1		YMONTH YMONTH*INCOME1			X1T1
	4	exp_month		SHARE_Drug							PANID 			YMONTH INCOME1		YMONTH YMONTH*INCOME1			X1T1
	5	exp_month		SHARE_Limited_assortment_stores		PANID 			YMONTH INCOME1		YMONTH YMONTH*INCOME1			X1T1
	6	exp_month		SHARE_Mass							PANID 			YMONTH INCOME1		YMONTH YMONTH*INCOME1			X1T1
	7	exp_month		SHARE_Super_stores					PANID 			YMONTH INCOME1		YMONTH YMONTH*INCOME1			X1T1
	8	exp_month		SHARE_Wholesale_clubs				PANID 			YMONTH INCOME1		YMONTH YMONTH*INCOME1			X1T1
	9	exp_month		SHARE_Other							PANID 			YMONTH INCOME1		YMONTH YMONTH*INCOME1			X1T1
	;
run;

%multiple_glm(model_run,share4);

data share4;
	set share4;
	where StdErr ^= .;
	level = 'Monthly';
	trend = 'Year-month*demo';
	Dependent = substr(Dependent,7);
run;

*******************;
* Output delivery *;
*******************;
*** Plot the trend ***; 
ods pdf file = 'E:\Users\ccv103\Desktop\tmp_graph.pdf';

*** Plot the year-month trend ***;
data tmpplot;
	set est2(keep=Parameter Dependent Estimate level);
	time = input(cats('01',substr(Parameter,9)),date9.);
	FORMAT time date9.;
proc sort; by level Dependent time;
run;
title 'Trend of average monthly expenditure';
proc sgplot data=tmpplot(where=( Dependent = 'BSKTDOL' & level = 'Monthly'));
	series x=time y=Estimate;
run;
title 'Trend of average weekly expenditure';
proc sgplot data=tmpplot(where=( Dependent = 'BSKTDOL' & level = 'Weekly'));
	series x=time y=Estimate;
run;

* Average trend of expenditure share;
data tmp;
	set share2(keep=Parameter Dependent Estimate level);
	time = input(cats('01',substr(Parameter,9)),date9.);
	FORMAT time date9.;
proc sort; by level time Dependent;
run;

title 'Average trend of expenditure share';
proc sgplot data=tmp;
	series x=time y=estimate / group = Dependent;
run;

* Trend of expenditure by demo;
data tmp1;
	set share4(keep=Parameter Dependent Estimate level);
	time= input(cats('01',scan(Parameter,2,' ')),date9.);
	Income = scan(Parameter,3,' ');
	if Income =' ' then Income = 'Ref';
	FORMAT time date9.;
proc sort; by Dependent time Income;
run;

proc transpose data=tmp1 out=tmpplot;
	by Dependent time;
	id Income;
	var Estimate;
run;

data tmpplot;
	set tmpplot;
	Low = Ref + Low;
	Med = Ref + Med;
run;

title 'Trend of expenditure share for low-income';
proc sgplot data=tmpplot;
	series x=time y=Low / group = Dependent;
run;

title 'Trend of expenditure share for Med-income';
proc sgplot data=tmpplot;
	series x=time y=Med / group = Dependent;
run;

title 'Trend of expenditure share for High-income';
proc sgplot data=tmpplot;
	series x=time y=Ref / group = Dependent;
run;

ods pdf close;
	
*** Output estimates ***;
ods RTF file = 'E:\Users\ccv103\Desktop\tmp_regression.rtf';
title 'Expenditure regression with recession dummy only';
proc tabulate data=est1 noseps;      
  class Parameter level Dependent;                      
  var Estimate StdErr;     
  table Parameter=''*(estimate =' '*sum=' '                          
                     stderr=' '*sum=' '*F=stderrf.),                
         level=' '*Dependent=' ' / box=[label="Parameter"] rts=15 row=float misstext=' '; 
run;

title 'Expenditure regression with year-month dummies';
proc tabulate data=est2 noseps;      
  class Parameter level Dependent;                      
  var Estimate StdErr;     
  table Parameter=''*(estimate =' '*sum=' '                          
                     stderr=' '*sum=' '*F=stderrf.),                
         level=' '*Dependent=' ' / box=[label="Parameter"] rts=15 row=float misstext=' '; 
run;

title 'Expenditure regression with recession dummy and demographics';
proc tabulate data=est3 noseps;      
  class Parameter level Dependent;                      
  var Estimate StdErr;     
  table Parameter=''*(estimate =' '*sum=' '                          
                     stderr=' '*sum=' '*F=stderrf.),                
         level=' '*Dependent=' ' / box=[label="Parameter"] rts=15 row=float misstext=' '; 
run;

title 'Expenditure regression with year-month dummies and demongraphics';
proc tabulate data=est4 noseps;      
  class Parameter level Dependent;                      
  var Estimate StdErr;     
  table Parameter=''*(estimate =' '*sum=' '                          
                     stderr=' '*sum=' '*F=stderrf.),                
         level=' '*Dependent=' ' / box=[label="Parameter"] rts=15 row=float misstext=' '; 
run;

title 'Expenditure share regression with recession dummy only';
proc tabulate data=share1 noseps;      
  class Parameter Dependent;                      
  var Estimate StdErr;     
  table Parameter=''*(estimate =' '*sum=' '*F=d5.4                          
                     stderr=' '*sum=' '*F=stderrfh.),                
        Dependent=' ' / box=[label="Parameter"] rts=15 row=float misstext=' '; 
run;

title 'Expenditure share regression with year-month dummy only';
proc tabulate data=share2 noseps;      
  class Parameter Dependent;                      
  var Estimate StdErr;     
  table Parameter=''*(estimate =' '*sum=' '*F=d5.4                          
                     stderr=' '*sum=' '*F=stderrfh.),                
        Dependent=' ' / box=[label="Parameter"] rts=15 row=float misstext=' '; 
run;

title 'Expenditure share regression with recession dummy and demo';
proc tabulate data=share3 noseps;      
  class Parameter Dependent;                      
  var Estimate StdErr;     
  table Parameter=''*(estimate =' '*sum=' '*F=d5.4                          
                     stderr=' '*sum=' '*F=stderrfh.),                
        Dependent=' ' / box=[label="Parameter"] rts=15 row=float misstext=' '; 
run;

title 'Expenditure share regression with year-month dummies and demo';
proc tabulate data=share4 noseps;      
  class Parameter Dependent;                      
  var Estimate StdErr;     
  table Parameter=''*(estimate =' '*sum=' '*F=d5.4                          
                     stderr=' '*sum=' '*F=stderrfh.),                
        Dependent=' ' / box=[label="Parameter"] rts=15 row=float misstext=' '; 
run;

ods RTF close;

*********************;
* Basket regression *;
*********************;
* Merge monthly expenditure with weekly basket; 
proc sql noprint;
	create table tmp as
	select A.*,B.BSKTDOL as month_bsktdol
	from exp_week as A left join exp_month as B
	on A.PANID = B.PANID and A.YMONTH = B.YMONTH;
	
	create table tmp1 as 
	select A.*,B.std_bsktdol as month_std
	from tmp as A left join vol_month as B
	on A.PANID = B.PANID and A.recession = B.recession
	order by PANID,WEEK;
	
	create table basket as 
	select A.*, B.std_bsktdol as week_std
	from tmp1 as A left join vol_week as B
	on A.PANID = B.PANID and A.recession = B.recession
	order by PANID,WEEK;	
quit;

* Restrict to the basket conditional on shopping; 
data basket;
	set basket;
	where DOL > 0;
	ln_month_bsktdol= log(month_bsktdol);
	ln_numcat		= log(num_cat);
	r_num_ed		= num_edible / num_cat;
	r_quant_ed 		= quant_ed/QUANT;
	* Add the lag of last purchase; 
	quant_last		= lag1(QUANT);
	num_cat_last	= lag1(num_cat);
	ln_numcat_last	= lag1(ln_numcat);
	r_num_ed_last	= lag1(r_num_ed);
	r_quant_ed_last	= lag1(r_quant_ed);
	by PANID;
	if first.PANID then do;
		quant_last 		= .;
		numcat_last 	= .;
		ln_numcat_last 	= .;
		r_num_ed_last	= .;
		r_quant_ed_last	= .;
	end;
run;

* Basket quantity against monthly expenditure and volatility;
proc sort data=basket; by PANID MONTH; run;	
proc glm data=basket(where=(QUANT>0)) order=data;
	absorb PANID MONTH;
	class recession INCOME1;
	model QUANT = recession recession*INCOME1 / solution;
run;

proc glm data=basket order=data;
	absorb PANID MONTH;
	class INCOME1;
	model QUANT = ln_month_bsktdol month_std month_std*INCOME1 quant_last/ solution;
run;

proc glm data=basket order=data;
	absorb PANID MONTH;
	class FAMSIZE1;
	model QUANT = ln_month_bsktdol month_std month_std*FAMSIZE1 quant_last/ solution;
run;

* Basket variety;
proc glm data=basket(where = (num_cat>0)) order = data;
	absorb PANID MONTH;
	class INCOME1;
	model num_cat = ln_month_bsktdol month_std month_std*INCOME1 numcat_last / solution;
run;

proc glm data=basket(where = (num_cat>0)) order = data;
	absorb PANID MONTH;
	class INCOME1;
	model ln_numcat = ln_month_bsktdol month_std month_std*INCOME1 ln_numcat_last / solution;
run;
	
proc glm data=basket(where = (num_cat>0)) order = data;
	absorb PANID MONTH;
	class INCOME1;
	model r_quant_ed = ln_month_bsktdol month_std month_std*INCOME1 r_quant_ed_last / solution;
run;

proc glm data=basket(where = (num_cat>0)) order = data;
	absorb PANID MONTH;
	class INCOME1;
	model r_num_ed = ln_month_bsktdol month_std month_std*INCOME1 r_num_ed_last / solution;
run;
	
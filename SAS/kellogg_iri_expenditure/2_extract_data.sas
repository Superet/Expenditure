libname mylib "E:\Users\ccv103\Documents\Research\Store switching\SAS_temp\RetailFormat";
libname ORIG "E:\Users\ccv103\Documents\Research\Store switching\SAS_temp\Orig";
%let dirname= E:\Users\Projects\Project_Marketing_IRI\kelloggirir;

******************;
* Data selection *;
******************;
/*
 Household selection ;
 Criterion: 
 1. HH must have stayed at least 6 months before and during the recession
 2. HH must has the spell between trips less than 30 days
 3. HH must have non-missing demographics 
 4. HH's weekly expenditure is not too extreme
*/

* Set parameters;
%let min_stay 		= 300;
%let max_spell 		= 30;
%let edible_dpt		= 5;

proc sort data=ORIG.trip_sum; by PANID DATE; run;
proc means data=ORIG.trip_sum noprint;
	by PANID DATE;
	var BSKTDOL;
	output out=tmp MAX=;
run;

proc sort data=tmp; by PANID DATE; run;
data tmp;
	set tmp;
	recession = 0;
	if DATE >= MDY(12,1,2007) then recession = 1;
	interpurchase_day = dif(DATE);
	by PANID;
	if first.PANID then interpurchase_day = .;
run;

proc sql noprint;
	create table tmp2 as 
	select PANID, recession, count(DATE) as num_trip, max(DATE) - min(DATE) as stay_length,
			max(interpurchase_day) as max_spell, 
			max(BSKTDOL) as max_BSKTDOL
	from tmp 
	group by PANID,recession;

	create table tmp1 as 
	select PANID, mean(recession) as mean_rec , min(num_trip) as recession_trip, min(stay_length) as stay_length,
			max(max_spell) as max_spell, 
			max(max_BSKTDOL) as max_BSKTDOL
	from tmp2 
	group by PANID;

	create table pan_char(drop=old_PANID) as 
	select A.*, B.*
	from tmp1 as A left join ORIG.demo(rename=(PANID=old_PANID)) as B
	on A.PANID = B.old_PANID
	order by PANID;
quit;

proc means data=pan_char N MEAN MEDIAN MIN MAX P1 P25 P75 P99;
	var recession_trip -- max_BSKTDOL;
run;
proc means data=pan_char noprint;
	var max_BSKTDOL;
	output out=tmp P99=;
run;

data tmp;
	set tmp;
	call symput ('BSKTDOL_cut',max_BSKTDOL);
run;

data pan_char;
	set pan_char;
	keep = 1;
	if (mean_rec=0 | mean_rec=1) then keep = 0;
	if stay_length < &min_stay then keep = 0;
	if max_spell >&max_spell then keep = 0;
	if max_BSKTDOL > &BSKTDOL_cut then keep = 0;
	array demovar{14} INCOME FAMSIZE RACE CHILDREN FMLE_AGE FMLE_EDU FMLE_HRS FMLE_OCC 
					MALE_AGE MALE_EDU MALE_HRS MALE_OCC M_STATUS RENTOWN;
	do i = 1 to 14;
		if demovar[i] = 99 then keep = 0;
	end;
run;
proc freq data=pan_char; tables keep; run;

data keeppan;
	set pan_char;
	if keep = 1;
	keep PANID;
run;

****************;
* Combine data *;
****************;
/** Group convenience and drug ;
data tmp;
	set ORIG.trip_det;
	if FORMAT = 'Convenience' | FORMAT = 'Drug' then FORMAT = 'Convenience_Drug';
run;*/

* Expenditure on dry goods in each format;
proc sql noprint;
	create table tmp1 as 
	select *
	from ORIG.trip_det
	where PANID in (select PANID from keeppan);

	create table tmp2 as 
	select A.*, B.CATEGORY, B.BRAND,B.vol_unit, B.DEPARTMENT, (B.DEPARTMENT<=&edible_dpt) as Edible, (B.DEPARTMENT>&edible_dpt) as NonEdible
	from tmp1 as A left join mylib.upc_unit as B
	on A.UPC = B.UPC;
	
	create table trip_long as
	select MKT,PANID, WEEK, FORMAT, sum(DOL) as DOL, sum(UNITS*vol_unit) as QUANT
	from tmp2 
	group by MKT, PANID, WEEK, FORMAT;

	create table tmp3 as
	select MKT,PANID,WEEK,DEPARTMENT,Edible,NonEdible,CATEGORY, sum(UNITS*vol_unit) as QUANT
	from tmp2 
	group by MKT,PANID,WEEK,DEPARTMENT,Edible,NonEdible,CATEGORY;
		
	create table tmp_cat as 
	select MKT,PANID,WEEK, count(distinct CATEGORY) as num_cat, sum(Edible) as num_edible, sum(NonEdible) as num_nonedible,
			sum(QUANT*Edible) as quant_ed, sum(QUANT*NonEdible) as quant_noned
	from tmp3
	group by MKT,PANID,WEEK
	order by MKT,PANID,WEEK;
quit;

proc transpose data=trip_long out=tmp_dol(drop=_name_) prefix=DOL_;
	by MKT PANID WEEK;
	id FORMAT;
	var DOL;
run;
proc transpose data=trip_long out=tmp_quant(drop=_name_) prefix=QUANT_;
	by MKT PANID WEEK;
	id FORMAT;
	var QUANT;
run;
proc sort data=mylib.fm_price; by MKT WEEK FORMAT; run;
proc transpose data=mylib.fm_price out=tmp_price(drop=_name_) prefix=PRICE_;
	by MKT WEEK;
	id FORMAT;
	var unit_price;
run;

data tmp_trip;
	merge tmp_dol tmp_quant tmp_cat;
	by MKT PANID WEEK;
	array allnum{*} _numeric_;
	array dolvar{*} DOL_:;
	array qntvar{*} QUANT_:;
/*	array dolvar{*} DOL_Convenience DOL_Conventional_supermarkets DOL_Dollar_variety 
					DOL_Drug DOL_Limited_assortment_stores DOL_Mass DOL_Super_stores DOL_Wholesale_clubs;
	array qntvar{*} QUANT_Convenience QUANT_Conventional_supermarkets QUANT_QUANTlar_variety 
					QUANT_Drug QUANT_Limited_assortment_stores QUANT_Mass QUANT_Super_stores QUANT_Wholesale_clubs; */
	do i = 1 to dim(allnum);
		if allnum[i] = . then allnum[i] = 0;
	end;
	DOL = 0;
	QUANT = 0;
	do i = 1 to dim(dolvar);
		DOL = DOL + dolvar[i];
		QUANT = QUANT + qntvar[i];
	end;
	drop i;
run;

/*proc sort data=tmp_trip; by MKT WEEK; run;
data tmp_trip;
	merge tmp_trip tmp_price;
	by MKT WEEK;
run;*/

* Merge together with total basket spending;
proc sql noprint;
	create table tmp1 as 
	select PANID, WEEK, sum(BSKTDOL) as BSKTDOL
	from ORIG.trip_sum
	where PANID in (select PANID from keeppan)
	group by PANID, WEEK;
	
	create table tmp2 as 
	select A.*,B.MKT
	from tmp1 as A left join ORIG.demo as B
	on A.PANID=B.PANID;

	create table mylib.pan_week(drop=old_PANID old_WEEK) as 
	select *
	from tmp2 as A left join tmp_trip(rename=(PANID=old_PANID WEEK=old_WEEK) drop=MKT) as B
	on A.PANID=B.old_PANID and A.WEEK = B.old_WEEK
	order by PANID, WEEK;
quit;

* Correct the baskt dollar;
data mylib.pan_week;
	set mylib.pan_week;
	if BSKTDOL <= DOL then BSKTDOL = DOL + 1;
	DATE = (WEEK - 400)*7 + 31900 - 21916;
	FORMAT DATE YYMMDD10.;
run;

proc datasets noprint; delete tmp tmp1 tmp2 tmp_dol tmp_quant tmp_price tmp_cat; run;

*** Monthly expenditure ***;
data tmp;
	set mylib.pan_week;
	YMONTH = mdy(month(DATE),1,year(DATE));
	MONTH = month(DATE);
	FORMAT YMONTH monyy7.;
run;

proc sql noprint;
	create table mylib.pan_month as 
	select PANID, MKT, YMONTH, MONTH, sum(BSKTDOL) as BSKTDOL, sum(DOL) as DOL, sum(QUANT) as QUANT,
		sum(DOL_Convenience) as DOL_Convenience, sum(DOL_Conventional_supermarkets) as DOL_Conventional_supermarkets,
		sum(DOL_Dollar_variety) as DOL_Dollar_variety, sum(DOL_Drug) as DOL_Drug,
		sum(DOL_Limited_assortment_stores) as DOL_Limited_assortment_stores, sum(DOL_Mass) as DOL_Mass,
		sum(DOL_Natural_gourmet_food_stor) as DOL_Natural_gourmet_food_stor, sum(DOL_Other) as DOL_Other,
		sum(DOL_Super_stores) as DOL_Super_stores,sum(DOL_Wholesale_clubs) as DOL_Wholesale_clubs
	from tmp
	group by MKT,PANID, YMONTH, MONTH
	order by PANID,YMONTH;
quit;


**********************;
* Charaterize basket *;
**********************;
data pan_week; 
	set mylib.pan_week;
	where MKT = 10 & DOL > 0;
run;

proc univariate data=pan_week noprint; 
	histogram QUANT; 
	histogram num_cat; 
	histogram num_edible;
	histogram num_nonedible;
	histogram quant_ed;
run;
proc means data=pan_week N MEAN P1 P25 P50 P75 P99;
var QUANT num_cat num_edible num_nonedible quant_ed quant_noned; 
run;

%let quant_cut = 25;
data pan_week;
	set pan_week;
	edible_pct = quant_ed/(quant_ed + quant_noned);
	if edible_pct <=.8 then variety = 'Non-edible';
	else if edible_pct < 1 then variety = 'Both';
	else variety = 'Edible';
	if QUANT <=20 then size = 'Small';
	else if QUANT <= 50 then size = 'Med';
	else size = 'Large';
run;

data tmp;
	set pan_week(keep = PANID WEEK variety num_cat size QUANT 
					DOL_Convenience DOL_Conventional_supermarkets DOL_Dollar_variety DOL_Drug 
					DOL_Limited_assortment_stores DOL_Mass DOL_Natural_gourmet_food_stor
					DOL_Other DOL_Super_stores DOL_Wholesale_clubs);
run;
proc sort data=tmp; by PANID WEEK;
proc transpose data=tmp out=tmp1;
	by PANID WEEK variety num_cat size QUANT;
run;
data tmp1;
	set tmp1;
	where col1 > 0 ;
	rename _NAME_ = FORMAT 
		   COL1 = DOL;
	grp = variety||size;
run;

*** Verify basket classification: frequency tables of format trips ***;
proc sort data=tmp1; by variety size;run;
proc freq data=tmp1;
	by grp;
	tables FORMAT/out = tmp2;
run;

proc sort data=tmp2; by grp percent; run;
proc gchart data=tmp2;
	hbar grp/ discrete subgroup=FORMAT freq=percent type=percent;
run;

******************;
* Export to .csv *;
******************;
* Export a subset *;
%let sel_mkt = 10;
data tmp;
	set mylib.pan_week(drop=QUANT_:);
/*	where MKT = &sel_mkt;*/
run;

PROC EXPORT DATA= tmp
            OUTFILE= "\\tsclient\Resear1\Store switching\processed data\0_subset.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;

PROC EXPORT DATA=mylib.format_mkt_static
			OUTFILE= "\\tsclient\Resear1\Store switching\processed data\0_format_mkt_static.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;

PROC EXPORT DATA=mylib.fm_price
			OUTFILE= "\\tsclient\Resear1\Store switching\processed data\0_format_price.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;

PROC EXPORT DATA=mylib.fmt_cat
			OUTFILE= "\\tsclient\Resear1\Store switching\processed data\0_format_category.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;











	




libname mylib "E:\Users\ccv103\Documents\Research\Store switching\SAS_temp\RetailFormat";
libname ORIG "E:\Users\ccv103\Documents\Research\Store switching\SAS_temp\Orig";
%let dirname= E:\Users\Projects\Project_Marketing_IRI\kelloggirir;

*************;
* Read data *; 
*************;
* A function of reading data;
%macro readfile(dirlist,outname);
	proc sql noprint;
		select count(fname) into :num_files
		from &dirlist;
	quit;

	%do j=1 %to &num_files;
		proc sql noprint;
		select fname into :fname
		from &dirlist
		where n=&j;
		quit;
	
		%let fpath=&dirname\&fname;
		%put &fpath;
		proc import out=tmp&j
			datafile="&fpath"
			dbms=dlm replace;
			delimiter='|';
			getnames=no;
			datarow=1;
		run;
	%end;

	%let tmpend=%sysfunc(catx(%str(),tmp,&num_files));
	%put &tmpend;
	data &outname;
		set tmp1-&tmpend;
	run;

	* Delete separate files;
	proc datasets noprint; delete tmp1-&tmpend;run;	
%mend readfile;		

*** Read store format data ***;
PROC IMPORT OUT= WORK.storefm
            DATAFILE= "\\tsclient\Resear1\Store switching\processed data\Chain_format.txt"  
            DBMS=tab REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
RUN;

data ORIG.storefm;
	set storefm(drop=VARIABLE rename=(CODE=CHAIN DSC=CHAIN_NAME));
	if FORMAT in ('Warehouse stores','Supercenters','Conv. mass merchandiser','General merchandise') then FORMAT='Mass';
	*if FORMAT = 'Specialty stores' then FORMAT = 'Other';
proc sort; by CHAIN;
run;

*** Read duplicated trip file ***;
PROC IMPORT OUT= WORK.duptrip
            DATAFILE= "\\tsclient\Resear1\Store switching\processed data\duplicated trip.csv"  
            DBMS=csv REPLACE;
     GETNAMES=YES;
RUN;

*** Read UPC ***;
PROC IMPORT OUT= WORK.dict 
            DATAFILE= "E:\Users\Projects\Project_Marketing_IRI\kelloggirir\DICT.txt" 
            DBMS=dlm REPLACE;
	 delimiter='|';
     GETNAMES=NO;
     DATAROW=1; 
RUN;

data ORIG.dict;
	set dict(drop=VAR10);
	rename VAR1=UPC VAR2=UPCDSC VAR3=CATEGORY VAR4=TYPE VAR5=VENDOR VAR6=BRAND VAR7=DEPARTMENT VAR8=VOL_EQ VAR9=UNIT_MEASURE;
proc sort; by UPC;
run;

*** Read demographics ***;
PROC IMPORT OUT= WORK.demo 
            DATAFILE= "E:\Users\Projects\Project_Marketing_IRI\kelloggirir\DEMO.txt" 
            DBMS=dlm REPLACE;
	 delimiter='|';
     GETNAMES=NO;
     DATAROW=1; 
RUN;

data ORIG.demo;
	set demo(drop=VAR25);
	rename VAR1=PANID VAR2=FAMSIZE VAR3=INCOME VAR4=RACE VAR5=CHILDREN
	 	VAR6=FMLE_AGE VAR7=FMLE_EDU VAR8=FMLE_HRS VAR9=FMLE_OCC	VAR10=MALE_AGE 
		VAR11=MALE_EDU VAR12=MALE_HRS VAR13=MALE_OCC VAR14=M_STATUS VAR15=RENTOWN 
		VAR16=NUM_CATS VAR17=NUM_DOGS VAR18=REGION VAR19=MKT VAR20=PROJ09 VAR21=PROJ08 
		VAR22=PROJ07 VAR23=PROJ06 VAR24=ZIPCODE;
run;
	
*** Read demo code ***; 
PROC IMPORT OUT=democode
            DATAFILE="E:\Users\Projects\Project_Marketing_IRI\kelloggirir\Academic Household File.xls"
            DBMS=EXCELCS REPLACE;
   SHEET='Demo Code';
RUN;

data ORIG.democode;
	set democode (firstobs=3 rename=(DEMO_CODE=VARIABLE F2=START F3=LABEL));
	if START=. then delete;
run;

*** Format demographics ***;
data demo_fmt(keep = fmtname type start label);
	retain type 'N';
	set ORIG.democode end=lastrec;
	fmtname = variable;
	output;
	/*if lastrec then do;
		hlo ='O'; 
		label = 'Missing';
		output;
	end;*/
run;

proc format library=ORIG cntlin=demo_fmt; run;

/*options  fmtsearch=(ORIG);
data ORIG.demo;
	set ORIG.demo;
	format FAMSIZE FAMSIZE.;
	format INCOME INCOME.;
	format RACE RACE.;
	format CHILDREN CHILDREN.;
	format FMLE_AGE FMLE_AGE.;
	format FMLE_EDU FMLE_EDU.;
	format FMLE_HRS FMLE_HRS.;
	format FMLE_OCC FMLE_OCC.;
	format MALE_AGE MALE_AGE.;
	format MALE_EDU MALE_EDU.;
	format MALE_HRS MALE_HRS.;
	format MALE_OCC MALE_OCC.;
	format M_STATUS M_STATUS.;
	format RENTOWN RENTOWN.;
	format NUM_CATS NUM_CATS.;
	format NUM_DOGS NUM_DOGS.;
	format REGION REGION.;
	format MKT MKT.;
run;
*/
*** Read trip summary data ***;
filename dirlist pipe "dir /B &dirname\Trip_summary*.txt";
data dirlist;
	length fname $30;
	infile dirlist length=reclen;
	input fname $varying30. reclen;
	n=_n_;
	*call symput ('num_files',_n_);
run;
proc print data=dirlist;run;

%readfile(dirlist,trip_sum);

data tmp;
	set trip_sum(drop=VAR9 rename=(VAR1=PANID VAR2=CHAIN VAR3=DATE VAR4=TRIP_ID VAR5=WEEK VAR6=BSKTDOL VAR7=TRIP_MISSION VAR8=TRIP_TYPE));
	DATE = INPUT(PUT(DATE,8.),YYMMDD8.);
	FORMAT DATE YYMMDD10.;
	pdc = cats(PANID,"*",DATE,"*", CHAIN);
run;

* Put together the multiple trips in the same store on the same day;
proc sql noprint;
	create table tmp1 as
	select PANID, WEEK, DATE, CHAIN, pdc,sum(BSKTDOL) as BSKTDOL
	from tmp
	group by PANID,WEEK,DATE,CHAIN, pdc
	order by PANID,DATE,CHAIN;
	
	* drop duplicated trip; 
	create table tmp2 as 
	select *
	from tmp1
	where pdc not in (select pdc from duptrip);
	
	* Drop specialty stores;
	create table trip_sum as
	select A.*,B.FORMAT
	from tmp2 as A left join ORIG.storefm(drop=CHAIN_NAME) as B
	on A.CHAIN = B.CHAIN
	where (FORMAT is  NOT missing) and (FORMAT ne 'Specialty stores');
quit;

* Drop trips with basket size equal to .01; 
data ORIG.trip_sum; set trip_sum; if BSKTDOL>.01; drop pdc; run;

*** Read trip detail data ***;
filename dirlist pipe "dir /B &dirname\Trip_detail*.txt";
data dirlist;
	length fname $30;
	infile dirlist length=reclen;
	input fname $varying30. reclen;
	n=_n_;
	*call symput ('num_files',_n_);
run;
proc print data=dirlist;run;

%readfile(dirlist,trip_det);

data tmp;
	set trip_det(drop=VAR13 rename=(VAR1=PANID VAR2=CHAIN VAR3=DATE VAR4=TRIP_ID VAR5=WEEK 
				VAR6=UPC VAR7=DOL VAR8=UNITS VAR9=COUPON VAR10=DISPLAY VAR11=FEATUR VAR12=PRICEOFF));
	DATE = INPUT(PUT(DATE,8.),YYMMDD8.);
	FORMAT DATE YYMMDD10.;
	if UNITS<1 then UNITS=1;
run;

proc sql noprint;
* Drop specialty store;
	create table tmp1 as
	select A.*,B.FORMAT
	from tmp as A left join ORIG.storefm(drop=CHAIN_NAME) as B
	on A.CHAIN = B.CHAIN
	where (FORMAT is  NOT missing) and (FORMAT ne 'Specialty stores');
	
	create table tmp2 as 
	select A.*, B.MKT
	from tmp1 as A left join ORIG.demo as B
	on A.PANID = B.PANID
	order by PANID, DATE;
	
	create table ORIG.trip_det as 
	select A.*, B.DEPARTMENT, B.CATEGORY, B.BRAND, B.VOL_EQ
	from tmp2 as A left join ORIG.dict as B
	on A.UPC = B.UPC
	order by PANID, DATE; 
quit;

proc datasets library=work kill noprint;run;

*****************************************;
* Cross-category universal volume units *;
*****************************************;
proc sql noprint;
	create table tmp as 
	select CATEGORY, sum(DOL) as DOL, sum(UNITS*VOL_EQ) as quant
	from ORIG.trip_det
	where year(DATE)=2006
	group by CATEGORY
	order by CATEGORY;	

	select sum(dol) into :base_dol from tmp;
quit;	

proc sort data=ORIG.dict(keep= CATEGORY UNIT_MEASURE) nodupkey out=tmp1;
	by CATEGORY UNIT_MEASURE;
RUN;

data cat_vol;
	merge tmp tmp1;
	by CATEGORY;
	Base_dol = &base_dol;
	univ_unit_measure = dol/base_dol;
run;

*** Size index within category package distribution ***;
proc sort data=ORIG.dict; by CATEGORY; run;
proc means data=ORIG.dict noprint;
	by CATEGORY;
	var VOL_EQ;
	output out=tmp MEDIAN=vol_med;
run;

data tmp;
	merge tmp cat_vol;
	by CATEGORY;
run;

proc sql noprint;
	create table mylib.upc_unit as 
	select A.*, B.univ_unit_measure, VOL_EQ*B.univ_unit_measure as vol_unit1,
				B.vol_med, 100 + 100*(VOL_EQ-vol_med)/vol_med as size_index, 
				(100 + 100*(VOL_EQ-vol_med)/vol_med)*B.univ_unit_measure as vol_unit
	from ORIG.dict as A left join tmp as B 
	on A.CATEGORY = B.CATEGORY
	order by CATEGORY, BRAND;
quit;

proc datasets noprint; delete tmp tmp1 unit1; run;

*******************;
* Format offering *;
*******************;
* 2006 scanned units as weights;
proc sql noprint;
	create table trip_det_unit as 
	select A.*,B.size_index, B.vol_unit, B.univ_unit_measure
	from ORIG.trip_det as A left join mylib.upc_unit as B
	on A.UPC = B.UPC;

	create table tmp as 
	select FORMAT, MKT, DEPARTMENT,CATEGORY, BRAND, UPC, vol_unit, sum(DOL) as DOL, mean(DOL/UNITS) as PRICE, sum(UNITS) as UNITS
	from trip_det_unit
	where year(DATE) = 2006
	group by FORMAT, MKT, DEPARTMENT,CATEGORY, BRAND, UPC, vol_unit; 
	
/*	create table tmp1 as 
	select *, sum(DOL) as cat_dol, DOL/sum(DOL) as within_weight
	from tmp 
	group by FORMAT, MKT, CATEGORY;
	
	create table weight_06rev as 
	select *, sum(cat_dol) as all_dol, cat_dol/sum(cat_dol) as cat_weight
	from tmp1 
	group by FORMAT, MKT
	order by FORMAT, MKT, CATEGORY;*/
	
	create table tmp1 as 
	select *, sum(UNITS) as cat_unit, UNITS/sum(UNITS) as within_weight
	from tmp 
	group by FORMAT, MKT, DEPARTMENT, CATEGORY;
	
	create table weight_06unit as 
	select *, sum(cat_unit) as all_unit, cat_unit/sum(cat_unit) as cat_weight
	from tmp1 
	group by FORMAT, MKT
	order by FORMAT, MKT, DEPARTMENT, CATEGORY;
quit;

*** A "typical" category offering in each format ***;
proc sql noprint;
	create table tmp as
	select MKT,FORMAT,CHAIN,DEPARTMENT, CATEGORY,BRAND,UPC,vol_unit,size_index,mean(DOL/UNITS) as PRICE
	from trip_det_unit
	group by MKT,FORMAT,CHAIN,DEPARTMENT,CATEGORY,BRAND,UPC,vol_unit,size_index
	order by MKT, FORMAT, CHAIN,DEPARTMENT,CATEGORY;
	
	create table tmp1 as 
	select A.*, B.within_weight,B.cat_weight
	from tmp as A inner join weight_06unit as B
	on A.MKT=B.MKT and A.FORMAT=B.FORMAT and A.UPC=B.UPC;
	
	create table chain_cat_quant as
	select MKT,FORMAT,CHAIN,DEPARTMENT,CATEGORY,
			sum(price*within_weight)/sum(within_weight) as cat_price, 
			sum(vol_unit*within_weight)/sum(within_weight) as cat_unit, 
			sum(size_index*within_weight)/sum(within_weight) as cat_size_index,
			mean(cat_weight) as cat_weight
	from tmp1
	group by MKT,FORMAT,CHAIN,DEPARTMENT,CATEGORY;

	create table chain_quant as 
	select MKT,FORMAT,CHAIN,sum(cat_price*cat_weight)/sum(cat_weight) as price, sum(cat_unit*cat_weight)/sum(cat_weight) as unit, 
		 	sum(cat_size_index*cat_weight)/sum(cat_weight) as size_index
	from chain_cat_quant
	group by MKT,FORMAT,CHAIN
	order by MKT,FORMAT,CHAIN;
quit;

data tmp2;
	set chain_cat_quant;
	where MKT = 10;
run;

proc sgpanel data=tmp2;
	panelby FORMAT/columns = 4;
	histogram cat_unit/binwidth = 10;
	colaxis values=(0 to 1000 by 100) label = 'Unit';
run;


******************************************************;
* Format attributes: weekly prices, assortment, size *;
******************************************************;
*** Weekly prices ***;
proc sql noprint;
	create table fm_price as 
	select MKT, FORMAT, WEEK, sum(DOL) as REVENUE, sum(UNITS*vol_unit) as QAUNT, sum(DOL)/sum(UNITS*vol_unit) as unit_price
	from trip_det_unit
	group by MKT, FORMAT, WEEK
	order by MKT, FORMAT, WEEK;
quit;

data mylib.fm_price;
	set fm_price;
	WEEK_DATE = (WEEK - 400)*7 + 31900 - 21916;
	format WEEK_DATE YYMMDD10.;

* Plot price sequence in Chicago;	
data tmp;
	set mylib.fm_price;
	if MKT = 10; * Chicago;
run;

proc sgplot data=tmp;
	title 'Unit price series in Chicago';
	series x=WEEK_DATE y=unit_price/group=FORMAT;
run;

*** Product assortment ***;
proc sort data=weight_06unit nodupkey out=tmp; by MKT FORMAT DEPARTMENT CATEGORY;run;
proc sql noprint;
	create table tmp1 as 
	select MKT,FORMAT,CHAIN,DEPARTMENT,CATEGORY,count(distinct BRAND) as num_brand,count(distinct UPC) as num_upc
	from trip_det_unit
	group by MKT,FORMAT,CHAIN,DEPARTMENT,CATEGORY;

	create table tmp2 as
	select A.*,B.cat_weight
	from tmp1 as A left join tmp as B
	on A.MKT=B.MKT and A.FORMAT=B.FORMAT and A.CATEGORY=B.CATEGORY;

	create table chain_assort as
	select MKT, FORMAT, CHAIN, count(distinct CATEGORY) as num_cat,sum(num_brand) as allnum_brand, 
			sum(num_brand*cat_weight)/sum(cat_weight) as avgnum_brand,
			sum(num_upc) as allnum_upc, sum(num_upc*cat_weight)/sum(cat_weight) as avgnum_upc
	from tmp2
	group by MKT, FORMAT,CHAIN
	order by MKT, FORMAT;
quit;

proc means data=chain_assort noprint;
	by MKT FORMAT;
	var num_cat allnum_brand avgnum_brand allnum_upc avgnum_upc;
	output out=fm_assort MEDIAN=;
run;

*** Merge data ***;
proc sort data=fm_quant; by MKT FORMAT; run;
proc sort data=fm_assort; by MKT FORMAT; run;
proc sort data=mylib.fm_price; by MKT FORMAT; run;
proc means data=mylib.fm_price noprint;
	by MKT FORMAT;
	var unit_price;
	output out=tmp(drop=_TYPE_ _FREQ_) MEDIAN=;
run; 

data mylib.format_mkt_static;
	merge fm_quant fm_assort(drop=_TYPE_ _FREQ_) tmp;
	by MKT FORMAT;
run;

proc sql noprint;
	create table chain_cat as 
	select A.*, B.num_brand, B.num_upc
	from chain_cat_quant as A inner join tmp2 as B
	on A.MKT=B.MKT and A.FORMAT=B.FORMAT and A.CHAIN=B.CHAIN and A.CATEGORY=B.CATEGORY
	order by MKT, FORMAT, DEPARTMENT,CATEGORY;
quit;

proc means data=chain_cat noprint;
	by MKT FORMAT DEPARTMENT CATEGORY;
	var cat_price cat_unit cat_size_index cat_weight num_brand num_upc;
	output out=mylib.fmt_cat MEDIAN=;
run;




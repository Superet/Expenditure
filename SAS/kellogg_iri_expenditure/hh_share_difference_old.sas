libname mylib "E:\Users\ccv103\Documents\Research\Store switching\SAS_temp\RetailFormat";
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

data storefm;
	set storefm(drop=VARIABLE rename=(CODE=CHAIN DSC=CHAIN_NAME));
	if FORMAT in ('Warehouse stores','Supercenters','Conv. mass merchandiser','General merchandise') then FORMAT='Mass';
	if FORMAT = 'Specialty stores' then FORMAT = 'Other';
proc sort; by CHAIN;
run;

*** Read UPC ***;
PROC IMPORT OUT= WORK.dict 
            DATAFILE= "E:\Users\Projects\Project_Marketing_IRI\kelloggirir\DICT.txt" 
            DBMS=dlm REPLACE;
	 delimiter='|';
     GETNAMES=NO;
     DATAROW=1; 
RUN;

data dict;
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

data demo;
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

data democode;
	set democode (firstobs=3 rename=(DEMO_CODE=VARIABLE F2=START F3=LABEL));
	if START=. then delete;
run;

*** Format demographics ***;
data demo_fmt(keep = fmtname type start label);
	retain type 'N';
	set democode end=lastrec;
	fmtname = variable;
	output;
	/*if lastrec then do;
		hlo ='O'; 
		label = 'Missing';
		output;
	end;*/
run;

proc format cntlin=demo_fmt; run;

data demo;
	set demo;
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
	
*** Read trip code;
PROC IMPORT OUT=tripcode
            DATAFILE="E:\Users\Projects\Project_Marketing_IRI\kelloggirir\Academic Household File.xls"
            DBMS=EXCELCS REPLACE;
   SHEET='Trip Code';
RUN;

data tripcode;
	set tripcode(firstobs=3 drop=F4-F7);
	if F2=. then delete;
	rename F2=CODE F3=TRIP_DSC;
run;

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
run;

* Put together the multiple trips in the same store on the same day;
proc sql noprint;
	create table trip_sum as
	select PANID, CHAIN, DATE, sum(BSKTDOL) as BSKTDOL
	from tmp
	group by PANID,CHAIN,DATE
	order by PANID,DATE,CHAIN;
quit;

* Check any systematic errors of basket dollars of trips;
proc freq data=trip_sum;
	where BSKTDOL<0.05;
	tables bsktdol;
run;

* Drop trips with basket size equal to .01; 
data tmp; set trip_sum; if BSKTDOL>.01; run;

* Merge with chain format;
proc sql noprint;
	create table trip_sum as
	select A.*,B.FORMAT 
	from tmp as A left join storefm as B
	on A.CHAIN=B.CHAIN
	where FORMAT is not missing
	order by PANID, DATE;
quit;

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
run;

proc sql noprint;
	create table trip_det as 
	select A.*,B.FORMAT
	from tmp as A left join storefm as B
	on A.CHAIN = B.CHAIN
	order by PANID, DATE;
quit;

*****************************************************;
* Split consumers based on overal expenditure share *;
*****************************************************;
* Compute expenditure shares in each format;
proc sql noprint;
	create table tmp as
	select PANID,FORMAT, sum(BSKTDOL) as BSKTDOL
	from trip_sum
	group by PANID,FORMAT;
	
	create table tmp1 as
	select *, sum(BSKTDOL) as all_dol
	from tmp
	group by PANID;
quit;

data pan_share;
	set tmp1;
	share = BSKTDOL/all_dol;
run;

* Median split the households based on their expenditure share of wholesale clubs;
data pan_share_wholesale;
	set pan_share;
	if FORMAT = 'Wholesale clubs';
run;

* Demographic differences between low and high wholesale club shoppers;
data pan_share_wholesale;
	merge pan_share_wholesale demo;
	by PANID;
	*if share = . then share = 0;
run;

proc means data=pan_share_wholesale noprint;
	var share;
	output out=tmp1 P25=Q25 MEDIAN=Q50 P75=Q75;
run;

data tmp1;
	set tmp1;
	call symput('Q75',Q75);
run;
%put &Q75;

data pan_share_wholesale;
	set pan_share_wholesale;
	wc_share = 'No';
	if share>0 & share<&Q75 then wc_share ='Low';
	if share >= &Q75 then wc_share = 'High';
run;

ods graphics on;
proc freq data=pan_share_wholesale;
	tables wc_share*INCOME / norow chisq plots=all; 
run;

*************************************************************;
* Tendency to purchase large-size products in the recession *;
*************************************************************;
*** Product size index ***;
proc sort data=dict; by CATEGORY; run;
proc means data=dict noprint;
	by CATEGORY;
	var VOL_EQ;
	output out=tmp MEDIAN=vol_med;
run;

proc sql noprint;
	create table upc_size as
	select A.CATEGORY,A.BRAND,A.UPC,A.VOL_EQ,B.vol_med, 100 + 100*(VOL_EQ-vol_med)/vol_med as size_index
	from dict as A left join tmp as B
	on A.CATEGORY = B.CATEGORY;
quit;

*** Choose those who stayed before and during the recession ***;
data tmp;
	set trip_sum;
	recession = 0;
	if year(DATE)>=2008 then recession = 1;
run;

proc sql noprint;
	create table tmp1 as
	select PANID, sum(1-recession) as num_pre, sum(recession) as num_rec
	from tmp
	group by PANID;
	
	create table keeppan as
	select PANID
	from tmp1
	where num_pre >= 10 & num_rec>=10;
quit;

proc sql noprint;
	create table trip_sum_rec as 
	select *
	from trip_sum 
	where PANID in (select PANID from keeppan);

	create table trip_det_rec as
	select *
	from trip_det
	where PANID in (select PANID from keeppan);
quit;

*** Frequency of purchasing large-size product ***;
data trip_det_rec;
	set trip_det_rec;
	recession = 0;
	if year(DATE)>=2008 then recession = 1;
run;

%let large_quantile = 100;
proc sql noprint;
	create table trip_det_rec_size as
	select A.*,B.CATEGORY,B.size_index, ifn(B.size_index>=&large_quantile,1,0) as largesize
	from trip_det_rec as A left join upc_size as B
	on A.UPC = B.UPC;
quit;

proc sql noprint;	
	create table tmp as 
	select PANID,CATEGORY,recession, sum(largesize*DOL)/sum(DOL) as share_largesize
	from trip_det_rec_size
	group by PANID,CATEGORY,recession
	order by PANID,CATEGORY,recession;
	
	create table pan_cat_size as
	select A.*, B.INCOME, B.FAMSIZE
	from tmp as A left join demo as B
	on A.PANID = B.PANID;
quit;

proc freq data=pan_cat_size; tables INCOME; run;
proc contents data=pan_cat_size;run;
data pan_cat_size1;
	set pan_cat_size;
	if INCOME^=13;
	INCOME1 = 1;
	if INCOME>7 & INCOME<=10 then INCOME1 = 2;
	if INCOME>10 then INCOME1 = 3;
run;

proc sort data=pan_cat_size1; by descending INCOME1 PANID CATEGORY descending recession;run;
proc glm data=pan_cat_size1 order = data;
	absorb PANID CATEGORY;
	class recession INCOME1;
	model share_largesize = recession recession*INCOME1/solution;
run;

*** Households favorite size before and during the recession ***;
proc sql noprint;
	create table tmp as 
	select PANID, CATEGORY, recession, size_index, count(distinct DATE) as num_purchase
	from trip_det_rec_size 
	group by PANID, CATEGORY, recession, size_index
	order by PANID, CATEGORY, recession, num_purchase descending;
	
	create table tmp1 as 
	select *, sum(num_purchase) as all_purchase, num_purchase/sum(num_purchase) as freq
	from tmp 
	group by PANID, CATEGORY, recession 
	order by PANID, CATEGORY, recession, num_purchase descending;
quit;

data tmp1;
	set tmp1;
	fav = 0 ;
	by PANID CATEGORY recession;
	if first.recession then fav= 1;
run;

data tmp2; set tmp1; if fav = 1 & all_purchase>=2;run;

proc sql noprint;
	create table tmp as
	select *,count(recession) as num_rec
	from tmp2	
	group by PANID,CATEGORY
	having num_rec=2;

	create table mydata as
	select A.*, B.INCOME, B.FAMSIZE
	from tmp as A left join demo as B
	on A.PANID = B.PANID
	order by PANID, CATEGORY;
quit;

data mydata;
	set mydata;
	if INCOME^=13;
	INCOME1 = 1;
	if INCOME>7 & INCOME<=10 then INCOME1 = 2;
	if INCOME>10 then INCOME1 = 3;
run;

proc sort data=mydata; by descending INCOME1 PANID CATEGORY descending recession;run;
proc glm data=mydata order = data;
	absorb PANID CATEGORY;
	class recession INCOME1;
	model size_index = recession recession*INCOME1/solution;
run;

* Shift to large size;
proc sort data=mydata; by PANID INCOME1 CATEGORY; run;
proc transpose data=mydata out=mydata1 prefix = size;
	by PANID INCOME1 CATEGORY;
	id recession;
	var size_index;
run;

data mydata1;
	set mydata1;
	shift_large = 0;
	if size1 > size0 then shift_large = 1;
run;

proc sort data=mydata1; by descending INCOME1; run;
proc logistic data=mydata1 order=data;
	class INCOME1(ref='1') CATEGORY/param=ref;
	model shift_large(event = '1') = INCOME1 CATEGORY;
run;

************************;
* Expenditure volitity *;
************************;
data trip_sum_rec;
	set trip_sum_rec;
	recession = 0;
	if year(DATE)>=2008 then recession = 1;
	WEEK = cats(year(DATE), '-',put(week(DATE),2.));
run;

proc sql noprint;
	create table tmp as
	select PANID, recession, WEEK, sum(BSKTDOL) as BSKTDOL
	from trip_sum_rec
	group by PANID,recession,WEEK
	order by PANID, recession;

	create table tmp1 as
	select PANID, recession, std(BSKTDOL) as std_dol
	from tmp
	group by PANID,recession
	order by PANID, recession;
	
	create table mydata as
	select A.*,B.INCOME, B.FAMSIZE
	from tmp1 as A left join demo as B
	on A.PANID = B.PANID;
quit;

data mydata;
	set mydata;
	if INCOME^=13;
	INCOME1 = 1;
	if INCOME>7 & INCOME<=10 then INCOME1 = 2;
	if INCOME>10 then INCOME1 = 3;
run;

proc format;
	value INCOMEfmt 
		3 = 'High'
		2 = 'Medium'
		 1 = 'Low';
run;

proc sort data=mydata; by PANID  descending recession;run;
proc glm data=mydata order = data;
	format INCOME1 INCOMEfmt.;
	absorb PANID;
	class recession INCOME1;
	model std_dol = recession recession*INCOME1/solution;
run;



	

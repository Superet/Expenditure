libname myiri "E:\Users\ccv103\Documents\Research\kelloggiri\processed data";
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
proc sort; by CHAIN;
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
	create table tmp1 as
	select PANID, CHAIN, DATE, sum(BSKTDOL) as BSKTDOL
	from tmp
	group by PANID,CHAIN,DATE
	order by PANID,DATE,CHAIN;
	
	create table trip_sum as 
	select A.*,B.CHAIN_NAME,B.FORMAT
	from tmp1 as A left join storefm as B
	on A.CHAIN=B.CHAIN;
quit;

* Check any systematic errors of basket dollars of trips;
proc freq data=trip_sum;
	where BSKTDOL<0.05;
	tables bsktdol;
run;

* Drop trips with basket size equal to .01; 
data trip_sum; set trip_sum; if BSKTDOL>.01; run;

* Select households;
%let min_day = 90;
data selyear;
	input year; 
	cards;
	2006
	2007
	;
run;

proc sql noprint;
	create table tmp as 
	select PANID, max(DATE)-min(DATE) as num_day
	from trip_sum
	where year(DATE) in (select year from selyear);
	
	create table keeppan as 
	select PANID from tmp
	where num_day>=&min_day;
	
	create table trip_sum_subset as
	select *
	from trip_sum
	where (PANID in (select PANID from keeppan)) and (year(DATE) in (select year from selyear))
	order by PANID;
quit;

proc sql; select count(distinct PANID) from keeppan; quit;
proc datasets noprint; delete tmp; run;

*** Consumer loyalty profile ***;
* A macro function classify loyal consumers based on the quantile of a given variable;
%macro class_loyal(dataset,var,pctval);
	proc sort data=&dataset; by PANID &var; run;
	data &dataset;
		set &dataset;
		loyal = .;
		by PANID;
		if last.PANID then prim = 1;
	run;

	* Plot the distribution of primary share/frequency;
	data tmp;
		set &dataset;
		if prim=1;
	run;
	proc univariate data=tmp; var &var; histogram &var; run;

	* Determine the threshold value of loyalty definition;
	proc univariate data=tmp;
		var &var;
		output out=tmp1 pctlpts=&pctval pctlpre=share;
	run;
	data tmp1;
		set tmp1;
		call symput('loyal_thresh',share&pctval);
	run;
	%put &loyal_thresh;

	* Flag loyal indicator in the dataset;
	data &dataset;
		set &dataset;
		loyal = 0;
		if &var >= &loyal_thresh then loyal = 1;
	run;
%mend class_loyal;

%let min_hh = 1000;

proc sql noprint;
	create table tmp as 
	select PANID, FORMAT, CHAIN_NAME, CHAIN, sum(BSKTDOL) as BSKTDOL, sum(DATE) as num_trip
	from trip_sum_subset
	group by PANID, FORMAT, CHAIN_NAME, CHAIN
	order by PANID, BSKTDOL;
	
	create table hh_chain_dol as
	select *, sum(BSKTDOL) as DOL, BSKTDOL/sum(BSKTDOL) as share, num_trip/sum(num_trip) as freq
	from tmp
	group by PANID
	order by PANID, share;
quit;

%class_loyal(hh_chain_dol,freq,50);

proc sort data=hh_chain_dol; by FORMAT CHAIN CHAIN_NAME; run;
proc means data=hh_chain_dol noprint;
	by FORMAT CHAIN CHAIN_NAME;
	var loyal;
	output out=tmp mean=;
run;

data keep_chain;
	set tmp;
	other = substr(CHAIN_NAME,1,5);
	if _FREQ_>=&min_hh & FORMAT ^= 'Convenience' & other^='Other';
	drop other;
run;

* Plot the distribution of loyal consumers by chain;
proc univariate data=keep_chain;
	var loyal;
	histogram loyal;
run;

proc sql noprint;
	create table tmp as 
	select *
	from hh_chain_dol 
	where CHAIN in (select CHAIN from keep_chain);
	
	create table tmp1 as
	select FORMAT,CHAIN_NAME, CHAIN, loyal, sum(BSKTDOL) as DOL
	from tmp
	group by FORMAT,CHAIN_NAME, CHAIN, loyal
	order by FORMAT,CHAIN_NAME, CHAIN, loyal;
quit;

proc transpose data=tmp1 out=tmp prefix=loyal;
	by FORMAT CHAIN_NAME CHAIN;
	id loyal;
	var DOL;
run;

proc sql noprint;
	create table chain_loyal as 
	select A.*,B.loyal0,B.loyal1 , B.loyal1/(B.loyal0 + B.loyal1) as loyal_contr
	from keep_chain as A full join tmp as B
	on A.CHAIN = B.CHAIN
	order by format, chain_name;
quit;

proc univariate data=chain_loyal;
	var loyal_contr;
	histogram loyal_contr;
run;

data chain_loyal;
	set chain_loyal(drop=CHAIN _TYPE_ rename=(_FREQ_=num_hh loyal=loyal_pct));
	label	num_hh		= "Number of households observed"
			loyal_pct 	= "Percentage of loyal consumers"
			loyal0		= "Sales from unloyal consumers"
			loyal1		= "Sales from loyal consumers"
			loyal_contr	= "Contribution of sales from loyal consumers";
run;
	
proc export data=chain_loyal
            OUTFILE= "\\tsclient\Resear1\Store switching\processed data\chain_loyal.csv" 
            DBMS=csv REPLACE;
     PUTNAMES=YES;
RUN;	


* Use first 20% households in each store as loyal indicator;	
proc sort data=hh_chain_dol; by CHAIN descending share; run;
proc univariate data=hh_chain_dol; 
	by CHAIN;
	var share;
	output out=tmp pctlpts=80 pctlpre=share;
run;

data tmp1;
	merge hh_chain_dol tmp;
	by CHAIN;
	loyal = 0;
	if share>= share80 then loyal = 1;
run;

proc sql noprint;
	create table tmp as 
	select FORMAT, CHAIN_NAME, CHAIN, loyal, sum(BSKTDOL) as DOL
	from tmp1 
	where CHAIN in (select CHAIN from keep_chain)
	group by FORMAT, CHAIN_NAME, CHAIN, loyal;
quit;

proc transpose data=tmp out=tmp1 prefix=loyal;
	by FORMAT CHAIN_NAME CHAIN;
	id loyal;
	var DOL;
run;	

proc sql noprint;
	create table chain_loyal as 
	select A.*,B.loyal0,B.loyal1 , B.loyal1/(B.loyal0 + B.loyal1) as loyal_contr
	from keep_chain as A full join tmp1 as B
	on A.CHAIN = B.CHAIN
	order by format, chain_name;
quit;

proc univariate data=chain_loyal;
	var loyal_contr;
	histogram loyal_contr;
run;

data chain_loyal;
	set chain_loyal(drop=CHAIN _TYPE_ rename=(_FREQ_=num_hh loyal=loyal_pct));
	label	num_hh		= "Number of households observed"
			loyal_pct 	= "Percentage of loyal consumers"
			loyal0		= "Sales from unloyal consumers"
			loyal1		= "Sales from loyal consumers"
			loyal_contr	= "Contribution of sales from loyal consumers";
run;

	
	




	

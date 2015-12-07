%let dirname= E:\Users\ccv103\Documents\Research\Nielsen\nielsen_extracts\HMS;
options compress=yes reuse=yes; 

* Read master product data;
%let fpath = &dirname\Master_Files\Latest\products.tsv;
%put &fpath;
PROC IMPORT OUT= WORK.hms_products
            DATAFILE= "&fpath" 
            DBMS=dlm REPLACE;
			DELIMITER ='09'x;
     		DATAROW=2; 
			GUESSINGROWS=5000;
RUN;

* Read product data from RMS; 
PROC IMPORT OUT= rms_products
            DATAFILE= "E:\Users\ccv103\Documents\Research\Nielsen\nielsen_extracts\RMS\Master_Files\Latest\products.tsv" 
            DBMS=dlm REPLACE;
			DELIMITER ='09'x;
     		DATAROW=2; 
			GUESSINGROWS=5000;
RUN;

* Check if product data from HMS and RMS are identical; 
proc contents data = hms_products; 
proc contents data = rms_products; 
run;

proc compare base = hms_products compare = rms_products; 
run;

*** Read UPC from IRI data ***;
%let dirname= E:\Users\Projects\Project_Marketing_IRI\kelloggirir;
PROC IMPORT OUT= dict 
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
proc contents data=dict;run;

*** Read new product attriubtes ***;
proc import out=upc
	datafile="E:\Users\Projects\Project_Marketing_IRI\kelloggirir\cleanattr.txt"
	dbms=dlm replace;
	delimiter='|';
	getnames=no;
	guessingrows=5000;
	datarow=2;
run;

data upc;
	set upc(keep= var1-var7 var233 var264 var270 var316);
	rename var233 = PROD_TYPE
		   var264 = SEASONAL
		   var270 = SERV_SIZE
		   var316 = STORE_LOC;
	SYS = input(var1,2.);
	GEN = input(var2,2.);
	VEN = input(var5,5.);   
	ITE = input(var6,5.);
	WFM = input(var3,6.);
	WLM = input(var4,6.);
	KEYCAT = input(var7,6.);
	format SYS z2. 
		   GEN z2. 
		   VEN z5.
		   ITE z5.;
	UPC = input( cats(put(sys,z2.),put(gen,z2.),put(ven,z5.),put(ite,z5.)),14.);
	upc12 = input( cats(put(sys,z2.),put(ven,z5.),put(ite,z5.)),12.);
	drop var1-var7;
	DFM = (WFM - 400)*7 + 31900 - 21916;
	DLM = (WLM - 400)*7 + 31900 - 21916;
	format DFM YYMMDD10. DLM YYMMDD10.;
	LIFE = WLM - WFM;
	Year = year(DFM);
proc sort; by UPC; 
run;
proc contents data=upc; run;

proc sql noprint;
	create table iri_products as 
	select A.*, B.upc12, WFM, WLM, DFM, DLM
	from dict as A inner join upc as B
	on A.upc = B.upc;
quit; 

* Match upc from two datasets; 
proc sql noprint;
	create table check_match as
	select A.upc, A.upc12 , B.upc as nls_upc, A.CATEGORY, B.product_module_descr, A.Brand, B.brand_descr, A.upcdsc,B.upc_descr
	from iri_products as A inner join (select * from rms_products where upc_ver_uc = 1) as B
	on A.upc12 = B.upc
	order by category;
quit; 

* Check the upc in IRI that do not have matches in Nielsen; 
proc sql noprint;
	create table check_unmatch as 
	select *
	from iri_products
	where upc not in (select upc from check_match); 
quit; 

data check_unmatch; 
	set check_unmatch; 
	year = year(DFM); 
proc freq; table year; 
run;

* Check how many new products from harbinger project do not have matches; 
libname prod 	"E:\Users\ccv103\Documents\Research\Store switching\SAS_temp\Prod";
proc sql noprint;
/*	95744 new UPC are matched*/
	create table newprod_match as 
	select *
	from PROD.new_full 
	where upc in (select upc from check_match);
	
/*	17285 new UPC can not be matched*/
	create table newprod_unmatch as 
	select *
	from PROD.new_full 
	where upc not in (select upc from check_match);
quit; 



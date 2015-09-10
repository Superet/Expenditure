LIBNAME mysqllib OLEDB
OLEDB_SERVICES=NO
Datasource="kdc01\kdwh02"
PROVIDER=SQLOLEDB.1
Properties=('Initial Catalog'=USRDB_ccv103
                  'Integrated Security'=SSPI)
SCHEMA=DBO
BULKLOAD=YES
bl_options='ROWS_PER_BATCH = 500000';

libname mylib "E:\Users\ccv103\Documents\Research\Store switching\SAS_temp\Nielsen";
libname mysel "E:\Users\ccv103\Documents\Research\Store switching\SAS_temp\Nielsen\Work";
libname mainex "E:\Users\ccv103\Documents\Research\Store switching\SAS_temp\Nielsen\Work\main";

%let dirname= E:\Users\ccv103\Documents\Research\Nielsen\nielsen_extracts\HMS;
options compress=yes reuse=yes;

*************;
* Read data *;
*************;
* Read in codebook;
PROC IMPORT OUT= mysel.my_codebook
            DATAFILE= "\\tsclient\Resear1\Data\Nielsen\reference_documentation\code_book.csv" 
            DBMS=csv REPLACE;
     		GETNAMES=YES; 
RUN;

* Read in the data of deleted UPC;
PROC IMPORT OUT= my_delete_upc
            DATAFILE= "\\tsclient\Resear1\Store switching\processed data\deleted_upc.csv" 
            DBMS=csv REPLACE;
     		GETNAMES=YES; 
RUN;

* Read in CPI data; 
PROC IMPORT OUT= cpi
            DATAFILE= "\\tsclient\Resear1\Store switching\processed data\CPI.csv" 
            DBMS=csv REPLACE;
     		GETNAMES=YES; 
RUN;

* Classify product department by food and non-edible;
data products;
	set mylib.products;
	* Classify food and non-edible departments; 
	Dpmt_type = "Non-edible";
	if my_magnet = 1 then Dpmt_type = "Magnet";
	else if	department_descr in ("ALCOHOLIC BEVERAGES","DAIRY","DELI","DRY GROCERY","FRESH PRODUCE","FROZEN FOODS","PACKAGED MEAT") 
	then do;
		Dpmt_type = "Food";
		Food_module = product_module_descr;
	end;
	else Nonedible_module = product_module_descr;
	* Classify storable and perishable departments; 
	Storable = 1;
	if department_descr in ("DAIRY","DELI","FRESH PRODUCE") then Storable = 0;
	if department_descr = "DRY GROCERY" and substr(product_module_descr,1,6) = "BAKERY" then Storable = 0; 
	if my_magnet = 1 then Storable = 0; 
	if Storable = 1 then storable_module = product_module_descr;
	if Storable = 0 then nonstorable_module = product_module_descr;
	if my_magnet = 1 then magnet_module = product_module_descr;
run;

* NOTE: deleted UPCs are already dropped in mylib.products; 
* Read master product data;
%let fpath = &dirname\Master_Files\Latest\products.tsv;
%put &fpath;
PROC IMPORT OUT= WORK.orig_products
            DATAFILE= "&fpath" 
            DBMS=dlm REPLACE;
			DELIMITER ='09'x;
     		DATAROW=2; 
			GUESSINGROWS=5000;
RUN;

* Find the UPC for deleted product modules;
proc sql noprint;
	create table tmp1 as 
	select  A.product_module_descr, B.upc, B.upc_ver_uc, cats(B.upc, "-", B.upc_ver_uc) as upcv, B.upc_descr
	from (select * from my_delete_upc where upc_descr^="") as A inner join orig_products as B
	on A.upc_descr=B.upc_descr;
	
	create table tmp2 as 
	select A.product_module_descr, B.upc, B.upc_ver_uc, cats(B.upc, "-", B.upc_ver_uc) as upcv, B.upc_descr
	from (select * from my_delete_upc where upc_descr="") as A inner join orig_products as B
	on A.product_module_descr=B.product_module_descr;
quit;

data mysel.my_delete_upc; set tmp1 tmp2; run;
proc datasets noprint; delete orig_products; run;

* Define the retailer scope; 
data mysel.my_channel;
	input channel_type $40.	channel_group $10.;
	datalines;
	Grocery								Grocery
	Drug Store							Grocery
	Dollar Store						Grocery
	Convenience Store					Grocery
	Discount Store						Grocery
	Warehouse Club						Grocery
	;
run;

* Read master retailer data;
%let fpath = &dirname\Master_Files\Latest\retailers.tsv;
%put &fpath;
PROC IMPORT OUT= WORK.orig_retailers
            DATAFILE= "&fpath" 
            DBMS=dlm REPLACE;
			DELIMITER ='09'x;
     		DATAROW=2; 
			GUESSINGROWS=5000;
RUN;

* NOTE: we extract purchase data from only grocery stores;  
proc sql noprint;
	create table mysel.retailers as
	select A.*, B.channel_group
	from orig_retailers as A inner join mysel.my_channel as B
	on A.channel_type = B.channel_type;
quit;

data mysqllib.tmp_products (bulkload = YES); 
	set products; 
run;

* Summarize a shopping basket in bi-weekly basis; 
proc sql;
	connect to oledb as mydb
	(OLEDB_SERVICES=NO Datasource="kdc01\kdwh02" PROVIDER=SQLOLEDB.1 
	Properties=('Initial Catalog'=USRDB_ccv103 'Integrated Security'=SSPI)
	SCHEMA=DBO);	
	
	EXECUTE(		
		/*Append product attributes to purchase data; */
		select C.*, D.size_index, D.size, D.Dpmt_type, D.Food_module, D.Nonedible_module, 
				D.storable, D.storable_module, D.nonstorable_module, D.product_module_descr, D.my_magnet, D.weight
		into dbo.tmp 		
		from dbo.purchases as C inner join 
		(select A.*, B.weight from dbo.tmp_products as A inner join dbo.basket as B on A.product_module_descr = B.product_module_descr) as D
		on C.upcv = D.upcv;
		
		/*Summarize basket size over two weeks; */
		select household_code, min(year) as year, biweek, sum(total_price_paid - coupon_value) as dol_purchases,
				count(distinct purchase_date) as num_day, count(distinct trip_code_uc) as num_trip, 
				sum((total_price_paid - coupon_value)*my_magnet) as dol_magnet, 
				count(distinct Food_module) as num_food_module, count(distinct Nonedible_module) as num_noned_module,
				count(distinct storable_module) as num_storable_module, count(distinct nonstorable_module) as num_nonstr_module,
				sum(CASE WHEN Dpmt_type = 'Food' THEN size_index*quantity*weight ELSE 0 END) as food_quant, 
				sum(CASE WHEN Dpmt_type = 'Non-edible' THEN size_index*quantity*weight ELSE 0 END) as nonedible_quant, 
				sum(storable*size_index*quantity*weight) as storable_quant, 
				sum((1-storable)*size_index*quantity*weight) as nonstr_quant
		into dbo.biweek_basket 
		from (select A.*, B.household_code, B.purchase_date from 
			  dbo.tmp as A inner join dbo.trips as B on A.trip_code_uc = B.trip_code_uc) as t
		group by household_code, biweek
		order by household_code, biweek;
		
		drop table dbo.tmp; 
	) by mydb; 
	disconnect from mydb;
quit;

proc datasets noprint; delete tmp tmp1 tmp2; run;

************************;
* Household selection * ;
************************;
/*Selection critetion:
1. The maximum inter-purchase gap is less than 30 days; 
2. Households have at least 6 observations in two continual years; */

%let gap_thresh = 30;
%let min_n		= 6;

* Compute the inter-purchase date;
proc sql; 
	connect to oledb as mydb
	(OLEDB_SERVICES=NO Datasource="kdc01\kdwh02" PROVIDER=SQLOLEDB.1 
	Properties=('Initial Catalog'=USRDB_ccv103 'Integrated Security'=SSPI)
	SCHEMA=DBO);
	
	create table tmp as 
	select * from connection to mydb (
		select *, recession = CASE WHEN year >= 2008 THEN 1 ELSE 0 END,
			DATEDIFF(day, LAG(purchase_date, 1) OVER (PARTITION BY household_code ORDER BY purchase_date), purchase_date) as inter_date
		from dbo.trips;
	);
	disconnect from mydb;
quit;

* Find the maximum interpurchase date for each household and plot the histogram;
* NOTE: households may move to a different geo market; 
proc sql noprint;
	create table tmp1 as 
	select scantrack_market_descr, household_code, max(inter_date) as gap, 
		sum(biweek*(1-recession) > 0) as n1, sum(biweek*recession > 0) as n2, 
		sum(biweek*(year=2004)) as nyr_2004, sum(biweek*(year=2005)) as nyr_2005, sum(biweek*(year=2006)) as nyr_2006, 
		sum(biweek*(year=2007)) as nyr_2007, sum(biweek*(year=2008)) as nyr_2008, sum(biweek*(year=2009)) as nyr_2009, 
		sum(biweek*(year=2010)) as nyr_2010, sum(biweek*(year=2011)) as nyr_2011, sum(biweek*(year=2012)) as nyr_2012
	from tmp 
	group by scantrack_market_descr, household_code;
	
	create table tmp2 as 
	select household_code, max(gap) as gap, count(distinct scantrack_market_descr) as n_mkt, sum(n1) as n1, sum(n2) as n2, 
		sum(nyr_2004) as nyr_2004, sum(nyr_2005) as nyr_2005, sum(nyr_2006) as nyr_2006, sum(nyr_2007) as nyr_2007, 
		sum(nyr_2008) as nyr_2008, sum(nyr_2009) as nyr_2009, sum(nyr_2010) as nyr_2010, sum(nyr_2011) as nyr_2011, 
		sum(nyr_2012) as nyr_2012
	from tmp1
	group by household_code;
quit;

* Compute the number of years of active stay; 
data tmp2; 
	set tmp2; 
	n_all = 0; 
	array ny{*} nyr:;
	do i = 1 to dim(ny);
		if ny[i] > &min_n then n_all + 1; 
	end;
run; 

* Check how many households will be dropped; 
proc univariate data=tmp2 noprint; histogram gap; run;
proc sql; 
	select count(gap) from tmp2 where gap <= &gap_thresh; 
	select sum(n_mkt>1) from tmp2;
	select sum(n1 > &min_n and n2 > &min_n) from tmp2;
	select sum(n_all >=2) from tmp2; 
quit;

* Drop the households who do not satisfy our selection criterion; 
proc sql noprint;
	create table mysel.keeppan as 
	select household_code 
	from tmp2 
	where gap <= &gap_thresh and n_mkt=1 and (n_all >= 2)
	order by household_code;
quit;

data mysqllib.keeppan; set mysel.keeppan; run;

*******************************************;
* Construct stone price index every biweek *;
*******************************************;
* COMPUTE STONE PRICE INDEX FROM RETAIL FORMAT-LEVEL PRICE INDEX;
proc sql;
	connect to oledb as mydb
	(OLEDB_SERVICES=NO Datasource="kdc01\kdwh02" PROVIDER=SQLOLEDB.1 
	Properties=('Initial Catalog'=USRDB_ccv103 'Integrated Security'=SSPI)
	SCHEMA=DBO);
	
	* Expenditure share of retail formats from all households in a scantrack_market_descr each biweek;
	create table tmp as 
	select * from connection to mydb (
		select scantrack_market_descr, biweek, channel_type, min(year) as year, 
				sum(total_spent) as revenue 
		from dbo.trips
		group by  scantrack_market_descr, biweek, channel_type
	);	
	
	* Compute stone price index;
	create table mysel.biweek_stone_price as 
	select  scantrack_market_descr, min(year) as year, biweek, 
			sum(revenue*log(bsk_price_paid))/sum(revenue) as stone_price
	from (select A.*, B.revenue from mylib.format_biweek_price as A inner join tmp as B 
		  on A.scantrack_market_descr=B.scantrack_market_descr 
			 and A.biweek=B.biweek and A.channel_type=B.channel_type)
	group by  scantrack_market_descr, biweek
	order by  scantrack_market_descr, biweek;
	
	disconnect from mydb;
quit;

proc datasets noprint; delete tmp tmp1 tmp2; run;

**************************************************************;
* Construct the data for the process of biweekly expenditure *;
**************************************************************;
*---------------------------*;
* Fill in nonshopping weeks *;
proc sql noprint;
	* Aggregate daily trip data to biweekly expenditure data; 
	create table tmp as 
	select household_code, min(year) as year, biweek, sum(net_spent) as dol, sum(dolp) as dolp, sum(dol_nonscan) as dol_nonscan
	from (select * from mysqllib.trips where household_code in (select household_code from mysel.keeppan))
	group by household_code, biweek;
	
	create table tmp1 as 
	select household_code, min(biweek) as min_n, max(biweek) as max_n
	from tmp
	group by household_code
	order by household_code; 
	
	create table biweek_year as 
	select biweek, min(year) as year
	from tmp
	group by biweek
	order by biweek; 
quit;

data tmp1(drop=min_n max_n);
	set tmp1;
	by household_code;
	biweek = min_n - 1;
	if first.household_code then do; 	
		do while(biweek<max_n);
			biweek + 1; output;
		end;
	end;
run;

* Merge together expenditure with complete biweek index; 	
proc sql noprint;
	create table tmp_exp as 
	select B.*, A.dol, A.dolp, A.dol_nonscan
	from tmp as A right join tmp1 as B
	on A.household_code=B.household_code and A.biweek=B.biweek
	order by biweek;
quit;

* Add year vairable; 
data tmp_exp; 
	merge tmp_exp biweek_year;
	by biweek;
proc sort; by household_code biweek; 
run; 

* Check the non-shopping weeks; 
proc sql; select nmiss(dol)/count(biweek) from tmp2; quit;
proc sql noprint;
	create table tmp as 
	select household_code, nmiss(dol)/count(biweek) as n_missbiweek
	from tmp_exp
	group by household_code
	order by household_code;
quit;
proc univariate data = tmp; histogram n_missbiweek; run;

proc datasets noprint; delete tmp tmp1 tmp2; run;

*---------------------------------------------------------------------*;
* Compute household expenditure share across retail format each biweek ;
proc sql noprint;
	connect to oledb as mydb
	(OLEDB_SERVICES=NO Datasource="kdc01\kdwh02" PROVIDER=SQLOLEDB.1 
	Properties=('Initial Catalog'=USRDB_ccv103 'Integrated Security'=SSPI)
	SCHEMA=DBO);
	
	create table tmp as
	select * from connection to mydb (
		select scantrack_market_descr, household_code, min(year) as year, biweek, channel_type, sum(net_spent) as DOL, 
				sum(dolp) as dolpurch
		from (select * from dbo.trips where household_code in (select household_code from dbo.keeppan)) as t
		group by scantrack_market_descr, household_code, biweek, channel_type
		order by scantrack_market_descr, household_code, biweek, channel_type;
	); 
	disconnect from mydb; 
quit;

proc transpose data = tmp out = tmp_share prefix=DOL_;
	by scantrack_market_descr household_code biweek;
	id channel_type;
	var DOL;
run;

data tmp_share; 
	set tmp_share(drop=_NAME_);
	array dolvar{*} DOL:;
	do i = 1 to dim(dolvar);
		if dolvar[i]=. then dolvar[i] = 0;
	end;
	drop i;
run;

*---------------------------------------*;
* Constrct data of household demographics;
data panelists; set mainex.panelists; run;
proc sort data=panelists; by household_code panel_year; run;

* Adjust the income reporting; 
data tmp; 
	set panelists(keep = household_code panel_year household_income);
	income_year = panel_year - 2;
	drop panel_year; 
run;

* NOTE: we lose data from the last two years because of the lag income reporting;
proc sql noprint;
	create table tmp1 as 
	select A.*, B.household_income as income_real
	from panelists as A inner join tmp as B
	on A.household_code = B.household_code and A.panel_year = B.income_year
	order by household_code, panel_year;
quit;

* The initial income status of households; 
data tmp2(keep = household_code first_income); 
	set tmp1;
	by household_code panel_year; 
	first_income = income_real; 
	if not first.household_code then delete;
run; 

data tmp;
	merge tmp1(keep = household_code panel_year scantrack_market_descr 
			income_real household_size type_of_residence male_head_employment female_head_employment age_and_presence_of_children)
		  tmp2; 
	by household_code; 		
	famsize = "Three+";
	if household_size = 2 then famsize = "Two";
	if household_size = 1 then famsize = "Single";
	condo = 0; 
	if mod(type_of_residence, 2)=0 then condo = 1;
	employed = 0;
	if (male_head_employment <= 3 | female_head_employment <=3) then employed = 1;
	NumChild = 0;
	if age_and_presence_of_children <= 3 then NumChild = 1;
	if age_and_presence_of_children >3 and age_and_presence_of_children <= 7 then NumChild = 2;
run;

* Add the middle value of income range; 
proc sql noprint;
	create table panelists as
	select A.*, B.mid_range as income_midvalue
	from tmp as A left join (select code_value, description, mid_range from mysel.my_codebook where variable='household_income') as B
	on A.income_real = B.code_value
	order by household_code;	
quit;

proc datasets noprint; delete tmp tmp1 tmp2; run; 

*-------------------------------------------------------------------------*;
* Merge household demongraphics and price index to biweekly expenditure data; 
proc sql noprint;
	* Combine expenditure share, and basket data together; 
	create table tmp1(drop = old_household_code old_biweek dol_purchases year) as 
	select * 
	from tmp_share as A left join mysqllib.biweek_basket(rename = (household_code = old_household_code biweek = old_biweek)) as B 
	on A.household_code = B.old_household_code and A.biweek = B.old_biweek; 

	* Combine complete shopping expenditure and expenditure share; 
	create table tmp2(drop = old_household_code old_biweek)  as 
	select *
	from tmp_exp as A left join tmp1(rename = (household_code = old_household_code biweek = old_biweek)) as B
	on A.household_code = B.old_household_code and A.biweek = B.old_biweek;
	
	* Merge in household demographics; 
	create table tmp3(drop = household_code_old panel_year_old) as 
	select *
	from tmp2(drop = scantrack_market_descr) as A inner join 
	panelists(rename=(household_code = household_code_old panel_year=panel_year_old)) as B
	on A.household_code = B.household_code_old and A.year = B.panel_year_old
	order by household_code, biweek;
	
	* Combine price index; 
	create table hh_biweek_exp_merge as 
	select A.*, B.stone_price
	from tmp3 as A left join mysel.biweek_stone_price as B
	on A.scantrack_market_descr = B.scantrack_market_descr and A.biweek = B.biweek;
quit; 

proc datasets noprint; delete tmp1 tmp2 tmp3; run;

*---------------------------------*;
* Adjust all dollar measure for CPI; 
proc sql noprint;
	select Annual into: price_base from cpi where year = 2004; 
	
	create table tmp as 
	select A.*, B.Annual/&price_base as cpi
	from hh_biweek_exp_merge as A left join cpi as B
	on A.year = B.year
	order by household_code, year, biweek;
quit;

* Extract size index and weight of non-scannables; 
proc sql noprint; 
	connect to oledb as mydb
	(OLEDB_SERVICES=NO Datasource="kdc01\kdwh02" PROVIDER=SQLOLEDB.1 
	Properties=('Initial Catalog'=USRDB_ccv103 'Integrated Security'=SSPI)
	SCHEMA=DBO);
	
	select weight into :nonscan_weight from mylib.basket where product_module_descr = "Nonscan"; 
	select * into :nonscan_size from connection to mydb (
		select avg(dol_nonscan)
		from (select household_code, biweek, sum(dol_nonscan) as dol_nonscan
				from dbo.trips group by household_code, biweek) as t; 
	); 
	disconnect from mydb; 
quit; 
%put &nonscan_size; 
%put &nonscan_weight; 

* Adjust CPI, set unmatched values to zero; 
* Also, compute overall basket quantity; 
data hh_biweek_exp_merge (drop=i j _LABEL_ age_and_presence_of_children condo female_head_employment male_head_employment type_of_residence);
	set tmp;
	array dolvar{*} dol: stone_price;
	do i=1 to dim(dolvar);
		dolvar[i] = dolvar[i]/cpi;
	end;
	array valvar{*} dol: num: food_quant nonedible_quant storable_quant nontr_quant stone_price; 
	do j = 1 to dim(valvar); 
		if valvar[j] = . then valvar[j] = 0; 
	end;
	recession = 0; 
	if year >= 2008 then recession = 1; 
	Q = food_quant + nonedible_quant + dol_nonscan/&nonscan_size* &nonscan_weight; 
run;	

* Check if the data has duplicated by key;
* NOTE: the data is not unique in keys because some households move; 
proc sort data=hh_biweek_exp_merge nodupkeys out=tmp; 
by household_code biweek;
run;

proc contents data = hh_biweek_exp_merge; run;

* Save data; 
data mylib.panelists; set panelists; run;
data mylib.hh_biweek_exp_merge; set hh_biweek_exp_merge; run;

proc datasets noprint; delete tmp tmp_share tmp_exp; run;

****************************;
* Export the data to excel *;
****************************;
PROC EXPORT DATA=mylib.panelists
			OUTFILE= "\\tsclient\Resear1\Store switching\processed data\2_panelists.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;

PROC EXPORT DATA=mylib.hh_biweek_exp_merge
			OUTFILE= "\\tsclient\Resear1\Store switching\processed data\2_hh_biweek_exp_merge.csv" 
            DBMS=CSV REPLACE;
     PUTNAMES=YES;
RUN;

/*
proc datasets library=mysel kill noprint;
proc datasets library=mainex kill noprint;
run;
*/


/*This scipt is to stack all the UPC weight files together. 
Each file is the volumn sold in 2006 for a product group. 
*/
libname mymaster "/sscc/home/c/ccv103/rms_attributes/master";
libname	grpyear "/sscc/home/c/ccv103/rms_attributes/group_year";
libname	newprod "/sscc/home/c/ccv103/rms_attributes/newprod";
%let dirname= /sscc/home/c/ccv103/rms_attributes;

* Read file names; 
%let subdir	= &dirname/master/upc_weights;
%put &subdir; 
filename dirlist pipe "ls &subdir";
data dirlist;
	length fname $30 filep $300.;
	infile dirlist length=reclen;
	input fname $varying30. reclen;
	filep = cats("&subdir", "/",fname); 
run;	

* Stack files; 
data upc_weight;
	set dirlist(keep = filep);
	infile dummy filevar = filep delimiter='09'x MISSOVER DSD LRECL=32767 FIRSTOBS=2 end = done;
	do while(not done);
		input store_code_uc upc upc_ver_uc product_module_code vol_weight;
		output; 
	end;
	drop filep;
run;
proc contents data = upc_weight; run;

data mymaster.upc_weight; 
	set upc_weight; 
run;



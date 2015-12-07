libname mymaster "/sscc/home/c/ccv103/rms_attributes/master";
libname	grpyear "/sscc/home/c/ccv103/rms_attributes/group_year";
libname	newprod "/sscc/home/c/ccv103/rms_attributes/newprod";
%let dirname= /sscc/home/c/ccv103/rms_attributes;

%macro assign_value(single); 
	%global sel_group sel_year; 
	%if (&single = 1) %then %do; 
		%let sel_group 	= %sysget(sel_group); 
		%let sel_year	= %sysget(sel_year); 
	%end; 
	%else %do; 
		%let grp_idx =  %sysget(PBS_ARRAY_INDEX); 
		%let sel_year	= %sysget(sel_year); 
		proc sql noprint;
			select product_group_code into: sel_group from mymaster.prod_group_index
			where n = &grp_idx; 
		quit; 				
	%end; 
%mend assing_value; 

%assign_value(%sysget(single));
%put &sel_group; 
%put &sel_year;
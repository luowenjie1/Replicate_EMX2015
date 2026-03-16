
****This STATA do file creates the sectoral moments in Table 1***

clear
clear matrix

set mem 800m
*************1. Plant-product Level Data*****************
cd "C:\Users\yx42\Documents\Research Projects\Competition and Markup\Revision\data_moments"

use "C:\Users\yx42\Documents\Research Projects\Competition and Markup\Revision\raw data\prod93.dta"
append using "C:\Users\yx42\Documents\Research Projects\Competition and Markup\Revision\raw data\prod92.dta"
append using "C:\Users\yx42\Documents\Research Projects\Competition and Markup\Revision\raw data\prod91.dta"
append using "C:\Users\yx42\Documents\Research Projects\Competition and Markup\Revision\raw data\prod89.dta"

***drop negative or missing values
drop if isalamt<=0
drop if psales<=0
drop if stno==.

**last two-digit indicating "Other"
g last2=prodno-floor(prodno/100)*100
g dlast=last2>=90

unique stno

table dlast year, c(sum psales)

egen last=sum(dlast), by(year stno)   /***last indicate A PLANT which has a product that's dropped due to classification 90-99***/
drop if last2>=90
drop if prodno<=1000000 /****agriculture and mining****/

g dsales=isalamt
g esales=esalamt

drop last2 dlast

sort year stno

collapse (sum) dsales esales psales last, by (stno year prodno) 
save product.dta, replace

sort prodno

*************************start to generate sectoral_moments*********************
#delimit;

sort stno year prodno;
egen nprod=count(prodno), by(year stno);
collapse (sum) dsales esales psales (mean) nprod last , by(year stno prodno);
egen prodplant=group(prodno stno);
xtset prodplant year;

egen indsales=sum(dsales), by(year prodno);
egen countcomp=count(stno), by(year prodno);

g share = dsales/indsales;

sort year prodno stno;
g share2=share^2;
egen HH=sum(share2), by(prodno year);
g iHH=HH^(-1);
egen smax = max(share), by(prodno year);

merge m:1 stno year using plant.dta;
tab _merge;
drop if _merge==2;

/*preserve;
keep if _merge==3;
save plant_product.dta, replace;
restore;*/

/***split wage bill based on sales share within firm***/
egen fsales=sum(psales), by (stno year);
gen pshare=dsales/fsales;
sum pshare;
gen pvl=pshare*vl;

/**Unconditional Share**/
sum share, d; 

/* concentration statistics excl. imports */
collapse (mean) countcomp indsales HH iHH smax (sum) pvl eindsales=esales (p90) share90=share (p10) share10=share (p75) share75=share (p25) share25=share  (sd) sdshare=share, by(prodno year);
replace sdshare=0 if countcomp==1;
gen dshare9010=share90-share10;
gen dshare7525=share75-share25;

**within sector moments;
tabstat countcomp iHH smax dshare9010 dshare7525 sdshare, stat(mean p50 sd count );
**cross sector moments directly generated in MATLAB

sort prodno year;
tab year;
gen gindsales=indsales;

**pool all years together;
gsort -gindsales;
gen ranks=[_n];
egen ts=sum(gindsales);
egen tl=sum(pvl);
egen tobs=count(year);

gen ranks1=1 if ranks/tobs<=0.01;
gen ranks5=1 if ranks/tobs<=0.05;

gen s1=gindsales*ranks1/ts;
gen s5=gindsales*ranks5/ts;

gen lv1=pvl*ranks1/tl;
gen lv5=pvl*ranks5/tl;

collapse (sum)  s1 s5 lv1 lv5;
tabstat s1 s5 lv1 lv5;






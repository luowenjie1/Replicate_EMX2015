clear

clear matrix

set mem 800m
*************1. Firm-product Level Data*****************
cd "C:\Users\yx42\Documents\Research Projects\Competition and Markup\Revision\data_moments\firm moments"

use "C:\Users\yx42\Documents\Research Projects\Competition and Markup\Revision\raw data\prod93.dta"
append using "C:\Users\yx42\Documents\Research Projects\Competition and Markup\Revision\raw data\prod92.dta"
append using "C:\Users\yx42\Documents\Research Projects\Competition and Markup\Revision\raw data\prod91.dta"
append using "C:\Users\yx42\Documents\Research Projects\Competition and Markup\Revision\raw data\prod89.dta"

**replicate cleaning we were doing before
drop if isalamt<=0
drop if psales<=0

**last two-digit in 90
g last2=prodno-floor(prodno/100)*100
g dlast=last2>=90

unique stno

sort year stno
merge m:1 year stno using  stno-bran.dta
keep if _merge==3

table dlast year, c(sum psales)

egen last=sum(dlast), by(year bran)   
drop if last2>=90
g dsales=isalamt
g esales=esalamt

drop last2 dlast

sort year bran
drop if bran==.
collapse (sum) dsales esales psales last, by (bran year prodno) 

sort prodno

rename bran stno
save product.dta, replace

**************2. Creat firm level data***************
use  plant.dta
sort year stno
merge m:1 year stno using  stno-bran.dta
keep if _merge==3
drop _merge
collapse (sum) tbusincome sales vm vl ql vk ve vs voth (mean) sic2, by (bran year) 
rename bran stno
save firm.dta, replace

use product.dta, clear

*************3. Rest of the code same as plant-level moments
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
sum share, d; 

merge m:1 stno year using  firm.dta;
tab _merge;
drop if _merge==2;

preserve;
keep if _merge==3;
save firm_product.dta, replace;
restore;

***split wage bill based on sales share within firm;
egen fsales=sum(psales), by (stno year);
gen pshare=dsales/fsales;
sum pshare;

gen pvl=pshare*vl;

/* concentration statistics excl. imports */
collapse (mean) countcomp indsales HH iHH smax (sum) pvl eindsales=esales (p90) share90=share (p10) share10=share (p75) share75=share (p25) share25=share  (sd) sdshare=share, by(prodno year);
replace sdshare=0 if countcomp==1;
gen dshare9010=share90-share10;
gen dshare7525=share75-share25;

**within sector moments;
tabstat countcomp iHH smax dshare9010 dshare7525 sdshare, stat(mean p50 sd count );
**cross sector moments;
tabstat iHH smax countcomp, stat(p10 p25 p50 p75 p90 p95);

**pool all years together;
gsort -indsales;
gen ranks=[_n];
egen ts=sum(indsales);
egen tl=sum(pvl);
egen tobs=count(year);

gen ranks1=1 if ranks/tobs<=0.01;
gen ranks5=1 if ranks/tobs<=0.05;

gen s1=indsales*ranks1/ts;
gen s5=indsales*ranks5/ts;

gen lv1=pvl*ranks1/tl;
gen lv5=pvl*ranks5/tl;


collapse (sum)  s1 s5 lv1 lv5;



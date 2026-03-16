/*
-------------------------------------------------------------------------------
sectoral_moments.do
Purpose: Construct sector-level concentration/distribution moments (Table 1).

[DATA AVAILABILITY NOTICE]
This script uses confidential raw files under `../../raw data/`:
- prod89.dta, prod91.dta, prod92.dta, prod93.dta
and an internal file `plant.dta` in `../../data_moments/`.
Without these files, the script cannot be executed end-to-end.
-------------------------------------------------------------------------------
*/

clear
clear matrix
set mem 800m

* 1) Build pooled plant-product panel (1989/1991/1992/1993).
cd "../../data_moments"

use "../../raw data/prod93.dta"
append using "../../raw data/prod92.dta"
append using "../../raw data/prod91.dta"
append using "../../raw data/prod89.dta"

* Remove invalid sales and missing plant IDs.
drop if isalamt<=0
drop if psales<=0
drop if stno==.

* Drop "other" product codes ending in 90-99 and non-manufacturing codes.
g last2 = prodno-floor(prodno/100)*100
g dlast = last2>=90

unique stno
table dlast year, c(sum psales)

egen last = sum(dlast), by(year stno)  /* plant ever has product in 90-99 */
drop if last2>=90
drop if prodno<=1000000               /* agriculture and mining */

* Rename domestic/export sales variables for later consistency.
g dsales = isalamt
g esales = esalamt
drop last2 dlast

* Save cleaned product panel.
sort year stno
collapse (sum) dsales esales psales last, by(stno year prodno)
save product.dta, replace
sort prodno

* 2) Compute sector-year concentration moments.
#delimit ;

sort stno year prodno;
egen nprod = count(prodno), by(year stno);
collapse (sum) dsales esales psales (mean) nprod last, by(year stno prodno);
egen prodplant = group(prodno stno);
xtset prodplant year;

egen indsales = sum(dsales), by(year prodno);
egen countcomp = count(stno), by(year prodno);

g share = dsales/indsales;

* HHI and top firm share by product-year.
sort year prodno stno;
g share2 = share^2;
egen HH = sum(share2), by(prodno year);
g iHH = HH^(-1);
egen smax = max(share), by(prodno year);

* Bring in plant-level covariates and keep matched sample.
merge m:1 stno year using plant.dta;
tab _merge;
drop if _merge==2;

/*
(Optional in original code) Save matched product-plant file for later scripts.
preserve;
keep if _merge==3;
save plant_product.dta, replace;
restore;
*/

* Split plant wage bill across products using product sales weights.
egen fsales = sum(psales), by(stno year);
gen pshare = dsales/fsales;
sum pshare;
gen pvl = pshare*vl;

* Unconditional distribution of firm market shares.
sum share, d;

/*
Collapse to product-year moments:
- countcomp: number of producers
- iHH: inverse HHI
- smax: largest firm share
- dshare9010 / dshare7525 / sdshare: dispersion of shares
- pvl / eindsales aggregates used for top-sector calculations
*/
collapse (mean) countcomp indsales HH iHH smax (sum) pvl eindsales = esales (p90) share90 = share (p10) share10 = share (p75) share75 = share (p25) share25 = share (sd) sdshare = share, by(prodno year);

replace sdshare = 0 if countcomp==1;
gen dshare9010 = share90-share10;
gen dshare7525 = share75-share25;

* Within-sector summary moments (cross-sector moments are computed in MATLAB).
tabstat countcomp iHH smax dshare9010 dshare7525 sdshare, stat(mean p50 sd count);

* 3) Pool all years to compute top 1% / 5% sector shares.
sort prodno year;
tab year;
gen gindsales = indsales;

gsort -gindsales;
gen ranks = [_n];
egen ts = sum(gindsales);
egen tl = sum(pvl);
egen tobs = count(year);

gen ranks1 = 1 if ranks/tobs<=0.01;
gen ranks5 = 1 if ranks/tobs<=0.05;

gen s1 = gindsales*ranks1/ts;
gen s5 = gindsales*ranks5/ts;

gen lv1 = pvl*ranks1/tl;
gen lv5 = pvl*ranks5/tl;

collapse (sum) s1 s5 lv1 lv5;
tabstat s1 s5 lv1 lv5;

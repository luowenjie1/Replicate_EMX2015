/*
-------------------------------------------------------------------------------
plant_moments.do
Purpose: Reproduce plant-level concentration/labor moments used in Table 1.

[DATA AVAILABILITY NOTICE]
This script depends on confidential micro data file `plant_product.dta`
located under `../../data_moments/`. If this file is missing in your copy,
you can still read the workflow but cannot replicate numeric outputs.
-------------------------------------------------------------------------------
*/

clear
clear matrix

cd "../../data_moments"

* Load cleaned plant-product panel (constructed in sectoral_moments.do).
use plant_product.dta, clear

* Keep observations with valid sales and labor-cost variables.
drop if sales<=0 | sales==. | vl==.

* Construct value added (sales minus intermediate inputs and other costs).
gen va = sales - vm - ve - vs - voth

* Drop impossible/negative value-added observations.
drop if va<0

* Collapse to plant-year level: one observation per plant-year.
collapse (mean) va dsales sales vl vm, by(year stno)

* Domestic-sales share at the plant-year level.
gen dshare = dsales/sales
drop if dshare>1

#delimit ;

/*
Compute top 1% / 5% plant shares year-by-year:
- s1/s5: share of domestic sales by top plants.
- lv1/lv5: share of wage bill by the same top plants.
*/
sort stno year;
tab year;
egen yobs = count(year), by(year);

gsort year -dsales;
by year: gen ranks = [_n];
egen ts = sum(dsales), by(year);

gen ranks1 = 1 if ranks/yobs<=0.01;
gen ranks5 = 1 if ranks/yobs<=0.05;

gen s1 = dsales*ranks1/ts;
gen s5 = dsales*ranks5/ts;

egen tl = sum(vl), by(year);

gen lv1 = vl*ranks1/tl;
gen lv5 = vl*ranks5/tl;

collapse (sum) s1 s5 lv1 lv5, by(year);

* Report annual moments used for Table 1 comparison.
tabstat s1 s5 lv1 lv5, by(year) stat(mean);

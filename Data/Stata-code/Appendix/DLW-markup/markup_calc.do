/*
-------------------------------------------------------------------------------
markup_calc.do (Appendix DLW markup procedure)
Purpose: Compute plant-level labor-input markup following De Loecker-Warzynski.

[DATA AVAILABILITY NOTICE]
Requires confidential datasets not shipped in this repo:
- ../../../raw data/taiwan.dta
- production intermediate files generated from restricted plant data.
Without these files, replication is not possible.
-------------------------------------------------------------------------------
*/

cd "../../../production"
clear
clear matrix
set mem 1g

* Load sector-specific production-function coefficients estimated in MATLAB.
insheet sic2 al ak am al2 ak2 am2 akl alm akm using "../../../production/prod_est_g.csv"
save prod_est_g.dta, replace

* Load confidential plant panel and construct capital stock recursively.
use "../../../raw data/taiwan.dta", clear
drop _merge
egen miny = min(year), by(stno)
sort stno year
by stno: replace vk = vk[_n-1]+invest if year>=89 & miny<89

* Keep replication years and build 2-digit sector code.
keep if year>=89
gen sic2 = floor(did/100)
order stno year sic2 tbusincome vm vl ql vk ve vs

* Construct key logs.
gen va = sales-vm-ve-voth
gen lnv = log(va)
gen lnq = log(sales)
gen lnl = log(vl)
gen lnk = log(vk)
gen lnm = log(vm)

egen sicy = group(sic2 year)
drop if stno==.
keep if sic2<=31 & sic2>=10

* Merge in estimated production-function coefficients.
merge m:1 sic2 using prod_est_g.dta

* Flexible polynomial terms.
gen lnk2 = lnk^2
gen lnl2 = lnl^2
gen lnm2 = lnm^2
gen lnk3 = lnk^3
gen lnl3 = lnl^3
gen lnm3 = lnm^3

gen lnkl = lnk*lnl
gen lnkm = lnk*lnm
gen lnlm = lnm*lnl

gen lnk2l = lnk2*lnl
gen lnk2m = lnk2*lnm
gen lnl2m = lnl2*lnm
gen lnl2k = lnl2*lnk
gen lnm2k = lnm2*lnk
gen lnm2l = lnm2*lnl

* Plant-year labor elasticity and markup: mu = elasticity / revenue share.
gen e_l = al+2*al2*lnl+alm*lnm+akl*lnk
tabstat e_l, by(sic2) stat(p25 p50 p75)

gen s_l = vl/exp(lnq)
sum s_l, d
gen mkup_l = e_l/s_l

sum mkup_l, d
replace mkup_l = . if mkup_l>=10 | mkup_l<=1
tabstat mkup_l, stat(p5 p25 p50 p75 p95)

* Save cleaned plant-level markup file for sector aggregation script.
sort stno year
drop _merge
rename mkup_l mkup_lg
order stno year mkup_lg sales vm vl vk ve voth vs did
keep stno year mkup_lg sales vm vl vk ve voth vs did
save plant_mkup_g.dta, replace

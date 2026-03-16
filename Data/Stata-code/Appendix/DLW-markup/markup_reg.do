/*
-------------------------------------------------------------------------------
markup_reg.do (Appendix DLW markup aggregation)
Purpose: Merge plant-product data with DLW plant markups and build sector markups.

[DATA AVAILABILITY NOTICE]
Depends on confidential/restricted files:
- ../../../data_moments/plant_product.dta
- ../../../production/plant_mkup_g.dta
If missing, script cannot reproduce Appendix results.
-------------------------------------------------------------------------------
*/

clear
clear matrix
set mem 1g

cd "../../../data_moments"

* Plant-product panel from main data construction.
use plant_product.dta, clear
drop _merge

* Merge plant-level markups from markup_calc.do.
merge m:1 stno year using "../../../production/plant_mkup_g.dta"
table _merge year, c(sum sales)

keep if _merge==3
drop _merge

* Keep reliable product-to-plant matching sample.
egen fsales = sum(psales), by(stno year)
drop if last>0
keep if abs((fsales-sales)/sales)<.20

sum mkup_lg, d

* Weighted mean inverse markup by year.
gen imkup_lg = 1/mkup_lg
tabstat imkup_lg [aweight = psales], by(year) stat(mean)

* Predict inverse markup from import share control.
gen lnmkup = log(mkup_lg)
gen lnshare = log(tshare)
reg imkup_lg tshare i.year
predict imkup, xb

* Aggregate product-year sector markups (harmonic weighting).
replace imkup = . if imkup<0.1
gen mk_missing = 0 if imkup==.
replace mk_missing = 1 if imkup~=.

egen norm_tshare = sum(tshare*mk_missing), by(prodno year)
gen norm_share = tshare*mk_missing/norm_tshare

gen imkup_sh = imkup*norm_share
collapse (sum) isec_markup = imkup_sh (sd) sd_share = norm_share, by(prodno year)
gen sec_markup = 1/isec_markup

* Dispersion-mean relationship reported in appendix discussion.
reg sec_markup sd_share, r

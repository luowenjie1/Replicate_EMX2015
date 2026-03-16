/*
-------------------------------------------------------------------------------
markup_moments.do
Purpose: Generate markup-related moments used in Table 1 / Table 2.

[DATA AVAILABILITY NOTICE]
This script depends on confidential file `plant_product.dta` under
`../../data_moments/` (built from restricted micro data). Missing file implies
results cannot be replicated numerically.
-------------------------------------------------------------------------------
*/

clear
clear matrix
set mem 1g

cd "../../data_moments"

use plant_product.dta, clear
drop _merge

* Construct value added and keep economically meaningful sample.
gen va = sales - vm - ve - vs - voth
keep if va>0

* Ensure product-level totals align with plant-level totals.
drop if last>0
egen fsales_check = sum(psales), by(stno year)
keep if abs((fsales-sales)/sales)<.20

* -----------------------------------------------------------------------------
* 1) Labor share regression for NON-exporters (used in text/table discussion).
* -----------------------------------------------------------------------------
preserve
collapse (sum) dsales esales (mean) vl va tshare [aweight = psales], by(year stno)
drop if esales>0

gen lshare = vl/va
sum lshare, d
keep if lshare>0 & lshare<1
reg lshare tshare i.year, r
restore

* -----------------------------------------------------------------------------
* 2) Model-implied markups at product-year level.
*    imarkup = marginal cost share inferred from model parameters.
* -----------------------------------------------------------------------------
gen gamma = 10.5
gen theta = 1.24
gen imarkup = (gamma-1)/gamma - (1/theta-1/gamma)*tshare
gen markup = 1/imarkup
sum markup, d

gen logmarkup = log(markup)
sum logmarkup

* Build normalized product-year shares over non-missing markups.
gen mk_missing = 0 if markup==.
replace mk_missing = 1 if markup~=.
egen norm_tshare = sum(tshare*mk_missing), by(prodno year)
gen norm_share = tshare*mk_missing/norm_tshare

* Sector-level harmonic aggregation of markups.
gen imkup_sh = (1/markup)*norm_share
collapse (sum) isec_markup = imkup_sh (sd) sd_share = norm_share, by(prodno year)
gen sec_markup = 1/isec_markup
sum sec_markup, d

gen logsec_markup = log(sec_markup)
sum logsec_markup

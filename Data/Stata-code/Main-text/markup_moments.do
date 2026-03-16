****This STATA do file creates the markup moments in Table 1/2***

clear

clear matrix

set mem 1g

*************Plant-product Level Data*****************
cd "../../data_moments"

use plant_product.dta, clear
drop _merge

************Construct labor share*************
gen va = sales - vm - ve -vs - voth
keep if va>0

************Drop if product-level and plant level data doesn't match******
drop if last>0
egen fsales_check = sum(psales), by(stno year)
keep if abs((fsales-sales)/sales)<.20

*************Collapse on firm year
preserve
collapse (sum) dsales esales (mean) vl va tshare [aweight = psales] , by(year stno)

**keep non-exporters
drop if esales>0
**************Regression
gen lshare = vl/va
sum lshare, d
keep if lshare>0&lshare<1
reg lshare tshare i.year, r
restore

**projected markup from model (in Table 2)
gen gamma = 10.5
gen theta = 1.24
gen imarkup = (gamma-1)/gamma-(1/theta-1/gamma)*tshare
gen markup = 1/imarkup
sum markup, d
gen logmarkup = log(markup)
sum logmarkup

gen mk_missing = 0 if markup==.
replace mk_missing = 1 if markup~=.
egen norm_tshare = sum(tshare*mk_missing), by (prodno year)
gen norm_share = tshare*mk_missing/norm_tshare
gen imkup_sh = (1/markup)*norm_share
collapse (sum) isec_markup = imkup_sh (sd) sd_share = norm_share, by (prodno year)
gen sec_markup = 1/isec_markup
sum sec_markup, d
gen logsec_markup = log(sec_markup)
sum logsec_markup




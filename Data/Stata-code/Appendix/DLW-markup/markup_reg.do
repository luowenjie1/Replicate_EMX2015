clear

clear matrix

set mem 1g

*************1. Plant-product Level Data*****************
cd "../../../data_moments"

use plant_product.dta, clear
drop _merge
**********************2. Merge in Plant Level Data of Markups and "Productivity"*****************
/**only losing around 10% of sales**/
merge m:1 stno year using "../../../production/plant_mkup_g.dta"
table _merge year, c(sum sales)

keep if _merge==3
drop _merge

egen fsales = sum(psales), by(stno year) /*the plant-level data dropping was problematic*/
drop if last>0              /** drop  firms selling in 9x industries (other) **/
keep if abs((fsales-sales)/sales)<.20 /* drop guys w/ bad match of plant sales */

sum mkup_lg, d

gen imkup_lg = 1/mkup_lg
tabstat imkup_lg [aweight = psales], by (year) stat (mean)

*********************Construct sectoral level variables****************
gen lnmkup = log(mkup_lg)
gen lnshare = log(tshare)
reg imkup_lg tshare i.year
**projected markup
predict imkup, xb

/***sectoral markups****/
replace imkup = . if imkup<0.1
gen mk_missing = 0 if imkup==.
replace mk_missing = 1 if imkup~=.

egen norm_tshare = sum(tshare*mk_missing), by (prodno year)
gen norm_share = tshare*mk_missing/norm_tshare

gen imkup_sh = imkup*norm_share
collapse (sum) isec_markup = imkup_sh (sd) sd_share = norm_share, by (prodno year)
gen sec_markup = 1/isec_markup

reg sec_markup sd_share, r

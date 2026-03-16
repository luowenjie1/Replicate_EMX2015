
cd "C:\Users\yx42\Documents\Research Projects\Competition and Markup\Revision\production"
clear

clear matrix

set mem 1g

insheet sic2 al ak am al2 ak2 am2 akl alm akm using "C:\Users\yx42\Documents\Research Projects\Competition and Markup\Revision\production\prod_est_g.csv"

save prod_est_g.dta, replace

use "C:\Users\yx42\Documents\Research Projects\Competition and Markup\Revision\raw data\taiwan.dta", clear
drop _merge
*construct capital stock measure for 2000-2004 using Inventory method
egen miny=min(year), by (stno)
sort stno year
by stno: replace vk=vk[_n-1]+invest if year>=89 &miny<89

*tab year, gen (dy)
keep if year>=89
gen sic2=floor(did/100)
tab sic2 year
order stno year sic2 tbusincome vm vl ql vk ve vs

*construct value-added
gen va=sales-vm-ve-voth

gen lnv=log(va)
gen lnq=log(sales)
gen lnl=log(vl)
gen lnk=log(vk)
gen lnm=log(vm)

egen sicy=group(sic2 year)

drop if stno==.
**only keep industries 10 - 31
keep if sic2<=31&sic2>=10

merge m:1 sic2 using prod_est_g.dta

*measured "TFP" - which is in fact contaminated by markups
gen lnk2=lnk^2
gen lnl2=lnl^2
gen lnm2=lnm^2
 gen lnk3=lnk^3
 gen lnl3=lnl^3
 gen lnm3=lnm^3
 
gen lnkl=lnk*lnl
gen lnkm=lnk*lnm
gen lnlm=lnm*lnl

 gen lnk2l=lnk2*lnl
 gen lnk2m=lnk2*lnm
 gen lnl2m=lnl2*lnm
 gen lnl2k=lnl2*lnk
 gen lnm2k=lnm2*lnk
 gen lnm2l=lnm2*lnl

  *calculate input elas. for each plant-year
gen e_l=al+2*al2*lnl+alm*lnm+akl*lnk

tabstat e_l, by (sic2) stat (p25 p50 p75)


gen s_l=vl/exp(lnq)
sum s_l, d
gen mkup_l=e_l/s_l

sum mkup_l, d
replace mkup_l=. if mkup_l>=10|mkup_l<=1
tabstat mkup_l, stat (p5 p25 p50 p75 p95)

sort stno year
drop _merge

rename mkup_l mkup_lg
order stno year mkup_lg sales vm vl vk ve voth vs did
keep stno year  mkup_lg sales vm vl vk ve voth vs did
save plant_mkup_g.dta, replace

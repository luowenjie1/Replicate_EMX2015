/*
-------------------------------------------------------------------------------
trade_moments.do
Purpose: Construct import/export moments in Table 1.

[DATA AVAILABILITY NOTICE]
This script requires confidential trade files under `../../import/`:
- import_final_89b.dta, import_final_91b.dta,
  import_final_92b.dta, import_final_93b.dta.
Without these files, the script cannot reproduce the table moments.
-------------------------------------------------------------------------------
*/

clear

cd "../../import"

* 1) Pool annual import-export sector panel.
use import_final_89b.dta
append using import_final_91b.dta
append using import_final_92b.dta
append using import_final_93b.dta

* Remove "other" products and excluded industries.
g last2 = prodno-floor(prodno/100)*100
drop if last2>=90
g first2 = floor(prodno/100000)
drop if first2==17 | first2==23

tab year

* 2) Persistence check of import penetration.
sort prodno year
by prodno: gen limport_pen = import_pen[_n-1]
corr import_pen limport_pen

* 3) Build aggregate yearly totals and sector shares.
gen tindsales = import + dsalesc
egen tds = sum(dsalesc), by(year)
egen tes = sum(esalesc), by(year)
egen ts = sum(tindsales), by(year)
egen ti = sum(import), by(year)

tabstat tds tes ts ti, by(year) stat(mean)

gen sj_d = dsalesc/tds
gen sj_i = import/ti
gen sj_e = esalesc/tes
gen sj_t = tindsales/ts

* Co-movement of import share and domestic-sales share.
reg sj_i sj_d i.year

* 4) Import dispersion index and normalized Grubel-Lloyd index.
gen lambdaj = (1-import_pen)
gen temp = (lambdaj)*sj_t
egen tlambda = sum(temp), by(year)
gen sl = 1-lambdaj*(1-lambdaj)/(tlambda*(1-tlambda))

gen gln = 1-abs(sj_e-sj_i)/(sj_e+sj_i)

* Report key trade moments.
tabstat sl gln import_pen [aw = tindsales], by(year) stat(mean)

save import_finalb.dta, replace

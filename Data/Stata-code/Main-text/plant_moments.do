
****This STATA do file creates the plant moments in Table 1***

cd "C:\Users\yx42\Documents\Research Projects\Competition and Markup\Revision\data_moments"

*****************************construct plant level moments************************************
use plant_product.dta, clear

drop if sales<=0|sales==.|vl==.

g va = sales - vm - ve - vs -voth

drop if va<0

*********************************************************

collapse (mean) va dsales sales vl vm, by (year stno)

gen dshare=dsales/sales
drop if dshare>1

#delimit;

/*calculate year by year*/
sort stno year;
tab year;
egen yobs=count(year), by (year);

gsort year -dsales;
by year: gen ranks=[_n]; 
egen ts=sum(dsales), by (year);

gen ranks1=1 if ranks/yobs<=0.01;
gen ranks5=1 if ranks/yobs<=0.05;

gen s1=dsales*ranks1/ts;
gen s5=dsales*ranks5/ts;

egen tl=sum(vl), by (year);

gen lv1=vl*ranks1/tl;
gen lv5=vl*ranks5/tl;

collapse (sum) s1 s5 lv1 lv5,  by (year);

tabstat s1 s5 lv1 lv5, by (year) stat (mean);

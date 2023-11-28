mata: mata clear
capture program drop _all
do gcomp_effects.ado

import delimited rhc.csv, clear

replace dschdte = "" if dschdte == "NA"
replace dthdte = "" if dthdte == "NA"
destring dschdte, replace
destring dthdte, replace
gen len_stay = dschdte - sadmdte
replace len_stay = dthdte - sadmdte if len_stay == .

gen rhc = swang1 == "RHC"
drop v1 dthdte lstctdte dschdte death t3d30 dth30 surv2md1 sadmdte ptid adld3p urin1 cat2 swang1


foreach var of varlist _all {
	if substr("`: type `var''", 1, 3) == "str" {
		tempvar x
		encode `var', gen(`x')
		drop `var'
		rename `x' `var'
		replace `var' = `var' - 1
	}
}

keep len_stay rhc age sex income cat1
mkspline age_spl = age, cubic nknots(5)

timer clear

* Model 1
timer on 1
qui regress len_stay b0.rhc b0.cat1 b0.sex c.age b0.income
gcomp_effects rhc, values(0 1) ate
return list
timer off 1

* Model 2
timer on 2
qui regress len_stay b0.rhc##b0.cat1 b0.sex c.age b0.income
gcomp_effects rhc, values(0 1) ate 
return list
timer off 2

* Model 3
timer on 3
qui regress len_stay b0.rhc##b0.cat1 b0.sex c.age_spl* b0.income
gcomp_effects rhc, values(0 1) ate 
return list
timer off 3

* Model 4
timer on 4
qui regress len_stay b0.rhc##b0.cat1 b0.rhc##b0.sex c.age_spl* b0.income
gcomp_effects rhc, values(0 1) ate 
return list
timer off 4

timer list
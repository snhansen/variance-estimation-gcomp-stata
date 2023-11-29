cd "C:\Users\au234616\Dropbox\work\gcomp"
mata: mata clear
capture program drop _all
do gcomp_effects_terza.ado
do sim_outcome.do

local n_obs = 500
local alpha = 0
local beta_x = -1
local beta_b = 0
local beta_xb = -5
local beta_c = 0
local beta_xc = 0
local beta_bc = 2
local beta_xbc = 0
local mean_b = 0
local mean_c = 0
local logit_xc = 3
local ate = `beta_x'+`beta_xb'*`mean_b'
local rep_n = 10000
local filename "ate_md_n`n_obs'_logit`logit_xc'_terza"
local seed = floor((100000)*runiform())
set seed `seed'


cd data
cap postclose simdata

capture confirm file "`filename'.dta"
local simdata = _rc
disp `simdata'
if `simdata' == 0 {
	capture confirm file "`filename'_result.dta"
	local resdata = _rc
	if `resdata' > 0 {
		shell rename "`filename'.dta" "`filename'_result.dta"
		use "`filename'_result.dta", clear
		local n_current = _N
	}
	else {
		use "`filename'_result.dta", clear
		append using "`filename'"
		local n_current = _N
		save "`filename'_result.dta", replace
	}
}

else {
	local n_current = 0
}


postfile simdata true_ate double(ate_est ate_se1 ate_se2) logit_xc n_obs seed using "`filename'", replace every(5) 

forvalues i=1/`=`rep_n'-`n_current'' {
	disp `i' " / " `=`rep_n'-`n_current''
	qui sim_cont_depend `n_obs' `mean_b' `mean_c' `logit_xc' `alpha' `beta_x' `beta_b' `beta_c' `beta_xb' `beta_xc' `beta_bc' `beta_xbc'
	qui glm y b0.x##c.b c.c
	qui gcomp_effects x, values(0 1) ate notable
	local ate_se1 = r(ate_se)
	local ate_se2 = r(ate_se2)
	local ate_est = r(ate)
	post simdata (`ate') (`ate_est') (`ate_se1') (`ate_se2') (`logit_xc') (`n_obs') (`seed')
}
postclose simdata

capture confirm file "`filename'_result.dta"
if _rc == 0 {
	use "`filename'_result.dta", clear
	append using "`filename'"
	save "`filename'_result.dta", replace
	erase "`filename'.dta"
}
else {
	shell rename "`filename'.dta" "`filename'_result.dta"
}

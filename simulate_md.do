mata: mata clear
capture program drop get_bvar
capture program drop gcomp_effects
capture program drop output_line
capture program drop sim_cont_outcome
capture program drop sim_bin_outcome
do gcomp_effects.ado
do get_bvar.ado
do sim_outcome.do

local n_obs = `1'
local alpha = 0
local beta_x = 1
local beta_c = -1
local beta_xc = 3
local mean_c = 2
local mean_po1 = `alpha' + `beta_x' + (`beta_c' + `beta_xc')*`mean_c'
local mean_po0 = `alpha' + `beta_c'*`mean_c'
local ate = `mean_po1' - `mean_po0'
local bstrap_n = `2'
local rep_n = 10000
local filename "ate_md_n`n_obs'_k`bstrap_n'"
local seed = floor((100000)*runiform())
set seed `seed'

cd data
cap postclose simdata

capture confirm file "`filename'.dta"
local simdata = _rc
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

postfile simdata true_ate double(ate_est ate_se1 ate_se2 ate_p25 ate_p975 t1 t2) bstrap_n n_obs seed using "`filename'", replace every(5) 

forvalues i=1/`=`rep_n'-`n_current'' {
	disp `i' " / " `=`rep_n'-`n_current''
	
	qui sim_cont_outcome `n_obs' `mean_c' `alpha' `beta_x' `beta_c' `beta_xc'
	qui glm y b0.x##c.c
	
	timer clear

	timer on 1
	qui gcomp_effects x, values(0 1) ate
	timer off 1
	local ate_se1 = r(ate_se)
	local ate_est = r(ate)
	
	timer on 2
	get_bvar x, values(0 1) reps(`bstrap_n')
	timer off 2
	local ate_se2 = r(ate_se)
	local ate_p25 = r(ate_p25)
	local ate_p975 = r(ate_p975)
	qui timer list

	post simdata (`ate') (`ate_est') (`ate_se1') (`ate_se2') (`ate_p25') (`ate_p975') (`r(t1)') (`r(t2)') (`bstrap_n') (`n_obs') (`seed')
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
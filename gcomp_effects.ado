* Estimation of potential outcome means and average treatment effects and their asymptotic (co-)variance matrix.
* This requires that a GLM has been fit (the Q-model) and is available in the return list.

program define gcomp_effects, rclass
	syntax varname, VALues(numlist min=1) [ATE RTE SUBgrp(varlist min=1 max=1) Level(real 95) NOTABle]
	local exposure `varlist'
	local ate_indicator = "`ate'" == "ate"
	local rte_indicator = "`rte'" == "rte"
	local subgrp_indicator = "`subgrp'" != ""
	local level = `level'/100
	
	mata {
		ate_indicator = st_local("ate") == "ate"
		rte_indicator = st_local("rte") == "rte"
		subgrp_indicator = st_local("subgrp") != ""
	}
	
	* Checking that a GLM is available in the return list.
	if ! inlist(e(cmd), "glm", "regress", "logit", "logistic", "poisson", "nbreg", "cloglog") {
 		disp as error "No GLM available in memory."
		exit
 	}

	* Need at least two exposure values to calculate treatment effects.
	if `: list sizeof values' == 1 & (`ate_indicator' | `rte_indicator') {
		disp as error "Unable to calculate treatment effects with only one exposure value specified."
		exit
	}
	
	* Checking whether the subgroup variable has any observations with a value of 1.
	if `subgrp_indicator' {
		quietly count if `subgrp' == 1
		if ! r(N) {
			disp as error "There are no observations with a value of 1 in the subgroup."
			exit
		}
	}

	* We take advantage of Stata's in-built functions to calculate relevant
	* quantities for a GLM. These are the glim_v* and glim_l* commands.
	if e(cmd) == "regress" {
		local varfunc = "glim_v1"
		local linkfunc = "glim_l01"
		local SGLM_m = 1
	}
	else if e(cmd) == "poisson" {
		local varfunc = "glim_v3"
		local linkfunc = "glim_l03"
		local SGLM_m = 1
	}
	else if e(cmd) == "nbreg" {
		local varfunc = "glim_v6"
		local linkfunc = "glim_l03"
		local SGLM_m = 1
	}
	else if e(cmd) == "probit" {
		local varfunc = "glim_v2"
		local linkfunc = "glim_l08"
		local SGLM_m = 1
	}
	else if inlist(e(cmd), "logit", "logistic") {
		local varfunc = "glim_v2"
		local linkfunc = "glim_l02"
		local SGLM_m = 1
	}
	else if e(cmd) == "cloglog" {
		local varfunc = "glim_v2"
		local linkfunc = "glim_l07"
		local SGLM_m = 1
	}
	else {
		local varfunc = "`e(varfunc)'"
		local linkfunc = "`e(link)'"
		local SGLM_m = 1
	}
	
	tempvar eta mu dmu d2mu v dv m
	quietly predict double `eta', xb // predict seems to drop all global macros, so we need to set SGLM_m (binomial denominator) after each predict
	 
	global SGLM_m = `SGLM_m'
	`linkfunc' 1 `eta' `mu'
	`linkfunc' 2 `eta' `mu' `dmu'
	`linkfunc' 3 `eta' `mu' `d2mu'
	`varfunc' 1 `eta' `mu' `v'
	`varfunc' 2 `eta' `mu' `dv'

	foreach var in `dmu' `d2mu' `v' `dv' {
		capture gen double `var' = `var' // some times glim_v* returns a scalar and other times a variable. this ensures that `var' is always a variable
	}
	
	gen double `m' = (`d2mu'/`v' - `dv'*`dmu'*`dmu'/(`v'*`v'))*(`e(depvar)'-`mu') - `dmu'*`dmu'/`v' 
	
	local 2 = e(cmdline)
	gettoken 1 2: 2, parse(" ")
	gettoken 1 2: 2, parse(" ")
	gettoken parlist 2: 2, parse(",")
	if regexm("`parlist'", "if") {
		local parlist = strrtrim(usubstr("`parlist'", 1, `=ustrpos("`parlist'", "if") - 1'))
	}
	fvexpand `parlist'
	local parlist = r(varlist)

	*disp "`parlist'"
	tempvar touse
	gen `touse' = e(sample)
	
	* Reading the relevant quantities into Mata.
	mata {
		n_obs = `e(N)'
		x_vals = tokens("`values'")
		eb = st_matrix("e(b)")
		X_obs = st_data(., "`parlist'", "`touse'") , J(n_obs, 1, 1)
		y = st_data(., "`e(depvar)'", "`touse'")
		M_vec = st_data(., "`m'", "`touse'")
		v = st_data(., "`v'", "`touse'")
		dmu = st_data(., "`dmu'", "`touse'")
		mu_hat = st_data(., "`mu'", "`touse'")
		M_hat = (X_obs:*M_vec)'*X_obs / rows(X_obs)
		M_invhat = pinv(M_hat)
		beta_infl = -(X_obs:*dmu:/v)*M_invhat':*(y-mu_hat)
		X_int = asarray_create("real") // Array of design matrices after intervening on the exposure variable.
		po = J(n_obs, length(x_vals), .)
		dmu_int = J(n_obs, length(x_vals), .)
	}

	* Creating the design matrices after intervening on the exposure variable.
	* We exploit Stata's st_data() function which updates all relevant columns 
	* in the design matrix after changing the exposure variable, i.e. it handles
	* interaction terms automatically for us.
	local i = 1
	foreach x_val in `values' {
		tempvar x_bk
		rename `exposure' `x_bk'
		gen `exposure' = `x_val'
		tempvar eta mu dmu v
		quietly predict double `eta', xb
		global SGLM_m = `SGLM_m'
		`linkfunc' 1 `eta' `mu'
		`linkfunc' 2 `eta' `mu' `dmu'
		mata: asarray(X_int, `i', (st_data(., "`parlist'", "`touse'") , J(n_obs, 1, 1)))
		mata: po[, `i'] = st_data(., "`mu'", "`touse'")
		mata: dmu_int[, `i'] = st_data(., "`dmu'", "`touse'")
		drop `exposure'
		rename `x_bk' `exposure'
		local ++i
	}
	mac drop SGLM_m
	
	* Calculating the predicted potential outcomes, their means, and the asymptotic covariance matrix.
	* If specified, we also calculate the ATE and RTE.
	* If the subgroup option is used, then all calculations will be done for the specified subgroup.
	mata {
		if(!subgrp_indicator) 
		{
			po_infl = J(n_obs, length(x_vals), .)
			po_infl2_1 = J(n_obs, length(x_vals), .)
			po_infl2_2 = J(n_obs, length(x_vals), .)
			po_mean = J(1, length(x_vals), .)
			for (i=1; i<=length(x_vals); i++) {
				po_mean[i] = mean(po[, i])
				po_infl[, i] = po[, i] :- po_mean[, i] + beta_infl*mean(asarray(X_int, i):*dmu_int[, i])'
				po_infl2_1[, i] = po[, i] :- po_mean[, i]
				po_infl2_2[, i] = beta_infl*mean(asarray(X_int, i):*dmu_int[, i])'
			}
			po_var = po_infl'*po_infl/n_obs
			po_var2 = po_infl2_1'*po_infl2_1/n_obs + po_infl2_2'*po_infl2_2/n_obs // Terza's expression.
			po_se = sqrt(diagonal(po_var):/n_obs)
			po_se2 = sqrt(diagonal(po_var2):/n_obs)
			st_matrix("po_mean", po_mean')
			st_matrix("po_var", po_var)
			st_matrix("po_se", po_se)
			st_matrix("po_var2", po_var2)
			st_matrix("po_se2", po_se2)
			
			if (ate_indicator) 
			{
				ate = po_mean[2] - po_mean[1]
				ate_var = (-1, 1)*po_var[1::2, 1::2]*(-1 \ 1)
				ate_var2 = (-1, 1)*po_var2[1::2, 1::2]*(-1 \ 1)
				st_local("ate", strofreal(ate))
				st_local("ate_se", strofreal(sqrt(ate_var/n_obs)))
				st_local("ate_se2", strofreal(sqrt(ate_var2/n_obs)))
			}
			
			if (rte_indicator)
			{
				if (po_mean[1] <= 0 | po_mean[2] <= 0)
				{
					_error(3300, "Potential outcome means are not positive, so logarithm can't be applied.")   
				}
				log_rte = log(po_mean[2]/po_mean[1])
				log_rte_var = (-1/po_mean[1], 1/po_mean[2])*po_var[1::2, 1::2]*(-1/po_mean[1] \ 1/po_mean[2])
				st_local("log_rte", strofreal(log_rte))
				st_local("log_rte_se", strofreal(sqrt(log_rte_var/n_obs)))
			}
		}
		
		if(subgrp_indicator)
		{
			subgrp = st_data(., "`subgrp'")
			include = subgrp :== 1
			po_sub = J(n_obs, length(x_vals), .)
			po_sub_mean = J(1, length(x_vals), .)
			for (i=1; i<=length(x_vals); i++) {
				po_sub[, i] = po[, i]:*include
				po_sub_mean[i] = mean(select(po_sub[, i], include))
			}
			po_sub_infl = J(n_obs, length(x_vals), .)
			p_sub = mean(include)			
			for (i=1; i<=length(x_vals); i++) {
				po_sub_infl[, i] = po_sub[, i]:/p_sub :- J(n_obs, 1, po_sub_mean[, i]/p_sub):*include + beta_infl*mean(select(asarray(X_int, i):*dmu_int[, i], include))'
			}
			po_sub_var = po_sub_infl'*po_sub_infl/n_obs
			po_sub_se = sqrt(diagonal(po_sub_var):/n_obs)   
			st_matrix("po_sub_mean", po_sub_mean')
			st_matrix("po_sub_var", po_sub_var)
			st_matrix("po_sub_se", po_sub_se)
			
			if(ate_indicator)
			{
				ate_sub = po_sub_mean[2] - po_sub_mean[1]
				ate_sub_var = (-1, 1)*po_sub_var[1::2, 1::2]*(-1 \ 1)
				st_local("ate_sub", strofreal(ate_sub))
				st_local("ate_sub_se", strofreal(sqrt(ate_sub_var/n_obs)))
			}
			
			if (rte_indicator)
			{
				if (po_sub_mean[1] <= 0 | po_sub_mean[2] <= 0)
				{
					_error(3300, "Potential outcome means are not positive, so logarithm can't be applied.")   
				}
				log_rte_sub = log(po_sub_mean[2]/po_sub_mean[1])
				log_rte_sub_var = (-1/po_sub_mean[1], 1/po_sub_mean[2])*po_sub_var[1::2, 1::2]*(-1/po_sub_mean[1] \ 1/po_sub_mean[2])
				st_local("log_rte_sub", strofreal(log_rte_sub))
				st_local("log_rte_sub_se", strofreal(sqrt(log_rte_sub_var/n_obs)))
			}
		}
	}

	* Writing the output table.
	if "`notable'" != "notable" {
		local maxlen: label (`exposure') maxlength
		if `ate_indicator' | `rte_indicator' {
			local maxlen = max(`maxlen', length(`"`: label (`exposure') `: word 2 of `values''' vs `: label (`exposure') `: word 1 of `values'''"'))
		}
		
		local ndash = 13

		if c(linesize) >= 80 {
			local usablerows = c(linesize)-68
			if `maxlen' < `usablerows'-1 {
				local ndash = max(`maxlen' + 3, `ndash')
			}
			else {
				local ndash = max(`usablerows' + 2, `ndash')
			}
		}
		
		local abexposure = abbrev("`exposure'", `=`ndash'-1')
		local abdepvar = abbrev("`e(depvar)'", `=`ndash'-1')
		
		local i = 1
		foreach val of local values {
			if c(linesize) > 100 {
				local abname`i' = "`: label (`exposure') `val''"
			}
			else if c(linesize) >= 80 {
				local abname`i' = abbrev("`: label (`exposure') `val''", c(linesize)-68)
			}
			else {
				local abname`i' = abbrev("`: label (`exposure') `val''", 11)
			}
			local ++i
		}
		
		if `subgrp_indicator' {
			local ndash = max(`ndash', 15)
		}
		disp "{hline `ndash'}{c TT}{hline 64}"
		disp as text "{ralign `=`ndash'-1':`abdepvar'} {c |}	  Coef.   Std. Err.	  z	P>|z|	 [95% Conf. Interval]"
		disp "{hline `ndash'}{c +}{hline 64}"
		if !`subgrp_indicator' {
			disp as result "{lalign `ndash':POM}{c |}"
			disp as text "{ralign `=`ndash'-1':`abexposure'} {c |}"
			local i = 1
			foreach val of local values {
				output_line "`abname`i''" po_mean[`i',1] po_se[`i',1] `level' `=`ndash'-2' 0
				local ++i
			}
			
			if `ate_indicator' {
				disp "{hline `ndash'}{c +}{hline 64}"
				disp as result "{lalign `ndash':ATE}{c |}"
				disp as text "{ralign `=`ndash'-1':`abexposure'} {c |}"
				local vstext `"`abname2' vs `abname1'"' 
				if length("`vstext'") <= `=`ndash'-2' {
					output_line "`vstext'" `ate' `ate_se' `level' `=`ndash'-2' 0
				}
				else {
					disp as text "{ralign `=`ndash'-2': `abname2'}  {c |}"
					disp as text "{ralign `=`ndash'-2': vs}  {c |}"
					output_line "`abname1'" `ate' `ate_se' `level' `=`ndash'-2' 0
				}
				
			}
			if `rte_indicator' {
				disp "{hline `ndash'}{c +}{hline 64}"
				disp as result "{lalign `ndash':RTE}{c |}"
				disp as text "{ralign `=`ndash'-1':`abexposure'} {c |}"
				local vstext `"`abname2' vs `abname1'"' 
				if length("`vstext'") <= `=`ndash'-2' {
					output_line "`vstext'" `log_rte' `log_rte_se' `level' `=`ndash'-2' 1
				}
				else {
					disp as text "{ralign `=`ndash'-2': `abname2'}  {c |}"
					disp as text "{ralign `=`ndash'-2': vs}  {c |}"
					output_line "`abname1'" `log_rte' `log_rte_se' `level' `=`ndash'-2' 1
				}
			}
			disp "{hline `ndash'}{c BT}{hline 64}"
			
		}
		else {
			disp as result "{lalign `ndash':POM (subgroup)}{c |}"
			disp as text "{ralign `=`ndash'-1':`abexposure'} {c |}"
			local i = 1
			foreach val of local values {
				output_line "`abname`i''" po_sub_mean[`i',1] po_sub_se[`i',1] `level' `=`ndash'-2' 0
				local ++i
			}
			if `ate_indicator' {
				disp "{hline `ndash'}{c +}{hline 64}"
				disp as result "{lalign `ndash':ATE (subgroup)}{c |}"
				disp as text "{ralign `=`ndash'-1':`abexposure'} {c |}"
				local vstext `"`abname2' vs `abname1'"' 
				if length("`vstext'") <= `=`ndash'-2' {
					output_line "`vstext'" `ate_sub' `ate_sub_se' `level' `=`ndash'-2' 0
				}
				else {
					disp as text "{ralign `=`ndash'-2': `abname2'}  {c |}"
					disp as text "{ralign `=`ndash'-2': vs}  {c |}"
					output_line "`abname1'" `ate_sub' `ate_sub_se' `level' `=`ndash'-2' 0
				}
				
			}
			if `rte_indicator' {
				disp "{hline `ndash'}{c +}{hline 64}"
				disp as result "{lalign `ndash':RTE (subgroup)}{c |}"
				disp as text "{ralign `=`ndash'-1':`abexposure'} {c |}"
				local vstext `"`abname2' vs `abname1'"' 
				if length("`vstext'") <= `=`ndash'-2' {
					output_line "`vstext'" `log_rte_sub' `log_rte_sub_se' `level' `=`ndash'-2' 1
				}
				else {
					disp as text "{ralign `=`ndash'-2': `abname2'}  {c |}"
					disp as text "{ralign `=`ndash'-2': vs}  {c |}"
					output_line "`abname1'" `log_rte_sub' `log_rte_sub_se' `level' `=`ndash'-2' 1
				}
			}
			disp "{hline `ndash'}{c BT}{hline 64}"
		}
	}
	
	* Adding to the return list
	if ! `subgrp_indicator' {
		return matrix po_mean = po_mean
		return matrix po_var = po_var
		return matrix po_se = po_se
		if `ate_indicator' {
			return scalar ate = `ate'
			return scalar ate_se = `ate_se'
			return scalar ate_se2 = `ate_se2'
		}
		if `rte_indicator' {
			return scalar log_rte = `log_rte'
			return scalar log_rte_se = `log_rte_se'
		}
	}
	else {
		return matrix po_sub_mean = po_sub_mean
		return matrix po_sub_var = po_sub_var
		return matrix po_sub_se = po_sub_se
		if `ate_indicator' {
			return scalar ate_sub = `ate_sub'
			return scalar ate_sub_se = `ate_sub_se'
		}
		if `rte_indicator' {
			return scalar log_rte_sub = `log_rte_sub'
			return scalar log_rte_sub_se = `log_rte_sub_se'
		}
	}   
end

* Stata functions
program output_line 
	args vname est se lev ndash exp
	local z = `est'/`se'
	local q = -invnormal((1-`lev')/2)
	if `exp' == 1 {
		disp as text "{ralign `ndash':`vname'}  {c |}" /*
		*/ "  " as result %9.0g exp(`est') "  " %9.0g `se' " " %8.2f `z' /*
		*/ "   " %5.3f 2*(1-normal(abs(`z'))) "	" %9.0g exp(`est'-`q'*`se') "   " %9.0g exp(`est'+`q'*`se')
	}
	else {
		disp as text "{ralign `ndash':`vname'}  {c |}" /*
		*/ "  " as result %9.0g `est' "  " %9.0g `se' " " %8.2f `z' /*
		*/ "   " %5.3f 2*(1-normal(abs(`z'))) "	" %9.0g `est'-`q'*`se' "   " %9.0g `est'+`q'*`se'
   }
end

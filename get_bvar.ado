program define get_bvar, rclass 
	syntax varname, VALues(numlist min=1) [Reps(integer 100) Seed(int 0) Level(real 95)]
	local exposure `varlist'
	local regcmd `e(cmdline)'
	local n_obs `e(N)'
	preserve
	
	* Creating a list of parameters in the regression model.
	local colnames: colfullnames e(b)
	tempname params
	mat `params' = e(b)

	local i = 0
	local parlist
	
	if e(cmd) == "regress" {
		foreach name of local colnames {
			local ++i
			if "`name'" != "_cons" && `params'[1,`i'] != 0 {
				local parlist `parlist' `name'
			}
		}
	}
	else {
		foreach name of local colnames {
			local ++i
			tokenize `name', parse(":")
			if "`3'" != "_cons" && `params'[1,`i'] != 0 {
				local parlist `parlist' `3'
			}
		}
	}
	
	* Creating a list of variables used in the regression model.
	local varlist

	foreach name of local parlist {
		tokenize `name', parse("#")
		local i = 1
		while "``=`i'''" != "" {
			if "``i''" != "#" {
				local temp = substr("``i''", strpos("``i''", ".") + 1, .)
				local varlist: list varlist | temp
			}
			local ++i
		}
	}
	* Specifying the link function.
	if e(cmd) == "regress" {
		mata: link = "Identity"
	}
	else if inlist(e(cmd), "poisson", "nbreg") {
		mata: link = "Log"
	}
	else if e(cmd) == "probit" {
		mata: link = "Probit"
	}
	else if inlist(e(cmd), "logit", "logistic") {
		mata: link = "Logit"
	}
	else if e(cmd) == "cloglog" {
		mata: link = "Complementary log-log"
	}
	else {
		mata: link = "`e(linkt)'"
	}
	
	* Setting the seed if specified.
	if `seed' {
		set seed `seed'
	}
	
	* Reading relevant quantities into Mata.
	mata {
		exp_vals = strtoreal(tokens(st_local("values")))
		n_obs = `n_obs'
		n_reps = `reps'
		eb = st_matrix("e(b)")
		par_est = select(eb, eb) 		// keeping only non-zero entries.
		ev = st_matrix("e(V)")
		par_var = select(select(ev, eb), eb') 	// keeping only non-zero entries.
		level = 1-(1-`level'/100)/2
		st_local("level", strofreal(level))
	}
	
	mata: data = st_data(., "`varlist'")
	mata: data_obs = asarray_create("real")
	
	local i = 1
	tempvar exp_bk
	gen `exp_bk' = `exposure'
	foreach val of local values {
		cap drop `exposure'
		gen `exposure' = `val'
		mata: asarray(data_obs, `i', (st_data(., "`parlist'"), J(n_obs, 1, 1)))
		local ++i
	}
	qui replace `exposure' = `exp_bk'
	
	mata: po_obs = pomeans(data_obs, J(n_obs, 1, par_est), J(n_obs, 1, 1), 1, n_obs, link)
	mata: data_samp = replsample("`varlist' `e(depvar)'", n_obs*n_reps)
	
	clear
	getmata (`varlist' `e(depvar)') = data_samp
	gen include = 1
	mata: regpars = J(0, length(par_est), .)
	
	forvalues i = 1/`reps' {
		qui cap `=regexr("`regcmd'", ", vce\(robust\)", "")' if inrange(_n, `=(`i'-1)*`n_obs'+1', `=`i'*`n_obs''), vce(robust)
		if _rc {
			disp "Can't use bootstrap sample"
			mata: n_reps = n_reps - 1
			qui replace include = 0 if inrange(_n, `=(`i'-1)*`n_obs'+1', `=`i'*`n_obs'')
		}
		else {
			if `e(rank)'<4 | `e(converged)' == 0 {
				disp "Can't use bootstrap sample"
				mata: n_reps = n_reps - 1
				qui replace include = 0 if inrange(_n, `=(`i'-1)*`n_obs'+1', `=`i'*`n_obs'')
			}
			else {
				mata: regpars = regpars \ J(n_obs, 1, select(st_matrix("e(b)"), st_matrix("e(b)")))
			}
		}
	}

	mata: data_boot = asarray_create("real")
	
	local i = 1
	foreach val of local values {
		cap drop `exposure'
		gen `exposure' = `val'
		mata: asarray(data_boot, `i', (st_data(., "`parlist'", "include"), J(n_reps*n_obs, 1, 1)))
		local ++i
	}

	mata: po_boot = pomeans(data_boot, regpars, J(n_reps*n_obs, 1, 1), n_reps, n_obs, link)
	mata: po_var = diagonal(variance(po_boot))
	local i = 1
	foreach val in `values' {
		mata: st_local("var`val'", strofreal(po_var[`i', 1]))
		return scalar po`val'_se3 = sqrt(`var`val'')
		local ++i
	}
	
	mata: te = po_boot[, 2] - po_boot[, 1]
	mata: ate_p25 = mm_quantile(te, 1, 0.025)
	mata: ate_p975 = mm_quantile(te, 1, 0.975)
	mata: ate_var = variance(po_boot[, 2] - po_boot[, 1])
	mata: st_local("ate_p25", strofreal(ate_p25))
	mata: st_local("ate_p975", strofreal(ate_p975))
	mata: st_local("varate", strofreal(ate_var))
	return scalar ate_se = sqrt(`varate')
	return scalar ate_p25 = `ate_p25'
	return scalar ate_p975 = `ate_p975'
	
	restore
end
	

* Mata functions
mata:

// Sampling n observations with replacement for variables in varlist among observations marked 
// by touse in the current dataset. 
matrix replsample(string scalar varlist, real scalar n) {
	data_obs = st_data(., varlist)
	index = ceil(rows(data_obs)*runiform(n, 1))
	return(data_obs[index, ])
}

// Simulating parameters from the limiting multivariate normal distribution.
matrix simpar(real vector par_est, real matrix par_var, real scalar n) {	
	return((cholesky(par_var)*invnormal(uniform(cols(par_est), n)))' + J(n, 1, par_est))
}

matrix mu2(matrix x, string scalar link) {
	if(link=="Identity") {
		mu = x
	}
	else if(link=="Log") {
		mu = exp(x)
	}
	else if(link=="Logit") {
		mu = 1:/(1:+exp(-x))
	}
	else if(link=="Probit") {
		mu = normal(x)
	}
	else if(link=="Complementary log-log") {
		mu = 1:-exp(-exp(x))
	}
	else if(link=="Neg. Binomial") {
		mu = -1:/(1:-exp(-x))
	}
	else if(link=="Log-log") {
		mu = exp(-exp(-x))
	}
	else if(link=="Log complement") {
		mu = 1:-exp(x)
	}
	return(mu)
}

// Calcuate, for each bootstrap sample, the mean of the potential outcomes. 
// The touse vector should be a vector of zeros and ones which specifies among which subjects the mean is to be calculated (i.e. among the treated).
matrix pomeans(string scalar dataarray, real matrix beta, real vector touse, real scalar n_reps, real scalar n_obs, string scalar link) {
	k = asarray_elements(dataarray)
	po = J(n_reps, k, .)
	for(i=1; i<=k; i++) {
		lin_pred = rowsum(select(beta, touse) :* select(asarray(dataarray, i), touse))
		cond_exp = mu2(lin_pred, link)
		index2 = 0
		for(j=1; j<=n_reps; j++) {
			index1 = index2 + 1
			index2 = index2 + sum(touse[(j-1)*n_obs+1::j*n_obs])
			po[j, i] = mean(cond_exp[index1::index2])
		}
	}
	return(po)
}

end

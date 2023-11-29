program define sim_cont_outcome
	clear
	set obs `1'
	gen x = rbinomial(1, 0.5)
	gen b = rnormal(`2', 1)
	gen y = `3' + `4'*x + `5'*b + `6'*x*b + rnormal(0, 1)
end

program define sim_cont_depend
	clear
	set obs `1'
	gen b = rnormal(`2', 1)
	gen c = rnormal(`3', 1)
	gen p = 1/(1+exp(-`4'*c))
	replace p = 0.99999 if p == 1
	gen x = rbinomial(1, p)
	gen y = `5' + `6'*x + `7'*b + `8'*c + `9'*x*b + `10'*x*c + `11'*b*c + `12'*x*b*c + rnormal(0, 1)
end


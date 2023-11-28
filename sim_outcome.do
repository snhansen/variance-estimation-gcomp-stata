program define sim_cont_outcome
	clear
	set obs `1'
	gen x = rbinomial(1, 0.5)
	gen c = rnormal(`2', 1)
	gen y = `3' + `4'*x + `5'*c + `6'*x*c + rnormal(0, 1)
end
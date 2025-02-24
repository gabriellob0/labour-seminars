	*1.) mylogit_mle1 - defines error variable 
	*2.) mylogit_mle2 - definites likelihood functions for error-corrected maximum likelihood logit model
	*3.) me_correction - applies ME correction procedure in order to obtain the approriate error rate for #2 and #3
	*4.) build_table4 - Saves MLE results to a matrix called "results."
 
************************************************************************
*1.) Define error rate variable
capture program drop mylogit_mle1
program define mylogit_mle1
	if `1' > 1 | `1' < 0 {
	display as error "Error term must be between 0 and 1, inclusive"
	exit
	}
	capture confirm variable error // Check if variable exists
	if !_rc {
		drop error
	}
	gen error = `1'
end
************************************************************************
*2.) Define likelihood functions for error-corrected maximum likelihood logit model (note: if error = 0, model is equivalent to standard logit)
capture program drop mylogit_mle2
program define mylogit_mle2
args lnf Xb
	replace `lnf' = ln(invlogit(`Xb')*(1-error) + (1-invlogit(`Xb'))*(error)) if $ML_y1==1 // if dep var==1
	replace `lnf' =  ln(invlogit(`Xb')*(error) + (1-invlogit(`Xb'))*(1-error))   if $ML_y1==0
end	
************************************************************************
*3.) Define program to caclulate measurement error correction

capture program drop me_correction
program define me_correction

	qui count if wagegap==-5
	if r(N)== 0 {
		preserve
		clear
		set obs 1
		gen wagegap=-5
		tempfile add
		save `add'
		restore 
		append using `add'
	}


	qui sum mm1 if tag==1 & wagegap <= -2
	if r(sd) == 0 {
		capture glm mm1 wagegap if wagegap < 0 & tag == 1, link(logit) vce(robust) iterate(200)
		if _rc != 0 qui reg mm1 wagegap if wagegap <= 2 & tag == 1, vce(robust)
		local ec_method = 2
		if _rc != 0 local ec_method = 3
		}
	if r(sd) != 0 {
		qui reg mm1 wagegap if wagegap <= -2 & tag == 1, vce(robust)
		local ec_method = 1
		}
	local ec_b	= _b[wagegap] 
	
	if _b[wagegap] <= 0 {
		predict p1 if wagegap <= -2
		qui summ p1 if wagegap==-5
		local error = 1 - r(mean)
	}
	if _b[wagegap] > 0 {
		predict p1 if wagegap <= -2
		qui sum mm1 if wagegap <= -2
		replace p1 = r(mean) if wagegap<= -2
		qui summ p1 if wagegap == -5
		local error = 1 - r(mean)
	}

	global error = .
	
	if `error' >= 0 {
		gen mmf1 = (mm1 - `error')/(1-2*`error')
		global error = `error'
	}
	else {
		gen mmf1 = mm1
		local ec_method = 0
		global error = 0
	}
		
	global ec_method = `ec_method'
	
	end

************************************************************************

*This do file loads the RN and EUT time preference ML routines

capture program drop ml_rn_discount_flex
capture program drop ml_eut_discount_flex
capture program drop ml_rdu_discount_flex
capture program drop ml_rdu_discount_flex_LN
capture program drop ml_rdu_discount_mixed
capture program drop ml_rdu_discount_mixed_qh
capture program drop ml_rdu_discount_mixed_ewb
capture program drop ml_rdu_discount_mixed_hqh
capture program drop ml_rdu_discount_mixed_hwb
capture program drop ml_rdu_discount_mixed_qhwb


*ML model for risk neutral discounting (exponential and hyperbolic)
program define ml_rn_discount_flex

* Specify the arguments of this program 
if "$discount" == "qh" | "$discount" == "weibull" {
	args lnf beta LNdelta LNnoiseDR
}
else {
	args lnf LNdelta LNnoiseDR
}

* Declare the temporary variables to be used 
tempvar choice risk ssamount ssdelay llamount lldelay delta noiseDR dSS dLL PVss PVll pvDiff    
tempvar lnftemp

quietly {

*Read in the data

generate int `choice' = $ML_y1

generate int `risk' = $ML_y2 

generate double `ssamount' = $ML_y3 
generate double `ssdelay' = $ML_y4/365 

generate double `llamount' = $ML_y5 
generate double `lldelay' = $ML_y6/365 

*Transform parameters
generate double `delta' = exp(`LNdelta')
generate double `noiseDR' = exp(`LNnoiseDR')


*Generate time preference data (i.e. present value)
if "$discount" == "exp" {
*SS
generate double `dSS' = 1 
replace `dSS' = 1/((1 + `delta')^`ssdelay')
generate double `PVss' = `dSS'*`ssamount' 

*LL
generate double `dLL' = 1 
replace `dLL' = 1/((1 + `delta')^`lldelay')
generate double `PVll' = `dLL'*`llamount'
}
else if "$discount" == "mazur" {
*SS
generate double `dSS' = 1
replace `dSS' = 1/(1 + (`delta'*`ssdelay'))
generate double `PVss' = `dSS' * `ssamount'

*LL
generate double `dLL' = 1
replace `dLL' = 1/(1 + (`delta'*`lldelay'))
generate double `PVll' = `dLL' * `llamount'
}
else if "$discount" == "qh" {
*SS
generate double `dSS' = 1
replace `dSS' = `beta'/((1 + `delta')^`ssdelay') if `ssdelay' > 0
generate double `PVss' = `dSS' * `ssamount'

*LL
generate double `dLL' = 1
replace `dLL' = `beta'/((1 + `delta')^`lldelay') if `lldelay' > 0
generate double `PVll' = `dLL' * `llamount'
}
else if "$discount" == "weibull" {
*SS
generate double `dSS' = 1 
replace `dSS' = exp(-`delta'*(`ssdelay')^(1/`beta'))
generate double `PVss' = `dSS'*`ssamount' 

*LL
generate double `dLL' = 1 
replace `dLL' = exp(-`delta'*(`lldelay')^(1/`beta'))
generate double `PVll' = `dLL'*`llamount'
}


*Calculate the PV difference using Fechner errors
generate double `pvDiff' = (`PVss' - `PVll')/`noiseDR' 

* Evaluate the likelihood for time choices
replace `lnf' = ln($cdf(`pvDiff')) if `choice'==0 & risk == 0
replace `lnf' = ln($cdf(-`pvDiff')) if `choice'==1 & risk == 0

}


end


*ML model for EUT with CU error and discounting with Fechner error
program define ml_eut_discount_flex

*Remember to define the globals $cdf, $ufunc and $discount

* Specify the arguments of this program 
if "$discount" == "qh" | "$discount" == "weibull" {
	args lnf r beta delta noiseRA noiseDR
}
else {
    args lnf r delta noiseRA noiseDR
}

* Declare the temporary variables to be used 
tempvar choice p1l p2l p3l p1r p2r p3r prize1L prize2L prize3L prize1R prize2R prize3R ///
u1L u2L u3L u1R u2R u3R euL euR euDiff uHigh uLow uMax uMin
tempvar risk ssamount ssdelay llamount lldelay uSS uLL dSS dLL PVss PVll pvDiff    

quietly {

*Read in the data

generate int `choice' = $ML_y1

generate double `p1l' = $ML_y2 
generate double `p2l' = $ML_y3 
generate double `p3l' = $ML_y4 

generate double `p1r' = $ML_y5 
generate double `p2r' = $ML_y6
generate double `p3r' = $ML_y7

generate double `prize1L' = $ML_y8 
generate double `prize2L' = $ML_y9 
generate double `prize3L' = $ML_y10 

generate double `prize1R' = $ML_y11 
generate double `prize2R' = $ML_y12 
generate double `prize3R' = $ML_y13 

generate double `uMax' = $ML_y14
generate double `uMin' = $ML_y15

generate int `risk' = $ML_y16
generate double `ssamount' = $ML_y17
generate double `ssdelay' = $ML_y18/365
generate double `llamount' = $ML_y19
generate double `lldelay' = $ML_y20/365

* Construct the argument of the utility function 
if "$ufunc" == "power" {
*Left lottery
	generate double `u1L' = `prize1L'^`r' 
	generate double `u2L' = `prize2L'^`r' 
	generate double `u3L' = `prize3L'^`r'

	replace `u1L' = ln(`prize1L') if `r' == 0
	replace `u2L' = ln(`prize2L') if `r' == 0
	replace `u3L' = ln(`prize3L') if `r' == 0

	replace `u1L' = -`prize1L'^`r' if `r' < 0
	replace `u2L' = -`prize2L'^`r' if `r' < 0
	replace `u3L' = -`prize3L'^`r' if `r' < 0 


*Right lottery
	generate double `u1R' = `prize1R'^`r' 
	generate double `u2R' = `prize2R'^`r' 
	generate double `u3R' = `prize3R'^`r'

	replace `u1R' = ln(`prize1R') if `r' == 0
	replace `u2R' = ln(`prize2R') if `r' == 0
	replace `u3R' = ln(`prize3R') if `r' == 0

	replace `u1R' = -`prize1R'^`r' if `r' < 0
	replace `u2R' = -`prize2R'^`r' if `r' < 0
	replace `u3R' = -`prize3R'^`r' if `r' < 0 


*Capture discounting information
	generate double `uSS' = `ssamount'^`r'
	generate double `uLL' = `llamount'^`r'

	replace `uSS' = ln(`ssamount') if `r' == 0
	replace `uLL' = ln(`llamount') if `r' == 0

	replace `uSS' = -`ssamount'^`r' if `r' < 0
	replace `uLL' = -`llamount'^`r' if `r' < 0
	
}
else if "$ufunc" == "crra" {
*This model converges when I remove the ln(prize) stuff
*Left lottery
	generate double `u1L' = (`prize1L'^(1-`r'))/(1-`r') 
	generate double `u2L' = (`prize2L'^(1-`r'))/(1-`r') 
	generate double `u3L' = (`prize3L'^(1-`r'))/(1-`r')

	replace `u1L' = ln(`prize1L') if `r' == 1
	replace `u2L' = ln(`prize2L') if `r' == 1
	replace `u3L' = ln(`prize3L') if `r' == 1

*Right lottery
	generate double `u1R' = (`prize1R'^(1-`r'))/(1-`r') 
	generate double `u2R' = (`prize2R'^(1-`r'))/(1-`r') 
	generate double `u3R' = (`prize3R'^(1-`r'))/(1-`r')
	
	replace `u1R' = ln(`prize1R') if `r' == 1
	replace `u2R' = ln(`prize2R') if `r' == 1
	replace `u3R' = ln(`prize3R') if `r' == 1

*Capture discounting information
	generate double `uSS' = (`ssamount'^(1-`r'))/(1-`r')
	generate double `uLL' = (`llamount'^(1-`r'))/(1-`r')

	replace `uSS' = ln(`ssamount') if `r' == 1
	replace `uLL' = ln(`llamount') if `r' == 1

}
else if "$ufunc" == "cara" {
*Left lottery
	generate double `u1L' = (1-exp(-`r'*`prize1L'))/`r'
	generate double `u2L' = (1-exp(-`r'*`prize2L'))/`r'
	generate double `u3L' = (1-exp(-`r'*`prize3L'))/`r'

*Right lottery
	generate double `u1R' = (1-exp(-`r'*`prize1R'))/`r'
	generate double `u2R' = (1-exp(-`r'*`prize2R'))/`r'
	generate double `u3R' = (1-exp(-`r'*`prize3R'))/`r'

*Capture discounting information
	generate double `uSS' = (1-exp(-`r'*`ssamount'))/`r'
	generate double `uLL' = (1-exp(-`r'*`llamount'))/`r'
}


* Calculate EU of each lottery 
generate double `euL' = (`p1l'*`u1L')+(`p2l'*`u2L')+(`p3l'*`u3L')
generate double `euR' = (`p1r'*`u1R')+(`p2r'*`u2R')+(`p3r'*`u3R')

* Calculate the EU difference by using Contextual errors 
if "$ufunc" == "power" {
	generate double `uHigh' = `uMax'^`r'
	generate double `uLow' = `uMin'^`r'
				
	replace `uHigh' = ln(`uMax') if `r' == 0
	replace `uLow' = ln(`uMin') if `r' == 0
				
	replace `uHigh' = -`uMax'^`r' if `r' < 0
	replace `uLow' = -`uMin'^`r' if `r' < 0
}
else if "$ufunc" == "crra" {
	generate double `uHigh' = (`uMax'^(1-`r'))/(1-`r')
	generate double `uLow' = (`uMin'^(1-`r'))/(1-`r')
	
	replace `uHigh' = ln(`uMax') if `r' == 1
	replace `uLow' = ln(`uMin') if `r' == 1
}
else if "$ufunc" == "cara" {
	generate double `uHigh' = (1-exp(-`r'*`uMax'))/`r'
	generate double `uLow' = (1-exp(-`r'*`uMin'))/`r'
}


generate double `euDiff' = ((`euR' - `euL')/(`uHigh' - `uLow'))/`noiseRA'

* Evaluate the likelihood for risky choices
replace `lnf' = ln($cdf( `euDiff')) if `choice'==1 & `risk' == 1 
replace `lnf' = ln($cdf(-`euDiff')) if `choice'==0 & `risk' == 1
 
*Now focus on time preference data (i.e. present value)
if "$discount" == "exp" {
*SS
generate double `dSS' = 1
replace `dSS' = 1/((1 +`delta')^`ssdelay')
generate double `PVss' = `dSS' * `uSS'

*LL
generate double `dLL' = 1
replace `dLL' = 1/((1 +`delta')^`lldelay')
generate double `PVll' = `dLL' * `uLL'
}
else if "$discount" == "mazur" {
*SS
generate double `dSS' = 1
replace `dSS' = 1/(1 + (`delta'*`ssdelay'))
generate double `PVss' = `dSS' * `uSS'

*LL
generate double `dLL' = 1
replace `dLL' = 1/(1 + (`delta'*`lldelay'))
generate double `PVll' = `dLL' * `uLL'
}
else if "$discount" == "qh" {
*SS
generate double `dSS' = 1
replace `dSS' = `beta'/((1 + `delta')^`ssdelay') if `ssdelay' > 0
generate double `PVss' = `dSS' * `uSS'

*LL
generate double `dLL' = 1
replace `dLL' = `beta'/((1 + `delta')^`lldelay') if `lldelay' > 0
generate double `PVll' = `dLL' * `uLL'
}
else if "$discount" == "weibull" {
*SS
generate double `dSS' = 1 
replace `dSS' = exp(-`delta'*(`ssdelay')^(1/`beta'))
generate double `PVss' = `dSS'*`uSS' 

*LL
generate double `dLL' = 1 
replace `dLL' = exp(-`delta'*(`lldelay')^(1/`beta'))
generate double `PVll' = `dLL'*`uLL'
}


*Calculate the PV difference using Fechner errors
generate double `pvDiff' = (`PVss' - `PVll')/`noiseDR'


*Evaluate the likelihood for time choices
replace `lnf' = ln($cdf( `pvDiff')) if `choice'==0 & `risk' == 0
replace `lnf' = ln($cdf(-`pvDiff')) if `choice'==1 & `risk' == 0

}


end




*ML model for RDU and discounting with CU error
program define ml_rdu_discount_flex

if "$discount" == "qh" | "$discount" == "weibull" {
	if "$weigh" == "prelec2" {
		args lnf r phi eta beta delta noiseRA noiseDR
	}
	else {
		args lnf r gamma beta delta noiseRA noiseDR
	}
}
else {
	if "$weigh" == "prelec2" {
	    args lnf r phi eta delta noiseRA noiseDR
	}
	else {	
    args lnf r gamma delta noiseRA noiseDR
    }
}

* Declare the temporary variables to be used 
tempvar choice p1l p2l p3l p1r p2r p3r prize1L prize2L prize3L prize1R prize2R prize3R ///
u1L u2L u3L u1R u2R u3R euL euR euDiff uHigh uLow uMax uMin
tempvar p1l p2l p3l p1r p2r p3r
tempvar pw_p1l pw_p2l pw_p3l pw_p1r pw_p2r pw_p3r
tempvar dw_p1l dw_p2l dw_p3l dw_p1r dw_p2r dw_p3r
tempvar dw_check_l dw_check_r 
tempvar risk ssamount ssdelay llamount lldelay uSS uLL dSS dLL PVss PVll pvDiff    

quietly {

*Read in the data

generate int `choice' = $ML_y1

generate double `prize1L' = $ML_y8 
generate double `prize2L' = $ML_y9 
generate double `prize3L' = $ML_y10 

generate double `prize1R' = $ML_y11 
generate double `prize2R' = $ML_y12 
generate double `prize3R' = $ML_y13 

generate double `uMax' = $ML_y14
generate double `uMin' = $ML_y15

generate int `risk' = $ML_y16
generate double `ssamount' = $ML_y17
generate double `ssdelay' = $ML_y18/365
generate double `llamount' = $ML_y19
generate double `lldelay' = $ML_y20/365


*Generate the aggregate probabilities for the Left lottery
generate double `p3l' = $ML_y4
generate double `p2l' = $ML_y3
replace `p2l' = `p2l' + `p3l' if `p2l' > 0
generate double `p1l' = $ML_y2
replace `p1l' = `p1l' + `p2l' if `p1l' > 0
replace `p1l' = `p1l' + `p3l' if `p1l' > 0 & `p2l' == 0

*Generate the weighted probabilities for the Left lottery using TK (1992), Power or Prelec weighting functions
if "$weigh" == "tk" {
	forvalues i = 1/3 {
		generate double `pw_p`i'l' = (`p`i'l'^`gamma')/(`p`i'l'^`gamma' + ((1-`p`i'l')^`gamma'))^(1/`gamma')
		*This replace ensures that cumulative probabilities of 1 are set as 1
		replace `pw_p`i'l' = 1 if `pw_p`i'l' == .
	}
}
else if "$weigh" == "power" {
	forvalues i = 1/3 {
		generate double `pw_p`i'l' = (`p`i'l'^`gamma')
	}
}
else if "$weigh" == "prelec" {
	forvalues i = 1/3 {
		generate double `pw_p`i'l' = `p`i'l'
		replace `pw_p`i'l' = exp(-((-ln(`p`i'l'))^`gamma')) if `p`i'l' > 0 & `p`i'l' < 1
	}
}
else if "$weigh" == "prelec2" {
	forvalues i = 1/3 {
		generate double `pw_p`i'l' = `p`i'l'
		replace `pw_p`i'l' = exp((-`eta')*((-ln(`p`i'l'))^`phi')) if `p`i'l' > 0 & `p`i'l' < 1
	}
}

*Generate the decision weights for the Left lottery
generate double `dw_p3l' = 0
replace 		`dw_p3l' = `pw_p3l' if `pw_p3l' > 0
generate double `dw_p2l' = 0
replace 		`dw_p2l' = `pw_p2l' if `pw_p2l' > 0
replace			`dw_p2l' = `pw_p2l' - `pw_p3l' if `dw_p3l' > 0 & `pw_p2l' > 0
generate double `dw_p1l' = 0
replace 		`dw_p1l' = `pw_p1l' if `pw_p1l' > 0
replace			`dw_p1l' = `pw_p1l' - `pw_p2l' if `dw_p2l' > 0 & `pw_p1l' > 0
replace			`dw_p1l' = `pw_p1l' - `pw_p3l' if `dw_p3l' > 0 & `pw_p2l' == 0 & `pw_p1l' > 0

*Check the decision weights
generate double `dw_check_l' = 0
replace `dw_check_l' = `dw_p3l' + `dw_p2l' + `dw_p1l'

*Now, generate the aggregate probabilities for the Right lottery
generate double `p3r' = $ML_y7
generate double `p2r' = $ML_y6
replace `p2r' = `p2r' + `p3r' if `p2r' > 0
generate double `p1r' = $ML_y5
replace `p1r' = `p1r' + `p2r' if `p1r' > 0
replace `p1r' = `p1r' + `p3r' if `p1r' > 0 & `p2r' == 0

*Generate the weighted probabilities for the Right lottery
if "$weigh" == "tk" {
	forvalues i = 1/3 {
		generate double `pw_p`i'r' = (`p`i'r'^`gamma')/(`p`i'r'^`gamma' + ((1-`p`i'r')^`gamma'))^(1/`gamma')
		*This replace ensures that cumulative probabilities of 1 are set as 1
		replace `pw_p`i'r' = 1 if `pw_p`i'r' == .
	}
}
else if "$weigh" == "power" {
	forvalues i = 1/3 {
		generate double `pw_p`i'r' = (`p`i'r'^`gamma')
	}
}
else if "$weigh" == "prelec" {
	forvalues i = 1/3 {
		generate double `pw_p`i'r' = `p`i'r'
		replace `pw_p`i'r' = exp(-((-ln(`p`i'r'))^`gamma')) if `p`i'r' > 0 & `p`i'r' < 1
	}
}
else if "$weigh" == "prelec2" {
	forvalues i = 1/3 {
		generate double `pw_p`i'r' = `p`i'r'
		replace `pw_p`i'r' = exp((-`eta')*((-ln(`p`i'r'))^`phi')) if `p`i'r' > 0 & `p`i'r' < 1
	}
}



*Generate the decision weights for the Right lottery
generate double `dw_p3r' = 0
replace 		`dw_p3r' = `pw_p3r' if `pw_p3r' > 0
generate double `dw_p2r' = 0
replace 		`dw_p2r' = `pw_p2r' if `pw_p2r' > 0
replace			`dw_p2r' = `pw_p2r' - `pw_p3r' if `dw_p3r' > 0 & `pw_p2r' > 0
generate double `dw_p1r' = 0
replace 		`dw_p1r' = `pw_p1r' if `pw_p1r' > 0
replace			`dw_p1r' = `pw_p1r' - `pw_p2r' if `dw_p2r' > 0 & `pw_p1r' > 0
replace			`dw_p1r' = `pw_p1r' - `pw_p3r' if `dw_p3r' > 0 & `pw_p2r' == 0 & `pw_p1r' > 0

*Check the decision weights
generate double `dw_check_r' = 0
replace `dw_check_r' = `dw_p3r' + `dw_p2r' + `dw_p1r'

*Check the probabilities
/* noisily summ `dw_check_l' `dw_check_r' `p1l' `p2l' `p3l' `p1r' `p2r' `p3r' ///
`pw_p1l' `pw_p2l' `pw_p3l' `pw_p1r' `pw_p2r' `pw_p3r' `dw_p1l' `dw_p2l' `dw_p3l' ///
`dw_p1r' `dw_p2r' `dw_p3r' */


*Construct the argument of the utility function 
if "$ufunc" == "power" {
*Left lottery
	generate double `u1L' = `prize1L'^`r' 
	generate double `u2L' = `prize2L'^`r' 
	generate double `u3L' = `prize3L'^`r'

	replace `u1L' = ln(`prize1L') if `r' == 0
	replace `u2L' = ln(`prize2L') if `r' == 0
	replace `u3L' = ln(`prize3L') if `r' == 0

	replace `u1L' = -`prize1L'^`r' if `r' < 0
	replace `u2L' = -`prize2L'^`r' if `r' < 0
	replace `u3L' = -`prize3L'^`r' if `r' < 0 

*Right lottery
	generate double `u1R' = `prize1R'^`r' 
	generate double `u2R' = `prize2R'^`r' 
	generate double `u3R' = `prize3R'^`r'

	replace `u1R' = ln(`prize1R') if `r' == 0
	replace `u2R' = ln(`prize2R') if `r' == 0
	replace `u3R' = ln(`prize3R') if `r' == 0

	replace `u1R' = -`prize1R'^`r' if `r' < 0
	replace `u2R' = -`prize2R'^`r' if `r' < 0
	replace `u3R' = -`prize3R'^`r' if `r' < 0 
	
*Capture discounting information
	generate double `uSS' = `ssamount'^`r'
	generate double `uLL' = `llamount'^`r'

	replace `uSS' = ln(`ssamount') if `r' == 0
	replace `uLL' = ln(`llamount') if `r' == 0

	replace `uSS' = -`ssamount'^`r' if `r' < 0
	replace `uLL' = -`llamount'^`r' if `r' < 0
	
}
else if "$ufunc" == "crra" {
*Left lottery
	generate double `u1L' = (`prize1L'^(1-`r'))/(1-`r') 
	generate double `u2L' = (`prize2L'^(1-`r'))/(1-`r')
	generate double `u3L' = (`prize3L'^(1-`r'))/(1-`r')

	replace `u1L' = ln(`prize1L') if `r' == 1
	replace `u2L' = ln(`prize2L') if `r' == 1
	replace `u3L' = ln(`prize3L') if `r' == 1


*Right lottery
	generate double `u1R' = (`prize1R'^(1-`r'))/(1-`r') 
	generate double `u2R' = (`prize2R'^(1-`r'))/(1-`r')
	generate double `u3R' = (`prize3R'^(1-`r'))/(1-`r')

	replace `u1R' = ln(`prize1R') if `r' == 1
	replace `u2R' = ln(`prize2R') if `r' == 1
	replace `u3R' = ln(`prize3R') if `r' == 1
	
*Capture discounting information
	generate double `uSS' = (`ssamount'^(1-`r'))/(1-`r')
	generate double `uLL' = (`llamount'^(1-`r'))/(1-`r')

	replace `uSS' = ln(`ssamount') if `r' == 1
	replace `uLL' = ln(`llamount') if `r' == 1

}
else if "$ufunc" == "cara" {
*Left lottery
	generate double `u1L' = (1-exp(-`r'*`prize1L'))/`r'
	generate double `u2L' = (1-exp(-`r'*`prize2L'))/`r'
	generate double `u3L' = (1-exp(-`r'*`prize3L'))/`r'

*Right lottery
	generate double `u1R' = (1-exp(-`r'*`prize1R'))/`r'
	generate double `u2R' = (1-exp(-`r'*`prize2R'))/`r'
	generate double `u3R' = (1-exp(-`r'*`prize3R'))/`r'

*Capture discounting information
	generate double `uSS' = (1-exp(-`r'*`ssamount'))/`r'
	generate double `uLL' = (1-exp(-`r'*`llamount'))/`r'

}


*Calculate EU of each lottery 
generate double `euL' = (`dw_p1l'*`u1L')+(`dw_p2l'*`u2L')+(`dw_p3l'*`u3L')
generate double `euR' = (`dw_p1r'*`u1R')+(`dw_p2r'*`u2R')+(`dw_p3r'*`u3R')

*Calculate the EU difference by using Contextual errors 
if "$ufunc" == "power" {
	generate double `uHigh' = `uMax'^`r'
	generate double `uLow' = `uMin'^`r'
				
	replace `uHigh' = ln(`uMax') if `r' == 0
	replace `uLow' = ln(`uMin') if `r' == 0
				
	replace `uHigh' = -`uMax'^`r' if `r' < 0
	replace `uLow' = -`uMin'^`r' if `r' < 0
}
else if "$ufunc" == "crra" {
	generate double `uHigh' = (`uMax'^(1-`r'))/(1-`r')
	generate double `uLow' = (`uMin'^(1-`r'))/(1-`r')
				
	replace `uHigh' = ln(`uMax') if `r' == 1
	replace `uLow' = ln(`uMin') if `r' == 1
}
else if "$ufunc" == "cara" {				
	generate double `uHigh' = (1-exp(-`r'*`uMax'))/`r' 
	generate double `uLow' = (1-exp(-`r'*`uMin'))/`r'
}

generate double `euDiff' = ((`euR' - `euL')/(`uHigh' - `uLow'))/`noiseRA'

*Evaluate the likelihood for risky choices
replace `lnf' = ln($cdf( `euDiff')) if `choice'==1 & `risk' == 1 
replace `lnf' = ln($cdf(-`euDiff')) if `choice'==0 & `risk' == 1


*Now focus on time preference data (i.e. present value)
if "$discount" == "exp" {
	*SS
	generate double `dSS' = 1
	replace `dSS' = 1/((1 +`delta')^`ssdelay')
	generate double `PVss' = `dSS' * `uSS'

	*LL
	generate double `dLL' = 1
	replace `dLL' = 1/((1 +`delta')^`lldelay')
	generate double `PVll' = `dLL' * `uLL'
}
else if "$discount" == "mazur" {
	*SS
	generate double `dSS' = 1
	replace `dSS' = 1/(1 + (`delta'*`ssdelay'))
	generate double `PVss' = `dSS' * `uSS'

	*LL
	generate double `dLL' = 1
	replace `dLL' = 1/(1 + (`delta'*`lldelay'))
	generate double `PVll' = `dLL' * `uLL'
}
else if "$discount" == "qh" {
	*SS
	generate double `dSS' = 1
	replace `dSS' = `beta'/((1 + `delta')^`ssdelay') if `ssdelay' > 0
	generate double `PVss' = `dSS' * `uSS'

	*LL
	generate double `dLL' = 1
	replace `dLL' = `beta'/((1 + `delta')^`lldelay') if `lldelay' > 0
	generate double `PVll' = `dLL' * `uLL'
}
else if "$discount" == "weibull" {
	*SS
	generate double `dSS' = 1 
	replace `dSS' = exp(-`delta'*(`ssdelay')^(1/`beta'))
	generate double `PVss' = `dSS'*`uSS' 

	*LL
	generate double `dLL' = 1 
	replace `dLL' = exp(-`delta'*(`lldelay')^(1/`beta'))
	generate double `PVll' = `dLL'*`uLL'
}


*Calculate the PV difference using Fechner errors
generate double `pvDiff' = (`PVss' - `PVll')/`noiseDR'


*Evaluate the likelihood for time choices
replace `lnf' = ln($cdf( `pvDiff')) if `choice'==0 & `risk' == 0
replace `lnf' = ln($cdf(-`pvDiff')) if `choice'==1 & `risk' == 0


}


end

*ML model for RDU and discounting with CU error
program define ml_rdu_discount_flex_LN

if "$discount" == "qh" | "$discount" == "weibull" {
	if "$weigh" == "prelec2" {
		args lnf r phi eta beta LNdelta noiseRA noiseDR
	}
	else {
		args lnf r gamma beta LNdelta noiseRA noiseDR
	}
}
else {
	if "$weigh" == "prelec2" {
	    args lnf r phi eta LNdelta noiseRA noiseDR
	}
	else {	
    args lnf r gamma LNdelta noiseRA noiseDR
    }
}

* Declare the temporary variables to be used 
tempvar choice p1l p2l p3l p1r p2r p3r prize1L prize2L prize3L prize1R prize2R prize3R ///
u1L u2L u3L u1R u2R u3R euL euR euDiff uHigh uLow uMax uMin
tempvar p1l p2l p3l p1r p2r p3r
tempvar pw_p1l pw_p2l pw_p3l pw_p1r pw_p2r pw_p3r
tempvar dw_p1l dw_p2l dw_p3l dw_p1r dw_p2r dw_p3r
tempvar dw_check_l dw_check_r 
tempvar risk ssamount ssdelay llamount lldelay uSS uLL dSS dLL PVss PVll pvDiff    
tempvar delta

quietly {

*Read in the data

generate int `choice' = $ML_y1

generate double `prize1L' = $ML_y8 
generate double `prize2L' = $ML_y9 
generate double `prize3L' = $ML_y10 

generate double `prize1R' = $ML_y11 
generate double `prize2R' = $ML_y12 
generate double `prize3R' = $ML_y13 

generate double `uMax' = $ML_y14
generate double `uMin' = $ML_y15

generate int `risk' = $ML_y16
generate double `ssamount' = $ML_y17
generate double `ssdelay' = $ML_y18/365
generate double `llamount' = $ML_y19
generate double `lldelay' = $ML_y20/365

*Impose non-negativity constraint
generate double `delta' = exp(`LNdelta')

*Generate the aggregate probabilities for the Left lottery
generate double `p3l' = $ML_y4
generate double `p2l' = $ML_y3
replace `p2l' = `p2l' + `p3l' if `p2l' > 0
generate double `p1l' = $ML_y2
replace `p1l' = `p1l' + `p2l' if `p1l' > 0
replace `p1l' = `p1l' + `p3l' if `p1l' > 0 & `p2l' == 0

*Generate the weighted probabilities for the Left lottery using TK (1992), Power or Prelec weighting functions
if "$weigh" == "tk" {
	forvalues i = 1/3 {
		generate double `pw_p`i'l' = (`p`i'l'^`gamma')/(`p`i'l'^`gamma' + ((1-`p`i'l')^`gamma'))^(1/`gamma')
		*This replace ensures that cumulative probabilities of 1 are set as 1
		replace `pw_p`i'l' = 1 if `pw_p`i'l' == .
	}
}
else if "$weigh" == "power" {
	forvalues i = 1/3 {
		generate double `pw_p`i'l' = (`p`i'l'^`gamma')
	}
}
else if "$weigh" == "prelec" {
	forvalues i = 1/3 {
		generate double `pw_p`i'l' = `p`i'l'
		replace `pw_p`i'l' = exp(-((-ln(`p`i'l'))^`gamma')) if `p`i'l' > 0 & `p`i'l' < 1
	}
}
else if "$weigh" == "prelec2" {
	forvalues i = 1/3 {
		generate double `pw_p`i'l' = `p`i'l'
		replace `pw_p`i'l' = exp((-`eta')*((-ln(`p`i'l'))^`phi')) if `p`i'l' > 0 & `p`i'l' < 1
	}
}

*Generate the decision weights for the Left lottery
generate double `dw_p3l' = 0
replace 		`dw_p3l' = `pw_p3l' if `pw_p3l' > 0
generate double `dw_p2l' = 0
replace 		`dw_p2l' = `pw_p2l' if `pw_p2l' > 0
replace			`dw_p2l' = `pw_p2l' - `pw_p3l' if `dw_p3l' > 0 & `pw_p2l' > 0
generate double `dw_p1l' = 0
replace 		`dw_p1l' = `pw_p1l' if `pw_p1l' > 0
replace			`dw_p1l' = `pw_p1l' - `pw_p2l' if `dw_p2l' > 0 & `pw_p1l' > 0
replace			`dw_p1l' = `pw_p1l' - `pw_p3l' if `dw_p3l' > 0 & `pw_p2l' == 0 & `pw_p1l' > 0

*Check the decision weights
generate double `dw_check_l' = 0
replace `dw_check_l' = `dw_p3l' + `dw_p2l' + `dw_p1l'

*Now, generate the aggregate probabilities for the Right lottery
generate double `p3r' = $ML_y7
generate double `p2r' = $ML_y6
replace `p2r' = `p2r' + `p3r' if `p2r' > 0
generate double `p1r' = $ML_y5
replace `p1r' = `p1r' + `p2r' if `p1r' > 0
replace `p1r' = `p1r' + `p3r' if `p1r' > 0 & `p2r' == 0

*Generate the weighted probabilities for the Right lottery
if "$weigh" == "tk" {
	forvalues i = 1/3 {
		generate double `pw_p`i'r' = (`p`i'r'^`gamma')/(`p`i'r'^`gamma' + ((1-`p`i'r')^`gamma'))^(1/`gamma')
		*This replace ensures that cumulative probabilities of 1 are set as 1
		replace `pw_p`i'r' = 1 if `pw_p`i'r' == .
	}
}
else if "$weigh" == "power" {
	forvalues i = 1/3 {
		generate double `pw_p`i'r' = (`p`i'r'^`gamma')
	}
}
else if "$weigh" == "prelec" {
	forvalues i = 1/3 {
		generate double `pw_p`i'r' = `p`i'r'
		replace `pw_p`i'r' = exp(-((-ln(`p`i'r'))^`gamma')) if `p`i'r' > 0 & `p`i'r' < 1
	}
}
else if "$weigh" == "prelec2" {
	forvalues i = 1/3 {
		generate double `pw_p`i'r' = `p`i'r'
		replace `pw_p`i'r' = exp((-`eta')*((-ln(`p`i'r'))^`phi')) if `p`i'r' > 0 & `p`i'r' < 1
	}
}



*Generate the decision weights for the Right lottery
generate double `dw_p3r' = 0
replace 		`dw_p3r' = `pw_p3r' if `pw_p3r' > 0
generate double `dw_p2r' = 0
replace 		`dw_p2r' = `pw_p2r' if `pw_p2r' > 0
replace			`dw_p2r' = `pw_p2r' - `pw_p3r' if `dw_p3r' > 0 & `pw_p2r' > 0
generate double `dw_p1r' = 0
replace 		`dw_p1r' = `pw_p1r' if `pw_p1r' > 0
replace			`dw_p1r' = `pw_p1r' - `pw_p2r' if `dw_p2r' > 0 & `pw_p1r' > 0
replace			`dw_p1r' = `pw_p1r' - `pw_p3r' if `dw_p3r' > 0 & `pw_p2r' == 0 & `pw_p1r' > 0

*Check the decision weights
generate double `dw_check_r' = 0
replace `dw_check_r' = `dw_p3r' + `dw_p2r' + `dw_p1r'

*Check the probabilities
/* noisily summ `dw_check_l' `dw_check_r' `p1l' `p2l' `p3l' `p1r' `p2r' `p3r' ///
`pw_p1l' `pw_p2l' `pw_p3l' `pw_p1r' `pw_p2r' `pw_p3r' `dw_p1l' `dw_p2l' `dw_p3l' ///
`dw_p1r' `dw_p2r' `dw_p3r' */


*Construct the argument of the utility function 
if "$ufunc" == "power" {
*Left lottery
	generate double `u1L' = `prize1L'^`r' 
	generate double `u2L' = `prize2L'^`r' 
	generate double `u3L' = `prize3L'^`r'

	replace `u1L' = ln(`prize1L') if `r' == 0
	replace `u2L' = ln(`prize2L') if `r' == 0
	replace `u3L' = ln(`prize3L') if `r' == 0

	replace `u1L' = -`prize1L'^`r' if `r' < 0
	replace `u2L' = -`prize2L'^`r' if `r' < 0
	replace `u3L' = -`prize3L'^`r' if `r' < 0 

*Right lottery
	generate double `u1R' = `prize1R'^`r' 
	generate double `u2R' = `prize2R'^`r' 
	generate double `u3R' = `prize3R'^`r'

	replace `u1R' = ln(`prize1R') if `r' == 0
	replace `u2R' = ln(`prize2R') if `r' == 0
	replace `u3R' = ln(`prize3R') if `r' == 0

	replace `u1R' = -`prize1R'^`r' if `r' < 0
	replace `u2R' = -`prize2R'^`r' if `r' < 0
	replace `u3R' = -`prize3R'^`r' if `r' < 0 
	
*Capture discounting information
	generate double `uSS' = `ssamount'^`r'
	generate double `uLL' = `llamount'^`r'

	replace `uSS' = ln(`ssamount') if `r' == 0
	replace `uLL' = ln(`llamount') if `r' == 0

	replace `uSS' = -`ssamount'^`r' if `r' < 0
	replace `uLL' = -`llamount'^`r' if `r' < 0
	
}
else if "$ufunc" == "crra" {
*Left lottery
	generate double `u1L' = (`prize1L'^(1-`r'))/(1-`r') 
	generate double `u2L' = (`prize2L'^(1-`r'))/(1-`r')
	generate double `u3L' = (`prize3L'^(1-`r'))/(1-`r')

	replace `u1L' = ln(`prize1L') if `r' == 1
	replace `u2L' = ln(`prize2L') if `r' == 1
	replace `u3L' = ln(`prize3L') if `r' == 1


*Right lottery
	generate double `u1R' = (`prize1R'^(1-`r'))/(1-`r') 
	generate double `u2R' = (`prize2R'^(1-`r'))/(1-`r')
	generate double `u3R' = (`prize3R'^(1-`r'))/(1-`r')

	replace `u1R' = ln(`prize1R') if `r' == 1
	replace `u2R' = ln(`prize2R') if `r' == 1
	replace `u3R' = ln(`prize3R') if `r' == 1
	
*Capture discounting information
	generate double `uSS' = (`ssamount'^(1-`r'))/(1-`r')
	generate double `uLL' = (`llamount'^(1-`r'))/(1-`r')

	replace `uSS' = ln(`ssamount') if `r' == 1
	replace `uLL' = ln(`llamount') if `r' == 1

}
else if "$ufunc" == "cara" {
*Left lottery
	generate double `u1L' = (1-exp(-`r'*`prize1L'))/`r'
	generate double `u2L' = (1-exp(-`r'*`prize2L'))/`r'
	generate double `u3L' = (1-exp(-`r'*`prize3L'))/`r'

*Right lottery
	generate double `u1R' = (1-exp(-`r'*`prize1R'))/`r'
	generate double `u2R' = (1-exp(-`r'*`prize2R'))/`r'
	generate double `u3R' = (1-exp(-`r'*`prize3R'))/`r'

*Capture discounting information
	generate double `uSS' = (1-exp(-`r'*`ssamount'))/`r'
	generate double `uLL' = (1-exp(-`r'*`llamount'))/`r'

}


*Calculate EU of each lottery 
generate double `euL' = (`dw_p1l'*`u1L')+(`dw_p2l'*`u2L')+(`dw_p3l'*`u3L')
generate double `euR' = (`dw_p1r'*`u1R')+(`dw_p2r'*`u2R')+(`dw_p3r'*`u3R')

*Calculate the EU difference by using Contextual errors 
if "$ufunc" == "power" {
	generate double `uHigh' = `uMax'^`r'
	generate double `uLow' = `uMin'^`r'
				
	replace `uHigh' = ln(`uMax') if `r' == 0
	replace `uLow' = ln(`uMin') if `r' == 0
				
	replace `uHigh' = -`uMax'^`r' if `r' < 0
	replace `uLow' = -`uMin'^`r' if `r' < 0
}
else if "$ufunc" == "crra" {
	generate double `uHigh' = (`uMax'^(1-`r'))/(1-`r')
	generate double `uLow' = (`uMin'^(1-`r'))/(1-`r')
				
	replace `uHigh' = ln(`uMax') if `r' == 1
	replace `uLow' = ln(`uMin') if `r' == 1
}
else if "$ufunc" == "cara" {				
	generate double `uHigh' = (1-exp(-`r'*`uMax'))/`r' 
	generate double `uLow' = (1-exp(-`r'*`uMin'))/`r'
}

generate double `euDiff' = ((`euR' - `euL')/(`uHigh' - `uLow'))/`noiseRA'

*Evaluate the likelihood for risky choices
replace `lnf' = ln($cdf( `euDiff')) if `choice'==1 & `risk' == 1 
replace `lnf' = ln($cdf(-`euDiff')) if `choice'==0 & `risk' == 1


*Now focus on time preference data (i.e. present value)
if "$discount" == "exp" {
	*SS
	generate double `dSS' = 1
	replace `dSS' = 1/((1 +`delta')^`ssdelay')
	generate double `PVss' = `dSS' * `uSS'

	*LL
	generate double `dLL' = 1
	replace `dLL' = 1/((1 +`delta')^`lldelay')
	generate double `PVll' = `dLL' * `uLL'
}
else if "$discount" == "mazur" {
	*SS
	generate double `dSS' = 1
	replace `dSS' = 1/(1 + (`delta'*`ssdelay'))
	generate double `PVss' = `dSS' * `uSS'

	*LL
	generate double `dLL' = 1
	replace `dLL' = 1/(1 + (`delta'*`lldelay'))
	generate double `PVll' = `dLL' * `uLL'
}
else if "$discount" == "qh" {
	*SS
	generate double `dSS' = 1
	replace `dSS' = `beta'/((1 + `delta')^`ssdelay') if `ssdelay' > 0
	generate double `PVss' = `dSS' * `uSS'

	*LL
	generate double `dLL' = 1
	replace `dLL' = `beta'/((1 + `delta')^`lldelay') if `lldelay' > 0
	generate double `PVll' = `dLL' * `uLL'
}
else if "$discount" == "weibull" {
	*SS
	generate double `dSS' = 1 
	replace `dSS' = exp(-`delta'*(`ssdelay')^(1/`beta'))
	generate double `PVss' = `dSS'*`uSS' 

	*LL
	generate double `dLL' = 1 
	replace `dLL' = exp(-`delta'*(`lldelay')^(1/`beta'))
	generate double `PVll' = `dLL'*`uLL'
}


*Calculate the PV difference using Fechner errors
generate double `pvDiff' = (`PVss' - `PVll')/`noiseDR'


*Evaluate the likelihood for time choices
replace `lnf' = ln($cdf( `pvDiff')) if `choice'==0 & `risk' == 0
replace `lnf' = ln($cdf(-`pvDiff')) if `choice'==1 & `risk' == 0


}


end


*Mixture model of Exponential and Hyperbolic Discounting (RDU utility function)
program define ml_rdu_discount_mixed

*Remember to define the globals $cdf, $ufunc

* Specify the arguments of this program 
if "$weigh" == "prelec2" {
    args lnf r phi eta deltaE deltaH noiseRA noiseDR kappa
}
else {
   args lnf r gamma deltaE deltaH noiseRA noiseDR kappa
}
* Declare the temporary variables to be used 
tempvar choice p1l p2l p3l p1r p2r p3r prize1L prize2L prize3L prize1R prize2R prize3R ///
u1L u2L u3L u1R u2R u3R euL euR euDiff uHigh uLow uMax uMin
tempvar p1l p2l p3l p1r p2r p3r
tempvar pw_p1l pw_p2l pw_p3l pw_p1r pw_p2r pw_p3r
tempvar dw_p1l dw_p2l dw_p3l dw_p1r dw_p2r dw_p3r
tempvar dw_check_l dw_check_r 
tempvar risk ssamount ssdelay llamount lldelay uSS uLL dSS dLL PVss PVll pvDiff
*Mixture model variables
tempvar lnfED lnfHD f1 f2 p1 p2    

quietly {

*Read in the data

generate int `choice' = $ML_y1

generate double `prize1L' = $ML_y8 
generate double `prize2L' = $ML_y9 
generate double `prize3L' = $ML_y10 

generate double `prize1R' = $ML_y11 
generate double `prize2R' = $ML_y12 
generate double `prize3R' = $ML_y13 

generate double `uMax' = $ML_y14
generate double `uMin' = $ML_y15

generate int `risk' = $ML_y16
generate double `ssamount' = $ML_y17
generate double `ssdelay' = $ML_y18/365
generate double `llamount' = $ML_y19
generate double `lldelay' = $ML_y20/365

*Generate the aggregate probabilities for the Left lottery
generate double `p3l' = $ML_y4
generate double `p2l' = $ML_y3
replace `p2l' = `p2l' + `p3l' if `p2l' > 0
generate double `p1l' = $ML_y2
replace `p1l' = `p1l' + `p2l' if `p1l' > 0
replace `p1l' = `p1l' + `p3l' if `p1l' > 0 & `p2l' == 0

*Generate the weighted probabilities for the Left lottery using TK (1992), Power or Prelec weighting functions
if "$weigh" == "tk" {
	forvalues i = 1/3 {
		generate double `pw_p`i'l' = (`p`i'l'^`gamma')/(`p`i'l'^`gamma' + ((1-`p`i'l')^`gamma'))^(1/`gamma')
		*This replace ensures that cumulative probabilities of 1 are set as 1
		replace `pw_p`i'l' = 1 if `pw_p`i'l' == .
	}
}
else if "$weigh" == "power" {
	forvalues i = 1/3 {
		generate double `pw_p`i'l' = (`p`i'l'^`gamma')
	}
}
else if "$weigh" == "prelec" {
	forvalues i = 1/3 {
		generate double `pw_p`i'l' = `p`i'l'
		replace `pw_p`i'l' = exp(-((-ln(`p`i'l'))^`gamma')) if `p`i'l' > 0 & `p`i'l' < 1
	}
}
else if "$weigh" == "prelec2" {
	forvalues i = 1/3 {
		generate double `pw_p`i'l' = `p`i'l'
		replace `pw_p`i'l' = exp((-`eta')*((-ln(`p`i'l'))^`phi')) if `p`i'l' > 0 & `p`i'l' < 1
	}
}


*Generate the decision weights for the Left lottery
generate double `dw_p3l' = 0
replace 		`dw_p3l' = `pw_p3l' if `pw_p3l' > 0
generate double `dw_p2l' = 0
replace 		`dw_p2l' = `pw_p2l' if `pw_p2l' > 0
replace			`dw_p2l' = `pw_p2l' - `pw_p3l' if `dw_p3l' > 0 & `pw_p2l' > 0
generate double `dw_p1l' = 0
replace 		`dw_p1l' = `pw_p1l' if `pw_p1l' > 0
replace			`dw_p1l' = `pw_p1l' - `pw_p2l' if `dw_p2l' > 0 & `pw_p1l' > 0
replace			`dw_p1l' = `pw_p1l' - `pw_p3l' if `dw_p3l' > 0 & `pw_p2l' == 0 & `pw_p1l' > 0

*Check the decision weights
generate double `dw_check_l' = 0
replace `dw_check_l' = `dw_p3l' + `dw_p2l' + `dw_p1l'

*Now, generate the aggregate probabilities for the Right lottery
generate double `p3r' = $ML_y7
generate double `p2r' = $ML_y6
replace `p2r' = `p2r' + `p3r' if `p2r' > 0
generate double `p1r' = $ML_y5
replace `p1r' = `p1r' + `p2r' if `p1r' > 0
replace `p1r' = `p1r' + `p3r' if `p1r' > 0 & `p2r' == 0

*Generate the weighted probabilities for the Right lottery
if "$weigh" == "tk" {
	forvalues i = 1/3 {
		generate double `pw_p`i'r' = (`p`i'r'^`gamma')/(`p`i'r'^`gamma' + ((1-`p`i'r')^`gamma'))^(1/`gamma')
		*This replace ensures that cumulative probabilities of 1 are set as 1
		replace `pw_p`i'r' = 1 if `pw_p`i'r' == .
	}
}
else if "$weigh" == "power" {
	forvalues i = 1/3 {
		generate double `pw_p`i'r' = (`p`i'r'^`gamma')
	}
}
else if "$weigh" == "prelec" {
	forvalues i = 1/3 {
		generate double `pw_p`i'r' = `p`i'r'
		replace `pw_p`i'r' = exp(-((-ln(`p`i'r'))^`gamma')) if `p`i'r' > 0 & `p`i'r' < 1
	}
}
else if "$weigh" == "prelec2" {
	forvalues i = 1/3 {
		generate double `pw_p`i'r' = `p`i'r'
		replace `pw_p`i'r' = exp((-`eta')*((-ln(`p`i'r'))^`phi')) if `p`i'r' > 0 & `p`i'r' < 1
	}
}


*Generate the decision weights for the Right lottery
generate double `dw_p3r' = 0
replace 		`dw_p3r' = `pw_p3r' if `pw_p3r' > 0
generate double `dw_p2r' = 0
replace 		`dw_p2r' = `pw_p2r' if `pw_p2r' > 0
replace			`dw_p2r' = `pw_p2r' - `pw_p3r' if `dw_p3r' > 0 & `pw_p2r' > 0
generate double `dw_p1r' = 0
replace 		`dw_p1r' = `pw_p1r' if `pw_p1r' > 0
replace			`dw_p1r' = `pw_p1r' - `pw_p2r' if `dw_p2r' > 0 & `pw_p1r' > 0
replace			`dw_p1r' = `pw_p1r' - `pw_p3r' if `dw_p3r' > 0 & `pw_p2r' == 0 & `pw_p1r' > 0

*Check the decision weights
generate double `dw_check_r' = 0
replace `dw_check_r' = `dw_p3r' + `dw_p2r' + `dw_p1r'

*Check the probabilities
/* noisily summ `dw_check_l' `dw_check_r' `p1l' `p2l' `p3l' `p1r' `p2r' `p3r' ///
`pw_p1l' `pw_p2l' `pw_p3l' `pw_p1r' `pw_p2r' `pw_p3r' `dw_p1l' `dw_p2l' `dw_p3l' ///
`dw_p1r' `dw_p2r' `dw_p3r' */


*Construct the argument of the utility function 
if "$ufunc" == "power" {
*Left lottery
	generate double `u1L' = `prize1L'^`r' 
	generate double `u2L' = `prize2L'^`r' 
	generate double `u3L' = `prize3L'^`r'

	replace `u1L' = ln(`prize1L') if `r' == 0
	replace `u2L' = ln(`prize2L') if `r' == 0
	replace `u3L' = ln(`prize3L') if `r' == 0

	replace `u1L' = -`prize1L'^`r' if `r' < 0
	replace `u2L' = -`prize2L'^`r' if `r' < 0
	replace `u3L' = -`prize3L'^`r' if `r' < 0 

*Right lottery
	generate double `u1R' = `prize1R'^`r' 
	generate double `u2R' = `prize2R'^`r' 
	generate double `u3R' = `prize3R'^`r'

	replace `u1R' = ln(`prize1R') if `r' == 0
	replace `u2R' = ln(`prize2R') if `r' == 0
	replace `u3R' = ln(`prize3R') if `r' == 0

	replace `u1R' = -`prize1R'^`r' if `r' < 0
	replace `u2R' = -`prize2R'^`r' if `r' < 0
	replace `u3R' = -`prize3R'^`r' if `r' < 0 
	
*Capture discounting information
	generate double `uSS' = `ssamount'^`r'
	generate double `uLL' = `llamount'^`r'

	replace `uSS' = ln(`ssamount') if `r' == 0
	replace `uLL' = ln(`llamount') if `r' == 0

	replace `uSS' = -`ssamount'^`r' if `r' < 0
	replace `uLL' = -`llamount'^`r' if `r' < 0
	
}
else if "$ufunc" == "crra" {
*Left lottery
	generate double `u1L' = (`prize1L'^(1-`r'))/(1-`r') 
	generate double `u2L' = (`prize2L'^(1-`r'))/(1-`r')
	generate double `u3L' = (`prize3L'^(1-`r'))/(1-`r')

	replace `u1L' = ln(`prize1L') if `r' == 1
	replace `u2L' = ln(`prize2L') if `r' == 1
	replace `u3L' = ln(`prize3L') if `r' == 1


*Right lottery
	generate double `u1R' = (`prize1R'^(1-`r'))/(1-`r') 
	generate double `u2R' = (`prize2R'^(1-`r'))/(1-`r')
	generate double `u3R' = (`prize3R'^(1-`r'))/(1-`r')

	replace `u1R' = ln(`prize1R') if `r' == 1
	replace `u2R' = ln(`prize2R') if `r' == 1
	replace `u3R' = ln(`prize3R') if `r' == 1
	
*Capture discounting information
	generate double `uSS' = (`ssamount'^(1-`r'))/(1-`r')
	generate double `uLL' = (`llamount'^(1-`r'))/(1-`r')

	replace `uSS' = ln(`ssamount') if `r' == 1
	replace `uLL' = ln(`llamount') if `r' == 1

}
else if "$ufunc" == "cara" {
*Left lottery
	generate double `u1L' = (1-exp(-`r'*`prize1L'))/`r'
	generate double `u2L' = (1-exp(-`r'*`prize2L'))/`r'
	generate double `u3L' = (1-exp(-`r'*`prize3L'))/`r'

*Right lottery
	generate double `u1R' = (1-exp(-`r'*`prize1R'))/`r'
	generate double `u2R' = (1-exp(-`r'*`prize2R'))/`r'
	generate double `u3R' = (1-exp(-`r'*`prize3R'))/`r'

*Capture discounting information
	generate double `uSS' = (1-exp(-`r'*`ssamount'))/`r'
	generate double `uLL' = (1-exp(-`r'*`llamount'))/`r'

}

*Calculate EU of each lottery 
generate double `euL' = (`dw_p1l'*`u1L')+(`dw_p2l'*`u2L')+(`dw_p3l'*`u3L')
generate double `euR' = (`dw_p1r'*`u1R')+(`dw_p2r'*`u2R')+(`dw_p3r'*`u3R')

*Calculate the EU difference by using Contextual errors 
if "$ufunc" == "power" {
	generate double `uHigh' = `uMax'^`r'
	generate double `uLow' = `uMin'^`r'
				
	replace `uHigh' = ln(`uMax') if `r' == 0
	replace `uLow' = ln(`uMin') if `r' == 0
				
	replace `uHigh' = -`uMax'^`r' if `r' < 0
	replace `uLow' = -`uMin'^`r' if `r' < 0
}
else if "$ufunc" == "crra" {
	generate double `uHigh' = (`uMax'^(1-`r'))/(1-`r')
	generate double `uLow' = (`uMin'^(1-`r'))/(1-`r')
				
	replace `uHigh' = ln(`uMax') if `r' == 1
	replace `uLow' = ln(`uMin') if `r' == 1
}
else if "$ufunc" == "cara" {				
	generate double `uHigh' = (1-exp(-`r'*`uMax'))/`r' 
	generate double `uLow' = (1-exp(-`r'*`uMin'))/`r'
}

generate double `euDiff' = ((`euR' - `euL')/(`uHigh' - `uLow'))/`noiseRA'

*Evaluate the likelihood for risky choices
replace `lnf' = ln($cdf( `euDiff')) if `choice'==1 & `risk' == 1 
replace `lnf' = ln($cdf(-`euDiff')) if `choice'==0 & `risk' == 1

 
*Now focus on time preference data (i.e. present value)
*Exponential
*SS
generate double `dSS' = 1
replace `dSS' = 1/((1 +`deltaE')^`ssdelay')
generate double `PVss' = `dSS' * `uSS'

*LL
generate double `dLL' = 1
replace `dLL' = 1/((1 +`deltaE')^`lldelay')
generate double `PVll' = `dLL' * `uLL'

*Calculate the PV difference using Fechner errors
generate double `pvDiff' = (`PVss' - `PVll')/`noiseDR' if `risk' == 0

*Evaluate the likelihood for time choices (Exponential)
generate double `lnfED' = ln($cdf( `pvDiff')) if `choice'==0 & `risk' == 0 
replace `lnfED' = ln($cdf(-`pvDiff')) if `choice'==1 & `risk' == 0



*Hyperbolic
*SS
replace `dSS' = 1
replace `dSS' = 1/(1 + (`deltaH'*`ssdelay'))
replace `PVss' = `dSS' * `uSS'

*LL
replace `dLL' = 1
replace `dLL' = 1/(1 + (`deltaH'*`lldelay'))
replace `PVll' = `dLL' * `uLL'

*Calculate the PV difference using Fechner errors
replace `pvDiff' = (`PVss' - `PVll')/`noiseDR' if `risk' == 0

*Evaluate the likelihood for time choices (Hyperbolic)
generate double `lnfHD' = ln($cdf( `pvDiff')) if `choice'==0 & `risk' == 0 
replace `lnfHD' = ln($cdf(-`pvDiff')) if `choice'==1 & `risk' == 0

*Calculate the grand likelihood for the discounting choices
generate double `f1' = exp(`lnfED') if `risk' == 0
generate double `f2' = exp(`lnfHD') if `risk' == 0
generate double `p1' = 1/(1+exp(`kappa')) if `risk' == 0
generate double `p2' = exp(`kappa')/(1+exp(`kappa')) if `risk' == 0
replace	`lnf' = ln((`p1' * `f1') + (`p2' * `f2')) if `risk' == 0

}
end


*Mixture model of Exponential and Quasi-Hyperbolic Discounting (RDU utility function)
program define ml_rdu_discount_mixed_qh

*Remember to define the globals $cdf, $ufunc

* Specify the arguments of this program 
if "$weigh" == "prelec2" {
    args lnf r phi eta deltaE betaQH deltaQH noiseRA noiseDR kappa
}
else {
    args lnf r gamma deltaE betaQH deltaQH noiseRA noiseDR kappa
}

* Declare the temporary variables to be used 
tempvar choice p1l p2l p3l p1r p2r p3r prize1L prize2L prize3L prize1R prize2R prize3R ///
u1L u2L u3L u1R u2R u3R euL euR euDiff uHigh uLow uMax uMin
tempvar p1l p2l p3l p1r p2r p3r
tempvar pw_p1l pw_p2l pw_p3l pw_p1r pw_p2r pw_p3r
tempvar dw_p1l dw_p2l dw_p3l dw_p1r dw_p2r dw_p3r
tempvar dw_check_l dw_check_r 
tempvar risk ssamount ssdelay llamount lldelay uSS uLL dSS dLL PVss PVll pvDiff
*Mixture model variables
tempvar lnfED lnfQHD f1 f2 p1 p2    

quietly {

*Read in the data

generate int `choice' = $ML_y1

generate double `prize1L' = $ML_y8 
generate double `prize2L' = $ML_y9 
generate double `prize3L' = $ML_y10 

generate double `prize1R' = $ML_y11 
generate double `prize2R' = $ML_y12 
generate double `prize3R' = $ML_y13 

generate double `uMax' = $ML_y14
generate double `uMin' = $ML_y15

generate int `risk' = $ML_y16
generate double `ssamount' = $ML_y17
generate double `ssdelay' = $ML_y18/365
generate double `llamount' = $ML_y19
generate double `lldelay' = $ML_y20/365

*Generate the aggregate probabilities for the Left lottery
generate double `p3l' = $ML_y4
generate double `p2l' = $ML_y3
replace `p2l' = `p2l' + `p3l' if `p2l' > 0
generate double `p1l' = $ML_y2
replace `p1l' = `p1l' + `p2l' if `p1l' > 0
replace `p1l' = `p1l' + `p3l' if `p1l' > 0 & `p2l' == 0

*Generate the weighted probabilities for the Left lottery using TK (1992), Power or Prelec weighting functions
if "$weigh" == "tk" {
	forvalues i = 1/3 {
		generate double `pw_p`i'l' = (`p`i'l'^`gamma')/(`p`i'l'^`gamma' + ((1-`p`i'l')^`gamma'))^(1/`gamma')
		*This replace ensures that cumulative probabilities of 1 are set as 1
		replace `pw_p`i'l' = 1 if `pw_p`i'l' == .
	}
}
else if "$weigh" == "power" {
	forvalues i = 1/3 {
		generate double `pw_p`i'l' = (`p`i'l'^`gamma')
	}
}
else if "$weigh" == "prelec" {
	forvalues i = 1/3 {
		generate double `pw_p`i'l' = `p`i'l'
		replace `pw_p`i'l' = exp(-((-ln(`p`i'l'))^`gamma')) if `p`i'l' > 0 & `p`i'l' < 1
	}
}
else if "$weigh" == "prelec2" {
	forvalues i = 1/3 {
		generate double `pw_p`i'l' = `p`i'l'
		replace `pw_p`i'l' = exp((-`eta')*((-ln(`p`i'l'))^`phi')) if `p`i'l' > 0 & `p`i'l' < 1
	}
}


*Generate the decision weights for the Left lottery
generate double `dw_p3l' = 0
replace 		`dw_p3l' = `pw_p3l' if `pw_p3l' > 0
generate double `dw_p2l' = 0
replace 		`dw_p2l' = `pw_p2l' if `pw_p2l' > 0
replace			`dw_p2l' = `pw_p2l' - `pw_p3l' if `dw_p3l' > 0 & `pw_p2l' > 0
generate double `dw_p1l' = 0
replace 		`dw_p1l' = `pw_p1l' if `pw_p1l' > 0
replace			`dw_p1l' = `pw_p1l' - `pw_p2l' if `dw_p2l' > 0 & `pw_p1l' > 0
replace			`dw_p1l' = `pw_p1l' - `pw_p3l' if `dw_p3l' > 0 & `pw_p2l' == 0 & `pw_p1l' > 0

*Check the decision weights
generate double `dw_check_l' = 0
replace `dw_check_l' = `dw_p3l' + `dw_p2l' + `dw_p1l'

*Now, generate the aggregate probabilities for the Right lottery
generate double `p3r' = $ML_y7
generate double `p2r' = $ML_y6
replace `p2r' = `p2r' + `p3r' if `p2r' > 0
generate double `p1r' = $ML_y5
replace `p1r' = `p1r' + `p2r' if `p1r' > 0
replace `p1r' = `p1r' + `p3r' if `p1r' > 0 & `p2r' == 0

*Generate the weighted probabilities for the Right lottery
if "$weigh" == "tk" {
	forvalues i = 1/3 {
		generate double `pw_p`i'r' = (`p`i'r'^`gamma')/(`p`i'r'^`gamma' + ((1-`p`i'r')^`gamma'))^(1/`gamma')
		*This replace ensures that cumulative probabilities of 1 are set as 1
		replace `pw_p`i'r' = 1 if `pw_p`i'r' == .
	}
}
else if "$weigh" == "power" {
	forvalues i = 1/3 {
		generate double `pw_p`i'r' = (`p`i'r'^`gamma')
	}
}
else if "$weigh" == "prelec" {
	forvalues i = 1/3 {
		generate double `pw_p`i'r' = `p`i'r'
		replace `pw_p`i'r' = exp(-((-ln(`p`i'r'))^`gamma')) if `p`i'r' > 0 & `p`i'r' < 1
	}
}
else if "$weigh" == "prelec2" {
	forvalues i = 1/3 {
		generate double `pw_p`i'r' = `p`i'r'
		replace `pw_p`i'r' = exp((-`eta')*((-ln(`p`i'r'))^`phi')) if `p`i'r' > 0 & `p`i'r' < 1
	}
}



*Generate the decision weights for the Right lottery
generate double `dw_p3r' = 0
replace 		`dw_p3r' = `pw_p3r' if `pw_p3r' > 0
generate double `dw_p2r' = 0
replace 		`dw_p2r' = `pw_p2r' if `pw_p2r' > 0
replace			`dw_p2r' = `pw_p2r' - `pw_p3r' if `dw_p3r' > 0 & `pw_p2r' > 0
generate double `dw_p1r' = 0
replace 		`dw_p1r' = `pw_p1r' if `pw_p1r' > 0
replace			`dw_p1r' = `pw_p1r' - `pw_p2r' if `dw_p2r' > 0 & `pw_p1r' > 0
replace			`dw_p1r' = `pw_p1r' - `pw_p3r' if `dw_p3r' > 0 & `pw_p2r' == 0 & `pw_p1r' > 0

*Check the decision weights
generate double `dw_check_r' = 0
replace `dw_check_r' = `dw_p3r' + `dw_p2r' + `dw_p1r'

*Check the probabilities
/* noisily summ `dw_check_l' `dw_check_r' `p1l' `p2l' `p3l' `p1r' `p2r' `p3r' ///
`pw_p1l' `pw_p2l' `pw_p3l' `pw_p1r' `pw_p2r' `pw_p3r' `dw_p1l' `dw_p2l' `dw_p3l' ///
`dw_p1r' `dw_p2r' `dw_p3r' */


*Construct the argument of the utility function 
if "$ufunc" == "power" {
*Left lottery
	generate double `u1L' = `prize1L'^`r' 
	generate double `u2L' = `prize2L'^`r' 
	generate double `u3L' = `prize3L'^`r'

	replace `u1L' = ln(`prize1L') if `r' == 0
	replace `u2L' = ln(`prize2L') if `r' == 0
	replace `u3L' = ln(`prize3L') if `r' == 0

	replace `u1L' = -`prize1L'^`r' if `r' < 0
	replace `u2L' = -`prize2L'^`r' if `r' < 0
	replace `u3L' = -`prize3L'^`r' if `r' < 0 

*Right lottery
	generate double `u1R' = `prize1R'^`r' 
	generate double `u2R' = `prize2R'^`r' 
	generate double `u3R' = `prize3R'^`r'

	replace `u1R' = ln(`prize1R') if `r' == 0
	replace `u2R' = ln(`prize2R') if `r' == 0
	replace `u3R' = ln(`prize3R') if `r' == 0

	replace `u1R' = -`prize1R'^`r' if `r' < 0
	replace `u2R' = -`prize2R'^`r' if `r' < 0
	replace `u3R' = -`prize3R'^`r' if `r' < 0 
	
*Capture discounting information
	generate double `uSS' = `ssamount'^`r'
	generate double `uLL' = `llamount'^`r'

	replace `uSS' = ln(`ssamount') if `r' == 0
	replace `uLL' = ln(`llamount') if `r' == 0

	replace `uSS' = -`ssamount'^`r' if `r' < 0
	replace `uLL' = -`llamount'^`r' if `r' < 0
	
}
else if "$ufunc" == "crra" {
*Left lottery
	generate double `u1L' = (`prize1L'^(1-`r'))/(1-`r') 
	generate double `u2L' = (`prize2L'^(1-`r'))/(1-`r')
	generate double `u3L' = (`prize3L'^(1-`r'))/(1-`r')

	replace `u1L' = ln(`prize1L') if `r' == 1
	replace `u2L' = ln(`prize2L') if `r' == 1
	replace `u3L' = ln(`prize3L') if `r' == 1


*Right lottery
	generate double `u1R' = (`prize1R'^(1-`r'))/(1-`r') 
	generate double `u2R' = (`prize2R'^(1-`r'))/(1-`r')
	generate double `u3R' = (`prize3R'^(1-`r'))/(1-`r')

	replace `u1R' = ln(`prize1R') if `r' == 1
	replace `u2R' = ln(`prize2R') if `r' == 1
	replace `u3R' = ln(`prize3R') if `r' == 1
	
*Capture discounting information
	generate double `uSS' = (`ssamount'^(1-`r'))/(1-`r')
	generate double `uLL' = (`llamount'^(1-`r'))/(1-`r')

	replace `uSS' = ln(`ssamount') if `r' == 1
	replace `uLL' = ln(`llamount') if `r' == 1

}
else if "$ufunc" == "cara" {
*Left lottery
	generate double `u1L' = (1-exp(-`r'*`prize1L'))/`r'
	generate double `u2L' = (1-exp(-`r'*`prize2L'))/`r'
	generate double `u3L' = (1-exp(-`r'*`prize3L'))/`r'

*Right lottery
	generate double `u1R' = (1-exp(-`r'*`prize1R'))/`r'
	generate double `u2R' = (1-exp(-`r'*`prize2R'))/`r'
	generate double `u3R' = (1-exp(-`r'*`prize3R'))/`r'

*Capture discounting information
	generate double `uSS' = (1-exp(-`r'*`ssamount'))/`r'
	generate double `uLL' = (1-exp(-`r'*`llamount'))/`r'

}

*Calculate EU of each lottery 
generate double `euL' = (`dw_p1l'*`u1L')+(`dw_p2l'*`u2L')+(`dw_p3l'*`u3L')
generate double `euR' = (`dw_p1r'*`u1R')+(`dw_p2r'*`u2R')+(`dw_p3r'*`u3R')

*Calculate the EU difference by using Contextual errors 
if "$ufunc" == "power" {
	generate double `uHigh' = `uMax'^`r'
	generate double `uLow' = `uMin'^`r'
				
	replace `uHigh' = ln(`uMax') if `r' == 0
	replace `uLow' = ln(`uMin') if `r' == 0
				
	replace `uHigh' = -`uMax'^`r' if `r' < 0
	replace `uLow' = -`uMin'^`r' if `r' < 0
}
else if "$ufunc" == "crra" {
	generate double `uHigh' = (`uMax'^(1-`r'))/(1-`r')
	generate double `uLow' = (`uMin'^(1-`r'))/(1-`r')
				
	replace `uHigh' = ln(`uMax') if `r' == 1
	replace `uLow' = ln(`uMin') if `r' == 1
}
else if "$ufunc" == "cara" {				
	generate double `uHigh' = (1-exp(-`r'*`uMax'))/`r' 
	generate double `uLow' = (1-exp(-`r'*`uMin'))/`r'
}

generate double `euDiff' = ((`euR' - `euL')/(`uHigh' - `uLow'))/`noiseRA'

*Evaluate the likelihood for risky choices
replace `lnf' = ln($cdf( `euDiff')) if `choice'==1 & `risk' == 1 
replace `lnf' = ln($cdf(-`euDiff')) if `choice'==0 & `risk' == 1

 
*Now focus on time preference data (i.e. present value)
*Exponential
*SS
generate double `dSS' = 1
replace `dSS' = 1/((1 +`deltaE')^`ssdelay')
generate double `PVss' = `dSS' * `uSS'

*LL
generate double `dLL' = 1
replace `dLL' = 1/((1 +`deltaE')^`lldelay')
generate double `PVll' = `dLL' * `uLL'

*Calculate the PV difference using Fechner errors
generate double `pvDiff' = (`PVss' - `PVll')/`noiseDR' if `risk' == 0

*Evaluate the likelihood for time choices (Exponential)
generate double `lnfED' = ln($cdf( `pvDiff')) if `choice'==0 & `risk' == 0 
replace `lnfED' = ln($cdf(-`pvDiff')) if `choice'==1 & `risk' == 0



*Quasi-Hyperbolic
*SS
replace `dSS' = 1
replace `dSS' = `betaQH'/((1 + `deltaQH')^`ssdelay') if `ssdelay' > 0
replace `PVss' = `dSS' * `uSS'

*LL
replace `dLL' = 1
replace `dLL' = `betaQH'/((1 + `deltaQH')^`lldelay') if `lldelay' > 0
replace `PVll' = `dLL' * `uLL'


*Calculate the PV difference using Fechner errors
replace `pvDiff' = (`PVss' - `PVll')/`noiseDR' if `risk' == 0

*Evaluate the likelihood for time choices (Quasi-Hyperbolic)
generate double `lnfQHD' = ln($cdf( `pvDiff')) if `choice'==0 & `risk' == 0 
replace `lnfQHD' = ln($cdf(-`pvDiff')) if `choice'==1 & `risk' == 0

*Calculate the grand likelihood for the discounting choices
generate double `f1' = exp(`lnfED') if `risk' == 0
generate double `f2' = exp(`lnfQHD') if `risk' == 0
generate double `p1' = 1/(1+exp(`kappa')) if `risk' == 0
generate double `p2' = exp(`kappa')/(1+exp(`kappa')) if `risk' == 0
replace	`lnf' = ln((`p1' * `f1') + (`p2' * `f2')) if `risk' == 0

}
end

*Mixture model of Exponential and Weibull Discounting (RDU utility function)
program define ml_rdu_discount_mixed_ewb

*Remember to define the globals $cdf, $ufunc

* Specify the arguments of this program 
if "$weigh" == "prelec2" {
    args lnf r phi eta deltaE deltaWB betaWB noiseRA noiseDR kappa
}
else {
    args lnf r gamma deltaE deltaWB betaWB noiseRA noiseDR kappa
}

* Declare the temporary variables to be used 
tempvar choice p1l p2l p3l p1r p2r p3r prize1L prize2L prize3L prize1R prize2R prize3R ///
u1L u2L u3L u1R u2R u3R euL euR euDiff uHigh uLow uMax uMin
tempvar p1l p2l p3l p1r p2r p3r
tempvar pw_p1l pw_p2l pw_p3l pw_p1r pw_p2r pw_p3r
tempvar dw_p1l dw_p2l dw_p3l dw_p1r dw_p2r dw_p3r
tempvar dw_check_l dw_check_r 
tempvar risk ssamount ssdelay llamount lldelay uSS uLL dSS dLL PVss PVll pvDiff
*Mixture model variables
tempvar lnfED lnfWB f1 f2 p1 p2    

quietly {

*Read in the data

generate int `choice' = $ML_y1

generate double `prize1L' = $ML_y8 
generate double `prize2L' = $ML_y9 
generate double `prize3L' = $ML_y10 

generate double `prize1R' = $ML_y11 
generate double `prize2R' = $ML_y12 
generate double `prize3R' = $ML_y13 

generate double `uMax' = $ML_y14
generate double `uMin' = $ML_y15

generate int `risk' = $ML_y16
generate double `ssamount' = $ML_y17
generate double `ssdelay' = $ML_y18/365
generate double `llamount' = $ML_y19
generate double `lldelay' = $ML_y20/365

*Generate the aggregate probabilities for the Left lottery
generate double `p3l' = $ML_y4
generate double `p2l' = $ML_y3
replace `p2l' = `p2l' + `p3l' if `p2l' > 0
generate double `p1l' = $ML_y2
replace `p1l' = `p1l' + `p2l' if `p1l' > 0
replace `p1l' = `p1l' + `p3l' if `p1l' > 0 & `p2l' == 0

*Generate the weighted probabilities for the Left lottery using TK (1992), Power or Prelec weighting functions
if "$weigh" == "tk" {
	forvalues i = 1/3 {
		generate double `pw_p`i'l' = (`p`i'l'^`gamma')/(`p`i'l'^`gamma' + ((1-`p`i'l')^`gamma'))^(1/`gamma')
		*This replace ensures that cumulative probabilities of 1 are set as 1
		replace `pw_p`i'l' = 1 if `pw_p`i'l' == .
	}
}
else if "$weigh" == "power" {
	forvalues i = 1/3 {
		generate double `pw_p`i'l' = (`p`i'l'^`gamma')
	}
}
else if "$weigh" == "prelec" {
	forvalues i = 1/3 {
		generate double `pw_p`i'l' = `p`i'l'
		replace `pw_p`i'l' = exp(-((-ln(`p`i'l'))^`gamma')) if `p`i'l' > 0 & `p`i'l' < 1
	}
}
else if "$weigh" == "prelec2" {
	forvalues i = 1/3 {
		generate double `pw_p`i'l' = `p`i'l'
		replace `pw_p`i'l' = exp((-`eta')*((-ln(`p`i'l'))^`phi')) if `p`i'l' > 0 & `p`i'l' < 1
	}
}


*Generate the decision weights for the Left lottery
generate double `dw_p3l' = 0
replace 		`dw_p3l' = `pw_p3l' if `pw_p3l' > 0
generate double `dw_p2l' = 0
replace 		`dw_p2l' = `pw_p2l' if `pw_p2l' > 0
replace			`dw_p2l' = `pw_p2l' - `pw_p3l' if `dw_p3l' > 0 & `pw_p2l' > 0
generate double `dw_p1l' = 0
replace 		`dw_p1l' = `pw_p1l' if `pw_p1l' > 0
replace			`dw_p1l' = `pw_p1l' - `pw_p2l' if `dw_p2l' > 0 & `pw_p1l' > 0
replace			`dw_p1l' = `pw_p1l' - `pw_p3l' if `dw_p3l' > 0 & `pw_p2l' == 0 & `pw_p1l' > 0

*Check the decision weights
generate double `dw_check_l' = 0
replace `dw_check_l' = `dw_p3l' + `dw_p2l' + `dw_p1l'

*Now, generate the aggregate probabilities for the Right lottery
generate double `p3r' = $ML_y7
generate double `p2r' = $ML_y6
replace `p2r' = `p2r' + `p3r' if `p2r' > 0
generate double `p1r' = $ML_y5
replace `p1r' = `p1r' + `p2r' if `p1r' > 0
replace `p1r' = `p1r' + `p3r' if `p1r' > 0 & `p2r' == 0

*Generate the weighted probabilities for the Right lottery
if "$weigh" == "tk" {
	forvalues i = 1/3 {
		generate double `pw_p`i'r' = (`p`i'r'^`gamma')/(`p`i'r'^`gamma' + ((1-`p`i'r')^`gamma'))^(1/`gamma')
		*This replace ensures that cumulative probabilities of 1 are set as 1
		replace `pw_p`i'r' = 1 if `pw_p`i'r' == .
	}
}
else if "$weigh" == "power" {
	forvalues i = 1/3 {
		generate double `pw_p`i'r' = (`p`i'r'^`gamma')
	}
}
else if "$weigh" == "prelec" {
	forvalues i = 1/3 {
		generate double `pw_p`i'r' = `p`i'r'
		replace `pw_p`i'r' = exp(-((-ln(`p`i'r'))^`gamma')) if `p`i'r' > 0 & `p`i'r' < 1
	}
}
else if "$weigh" == "prelec2" {
	forvalues i = 1/3 {
		generate double `pw_p`i'r' = `p`i'r'
		replace `pw_p`i'r' = exp((-`eta')*((-ln(`p`i'r'))^`phi')) if `p`i'r' > 0 & `p`i'r' < 1
	}
}


*Generate the decision weights for the Right lottery
generate double `dw_p3r' = 0
replace 		`dw_p3r' = `pw_p3r' if `pw_p3r' > 0
generate double `dw_p2r' = 0
replace 		`dw_p2r' = `pw_p2r' if `pw_p2r' > 0
replace			`dw_p2r' = `pw_p2r' - `pw_p3r' if `dw_p3r' > 0 & `pw_p2r' > 0
generate double `dw_p1r' = 0
replace 		`dw_p1r' = `pw_p1r' if `pw_p1r' > 0
replace			`dw_p1r' = `pw_p1r' - `pw_p2r' if `dw_p2r' > 0 & `pw_p1r' > 0
replace			`dw_p1r' = `pw_p1r' - `pw_p3r' if `dw_p3r' > 0 & `pw_p2r' == 0 & `pw_p1r' > 0

*Check the decision weights
generate double `dw_check_r' = 0
replace `dw_check_r' = `dw_p3r' + `dw_p2r' + `dw_p1r'

*Check the probabilities
/* noisily summ `dw_check_l' `dw_check_r' `p1l' `p2l' `p3l' `p1r' `p2r' `p3r' ///
`pw_p1l' `pw_p2l' `pw_p3l' `pw_p1r' `pw_p2r' `pw_p3r' `dw_p1l' `dw_p2l' `dw_p3l' ///
`dw_p1r' `dw_p2r' `dw_p3r' */


*Construct the argument of the utility function 
if "$ufunc" == "power" {
*Left lottery
	generate double `u1L' = `prize1L'^`r' 
	generate double `u2L' = `prize2L'^`r' 
	generate double `u3L' = `prize3L'^`r'

	replace `u1L' = ln(`prize1L') if `r' == 0
	replace `u2L' = ln(`prize2L') if `r' == 0
	replace `u3L' = ln(`prize3L') if `r' == 0

	replace `u1L' = -`prize1L'^`r' if `r' < 0
	replace `u2L' = -`prize2L'^`r' if `r' < 0
	replace `u3L' = -`prize3L'^`r' if `r' < 0 

*Right lottery
	generate double `u1R' = `prize1R'^`r' 
	generate double `u2R' = `prize2R'^`r' 
	generate double `u3R' = `prize3R'^`r'

	replace `u1R' = ln(`prize1R') if `r' == 0
	replace `u2R' = ln(`prize2R') if `r' == 0
	replace `u3R' = ln(`prize3R') if `r' == 0

	replace `u1R' = -`prize1R'^`r' if `r' < 0
	replace `u2R' = -`prize2R'^`r' if `r' < 0
	replace `u3R' = -`prize3R'^`r' if `r' < 0 
	
*Capture discounting information
	generate double `uSS' = `ssamount'^`r'
	generate double `uLL' = `llamount'^`r'

	replace `uSS' = ln(`ssamount') if `r' == 0
	replace `uLL' = ln(`llamount') if `r' == 0

	replace `uSS' = -`ssamount'^`r' if `r' < 0
	replace `uLL' = -`llamount'^`r' if `r' < 0
	
}
else if "$ufunc" == "crra" {
*Left lottery
	generate double `u1L' = (`prize1L'^(1-`r'))/(1-`r') 
	generate double `u2L' = (`prize2L'^(1-`r'))/(1-`r')
	generate double `u3L' = (`prize3L'^(1-`r'))/(1-`r')

	replace `u1L' = ln(`prize1L') if `r' == 1
	replace `u2L' = ln(`prize2L') if `r' == 1
	replace `u3L' = ln(`prize3L') if `r' == 1


*Right lottery
	generate double `u1R' = (`prize1R'^(1-`r'))/(1-`r') 
	generate double `u2R' = (`prize2R'^(1-`r'))/(1-`r')
	generate double `u3R' = (`prize3R'^(1-`r'))/(1-`r')

	replace `u1R' = ln(`prize1R') if `r' == 1
	replace `u2R' = ln(`prize2R') if `r' == 1
	replace `u3R' = ln(`prize3R') if `r' == 1
	
*Capture discounting information
	generate double `uSS' = (`ssamount'^(1-`r'))/(1-`r')
	generate double `uLL' = (`llamount'^(1-`r'))/(1-`r')

	replace `uSS' = ln(`ssamount') if `r' == 1
	replace `uLL' = ln(`llamount') if `r' == 1

}
else if "$ufunc" == "cara" {
*Left lottery
	generate double `u1L' = (1-exp(-`r'*`prize1L'))/`r'
	generate double `u2L' = (1-exp(-`r'*`prize2L'))/`r'
	generate double `u3L' = (1-exp(-`r'*`prize3L'))/`r'

*Right lottery
	generate double `u1R' = (1-exp(-`r'*`prize1R'))/`r'
	generate double `u2R' = (1-exp(-`r'*`prize2R'))/`r'
	generate double `u3R' = (1-exp(-`r'*`prize3R'))/`r'

*Capture discounting information
	generate double `uSS' = (1-exp(-`r'*`ssamount'))/`r'
	generate double `uLL' = (1-exp(-`r'*`llamount'))/`r'

}

*Calculate EU of each lottery 
generate double `euL' = (`dw_p1l'*`u1L')+(`dw_p2l'*`u2L')+(`dw_p3l'*`u3L')
generate double `euR' = (`dw_p1r'*`u1R')+(`dw_p2r'*`u2R')+(`dw_p3r'*`u3R')

*Calculate the EU difference by using Contextual errors 
if "$ufunc" == "power" {
	generate double `uHigh' = `uMax'^`r'
	generate double `uLow' = `uMin'^`r'
				
	replace `uHigh' = ln(`uMax') if `r' == 0
	replace `uLow' = ln(`uMin') if `r' == 0
				
	replace `uHigh' = -`uMax'^`r' if `r' < 0
	replace `uLow' = -`uMin'^`r' if `r' < 0
}
else if "$ufunc" == "crra" {
	generate double `uHigh' = (`uMax'^(1-`r'))/(1-`r')
	generate double `uLow' = (`uMin'^(1-`r'))/(1-`r')
				
	replace `uHigh' = ln(`uMax') if `r' == 1
	replace `uLow' = ln(`uMin') if `r' == 1
}
else if "$ufunc" == "cara" {				
	generate double `uHigh' = (1-exp(-`r'*`uMax'))/`r' 
	generate double `uLow' = (1-exp(-`r'*`uMin'))/`r'
}

generate double `euDiff' = ((`euR' - `euL')/(`uHigh' - `uLow'))/`noiseRA'

*Evaluate the likelihood for risky choices
replace `lnf' = ln($cdf( `euDiff')) if `choice'==1 & `risk' == 1 
replace `lnf' = ln($cdf(-`euDiff')) if `choice'==0 & `risk' == 1

 
*Now focus on time preference data (i.e. present value)
*Exponential
*SS
generate double `dSS' = 1
replace `dSS' = 1/((1 +`deltaE')^`ssdelay')
generate double `PVss' = `dSS' * `uSS'

*LL
generate double `dLL' = 1
replace `dLL' = 1/((1 +`deltaE')^`lldelay')
generate double `PVll' = `dLL' * `uLL'

*Calculate the PV difference using Fechner errors
generate double `pvDiff' = (`PVss' - `PVll')/`noiseDR' if `risk' == 0

*Evaluate the likelihood for time choices (Exponential)
generate double `lnfED' = ln($cdf( `pvDiff')) if `choice'==0 & `risk' == 0 
replace `lnfED' = ln($cdf(-`pvDiff')) if `choice'==1 & `risk' == 0



*Weibull
*SS
replace `dSS' = 1 
replace `dSS' = exp(-`deltaWB'*(`ssdelay')^(1/`betaWB'))
replace `PVss' = `dSS'*`uSS' 

*LL
replace `dLL' = 1 
replace `dLL' = exp(-`deltaWB'*(`lldelay')^(1/`betaWB'))
replace `PVll' = `dLL'*`uLL'


*Calculate the PV difference using Fechner errors
replace `pvDiff' = (`PVss' - `PVll')/`noiseDR' if `risk' == 0

*Evaluate the likelihood for time choices (Weibull)
generate double `lnfWB' = ln($cdf( `pvDiff')) if `choice'==0 & `risk' == 0 
replace `lnfWB' = ln($cdf(-`pvDiff')) if `choice'==1 & `risk' == 0

*Calculate the grand likelihood for the discounting choices
generate double `f1' = exp(`lnfED') if `risk' == 0
generate double `f2' = exp(`lnfWB') if `risk' == 0
generate double `p1' = 1/(1+exp(`kappa')) if `risk' == 0
generate double `p2' = exp(`kappa')/(1+exp(`kappa')) if `risk' == 0
replace	`lnf' = ln((`p1' * `f1') + (`p2' * `f2')) if `risk' == 0

}
end


*Mixture model of Hyperbolic and Quasi-Hyperbolic Discounting (RDU utility function)
program define ml_rdu_discount_mixed_hqh

*Remember to define the globals $cdf, $ufunc

* Specify the arguments of this program 
if "$weigh" == "prelec2" {
    args lnf r phi eta deltaH betaQH deltaQH noiseRA noiseDR kappa
}
else {
    args lnf r gamma deltaH betaQH deltaQH noiseRA noiseDR kappa
}

* Declare the temporary variables to be used 
tempvar choice p1l p2l p3l p1r p2r p3r prize1L prize2L prize3L prize1R prize2R prize3R ///
u1L u2L u3L u1R u2R u3R euL euR euDiff uHigh uLow uMax uMin
tempvar p1l p2l p3l p1r p2r p3r
tempvar pw_p1l pw_p2l pw_p3l pw_p1r pw_p2r pw_p3r
tempvar dw_p1l dw_p2l dw_p3l dw_p1r dw_p2r dw_p3r
tempvar dw_check_l dw_check_r 
tempvar risk ssamount ssdelay llamount lldelay uSS uLL dSS dLL PVss PVll pvDiff
*Mixture model variables
tempvar lnfHD lnfQHD f1 f2 p1 p2    

quietly {

*Read in the data

generate int `choice' = $ML_y1

generate double `prize1L' = $ML_y8 
generate double `prize2L' = $ML_y9 
generate double `prize3L' = $ML_y10 

generate double `prize1R' = $ML_y11 
generate double `prize2R' = $ML_y12 
generate double `prize3R' = $ML_y13 

generate double `uMax' = $ML_y14
generate double `uMin' = $ML_y15

generate int `risk' = $ML_y16
generate double `ssamount' = $ML_y17
generate double `ssdelay' = $ML_y18/365
generate double `llamount' = $ML_y19
generate double `lldelay' = $ML_y20/365

*Generate the aggregate probabilities for the Left lottery
generate double `p3l' = $ML_y4
generate double `p2l' = $ML_y3
replace `p2l' = `p2l' + `p3l' if `p2l' > 0
generate double `p1l' = $ML_y2
replace `p1l' = `p1l' + `p2l' if `p1l' > 0
replace `p1l' = `p1l' + `p3l' if `p1l' > 0 & `p2l' == 0

*Generate the weighted probabilities for the Left lottery using TK (1992), Power or Prelec weighting functions
if "$weigh" == "tk" {
	forvalues i = 1/3 {
		generate double `pw_p`i'l' = (`p`i'l'^`gamma')/(`p`i'l'^`gamma' + ((1-`p`i'l')^`gamma'))^(1/`gamma')
		*This replace ensures that cumulative probabilities of 1 are set as 1
		replace `pw_p`i'l' = 1 if `pw_p`i'l' == .
	}
}
else if "$weigh" == "power" {
	forvalues i = 1/3 {
		generate double `pw_p`i'l' = (`p`i'l'^`gamma')
	}
}
else if "$weigh" == "prelec" {
	forvalues i = 1/3 {
		generate double `pw_p`i'l' = `p`i'l'
		replace `pw_p`i'l' = exp(-((-ln(`p`i'l'))^`gamma')) if `p`i'l' > 0 & `p`i'l' < 1
	}
}
else if "$weigh" == "prelec2" {
	forvalues i = 1/3 {
		generate double `pw_p`i'l' = `p`i'l'
		replace `pw_p`i'l' = exp((-`eta')*((-ln(`p`i'l'))^`phi')) if `p`i'l' > 0 & `p`i'l' < 1
	}
}


*Generate the decision weights for the Left lottery
generate double `dw_p3l' = 0
replace 		`dw_p3l' = `pw_p3l' if `pw_p3l' > 0
generate double `dw_p2l' = 0
replace 		`dw_p2l' = `pw_p2l' if `pw_p2l' > 0
replace			`dw_p2l' = `pw_p2l' - `pw_p3l' if `dw_p3l' > 0 & `pw_p2l' > 0
generate double `dw_p1l' = 0
replace 		`dw_p1l' = `pw_p1l' if `pw_p1l' > 0
replace			`dw_p1l' = `pw_p1l' - `pw_p2l' if `dw_p2l' > 0 & `pw_p1l' > 0
replace			`dw_p1l' = `pw_p1l' - `pw_p3l' if `dw_p3l' > 0 & `pw_p2l' == 0 & `pw_p1l' > 0

*Check the decision weights
generate double `dw_check_l' = 0
replace `dw_check_l' = `dw_p3l' + `dw_p2l' + `dw_p1l'

*Now, generate the aggregate probabilities for the Right lottery
generate double `p3r' = $ML_y7
generate double `p2r' = $ML_y6
replace `p2r' = `p2r' + `p3r' if `p2r' > 0
generate double `p1r' = $ML_y5
replace `p1r' = `p1r' + `p2r' if `p1r' > 0
replace `p1r' = `p1r' + `p3r' if `p1r' > 0 & `p2r' == 0

*Generate the weighted probabilities for the Right lottery
if "$weigh" == "tk" {
	forvalues i = 1/3 {
		generate double `pw_p`i'r' = (`p`i'r'^`gamma')/(`p`i'r'^`gamma' + ((1-`p`i'r')^`gamma'))^(1/`gamma')
		*This replace ensures that cumulative probabilities of 1 are set as 1
		replace `pw_p`i'r' = 1 if `pw_p`i'r' == .
	}
}
else if "$weigh" == "power" {
	forvalues i = 1/3 {
		generate double `pw_p`i'r' = (`p`i'r'^`gamma')
	}
}
else if "$weigh" == "prelec" {
	forvalues i = 1/3 {
		generate double `pw_p`i'r' = `p`i'r'
		replace `pw_p`i'r' = exp(-((-ln(`p`i'r'))^`gamma')) if `p`i'r' > 0 & `p`i'r' < 1
	}
}
else if "$weigh" == "prelec2" {
	forvalues i = 1/3 {
		generate double `pw_p`i'r' = `p`i'r'
		replace `pw_p`i'r' = exp((-`eta')*((-ln(`p`i'r'))^`phi')) if `p`i'r' > 0 & `p`i'r' < 1
	}
}


*Generate the decision weights for the Right lottery
generate double `dw_p3r' = 0
replace 		`dw_p3r' = `pw_p3r' if `pw_p3r' > 0
generate double `dw_p2r' = 0
replace 		`dw_p2r' = `pw_p2r' if `pw_p2r' > 0
replace			`dw_p2r' = `pw_p2r' - `pw_p3r' if `dw_p3r' > 0 & `pw_p2r' > 0
generate double `dw_p1r' = 0
replace 		`dw_p1r' = `pw_p1r' if `pw_p1r' > 0
replace			`dw_p1r' = `pw_p1r' - `pw_p2r' if `dw_p2r' > 0 & `pw_p1r' > 0
replace			`dw_p1r' = `pw_p1r' - `pw_p3r' if `dw_p3r' > 0 & `pw_p2r' == 0 & `pw_p1r' > 0

*Check the decision weights
generate double `dw_check_r' = 0
replace `dw_check_r' = `dw_p3r' + `dw_p2r' + `dw_p1r'

*Check the probabilities
/* noisily summ `dw_check_l' `dw_check_r' `p1l' `p2l' `p3l' `p1r' `p2r' `p3r' ///
`pw_p1l' `pw_p2l' `pw_p3l' `pw_p1r' `pw_p2r' `pw_p3r' `dw_p1l' `dw_p2l' `dw_p3l' ///
`dw_p1r' `dw_p2r' `dw_p3r' */


*Construct the argument of the utility function 
if "$ufunc" == "power" {
*Left lottery
	generate double `u1L' = `prize1L'^`r' 
	generate double `u2L' = `prize2L'^`r' 
	generate double `u3L' = `prize3L'^`r'

	replace `u1L' = ln(`prize1L') if `r' == 0
	replace `u2L' = ln(`prize2L') if `r' == 0
	replace `u3L' = ln(`prize3L') if `r' == 0

	replace `u1L' = -`prize1L'^`r' if `r' < 0
	replace `u2L' = -`prize2L'^`r' if `r' < 0
	replace `u3L' = -`prize3L'^`r' if `r' < 0 

*Right lottery
	generate double `u1R' = `prize1R'^`r' 
	generate double `u2R' = `prize2R'^`r' 
	generate double `u3R' = `prize3R'^`r'

	replace `u1R' = ln(`prize1R') if `r' == 0
	replace `u2R' = ln(`prize2R') if `r' == 0
	replace `u3R' = ln(`prize3R') if `r' == 0

	replace `u1R' = -`prize1R'^`r' if `r' < 0
	replace `u2R' = -`prize2R'^`r' if `r' < 0
	replace `u3R' = -`prize3R'^`r' if `r' < 0 
	
*Capture discounting information
	generate double `uSS' = `ssamount'^`r'
	generate double `uLL' = `llamount'^`r'

	replace `uSS' = ln(`ssamount') if `r' == 0
	replace `uLL' = ln(`llamount') if `r' == 0

	replace `uSS' = -`ssamount'^`r' if `r' < 0
	replace `uLL' = -`llamount'^`r' if `r' < 0
	
}
else if "$ufunc" == "crra" {
*Left lottery
	generate double `u1L' = (`prize1L'^(1-`r'))/(1-`r') 
	generate double `u2L' = (`prize2L'^(1-`r'))/(1-`r')
	generate double `u3L' = (`prize3L'^(1-`r'))/(1-`r')

	replace `u1L' = ln(`prize1L') if `r' == 1
	replace `u2L' = ln(`prize2L') if `r' == 1
	replace `u3L' = ln(`prize3L') if `r' == 1


*Right lottery
	generate double `u1R' = (`prize1R'^(1-`r'))/(1-`r') 
	generate double `u2R' = (`prize2R'^(1-`r'))/(1-`r')
	generate double `u3R' = (`prize3R'^(1-`r'))/(1-`r')

	replace `u1R' = ln(`prize1R') if `r' == 1
	replace `u2R' = ln(`prize2R') if `r' == 1
	replace `u3R' = ln(`prize3R') if `r' == 1
	
*Capture discounting information
	generate double `uSS' = (`ssamount'^(1-`r'))/(1-`r')
	generate double `uLL' = (`llamount'^(1-`r'))/(1-`r')

	replace `uSS' = ln(`ssamount') if `r' == 1
	replace `uLL' = ln(`llamount') if `r' == 1

}
else if "$ufunc" == "cara" {
*Left lottery
	generate double `u1L' = (1-exp(-`r'*`prize1L'))/`r'
	generate double `u2L' = (1-exp(-`r'*`prize2L'))/`r'
	generate double `u3L' = (1-exp(-`r'*`prize3L'))/`r'

*Right lottery
	generate double `u1R' = (1-exp(-`r'*`prize1R'))/`r'
	generate double `u2R' = (1-exp(-`r'*`prize2R'))/`r'
	generate double `u3R' = (1-exp(-`r'*`prize3R'))/`r'

*Capture discounting information
	generate double `uSS' = (1-exp(-`r'*`ssamount'))/`r'
	generate double `uLL' = (1-exp(-`r'*`llamount'))/`r'

}

*Calculate EU of each lottery 
generate double `euL' = (`dw_p1l'*`u1L')+(`dw_p2l'*`u2L')+(`dw_p3l'*`u3L')
generate double `euR' = (`dw_p1r'*`u1R')+(`dw_p2r'*`u2R')+(`dw_p3r'*`u3R')

*Calculate the EU difference by using Contextual errors 
if "$ufunc" == "power" {
	generate double `uHigh' = `uMax'^`r'
	generate double `uLow' = `uMin'^`r'
				
	replace `uHigh' = ln(`uMax') if `r' == 0
	replace `uLow' = ln(`uMin') if `r' == 0
				
	replace `uHigh' = -`uMax'^`r' if `r' < 0
	replace `uLow' = -`uMin'^`r' if `r' < 0
}
else if "$ufunc" == "crra" {
	generate double `uHigh' = (`uMax'^(1-`r'))/(1-`r')
	generate double `uLow' = (`uMin'^(1-`r'))/(1-`r')
				
	replace `uHigh' = ln(`uMax') if `r' == 1
	replace `uLow' = ln(`uMin') if `r' == 1
}
else if "$ufunc" == "cara" {				
	generate double `uHigh' = (1-exp(-`r'*`uMax'))/`r' 
	generate double `uLow' = (1-exp(-`r'*`uMin'))/`r'
}

generate double `euDiff' = ((`euR' - `euL')/(`uHigh' - `uLow'))/`noiseRA'

*Evaluate the likelihood for risky choices
replace `lnf' = ln($cdf( `euDiff')) if `choice'==1 & `risk' == 1 
replace `lnf' = ln($cdf(-`euDiff')) if `choice'==0 & `risk' == 1

 
*Now focus on time preference data (i.e. present value)
*Hyperbolic
*SS
generate double `dSS' = 1
replace `dSS' = 1/(1 + (`deltaH'*`ssdelay'))
generate double `PVss' = `dSS' * `uSS'

*LL
generate double `dLL' = 1
replace `dLL' = 1/(1 + (`deltaH'*`lldelay'))
generate double `PVll' = `dLL' * `uLL'


*Calculate the PV difference using Fechner errors
generate double `pvDiff' = (`PVss' - `PVll')/`noiseDR' if `risk' == 0

*Evaluate the likelihood for time choices (Hyperbolic)
generate double `lnfHD' = ln($cdf( `pvDiff')) if `choice'==0 & `risk' == 0 
replace `lnfHD' = ln($cdf(-`pvDiff')) if `choice'==1 & `risk' == 0



*Quasi-Hyperbolic
*SS
replace `dSS' = 1
replace `dSS' = `betaQH'/((1 + `deltaQH')^`ssdelay') if `ssdelay' > 0
replace `PVss' = `dSS' * `uSS'

*LL
replace `dLL' = 1
replace `dLL' = `betaQH'/((1 + `deltaQH')^`lldelay') if `lldelay' > 0
replace `PVll' = `dLL' * `uLL'


*Calculate the PV difference using Fechner errors
replace `pvDiff' = (`PVss' - `PVll')/`noiseDR' if `risk' == 0

*Evaluate the likelihood for time choices (Quasi-Hyperbolic)
generate double `lnfQHD' = ln($cdf( `pvDiff')) if `choice'==0 & `risk' == 0 
replace `lnfQHD' = ln($cdf(-`pvDiff')) if `choice'==1 & `risk' == 0

*Calculate the grand likelihood for the discounting choices
generate double `f1' = exp(`lnfHD') if `risk' == 0
generate double `f2' = exp(`lnfQHD') if `risk' == 0
generate double `p1' = 1/(1+exp(`kappa')) if `risk' == 0
generate double `p2' = exp(`kappa')/(1+exp(`kappa')) if `risk' == 0
replace	`lnf' = ln((`p1' * `f1') + (`p2' * `f2')) if `risk' == 0

}
end


*Mixture model of Hyperbolic and Weibull Discounting (RDU utility function)
program define ml_rdu_discount_mixed_hwb

*Remember to define the globals $cdf, $ufunc

* Specify the arguments of this program 
if "$weigh" == "prelec2" {
    args lnf r phi eta deltaH deltaWB betaWB noiseRA noiseDR kappa
}
else {
    args lnf r gamma deltaH deltaWB betaWB noiseRA noiseDR kappa
}

* Declare the temporary variables to be used 
tempvar choice p1l p2l p3l p1r p2r p3r prize1L prize2L prize3L prize1R prize2R prize3R ///
u1L u2L u3L u1R u2R u3R euL euR euDiff uHigh uLow uMax uMin
tempvar p1l p2l p3l p1r p2r p3r
tempvar pw_p1l pw_p2l pw_p3l pw_p1r pw_p2r pw_p3r
tempvar dw_p1l dw_p2l dw_p3l dw_p1r dw_p2r dw_p3r
tempvar dw_check_l dw_check_r 
tempvar risk ssamount ssdelay llamount lldelay uSS uLL dSS dLL PVss PVll pvDiff
*Mixture model variables
tempvar lnfHD lnfWB f1 f2 p1 p2    

quietly {

*Read in the data

generate int `choice' = $ML_y1

generate double `prize1L' = $ML_y8 
generate double `prize2L' = $ML_y9 
generate double `prize3L' = $ML_y10 

generate double `prize1R' = $ML_y11 
generate double `prize2R' = $ML_y12 
generate double `prize3R' = $ML_y13 

generate double `uMax' = $ML_y14
generate double `uMin' = $ML_y15

generate int `risk' = $ML_y16
generate double `ssamount' = $ML_y17
generate double `ssdelay' = $ML_y18/365
generate double `llamount' = $ML_y19
generate double `lldelay' = $ML_y20/365

*Generate the aggregate probabilities for the Left lottery
generate double `p3l' = $ML_y4
generate double `p2l' = $ML_y3
replace `p2l' = `p2l' + `p3l' if `p2l' > 0
generate double `p1l' = $ML_y2
replace `p1l' = `p1l' + `p2l' if `p1l' > 0
replace `p1l' = `p1l' + `p3l' if `p1l' > 0 & `p2l' == 0

*Generate the weighted probabilities for the Left lottery using TK (1992), Power or Prelec weighting functions
if "$weigh" == "tk" {
	forvalues i = 1/3 {
		generate double `pw_p`i'l' = (`p`i'l'^`gamma')/(`p`i'l'^`gamma' + ((1-`p`i'l')^`gamma'))^(1/`gamma')
		*This replace ensures that cumulative probabilities of 1 are set as 1
		replace `pw_p`i'l' = 1 if `pw_p`i'l' == .
	}
}
else if "$weigh" == "power" {
	forvalues i = 1/3 {
		generate double `pw_p`i'l' = (`p`i'l'^`gamma')
	}
}
else if "$weigh" == "prelec" {
	forvalues i = 1/3 {
		generate double `pw_p`i'l' = `p`i'l'
		replace `pw_p`i'l' = exp(-((-ln(`p`i'l'))^`gamma')) if `p`i'l' > 0 & `p`i'l' < 1
	}
}
else if "$weigh" == "prelec2" {
	forvalues i = 1/3 {
		generate double `pw_p`i'l' = `p`i'l'
		replace `pw_p`i'l' = exp((-`eta')*((-ln(`p`i'l'))^`phi')) if `p`i'l' > 0 & `p`i'l' < 1
	}
}


*Generate the decision weights for the Left lottery
generate double `dw_p3l' = 0
replace 		`dw_p3l' = `pw_p3l' if `pw_p3l' > 0
generate double `dw_p2l' = 0
replace 		`dw_p2l' = `pw_p2l' if `pw_p2l' > 0
replace			`dw_p2l' = `pw_p2l' - `pw_p3l' if `dw_p3l' > 0 & `pw_p2l' > 0
generate double `dw_p1l' = 0
replace 		`dw_p1l' = `pw_p1l' if `pw_p1l' > 0
replace			`dw_p1l' = `pw_p1l' - `pw_p2l' if `dw_p2l' > 0 & `pw_p1l' > 0
replace			`dw_p1l' = `pw_p1l' - `pw_p3l' if `dw_p3l' > 0 & `pw_p2l' == 0 & `pw_p1l' > 0

*Check the decision weights
generate double `dw_check_l' = 0
replace `dw_check_l' = `dw_p3l' + `dw_p2l' + `dw_p1l'

*Now, generate the aggregate probabilities for the Right lottery
generate double `p3r' = $ML_y7
generate double `p2r' = $ML_y6
replace `p2r' = `p2r' + `p3r' if `p2r' > 0
generate double `p1r' = $ML_y5
replace `p1r' = `p1r' + `p2r' if `p1r' > 0
replace `p1r' = `p1r' + `p3r' if `p1r' > 0 & `p2r' == 0

*Generate the weighted probabilities for the Right lottery
if "$weigh" == "tk" {
	forvalues i = 1/3 {
		generate double `pw_p`i'r' = (`p`i'r'^`gamma')/(`p`i'r'^`gamma' + ((1-`p`i'r')^`gamma'))^(1/`gamma')
		*This replace ensures that cumulative probabilities of 1 are set as 1
		replace `pw_p`i'r' = 1 if `pw_p`i'r' == .
	}
}
else if "$weigh" == "power" {
	forvalues i = 1/3 {
		generate double `pw_p`i'r' = (`p`i'r'^`gamma')
	}
}
else if "$weigh" == "prelec" {
	forvalues i = 1/3 {
		generate double `pw_p`i'r' = `p`i'r'
		replace `pw_p`i'r' = exp(-((-ln(`p`i'r'))^`gamma')) if `p`i'r' > 0 & `p`i'r' < 1
	}
}
else if "$weigh" == "prelec2" {
	forvalues i = 1/3 {
		generate double `pw_p`i'r' = `p`i'r'
		replace `pw_p`i'r' = exp((-`eta')*((-ln(`p`i'r'))^`phi')) if `p`i'r' > 0 & `p`i'r' < 1
	}
}


*Generate the decision weights for the Right lottery
generate double `dw_p3r' = 0
replace 		`dw_p3r' = `pw_p3r' if `pw_p3r' > 0
generate double `dw_p2r' = 0
replace 		`dw_p2r' = `pw_p2r' if `pw_p2r' > 0
replace			`dw_p2r' = `pw_p2r' - `pw_p3r' if `dw_p3r' > 0 & `pw_p2r' > 0
generate double `dw_p1r' = 0
replace 		`dw_p1r' = `pw_p1r' if `pw_p1r' > 0
replace			`dw_p1r' = `pw_p1r' - `pw_p2r' if `dw_p2r' > 0 & `pw_p1r' > 0
replace			`dw_p1r' = `pw_p1r' - `pw_p3r' if `dw_p3r' > 0 & `pw_p2r' == 0 & `pw_p1r' > 0

*Check the decision weights
generate double `dw_check_r' = 0
replace `dw_check_r' = `dw_p3r' + `dw_p2r' + `dw_p1r'

*Check the probabilities
/* noisily summ `dw_check_l' `dw_check_r' `p1l' `p2l' `p3l' `p1r' `p2r' `p3r' ///
`pw_p1l' `pw_p2l' `pw_p3l' `pw_p1r' `pw_p2r' `pw_p3r' `dw_p1l' `dw_p2l' `dw_p3l' ///
`dw_p1r' `dw_p2r' `dw_p3r' */


*Construct the argument of the utility function 
if "$ufunc" == "power" {
*Left lottery
	generate double `u1L' = `prize1L'^`r' 
	generate double `u2L' = `prize2L'^`r' 
	generate double `u3L' = `prize3L'^`r'

	replace `u1L' = ln(`prize1L') if `r' == 0
	replace `u2L' = ln(`prize2L') if `r' == 0
	replace `u3L' = ln(`prize3L') if `r' == 0

	replace `u1L' = -`prize1L'^`r' if `r' < 0
	replace `u2L' = -`prize2L'^`r' if `r' < 0
	replace `u3L' = -`prize3L'^`r' if `r' < 0 

*Right lottery
	generate double `u1R' = `prize1R'^`r' 
	generate double `u2R' = `prize2R'^`r' 
	generate double `u3R' = `prize3R'^`r'

	replace `u1R' = ln(`prize1R') if `r' == 0
	replace `u2R' = ln(`prize2R') if `r' == 0
	replace `u3R' = ln(`prize3R') if `r' == 0

	replace `u1R' = -`prize1R'^`r' if `r' < 0
	replace `u2R' = -`prize2R'^`r' if `r' < 0
	replace `u3R' = -`prize3R'^`r' if `r' < 0 
	
*Capture discounting information
	generate double `uSS' = `ssamount'^`r'
	generate double `uLL' = `llamount'^`r'

	replace `uSS' = ln(`ssamount') if `r' == 0
	replace `uLL' = ln(`llamount') if `r' == 0

	replace `uSS' = -`ssamount'^`r' if `r' < 0
	replace `uLL' = -`llamount'^`r' if `r' < 0
	
}
else if "$ufunc" == "crra" {
*Left lottery
	generate double `u1L' = (`prize1L'^(1-`r'))/(1-`r') 
	generate double `u2L' = (`prize2L'^(1-`r'))/(1-`r')
	generate double `u3L' = (`prize3L'^(1-`r'))/(1-`r')

	replace `u1L' = ln(`prize1L') if `r' == 1
	replace `u2L' = ln(`prize2L') if `r' == 1
	replace `u3L' = ln(`prize3L') if `r' == 1


*Right lottery
	generate double `u1R' = (`prize1R'^(1-`r'))/(1-`r') 
	generate double `u2R' = (`prize2R'^(1-`r'))/(1-`r')
	generate double `u3R' = (`prize3R'^(1-`r'))/(1-`r')

	replace `u1R' = ln(`prize1R') if `r' == 1
	replace `u2R' = ln(`prize2R') if `r' == 1
	replace `u3R' = ln(`prize3R') if `r' == 1
	
*Capture discounting information
	generate double `uSS' = (`ssamount'^(1-`r'))/(1-`r')
	generate double `uLL' = (`llamount'^(1-`r'))/(1-`r')

	replace `uSS' = ln(`ssamount') if `r' == 1
	replace `uLL' = ln(`llamount') if `r' == 1

}
else if "$ufunc" == "cara" {
*Left lottery
	generate double `u1L' = (1-exp(-`r'*`prize1L'))/`r'
	generate double `u2L' = (1-exp(-`r'*`prize2L'))/`r'
	generate double `u3L' = (1-exp(-`r'*`prize3L'))/`r'

*Right lottery
	generate double `u1R' = (1-exp(-`r'*`prize1R'))/`r'
	generate double `u2R' = (1-exp(-`r'*`prize2R'))/`r'
	generate double `u3R' = (1-exp(-`r'*`prize3R'))/`r'

*Capture discounting information
	generate double `uSS' = (1-exp(-`r'*`ssamount'))/`r'
	generate double `uLL' = (1-exp(-`r'*`llamount'))/`r'

}

*Calculate EU of each lottery 
generate double `euL' = (`dw_p1l'*`u1L')+(`dw_p2l'*`u2L')+(`dw_p3l'*`u3L')
generate double `euR' = (`dw_p1r'*`u1R')+(`dw_p2r'*`u2R')+(`dw_p3r'*`u3R')

*Calculate the EU difference by using Contextual errors 
if "$ufunc" == "power" {
	generate double `uHigh' = `uMax'^`r'
	generate double `uLow' = `uMin'^`r'
				
	replace `uHigh' = ln(`uMax') if `r' == 0
	replace `uLow' = ln(`uMin') if `r' == 0
				
	replace `uHigh' = -`uMax'^`r' if `r' < 0
	replace `uLow' = -`uMin'^`r' if `r' < 0
}
else if "$ufunc" == "crra" {
	generate double `uHigh' = (`uMax'^(1-`r'))/(1-`r')
	generate double `uLow' = (`uMin'^(1-`r'))/(1-`r')
				
	replace `uHigh' = ln(`uMax') if `r' == 1
	replace `uLow' = ln(`uMin') if `r' == 1
}
else if "$ufunc" == "cara" {				
	generate double `uHigh' = (1-exp(-`r'*`uMax'))/`r' 
	generate double `uLow' = (1-exp(-`r'*`uMin'))/`r'
}

generate double `euDiff' = ((`euR' - `euL')/(`uHigh' - `uLow'))/`noiseRA'

*Evaluate the likelihood for risky choices
replace `lnf' = ln($cdf( `euDiff')) if `choice'==1 & `risk' == 1 
replace `lnf' = ln($cdf(-`euDiff')) if `choice'==0 & `risk' == 1

 
*Now focus on time preference data (i.e. present value)
*Hyperbolic
*SS
generate double `dSS' = 1
replace `dSS' = 1/(1 + (`deltaH'*`ssdelay'))
generate double `PVss' = `dSS' * `uSS'

*LL
generate double `dLL' = 1
replace `dLL' = 1/(1 + (`deltaH'*`lldelay'))
generate double `PVll' = `dLL' * `uLL'


*Calculate the PV difference using Fechner errors
generate double `pvDiff' = (`PVss' - `PVll')/`noiseDR' if `risk' == 0

*Evaluate the likelihood for time choices (Hyperbolic)
generate double `lnfHD' = ln($cdf( `pvDiff')) if `choice'==0 & `risk' == 0 
replace `lnfHD' = ln($cdf(-`pvDiff')) if `choice'==1 & `risk' == 0



*Weibull
*SS
replace `dSS' = 1 
replace `dSS' = exp(-`deltaWB'*(`ssdelay')^(1/`betaWB'))
replace `PVss' = `dSS'*`uSS' 

*LL
replace `dLL' = 1 
replace `dLL' = exp(-`deltaWB'*(`lldelay')^(1/`betaWB'))
replace `PVll' = `dLL'*`uLL'


*Calculate the PV difference using Fechner errors
replace `pvDiff' = (`PVss' - `PVll')/`noiseDR' if `risk' == 0

*Evaluate the likelihood for time choices (Quasi-Hyperbolic)
generate double `lnfWB' = ln($cdf( `pvDiff')) if `choice'==0 & `risk' == 0 
replace `lnfWB' = ln($cdf(-`pvDiff')) if `choice'==1 & `risk' == 0

*Calculate the grand likelihood for the discounting choices
generate double `f1' = exp(`lnfHD') if `risk' == 0
generate double `f2' = exp(`lnfWB') if `risk' == 0
generate double `p1' = 1/(1+exp(`kappa')) if `risk' == 0
generate double `p2' = exp(`kappa')/(1+exp(`kappa')) if `risk' == 0
replace	`lnf' = ln((`p1' * `f1') + (`p2' * `f2')) if `risk' == 0

}
end


*Mixture model of Quasi-Hyperbolic and Weibull Discounting (RDU utility function)
program define ml_rdu_discount_mixed_qhwb

*Remember to define the globals $cdf, $ufunc

* Specify the arguments of this program 
if "$weigh" == "prelec2" {
    args lnf r phi eta betaQH deltaQH deltaWB betaWB noiseRA noiseDR kappa
}
else {
    args lnf r gamma betaQH deltaQH deltaWB betaWB noiseRA noiseDR kappa
}

* Declare the temporary variables to be used 
tempvar choice p1l p2l p3l p1r p2r p3r prize1L prize2L prize3L prize1R prize2R prize3R ///
u1L u2L u3L u1R u2R u3R euL euR euDiff uHigh uLow uMax uMin
tempvar p1l p2l p3l p1r p2r p3r
tempvar pw_p1l pw_p2l pw_p3l pw_p1r pw_p2r pw_p3r
tempvar dw_p1l dw_p2l dw_p3l dw_p1r dw_p2r dw_p3r
tempvar dw_check_l dw_check_r 
tempvar risk ssamount ssdelay llamount lldelay uSS uLL dSS dLL PVss PVll pvDiff
*Mixture model variables
tempvar lnfQHD lnfWB f1 f2 p1 p2    

quietly {

*Read in the data

generate int `choice' = $ML_y1

generate double `prize1L' = $ML_y8 
generate double `prize2L' = $ML_y9 
generate double `prize3L' = $ML_y10 

generate double `prize1R' = $ML_y11 
generate double `prize2R' = $ML_y12 
generate double `prize3R' = $ML_y13 

generate double `uMax' = $ML_y14
generate double `uMin' = $ML_y15

generate int `risk' = $ML_y16
generate double `ssamount' = $ML_y17
generate double `ssdelay' = $ML_y18/365
generate double `llamount' = $ML_y19
generate double `lldelay' = $ML_y20/365

*Generate the aggregate probabilities for the Left lottery
generate double `p3l' = $ML_y4
generate double `p2l' = $ML_y3
replace `p2l' = `p2l' + `p3l' if `p2l' > 0
generate double `p1l' = $ML_y2
replace `p1l' = `p1l' + `p2l' if `p1l' > 0
replace `p1l' = `p1l' + `p3l' if `p1l' > 0 & `p2l' == 0

*Generate the weighted probabilities for the Left lottery using TK (1992), Power or Prelec weighting functions
if "$weigh" == "tk" {
	forvalues i = 1/3 {
		generate double `pw_p`i'l' = (`p`i'l'^`gamma')/(`p`i'l'^`gamma' + ((1-`p`i'l')^`gamma'))^(1/`gamma')
		*This replace ensures that cumulative probabilities of 1 are set as 1
		replace `pw_p`i'l' = 1 if `pw_p`i'l' == .
	}
}
else if "$weigh" == "power" {
	forvalues i = 1/3 {
		generate double `pw_p`i'l' = (`p`i'l'^`gamma')
	}
}
else if "$weigh" == "prelec" {
	forvalues i = 1/3 {
		generate double `pw_p`i'l' = `p`i'l'
		replace `pw_p`i'l' = exp(-((-ln(`p`i'l'))^`gamma')) if `p`i'l' > 0 & `p`i'l' < 1
	}
}
else if "$weigh" == "prelec2" {
	forvalues i = 1/3 {
		generate double `pw_p`i'l' = `p`i'l'
		replace `pw_p`i'l' = exp((-`eta')*((-ln(`p`i'l'))^`phi')) if `p`i'l' > 0 & `p`i'l' < 1
	}
}


*Generate the decision weights for the Left lottery
generate double `dw_p3l' = 0
replace 		`dw_p3l' = `pw_p3l' if `pw_p3l' > 0
generate double `dw_p2l' = 0
replace 		`dw_p2l' = `pw_p2l' if `pw_p2l' > 0
replace			`dw_p2l' = `pw_p2l' - `pw_p3l' if `dw_p3l' > 0 & `pw_p2l' > 0
generate double `dw_p1l' = 0
replace 		`dw_p1l' = `pw_p1l' if `pw_p1l' > 0
replace			`dw_p1l' = `pw_p1l' - `pw_p2l' if `dw_p2l' > 0 & `pw_p1l' > 0
replace			`dw_p1l' = `pw_p1l' - `pw_p3l' if `dw_p3l' > 0 & `pw_p2l' == 0 & `pw_p1l' > 0

*Check the decision weights
generate double `dw_check_l' = 0
replace `dw_check_l' = `dw_p3l' + `dw_p2l' + `dw_p1l'

*Now, generate the aggregate probabilities for the Right lottery
generate double `p3r' = $ML_y7
generate double `p2r' = $ML_y6
replace `p2r' = `p2r' + `p3r' if `p2r' > 0
generate double `p1r' = $ML_y5
replace `p1r' = `p1r' + `p2r' if `p1r' > 0
replace `p1r' = `p1r' + `p3r' if `p1r' > 0 & `p2r' == 0

*Generate the weighted probabilities for the Right lottery
if "$weigh" == "tk" {
	forvalues i = 1/3 {
		generate double `pw_p`i'r' = (`p`i'r'^`gamma')/(`p`i'r'^`gamma' + ((1-`p`i'r')^`gamma'))^(1/`gamma')
		*This replace ensures that cumulative probabilities of 1 are set as 1
		replace `pw_p`i'r' = 1 if `pw_p`i'r' == .
	}
}
else if "$weigh" == "power" {
	forvalues i = 1/3 {
		generate double `pw_p`i'r' = (`p`i'r'^`gamma')
	}
}
else if "$weigh" == "prelec" {
	forvalues i = 1/3 {
		generate double `pw_p`i'r' = `p`i'r'
		replace `pw_p`i'r' = exp(-((-ln(`p`i'r'))^`gamma')) if `p`i'r' > 0 & `p`i'r' < 1
	}
}
else if "$weigh" == "prelec2" {
	forvalues i = 1/3 {
		generate double `pw_p`i'r' = `p`i'r'
		replace `pw_p`i'r' = exp((-`eta')*((-ln(`p`i'r'))^`phi')) if `p`i'r' > 0 & `p`i'r' < 1
	}
}


*Generate the decision weights for the Right lottery
generate double `dw_p3r' = 0
replace 		`dw_p3r' = `pw_p3r' if `pw_p3r' > 0
generate double `dw_p2r' = 0
replace 		`dw_p2r' = `pw_p2r' if `pw_p2r' > 0
replace			`dw_p2r' = `pw_p2r' - `pw_p3r' if `dw_p3r' > 0 & `pw_p2r' > 0
generate double `dw_p1r' = 0
replace 		`dw_p1r' = `pw_p1r' if `pw_p1r' > 0
replace			`dw_p1r' = `pw_p1r' - `pw_p2r' if `dw_p2r' > 0 & `pw_p1r' > 0
replace			`dw_p1r' = `pw_p1r' - `pw_p3r' if `dw_p3r' > 0 & `pw_p2r' == 0 & `pw_p1r' > 0

*Check the decision weights
generate double `dw_check_r' = 0
replace `dw_check_r' = `dw_p3r' + `dw_p2r' + `dw_p1r'

*Check the probabilities
/* noisily summ `dw_check_l' `dw_check_r' `p1l' `p2l' `p3l' `p1r' `p2r' `p3r' ///
`pw_p1l' `pw_p2l' `pw_p3l' `pw_p1r' `pw_p2r' `pw_p3r' `dw_p1l' `dw_p2l' `dw_p3l' ///
`dw_p1r' `dw_p2r' `dw_p3r' */


*Construct the argument of the utility function 
if "$ufunc" == "power" {
*Left lottery
	generate double `u1L' = `prize1L'^`r' 
	generate double `u2L' = `prize2L'^`r' 
	generate double `u3L' = `prize3L'^`r'

	replace `u1L' = ln(`prize1L') if `r' == 0
	replace `u2L' = ln(`prize2L') if `r' == 0
	replace `u3L' = ln(`prize3L') if `r' == 0

	replace `u1L' = -`prize1L'^`r' if `r' < 0
	replace `u2L' = -`prize2L'^`r' if `r' < 0
	replace `u3L' = -`prize3L'^`r' if `r' < 0 

*Right lottery
	generate double `u1R' = `prize1R'^`r' 
	generate double `u2R' = `prize2R'^`r' 
	generate double `u3R' = `prize3R'^`r'

	replace `u1R' = ln(`prize1R') if `r' == 0
	replace `u2R' = ln(`prize2R') if `r' == 0
	replace `u3R' = ln(`prize3R') if `r' == 0

	replace `u1R' = -`prize1R'^`r' if `r' < 0
	replace `u2R' = -`prize2R'^`r' if `r' < 0
	replace `u3R' = -`prize3R'^`r' if `r' < 0 
	
*Capture discounting information
	generate double `uSS' = `ssamount'^`r'
	generate double `uLL' = `llamount'^`r'

	replace `uSS' = ln(`ssamount') if `r' == 0
	replace `uLL' = ln(`llamount') if `r' == 0

	replace `uSS' = -`ssamount'^`r' if `r' < 0
	replace `uLL' = -`llamount'^`r' if `r' < 0
	
}
else if "$ufunc" == "crra" {
*Left lottery
	generate double `u1L' = (`prize1L'^(1-`r'))/(1-`r') 
	generate double `u2L' = (`prize2L'^(1-`r'))/(1-`r')
	generate double `u3L' = (`prize3L'^(1-`r'))/(1-`r')

	replace `u1L' = ln(`prize1L') if `r' == 1
	replace `u2L' = ln(`prize2L') if `r' == 1
	replace `u3L' = ln(`prize3L') if `r' == 1


*Right lottery
	generate double `u1R' = (`prize1R'^(1-`r'))/(1-`r') 
	generate double `u2R' = (`prize2R'^(1-`r'))/(1-`r')
	generate double `u3R' = (`prize3R'^(1-`r'))/(1-`r')

	replace `u1R' = ln(`prize1R') if `r' == 1
	replace `u2R' = ln(`prize2R') if `r' == 1
	replace `u3R' = ln(`prize3R') if `r' == 1
	
*Capture discounting information
	generate double `uSS' = (`ssamount'^(1-`r'))/(1-`r')
	generate double `uLL' = (`llamount'^(1-`r'))/(1-`r')

	replace `uSS' = ln(`ssamount') if `r' == 1
	replace `uLL' = ln(`llamount') if `r' == 1

}
else if "$ufunc" == "cara" {
*Left lottery
	generate double `u1L' = (1-exp(-`r'*`prize1L'))/`r'
	generate double `u2L' = (1-exp(-`r'*`prize2L'))/`r'
	generate double `u3L' = (1-exp(-`r'*`prize3L'))/`r'

*Right lottery
	generate double `u1R' = (1-exp(-`r'*`prize1R'))/`r'
	generate double `u2R' = (1-exp(-`r'*`prize2R'))/`r'
	generate double `u3R' = (1-exp(-`r'*`prize3R'))/`r'

*Capture discounting information
	generate double `uSS' = (1-exp(-`r'*`ssamount'))/`r'
	generate double `uLL' = (1-exp(-`r'*`llamount'))/`r'

}

*Calculate EU of each lottery 
generate double `euL' = (`dw_p1l'*`u1L')+(`dw_p2l'*`u2L')+(`dw_p3l'*`u3L')
generate double `euR' = (`dw_p1r'*`u1R')+(`dw_p2r'*`u2R')+(`dw_p3r'*`u3R')

*Calculate the EU difference by using Contextual errors 
if "$ufunc" == "power" {
	generate double `uHigh' = `uMax'^`r'
	generate double `uLow' = `uMin'^`r'
				
	replace `uHigh' = ln(`uMax') if `r' == 0
	replace `uLow' = ln(`uMin') if `r' == 0
				
	replace `uHigh' = -`uMax'^`r' if `r' < 0
	replace `uLow' = -`uMin'^`r' if `r' < 0
}
else if "$ufunc" == "crra" {
	generate double `uHigh' = (`uMax'^(1-`r'))/(1-`r')
	generate double `uLow' = (`uMin'^(1-`r'))/(1-`r')
				
	replace `uHigh' = ln(`uMax') if `r' == 1
	replace `uLow' = ln(`uMin') if `r' == 1
}
else if "$ufunc" == "cara" {				
	generate double `uHigh' = (1-exp(-`r'*`uMax'))/`r' 
	generate double `uLow' = (1-exp(-`r'*`uMin'))/`r'
}

generate double `euDiff' = ((`euR' - `euL')/(`uHigh' - `uLow'))/`noiseRA'

*Evaluate the likelihood for risky choices
replace `lnf' = ln($cdf( `euDiff')) if `choice'==1 & `risk' == 1 
replace `lnf' = ln($cdf(-`euDiff')) if `choice'==0 & `risk' == 1

 
*Now focus on time preference data (i.e. present value)
*Quasi-Hyperbolic
*SS
generate double `dSS' = 1
replace `dSS' = `betaQH'/((1 + `deltaQH')^`ssdelay') if `ssdelay' > 0
generate double `PVss' = `dSS' * `uSS'

*LL
generate double `dLL' = 1
replace `dLL' = `betaQH'/((1 + `deltaQH')^`lldelay') if `lldelay' > 0
generate double `PVll' = `dLL' * `uLL'


*Calculate the PV difference using Fechner errors
generate double `pvDiff' = (`PVss' - `PVll')/`noiseDR' if `risk' == 0

*Evaluate the likelihood for time choices (Hyperbolic)
generate double `lnfQHD' = ln($cdf( `pvDiff')) if `choice'==0 & `risk' == 0 
replace `lnfQHD' = ln($cdf(-`pvDiff')) if `choice'==1 & `risk' == 0



*Weibull
*SS
replace `dSS' = 1 
replace `dSS' = exp(-`deltaWB'*(`ssdelay')^(1/`betaWB'))
replace `PVss' = `dSS'*`uSS' 

*LL
replace `dLL' = 1 
replace `dLL' = exp(-`deltaWB'*(`lldelay')^(1/`betaWB'))
replace `PVll' = `dLL'*`uLL'


*Calculate the PV difference using Fechner errors
replace `pvDiff' = (`PVss' - `PVll')/`noiseDR' if `risk' == 0

*Evaluate the likelihood for time choices (Quasi-Hyperbolic)
generate double `lnfWB' = ln($cdf( `pvDiff')) if `choice'==0 & `risk' == 0 
replace `lnfWB' = ln($cdf(-`pvDiff')) if `choice'==1 & `risk' == 0

*Calculate the grand likelihood for the discounting choices
generate double `f1' = exp(`lnfQHD') if `risk' == 0
generate double `f2' = exp(`lnfWB') if `risk' == 0
generate double `p1' = 1/(1+exp(`kappa')) if `risk' == 0
generate double `p2' = exp(`kappa')/(1+exp(`kappa')) if `risk' == 0
replace	`lnf' = ln((`p1' * `f1') + (`p2' * `f2')) if `risk' == 0

}
end


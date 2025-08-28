*This do file loads the EUT and RDEU ML routines

capture program drop ml_eut_flex
capture program drop ml_rdu_flex
capture program drop ml_rdu_prelec2_LN

program define ml_eut_flex
*Remember to define the globals $error, $cdf and $ufunc prior to estimation

*Specify the arguments of this program 
if "$ufunc" == "expo" {
		if "$error" == "no" {
			args lnf r alpha
		}
		else {
			args lnf r alpha noise
		}
}
else if "$ufunc" == "expo2" {
		if "$error" == "no" {
			args lnf r alpha
		}
		else {
			args lnf r alpha noise
		}
}
else {
	if "$error" != "no" {
	args lnf r noise
	}
	else {
	args lnf r
	}

}

*Declare the temporary variables to be used 
tempvar choice p1l p2l p3l p1r p2r p3r prize1L prize2L prize3L prize1R prize2R prize3R ///
u1L u2L u3L u1R u2R u3R euL euR euDiff uHigh uLow uMax uMin


quietly {


* Read in the data 

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


}
else if "$ufunc" == "expo" {
*Left lottery
	generate double `u1L' = (1-exp(-`alpha'*((`prize1L')^(1-`r'))))/`alpha'
	generate double `u2L' = (1-exp(-`alpha'*((`prize2L')^(1-`r'))))/`alpha'
	generate double `u3L' = (1-exp(-`alpha'*((`prize3L')^(1-`r'))))/`alpha'

*Right lottery
	generate double `u1R' = (1-exp(-`alpha'*((`prize1R')^(1-`r'))))/`alpha'
	generate double `u2R' = (1-exp(-`alpha'*((`prize2R')^(1-`r'))))/`alpha'
	generate double `u3R' = (1-exp(-`alpha'*((`prize3R')^(1-`r'))))/`alpha'

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
}
else if "$ufunc" == "expo2" {
*Left lottery
	generate double `u1L' = -exp(-`alpha'*((`prize1L')^`r'))
	generate double `u2L' = -exp(-`alpha'*((`prize2L')^`r'))
	generate double `u3L' = -exp(-`alpha'*((`prize3L')^`r'))

*Right lottery
	generate double `u1R' = -exp(-`alpha'*((`prize1R')^`r'))
	generate double `u2R' = -exp(-`alpha'*((`prize2R')^`r'))
	generate double `u3R' = -exp(-`alpha'*((`prize3R')^`r'))
}

* Calculate EU of each lottery 
generate double `euL' = (`p1l'*`u1L')+(`p2l'*`u2L')+(`p3l'*`u3L')
generate double `euR' = (`p1r'*`u1R')+(`p2r'*`u2R')+(`p3r'*`u3R')


* Calculate the EU difference by using Fechner, Luce or Contextual errors 
if "$error" == "fechner" {
	generate double `euDiff' = (`euR'-`euL')/(`noise')
}
else if "$error" == "luce" {
	generate double `euDiff' = (`euR'^(1/`noise'))/((`euR'^(1/`noise'))+(`euL'^(1/`noise')))
}
else if "$error" == "context" {
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
			else if "$ufunc" == "expo" {
				generate double `uHigh' = (1-exp(-`alpha'*((`uMax')^(1-`r'))))/`alpha'
				generate double `uLow' = (1-exp(-`alpha'*((`uMin')^(1-`r'))))/`alpha'
			}
			else if "$ufunc" == "cara" {
				generate double `uHigh' = (1-exp(-`r'*`uMax'))/`r'
				generate double `uLow' = (1-exp(-`r'*`uMin'))/`r'
			}
			else if "$ufunc" == "expo2" {
				generate double `uHigh' = -exp(-`alpha'*((`uMax')^`r'))
				generate double `uLow' = -exp(-`alpha'*((`uMin')^`r'))
			}


	generate double `euDiff' = ((`euR' - `euL')/(`uHigh' - `uLow'))/`noise'
}
else if "$error" == "no" {
	generate double `euDiff' = `euR' - `euL'

}	



* Evaluate the likelihood 
if "$error" == "luce" {
	replace `lnf' = ln(`euDiff') if `choice'==1 
	replace `lnf' = ln(1-`euDiff') if `choice'==0
}
else { 
	replace `lnf' = ln($cdf( `euDiff')) if `choice'==1 
	replace `lnf' = ln($cdf(-`euDiff')) if `choice'==0
}

}
end



program define ml_rdu_flex
*Remember to define the globals $error, $cdf, $ufunc and $weigh prior to estimation

*Specify the arguments of this program 
	if "$error" == "no" {
		args lnf r gamma
	}
	else {
		if "$weigh" == "prelec2" {
			if "$ufunc" == "expo2" {
				args lnf r alpha phi eta noise
			}
			else {
				args lnf r phi eta noise
			}	
		}
		else {
			if "$ufunc" == "expo2" {
				args lnf r alpha gamma noise
			}
			else {
				args lnf r gamma noise
			}	
		}
	}
	
*Declare the temporary variables to be used 
tempvar choice prize1L prize2L prize3L prize1R prize2R prize3R 
tempvar u1L u2L u3L u1R u2R u3R euL euR euDiff uMax uMin uHigh uLow
tempvar p1l p2l p3l p1r p2r p3r
tempvar pw_p1l pw_p2l pw_p3l pw_p1r pw_p2r pw_p3r
tempvar dw_p1l dw_p2l dw_p3l dw_p1r dw_p2r dw_p3r
tempvar dw_check_l dw_check_r 


quietly {

*Initialize the data 
generate int `choice' = $ML_y1

generate double `prize1L' = $ML_y8 
generate double `prize2L' = $ML_y9 
generate double `prize3L' = $ML_y10 

generate double `prize1R' = $ML_y11 
generate double `prize2R' = $ML_y12 
generate double `prize3R' = $ML_y13 

generate double `uMax' = $ML_y14
generate double `uMin' = $ML_y15

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
}
else if "$ufunc" == "expo2" {
*Left lottery
	generate double `u1L' = -exp(-`alpha'*((`prize1L')^`r'))
	generate double `u2L' = -exp(-`alpha'*((`prize2L')^`r'))
	generate double `u3L' = -exp(-`alpha'*((`prize3L')^`r'))

*Right lottery
	generate double `u1R' = -exp(-`alpha'*((`prize1R')^`r'))
	generate double `u2R' = -exp(-`alpha'*((`prize2R')^`r'))
	generate double `u3R' = -exp(-`alpha'*((`prize3R')^`r'))
}



*Calculate RDU of each lottery 
generate double `euL' = (`dw_p1l'*`u1L')+(`dw_p2l'*`u2L')+(`dw_p3l'*`u3L')
generate double `euR' = (`dw_p1r'*`u1R')+(`dw_p2r'*`u2R')+(`dw_p3r'*`u3R')

* Calculate the RDU difference by using Fechner, Luce or Contextual errors 
if "$error" == "fechner" {
	generate double `euDiff' = (`euR'-`euL')/(`noise')
}
else if "$error" == "luce" {
	generate double `euDiff' = (`euR'^(1/`noise'))/((`euR'^(1/`noise'))+(`euL'^(1/`noise')))
}
else if "$error" == "context" {
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
			else if "$ufunc" == "expo2" {
				generate double `uHigh' = -exp(-`alpha'*((`uMax')^`r'))
				generate double `uLow' = -exp(-`alpha'*((`uMin')^`r'))
			}


	generate double `euDiff' = ((`euR' - `euL')/(`uHigh' - `uLow'))/`noise'
}
else if "$error" == "no" {
	generate double `euDiff' = `euR' - `euL'

}	


* Evaluate the likelihood 
if "$error" == "luce" {
	replace `lnf' = ln(`euDiff') if `choice'==1 
	replace `lnf' = ln(1-`euDiff') if `choice'==0
}
else { 
	replace `lnf' = ln($cdf(`euDiff')) if `choice'==1 
	replace `lnf' = ln($cdf(-`euDiff')) if `choice'==0
}

}
end


program define ml_rdu_prelec2_LN
*Remember to define the globals $error, $cdf, $ufunc and $weigh prior to estimation

*Specify the arguments of this program 
if "$ufunc" == "expo2" {
	args lnf r alpha LNphi LNeta noise
	}
else {
	args lnf r LNphi LNeta noise
	}	
	
*Declare the temporary variables to be used 
tempvar choice prize1L prize2L prize3L prize1R prize2R prize3R 
tempvar u1L u2L u3L u1R u2R u3R euL euR euDiff uMax uMin uHigh uLow
tempvar p1l p2l p3l p1r p2r p3r
tempvar pw_p1l pw_p2l pw_p3l pw_p1r pw_p2r pw_p3r
tempvar dw_p1l dw_p2l dw_p3l dw_p1r dw_p2r dw_p3r
tempvar dw_check_l dw_check_r 
tempvar phi eta

quietly {

*Initialize the data 
generate int `choice' = $ML_y1

generate double `prize1L' = $ML_y8 
generate double `prize2L' = $ML_y9 
generate double `prize3L' = $ML_y10 

generate double `prize1R' = $ML_y11 
generate double `prize2R' = $ML_y12 
generate double `prize3R' = $ML_y13 

generate double `uMax' = $ML_y14
generate double `uMin' = $ML_y15

*Impose non-negativity constraints
generate double `phi' = exp(`LNphi')
generate double `eta' = exp(`LNeta')

*Generate the aggregate probabilities for the Left lottery
generate double `p3l' = $ML_y4
generate double `p2l' = $ML_y3
replace `p2l' = `p2l' + `p3l' if `p2l' > 0
generate double `p1l' = $ML_y2
replace `p1l' = `p1l' + `p2l' if `p1l' > 0
replace `p1l' = `p1l' + `p3l' if `p1l' > 0 & `p2l' == 0

*Generate the weighted probabilities for the Left lottery using Prelec weighting functions
	forvalues i = 1/3 {
		generate double `pw_p`i'l' = `p`i'l'
		replace `pw_p`i'l' = exp((-`eta')*((-ln(`p`i'l'))^`phi')) if `p`i'l' > 0 & `p`i'l' < 1
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
	forvalues i = 1/3 {
		generate double `pw_p`i'r' = `p`i'r'
		replace `pw_p`i'r' = exp((-`eta')*((-ln(`p`i'r'))^`phi')) if `p`i'r' > 0 & `p`i'r' < 1
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
}
else if "$ufunc" == "expo2" {
*Left lottery
	generate double `u1L' = -exp(-`alpha'*((`prize1L')^`r'))
	generate double `u2L' = -exp(-`alpha'*((`prize2L')^`r'))
	generate double `u3L' = -exp(-`alpha'*((`prize3L')^`r'))

*Right lottery
	generate double `u1R' = -exp(-`alpha'*((`prize1R')^`r'))
	generate double `u2R' = -exp(-`alpha'*((`prize2R')^`r'))
	generate double `u3R' = -exp(-`alpha'*((`prize3R')^`r'))
}



*Calculate RDU of each lottery 
generate double `euL' = (`dw_p1l'*`u1L')+(`dw_p2l'*`u2L')+(`dw_p3l'*`u3L')
generate double `euR' = (`dw_p1r'*`u1R')+(`dw_p2r'*`u2R')+(`dw_p3r'*`u3R')

* Calculate the RDU difference by using Fechner, Luce or Contextual errors 
if "$error" == "fechner" {
	generate double `euDiff' = (`euR'-`euL')/(`noise')
}
else if "$error" == "luce" {
	generate double `euDiff' = (`euR'^(1/`noise'))/((`euR'^(1/`noise'))+(`euL'^(1/`noise')))
}
else if "$error" == "context" {
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
			else if "$ufunc" == "expo2" {
				generate double `uHigh' = -exp(-`alpha'*((`uMax')^`r'))
				generate double `uLow' = -exp(-`alpha'*((`uMin')^`r'))
			}


	generate double `euDiff' = ((`euR' - `euL')/(`uHigh' - `uLow'))/`noise'
}
else if "$error" == "no" {
	generate double `euDiff' = `euR' - `euL'

}	


* Evaluate the likelihood 
if "$error" == "luce" {
	replace `lnf' = ln(`euDiff') if `choice'==1 
	replace `lnf' = ln(1-`euDiff') if `choice'==0
}
else { 
	replace `lnf' = ln($cdf(`euDiff')) if `choice'==1 
	replace `lnf' = ln($cdf(-`euDiff')) if `choice'==0
}

}
end

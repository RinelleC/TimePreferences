*********************************************************************
*   DO FILE: Margins                                                *
*   Estimates margins on time preferences 							*
*********************************************************************

* Start log file 
cap log close 
log using "$logfiles/Log_Time_3_Margins.txt", text replace

*********************************************************************
***     				Exponential Discounting      			  ***
*********************************************************************

* Delta Equation
estimates restore m1hetero
asdoc margins, over(wave) predict(equation(delta)) post ///
	replace save($stata_tables/Discounting_Exponential) label dec(5) ///
	title(Delta Estimates)

*----------------------------------------------*
*  Table of PVs under Exponential Discounting  *
*----------------------------------------------*

	* R300 in 14 days 
	estimates restore m1hetero
	asdoc margins, over(wave) expression(300*(1/((1+predict(equation(delta)))^(14/365)))) ///
		append save($stata_tables/Discounting_Exponential) label dec(0) ///
		title(PV for R300)
	
	* R400 in 14 days 
	estimates restore m1hetero
	asdoc margins, over(wave) expression(400*(1/((1+predict(equation(delta)))^(14/365)))) ///
		append save($stata_tables/Discounting_Exponential) label dec(0) ///
		title(PV for R400)

	* R500 in 14 days 
	estimates restore m1hetero
	asdoc margins, over(wave) expression(500*(1/((1+predict(equation(delta)))^(14/365)))) ///
		append save($stata_tables/Discounting_Exponential) label dec(0) ///
		title(PV for R500 and 14 days) ///
		saving($estimations/pv_E_500_14days, replace) post

					* Test for wave effects (R500 and 14 days)
					foreach i in 1 2 3 4 5 6 {
						foreach j in `ferest()' {
						test `i'.wave == `j'.wave
							if r(p) < 0.05 {
								di as error r(p) 
							}
						}
					}

	* 600 in 14 days 
	estimates restore m1hetero
	asdoc margins, over(wave) expression(600*(1/((1+predict(equation(delta)))^(14/365)))) ///
		append save($stata_tables/Discounting_Exponential) label dec(0) ///
		title(PV for R600)

*********************************************************************
***   				Quasi-Hyperbolic Discounting   				  ***
*********************************************************************
    
* Beta Equation 
estimates restore m2hetero
asdoc margins, over(wave) predict(equation(beta)) post ///
	replace save($stata_tables/Discounting_QuasiHyperbolic) label dec(5) ///
	title(Beta Estimates)

* Delta Equation
estimates restore m2hetero
asdoc margins, over(wave) predict(equation(delta)) post ///
	append save($stata_tables/Discounting_QuasiHyperbolic) label dec(5) ///
	title(Delta Estimates)

local beta "(predict(equation(beta)))"

*------------------------------------------------*
*  Table of Present Values under QH Discounting  *
*------------------------------------------------*

    * R300 in 14 days 
	estimates restore m2hetero
	asdoc margins, over(wave) expression(300*`beta'*(1/((1+predict(equation(delta)))^(14/365)))) post ///
		append save($stata_tables/Discounting_QuasiHyperbolic) label dec(0) ///
		title(PV for R300 and 14 days)
		
	* R400 in 14 days 
	estimates restore m2hetero
	asdoc margins, over(wave) expression(400*`beta'*(1/((1+predict(equation(delta)))^(14/365)))) post ///
		append save($stata_tables/Discounting_QuasiHyperbolic) label dec(0) ///
		title(PV for R400 and 14 days)

	* R500 in 14 days 
	estimates restore m2hetero
	asdoc margins, over(wave) expression(500*`beta'*(1/((1+predict(equation(delta)))^(14/365)))) post ///
		append save($stata_tables/Discounting_QuasiHyperbolic) label dec(0) ///
		title(PV for R500 and 14 days) ///
	    saving($estimations/pv_QH_500_14days, replace)

						* Test for wave effects
						foreach i in 1 2 3 4 5 6 {
							foreach j in `ferest()' {
							test `i'.wave == `j'.wave
								if r(p) < 0.05 {
											di as error r(p) 
								}
							}
						}

	* 600 in 14 days 
	estimates restore m2hetero
	asdoc margins, over(wave) expression(600*`beta'*(1/((1+predict(equation(delta)))^(14/365)))) post ///
		append save($stata_tables/Discounting_QuasiHyperbolic) label dec(0) ///
		title(PV for R600 and 14 days)

*********************************************************************
***     				Hyperbolic Discounting      			  ***
*********************************************************************
    
* Delta Equation
estimates restore m3hetero
asdoc margins, over(wave) predict(equation(delta)) post ///
	replace save($stata_tables/Discounting_Hyperbolic) label dec(5) ///
	title(Delta Estimates)

*---------------------------------------------*
*  Table of PVs under Hyperbolic Discounting  *
*---------------------------------------------*

	* R300 in 14 days 
	estimates restore m3hetero
	asdoc margins, over(wave) expression(300*(1/(1+predict(equation(delta))*(14/365)))) ///
		append save($stata_tables/Discounting_Hyperbolic) label dec(0) ///
		title(PV for R300)
	
	* R400 in 14 days 
	estimates restore m3hetero
	asdoc margins, over(wave) expression(400*(1/(1+predict(equation(delta))*(14/365)))) ///
		append save($stata_tables/Discounting_Hyperbolic) label dec(0) ///
		title(PV for R400)

	* R500 in 14 days 
	estimates restore m3hetero
	asdoc margins, over(wave) expression(500*(1/(1+predict(equation(delta))*(14/365)))) ///
		append save($stata_tables/Discounting_Hyperbolic) label dec(0) ///
		title(PV for R500 and 14 days) ///
		saving($estimations/pv_H_500_14days, replace) post

					* Test for wave effects (R500 and 14 days)
					foreach i in 1 2 3 4 5 6 {
						foreach j in `ferest()' {
						test `i'.wave == `j'.wave
							if r(p) < 0.05 {
								di as error r(p) 
							}
						}
					}

	* 600 in 14 days 
	estimates restore m3hetero
	asdoc margins, over(wave) expression(600*(1/(1+predict(equation(delta))*(14/365)))) ///
		append save($stata_tables/Discounting_Hyperbolic) label dec(0) ///
		title(PV for R600)

*********************************************************************
***   					Weibull Discounting   					  ***
*********************************************************************
    
* Beta Equation 
estimates restore m4hetero
asdoc margins, over(wave) predict(equation(beta)) post ///
	replace save($stata_tables/Discounting_Weibull) label dec(5) ///
	title(Beta Estimates)

* Delta Equation
estimates restore m4hetero
asdoc margins, over(wave) predict(equation(delta)) post ///
	append save($stata_tables/Discounting_Weibull) label dec(5) ///
	title(Delta Estimates)

*-----------------------------------------------*
*  Table of Present Values under W Discounting  *
*-----------------------------------------------*

local beta "(predict(equation(beta)))"

    * R300 in 14 days 
	estimates restore m4hetero
	asdoc margins, over(wave) expression(300*exp(-predict(equation(delta))*((14/365)^(1/`beta')))) post ///
		append save($stata_tables/Discounting_Weibull) label dec(0) ///
		title(PV for R300 and 14 days)
		
	* R400 in 14 days 
	estimates restore m4hetero
	asdoc margins, over(wave) expression(400*exp(-predict(equation(delta))*((14/365)^(1/`beta')))) post ///
		append save($stata_tables/Discounting_Weibull) label dec(0) ///
		title(PV for R400 and 14 days)

	* R500 in 14 days 
	estimates restore m4hetero
	asdoc margins, over(wave) expression(500*exp(-predict(equation(delta))*((14/365)^(1/`beta')))) post ///
		append save($stata_tables/Discounting_Weibull) label dec(0) ///
		title(PV for R500 and 14 days) ///
	    saving($estimations/pv_W_500_14days, replace)

						* Test for wave effects
						foreach i in 1 2 3 4 5 6 {
							foreach j in `ferest()' {
							test `i'.wave == `j'.wave
								if r(p) < 0.05 {
											di as error r(p) 
								}
							}
						}

	* 600 in 14 days 
	estimates restore m4hetero
	asdoc margins, over(wave) expression(600*exp(-predict(equation(delta))*((14/365)^(1/`beta')))) post ///
		append save($stata_tables/Discounting_Weibull) label dec(0) ///
		title(PV for R600 and 14 days)

log close 

*********************************************************************
***   				Race Group Margins and Tests    			  ***
*********************************************************************

if "$doMARGINSRACE" == "y" {

	cap log close 
	log using "$logfiles/Log_Time_3_MarginsRace.txt", text replace

	*---------------------*
	*  Exponential        *
	*---------------------*

	* R500 in 14 days 
	estimates restore m1hetero
	asdoc margins, over(race wave) expression(500*(1/((1+predict(equation(delta)))^(14/365)))) ///
		replace save($stata_tables/Racial_Groups) label dec(0) ///
		title(Exponential - PV for R500 and 14 days) ///
		saving($explanatory/Race_Exponential, replace) post

	* Test whether the margin for race differs between two waves
	di _newline(1) "Test whether race differs between two waves"
	* Race = 1
	forvalues w1 = 1/5 {
		local w2start = `w1' + 1
		forvalues w2 = `w2start'/6 {
			test 1.race#`w1'.wave == 1.race#`w2'.wave
			}
		}
	* Race = 2
	forvalues w1 = 1/5 {
		local w2start = `w1' + 1
		forvalues w2 = `w2start'/6 {
			test 2.race#`w1'.wave == 2.race#`w2'.wave
			}
		}
	*

	*---------------------*
	*  Quasi-Hyperbolic   *
	*---------------------*

	local beta "(predict(equation(beta)))"

	* R500 in 14 days 
	estimates restore m2hetero
	asdoc margins, over(race wave) expression(600*`beta'*(1/((1+predict(equation(delta)))^(14/365)))) ///
		append save($stata_tables/Racial_Groups) label dec(0) ///
		title(Quasi-Hyperbolic - PV for R500 and 14 days) ///
		saving($explanatory/Race_QuasiHyperbolic, replace) post

	*---------------------*
	*  Hyperbolic 	      *
	*---------------------*

	* R500 in 14 days 
	estimates restore m3hetero
	asdoc margins, over(race wave) expression(500*(1/(1+predict(equation(delta))*(14/365)))) ///
		append save($stata_tables/Racial_Groups) label dec(0) ///
		title(Hyperbolic - PV for R500 and 14 days) ///
		saving($explanatory/Race_Hyperbolic, replace) post

	*---------------------*
	*  Weibull  	      *
	*---------------------*

	local beta "(predict(equation(beta)))"

	* R500 in 14 days 
	estimates restore m4hetero
	asdoc margins, over(race wave) expression(600*exp(-predict(equation(delta))*((14/365)^(1/`beta')))) ///
		append save($stata_tables/Racial_Groups) label dec(0) ///
		title(Weibull - PV for R500 and 14 days) ///
		saving($explanatory/Race_Weibull, replace) post

log close 

}

*********************************************************************

di as error "End of Margins do-file" 
*****************************************************************************
* This do file builds on the analysis of Time Preferences in South Africa.  *
* It fits within section 8 of the Main do-file. Analysis assumes RDU.       * 
* This do file only conducts margins commands and are kept separate         *
* because of how long the estimations take.                                 * 
*                                                                           *
* Date first generated:             8 September 2025                        *
* Created by:                       Rinelle Chetty                          * 
*****************************************************************************

****************************************
***     Exponential Discounting      ***
****************************************
    
* Delta Equation
estimates restore m1hetero
margins, over(wave) predict(equation(delta)) post


* Table of Present Values under Exponential Discounting 
	* R300 in 14 days 
	estimates restore m1hetero
	margins, over(wave) expression(300*(1/((1+predict(equation(delta)))^(14/365)))) 
	* R400 in 14 days 
	estimates restore m1hetero
	margins, over(wave) expression(400*(1/((1+predict(equation(delta)))^(14/365))))
	* R500 in 14 days 
	estimates restore m1hetero
	margins, over(wave) expression(500*(1/((1+predict(equation(delta)))^(14/365)))) ///
		saving($estimations/pvExp50margin, replace) post
	* 600 in 14 days 
	estimates restore m1hetero
	margins, over(wave) expression(600*(1/((1+predict(equation(delta)))^(14/365))))


****************************************
***   Quasi-Hyperbolic Discounting   ***
****************************************
    
* Beta Equation 
estimates restore m3hetero
margins, over(wave) predict(equation(beta)) post

* Delta Equation
estimates restore m3hetero
margins, over(wave) predict(equation(delta)) post

estimates restore m3hetero

local beta "(predict(equation(beta)))"


* Table of Present Values under Quasi-Hyperbolic Discounting 
    * R300 in 14 days 
	estimates restore m1hetero
	margins, over(wave) expression(300*`beta'*(1/((1+predict(equation(delta)))^(14/365))))
	* R400 in 14 days 
	estimates restore m1hetero
	margins, over(wave) expression(400*`beta'*(1/((1+predict(equation(delta)))^(14/365))))
	* R500 in 14 days 
	estimates restore m1hetero
	margins, over(wave) expression(500*`beta'*(1/((1+predict(equation(delta)))^(14/365)))) ///
		saving($estimations/pvQH50margin, replace) post
	* 600 in 14 days 
	estimates restore m1hetero
	margins, over(wave) expression(600*`beta'*(1/((1+predict(equation(delta)))^(14/365))))


*******************************************************************************

di as error "End of Margins do-file" 
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
asdoc margins, over(wave) predict(equation(delta)) post ///
	replace save($stata_tables/Exp_PresentValues) label dec(5) ///
	title(Delta Estimates)

	*********************************************************
	* Table of Present Values under Exponential Discounting *
	*********************************************************

	* R300 in 14 days 
	estimates restore m1hetero
	asdoc margins, over(wave) expression(300*(1/((1+predict(equation(delta)))^(14/365)))) ///
		append save($stata_tables/Exp_PresentValues) label dec(2) ///
		title(PV for R300)
	
	* R400 in 14 days 
	estimates restore m1hetero
	asdoc margins, over(wave) expression(400*(1/((1+predict(equation(delta)))^(14/365)))) ///
		append save($stata_tables/Exp_PresentValues) label dec(2) ///
		title(PV for R400)

	* R500 in 14 days 
	estimates restore m1hetero
	asdoc margins, over(wave) expression(500*(1/((1+predict(equation(delta)))^(14/365)))) ///
		append save($stata_tables/Exp_PresentValues) label dec(2) ///
		title(PV for R500)
	
	* 600 in 14 days 
	estimates restore m1hetero
	asdoc margins, over(wave) expression(600*(1/((1+predict(equation(delta)))^(14/365)))) ///
		append save($stata_tables/Exp_PresentValues) label dec(2) ///
		title(PV for R600)

	****************************************************************
	* Comparing R500 Present Values across time horizons under EXP *
	****************************************************************

    * 7 days 
    estimates restore m1hetero
	asdoc margins, over(wave) expression(500*(1/((1+predict(equation(delta)))^(7/365))))  ///
		append save($stata_tables/Exp_PresentValues) label dec(2) ///
		title(PV for R500 and 7 days) /// 
		saving($estimations/pv_E_500_7days, replace) post
	
	* 14 days - conducted in AnalysisTime.do
   	estimates restore m1hetero
	asdoc margins, over(wave) expression(500*(1/((1+predict(equation(delta)))^(14/365))))  ///
		append save($stata_tables/Exp_PresentValues) label dec(2) ///
		title(PV for R500 and 14 days) ///     	
		saving($estimations/pv_E_500_14days, replace) post

    * 42 days 
    estimates restore m1hetero
	asdoc margins, over(wave) expression(500*(1/((1+predict(equation(delta)))^42/365))))  ///
		append save($stata_tables/Exp_PresentValues) label dec(2) ///
		title(PV for R500 and 42 days) ///         
		saving($estimations/pv_E_500_42days, replace) post

    * 48 days 
   	estimates restore m1hetero
	asdoc margins, over(wave) expression(500*(1/((1+predict(equation(delta)))^(48/365))))  ///
		append save($stata_tables/Exp_PresentValues) label dec(2) ///
		title(PV for R500 and 48 days) ///         
		saving($estimations/pv_E_500_48days, replace) post



****************************************
***   Quasi-Hyperbolic Discounting   ***
****************************************
    
* Beta Equation 
estimates restore m3hetero
asdoc margins, over(wave) predict(equation(beta)) post ///
	replace save($stata_tables/QH_PresentValues) label dec(5) ///
	title(Beta Estimates)

* Delta Equation
estimates restore m3hetero
asdoc margins, over(wave) predict(equation(delta)) post ///
	append save($stata_tables/QH_PresentValues) label dec(5) ///
	title(Delta Estimates)

local beta "(predict(equation(beta)))"

	**************************************************************
	* Table of Present Values under Quasi-Hyperbolic Discounting *
	**************************************************************

    * R300 in 14 days 
	estimates restore m3hetero
	asdoc margins, over(wave) expression(300*`beta'*(1/((1+predict(equation(delta)))^(14/365)))) post ///
		append save($stata_tables/QH_PresentValues) label dec(2) ///
		title(PV for R300 and 14 days)
		
	* R400 in 14 days 
	estimates restore m3hetero
	asdoc margins, over(wave) expression(400*`beta'*(1/((1+predict(equation(delta)))^(14/365)))) post ///
		append save($stata_tables/QH_PresentValues) label dec(2) ///
		title(PV for R400 and 14 days)

	* R500 in 14 days - conducted in AnalysisTime.do 
	estimates restore m3hetero
	asdoc margins, over(wave) expression(500*`beta'*(1/((1+predict(equation(delta)))^(14/365)))) post ///
		append save($stata_tables/QH_PresentValues) label dec(2) ///
		title(PV for R500 and 14 days)
	
	* 600 in 14 days 
	estimates restore m3hetero
	asdoc margins, over(wave) expression(600*`beta'*(1/((1+predict(equation(delta)))^(14/365)))) post ///
		append save($stata_tables/QH_PresentValues) label dec(2) ///
		title(PV for R600 and 14 days)

	***************************************************************
	* Comparing R500 Present Values across time horizons under QH *
	***************************************************************

    * 7 days 
    estimates restore m3hetero
	asdoc margins, over(wave) expression(500*`beta'*(1/((1+predict(equation(delta)))^(7/365)))) ///  
        append save($stata_tables/QH_PresentValues) label dec(2) ///
		title(PV for R500 and 7 days) ///
		saving($estimations/pv_QH_500_7days, replace) post
	
	* 14 days 
    estimates restore m3hetero
	asdoc margins, over(wave) expression(500*`beta'*(1/((1+predict(equation(delta)))^(14/365)))) ///  
        append save($stata_tables/QH_PresentValues) label dec(2) ///
		title(PV for R500 and 14 days) ///
        saving($estimations/pv_QH_500_14days, replace) post

    * 42 days 
    estimates restore m3hetero
	asdoc margins, over(wave) expression(500*`beta'*(1/((1+predict(equation(delta)))^(42/365)))) ///  
        append save($stata_tables/QH_PresentValues) label dec(2) ///
		title(PV for R500 and 42 days) ///
        saving($estimations/pv_QH_500_42days, replace) post

    * 48 days 
    estimates restore m3hetero
	asdoc margins, over(wave) expression(500*`beta'*(1/((1+predict(equation(delta)))^(48/365)))) ///  
        append save($stata_tables/QH_PresentValues) label dec(2) ///
		title(PV for R500 and 48 days) ///
		saving($estimations/pv_QH_500_48days, replace) post


*******************************************************************************

di as error "End of Margins do-file" 
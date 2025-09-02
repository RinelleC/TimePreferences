*****************************************************************************
* This is the Analysis do file to analyse Time Preferences in South Africa  *
* and fits within section 8 of the Main do-file. Analysis assumes RDU.      * 
*                                                                           *
* Date first generated:             28 August 2025                          *
* Created by:                       Rinelle Chetty                          * 
*****************************************************************************


*******************************************************************************
*** 	8.1 -- Homogenous Preferences               						***
*******************************************************************************

estimates clear
set more off

global cdf          "invlogit"
global maxtech      "nr"
global riskvars     "prob1L prob2L prob3L prob1R prob2R prob3R prize1L prize2L prize3L prize1R prize2R prize3R uMax uMin"
global timevars     "risk ssamount ssdelay llamount lldelay"
global demog        ""
global hetero       ""

global error "context"      // contextual error 
global weigh "prelec2"      // prelec2 weighting function 

* Set to "crra" if you only want CRRA functions with the "1-r" form; otherwise set to "power"
// global doCRRA "power"
global doCRRA "crra"

if "$doCRRA" == "power" {

	global ufunc "power"

    ***********************************************************
    ***       Power utility, exponential discounting        ***
    ***********************************************************

	set more off
	global discount "exp"
	ml model lf ml_rdu_discount_flex (r: choice $riskvars $timevars = $demog) ///
	(phi: $demog) (eta: $demog) (delta: $demog) ///
	(noiseRA: $hetero) (noiseDR: $hetero) if risk == 1 | time == 1, ///
	cluster(id) technique($maxtech) init(0.595 0.611 0.925 0.818 0.146 4.315, copy)
	ml maximize, difficult

	estimates store m1, title(Model - Prelec2Exp)

	esttab m1 using "$stata_tables/ml_model_homogenous.rtf" , replace ///
		label se b(%15.10g) ///
        mtitle("Homogenous Preferences A") ///
        title(Power Utility & Exponential Discounting)

    ***********************************************************
    ***     Power utility, quasi-hyperbolic discounting     ***
    ***********************************************************

	set more off
	global discount "qh"
	ml model lf ml_rdu_discount_flex (r: choice $riskvars $timevars = $demog) ///
	(phi: $demog) (eta: $demog) (beta: $demog) (delta: $demog) ///
    (noiseRA: $hetero) (noiseDR: $hetero) if risk == 1 | time == 1, ///
	cluster(id) technique($maxtech) init(0.596 0.611 0.925 0.978 0.733 0.146 4.094, copy)
	ml maximize, difficult

	estimates store m3, title(Model 3 - Prelec2QHyp)

	esttab m3 using "$stata_tables/ml_model_homogenous.rtf" , append ///
		label se b(%15.10g) ///
        mtitle("Homogenous Preferences B") ///
        title(Power Utility & Quasi-Hyperbolic Discounting)

	test [beta]_cons == 1
}

else if "$doCRRA" == "crra" {

	global ufunc "crra"

    ***********************************************************
    ***        CRRA utility, exponential discounting        ***
    ***********************************************************

	set more off
	global discount "exp"
	ml model lf ml_rdu_discount_flex (r: choice $riskvars $timevars = $demog) ///
	(phi: $demog) (eta: $demog) (delta: $demog) ///
	(noiseRA: $hetero) (noiseDR: $hetero) if risk == 1 | time == 1, ///
	cluster(id) technique($maxtech) init(0.4041356 0.6118054 0.9252319 2.268743 0.0372434 1.836531, copy)
	ml maximize, difficult

	estimates store m1, title(Model - Prelec2Exp)

	esttab m1 using "$stata_tables/ml_model_homogenous.rtf" , replace ///
		label se b(%15.10g) ///
        mtitle("Homogenous Preferences A") ///
        title(CRRA Utility & Exponential Discounting)

    ***********************************************************
    ***     CRRA utility, quasi-hyperbolic discounting      ***
    ***********************************************************
    
	set more off
	global discount "qh"
	ml model lf ml_rdu_discount_flex (r: choice $riskvars $timevars = $demog) ///
	(phi: $demog) (eta: $demog) (beta: $demog) (delta: $demog) ///
	(noiseRA: $hetero) (noiseDR: $hetero) if risk == 1 | time == 1, ///
	cluster(id) technique($maxtech) init(0.4034 0.6116 0.9257 0.9789 2.0822 0.0371 1.7374, copy)
	ml maximize, difficult

	estimates store m3, title(Model 3 - Prelec2QHyp)

	esttab m3 using "$stata_tables/ml_model_homogenous.rtf" , append ///
		label se b(%15.10g) ///
        mtitle("Homogenous Preferences B") ///
        title(CRRA Utility & Quasi-Hyperbolic Discounting)

	test [beta]_cons == 1
}


    * Evaluate QH Discount rate at different horizons, specified in days
    foreach x of numlist 7 14 42 84 {
        di as error "QH discount rate evaluated at `x' day horizon"
        nlcom (QHDiscountRate : ([beta]_cons/(1+[delta]_cons)^(`x'/365))^(-1/(`x'/365)) - 1)
    }	


* Homogeneous preferences across waves to test whether all waves are best characterised by QH
forvalues w = 1/6 {
	estimates restore m3
	di "" ""	
	di as error "Estimates for Wave #`w'"
	global discount "qh"
		
    ml model lf ml_rdu_discount_flex (r: choice $riskvars $timevars = $demog) ///
        (phi: $demog) (eta: $demog) (beta: $demog) (delta: $demog) ///
	    (noiseRA: $hetero) (noiseDR: $hetero) if risk == 1 | time == 1, ///
	    cluster(id) technique($maxtech) continue
	ml maximize, difficult

	* test for QH
	di as error "Test for QH in wave #`w'"
	test [beta]_cons == 1
}


* Export the estimates to .TSV for easy import in Excel
estout * using "$estimations/RDUDiscEstimates_Homogenous.tsv", ///
replace starlevels(* 0.10 ** 0.05 *** 0.01) cells( b(star label("Estimate") fmt(3)) se(label("Std Error") par(`"="("' `")""') fmt(3))  ) ///
stats(N ll, fmt(%5.0f %10.3f) labels(N "log-likelihood")) nobaselevels ///
varlabels(r:_cons "CRRA function parameter (r)" phi:_cons "PWF parameter (phi)" eta:_cons "PWF parameter (eta)" beta:_cons "Discounting parameter (beta)" delta:_cons "Discounting parameter (delta)" noiseRA:_cons "Risk error (mu)" noiseDR:_cons "Time error (nu)") ///
prehead("Table" "Discounting Function ML Estimates" @title) ///
title("Concave Utility, Homogenous Preferences") ///
legend eqlabels(none) mlabels(,titles) /// ///
postfoot("Results account for clustering at the individual level" "Standard errors in parentheses")
	


*******************************************************************************
*** 	8.2 -- Heterogenous Preferences               						***
*******************************************************************************

set more off
global cdf "invlogit"
global maxtech "nr"
global riskvars "prob1L prob2L prob3L prob1R prob2R prob3R prize1L prize2L prize3L prize1R prize2R prize3R uMax uMin"
global timevars "risk ssamount ssdelay llamount lldelay"
global hetero ""

global error "context"
global weigh "prelec2"


if "$doCRRA" == "power" {

	global ufunc "power"
    
    ***********************************************************
    ***       Power utility, exponential discounting        ***
    ***********************************************************

	set more off
	global discount "exp"

	forvalues l = 1/2 {
		
        if `l' == 1 { 
			estimates restore m1
			global demog "i.wave c.age i.male c.anxiety_total c.depression_total i.race i.race#i.wave i.male#i.wave"
			ml model lf ml_rdu_discount_flex (r: choice $riskvars $timevars = $demog) ///
			(phi: $demog) (eta: $demog) (delta: ) ///
			(noiseRA: $hetero) (noiseDR: $hetero) if risk == 1 | time == 1, ///
			cluster(id) technique($maxtech) continue
		}
		
        if `l' == 2 { 
			estimates restore m1hetero
			ml model lf ml_rdu_discount_flex (r: choice $riskvars $timevars = $demog) ///
			(phi: $demog) (eta: $demog) (delta: $demog) ///
			(noiseRA: $hetero) (noiseDR: $hetero) if risk == 1 | time == 1, ///
			cluster(id) technique($maxtech) continue
		}
		
        ml maximize, difficult
		estimates store m1hetero

		esttab m1hetero using "$stata_tables/ml_model_heterogenous.rtf" , replace ///
		label se b(%15.10g) ///
        mtitle("Heterogenous Preferences A") ///
        title(Power Utility & Exponential Discounting)

	}


    ***********************************************************
    ***     Power utility, quasi-hyperbolic discounting     ***
    ***********************************************************

	set more off
	global discount "qh"

	ml model lf ml_rdu_discount_flex (r: choice $riskvars $timevars = $demog) ///
	(phi: $demog) (eta: $demog) (beta: $demog) (delta: $demog) ///
	(noiseRA: $hetero) (noiseDR: $hetero) if risk == 1 | time == 1, ///
	cluster(id) technique($maxtech) continue
	ml maximize, difficult

	estimates store m3hetero, title(Model 3 - Prelec2QHyp)

	esttab m1hetero using "$stata_tables/ml_model_heterogenous.rtf" , append ///
		label se b(%15.10g) ///
        mtitle("Heterogenous Preferences B") ///
        title(Power Utility & Quasi-Hyperbolic Discounting)

	test [beta]_cons == 1
	
}

else if "$doCRRA" == "crra" {

	global ufunc "crra"

    ***********************************************************
    ***        CRRA utility, exponential discounting        ***
    ***********************************************************

	set more off
	global discount "exp"

	forvalues l = 1/2 {
	
    	if `l' == 1 { 
			estimates restore m1
			global demog "i.wave c.age i.male c.anxiety_total c.depression_total i.race i.race#i.wave i.male#i.wave"
			ml model lf ml_rdu_discount_flex (r: choice $riskvars $timevars = $demog) ///
			(phi: $demog) (eta: $demog) (delta: ) ///
			(noiseRA: $hetero) (noiseDR: $hetero) if risk == 1 | time == 1, ///
			cluster(id) technique($maxtech) continue
		}

		if `l' == 2 { 
			estimates restore m1hetero
			ml model lf ml_rdu_discount_flex (r: choice $riskvars $timevars = $demog) ///
			(phi: $demog) (eta: $demog) (delta: $demog) ///
			(noiseRA: $hetero) (noiseDR: $hetero) if risk == 1 | time == 1, ///
			cluster(id) technique($maxtech) continue
		}

		ml maximize, difficult tolerance(1e-04) ltolerance(0) nrtolerance(1e-05)
		
		estimates store m1hetero
	
		esttab m1hetero using "$stata_tables/ml_model_heterogenous.rtf" , replace ///
		label se b(%15.10g) ///
        mtitle("Heterogenous Preferences A") ///
        title(CRRA Utility & Exponential Discounting)

    }

	estimates save "$estimations/time_het_exp", replace


    ***********************************************************
    ***     CRRA utility, quasi-hyperbolic discounting      ***
    ***********************************************************
    
	set more off
	global discount "qh"
	estimates restore m3

	global demog "i.wave c.age i.male c.anxiety_total c.depression_total i.race i.race#i.wave i.male#i.wave"

	ml model lf ml_rdu_discount_flex (r: choice $riskvars $timevars = $demog) ///
	(phi: $demog) (eta: $demog) (beta: $demog) (delta: $demog) ///
	(noiseRA: $hetero) (noiseDR: $hetero) if risk == 1 | time == 1, ///
	cluster(id) technique($maxtech) continue
	ml maximize, difficult		

	estimates store m3hetero, title(Model 3 - Prelec2QHyp)
	
	esttab m3hetero using "$stata_tables/ml_model_heterogenous.rtf" , append ///
		label se b(%15.10g) ///
        mtitle("Heterogenous Preferences B") ///
        title(CRRA Utility & Quasi-Hyperbolic Discounting)

	estimates save "$estimations/time_het_qh", replace

	test [beta]_cons == 1
}
		

* Now estimate with simpler demographics and covid_scale
* Specify stripped down demographics
global demogX "i.wave c.age i.male i.race c.anxiety_total c.depression_total covid_scale_deaths covid_scale_deaths_sq"

estimates restore m3
ml model lf ml_rdu_discount_flex (r: choice $riskvars $timevars = $demogX) ///
    (phi: $demogX) (eta: $demogX) (beta: $demogX) (delta: $demogX) ///
    (noiseRA: $hetero) (noiseDR: $hetero) if risk == 1 | time == 1, ///
    cluster(id) technique($maxtech) continue
ml maximize, difficult

estimates store mX 

esttab mX using "$stata_tables/ml_model_coviddeaths.rtf" , replace ///
		label se b(%15.10g) ///
        mtitle("Heterogenous Preferences with Covid Deaths") ///
        title(? Utility & ? Discounting)

estimates save "$estimations/time_covid_scale_deaths", replace
test covid_scale_deaths covid_scale_deaths_sq, mtest(noadjust)



*******************************************************************************
*** 	8.3 -- Margins for Exponential & Quasi-Hyperbolic Discounting       ***
*******************************************************************************

    ****************************************
    ***     Exponential Discounting      ***
    ****************************************
    
    * Delta Equation
    estimates restore m1hetero
    margins, over(wave) predict(equation(delta)) post

    * Test for wave effects
    foreach i in 1 2 3 4 5 6 {
        foreach j in `ferest()' {
        test `i'.wave == `j'.wave
            if r(p) < 0.1 {
                di as error r(p) 
            }
        }
    }

    * Estimate present value of $50 margins for comparisons across waves and downstream figures
    estimates restore m1hetero
    margins, over(wave) expression(50*(1/((1+predict(equation(delta)))^(14/365)))) ///
    saving($estimations/pvExp50margin, replace) post

    * Test for wave effects
    foreach i in 1 2 3 4 5 6 {
        foreach j in `ferest()' {
        test `i'.wave == `j'.wave
            if r(p) < 0.1 {
                di as error r(p) 
            }
        }
    }


    ****************************************
    ***   Quasi-Hyperbolic Discounting   ***
    ****************************************
    
    * Beta Equation 
	estimates restore m3hetero
	margins, over(wave) predict(equation(beta)) post

    * Test for wave effects
    foreach i in 1 2 3 4 5 6 {
        foreach j in `ferest()' {
        test `i'.wave == `j'.wave
            if r(p) < 0.1 {
                di as error r(p) 
            }
        }
    }

    * Delta Equation
    estimates restore m3hetero
    margins, over(wave) predict(equation(delta)) post

    * Test for wave effects
    foreach i in 1 2 3 4 5 6 {
        foreach j in `ferest()' {
        test `i'.wave == `j'.wave
            if r(p) < 0.1 {
                di as error r(p) 
            }
        }
    }

    * Estimate present value of $50 margins for comparisons across waves and downstream figures
    estimates restore m3hetero

    local beta "(predict(equation(beta)))"

    margins, over(wave) expression(50*`beta'*(1/((1+predict(equation(delta)))^(14/365)))) ///
    saving($estimations/pvQH50margin, replace) post

    * Test for wave effects
    foreach i in 1 2 3 4 5 6 {
        foreach j in `ferest()' {
        test `i'.wave == `j'.wave
            if r(p) < 0.1 {
                di as error r(p) 
            }
        }
    }


* Export the estimates to .TSV 
estout m1hetero m3hetero using "$estimations/RDUDiscEstimates_Heterogenous.tsv", ///
replace starlevels(* 0.10 ** 0.05 *** 0.01) cells( (b(star label("Estimate") fmt(3)) se(label("Std error") fmt(3)) )) ///
stats(N ll, fmt(%5.0f %10.3f) labels(N "log-likelihood")) nobaselevels ///
varlabels("$varlabels") ///
prehead("Table" "Discounting Function ML Estimates" @title) ///
title("Concave Utility, Heterogenous(?) Preferences") ///
legend eqlabels("CRRA function parameter (r)" "PWF parameter (phi)" "PWF parameter (eta)" "Discounting parameter (delta)" "Risk error (mu)" "Time error (nu)" "Discounting parameter (beta)") mlabels(,titles) ///
postfoot("Results account for clustering at the individual level" "Standard errors in parentheses")


*******************************************************************************
*** 	8.4 -- Getting metrics for the graphs						        ***
*******************************************************************************

tab ssamount, m 			// R250 and R400 principal amounts 
tab llamount, m 

* Mean of LL rewards across all time choices is R394.11
su llamount if time == 1 	

* Mean of LL rewards for the principal of R250 is R307.71
su llamount if ssamount == 250
su lldelay if ssamount == 250 

* Mean of LL rewards for the principal of R400 is R479.81
su llamount if ssamount == 400
	*---> use R500 as the representative amount for graphs 
su lldelay  if ssamount == 400 
tab lldelay if ssamount == 400 


*******************************************************************************

di as error "End of Analysis do-file" 
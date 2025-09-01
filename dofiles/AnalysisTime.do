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

    ***********************************************************
    ***     CRRA utility, quasi-hyperbolic discounting     ***
    ***********************************************************
    
	set more off
	global discount "qh"
	ml model lf ml_rdu_discount_flex (r: choice $riskvars $timevars = $demog) ///
	(phi: $demog) (eta: $demog) (beta: $demog) (delta: $demog) ///
	(noiseRA: $hetero) (noiseDR: $hetero) if risk == 1 | time == 1, ///
	cluster(id) technique($maxtech) init(0.4034 0.6116 0.9257 0.9789 2.0822 0.0371 1.7374, copy)
	ml maximize, difficult

	estimates store m3, title(Model 3 - Prelec2QHyp)

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
	estout * using "$estimations/RDUDiscEstimates.tsv", ///
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





*******************************************************************************

di as error "End of Analysis do-file" 
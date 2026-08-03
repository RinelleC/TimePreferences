*****************************************************************************
* This is the Analysis do file to analyse Time Preferences in South Africa  *
* and fits within section 7 of the Main do-file. Analysis assumes RDU.      * 
*                                                                           *
* Date first generated:             28 August 2025                          *
* Created by:                       Rinelle Chetty                          * 
*****************************************************************************


*******************************************************************************
*** 	7.1 -- Homogenous Preferences               						***
*******************************************************************************

estimates clear
set more off

global cdf          "invlogit"
global riskvars     "prob1L prob2L prob3L prob1R prob2R prob3R prize1L prize2L prize3L prize1R prize2R prize3R uMax uMin"
global timevars     "risk ssamount ssdelay llamount lldelay"
global demog        ""
global hetero       ""

global error "context"      // contextual error 
global weigh "prelec2"      // prelec2 weighting function 
global ufunc "crra"

    *-------------------------------------------*
    *       Exponential Discounting             *
    *-------------------------------------------*

	global discount "exp"
	
    ml model lf ml_rdu_discount_flex (r: choice $riskvars $timevars = $demog) ///
        (phi: $demog) (eta: $demog) (delta: $demog) ///
        (noiseRA: $hetero) (noiseDR: $hetero) if risk == 1 | time == 1, ///
        cluster(id) technique(nr) init(0.351 0.489 0.8839 1.48568 0.1423 6.268, copy)
	ml maximize, difficult

	estimates store m1

	esttab m1 using "$stata_tables/ml_model_homogenous.rtf" , replace ///
		label se b(%15.3g) mtitle("Homogenous Preferences A") ///
        title(Exponential Discounting)

    * Evaluate EXP Discount rate at different horizons, specified in days
    foreach x of numlist 7 14 42 84 {
        di as error "EXP discount rate evaluated at `x' day horizon"
        nlcom (EXPDiscountRate : (1/(1+[delta]_cons)^(`x'/365))^(-1/(`x'/365)) - 1)
    }	

    *-------------------------------------------*
    *       Quasi-Hyperbolic Discounting        *
    *-------------------------------------------*
    
	global discount "qh"
	
    ml model lf ml_rdu_discount_flex (r: choice $riskvars $timevars = $demog) ///
	    (phi: $demog) (eta: $demog) (beta: $demog) (delta: $demog) ///
	    (noiseRA: $hetero) (noiseDR: $hetero) if risk == 1 | time == 1, ///
	    cluster(id) technique(nr) init(0.351 0.489 0.8839 0.978 1.48568 0.1423 6.268, copy)
	ml maximize, difficult

	estimates store m2

	esttab m2 using "$stata_tables/ml_model_homogenous.rtf" , append ///
		label se b(%15.3g) mtitle("Homogenous Preferences B") ///
        title(Quasi-Hyperbolic Discounting)

	test [beta]_cons == 1

    * Evaluate QH Discount rate at different horizons, specified in days
    foreach x of numlist 7 14 42 84 {
        di as error "QH discount rate evaluated at `x' day horizon"
        nlcom (QHDiscountRate : ([beta]_cons/(1+[delta]_cons)^(`x'/365))^(-1/(`x'/365)) - 1)
    }	

* Homogeneous preferences across waves to test whether all waves are best characterised by QH
forvalues w = 1/6 {
	estimates restore m2
	di "" ""	
	di as error "Estimates for Wave #`w'"
	global discount "qh"
		
    ml model lf ml_rdu_discount_flex (r: choice $riskvars $timevars = $demog) ///
        (phi: $demog) (eta: $demog) (beta: $demog) (delta: $demog) ///
	    (noiseRA: $hetero) (noiseDR: $hetero) if risk == 1 | time == 1, ///
	    cluster(id) technique(nr) continue
	ml maximize, difficult

	* test for QH
	di as error "Test for QH in wave #`w'"
	test [beta]_cons == 1
    }
	
    *-------------------------------------------*
    *       Hyperbolic Discounting              *
    *-------------------------------------------*

	global discount "mazur"
	
    ml model lf ml_rdu_discount_flex (r: choice $riskvars $timevars = $demog) ///
	    (phi: $demog) (eta: $demog) (delta: $demog) ///
	    (noiseRA: $hetero) (noiseDR: $hetero) if risk == 1 | time == 1, ///
	    cluster(id) technique(nr) init(0.351 0.489 0.8839 1.48568 0.1423 6.268, copy)
	ml maximize, difficult

	estimates store m3

	esttab m3 using "$stata_tables/ml_model_homogenous.rtf" , append ///
		label se b(%15.3g) mtitle("Homogenous Preferences C") ///
        title(Hyperbolic Discounting)


    * Evaluate Discount rate at different horizons, specified in days
    foreach x of numlist 7 14 42 84 {
        di as error "Hyperbolic discount rate evaluated at `x' day horizon"
        nlcom (EXPDiscountRate : (1/(1+[delta]_cons)^(`x'/365))^(-1/(`x'/365)) - 1)
    }	

    *-------------------------------------------*
    *       Weibull Discounting                 *
    *-------------------------------------------*

	global discount "weibull"

    ml model lf ml_rdu_discount_flex (r: choice $riskvars $timevars = $demog) ///
	    (phi: $demog) (eta: $demog) (beta: $demog) (delta: $demog) ///
	    (noiseRA: $hetero) (noiseDR: $hetero) if risk == 1 | time == 1, ///
	    cluster(id) technique(nr) init(0.351 0.489 0.8839 1 1.48568 0.1423 6.268, copy)
	ml maximize, difficult

	estimates store m4

	esttab m4 using "$stata_tables/ml_model_homogenous.rtf" , append ///
		label se b(%15.3g) mtitle("Homogenous Preferences D") ///
        title(Weibull Discounting)

    *------------------------------------------------------*
    *    Export all homogenous results to one TSV table    *
    *------------------------------------------------------*
    
estout m1 m2 m3 m4 using "$estimations/Allmodels_Homogenous.tsv", replace ///
    starlevels(* 0.10 ** 0.05 *** 0.01) ///
    cells( b(star label("Estimate") fmt(3)) se(label("Std Error") fmt(3)) ) ///
    stats(N ll, fmt(%5.0f %10.3f) labels(N "log-likelihood")) nobaselevels ///
    varlabels(r:_cons "CRRA function parameter (r)" phi:_cons "PWF parameter (phi)" eta:_cons "PWF parameter (eta)" beta:_cons "Discounting parameter (beta)" delta:_cons "Discounting parameter (delta)" noiseRA:_cons "Risk error (mu)" noiseDR:_cons "Time error (nu)") ///
    prehead("Table" "Discounting Function ML Estimates" @title) ///
    title("Concave Utility, Homogenous Preferences") ///
    legend eqlabels(none) mlabels(,titles) /// ///
    postfoot("Results account for clustering at the individual level" "Standard errors in parentheses")
            

*******************************************************************************
*** 	7.2 -- Heterogenous Preferences               						***
*******************************************************************************

set more off

global cdf "invlogit"
global riskvars "prob1L prob2L prob3L prob1R prob2R prob3R prize1L prize2L prize3L prize1R prize2R prize3R uMax uMin"
global timevars "risk ssamount ssdelay llamount lldelay"
global hetero ""

global error "context"
global weigh "prelec2"
global ufunc "crra"

    *-------------------------------------------*
    *       Exponential Discounting             *
    *-------------------------------------------*

	global discount "exp"

	forvalues l = 1/2 {
	
    	if `l' == 1 { 
			estimates restore m1
			global demog "i.wave c.age i.male c.anxiety_total c.depression_total i.race i.race#i.wave i.male#i.wave"
			ml model lf ml_rdu_discount_flex (r: choice $riskvars $timevars = $demog) ///
			(phi: $demog) (eta: $demog) (delta: ) ///
			(noiseRA: $hetero) (noiseDR: $hetero) if risk == 1 | time == 1, ///
			cluster(id) technique(nr) continue
		}

		if `l' == 2 { 
			estimates restore m1hetero
			ml model lf ml_rdu_discount_flex (r: choice $riskvars $timevars = $demog) ///
			(phi: $demog) (eta: $demog) (delta: $demog) ///
			(noiseRA: $hetero) (noiseDR: $hetero) if risk == 1 | time == 1, ///
			cluster(id) technique(nr) continue
		}

		ml maximize, difficult tolerance(1e-04) ltolerance(0) nrtolerance(1e-05)
		
		estimates store m1hetero
	
		esttab m1hetero using "$stata_tables/ml_model_heterogenous.rtf" , replace ///
		label se b(%15.3g) mtitle("Heterogenous Preferences A") ///
        title(Exponential Discounting)

    }

    *-------------------------------------------*
    *       Quasi-Hyperbolic Discounting        *
    *-------------------------------------------*
    
	global discount "qh"

	estimates restore m2

	global demog "i.wave c.age i.male c.anxiety_total c.depression_total i.race i.race#i.wave i.male#i.wave"

	ml model lf ml_rdu_discount_flex (r: choice $riskvars $timevars = $demog) ///
        (phi: $demog) (eta: $demog) (beta: $demog) (delta: $demog) ///
        (noiseRA: $hetero) (noiseDR: $hetero) if risk == 1 | time == 1, ///
        cluster(id) technique(nr) continue
	ml maximize, difficult		

	estimates store m2hetero
	
	esttab m2hetero using "$stata_tables/ml_model_heterogenous.rtf" , append ///
		label se b(%15.3g) mtitle("Heterogenous Preferences B") ///
        title(Quasi-Hyperbolic Discounting)

	test [beta]_cons == 1
		

* Now estimate with simpler demographics and covid_scale
* Specify stripped down demographics
global demogX "i.wave c.age i.male i.race c.anxiety_total c.depression_total covid_scale_deaths covid_scale_deaths_sq"

estimates restore m2
ml model lf ml_rdu_discount_flex (r: choice $riskvars $timevars = $demogX) ///
    (phi: $demogX) (eta: $demogX) (beta: $demogX) (delta: $demogX) ///
    (noiseRA: $hetero) (noiseDR: $hetero) if risk == 1 | time == 1, ///
    cluster(id) technique(nr) continue
ml maximize, difficult

estimates store m2X 

esttab m2X using "$stata_tables/ml_model_coviddeaths.rtf" , replace ///
		label se b(%15.3g) mtitle("Heterogenous Preferences with Covid Deaths") ///
        title(Quasi-Hyperbolic Discounting)

test covid_scale_deaths covid_scale_deaths_sq, mtest(noadjust)

    *-------------------------------------------*
    *       Hyperbolic Discounting              *
    *-------------------------------------------*
    
	global discount "mazur"

	estimates restore m3

	global demog "i.wave c.age i.male c.anxiety_total c.depression_total i.race i.race#i.wave i.male#i.wave"

	ml model lf ml_rdu_discount_flex (r: choice $riskvars $timevars = $demog) ///
        (phi: $demog) (eta: $demog) (delta: $demog) ///
        (noiseRA: $hetero) (noiseDR: $hetero) if risk == 1 | time == 1, ///
        cluster(id) technique(nr) continue
	ml maximize, difficult		

	estimates store m3hetero
	
	esttab m3hetero using "$stata_tables/ml_model_heterogenous.rtf" , append ///
		label se b(%15.3g) mtitle("Heterogenous Preferences C") ///
        title(Hyperbolic Discounting)

    *-------------------------------------------*
    *       Weibull Discounting                 *
    *-------------------------------------------*
    
	global discount "weibull"

	global demog    "i.wave c.age i.male c.anxiety_total c.depression_total i.race i.race#i.wave i.male#i.wave"

	forvalues l = 1/2 {

		if `l' == 1 {
			estimates restore m4
			ml model lf ml_rdu_discount_flex (r: choice $riskvars $timevars = $demog) ///
			(phi: $demog) (eta: $demog) (beta: ) (delta: ) ///
			(noiseRA: $hetero) (noiseDR: $hetero) if risk == 1 | time == 1, ///
			cluster(id) technique(bhhh 10 nr 100) continue
		    }

		if `l' == 2 {
			estimates restore m4hetero
			ml model lf ml_rdu_discount_flex (r: choice $riskvars $timevars = $demog) ///
			(phi: $demog) (eta: $demog) (beta: $demog) (delta: $demog) ///
			(noiseRA: $hetero) (noiseDR: $hetero) if risk == 1 | time == 1, ///
			cluster(id) technique(bhhh 10 nr 100) continue
		    }

		ml maximize, difficult iterate(200) tolerance(1e-04) ltolerance(0) nrtolerance(1e-05)

		estimates store m4hetero

		esttab m4hetero using "$stata_tables/ml_model_heterogenous.rtf" , append ///
		label se b(%15.3g) mtitle("Heterogenous Preferences D") ///
        title(Weibull Discounting)

	    }

    test [beta]_cons == 1

    *--------------------------------------------------------*
    *    Export all heterogenous results to one TSV table    *
    *--------------------------------------------------------*
    
estout m1hetero m2hetero m3hetero m4hetero using "$estimations/Allmodels_Heterogenous.tsv", ///
    replace starlevels(* 0.10 ** 0.05 *** 0.01) varlabels("$varlabels") ///
    cells( (b(star label("Estimate") fmt(3)) se(label("Std error") fmt(3)) )) ///
    stats(N ll, fmt(%5.0f %10.3f) labels(N "log-likelihood")) nobaselevels ///
    prehead("Table" "Discounting Function ML Estimates" @title) ///
    title("Concave Utility, Heterogenous Preferences") ///
    legend eqlabels("CRRA function parameter (r)" "PWF parameter (phi)" "PWF parameter (eta)" "Discounting parameter (delta)" "Risk error (mu)" "Time error (nu)" "Discounting parameter (beta)") mlabels(,titles) ///
    postfoot("Results account for clustering at the individual level" "Standard errors in parentheses")


*******************************************************************************
*** 	7.3 -- Margins for Discounting Models                               ***
*******************************************************************************

    *-------------------------------------------*
    *       Exponential Discounting             *
    *-------------------------------------------*
    
    * Delta Equation
    estimates restore m1hetero
    margins, over(wave) predict(equation(delta)) post ///
        saving($estimations/Exponential_Delta, replace)

    * Test for wave effects
    foreach i in 1 2 3 4 5 6 {
        foreach j in `ferest()' {
        test `i'.wave == `j'.wave
            if r(p) < 0.1 {
                di as error r(p) 
            }
        }
    }

    *-------------------------------------------*
    *       Quasi-Hyperbolic Discounting        *
    *-------------------------------------------*
    
    * Beta Equation 
	estimates restore m2hetero
	margins, over(wave) predict(equation(beta)) post ///
        saving($estimations/Quasi_Beta, replace)

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
    estimates restore m2hetero
    margins, over(wave) predict(equation(delta)) post /// 
        saving($estimations/Quasi_Delta, replace)

    * Test for wave effects
    foreach i in 1 2 3 4 5 6 {
        foreach j in `ferest()' {
        test `i'.wave == `j'.wave
            if r(p) < 0.1 {
                di as error r(p) 
            }
        }
    }

    *-------------------------------------------*
    *       Hyperbolic Discounting             *
    *-------------------------------------------*
    
    * Delta Equation
    estimates restore m3hetero
    margins, over(wave) predict(equation(delta)) post ///
        saving($estimations/Hyperbolic_Delta, replace)

    * Test for wave effects
    foreach i in 1 2 3 4 5 6 {
        foreach j in `ferest()' {
        test `i'.wave == `j'.wave
            if r(p) < 0.1 {
                di as error r(p) 
            }
        }
    }

    *-------------------------------------------*
    *       Weibull Discounting                 *
    *-------------------------------------------*
    
    * Beta Equation 
	estimates restore m4hetero
	margins, over(wave) predict(equation(beta)) post ///
        saving($estimations/Weibull_Beta, replace)

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
    estimates restore m4hetero
    margins, over(wave) predict(equation(delta)) post /// 
        saving($estimations/Weibull_Delta, replace)

    * Test for wave effects
    foreach i in 1 2 3 4 5 6 {
        foreach j in `ferest()' {
        test `i'.wave == `j'.wave
            if r(p) < 0.1 {
                di as error r(p) 
            }
        }
    }


*******************************************************************************
*** 	7.4 -- Getting metrics for the graphs						        ***
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

tab lldelay if time==1 & ssdelay==0
tab lldelay if time==1 & ssdelay==7
tab lldelay if time==1 & ssdelay==14 
gen delaydiff = lldelay - ssdelay 
tab delaydiff, m 


*******************************************************************************

di as error "End of Analysis do-file" 
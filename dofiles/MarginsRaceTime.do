*********************************************************************
*   DO FILE: Margins for Ethnic Groups and Tests of Significance    *
*   Estimates margins on time preferences 							*
*********************************************************************

* Start log file 
cap log close 
log using "$logfiles/Log_Time_3_MarginsRace.txt", text replace

*---------------------*
*  Exponential        *
*---------------------*

	* R500 in 14 days 
	estimates restore m1hetero
	margins, over(race wave) ///
		expression(500*(1/((1+predict(equation(delta)))^(14/365)))) post ///
		saving($explanatory/Race_Exponential, replace)

	* Test whether the margins for race differs between two waves under Exponential 
	* Race = 0 (Black) [Exp]
	forvalues w1 = 1/5 {
		local w2start = `w1' + 1
		forvalues w2 = `w2start'/6 {
			test 0.race#`w1'.wave == 0.race#`w2'.wave
			}
		}
	* Race = 1 (Asian/Indian) [Exp]
	forvalues w1 = 1/5 {
		local w2start = `w1' + 1
		forvalues w2 = `w2start'/6 {
			test 1.race#`w1'.wave == 1.race#`w2'.wave
			}
		}
	* Race = 2 (Coloured) [Exp]
	forvalues w1 = 1/5 {
		local w2start = `w1' + 1
		forvalues w2 = `w2start'/6 {
			test 2.race#`w1'.wave == 2.race#`w2'.wave
			}
		}
	* Race = 3 (White) [Exp]
	forvalues w1 = 1/5 {
		local w2start = `w1' + 1
		forvalues w2 = `w2start'/6 {
			test 3.race#`w1'.wave == 3.race#`w2'.wave
			}
		}

	* Test whether the margins differ between race groups within each wave under Exponential
	forvalues w = 1/6 {
		* Joint test: all four race groups equal in wave `w' [Exp]
		di as text _n "Wave `w': joint test of equality across race groups [Exp]"
		test (0.race#`w'.wave == 1.race#`w'.wave) ///
			 (0.race#`w'.wave == 2.race#`w'.wave) ///
			 (0.race#`w'.wave == 3.race#`w'.wave)

		* Pairwise tests between race groups in wave `w' [Exp]
		forvalues r1 = 0/2 {
			local r2start = `r1' + 1
			forvalues r2 = `r2start'/3 {
				test `r1'.race#`w'.wave == `r2'.race#`w'.wave
				}
			}
		}

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

	* Test whether the margins for race differs between two waves under QH 
	* Race = 0 (Black) [QH]
	forvalues w1 = 1/5 {
		local w2start = `w1' + 1
		forvalues w2 = `w2start'/6 {
			test 0.race#`w1'.wave == 0.race#`w2'.wave
			}
		}
	* Race = 1 (Asian/Indian) [QH]
	forvalues w1 = 1/5 {
		local w2start = `w1' + 1
		forvalues w2 = `w2start'/6 {
			test 1.race#`w1'.wave == 1.race#`w2'.wave
			}
		}
	* Race = 2 (Coloured) [QH]
	forvalues w1 = 1/5 {
		local w2start = `w1' + 1
		forvalues w2 = `w2start'/6 {
			test 2.race#`w1'.wave == 2.race#`w2'.wave
			}
		}
	* Race = 3 (White) [QH]
	forvalues w1 = 1/5 {
		local w2start = `w1' + 1
		forvalues w2 = `w2start'/6 {
			test 3.race#`w1'.wave == 3.race#`w2'.wave
			}
		}

	* Test whether the margins differ between race groups within each wave under Exponential
	forvalues w = 1/6 {
		* Joint test: all four race groups equal in wave `w' [Exp]
		di as text _n "Wave `w': joint test of equality across race groups [Exp]"
		test (0.race#`w'.wave == 1.race#`w'.wave) ///
			 (0.race#`w'.wave == 2.race#`w'.wave) ///
			 (0.race#`w'.wave == 3.race#`w'.wave)

		* Pairwise tests between race groups in wave `w' [Exp]
		forvalues r1 = 0/2 {
			local r2start = `r1' + 1
			forvalues r2 = `r2start'/3 {
				test `r1'.race#`w'.wave == `r2'.race#`w'.wave
				}
			}
		}

*---------------------*
*  Hyperbolic 	      *
*---------------------*

	* R500 in 14 days 
	estimates restore m3hetero
	asdoc margins, over(race wave) expression(500*(1/(1+predict(equation(delta))*(14/365)))) ///
		append save($stata_tables/Racial_Groups) label dec(0) ///
		title(Hyperbolic - PV for R500 and 14 days) ///
		saving($explanatory/Race_Hyperbolic, replace) post

	* Test whether the margins for race differs between two waves under Hyperbolic 
	* Race = 0 (Black) [H]
	forvalues w1 = 1/5 {
		local w2start = `w1' + 1
		forvalues w2 = `w2start'/6 {
			test 0.race#`w1'.wave == 0.race#`w2'.wave
			}
		}
	* Race = 1 (Asian/Indian) [H]
	forvalues w1 = 1/5 {
		local w2start = `w1' + 1
		forvalues w2 = `w2start'/6 {
			test 1.race#`w1'.wave == 1.race#`w2'.wave
			}
		}
	* Race = 2 (Coloured) [H]
	forvalues w1 = 1/5 {
		local w2start = `w1' + 1
		forvalues w2 = `w2start'/6 {
			test 2.race#`w1'.wave == 2.race#`w2'.wave
			}
		}
	* Race = 3 (White) [H]
	forvalues w1 = 1/5 {
		local w2start = `w1' + 1
		forvalues w2 = `w2start'/6 {
			test 3.race#`w1'.wave == 3.race#`w2'.wave
			}
		}

	* Test whether the margins differ between race groups within each wave under Exponential
	forvalues w = 1/6 {
		* Joint test: all four race groups equal in wave `w' [Exp]
		di as text _n "Wave `w': joint test of equality across race groups [Exp]"
		test (0.race#`w'.wave == 1.race#`w'.wave) ///
			 (0.race#`w'.wave == 2.race#`w'.wave) ///
			 (0.race#`w'.wave == 3.race#`w'.wave)

		* Pairwise tests between race groups in wave `w' [Exp]
		forvalues r1 = 0/2 {
			local r2start = `r1' + 1
			forvalues r2 = `r2start'/3 {
				test `r1'.race#`w'.wave == `r2'.race#`w'.wave
				}
			}
		}

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

	* Test whether the margins for race differs between two waves under Weibull 
	* Race = 0 (Black) [W]
	forvalues w1 = 1/5 {
		local w2start = `w1' + 1
		forvalues w2 = `w2start'/6 {
			test 0.race#`w1'.wave == 0.race#`w2'.wave
			}
		}
	* Race = 1 (Asian/Indian) [W]
	forvalues w1 = 1/5 {
		local w2start = `w1' + 1
		forvalues w2 = `w2start'/6 {
			test 1.race#`w1'.wave == 1.race#`w2'.wave
			}
		}
	* Race = 2 (Coloured) [W]
	forvalues w1 = 1/5 {
		local w2start = `w1' + 1
		forvalues w2 = `w2start'/6 {
			test 2.race#`w1'.wave == 2.race#`w2'.wave
			}
		}
	* Race = 3 (White) [W]
	forvalues w1 = 1/5 {
		local w2start = `w1' + 1
		forvalues w2 = `w2start'/6 {
			test 3.race#`w1'.wave == 3.race#`w2'.wave
			}
		}

	* Test whether the margins differ between race groups within each wave under Exponential
	forvalues w = 1/6 {
		* Joint test: all four race groups equal in wave `w' [Exp]
		di as text _n "Wave `w': joint test of equality across race groups [Exp]"
		test (0.race#`w'.wave == 1.race#`w'.wave) ///
			 (0.race#`w'.wave == 2.race#`w'.wave) ///
			 (0.race#`w'.wave == 3.race#`w'.wave)

		* Pairwise tests between race groups in wave `w' [Exp]
		forvalues r1 = 0/2 {
			local r2start = `r1' + 1
			forvalues r2 = `r2start'/3 {
				test `r1'.race#`w'.wave == `r2'.race#`w'.wave
				}
			}
		}


*********************************************************************

log close 
di as error "End of MarginsRace do-file" 
*****************************************************************************
* This is the Cleaning do file to clean and generate relevant variables     *
* for the Time Preference Task in South Africa. This fits in section 6      *
* of the Main do file titled Main.do.                                       *
*                                                                           *
* Date first generated:             28 August 2025                          *
* Created by:                       Rinelle Chetty                          * 
*****************************************************************************


*******************************************************************************
***		6.1 -- Drop USA & unnecessary variables, drop other tasks           ***
*******************************************************************************

* Keep only Risk and Time task observations 
drop if task == 3 | task == 4 // drop CA data and beliefs data 

* Drop some USA- and SA-specific variables and other unnecessary variables
drop expenditure_daily expenditure_daily_cat
drop prize4* prob4*
drop residence_type rsa_home_type 
drop rsa_householdN
drop rsa_covid_info_* rsa_covid_information rsa_information_confidence 
drop rsa_experts_confidence rsa_covid_one_month rsa_covid_peak_capetown rsa_covid_situation_one_month 
drop covid_*     
drop hyp_risk BIS2Attimpulse BIS2Motorimpulse 
drop householdN-BIS2Nonplanimpulse
drop race location_completed relationship financial_situation 
drop rsa_covid_personal_* 

* Drop depression and anxiety line items 
drop depression_pleasure depression_hopeless depression_sleep depression_energy 
drop depression_appetite depression_failure depression_concentration 
drop depression_movement depression_dead
drop anxiety_nervous anxiety_worry_uncontrolled anxiety_worry_multiple 
drop anxiety_relax anxiety_restless anxiety_irritable anxiety_afraid

* Keep the variables for SA analyses
/*keep subjectid wave taskorder showupfee period age male rsa_* ///
work anxiety_total anxiety_difficulty depression_total depression_difficulty ///
task qid choice risk time prize* prob* uMax uMin */ 

* Drop beliefs and CA task variables 
drop beliefs gid alpha beta text bin* token_allocation_bin* 
drop ca cashorterdays calongerdays camoneysmaller camoneylarger calsprob

* Drop all USA data
drop usa_*
drop if usa == 1
drop usa


*******************************************************************************
***		6.2 -- Generate, relabel, and rename variables 		                ***
*******************************************************************************

rename rsa_race race 
rename subjectid id 


*******************************************************************************
***		6.3 -- Calculate Covid Deaths                		                ***
*******************************************************************************

* Generate Covid cumulative deaths scale for complementary analyses 
capture: generate covid_scale_deaths = 0
replace covid_scale_deaths = 106031 if wave == 1
replace covid_scale_deaths = 127787 if wave == 2
replace covid_scale_deaths = 154093 if wave == 3
replace covid_scale_deaths = 183684 if wave == 4
replace covid_scale_deaths = 206257 if wave == 5
replace covid_scale_deaths = 229214 if wave == 6

* Normalize
summ covid_scale_deaths
local max = r(max)
replace covid_scale_deaths = covid_scale_deaths/`max'
summ covid_scale_deaths

* Add the squared term
capture: generate covid_scale_deaths_sq = covid_scale_deaths * covid_scale_deaths


*******************************************************************************
***		6.4 -- New Condensed Anxiety and Depression Variables               ***
*******************************************************************************

* Anxiety condensed (from 21 items to 4 items, according to GAD-7)
tab anxiety_total, m 
gen anxcat = . 
replace anxcat = 1 if (anxiety_total >= 0  & anxiety_total <= 4)    // minimal anxiety
replace anxcat = 2 if (anxiety_total >= 5  & anxiety_total <= 9)    // mild 
replace anxcat = 3 if (anxiety_total >= 10 & anxiety_total <= 14)   // moderate
replace anxcat = 4 if (anxiety_total >= 15)                         // severe anxiety 
* 
label define anxietylabel 1 "Minimal" 2 "Mild" 3 "Moderate" 4 "Severe"
label val anxcat anxietylabel 
numlabel, add 
tab anxcat, m  


* Depression condensed (from 27 categories to 5 categories)
tab depression_total, m  
gen depcat = . 
replace depcat = 1 if (depression_total >= 0  & depression_total <= 4)  // minimal
replace depcat = 2 if (depression_total >= 5  & depression_total <= 9)  // mild 
replace depcat = 3 if (depression_total >= 10 & depression_total <= 14) // moderate
replace depcat = 4 if (depression_total >= 15 & depression_total <= 19) // moderate severe
replace depcat = 5 if (depression_total >= 20 & depression_total <= 27) // severe 
tab depression_total depcat, m 
* 
label define depressionlabel 1 "Minimal Depression" 2 "Mild Depression" 3 "Moderate Depression" 4 "Moderate Severe Depression" 5 "Severe Depression"
label val depcat depressionlabel 
numlabel, add 
tab depcat 


* Depression condensed (from 27 categories to 3 categories)
tab depcat, m  
gen depcat2 = . 
replace depcat2 = 1 if (depcat == 1)                // no or minimal depression 
replace depcat2 = 2 if (depcat == 2 | depcat == 3)  // mild or moderate depression 
replace depcat2 = 3 if (depcat == 4 | depcat == 5)  // moderate severe or severe depression 
tab depcat depcat2, m
* 
label define depressionlabel2 1 "None or Minimal Depression" 2 "Mild or Moderate Depression" 3 "Moderate Severe or Severe Depression"
label val depcat2 depressionlabel2
numlabel, add 
tab depcat2 


*******************************************************************************

di as error "End of Cleaning do-file" 
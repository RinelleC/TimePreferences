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

di as error "End of Cleaning do-file" 
*****************************************************************************
* This is the main do file to analyse Beliefs in South Africa               *
* It relies on two MLfunctions files and three do-files.                    *
* The main folder for this analysis is titled "beliefs".					*
*                                                                           *
* Date first generated:             19 July 2023                            *
* Created by:                       Rinelle Chetty                          * 
*****************************************************************************

clear all 

*******************************************************************************
*** 	1. Set up global file paths and stata version 						***
*******************************************************************************

* Move to main folder - remove this before sharing 
global mainfolder "/Volumes/RinellePhD/SAbeliefs" 
cd $mainfolder

* Set paths of subfolders
global dofiles 		"dofiles"
global logfiles 	"logfiles"
global stata_tables	"stata_tables"
global graphs 		"stata_graphs" 
global figures 		"figures"
	global simple 		"$figures/simplefigures"
	global figuresSA 	"$figures/figuresSA"
	global figuresUSA 	"$figures/figuresUSA"

* Start log file 
cap log close 
//log using "$logfiles/Log_Beliefs_Analysis.txt", replace 
log using "$logfiles/Log_Beliefs_Analysis.smcl", replace 

* Drop any labels in memory
capture: label drop _all

* Specify version
capture: version 15.1
capture: version 16.1
capture: version 17
set more off

* Document what ran
about
capture: timer clear
timer on 1

*******************************************************************************
*** 	2. Set up globals for running separate do files or analyses later 	***
*******************************************************************************

* Global for installing packages
global doPACKAGES	"n"

* Global calculate covid deaths scale 
global doDEATHS		"y"

* Globals for cleaning & relabelling data
global doCLEAN		"y"

* Globals for estimations
global doBELIEFS	"y"

* Globals for moments (mean & standard deviations) estimations 
global doMOMENTS	"y"

* Global for figures
global doFIGURES	"n" 

* Global for appendix
//global doAPPENDIX	"n"


*******************************************************************************
*** 	3. Configurations, fonts and graph colours							***
*******************************************************************************

* Allow big machines to process (needs capture so SE can run; comment out if you want to avoid too many processors)
capture: set processors 2
capture: set processors 4
capture: set processors 32

* Need to ensure replications across computers with multiple processors
set seed 1234

* Graphics font
graph set window fontface "Garamond"

* Graphics scheme
set scheme s1color


*******************************************************************************
*** 	4. Install all packages						 						***
*******************************************************************************

if "$doPACKAGES" == "y"		{
	cap net install asdoc.pkg
	cap net install st0085_2.pkg
	cap net install gr0092.pkg
	cap net install grc1leg2, from("http://digital.cgdev.org/doc/stata/MO/Misc")
	cap ssc install dataex
	cap net install gr0075.pkg
	cap ssc install estout
	cap ssc install outreg2
	cap net install tabout.pkg
}


*******************************************************************************
*** 	5. Run code to define programs 										***
*******************************************************************************

* Get the ML functions for Risk and Beliefs installed
qui {
	do MLfunctionsRisk
	do MLfunctionsBeliefs
}


*******************************************************************************
*** 	6. Open the data, drop other tasks, calculate deaths scale 			***
*******************************************************************************

* Unzip the data file
unzipfile ExpData.zip, replace

* Open original data file 
use ExpData.dta                                                                             

* Keep only Risk and Belief task observations 
drop if task == 2 | task == 3 // drop time and CA data

* Calculate covid deaths 
if "$doDEATHS" == "y" {
	do "dofiles/CovidDeathsScale"
}


*******************************************************************************
*** 	7. Clean and relabel the data, drop USA data 				 		***
*******************************************************************************

* Clean and prepare the data 
do "dofiles/CleanBeliefs"

* Save final dataset 
save "beliefsdata.dta", replace 

* Delete un-zipped version of original data 
erase "ExpData.dta"

//use "beliefsdata.dta" 


*******************************************************************************
*** 	8. Conduct Beliefs Analyses and Generate Figures 					***
*******************************************************************************

if "$doBELIEFS" == "y" {
	do "dofiles/AnalysisBeliefs"
}

if "$doMOMENTS" == "y" {
	do "dofiles/MomentsBeliefs"
}

if "$doFIGURES" == "y" {
	do "dofiles/FiguresBeliefs" 
}


*********************************************************************************

* End Timer
timer off 1
timer list
local secs = r(t1)
local mins = `secs'/60
local hrs = `mins'/60
local secs_ = string(`secs', "%10.0f")
local mins_ = string(`mins', "%4.1f")
local hrs_ = string(`hrs', "%4.2f")
di "Calculations took `secs_' seconds, `mins_' minutes, or `hrs_' hours."

cap log close 

* END of MAIN.do
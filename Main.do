*****************************************************************************
* This is the main do file to analyse Time Preferences in South Africa      *
* It relies on two MLfunctions files and three do-files.                    *
*                                                                           *
* Date first generated:             28 August 2025                          *
* Created by:                       Rinelle Chetty                          * 
*****************************************************************************

clear all 

*******************************************************************************
*** 	1. Set up global file paths and stata version 						***
*******************************************************************************

* Move to main folder - remove this before sharing 
global mainfolder "/Volumes/RinellePhD/TimePreferences" 
cd $mainfolder

* Set paths of subfolders
global dofiles 		"dofiles"
global estimations	"estimates"
global stata_tables	"stata_tables"
global graphs 		"stata_graphs" 
global figures 		"figures"
	global simple 		"$figures/simplefigures"

* Start log file 
cap log close 
log using "Log_Time_Analysis.txt", text replace 	// text file

* Drop any labels in memory
capture: label drop _all

* Specify version
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
global doPACKAGES		"n"

* Globals for cleaning & relabelling data
global doCLEAN			"y"

* Globals for estimations
global doTIMEANALYSIS	"y"

* Globals for margins commands and present value table
global doMARGINS		"n"

* Global for figures
global doFIGURES		"y" 


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

* Get the ML functions for Risk and Time installed
qui {
	do MLfunctionsRisk
	do MLfunctionsTime
	}


*******************************************************************************
*** 	6. Open the data, calculate deaths, clean and relabel, drop vars	***
*******************************************************************************

* Unzip the data file
unzipfile ExpData.zip, replace

* Open original data file 
use ExpData.dta                                                                             

* Clean and prepare the data 
do "dofiles/CleanTime"

* Save final dataset 
save "timedata.dta", replace 

* Delete un-zipped version of original data 
erase "ExpData.dta"


*******************************************************************************
*** 	7. Conduct Time Analyses and Generate Figures 						***
*******************************************************************************

if "$doTIMEANALYSIS" == "y" {
	do "dofiles/AnalysisTime"
}

if "$doMARGINS" == "y" {
	do "dofiles/MarginsTime"
}

if "$doFIGURES" == "y" {
	do "dofiles/FiguresTime" 
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


*******************************************************************************
* END of MAIN.do
*******************************************************************************

di as error "END of TIME PREFERENCES" 
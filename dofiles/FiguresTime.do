*********************************************************************
*   DO FILE: Figures                                                *
*   Generates all the figures on covid deaths, and present values   *
*********************************************************************

* Start log file 
cap log close 
log using "$logfiles/Log_Time_5_Figures.txt", text replace

*********************************************************************
*********                  JHU SA Covid Data                *********
*********************************************************************

* Open JHU Data and set file paths 
use "$covid/jhu_data_rsa.dta", clear   
sort date                                       // format: mm/dd/yyyy 

* Daily infections and deaths 
generate confirmed_sa_daily 	= confirmed_sa[_n] - confirmed_sa[_n-1]	
generate deaths_sa_daily 	    = deaths_sa[_n] - deaths_sa[_n-1]

* Generate smoothed data 
lowess confirmed_sa_daily date, bwidth(0.2) generate(confirmed_sa_daily_s)  nograph
lowess deaths_sa_daily date,    bwidth(0.2) generate(deaths_sa_daily_s)     nograph

* Create deciles
xtile dec_c_sa = confirmed_sa_daily_s, nq(10)
xtile dec_d_sa = deaths_sa_daily_s, nq(10)

* Add labels 
mylabels 0(3000)12000,  myscale(@) format(%7.0fc) local(confirmed_sa)   // cases 
mylabels 0(70)280,      myscale(@) format(%2.0fc) local(deaths_sa)      // deaths

* See what the dates are to get better xlabel values
foreach m in jan feb mar apr may jun jul aug sep oct nov dec {
	di d(1`m'2020') " " _c
}
* Full display is: 21915 21946 21975 22006 22036 22067 22097 22128 22159 22189 22220 22250
* Select the months from February onwards
local months = "21946 21975 22006 22036 22067 22097 22128 22159 22189 22220 22250"

* Max values of infections and deaths 
su confirmed_sa_daily_s
di r(max)
su deaths_sa_daily_s 
di r(max) 

* Generate a more dense bar
expand 1000

* Generate day and month variables 
generate int day = day(date)
generate int month = month(date)

* Retain months we have experiments for and regenerate
drop if month<5

* Generate the bar legends
generate s = uniform()*12000
twoway  (bar s date if dec_c_sa == 1,   sort fcolor(blue*0.05) lcolor(blue*0.05)) || ///
        (bar s date if dec_c_sa == 2,   sort fcolor(blue*0.15) lcolor(blue*0.15)) || ///
        (bar s date if dec_c_sa == 3,   sort fcolor(blue*0.25) lcolor(blue*0.25)) || ///
        (bar s date if dec_c_sa == 4,   sort fcolor(blue*0.35) lcolor(blue*0.35)) || ///
        (bar s date if dec_c_sa == 5,   sort fcolor(blue*0.45) lcolor(blue*0.45)) || ///
        (bar s date if dec_c_sa == 6,   sort fcolor(blue*0.55) lcolor(blue*0.55)) || ///
        (bar s date if dec_c_sa == 7,   sort fcolor(blue*0.65) lcolor(blue*0.65)) || ///
        (bar s date if dec_c_sa == 8,   sort fcolor(blue*0.75) lcolor(blue*0.75)) || ///
        (bar s date if dec_c_sa == 9,   sort fcolor(blue*0.85) lcolor(blue*0.85)) || ///
        (bar s date if dec_c_sa == 10,  sort fcolor(blue*0.95) lcolor(blue*0.95)), ///
            legend(off) ytitle("") ///
            plotregion(lcolor(black) lwidth(thin)) ///
            ylabel(none, labcolor(white) angle(horizontal) tlcolor(white)) ///
            xtitle("") xlabel(none, nolabels noticks) fysize(7.5) ///
            saving($covid/c_sa_bar, replace)
 
replace s = uniform()*250
twoway  (bar s date if dec_d_sa == 1,   sort fcolor(red*0.05) lcolor(red*0.05)) || ///
        (bar s date if dec_d_sa == 2,   sort fcolor(red*0.15) lcolor(red*0.15)) || ///
        (bar s date if dec_d_sa == 3,   sort fcolor(red*0.25) lcolor(red*0.25)) || ///
        (bar s date if dec_d_sa == 4,   sort fcolor(red*0.35) lcolor(red*0.35)) || ///
        (bar s date if dec_d_sa == 5,   sort fcolor(red*0.45) lcolor(red*0.45)) || ///
        (bar s date if dec_d_sa == 6,   sort fcolor(red*0.55) lcolor(red*0.55)) || ///
        (bar s date if dec_d_sa == 7,   sort fcolor(red*0.65) lcolor(red*0.65)) || ///
        (bar s date if dec_d_sa == 8,   sort fcolor(red*0.75) lcolor(red*0.75)) || ///
        (bar s date if dec_d_sa == 9,   sort fcolor(red*0.85) lcolor(red*0.85)) || ///
        (bar s date if dec_d_sa == 10,  sort fcolor(red*0.95) lcolor(red*0.95)), ///
            legend(off) ytitle("")  ///
            plotregion(lcolor(black) lwidth(thin)) ///
            ylabel(none, labcolor(white) angle(horizontal) tlcolor(white)) ///
            xtitle("") xlabel(none, nolabels noticks) fysize(7.5) ///
            saving($covid/d_sa_bar, replace)

* Save data for regenerating the bar
keep s date dec_c_sa dec_d_sa

*********************************************************************
********        South African Time - Present Values          ********
*********************************************************************

* Set size of LL reward
local LL "500"

* Set ylabel scaling for all subsequent graphs
local ylabel "400(10)500"

* Reset
mylabels 400(10)500, myscale(@) prefix(R) format(%4.2f) local(ylabel)

* Combine all models starting with exponential
use "$estimations/pv_E_500_14days", clear
* add quasi-hyperbolic
append using "$estimations/pv_QH_500_14days", generate(_by2)
* add hyperbolic
append using "$estimations/pv_H_500_14days", generate(_by3)
* add weibull
append using "$estimations/pv_W_500_14days", generate(_by4)

* set up wave timeline axis 
rename _by1 by1
generate _by1 = date("5/29/2020", "MDY")
format _by1 %td
replace _by1 = date("6/30/2020", "MDY")     if by1 == 2
replace _by1 = date("7/31/2020", "MDY")     if by1 == 3
replace _by1 = date("8/31/2020", "MDY")     if by1 == 4
replace _by1 = date("9/29/2020", "MDY")     if by1 == 5
replace _by1 = date("10/29/2020", "MDY")    if by1 == 6
drop by1
order _by2 _by3 _by4, after(_by1)
sort _by1
save "$estimations/pv500margin", replace

* Set graph colours
local exp_colour    "navy*.7"
local qh_colour     "purple*.7"
local hyp_colour    "forest_green*.7"
local wei_colour    "orange*.7"

* Now plot the combined margins dataset
marginsplot using "$estimations/pv500margin", l1title("Rand", orientation(horizontal)) ///
    ytitle("") ylabel(, angle(horizontal)) title("") ///
    xlabel("", format(%tdm)) xtitle("") ///
    plot1opts(lwidth(thick) lpattern(solid) lcolor(`exp_colour') mcolor(`exp_colour')) ci1opts(lcolor(`exp_colour')) ///
    plot2opts(lwidth(thick) lpattern(dash)  lcolor(`qh_colour')  mcolor(`qh_colour'))  ci2opts(lcolor(`qh_colur'))   ///
    plot3opts(lwidth(thick) lpattern(solid) lcolor(`hyp_colour') mcolor(`hyp_colour')) ci3opts(lcolor(`hyp_colur'))  ///
    plot4opts(lwidth(thick) lpattern(dash)  lcolor(`wei_colour') mcolor(`wei_colour')) ci4opts(lcolor(`wei_colur'))  ///
    legend(order(5 "Exponential" 6 "QH" 7 "Hyperbolic" 8 "Weibull") size(small) cols(1) ring(0) pos(5) nobox) ///
    saving("$timepref/presentvalue", replace)

* Caption
local caption ""The circles represent point estimates with 95% confidence intervals. The solid blue line shows time preferences under Exponential discounting, the dashed purple" "line QH discounting, the solid green line Hyperbolic discounting, and the dashed yellow line Weibull discounting. Daily national COVID-19 infection rate (blue)" "and death rate(red) in South Africa are indicated in the horizontal bars.""

* Combine the graphs and export
gr combine "$timepref/presentvalue.gph" "$covid/c_sa_bar.gph" "$covid/d_sa_bar.gph", ///
    cols(1) imargin(zero) xcommon ///
    title("Discounting Behaviour", size(vlarge)) ///
    subtitle("Based on the Present Value of a R500 Reward Received in 14 Days", ///
	size(medium) margin(medsmall)) caption(`caption', size(vsmall))
graph export "$timepref/discountingbehaviour.pdf", replace

*********************************************************************
*****        Racial Groups - Present Values - Exponential       *****
*********************************************************************

* Load the race x wave margins saved in MarginsTables.do
use "$explanatory/Race_Exponential", clear

rename _by1 race
rename _by2 wave
label define racesa 0 "Black/African" 1 "Asian/Indian" 2 "Coloured" 3 "White", replace
label values race racesa

* Truncate the CIs for display only (Asian/Indian wave 1 has a very wide CI)
local ymax = 510
generate double ci_lb_plot = max(_ci_lb, 0)
generate double ci_ub_plot = min(_ci_ub, `ymax')

* Colors for each race group
local c0 "cranberry*.7"
local c1 "midblue*.7"
local c2 "midgreen*.7"
local c3 "dkorange*.7"

* full graph
twoway 	(rcap ci_lb_plot ci_ub_plot wave if race == 0, lcolor(`c0'%50)) ///
		(rcap ci_lb_plot ci_ub_plot wave if race == 1, lcolor(`c1'%50)) ///
		(rcap ci_lb_plot ci_ub_plot wave if race == 2, lcolor(`c2'%50)) ///
		(rcap ci_lb_plot ci_ub_plot wave if race == 3, lcolor(`c3'%50)) ///
		(connected _margin wave if race == 0, lcolor(`c0') mcolor(`c0') lwidth(medthick) msymbol(O)) ///
		(connected _margin wave if race == 1, lcolor(`c1') mcolor(`c1') lwidth(medthick) msymbol(O)) ///
		(connected _margin wave if race == 2, lcolor(`c2') mcolor(`c2') lwidth(medthick) msymbol(O)) ///
		(connected _margin wave if race == 3, lcolor(`c3') mcolor(`c3') lwidth(medthick) msymbol(O)), ///
		xlabel(1 "Wave 1" 2 "Wave 2" 3 "Wave 3" 4 "Wave 4" 5 "Wave 5" 6 "Wave 6", labgap(small)) ///
		xtitle("") ///
        ylabel(, angle(horizontal) labgap(small)) ///
		ytitle("Rand", orientation(horizontal) margin(r=3)) ///
		title("Discounting Behaviour by Ethnic Group" "Exponential Discounting", linegap(2) margin(medium) size(vlarge) color(black)) ///
		legend(order(5 "Black/African" 6 "Asian/Indian" 7 "Coloured" 8 "White") ///
			   rows(1) position(6) region(lcolor(black)) size(small) symxsize(*.6) keygap(*.6) colgap(*2)) ///
		caption("Point estimates represented by the circles with 95% confidence intervals. Estimates show the present values for each ethnic group," "under Exponential Discounting.", size(vsmall)) ///
		graphregion(fcolor(white) color(white)) scheme(s1color) xsize(7) ysize(5) 
graph export "$explanatory/Race_Exponential.png", replace

* plain graph
twoway 	(rcap ci_lb_plot ci_ub_plot wave if race == 0, lcolor(`c0'%50)) ///
		(rcap ci_lb_plot ci_ub_plot wave if race == 1, lcolor(`c1'%50)) ///
		(rcap ci_lb_plot ci_ub_plot wave if race == 2, lcolor(`c2'%50)) ///
		(rcap ci_lb_plot ci_ub_plot wave if race == 3, lcolor(`c3'%50)) ///
		(connected _margin wave if race == 0, lcolor(`c0') mcolor(`c0') lwidth(medthick) msymbol(O)) ///
		(connected _margin wave if race == 1, lcolor(`c1') mcolor(`c1') lwidth(medthick) msymbol(O)) ///
		(connected _margin wave if race == 2, lcolor(`c2') mcolor(`c2') lwidth(medthick) msymbol(O)) ///
		(connected _margin wave if race == 3, lcolor(`c3') mcolor(`c3') lwidth(medthick) msymbol(O)), ///
		xlabel(1 "Wave 1" 2 "Wave 2" 3 "Wave 3" 4 "Wave 4" 5 "Wave 5" 6 "Wave 6", labgap(small)) ///
		xtitle("") ///
        ylabel(, angle(horizontal) labgap(small)) ///
		ytitle("Rand", orientation(horizontal) margin(r=3)) ///
        title("Exponential", box ring(0) pos(1) fcolor(khaki) color(black) size(medsmall)) ///
		legend(order(5 "Black/African" 6 "Asian/Indian" 7 "Coloured" 8 "White") ///
			   rows(1) position(6) region(lcolor(black)) size(small) symxsize(*.6) keygap(*.6) colgap(*2)) ///
		graphregion(fcolor(white) color(white)) scheme(s1color) xsize(7) ysize(5) 
graph save "$explanatory/Race_Exponential_plain", replace 

*********************************************************************
*****        Racial Groups - Present Values - Hyperbolic        *****
*********************************************************************

* Load the race x wave margins saved in MarginsTables.do
use "$explanatory/Race_Hyperbolic", clear

rename _by1 race
rename _by2 wave
label define racesa 0 "Black/African" 1 "Asian/Indian" 2 "Coloured" 3 "White", replace
label values race racesa

* Truncate the CIs for display only (Asian/Indian wave 1 has a very wide CI)
local ymax = 510
generate double ci_lb_plot = max(_ci_lb, 0)
generate double ci_ub_plot = min(_ci_ub, `ymax')

* Colors for each race group
local c0 "cranberry*.7"
local c1 "midblue*.7"
local c2 "midgreen*.7"
local c3 "dkorange*.7"

* full graph
twoway 	(rcap ci_lb_plot ci_ub_plot wave if race == 0, lcolor(`c0'%50)) ///
		(rcap ci_lb_plot ci_ub_plot wave if race == 1, lcolor(`c1'%50)) ///
		(rcap ci_lb_plot ci_ub_plot wave if race == 2, lcolor(`c2'%50)) ///
		(rcap ci_lb_plot ci_ub_plot wave if race == 3, lcolor(`c3'%50)) ///
		(connected _margin wave if race == 0, lcolor(`c0') mcolor(`c0') lwidth(medthick) msymbol(O)) ///
		(connected _margin wave if race == 1, lcolor(`c1') mcolor(`c1') lwidth(medthick) msymbol(O)) ///
		(connected _margin wave if race == 2, lcolor(`c2') mcolor(`c2') lwidth(medthick) msymbol(O)) ///
		(connected _margin wave if race == 3, lcolor(`c3') mcolor(`c3') lwidth(medthick) msymbol(O)), ///
		xlabel(1 "Wave 1" 2 "Wave 2" 3 "Wave 3" 4 "Wave 4" 5 "Wave 5" 6 "Wave 6", labgap(small)) ///
		xtitle("") ///
        ylabel(, angle(horizontal) labgap(small)) ///
		ytitle("Rand", orientation(horizontal) margin(r=3)) ///
		title("Discounting Behaviour by Ethnic Group" "Hyperbolic Discounting", linegap(2) margin(medium) size(vlarge) color(black)) ///
		legend(order(5 "Black/African" 6 "Asian/Indian" 7 "Coloured" 8 "White") ///
			   rows(1) position(6) region(lcolor(black)) size(small) symxsize(*.6) keygap(*.6) colgap(*2)) ///
		caption("Point estimates represented by the circles with 95% confidence intervals. Estimates show the present values for each ethnic group," "under Hyperbolic Discounting.", size(vsmall)) ///
		graphregion(fcolor(white) color(white)) scheme(s1color) xsize(7) ysize(5) 
graph export "$explanatory/Race_Hyperbolic.png", replace

* plain graph
twoway 	(rcap ci_lb_plot ci_ub_plot wave if race == 0, lcolor(`c0'%50)) ///
		(rcap ci_lb_plot ci_ub_plot wave if race == 1, lcolor(`c1'%50)) ///
		(rcap ci_lb_plot ci_ub_plot wave if race == 2, lcolor(`c2'%50)) ///
		(rcap ci_lb_plot ci_ub_plot wave if race == 3, lcolor(`c3'%50)) ///
		(connected _margin wave if race == 0, lcolor(`c0') mcolor(`c0') lwidth(medthick) msymbol(O)) ///
		(connected _margin wave if race == 1, lcolor(`c1') mcolor(`c1') lwidth(medthick) msymbol(O)) ///
		(connected _margin wave if race == 2, lcolor(`c2') mcolor(`c2') lwidth(medthick) msymbol(O)) ///
		(connected _margin wave if race == 3, lcolor(`c3') mcolor(`c3') lwidth(medthick) msymbol(O)), ///
		xlabel(1 "Wave 1" 2 "Wave 2" 3 "Wave 3" 4 "Wave 4" 5 "Wave 5" 6 "Wave 6", labgap(small)) ///
		xtitle("") ///
        ylabel(, angle(horizontal) labgap(small)) ///
		ytitle("Rand", orientation(horizontal) margin(r=3)) ///
        title("Hyperbolic", box ring(0) pos(1) fcolor(khaki) color(black) size(medsmall)) ///
		legend(order(5 "Black/African" 6 "Asian/Indian" 7 "Coloured" 8 "White") ///
			   rows(1) position(6) region(lcolor(black)) size(small) symxsize(*.6) keygap(*.6) colgap(*2)) ///
		graphregion(fcolor(white) color(white)) scheme(s1color) xsize(7) ysize(5) 
graph save "$explanatory/Race_Hyperbolic_plain", replace 

*********************************************************************
*****        Racial Groups - Present Values - QH                *****
*********************************************************************

* Load the race x wave margins saved in MarginsTables.do
use "$explanatory/Race_QuasiHyperbolic", clear

rename _by1 race
rename _by2 wave
label define racesa 0 "Black/African" 1 "Asian/Indian" 2 "Coloured" 3 "White", replace
label values race racesa

* Truncate the CIs for display only (Asian/Indian wave 1 has a very wide CI)
local ymax = 510
generate double ci_lb_plot = max(_ci_lb, 0)
generate double ci_ub_plot = min(_ci_ub, `ymax')

* Colors for each race group
local c0 "cranberry*.7"
local c1 "midblue*.7"
local c2 "midgreen*.7"
local c3 "dkorange*.7"

* full graph
twoway 	(rcap ci_lb_plot ci_ub_plot wave if race == 0, lcolor(`c0'%50)) ///
		(rcap ci_lb_plot ci_ub_plot wave if race == 1, lcolor(`c1'%50)) ///
		(rcap ci_lb_plot ci_ub_plot wave if race == 2, lcolor(`c2'%50)) ///
		(rcap ci_lb_plot ci_ub_plot wave if race == 3, lcolor(`c3'%50)) ///
		(connected _margin wave if race == 0, lcolor(`c0') mcolor(`c0') lwidth(medthick) msymbol(O)) ///
		(connected _margin wave if race == 1, lcolor(`c1') mcolor(`c1') lwidth(medthick) msymbol(O)) ///
		(connected _margin wave if race == 2, lcolor(`c2') mcolor(`c2') lwidth(medthick) msymbol(O)) ///
		(connected _margin wave if race == 3, lcolor(`c3') mcolor(`c3') lwidth(medthick) msymbol(O)), ///
		xlabel(1 "Wave 1" 2 "Wave 2" 3 "Wave 3" 4 "Wave 4" 5 "Wave 5" 6 "Wave 6", labgap(small)) ///
		xtitle("") ///
        ylabel(, angle(horizontal) labgap(small)) ///
		ytitle("Rand", orientation(horizontal) margin(r=3)) ///
		title("Discounting Behaviour by Ethnic Group" "QH Discounting", linegap(2) margin(medium) size(vlarge) color(black)) ///
		legend(order(5 "Black/African" 6 "Asian/Indian" 7 "Coloured" 8 "White") ///
			   rows(1) position(6) region(lcolor(black)) size(small) symxsize(*.6) keygap(*.6) colgap(*2)) ///
		caption("Point estimates represented by the circles with 95% confidence intervals. Estimates show the present values for each ethnic group," "under QH Discounting.", size(vsmall)) ///
		graphregion(fcolor(white) color(white)) scheme(s1color) xsize(7) ysize(5) 
graph export "$explanatory/Race_QuasiHyperbolic.png", replace

* plain graph
twoway 	(rcap ci_lb_plot ci_ub_plot wave if race == 0, lcolor(`c0'%50)) ///
		(rcap ci_lb_plot ci_ub_plot wave if race == 1, lcolor(`c1'%50)) ///
		(rcap ci_lb_plot ci_ub_plot wave if race == 2, lcolor(`c2'%50)) ///
		(rcap ci_lb_plot ci_ub_plot wave if race == 3, lcolor(`c3'%50)) ///
		(connected _margin wave if race == 0, lcolor(`c0') mcolor(`c0') lwidth(medthick) msymbol(O)) ///
		(connected _margin wave if race == 1, lcolor(`c1') mcolor(`c1') lwidth(medthick) msymbol(O)) ///
		(connected _margin wave if race == 2, lcolor(`c2') mcolor(`c2') lwidth(medthick) msymbol(O)) ///
		(connected _margin wave if race == 3, lcolor(`c3') mcolor(`c3') lwidth(medthick) msymbol(O)), ///
		xlabel(1 "Wave 1" 2 "Wave 2" 3 "Wave 3" 4 "Wave 4" 5 "Wave 5" 6 "Wave 6", labgap(small)) ///
		xtitle("") ///
        ylabel(, angle(horizontal) labgap(small)) ///
		ytitle("Rand", orientation(horizontal) margin(r=3)) ///
        title("Quasi-Hyperbolic", box ring(0) pos(1) fcolor(khaki) color(black) size(medsmall)) ///
		legend(order(5 "Black/African" 6 "Asian/Indian" 7 "Coloured" 8 "White") ///
			   rows(1) position(6) region(lcolor(black)) size(small) symxsize(*.6) keygap(*.6) colgap(*2)) ///
		graphregion(fcolor(white) color(white)) scheme(s1color) xsize(7) ysize(5) 
graph save "$explanatory/Race_QuasiHyperbolic_plain", replace 

*********************************************************************
*****        Racial Groups - Present Values - Weibull           *****
*********************************************************************

* Load the race x wave margins saved in MarginsTables.do
use "$explanatory/Race_Weibull", clear

rename _by1 race
rename _by2 wave
label define racesa 0 "Black/African" 1 "Asian/Indian" 2 "Coloured" 3 "White", replace
label values race racesa

* Truncate the CIs for display only 
local ymax = 510
generate double ci_lb_plot = max(_ci_lb, 0)
generate double ci_ub_plot = min(_ci_ub, `ymax')

* Colors for each race group
local c0 "cranberry*.7"
local c1 "midblue*.7"
local c2 "midgreen*.7"
local c3 "dkorange*.7"

* full graph
twoway 	(rcap ci_lb_plot ci_ub_plot wave if race == 0, lcolor(`c0'%50)) ///
		(rcap ci_lb_plot ci_ub_plot wave if race == 1, lcolor(`c1'%50)) ///
		(rcap ci_lb_plot ci_ub_plot wave if race == 2, lcolor(`c2'%50)) ///
		(rcap ci_lb_plot ci_ub_plot wave if race == 3, lcolor(`c3'%50)) ///
		(connected _margin wave if race == 0, lcolor(`c0') mcolor(`c0') lwidth(medthick) msymbol(O)) ///
		(connected _margin wave if race == 1, lcolor(`c1') mcolor(`c1') lwidth(medthick) msymbol(O)) ///
		(connected _margin wave if race == 2, lcolor(`c2') mcolor(`c2') lwidth(medthick) msymbol(O)) ///
		(connected _margin wave if race == 3, lcolor(`c3') mcolor(`c3') lwidth(medthick) msymbol(O)), ///
		xlabel(1 "Wave 1" 2 "Wave 2" 3 "Wave 3" 4 "Wave 4" 5 "Wave 5" 6 "Wave 6", labgap(small)) ///
		xtitle("") ///
        ylabel(, angle(horizontal) labgap(small)) ///
		ytitle("Rand", orientation(horizontal) margin(r=3)) ///
		title("Discounting Behaviour by Ethnic Group" "Weibull Discounting", linegap(2) margin(medium) size(vlarge) color(black)) ///
		legend(order(5 "Black/African" 6 "Asian/Indian" 7 "Coloured" 8 "White") ///
			   rows(1) position(6) region(lcolor(black)) size(small) symxsize(*.6) keygap(*.6) colgap(*2)) ///
		caption("Point estimates represented by the circles with 95% confidence intervals. Estimates show the present values for each ethnic group," "under Weibull Discounting.", size(vsmall)) ///
		graphregion(fcolor(white) color(white)) scheme(s1color) xsize(7) ysize(5) 
graph export "$explanatory/Race_Weibull.png", replace

* plain graph
twoway 	(rcap ci_lb_plot ci_ub_plot wave if race == 0, lcolor(`c0'%50)) ///
		(rcap ci_lb_plot ci_ub_plot wave if race == 1, lcolor(`c1'%50)) ///
		(rcap ci_lb_plot ci_ub_plot wave if race == 2, lcolor(`c2'%50)) ///
		(rcap ci_lb_plot ci_ub_plot wave if race == 3, lcolor(`c3'%50)) ///
		(connected _margin wave if race == 0, lcolor(`c0') mcolor(`c0') lwidth(medthick) msymbol(O)) ///
		(connected _margin wave if race == 1, lcolor(`c1') mcolor(`c1') lwidth(medthick) msymbol(O)) ///
		(connected _margin wave if race == 2, lcolor(`c2') mcolor(`c2') lwidth(medthick) msymbol(O)) ///
		(connected _margin wave if race == 3, lcolor(`c3') mcolor(`c3') lwidth(medthick) msymbol(O)), ///
		xlabel(1 "Wave 1" 2 "Wave 2" 3 "Wave 3" 4 "Wave 4" 5 "Wave 5" 6 "Wave 6", labgap(small)) ///
		xtitle("") ///
        ylabel(, angle(horizontal) labgap(small)) ///
		ytitle("Rand", orientation(horizontal) margin(r=3)) ///
        title("Weibull", box ring(0) pos(1) fcolor(khaki) color(black) size(medsmall)) ///
		legend(order(5 "Black/African" 6 "Asian/Indian" 7 "Coloured" 8 "White") ///
			   rows(1) position(6) region(lcolor(black)) size(small) symxsize(*.6) keygap(*.6) colgap(*2)) ///
		graphregion(fcolor(white) color(white)) scheme(s1color) xsize(7) ysize(5) 
graph save "$explanatory/Race_Weibull_plain", replace 

*********************************************************************
*****                   Racial Groups - Combining               *****
*********************************************************************

graph use "$explanatory/Race_Exponential_plain.gph",        name(exponential, replace)
graph use "$explanatory/Race_Hyperbolic_plain.gph",         name(hyperbolic, replace)
//graph use "$explanatory/Race_QuasiHyperbolic_plain.gph",    name(quasihyperbolic, replace)
//graph use "$explanatory/Race_Weibull_plain.gph",            name(weib, replace)

* caption 
local caption "The circles represent point estimates. 95% confidence intervals. Estimates show the present values for each ethnic group. "

grc1leg2 	exponential hyperbolic, ///
			cols(2) legendfrom(exponential) ///
			b1title(`b1title', size(small)) ///
			l1title("Rand", size(small) orientation(horizontal)) ///
			position(6) ring(2) graphregion(fcolor(white) color(white)) xtob1title ytol1title ///
			title("Discounting Behaviour By Ethnic Group", size(vlarge) color(black) span justification(right)) ///
			caption(`caption', size(vsmall) margin(small)) 

graph export "$explanatory/RaceCombined.png", replace 

erase     "$explanatory/Race_Exponential_plain.gph"
erase     "$explanatory/Race_Hyperbolic_plain.gph"
erase     "$explanatory/Race_QuasiHyperbolic_plain.gph"
erase     "$explanatory/Race_Weibull_plain.gph"

*********************************************************************

log close  
di as error "End of Figures do-file" 
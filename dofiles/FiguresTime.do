*****************************************************************************
* This is the Analysis do file to analyse Time Preferences in South Africa  *
* and fits within section 8 of the Main do-file.                            * 
*                                                                           *
* Date first generated:             1 September 2025                        *
* Created by:                       Rinelle Chetty                          * 
*****************************************************************************


*********************************************************************************
*******************           JHU SA Covid Data             *********************
*********************************************************************************

clear all 
use "$figures/jhu_data_rsa.dta", clear   
sort date                                               // format: mm/dd/yyyy 

* Daily infections and deaths 
generate confirmed_sa_daily 	= confirmed_sa[_n] - confirmed_sa[_n-1]	
generate deaths_sa_daily 	= deaths_sa[_n] - deaths_sa[_n-1]

* Generate smoothed data 
lowess confirmed_sa_daily date, bwidth(0.2) generate(confirmed_sa_daily_s) nograph
lowess deaths_sa_daily date, bwidth(0.2) generate(deaths_sa_daily_s) nograph

* Create deciles
xtile dec_c_sa = confirmed_sa_daily_s, nq(10)
tab dec_c_sa
xtile dec_d_sa = deaths_sa_daily_s, nq(10)
tab dec_d_sa, m

* Add labels 
mylabels 0(3000)12000,  myscale(@) format(%7.0fc) local(confirmed_sa)   // cases 
mylabels 0(70)280,      myscale(@) format(%2.0fc) local(deaths_sa)      // deaths
//mylabels 0(200)600, myscale(@) format(%1.0fc) local(deaths_sa)

* Identify the six waves 
local wave1 = td(29-05-2020)
local wave2 = td(30-06-2020)
local wave3 = td(31-07-2020)
local wave4 = td(31-08-2020)
local wave5 = td(29-09-2020)
local wave6 = td(29-10-2020)

* See what the dates are to get better xlabel values
foreach m in jan feb mar apr may jun jul aug sep oct nov dec {
	di d(1`m'2020') " " _c
}
* the full display is: 21915 21946 21975 22006 22036 22067 22097 22128 22159 22189 22220 22250
* so just hand-pick the months from february on, not january 2021
local months = "21946 21975 22006 22036 22067 22097 22128 22159 22189 22220 22250"

* Max values of infections and deaths 
su confirmed_sa_daily_s
di r(max)
su deaths_sa_daily_s 
di r(max) 

* Main picture 1 - Infections curve
twoway (line confirmed_sa_daily_s date, lcolor(blue) lwidth(thick)), ytitle("") ///
ytitle(, size(huge) margin(10-pt)) ///
ylabel(`confirmed_sa', angle(horizontal) labgap(small)) ///
xtitle("") xlabel(`months', format(%tdm) labsize(small) labgap(medium)) ///
xline(`wave1' `wave2' `wave3' `wave4' `wave5' `wave6', lcolor(gs14) lwidth(thick)) ///
plotregion(lcolor(black) lwidth(thin)) ///
title("Infections", size(vlarge) margin(small)) saving(c_sa, replace)

* Main picture 2 - Deaths curve
twoway (line deaths_sa_daily_s date, lcolor(red) lwidth(thick)), ytitle("") ///
ytitle(, size(huge) margin(10-pt)) ///
ylabel(`deaths_sa', angle(horizontal) labgap(small)) ///
xtitle("") xlabel(`months', format(%tdm) labsize(small) labgap(medium)) ///
xline(`wave1' `wave2' `wave3' `wave4' `wave5' `wave6', lcolor(gs14) lwidth(thick)) ///
plotregion(lcolor(black) lwidth(thin)) ///
title("Deaths", size(vlarge) margin(small)) saving(d_sa, replace)

* Generate a more dense bar
expand 1000

* Graph of Infections and Deaths - but no bars 
gr combine c_sa.gph d_sa.gph, cols(2) imargin(zero) ///
title("Infections and Deaths", size(vlarge)) ///
subtitle("South Africa") scheme(s1color) saving(figure_A, replace)

* Generate the bar legends and collate
generate s = uniform()*12000
twoway (bar s date if dec_c_sa == 1, sort fcolor(blue*0.05) lcolor(blue*0.05)) || ///
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
ylabel(`confirmed_sa', labcolor(white) angle(horizontal) tlcolor(white)) ///
xtitle("") xlabel(none, nolabels noticks) fysize(7.5) ///
saving(c_sa_bar, replace)
graph export bar_blue.pdf, replace
 
replace s = uniform()*250
twoway (bar s date if dec_d_sa == 1, sort fcolor(red*0.05) lcolor(red*0.05)) || ///
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
ylabel(`deaths_sa', labcolor(white) angle(horizontal) tlcolor(white)) ///
xtitle("") xlabel(none, nolabels noticks) fysize(7.5) ///
saving(d_sa_bar, replace)
graph export bar_red.pdf, replace

* figure 1 - for presentations (white background)
gr combine c_sa.gph d_sa.gph c_sa_bar.gph d_sa_bar.gph, cols(2) imargin(medium) ///
title("Infections and Deaths", size(vlarge)) ///
scheme(s1color) /// 
subtitle("South Africa") saving(figure1, replace) 
gr export figure1.pdf, replace 

* figure 1 - for PhD draft
gr combine c_sa.gph d_sa.gph c_sa_bar.gph d_sa_bar.gph, cols(2) imargin(medium) ///
title("Daily Infections and Deaths in South Africa", size(vlarge)) ///
scheme(s1color) /// 
subtitle("Vertical grey lines show dates of each experimental wave", margin(small)) ///
note("Smoothed data" "Data source: Johns Hopkins Coronavirus Resource Center") ///
saving(figure1_phd, replace) 
gr export figure1_phd.pdf, replace 

* Save data for regenerating the bar
keep s date dec_c_sa dec_d_sa
save legend, replace


*********************************************************************************
**************           South African Time Preferences            **************
*********************************************************************************

* Set size of LL reward
local LL "50"

* Set ylabel scaling for all subsequent graphs
local ylabel "44(1)50"

* Reset
mylabels 44(1)50, myscale(@) prefix(R) format(%4.2f) local(ylabel)

* Get time preference estimates linked to JHU data
use legend, clear
generate int day = day(date)
generate int month = month(date)

* Retain months we have experiments for and regenerate
drop if month<5

* Generate the bar legends
replace s = uniform()*50
twoway  (bar s date if dec_c_us == 1, sort fcolor(blue*0.05) lcolor(blue*0.05)) || ///
        (bar s date if dec_c_us == 2, sort fcolor(blue*0.15) lcolor(blue*0.15)) || ///
        (bar s date if dec_c_us == 3, sort fcolor(blue*0.25) lcolor(blue*0.25)) || ///
        (bar s date if dec_c_us == 4, sort fcolor(blue*0.35) lcolor(blue*0.35)) || ///
        (bar s date if dec_c_us == 5, sort fcolor(blue*0.45) lcolor(blue*0.45)) || ///
        (bar s date if dec_c_us == 6, sort fcolor(blue*0.55) lcolor(blue*0.55)) || ///
        (bar s date if dec_c_us == 7, sort fcolor(blue*0.65) lcolor(blue*0.65)) || ///
        (bar s date if dec_c_us == 8, sort fcolor(blue*0.75) lcolor(blue*0.75)) || ///
        (bar s date if dec_c_us == 9, sort fcolor(blue*0.85) lcolor(blue*0.85)) || ///
        (bar s date if dec_c_us == 10, sort fcolor(blue*0.95) lcolor(blue*0.95)), ///
            legend(off) ytitle("") ///
            l1title("Dollars", color(white) orientation(horizontal)) ///
            ylabel(`ylabel', labcolor(white) angle(horizontal) tlcolor(white)) ///
            xtitle("") xlabel(none, nolabels noticks) fysize(7.5) ///
            saving(c_bar, replace)

twoway  (bar s date if dec_d_us == 1, sort fcolor(red*0.05) lcolor(red*0.05)) || ///
        (bar s date if dec_d_us == 2, sort fcolor(red*0.15) lcolor(red*0.15)) || ///
        (bar s date if dec_d_us == 3, sort fcolor(red*0.25) lcolor(red*0.25)) || ///
        (bar s date if dec_d_us == 4, sort fcolor(red*0.35) lcolor(red*0.35)) || ///
        (bar s date if dec_d_us == 5, sort fcolor(red*0.45) lcolor(red*0.45)) || ///
        (bar s date if dec_d_us == 6, sort fcolor(red*0.55) lcolor(red*0.55)) || ///
        (bar s date if dec_d_us == 7, sort fcolor(red*0.65) lcolor(red*0.65)) || ///
        (bar s date if dec_d_us == 8, sort fcolor(red*0.75) lcolor(red*0.75)) || ///
        (bar s date if dec_d_us == 9, sort fcolor(red*0.85) lcolor(red*0.85)) || ///
        (bar s date if dec_d_us == 10, sort fcolor(red*0.95) lcolor(red*0.95)), ///
            legend(off) ytitle("") ///
            l1title("Dollars", color(white) orientation(horizontal)) ///
            ylabel(`ylabel', labcolor(white) angle(horizontal) tlcolor(white)) ////
            xtitle("") xlabel(none, noticks nolabels) fysize(7.5) ///
            saving(d_bar, replace)


* Get the combined margins.dta files into the same format as the infection and death bars above
use pvExp50margin, clear
append using pvQH50margin, generate(_by2)
rename _by1 by1
generate _by1 = date("5/29/2020", "MDY")
format _by1 %td
replace _by1 = date("6/30/2020", "MDY") if by1 == 2
replace _by1 = date("7/31/2020", "MDY") if by1 == 3
replace _by1 = date("8/31/2020", "MDY") if by1 == 4
replace _by1 = date("9/29/2020", "MDY") if by1 == 5
replace _by1 = date("10/29/2020", "MDY") if by1 == 6
drop by1
order _by2, after(_by1)
sort _by1
save pvExpQH50margin, replace


* Now plot the combined margins dataset
local exp_color "black"
local qh_color "dkorange"
local wide "thick"
marginsplot using pvExpQH50margin, l1title("Dollars", orientation(horizontal)) ///
ytitle("") title("") xlabel("", format(%tdm)) xtitle("") ///
plot1opts(lwidth(`wide') lcolor(`exp_color') mcolor(`exp_color')) ///
ci1opts(lcolor(`exp_color')) ///
plot2opts(lwidth(`wide') lpattern(dash) lcolor(`qh_color') mcolor(`qh_color')) ///
ci2opts(lcolor(`qh_color')) ///
legend(order(3 "Exponential" 4 "Quasi-Hyperbolic") size(medlarge) cols(1) ring(0) pos(2) nobox) ///
yline(50, lcolor(gs14) lwidth(vthick)) ///
ylabel(`ylabel', angle(horizontal)) saving(pv_us, replace)

* caption
local caption ""Point estimates in circles, and 95% confidence intervals above and below in bars" "Solid black line is for Exponential discounting, and dashed orange line is for Quasi-Hyperbolic discounting" "Solid, thick green line is for Exponential pre-pandemic, and dashed, thick green line is for Quasi-Hyperbolic pre-pandemic" "Daily national infection rate (blue) and death rate (red) indicated in bars at bottom: see Figure 1 for legend""

*Combine the graphs
gr combine pv_us.gph c_us_bar.gph d_us_bar.gph, cols(1) imargin(zero) xcommon ///
title("Figure 7: Discounting Behavior", size(vlarge)) ///
subtitle("Present value of $50 reward in 2 weeks" "Green bars show pre-pandemic estimates from the same population", ///
	size(medium) margin(medsmall)) caption(`caption', size(vsmall))  saving(figure7, replace)


*******************************************************************************

di as error "End of Figures do-file" 
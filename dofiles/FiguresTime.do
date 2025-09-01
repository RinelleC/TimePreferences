*****************************************************************************
* This is the Analysis do file to analyse Time Preferences in South Africa  *
* and fits within section 8 of the Main do-file.                            * 
*                                                                           *
* Date first generated:             1 September 2025                        *
* Created by:                       Rinelle Chetty                          * 
*****************************************************************************

* Set size of LL reward
local LL "50"

* Set ylabel scaling for all subsequent graphs
local ylabel "44(1)50"

* Reset
mylabels 44(1)50, myscale(@) prefix($) format(%4.2f) local(ylabel)


* Get time preference estimates linked to JHU data
use legend, clear
generate int day = day(date)
generate int month = month(date)


* Retain months we have experiments for and regenerate
drop if month<5

* Generate the bar legends
replace s = uniform()*50
twoway (bar s date if dec_c_us == 1, sort fcolor(blue*0.05) lcolor(blue*0.05)) || ///
(bar s date if dec_c_us == 2, sort fcolor(blue*0.15) lcolor(blue*0.15)) || ///
(bar s date if dec_c_us == 3, sort fcolor(blue*0.25) lcolor(blue*0.25)) || ///
(bar s date if dec_c_us == 4, sort fcolor(blue*0.35) lcolor(blue*0.35)) || ///
(bar s date if dec_c_us == 5, sort fcolor(blue*0.45) lcolor(blue*0.45)) || ///
(bar s date if dec_c_us == 6, sort fcolor(blue*0.55) lcolor(blue*0.55)) || ///
(bar s date if dec_c_us == 7, sort fcolor(blue*0.65) lcolor(blue*0.65)) || ///
(bar s date if dec_c_us == 8, sort fcolor(blue*0.75) lcolor(blue*0.75)) || ///
(bar s date if dec_c_us == 9, sort fcolor(blue*0.85) lcolor(blue*0.85)) || ///
(bar s date if dec_c_us == 10, sort fcolor(blue*0.95) lcolor(blue*0.95)), ///
legend(off) ytitle("") l1title("Dollars", color(white) orientation(horizontal)) ///
ylabel(`ylabel', labcolor(white) angle(horizontal) tlcolor(white)) ///
xtitle("") xlabel(none, nolabels noticks) fysize(7.5) ///
saving(c_us_bar, replace)

twoway (bar s date if dec_d_us == 1, sort fcolor(red*0.05) lcolor(red*0.05)) || ///
(bar s date if dec_d_us == 2, sort fcolor(red*0.15) lcolor(red*0.15)) || ///
(bar s date if dec_d_us == 3, sort fcolor(red*0.25) lcolor(red*0.25)) || ///
(bar s date if dec_d_us == 4, sort fcolor(red*0.35) lcolor(red*0.35)) || ///
(bar s date if dec_d_us == 5, sort fcolor(red*0.45) lcolor(red*0.45)) || ///
(bar s date if dec_d_us == 6, sort fcolor(red*0.55) lcolor(red*0.55)) || ///
(bar s date if dec_d_us == 7, sort fcolor(red*0.65) lcolor(red*0.65)) || ///
(bar s date if dec_d_us == 8, sort fcolor(red*0.75) lcolor(red*0.75)) || ///
(bar s date if dec_d_us == 9, sort fcolor(red*0.85) lcolor(red*0.85)) || ///
(bar s date if dec_d_us == 10, sort fcolor(red*0.95) lcolor(red*0.95)), ///
legend(off) ytitle("") l1title("Dollars", color(white) orientation(horizontal)) ///
ylabel(`ylabel', labcolor(white) angle(horizontal) tlcolor(white)) ////
xtitle("") xlabel(none, noticks nolabels) fysize(7.5) ///
saving(d_us_bar, replace)

* Estimated pre-COVID
local exp_pre = 48.03
local qh_pre = 47.07


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
local exp_pre = 48.03
local qh_pre = 47.07
marginsplot using pvExpQH50margin, l1title("Dollars", orientation(horizontal)) ///
ytitle("") title("") xlabel("", format(%tdm)) xtitle("") ///
plot1opts(lwidth(`wide') lcolor(`exp_color') mcolor(`exp_color')) ///
ci1opts(lcolor(`exp_color')) ///
plot2opts(lwidth(`wide') lpattern(dash) lcolor(`qh_color') mcolor(`qh_color')) ///
ci2opts(lcolor(`qh_color')) ///
legend(order(3 "Exponential" 4 "Quasi-Hyperbolic") size(medlarge) cols(1) ring(0) pos(2) nobox) ///
yline(50, lcolor(gs14) lwidth(vthick)) ///
yline(`exp_pre', lcolor(green%30) lwidth(vvvthick) lpattern(solid)) ///
yline(`qh_pre', lcolor(green%30) lwidth(vvvthick) lpattern(dash)) ///
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
*******************************************************************************
***   MIXTURE MODEL OF EXPONENTIAL AND HYPERBOLIC DISCOUNTING               ***
***   (RDU utility, joint estimation on the risk and time choices)          ***
***                                                                         ***
***   This is PART B for the time preferences chapter: the estimation code  ***
***   that drives  ml_rdu_discount_mixed , which is already in              ***
***   MLfunctionsTime.do. It follows the same sequence as Part B of         ***
***   MixtureEUTRDU.do in the risk chapter:                                 ***
***                                                                         ***
***     1. fit each single model for starting values                        ***
***     2. fit the mixture from those values                                ***
***     3. sanity checks, including other starting values                   ***
***     4. export the table                                                 ***
***     5. compare the mixture with the single models (no p-value)          ***
***     6. an interior test on the two discounting parameters               ***
***     7. wave effects on the mixing probability                           ***
***                                                                         ***
***   It is self-contained: nothing here depends on estimates left in       ***
***   memory by another do-file.                                            ***
***                                                                         ***
***   FOUR THINGS TO CHECK BEFORE THE FIRST RUN (all marked CHECK below):   ***
***     a. the name of the data file                                        ***
***     b. the names of the time variables in  global timevars              ***
***     c. the starting values for the two single-model fits                ***
***     d. the output folders used by  esttab  and  margins                 ***
*******************************************************************************

* Start log file 
cap log close 
log using "Log_Time_3_Mixture.txt", text replace

set more off

* CHECK (a): the data file, with risk and time choices stacked in one file.
use "timedata.dta", clear

/* Only needed to run this file on its own, outside the main do-file.
   Uncomment, and drop again once the file is wired into the main run.

clear all
set more off
global stata_tables "stata_tables"
global margins      "margins"
qui do "MLfunctionsTime.do"
use "timedata.dta", clear
*/


*-------------------------------------------------------------------*
*   Settings                                                        *
*-------------------------------------------------------------------*

/* ml_rdu_discount_mixed reads $cdf, $ufunc and $weigh. It does not read
   $discount (the two discounting functions are fixed: exponential and Mazur
   hyperbolic) and it does not read $error (risk choices always get contextual
   errors, time choices always get Fechner errors).

   With $weigh set to prelec2 the equations must be handed to -ml model- in
   this order:
       r  phi  eta  deltaE  deltaH  noiseRA  noiseDR  kappa
   Every parameter is in LEVELS. There is no log version of this program, so
   nothing keeps phi, eta, the two deltas or the two noise terms positive.
   The sanity checks below test for that.

   The dependent variables must be in this order, because the program reads
   them by position ($ML_y1 to $ML_y20):
       choice, six probabilities, six prizes, uMax, uMin,
       risk, ssamount, ssdelay, llamount, lldelay
   Delays are in days: the program divides them by 365. A time choice is coded
   choice = 0 for smaller-sooner and choice = 1 for larger-later. */

global cdf          "invlogit"
global maxtech      "nr"
global ufunc        "crra"
global weigh        "prelec2"
global riskvars     "prob1L prob2L prob3L prob1R prob2R prob3R prize1L prize2L prize3L prize1R prize2R prize3R uMax uMin"

* CHECK (b): risk indicator, SS amount, SS delay, LL amount, LL delay,
* in exactly that order.
global timevars     "risk ssamount ssdelay llamount lldelay"

/* Covariates. Three separate lists, all empty for the homogeneous model:
     demog      goes on r, phi, eta, deltaE and deltaH
     hetero     goes on the two noise equations
     kappavars  goes on the mixing equation only
   demog is deliberately NOT put on kappa. In the risk file it was, which
   would have put the demographics on the mixing probability as well. */
global demog        ""
global hetero       ""
global kappavars    ""

* The individual time variable names, for the checks further down
local rsk : word 1 of $timevars
local lld : word 5 of $timevars


*-------------------------------------------------------------------*
*   Starting values                                                 *
*-------------------------------------------------------------------*

/* Always start the mixture from the two single-model fits: RDU with
   exponential discounting and RDU with hyperbolic discounting.

   Exponential and Mazur hyperbolic are not nested in each other, which makes
   this mixture easier to maximise than EUT and RDU. The difficulty here is
   different: over short horizons the two discount functions look alike, so
   the mixing probability can be weakly identified. See check 4 below.

   CHECK (c): the single fits are started from the values below, in the order
       r  phi  eta  delta  noiseRA  noiseDR
   The first three and noiseRA come from the RDU fit in the risk chapter.
   delta and noiseDR are guesses. If your existing RDU discounting do-file has
   init() values that are known to converge on these data, paste them here
   instead. Set a local to "" to let -ml- search for its own start. */

local startE "0.35 0.49 0.88 0.50 0.14 1.00"
local startH "0.35 0.49 0.88 0.50 0.14 1.00"

* RDU with exponential discounting on its own
global discount "exp"
ml model lf ml_rdu_discount_flex (r: choice $riskvars $timevars =) (phi:) (eta:) ///
    (delta:) (noiseRA: $hetero) (noiseDR: $hetero), cluster(id) technique($maxtech)
if "`startE'" != "" ml init `startE', copy
ml maximize, difficult nolog
di "RDU + exponential alone: converged = " e(converged) ", ll = " e(ll) ", N = " e(N)
estimates store mE, title(RDU with exponential discounting)
local llE = e(ll)
local NE  = e(N)
foreach p in r phi eta noiseRA noiseDR {
    local `p'E = [`p']_b[_cons]
}
local dE = [delta]_b[_cons]

* RDU with hyperbolic (Mazur) discounting on its own
global discount "mazur"
ml model lf ml_rdu_discount_flex (r: choice $riskvars $timevars =) (phi:) (eta:) ///
    (delta:) (noiseRA: $hetero) (noiseDR: $hetero), cluster(id) technique($maxtech)
if "`startH'" != "" ml init `startH', copy
ml maximize, difficult nolog
di "RDU + hyperbolic alone:  converged = " e(converged) ", ll = " e(ll) ", N = " e(N)
estimates store mH, title(RDU with hyperbolic discounting)
local llH = e(ll)
foreach p in r phi eta noiseRA noiseDR {
    local `p'H = [`p']_b[_cons]
}
local dH = [delta]_b[_cons]

/* The mixture has one r, one weighting function and one noise term of each
   kind, shared by the two discounting rules. Start those from whichever
   single model fits better. deltaE comes from the exponential fit and deltaH
   from the hyperbolic fit. kappa = 0 is a 50:50 mixture. */
local best = cond(`llH' >= `llE', "H", "E")
foreach p in r phi eta noiseRA noiseDR {
    local `p'0 = ``p'`best''
}

di "shared start values taken from single model `best'"
di "start values: `r0' `phi0' `eta0' `dE' `dH' `noiseRA0' `noiseDR0' 0"


*-------------------------------------------------------------------*
*   The mixture                                                     *
*-------------------------------------------------------------------*

/* The weight on the exponential rule is invlogit(-kappa): kappa = 0 is a
   50:50 mixture and large positive kappa puts all the weight on hyperbolic.
   The mixture applies to the time choices only. The risk choices enter the
   likelihood through RDU alone and pin down r, phi, eta and noiseRA.

   Starting values are set by equation name, not with init(..., copy), so
   this still works when the covariate globals are filled in: the covariate
   coefficients simply start at zero. */

ml model lf ml_rdu_discount_mixed (r: choice $riskvars $timevars = $demog)      ///
    (phi: $demog) (eta: $demog) (deltaE: $demog) (deltaH: $demog)               ///
    (noiseRA: $hetero) (noiseDR: $hetero) (kappa: $kappavars),                  ///
    cluster(id) technique($maxtech)
ml init r:_cons=`r0' phi:_cons=`phi0' eta:_cons=`eta0' deltaE:_cons=`dE'        ///
    deltaH:_cons=`dH' noiseRA:_cons=`noiseRA0' noiseDR:_cons=`noiseDR0'         ///
    kappa:_cons=0
ml maximize, difficult iterate(150)

di as error "converged = " e(converged)
estimates store mEH, title(Exponential/Hyperbolic mixture)
local llM = e(ll)

* Share of time choices made by each rule. Only meaningful as written when
* there are no covariates on kappa.
nlcom (probExp: invlogit(-[kappa]_b[_cons])) (probHyp: invlogit([kappa]_b[_cons]))


*-------------------------------------------------------------------*
*   Sanity checks. Read these before believing the estimates.       *
*-------------------------------------------------------------------*

/* 1. Did it converge? A non-converged fit can report a log likelihood that is
      not comparable with anything else, so never rank models on e(ll) without
      checking this first. The three models must also use the same sample. */
if e(converged) != 1 {
    di as error "MIXTURE DID NOT CONVERGE: do not use these estimates"
}
if e(N) != `NE' {
    di as error "mixture N = " e(N) " but single-model N = `NE': samples differ"
}

/* 2. Are the parameters in range? Everything is estimated in levels, so
      nothing stops a weighting parameter, a discount rate or a noise term
      going to zero or below. A flat weighting function (phi near 0) means the
      RDU part has stopped using the probabilities, as in the risk chapter. */
foreach p in phi eta deltaE deltaH noiseRA noiseDR {
    if [`p']_b[_cons] <= 0 {
        di as error "`p' = " [`p']_b[_cons] " is not positive: estimate is not meaningful"
    }
}
if [phi]_b[_cons] < 0.05 {
    di as error "phi = " [phi]_b[_cons] " is close to 0: weighting function is flat"
}

/* 3. Has the mixture collapsed onto one rule? If the share of either rule is
      close to 0, the discount rate of that rule is barely identified and its
      standard error will be very large. */
local pE = invlogit(-[kappa]_b[_cons])
if `pE' < 0.02 | `pE' > 0.98 {
    di as error "share of exponential choices = " %6.4f `pE' ": the mixture is close to a single rule"
}

/* 4. Can the data tell the two rules apart? Compare the discount factors each
      rule implies at the shortest, median and longest LL delay in the data.
      If the two columns are nearly the same at every horizon, the two rules
      are making the same predictions and kappa is not well identified. */
quietly summarize `lld' if `rsk' == 0, detail
di "Discount factors implied by each rule at the horizons in the data:"
foreach d in `r(min)' `r(p50)' `r(max)' {
    di "  LL delay " %6.0f `d' " days:  exponential " %6.4f (1/((1 + [deltaE]_b[_cons])^(`d'/365))) ///
       "   hyperbolic " %6.4f (1/(1 + [deltaH]_b[_cons]*`d'/365))
}

/* 5. Try other starting values. All the sensible ones should land on the same
      log likelihood and the same deltas. If one lands somewhere better, check
      its e(converged) and e(N) before believing it.
        alt1  the two discount rates swapped
        alt2  start mostly hyperbolic
        alt3  start mostly exponential
        alt4  the two discount rates pulled apart */
local dElow  = `dE'/2
local dHhigh = `dH'*2
local alt1 "`r0' `phi0' `eta0' `dH' `dE' `noiseRA0' `noiseDR0' 0"
local alt2 "`r0' `phi0' `eta0' `dE' `dH' `noiseRA0' `noiseDR0' 1.5"
local alt3 "`r0' `phi0' `eta0' `dE' `dH' `noiseRA0' `noiseDR0' -1.5"
local alt4 "`r0' `phi0' `eta0' `dElow' `dHhigh' `noiseRA0' `noiseDR0' 0"

forvalues a = 1/4 {
    capture {
        ml model lf ml_rdu_discount_mixed (r: choice $riskvars $timevars =) (phi:) ///
            (eta:) (deltaE:) (deltaH:) (noiseRA:) (noiseDR:) (kappa:),             ///
            cluster(id) technique($maxtech)
        ml init `alt`a'', copy
        ml maximize, difficult nolog iterate(150)
    }
    if _rc di as error "alt`a' (`alt`a''): failed, rc = " _rc
    else    di "alt`a': ll = " %14.6f e(ll) ", converged = " e(converged)          ///
               ", deltaE = " %7.4f [deltaE]_b[_cons] ", deltaH = " %7.4f [deltaH]_b[_cons] ///
               ", kappa = " %7.4f [kappa]_b[_cons]
}

* The loop leaves the last alternative fit active, so go back to the mixture.
estimates restore mEH


*-------------------------------------------------------------------*
*   Table                                                           *
*-------------------------------------------------------------------*

/* The risk file refits the mixture with phi and eta in levels at this point.
   That step is not needed here: this program already estimates everything in
   levels, so the standard errors are on the parameters themselves. */

* CHECK (d): the folder must exist.
esttab mEH using "$stata_tables/mixture_exp_hyp.rtf",                    ///
    replace label se mtitle("All waves") b(%9.3f) se(%9.3f)              ///
    title(Mixture of exponential and hyperbolic discounting - homogenous preferences)


*-------------------------------------------------------------------*
*   Does the mixture beat each discounting model on its own?        *
*-------------------------------------------------------------------*

/* The mixture nests both single models: it adds one discount rate and kappa
   to each. Report the three log likelihoods, but do NOT read a likelihood
   ratio test on them as a test of whether a second discounting rule exists.
   Under the null that every choice is exponential, deltaH has dropped out of
   the likelihood and is unidentified (and the same for deltaE under the null
   that every choice is hyperbolic). This is Davies' problem, exactly as in the
   risk chapter, and the usual chi-squared critical values do not apply. The
   same goes for the z statistics that nlcom reports on probExp and probHyp,
   which test against a share of zero.

   What IS an ordinary test is the z on kappa in the mixture table: kappa = 0
   is a 50:50 split, which is in the interior of the parameter space.

   The two single models have the same number of parameters and are not
   nested in each other, so their log likelihoods can be compared directly. */

di "RDU + exponential alone: ll = " %12.3f `llE'
di "RDU + hyperbolic alone:  ll = " %12.3f `llH'
di "Mixture:                 ll = " %12.3f `llM'
di "2 x (mixture - exponential) = " %9.1f 2*(`llM' - `llE') "   (no chi-squared p-value)"
di "2 x (mixture - hyperbolic)  = " %9.1f 2*(`llM' - `llH') "   (no chi-squared p-value)"


*-------------------------------------------------------------------*
*   Do the two rules carry the same discounting parameter?          *
*-------------------------------------------------------------------*

/* The risk file asks here whether the two rules need different utility
   functions (rEUT = rRDU). That question does not arise in this model: the
   program already imposes one r and one noiseDR on both discounting rules.

   The nearest interior restriction is deltaE = deltaH. It does not make the
   two rules the same, because the functional forms still differ, so both
   components remain present and identified under the null and this is an
   ordinary Wald test. Read it narrowly: deltaE is an annual exponential rate
   and deltaH is Mazur's k, so equality only says the two rules discount a
   one-year delay by the same amount. */

test [deltaE]_cons = [deltaH]_cons


*-------------------------------------------------------------------*
*   Does the share of exponential choices move across waves?        *
*-------------------------------------------------------------------*

/* Putting covariates on the kappa equation lets the mixing probability vary.
   This comparison IS legitimate: the mixture is present under both
   hypotheses, so the identification problem above does not affect it.

   wave goes on kappa only, so that differences across waves show up in the
   mixing share and are not absorbed by the component parameters. Fill in
   global demog  above to add covariates to the component equations.

   The fit is started from the stored mixture, matched by coefficient name. */

estimates restore mEH
matrix b0 = e(b)

ml model lf ml_rdu_discount_mixed (r: choice $riskvars $timevars = $demog)      ///
    (phi: $demog) (eta: $demog) (deltaE: $demog) (deltaH: $demog)               ///
    (noiseRA: $hetero) (noiseDR: $hetero) (kappa: i.wave $kappavars),           ///
    cluster(id) technique($maxtech)
ml init b0, skip
ml maximize, difficult iterate(150)
di as error "converged = " e(converged)
estimates store mEHwave, title(Exponential/Hyperbolic mixture by wave)

* Joint test of the wave effects on the mixing probability. This is the
* headline result. The pairwise tests below are descriptive only.
testparm i.wave, equation(kappa)

* Share of time choices made by the exponential rule, by wave
quietly levelsof wave if e(sample), local(waves)
margins, over(wave) expression(invlogit(-predict(equation(kappa))))  ///
    saving("$margins/mixtureEXPshare", replace) post

/* Pairwise differences. With six waves there are fifteen of these, so a few
   will come in under 5% by chance. Do not report them as wave effects unless
   the joint test above rejects. */
foreach i in `waves' {
    foreach j in `ferest()' {
        test `i'.wave == `j'.wave
        if r(p) < 0.05 {
            di as error "waves `i' and `j' differ at 5%"
        }
    }
}

* -margins, post- has replaced the active estimates. Put the wave model back
* before doing anything else with it.
estimates restore mEHwave

*******************************************************************************

log close 
di as error "END of DO-FILE: EXPONENTIAL/HYPERBOLIC MIXTURE"

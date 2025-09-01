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
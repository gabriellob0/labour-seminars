clear all
set more off

qui do "mle_programs.do"

use Mas_Pallais_2017, clear

egen tag = tag (wagegap)
bys wagegap: egen mm1 = mean (chose_position1)
replace mm1 = . if tag != 1

me_correction
mylogit_mle1 ${error}

ml model lf mylogit_mle2 (chose_position1 = wagegap), vce(robust)
qui ml maximize, iterate(100)

save ps2, replace
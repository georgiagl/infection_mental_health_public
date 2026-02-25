*Cleaning data for Table 1 
*Author: Sl Cadogan

*encode variables

foreach var in gender ethnicity_4_na smoking_status_simp indexdate_age_cat obese_na cci_cat{
    encode `var', gen(`var'_encoded)
}


drop gender
rename gender_encoded gender
drop indexdate_age_cat
rename indexdate_age_cat_encoded indexdate_age_cat
drop ethnicity_4_na
rename ethnicity_4_na_encoded ethnicity_4_na
drop cci_cat
rename cci_cat_encoded cci_cat
drop obese_na
rename obese_na_encoded obese_na
drop smoking_status_simp
rename smoking_status_simp_encoded smoking_status_simp

*imd
tab imd, nolabel
encode imd, gen(imd_encoded)
drop imd
rename imd_encoded imd

*change missing to unknown for smoking, obesity and ethnicity
replace ethnicity_4_na=4 if ethnicity_4_na==.
replace smoking_status_simp=3 if smoking_status_simp==.
label define smoking_status_simp 1 "Current-or ex-smoker " 2 "Non-smoker" 3"Unknown"
lab values smoking_status_simp smoking_status_simp
tab obese_na, m
replace obese_na=3 if obese_na==.
label define obese_na 1"Not obese" 2"Obese" 3"Unknown"
lab values obese_na obese_na
replace imd=6 if imd==.
label define imd 1"1, least deprived" 2"2"3"3"4"4"5"5, most deprived"
lab values imd imd

* drop incorrect obese cat without missing data and rename correct obese
drop obese
rename obese_na obese

* change unknown categories from descriptives to missing
replace ethnicity_4_na=. if ethnicity_4_na==4
replace smoking_status_simp=. if smoking_status_simp==3
replace obese=. if obese==3
replace imd=. if imd==6
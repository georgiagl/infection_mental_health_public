

*Do file to create Table 1
*Table title: Baseline characteristics of samples by cohort and infection status
*Author: Sl Cadogan
*last modified: 16-09-25

***************************************************
clear 
*append data sets for descriptives

*drop matched sets without at least one exposure or unexposure
unique patid

tab exposed
sort setid
egen group_sum_exposed = total(exposed), by(setid)
browse setid patid exposed if group_sum_exposed==0
tab group_sum_exposed


gen unexposed=0
replace unexposed=1 if exposed==0
tab unexposed

sort setid
egen group_sum_unexposed = total(unexposed), by(setid)
browse setid patid exposed if group_sum_unexposed==0
tab group_sum_unexposed


*drop if no exposed in matched set
drop if group_sum_exposed==0


*drop if no unexposed in matched set
drop if group_sum_unexposed==0

unique patid
*save temp GE file
save temp_GE

clear


*drop matched sets without at least one exposure or unexposure

unique patid


tab exposed
sort setid
egen group_sum_exposed = total(exposed), by(setid)
browse setid patid exposed if group_sum_exposed==0
tab group_sum_exposed


gen unexposed=0
replace unexposed=1 if exposed==0
tab unexposed

sort setid
egen group_sum_unexposed = total(unexposed), by(setid)
browse setid patid exposed if group_sum_unexposed==0
tab group_sum_unexposed


*drop if no exposed in matched set
drop if group_sum_exposed==0

*drop if no unexposed in matched set
drop if group_sum_unexposed==0

unique patid

*save temp LRTI file
save temp_LRTI


**************
*SSTI
**************
clear


*drop matched sets without an exposure
unique patid

tab exposed
sort setid
egen group_sum_exposed = total(exposed), by(setid)
browse setid patid exposed if group_sum_exposed==0
tab group_sum_exposed


gen unexposed=0
replace unexposed=1 if exposed==0
tab unexposed

sort setid
egen group_sum_unexposed = total(unexposed), by(setid)
browse setid patid exposed if group_sum_unexposed==0
tab group_sum_unexposed


*drop if no exposed in matched set
drop if group_sum_exposed==0

*drop if no unexposed in matched set
drop if group_sum_unexposed==0

unique patid

*save temp SSTI file
save temp_SSTI

************************
*Sepsis
************************
clear

*drop matched sets without an exposure
unique patid


tab exposed
sort setid
egen group_sum_exposed = total(exposed), by(setid)
browse setid patid exposed if group_sum_exposed==0
tab group_sum_exposed


gen unexposed=0
replace unexposed=1 if exposed==0
tab unexposed

sort setid
egen group_sum_unexposed = total(unexposed), by(setid)
browse setid patid exposed if group_sum_unexposed==0
tab group_sum_unexposed

*drop if no exposed in matched set
drop if group_sum_exposed==0

*drop if no unexposed in matched set
drop if group_sum_unexposed==0
unique patid

*save temp Sepsis file
save temp_Sepsis


*********************
*UTI
*********************
clear

*drop matched sets without an exposure
unique patid

tab exposed
sort setid
egen group_sum_exposed = total(exposed), by(setid)
browse setid patid exposed if group_sum_exposed==0
tab group_sum_exposed


gen unexposed=0
replace unexposed=1 if exposed==0
tab unexposed

sort setid
egen group_sum_unexposed = total(unexposed), by(setid)
browse setid patid exposed if group_sum_unexposed==0
tab group_sum_unexposed

*drop if no exposed in matched set
drop if group_sum_exposed==0

*drop if no unexposed in matched set
drop if group_sum_unexposed==0

unique patid

*save temp UTI file
save temp_UTI

*************************
*Mening/encephalitis
**************************
clear

*drop matched sets without an exposure
unique patid

tab exposed
sort setid
egen group_sum_exposed = total(exposed), by(setid)
browse setid patid exposed if group_sum_exposed==0
tab group_sum_exposed


gen unexposed=0
replace unexposed=1 if exposed==0
tab unexposed

sort setid
egen group_sum_unexposed = total(unexposed), by(setid)
browse setid patid exposed if group_sum_unexposed==0
tab group_sum_unexposed

*drop if no exposed in matched set
drop if group_sum_exposed==0

*drop if no unexposed in matched set
drop if group_sum_unexposed==0
* o observations dropped

unique patid

*save temp mening/enceph file
save temp_mening_enceph


*append temp files
clear
use "temp_GE"
append using "temp_LRTI"
append using "temp_SSTI"
append using "temp_Sepsis"
append using "temp_UTI"
append using "temp_mening_enceph"

tab inf, nolabel



*run cleaning do



*to change order of infection/no infection in table
tab exposed
replace exposed=2 if exposed==0
tab exposed
label define exposed 1"with infection" 2"without_infection"
lab val exposed exposed
tab exposed


*change order of infections
encode inf, gen(inf_encoded)
tab inf_encoded
tab inf_encoded, nolabel
recode inf_encoded (6=5) (5=6), gen(inf_encoded_swapped)
tab inf_encoded_swapped
label define inf_encoded_swapped 1"GE" 2"LRTI" 3"SSTI" 4"UTI" 5"sepsis" 6"meningitis/encephalitis", replace
label values inf_encoded_swapped inf_encoded_swapped
drop inf_encoded
drop inf
rename inf_encoded_swapped inf


sort inf exposed

egen inf_exposed= group (inf exposed), label
tab inf_exposed



*check included variables
summ indexdate_age, detail
tab2 inf_exposed indexdate_age_cat, row
tab gender

tab ethnicity_4_na, m
rename ethnicity_4_na Ethnicity
tab smoking_status_simp, m
tab obese_na, m
rename indexdate_age Age

 
dtable, by(inf_exposed, tests testnotes nototal) ///
continuous(Age, statistics(median q1 q3)///
factor(indexdate_age_cat gender ethnicity smoking_status_simp obese cci_cat imd, statistics(fvfrequency fvproportion))  ///
nformat(%9.1f) sformat("(%9.2f)")///
export(baseline_table.xls, replace)

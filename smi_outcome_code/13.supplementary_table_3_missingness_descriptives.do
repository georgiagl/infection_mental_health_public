
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



generate missingness_model4 = !missing(exposed, imd, alcohol_abuse, smoking_status_simp,obese, indexdate_age,indexdate_age_cat, cci_cat, ethnicity_4_na, gender)
label define missingness_model4 0 "Excluded from model 4" 1 "In model 4 (complete case)"
label values missingness_model4 missingness_model4


egen grp = group(missingness_model4 inf), label

*v1
dtable, by(grp, tests testnotes nototal) sample(,statistic(frequency percent) place(seplabels)) continuous(Age, statistics(median q1 q3) factor(indexdate_age_cat gender, statistics(fvfrequency) test(fisher)format("(%s)")  export(missingdata.xls, replace)



tabstat indexdate_age, by(missingness_model4) statistics(mean sd n)

tabulate missingness_model4 gender, row

*Main cox regression analysis

*Notes: Indexdate is infection date for exposed or matching data for unexposed
*Enddate:  earliest of: end of registration, date of death, date of the last data collection from the practice, or the end of the study (2024-06-15)

*******************************************
*GE  
*******************************************
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



*clean data for cox regression analysis



drop if organic_psychosis_obsdate <=indexdate
drop if infection_censor_obsdate<=indexdate
drop if smi_incident_obsdate<=indexdate
drop if enddate<=indexdate

*Set the data for Cox regression: censoring variables (for SMI): organic_psychosis_obsdate,  infection censor: infection_censor_obsdate  
gen end_date=min(enddate, smi_incident_obsdate, infection_censor_obsdate, organic_psychosis_obsdate)
format end_date %d  
 
*generate smi_status variable
gen smi_status=1 if smi_incident_obsdate==end_date
tab smi_status

stset end_date, failure(smi_status=1) enter (indexdate) origin (indexdate) id(patid)scale(365.25)


*1. crude cox regression
stcox i.exposed ,strata (setid)

*2. +IMD
stcox i.exposed i.imd, strata (setid)

*3. + lifestlye and co-morbidities
stcox i.exposed i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat, strata (setid)

*4. Fully adjusted cox regression (adjusted for IMD, alcohol consumption, smoking status, obesity status, CCI score cat, + ethnicity)
stcox i.exposed i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat i.ethnicity_4_na, strata (setid)


*person time and rates
stptime, per(100000)


*Followup time for table 1
bysort exposed: summ _t, d


*age as effect modifier
*4. Fully adjusted cox regression (adjusted for IMD, alcohol consumption, smoking status, obesity status, CCI score cat, + ethnicity)
gen exposed_age = exposed*indexdate_age_cat

stcox i.exposed_age i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat i.ethnicity_4_na, strata (setid)


*sex as an effect modifier
*4. Fully adjusted cox regression (adjusted for IMD, alcohol consumption, smoking status, obesity status, CCI score cat, + ethnicity)
gen exposed_gender = exposed*gender


stcox i.exposed_gender i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat i.ethnicity_4_na, strata (setid)

*frailty as effect modifier
*keep aged 65-95 year olds
keep if indexdate_age>64
keep if indexdate_age<96

encode frailty_score, gen (frailty_score_encoded)
gen exposed_frailty=exposed*frailty_score_encoded

stcox i.exposed_frailty i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat i.ethnicity_4_na, strata (setid)



log close


*******************************************
*LRTI
***********************************************
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

*run cleaning dofile

*cox regression analysis
drop if organic_psychosis_obsdate <=indexdate
drop if infection_censor_obsdate<=indexdate
drop if smi_incident_obsdate<=indexdate
drop if enddate<=indexdate


*Set the data for Cox regression: censoring variables (for SMI): organic_psychosis_obsdate,  infection censor: infection_censor_obsdate  
gen end_date=min(enddate, smi_incident_obsdate,infection_censor_obsdate, organic_psychosis_obsdate)
format end_date %d  

*generate smi_status variable
gen smi_status=1 if smi_incident_obsdate==end_date
tab smi_status 

stset end_date, failure(smi_status=1) enter (indexdate) origin (indexdate) id(patid)scale(365.25)


*crude cox regression
stcox i.exposed ,strata (setid)

*2. +IMD
stcox i.exposed i.imd, strata (setid)

*3. + lifestlye and co-morbidities
stcox i.exposed i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat, strata (setid)

*4. Fully adjusted cox regression (adjusted for IMD, alcohol consumption, smoking status, obesity status, CCI score cat, + ethnicity)
stcox i.exposed i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat i.ethnicity_4_na, strata (setid)

*person time and rates
stptime, per(100000)

*age as effect modifier
*4. Fully adjusted cox regression (adjusted for IMD, alcohol consumption, smoking status, obesity status, CCI score cat, + ethnicity)
gen exposed_age = exposed*indexdate_age_cat

stcox i.exposed_age i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat i.ethnicity_4_na, strata (setid)

*sex as an effect modifier
*4. Fully adjusted cox regression (adjusted for IMD, alcohol consumption, smoking status, obesity status, CCI score cat, + ethnicity)
gen exposed_gender = exposed*gender


stcox i.exposed_gender i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat i.ethnicity_4_na, strata (setid)

*frailty as effect modifier
*keep aged 65-95 year olds
keep if indexdate_age>64
keep if indexdate_age<96

encode frailty_score, gen (frailty_score_encoded)
gen exposed_frailty=exposed*frailty_score_encoded

stcox i.exposed_frailty i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat i.ethnicity_4_na, strata (setid)

log close

********************************************************
*SSTI 
********************************************************
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

*run cleaning dofile

*cox regression analysis
drop if organic_psychosis_obsdate <=indexdate
drop if infection_censor_obsdate<=indexdate
drop if smi_incident_obsdate<=indexdate
drop if enddate<=indexdate


*Set the data for Cox regression: censoring variables (for SMI): organic_psychosis_obsdate,  infection censor: infection_censor_obsdate  
gen end_date=min(enddate, smi_incident_obsdate, infection_censor_obsdate, organic_psychosis_obsdate)
format end_date %d  
 
*generate smi_status variable
gen smi_status=1 if smi_incident_obsdate==end_date
tab smi_status

stset end_date, failure(smi_status=1) enter (indexdate) origin (indexdate) id(patid)scale(365.25)

*crude cox regression
stcox i.exposed ,strata (setid)

*2. +IMD
stcox i.exposed i.imd, strata (setid)

*3. + lifestlye and co-morbidities
stcox i.exposed i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat, strata (setid)

*4. Fully adjusted cox regression (adjusted for IMD, alcohol consumption, smoking status, obesity status, CCI score cat, + ethnicity)
stcox i.exposed i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat i.ethnicity_4_na, strata (setid)

*follow-up time
summ _t, detail

*person time and rates
stptime, per(100000)

*age as effect modifier
*4. Fully adjusted cox regression (adjusted for IMD, alcohol consumption, smoking status, obesity status, CCI score cat, + ethnicity)
gen exposed_age = exposed*indexdate_age_cat

stcox i.exposed_age i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat i.ethnicity_4_na, strata (setid)

*sex as an effect modifier
*4. Fully adjusted cox regression (adjusted for IMD, alcohol consumption, smoking status, obesity status, CCI score cat, + ethnicity)
gen exposed_gender = exposed*gender


stcox i.exposed_gender i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat i.ethnicity_4_na, strata (setid)

*frailty as effect modifier
*keep aged 65-95 year olds
keep if indexdate_age>64
keep if indexdate_age<96

encode frailty_score, gen (frailty_score_encoded)
gen exposed_frailty=exposed*frailty_score_encoded

stcox i.exposed_frailty i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat i.ethnicity_4_na, strata (setid)


*table one follow up time data
bysort exposed: summ _t, d

*table 1 consultations data
bysort exposed: summ n_cons_in_year_pre_index , d

log close


*********************************************
*meningitis/encephalitis (N=
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

unique patid

*run cleaning dofile

*cox regression analysis
*cox regression analysis
drop if organic_psychosis_obsdate <=indexdate
drop if infection_censor_obsdate<=indexdate
drop if smi_incident_obsdate<=indexdate
drop if enddate<=indexdate



*For setting the data for Cox regression, need to also add the censoring variables (for SMI): organic_psychosis_obsdate 

gen end_date=min(enddate, smi_incident_obsdate, infection_censor_obsdate, organic_psychosis_obsdate)
format end_date %d

*generate smi_status variable
gen smi_status=1 if smi_incident_obsdate==end_date
tab smi_status

stset end_date, failure(smi_status=1) enter (indexdate) origin (indexdate) id(patid)scale(365.25)


*crude cox regression
stcox i.exposed ,strata (setid)

*2. +IMD
stcox i.exposed i.imd, strata (setid)

*3. + lifestlye and co-morbidities
stcox i.exposed i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat, strata (setid)


*4. Fully adjusted cox regression (adjusted for IMD, alcohol consumption, smoking status, obesity status, CCI score cat, + ethnicity)
stcox i.exposed i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat i.ethnicity_4_na, strata (setid)


summ _t, d

*table one follow up time data
bysort exposed: summ _t, d



*table 1 consultations data
bysort exposed: summ n_cons_in_year_pre_index , d


*age as effect modifier
*4. Fully adjusted cox regression (adjusted for IMD, alcohol consumption, smoking status, obesity status, CCI score cat, + ethnicity)
gen exposed_age = exposed*indexdate_age_cat

stcox i.exposed_age i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat i.ethnicity_4_na, strata (setid)

*sex as an effect modifier
*4. Fully adjusted cox regression (adjusted for IMD, alcohol consumption, smoking status, obesity status, CCI score cat, + ethnicity)
gen exposed_gender = exposed*gender


stcox i.exposed_gender i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat i.ethnicity_4_na, strata (setid)

*frailty as effect modifier
*keep aged 65-95 year olds
keep if indexdate_age>64
keep if indexdate_age<96

encode frailty_score, gen (frailty_score_encoded)
gen exposed_frailty=exposed*frailty_score_encoded

stcox i.exposed_frailty i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat i.ethnicity_4_na, strata (setid)


log close



***********************************************
*SEPSIS
**************************************************

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


*run cleaning dofile

*cox regression analysis
drop if organic_psychosis_obsdate <=indexdate
drop if infection_censor_obsdate<=indexdate
drop if smi_incident_obsdate<=indexdate
drop if enddate<=indexdate

*For setting the data for Cox regression, need to also add the censoring variables (for SMI): organic_psychosis_obsdate 

gen end_date=min(enddate, smi_incident_obsdate, infection_censor_obsdate, organic_psychosis_obsdate)
format end_date %d

*generate smi_status variable
gen smi_status=1 if smi_incident_obsdate==end_date
tab smi_status

stset end_date, failure(smi_status=1) enter (indexdate) origin (indexdate) id(patid)scale(365.25)

summ _t, d

*crude cox regression
stcox i.exposed ,strata (setid)

*2. +IMD
stcox i.exposed i.imd, strata (setid)

*3. + lifestlye and co-morbidities
stcox i.exposed i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat, strata (setid)

*4. Fully adjusted cox regression (adjusted for IMD, alcohol consumption, smoking status, obesity status, CCI score cat, + ethnicity)
stcox i.exposed i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat i.ethnicity_4_na, strata (setid)


*person time and rates
stptime, per(100000)

*table one follow up time data
bysort exposed: summ _t, d

*table 1 consultations data
bysort exposed: summ n_cons_in_year_pre_index , d



*age as effect modifier
*4. Fully adjusted cox regression (adjusted for IMD, alcohol consumption, smoking status, obesity status, CCI score cat, + ethnicity)
gen exposed_age = exposed*indexdate_age_cat

stcox i.exposed_age i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat i.ethnicity_4_na, strata (setid)

*sex as an effect modifier
*4. Fully adjusted cox regression (adjusted for IMD, alcohol consumption, smoking status, obesity status, CCI score cat, + ethnicity)
gen exposed_gender = exposed*gender


stcox i.exposed_gender i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat i.ethnicity_4_na, strata (setid)

*frailty as effect modifier
*keep aged 65-95 year olds
keep if indexdate_age>64
keep if indexdate_age<96

encode frailty_score, gen (frailty_score_encoded)
gen exposed_frailty=exposed*frailty_score_encoded

stcox i.exposed_frailty i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat i.ethnicity_4_na, strata (setid)



log close


*********************
*UTI
*********************
*/
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


*run cleaning dofile


*cox regression analysis
drop if organic_psychosis_obsdate <=indexdate
drop if infection_censor_obsdate<=indexdate
drop if smi_incident_obsdate<=indexdate
drop if enddate<=indexdate



*For setting the data for Cox regression, need to also add the censoring variables (for SMI): organic_psychosis_obsdate 

gen end_date=min(enddate, smi_incident_obsdate, infection_censor_obsdate, organic_psychosis_obsdate)
format end_date %d

*generate smi_status variable
gen smi_status=1 if smi_incident_obsdate==end_date
tab smi_status

stset end_date, failure(smi_status=1) enter (indexdate) origin (indexdate) id(patid)scale(365.25)


*crude cox regression
stcox i.exposed ,strata (setid)

*2. +IMD
stcox i.exposed i.imd, strata (setid)

*3. + lifestlye and co-morbidities
stcox i.exposed i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat, strata (setid)

*4. Fully adjusted cox regression (adjusted for IMD, alcohol consumption, smoking status, obesity status, CCI score cat, + ethnicity)
stcox i.exposed i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat i.ethnicity_4_na, strata (setid)


*person time and rates
stptime, per(100000)

*overall follow-up 
summ _t, d

*table one follow up time data
bysort exposed: summ _t, d

*table 1 consultations data
bysort exposed: summ n_cons_in_year_pre_index , d


*age as effect modifier
*4. Fully adjusted cox regression (adjusted for IMD, alcohol consumption, smoking status, obesity status, CCI score cat, + ethnicity)
gen exposed_age = exposed*indexdate_age_cat

stcox i.exposed_age i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat i.ethnicity_4_na, strata (setid)

*sex as an effect modifier
*4. Fully adjusted cox regression (adjusted for IMD, alcohol consumption, smoking status, obesity status, CCI score cat, + ethnicity)
gen exposed_gender = exposed*gender


stcox i.exposed_gender i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat i.ethnicity_4_na, strata (setid)

*frailty as effect modifier
*keep aged 65-95 year olds
keep if indexdate_age>64
keep if indexdate_age<96

encode frailty_score, gen (frailty_score_encoded)
gen exposed_frailty=exposed*frailty_score_encoded

stcox i.exposed_frailty i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat i.ethnicity_4_na, strata (setid)



log close

********************************************







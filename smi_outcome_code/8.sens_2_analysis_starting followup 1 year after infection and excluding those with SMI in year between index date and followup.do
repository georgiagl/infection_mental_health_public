*Sensitivity analysis:
*Main analysis starting follow-up 1 year after initial infection, excluding those with SMI in the year between index date and start of follow-up



*Notes: MAIN ANALYSIS: Indexdate is infection date for exposed or matching data for unexposed
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


*Index date is date of infection_censor_date
*Create a new variable for entry after 1 year:
gen entry1yr = indexdate + 365.25
format entry1yr %d

*Set the endate for Cox regression: censoring variables (for SMI): organic_psychosis_obsdate,  infection censor: infection_censor_obsdate

*drop if censors are before indexdate
drop if organic_psychosis_obsdate <=indexdate
drop if infection_censor_obsdate<=indexdate
drop if smi_incident_obsdate<=indexdate

*enddate
gen end_date=min(enddate, smi_incident_obsdate, infection_censor_obsdate, organic_psychosis_obsdate)
format end_date %d  

*generate smi_status variable
gen smi_status=1 if smi_incident_obsdate==end_date
tab smi_status

*Drop those with SMI within first year :
drop if smi_status == 1 & smi_incident_obsdate >= indexdate & smi_incident_obsdate < entry1yr


*set data using new entry
stset end_date, failure(smi_status==1) enter(entry1yr) origin(indexdate) id(patid) scale(365.25)



*1. crude cox regression
stcox i.exposed ,strata (setid)

*2. +IMD
stcox i.exposed i.imd, strata (setid)

*3. + lifestlye and co-morbidities
stcox i.exposed i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat, strata (setid)

*4. Fully adjusted cox regression (adjusted for IMD, alcohol consumption, smoking status, obesity status, CCI score cat, + ethnicity)
stcox i.exposed i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat i.ethnicity_4_na, strata (setid)

*******************************
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

*COX regression

*Index date is date of infection_censor_date
*Create a new variable for entry after 1 year:
gen entry1yr = indexdate + 365.25
format entry1yr %d

*Set the endate for Cox regression: censoring variables (for SMI): organic_psychosis_obsdate,  infection censor: infection_censor_obsdate

*drop if censors are before indexdate
drop if organic_psychosis_obsdate <=indexdate
drop if infection_censor_obsdate<=indexdate
drop if smi_incident_obsdate<=indexdate

*enddate
gen end_date=min(enddate, smi_incident_obsdate, infection_censor_obsdate, organic_psychosis_obsdate)
format end_date %d  

*generate smi_status variable
gen smi_status=1 if smi_incident_obsdate==end_date
tab smi_status

*Drop those with SMI within first year :
drop if smi_status == 1 & smi_incident_obsdate >= indexdate & smi_incident_obsdate < entry1yr


*set data using new entry
stset end_date, failure(smi_status==1) enter(entry1yr) origin(indexdate) id(patid) scale(365.25)

*crude cox regression
stcox i.exposed ,strata (setid)


*2. +IMD
stcox i.exposed i.imd, strata (setid)


*3. + lifestlye and co-morbidities
stcox i.exposed i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat, strata (setid)


*4. Fully adjusted cox regression (adjusted for IMD, alcohol consumption, smoking status, obesity status, CCI score cat, + ethnicity)
stcox i.exposed i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat i.ethnicity_4_na, strata (setid)


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


*COX regression

*Index date is date of infection_censor_date
*Create a new variable for entry after 1 year:
gen entry1yr = indexdate + 365.25
format entry1yr %d

*Set the endate for Cox regression: censoring variables (for SMI): organic_psychosis_obsdate,  infection censor: infection_censor_obsdate

*drop if censors are before indexdate
drop if organic_psychosis_obsdate <=indexdate
drop if infection_censor_obsdate<=indexdate
drop if smi_incident_obsdate<=indexdate

*enddate
gen end_date=min(enddate, smi_incident_obsdate, infection_censor_obsdate, organic_psychosis_obsdate)
format end_date %d  

*generate smi_status variable
gen smi_status=1 if smi_incident_obsdate==end_date
tab smi_status

*Drop those with SMI within first year :
drop if smi_status == 1 & smi_incident_obsdate >= indexdate & smi_incident_obsdate < entry1yr


*set data using new entry
stset end_date, failure(smi_status==1) enter(entry1yr) origin(indexdate) id(patid) scale(365.25)

*crude cox regression
stcox i.exposed ,strata (setid)


*2. +IMD
stcox i.exposed i.imd, strata (setid)


*3. + lifestlye and co-morbidities
stcox i.exposed i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat, strata (setid)


*4. Fully adjusted cox regression (adjusted for IMD, alcohol consumption, smoking status, obesity status, CCI score cat, + ethnicity)
stcox i.exposed i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat i.ethnicity_4_na, strata (setid)



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

*COX regression

*Index date is date of infection_censor_date
*Create a new variable for entry after 1 year:
gen entry1yr = indexdate + 365.25
format entry1yr %d

*Set the endate for Cox regression: censoring variables (for SMI): organic_psychosis_obsdate,  infection censor: infection_censor_obsdate

*drop if censors are before indexdate
drop if organic_psychosis_obsdate <=indexdate
drop if infection_censor_obsdate<=indexdate
drop if smi_incident_obsdate<=indexdate

*enddate
gen end_date=min(enddate, smi_incident_obsdate, infection_censor_obsdate, organic_psychosis_obsdate)
format end_date %d  

*generate smi_status variable
gen smi_status=1 if smi_incident_obsdate==end_date
tab smi_status

*Drop those with SMI within first year :
drop if smi_status == 1 & smi_incident_obsdate >= indexdate & smi_incident_obsdate < entry1yr


*set data using new entry
stset end_date, failure(smi_status==1) enter(entry1yr) origin(indexdate) id(patid) scale(365.25)



*crude cox regression
stcox i.exposed ,strata (setid)


*2. +IMD
stcox i.exposed i.imd, strata (setid)


*3. + lifestlye and co-morbidities
stcox i.exposed i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat, strata (setid)



*4. Fully adjusted cox regression (adjusted for IMD, alcohol consumption, smoking status, obesity status, CCI score cat, + ethnicity)
stcox i.exposed i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat i.ethnicity_4_na, strata (setid)


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

*COX regression

*Index date is date of infection_censor_date
*Create a new variable for entry after 1 year:
gen entry1yr = indexdate + 365.25
format entry1yr %d

*Set the endate for Cox regression: censoring variables (for SMI): organic_psychosis_obsdate,  infection censor: infection_censor_obsdate

*drop if censors are before indexdate
drop if organic_psychosis_obsdate <=indexdate
drop if infection_censor_obsdate<=indexdate
drop if smi_incident_obsdate<=indexdate

*enddate
gen end_date=min(enddate, smi_incident_obsdate, infection_censor_obsdate, organic_psychosis_obsdate)
format end_date %d  

*generate smi_status variable
gen smi_status=1 if smi_incident_obsdate==end_date
tab smi_status

*Drop those with SMI within first year :
drop if smi_status == 1 & smi_incident_obsdate >= indexdate & smi_incident_obsdate < entry1yr


*set data using new entry
stset end_date, failure(smi_status==1) enter(entry1yr) origin(indexdate) id(patid) scale(365.25)


*crude cox regression
stcox i.exposed ,strata (setid)


*2. +IMD
stcox i.exposed i.imd, strata (setid)


*3. + lifestlye and co-morbidities
stcox i.exposed i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat, strata (setid)


*4. Fully adjusted cox regression (adjusted for IMD, alcohol consumption, smoking status, obesity status, CCI score cat, + ethnicity)
stcox i.exposed i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat i.ethnicity_4_na, strata (setid)





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



*COX regression

*Index date is date of infection_censor_date
*Create a new variable for entry after 1 year:
gen entry1yr = indexdate + 365.25
format entry1yr %d

*Set the endate for Cox regression: censoring variables (for SMI): organic_psychosis_obsdate,  infection censor: infection_censor_obsdate

*drop if censors are before indexdate
drop if organic_psychosis_obsdate <=indexdate
drop if infection_censor_obsdate<=indexdate
drop if smi_incident_obsdate<=indexdate

*enddate
gen end_date=min(enddate, smi_incident_obsdate, infection_censor_obsdate, organic_psychosis_obsdate)
format end_date %d  

*generate smi_status variable
gen smi_status=1 if smi_incident_obsdate==end_date
tab smi_status

*Drop those with SMI within first year :
drop if smi_status == 1 & smi_incident_obsdate >= indexdate & smi_incident_obsdate < entry1yr


*set data using new entry
stset end_date, failure(smi_status==1) enter(entry1yr) origin(indexdate) id(patid) scale(365.25)


*crude cox regression
stcox i.exposed ,strata (setid)


*2. +IMD
stcox i.exposed i.imd, strata (setid)


*3. + lifestlye and co-morbidities
stcox i.exposed i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat, strata (setid)


*4. Fully adjusted cox regression (adjusted for IMD, alcohol consumption, smoking status, obesity status, CCI score cat, + ethnicity)
stcox i.exposed i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat i.ethnicity_4_na, strata (setid)

log close
 
*Sensitivity analysis: Repeat analyses ending study on 1st March 2020, to account for impact of COVID pandemic



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

drop if organic_psychosis_obsdate <=indexdate
drop if infection_censor_obsdate<=indexdate
drop if smi_incident_obsdate<=indexdate
drop if enddate<=indexdate

*cox regression analysis

*new enddate of 01-05-2020
gen covid_end_date=d(01-05-2020)
format covid_end_date %d

*add covid censor var

*Set the data for Cox regression: censoring variables (for SMI): organic_psychosis_obsdate,  infection censor: infection_censor_obsdate  
gen end_date=min(enddate, covid_end_date, smi_incident_obsdate, infection_censor_obsdate, organic_psychosis_obsdate)
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

*new enddate of 01-05-2020
gen covid_end_date=d(01-05-2020)
format covid_end_date %d

*add covid censor var

*Set the data for Cox regression: censoring variables (for SMI): organic_psychosis_obsdate,  infection censor: infection_censor_obsdate  
gen end_date=min(enddate, covid_end_date, smi_incident_obsdate, infection_censor_obsdate, organic_psychosis_obsdate)
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


*new enddate of 01-05-2020
gen covid_end_date=d(01-05-2020)
format covid_end_date %d

*add covid censor var

*Set the data for Cox regression: censoring variables (for SMI): organic_psychosis_obsdate,  infection censor: infection_censor_obsdate  
gen end_date=min(enddate, covid_end_date, smi_incident_obsdate, infection_censor_obsdate, organic_psychosis_obsdate)
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



*********************************************
*meningitis/encephalitis 
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
drop if organic_psychosis_obsdate <=indexdate
drop if infection_censor_obsdate<=indexdate
drop if smi_incident_obsdate<=indexdate
drop if enddate<=indexdate


*new enddate of 01-05-2020
gen covid_end_date=d(01-05-2020)
format covid_end_date %d

*add covid censor var

*Set the data for Cox regression: censoring variables (for SMI): organic_psychosis_obsdate,  infection censor: infection_censor_obsdate  
gen end_date=min(enddate, covid_end_date, smi_incident_obsdate, infection_censor_obsdate, organic_psychosis_obsdate)
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

*new enddate of 01-05-2020
gen covid_end_date=d(01-05-2020)
format covid_end_date %d

*add covid censor var

*Set the data for Cox regression: censoring variables (for SMI): organic_psychosis_obsdate,  infection censor: infection_censor_obsdate  
gen end_date=min(enddate, covid_end_date, smi_incident_obsdate, infection_censor_obsdate, organic_psychosis_obsdate)
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


*new enddate of 01-05-2020
gen covid_end_date=d(01-05-2020)
format covid_end_date %d
drop if covid_end_date<=indexdate



*add covid censor var

*Set the data for Cox regression: censoring variables (for SMI): organic_psychosis_obsdate,  infection censor: infection_censor_obsdate  
gen end_date=min(enddate, covid_end_date, smi_incident_obsdate, infection_censor_obsdate, organic_psychosis_obsdate)
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

log close


*Secondary analysis

*1) restricting sample to HES linked records and examining severity of infection for LRTI, GE, UTI and SSTI infection cohorts
*******************************************
*GE  
*******************************************
clear 

log using "J:\EHR-Working\Charlotte_Sharon\Infections and SMI\logfiles\new models\GE_severity_HES.smcl", replace

pq use "Z:\GPRD_GOLD\Georgia\Georgia_extract_2024_12_13\processed_data\smi_cohort\20250902\GE.parquet"


*restrict to HES linked data
keep if hes_linkage_eligible== 1
unique patid
*1300352 with linked HES data


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
* 37,995 obsverations dropped


*drop if no unexposed in matched set
drop if group_sum_unexposed==0
*523 observations droped

unique patid
*N=1,262319 


*Check severe flag for analysis
tab severe_flag_composite


*create infection severity variable for analysis
gen infec_sev = .
replace infec_sev = 0 if severe_flag_composite == . & exposed == 0    // No infection
replace infec_sev = 1 if severe_flag_composite == 0 & exposed == 1    // Non-severe infection
replace infec_sev = 2 if severe_flag_composite == 1 & exposed == 1    // Severe infection

label define infec_sev_lbl 0 "No infection" 1 "Non-severe" 2 "Severe"
label values infec_sev infec_sev_lbl
tab infec_sev


*run cleaning dofile
run "J:\EHR-Working\Charlotte_Sharon\Infections and SMI\dofiles\1.cleaning.do" 

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
*4457

*set data for analysis
stset end_date, failure(smi_status=1) enter (indexdate) origin (indexdate) id(patid)scale(365.25)



*1. crude cox regression
stcox i.infec_sev ,strata (setid)


*2. +IMD
stcox i.infec_sev i.imd, strata (setid)

*3. + lifestlye and co-morbidities
stcox i.infec_sev i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat, strata (setid)


*4. Fully adjusted cox regression (adjusted for IMD, alcohol consumption, smoking status, obesity status, CCI score cat, + ethnicity)
stcox i.infec_sev i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat i.ethnicity_4_na, strata (setid)

log close


********************************************************
*SSTI 
********************************************************
log using "J:\EHR-Working\Charlotte_Sharon\Infections and SMI\logfiles\new models\SSTI_severity_HES.smcl", replace

clear

pq use "Z:\GPRD_GOLD\Georgia\Georgia_extract_2024_12_13\processed_data\smi_cohort\20250902\SSTI.parquet"

*restrict to HES linked data
keep if hes_linkage_eligible== 1
unique patid
*1617439

*drop matched sets without an exposure

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
* 44529 obsverations dropped

*drop if no unexposed in matched set
drop if group_sum_unexposed==0
*790 observations dropped

unique patid
*N=1572874


tab severe_flag_composite

*create new infection variable
gen infec_sev = .
replace infec_sev = 0 if severe_flag_composite == . & exposed == 0    // No infection
replace infec_sev = 1 if severe_flag_composite == 0 & exposed == 1    // Non-severe infection
replace infec_sev = 2 if severe_flag_composite == 1 & exposed == 1    // Severe infection

label define infec_sev_lbl 0 "No infection" 1 "Non-severe" 2 "Severe"
label values infec_sev infec_sev_lbl


*run cleaning dofile
run "J:\EHR-Working\Charlotte_Sharon\Infections and SMI\dofiles\1.cleaning.do" 

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
*5145

stset end_date, failure(smi_status=1) enter (indexdate) origin (indexdate) id(patid)scale(365.25)



*1. crude cox regression
stcox i.infec_sev ,strata (setid)


*2. +IMD
stcox i.infec_sev i.imd, strata (setid)

*3. + lifestlye and co-morbidities
stcox i.infec_sev i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat, strata (setid)


*4. Fully adjusted cox regression (adjusted for IMD, alcohol consumption, smoking status, obesity status, CCI score cat, + ethnicity)
stcox i.infec_sev i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat i.ethnicity_4_na, strata (setid)

log close



*******************************************
*LRTI
***********************************************
log using "J:\EHR-Working\Charlotte_Sharon\Infections and SMI\logfiles\new models\LRTI_severity_HES.smcl", replace

clear
pq use "Z:\GPRD_GOLD\Georgia\Georgia_extract_2024_12_13\processed_data\smi_cohort\20250902\LRTI.parquet"

*restrict to HES linked data
keep if hes_linkage_eligible== 1

unique patid
*N=1512440


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
* 243,902 obsverations dropped

*drop if no unexposed in matched set
drop if group_sum_unexposed==0
*1214 observations dropped

unique patid



gen infec_sev = .
replace infec_sev = 0 if severe_flag_composite == . & exposed == 0    // No infection
replace infec_sev = 1 if severe_flag_composite == 0 & exposed == 1    // Non-severe infection
replace infec_sev = 2 if severe_flag_composite == 1 & exposed == 1    // Severe infection

label define infec_sev_lbl 0 "No infection" 1 "Non-severe" 2 "Severe"
label values infec_sev infec_sev_lbl


*run cleaning dofile
run "J:\EHR-Working\Charlotte_Sharon\Infections and SMI\dofiles\1.cleaning.do" 


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
*4,999

stset end_date, failure(smi_status=1) enter (indexdate) origin (indexdate) id(patid)scale(365.25)


*1. crude cox regression
stcox i.infec_sev ,strata (setid)


*2. +IMD
stcox i.infec_sev i.imd, strata (setid)

*3. + lifestlye and co-morbidities
stcox i.infec_sev i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat, strata (setid)


*4. Fully adjusted cox regression (adjusted for IMD, alcohol consumption, smoking status, obesity status, CCI score cat, + ethnicity)
stcox i.infec_sev i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat i.ethnicity_4_na, strata (setid)


log close


*********************
*UTI
*********************
*/

log using "J:\EHR-Working\Charlotte_Sharon\Infections and SMI\logfiles\new models\UTI_severity_HES.smcl", replace

clear
pq use "Z:\GPRD_GOLD\Georgia\Georgia_extract_2024_12_13\processed_data\smi_cohort\20250902\UTI.parquet"

*restrict to HES linked data
keep if hes_linkage_eligible== 1
unique patid
*813337 with linked HES data


*drop matched sets without an exposure

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
*22241 obsverations dropped

*drop if no unexposed in matched set
drop if group_sum_unexposed==0
*370 observations dropped

unique patid
*N=791002




*create infection severity variable for analysis
gen infec_sev = .
replace infec_sev = 0 if severe_flag_composite == . & exposed == 0    // No infection
replace infec_sev = 1 if severe_flag_composite == 0 & exposed == 1    // Non-severe infection
replace infec_sev = 2 if severe_flag_composite == 1 & exposed == 1    // Severe infection

label define infec_sev_lbl 0 "No infection" 1 "Non-severe" 2 "Severe"
label values infec_sev infec_sev_lbl


*run cleaning dofile
run "J:\EHR-Working\Charlotte_Sharon\Infections and SMI\dofiles\1.cleaning.do" 

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
*2941

*set data for analysis
stset end_date, failure(smi_status=1) enter (indexdate) origin (indexdate) id(patid)scale(365.25)


*1. crude cox regression
stcox i.infec_sev ,strata (setid)


*2. +IMD
stcox i.infec_sev i.imd, strata (setid)

*3. + lifestlye and co-morbidities
stcox i.infec_sev i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat, strata (setid)


*4. Fully adjusted cox regression (adjusted for IMD, alcohol consumption, smoking status, obesity status, CCI score cat, + ethnicity)
stcox i.infec_sev i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat i.ethnicity_4_na, strata (setid)

log close

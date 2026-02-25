*secondary analysis - wanes over time
*SL Cadogan


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


*gen new censors for time periods
gen enddate_6mon=indexdate+182.5
gen enddate_1yr=indexdate+365
gen enddate_2yr=indexdate+(365*2)
gen enddate_3yr=indexdate+(365*3)
gen enddate_4yr=indexdate+(365*4)
gen enddate_5yr=indexdate+(365*5)


local censorvars enddate_6mon enddate_1yr enddate_2yr enddate_3yr enddate_4yr enddate_5yr

foreach cvar of local censorvars {
    // Create composite end date variable
    gen end_date_`cvar' = min(enddate, smi_incident_obsdate, infection_censor_obsdate, organic_psychosis_obsdate, `cvar')
    format end_date_`cvar' %d

    // Create event indicator variable
    gen smi_status_`cvar' = .
    replace smi_status_`cvar' = 1 if smi_incident_obsdate == end_date_`cvar'
    replace smi_status_`cvar' = 0 if smi_status_`cvar' == .

    // Set survival data
    stset end_date_`cvar', failure(smi_status_`cvar' == 1) origin(indexdate) id(patid) scale(365.25)

    // Run the fully adjusted Cox regression model and store results
    stcox i.exposed i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat i.ethnicity_4_na, strata(setid)
    estimates store adj_`cvar'
}


	
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


drop if organic_psychosis_obsdate <=indexdate
drop if infection_censor_obsdate<=indexdate
drop if smi_incident_obsdate<=indexdate
drop if enddate<=indexdate


*gen new censors for time periods
gen enddate_6mon=indexdate+182.5
gen enddate_1yr=indexdate+365
gen enddate_2yr=indexdate+(365*2)
gen enddate_3yr=indexdate+(365*3)
gen enddate_4yr=indexdate+(365*4)
gen enddate_5yr=indexdate+(365*5)


local censorvars enddate_6mon enddate_1yr enddate_2yr enddate_3yr enddate_4yr enddate_5yr

foreach cvar of local censorvars {
    // Create composite end date variable
    gen end_date_`cvar' = min(enddate, smi_incident_obsdate, infection_censor_obsdate, organic_psychosis_obsdate, `cvar')
    format end_date_`cvar' %d

    // Create event indicator variable
    gen smi_status_`cvar' = .
    replace smi_status_`cvar' = 1 if smi_incident_obsdate == end_date_`cvar'
    replace smi_status_`cvar' = 0 if smi_status_`cvar' == .

    // Set survival data
    stset end_date_`cvar', failure(smi_status_`cvar' == 1) origin(indexdate) id(patid) scale(365.25)

    // Run the fully adjusted Cox regression model and store results
    stcox i.exposed i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat i.ethnicity_4_na, strata(setid)
    estimates store adj_`cvar'
}

log close


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

drop if organic_psychosis_obsdate <=indexdate
drop if infection_censor_obsdate<=indexdate
drop if smi_incident_obsdate<=indexdate
drop if enddate<=indexdate


*gen new censors for time periods
gen enddate_6mon=indexdate+182.5
gen enddate_1yr=indexdate+365
gen enddate_2yr=indexdate+(365*2)
gen enddate_3yr=indexdate+(365*3)
gen enddate_4yr=indexdate+(365*4)
gen enddate_5yr=indexdate+(365*5)


local censorvars enddate_6mon enddate_1yr enddate_2yr enddate_3yr enddate_4yr enddate_5yr

foreach cvar of local censorvars {
    // Create composite end date variable
    gen end_date_`cvar' = min(enddate, smi_incident_obsdate, infection_censor_obsdate, organic_psychosis_obsdate, `cvar')
    format end_date_`cvar' %d

    // Create event indicator variable
    gen smi_status_`cvar' = .
    replace smi_status_`cvar' = 1 if smi_incident_obsdate == end_date_`cvar'
    replace smi_status_`cvar' = 0 if smi_status_`cvar' == .

    // Set survival data
    stset end_date_`cvar', failure(smi_status_`cvar' == 1) origin(indexdate) id(patid) scale(365.25)

    // Run the fully adjusted Cox regression model and store results
    stcox i.exposed i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat i.ethnicity_4_na, strata(setid)
    estimates store adj_`cvar'
}



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

drop if organic_psychosis_obsdate <=indexdate
drop if infection_censor_obsdate<=indexdate
drop if smi_incident_obsdate<=indexdate
drop if enddate<=indexdate


*gen new censors for time periods
gen enddate_6mon=indexdate+182.5
gen enddate_1yr=indexdate+365
gen enddate_2yr=indexdate+(365*2)
gen enddate_3yr=indexdate+(365*3)
gen enddate_4yr=indexdate+(365*4)
gen enddate_5yr=indexdate+(365*5)


local censorvars enddate_6mon enddate_1yr enddate_2yr enddate_3yr enddate_4yr enddate_5yr

foreach cvar of local censorvars {
    // Create composite end date variable
    gen end_date_`cvar' = min(enddate, smi_incident_obsdate, infection_censor_obsdate, organic_psychosis_obsdate, `cvar')
    format end_date_`cvar' %d

    // Create event indicator variable
    gen smi_status_`cvar' = .
    replace smi_status_`cvar' = 1 if smi_incident_obsdate == end_date_`cvar'
    replace smi_status_`cvar' = 0 if smi_status_`cvar' == .

    // Set survival data
    stset end_date_`cvar', failure(smi_status_`cvar' == 1) origin(indexdate) id(patid) scale(365.25)

    // Run the fully adjusted Cox regression model and store results
    stcox i.exposed i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat i.ethnicity_4_na, strata(setid)
    estimates store adj_`cvar'
}

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

drop if organic_psychosis_obsdate <=indexdate
drop if infection_censor_obsdate<=indexdate
drop if smi_incident_obsdate<=indexdate
drop if enddate<=indexdate


*gen new censors for time periods
gen enddate_6mon=indexdate+182.5
gen enddate_1yr=indexdate+365
gen enddate_2yr=indexdate+(365*2)
gen enddate_3yr=indexdate+(365*3)
gen enddate_4yr=indexdate+(365*4)
gen enddate_5yr=indexdate+(365*5)


local censorvars enddate_6mon enddate_1yr enddate_2yr enddate_3yr enddate_4yr enddate_5yr

foreach cvar of local censorvars {
    // Create composite end date variable
    gen end_date_`cvar' = min(enddate, smi_incident_obsdate, infection_censor_obsdate, organic_psychosis_obsdate, `cvar')
    format end_date_`cvar' %d

    // Create event indicator variable
    gen smi_status_`cvar' = .
    replace smi_status_`cvar' = 1 if smi_incident_obsdate == end_date_`cvar'
    replace smi_status_`cvar' = 0 if smi_status_`cvar' == .

    // Set survival data
    stset end_date_`cvar', failure(smi_status_`cvar' == 1) origin(indexdate) id(patid) scale(365.25)

    // Run the fully adjusted Cox regression model and store results
    stcox i.exposed i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat i.ethnicity_4_na, strata(setid)
    estimates store adj_`cvar'
}

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

drop if organic_psychosis_obsdate <=indexdate
drop if infection_censor_obsdate<=indexdate
drop if smi_incident_obsdate<=indexdate
drop if enddate<=indexdate


*gen new censors for time periods
gen enddate_6mon=indexdate+182.5
gen enddate_1yr=indexdate+365
gen enddate_2yr=indexdate+(365*2)
gen enddate_3yr=indexdate+(365*3)
gen enddate_4yr=indexdate+(365*4)
gen enddate_5yr=indexdate+(365*5)


local censorvars enddate_6mon enddate_1yr enddate_2yr enddate_3yr enddate_4yr enddate_5yr

foreach cvar of local censorvars {
    // Create composite end date variable
    gen end_date_`cvar' = min(enddate, smi_incident_obsdate, infection_censor_obsdate, organic_psychosis_obsdate, `cvar')
    format end_date_`cvar' %d

    // Create event indicator variable
    gen smi_status_`cvar' = .
    replace smi_status_`cvar' = 1 if smi_incident_obsdate == end_date_`cvar'
    replace smi_status_`cvar' = 0 if smi_status_`cvar' == .

    // Set survival data
    stset end_date_`cvar', failure(smi_status_`cvar' == 1) origin(indexdate) id(patid) scale(365.25)

    // Run the fully adjusted Cox regression model and store results
    stcox i.exposed i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat i.ethnicity_4_na, strata(setid)
    estimates store adj_`cvar'
}

log close
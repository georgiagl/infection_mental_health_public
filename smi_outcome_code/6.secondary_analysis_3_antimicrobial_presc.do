*Anti-microbial resistance
*variable


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

tab antimicro_presc
* gen new exposure variables for analysis

*exposure with prescription
gen exposed_antimicro_presc_7days=1 if antimicro_presc==1 &exposed==1
replace exposed_antimicro_presc_7days=0 if exposed_antimicro_presc_7days==.
tab exposed_antimicro_presc_7days

*exposure without prescription
gen exposed_antimicro_no_presc_7days=1 if antimicro_presc==0 &exposed==1
replace exposed_antimicro_no_presc_7days=0 if exposed_antimicro_no_presc_7days==.
tab exposed_antimicro_no_presc_7days



drop if organic_psychosis_obsdate <=indexdate
drop if infection_censor_obsdate<=indexdate
drop if smi_incident_obsdate<=indexdate

*Set the data for Cox regression: censoring variables (for SMI): organic_psychosis_obsdate,  infection censor: infection_censor_obsdate  
gen end_date=min(enddate, smi_incident_obsdate, infection_censor_obsdate, organic_psychosis_obsdate)
format end_date %d  
 
*generate smi_status variable
gen smi_status=1 if smi_incident_obsdate==end_date
tab smi_status


stset end_date, failure(smi_status=1) enter (indexdate) origin (indexdate) id(patid)scale(365.25)


*4. Fully adjusted cox regression (adjusted for IMD, alcohol consumption, smoking status, obesity status, CCI score cat, + ethnicity)
stcox i.exposed_antimicro_presc_7days i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat i.ethnicity_4_na

*EXPOSURE: with infection, but without antimicro_presc within 7 days 
stcox i.exposed_antimicro_no_presc_7days i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat i.ethnicity_4_na


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

tab antimicro_presc
* gen new exposure variables for analysis

*exposure with prescription
gen exposed_antimicro_presc_7days=1 if antimicro_presc==1 &exposed==1
replace exposed_antimicro_presc_7days=0 if exposed_antimicro_presc_7days==.
tab exposed_antimicro_presc_7days

*exposure without prescription
gen exposed_antimicro_no_presc_7days=1 if antimicro_presc==0 &exposed==1
replace exposed_antimicro_no_presc_7days=0 if exposed_antimicro_no_presc_7days==.
tab exposed_antimicro_no_presc_7days



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


*4. Fully adjusted cox regression (adjusted for IMD, alcohol consumption, smoking status, obesity status, CCI score cat, + ethnicity)
stcox i.exposed_antimicro_presc_7days i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat i.ethnicity_4_na

*EXPOSURE: with infection, but without antimicro_presc within 7 days 
stcox i.exposed_antimicro_no_presc_7days i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat i.ethnicity_4_na



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

tab antimicro_presc
* gen new exposure variables for analysis

*exposure with prescription
gen exposed_antimicro_presc_7days=1 if antimicro_presc==1 &exposed==1
replace exposed_antimicro_presc_7days=0 if exposed_antimicro_presc_7days==.
tab exposed_antimicro_presc_7days

*exposure without prescription
gen exposed_antimicro_no_presc_7days=1 if antimicro_presc==0 &exposed==1
replace exposed_antimicro_no_presc_7days=0 if exposed_antimicro_no_presc_7days==.
tab exposed_antimicro_no_presc_7days



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


*4. Fully adjusted cox regression (adjusted for IMD, alcohol consumption, smoking status, obesity status, CCI score cat, + ethnicity)
stcox i.exposed_antimicro_presc_7days i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat i.ethnicity_4_na

*EXPOSURE: with infection, but without antimicro_presc within 7 days 
stcox i.exposed_antimicro_no_presc_7days i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat i.ethnicity_4_na




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

tab antimicro_presc
* gen new exposure variables for analysis

*exposure with prescription
gen exposed_antimicro_presc_7days=1 if antimicro_presc==1 &exposed==1
replace exposed_antimicro_presc_7days=0 if exposed_antimicro_presc_7days==.
tab exposed_antimicro_presc_7days

*exposure without prescription
gen exposed_antimicro_no_presc_7days=1 if antimicro_presc==0 &exposed==1
replace exposed_antimicro_no_presc_7days=0 if exposed_antimicro_no_presc_7days==.
tab exposed_antimicro_no_presc_7days



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


*4. Fully adjusted cox regression (adjusted for IMD, alcohol consumption, smoking status, obesity status, CCI score cat, + ethnicity)
stcox i.exposed_antimicro_presc_7days i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat i.ethnicity_4_na

*EXPOSURE: with infection, but without antimicro_presc within 7 days 
stcox i.exposed_antimicro_no_presc_7days i.imd i.alcohol_abuse i.smoking_status_simp i.obese i.cci_cat i.ethnicity_4_na

log close
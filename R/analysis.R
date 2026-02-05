
# set up for main CMD analyses- with all censoring vars and "obsdate" used- this is "cmd_incident_obsdate" renamed in data cleaning
create_cohort_post_exclusion <- function(cohort_wide) {

  temp <- cohort_wide 
  
  #Set enddate and remove patients without follow-up time: have added CMD censoring vars
  # NB obsdate is CMD event date 
  # added infection_censor_obsdate which censors unexposed when they get an infection (these are different to the exposed who start off unexposed)
  temp <- temp |>
    ftransform(enddate=pmin(enddate, obsdate, infection_censor_obsdate, organic_cmd_obsdate, organic_dep_obsdate, dementia_obsdate, smi_incident_obsdate, smi_historic_obsdate, ocd_obsdate , ptsd_obsdate, cmd_historic_obsdate, reg_date_of_death_ons,  na.rm = TRUE)) |> 
    ftransform(event=if_else(obsdate==enddate, 1, 0, missing=0)) |> 
    fsubset(enddate>indexdate)
  
  #Remove matched sets that have become incomplete (i.e. no longer containing at least 1 exposed and 1 unexposed person)
  complete_sets <- temp |> 
    fgroup_by(setid) |> 
    get_vars("exposed") |> 
    fmean() |> 
    fsubset(exposed!=0 & exposed !=1) |> 
    pull(setid)
  
  cohort_post_exclusion <- temp |> 
    fsubset(setid %in% complete_sets) |> 
    ftransform(indexdate_copy=indexdate) #Make copy because variable will be deleted when used in SurvSplit
}

create_cohort_post_exclusion_smi <- function(cohort_wide) {
  
  temp <- cohort_wide 
  
  #Set enddate and remove patients without follow-up time: have added SMI censoring vars
  # NB obsdate is SMI event date 
  # added infection_censor_obsdate which censors unexposed when they get an infection (these are different to the exposed who start off unexposed)
  temp <- temp |>
    ftransform(enddate=pmin(enddate, obsdate, infection_censor_obsdate, organic_psychosis_obsdate, na.rm = TRUE)) |> 
    ftransform(event=if_else(obsdate==enddate, 1, 0, missing=0)) |> 
    fsubset(enddate>indexdate)
  
  #Remove matched sets that have become incomplete (i.e. no longer containing at least 1 exposed and 1 unexposed person)
  complete_sets <- temp |> 
    fgroup_by(setid) |> 
    get_vars("exposed") |> 
    fmean() |> 
    fsubset(exposed!=0 & exposed !=1) |> 
    pull(setid)
  
  cohort_post_exclusion <- temp |> 
    fsubset(setid %in% complete_sets) |> 
    ftransform(indexdate_copy=indexdate) #Make copy because variable will be deleted when used in SurvSplit
}



analysis_hrs <- function(outcome, cohort_post_exclusion, cohort_labels, exposure, model) {
  
  #Create survival object
  survival_object <- Surv(time = as.numeric(cohort_post_exclusion$indexdate),
                          time2 = as.numeric(cohort_post_exclusion$enddate),
                          origin = as.numeric(cohort_post_exclusion$indexdate),
                          event = cohort_post_exclusion$event)  
  
  survival::coxph(formula(paste("survival_object ~", exposure, model, " + strata(setid) + cluster(patid)")), 
                  data=cohort_post_exclusion, id=patid) |>
    broom::tidy(exponentiate=TRUE, conf.int=TRUE) |> 
    dplyr::mutate(cohort=cohort_labels,
                  exposure=exposure,
                  outcome=outcome,
                  model=names(model)
    )
}

analysis_rates <- function(outcome, cohort_post_exclusion, cohort_labels, model, exposure) {
  
  #Exclude rows with missing values for any variable
  vars <- str_split_1(model, " \\+ | \\* ") |> str_subset(".+") 
  vars <- vars[vars %in% names(cohort_post_exclusion)]
  temp <- cohort_post_exclusion |> filter(if_all(vars, \(x) !is.na(x)))
  
      survival_object <- Surv(time = as.numeric(temp$indexdate),
                            time2 = as.numeric(temp$enddate),
                            origin = as.numeric(temp$indexdate),
                            event = temp$event)

  fit <- survfit(formula(paste("survival_object ~", exposure)), data = temp)
  
  
  # Get person years
  results_rates <- pyears(formula(paste("survival_object ~", exposure)), 
                          data = temp, 
                          scale = 365.25) |>
    tidy() |> 
    mutate(rate=(event/pyears)*1000) |> 
    bind_cols(tibble(term=paste0(exposure, levels(as.factor(temp[[exposure]]))), 
                     cohort=cohort_labels,
                     outcome=outcome,
                     model=names(model)))  |> 
    # gives quintiles of FU time: 
  bind_cols(as.data.frame(quantile(fit), row.names=NULL))
  
  #Get unexposed person years in the same row
  results_rates |> 
    group_by(cohort, outcome, model) |>
    mutate(pyears_unexposed=pyears[1],
           n_unexposed=n[1],
           event_unexposed=event[1],
           rate_unexposed=rate[1],
           ratio=1/(pyears[2]/pyears[1]),
           ratio_text=paste0("1:", round(ratio))) |>
    ungroup() 
  
}


## FUNCTIONS FOR SENSITIVITY ANALS: 
# create cohort with no censoring at alternative MH diagnosis: 

create_cohort_post_exclusion_no_cmd_censor <- function(cohort_wide) {
  
  temp <- cohort_wide 
  
  #Set enddate and remove patients without follow-up time:  note: have removed CMD censoring vars
  # NB obsdate is CMD event date 
  temp <- temp |>
    ftransform(enddate=pmin(enddate, obsdate, infection_censor_obsdate, dementia_obsdate, cmd_historic_obsdate, reg_date_of_death_ons,  na.rm = TRUE)) |> 
    ftransform(event=if_else(obsdate==enddate, 1, 0, missing=0)) |> 
    fsubset(enddate>indexdate)
  
  #Remove matched sets that have become incomplete (i.e. no longer containing at least 1 exposed and 1 unexposed person)
  complete_sets <- temp |> 
    fgroup_by(setid) |> 
    get_vars("exposed") |> 
    fmean() |> 
    fsubset(exposed!=0 & exposed !=1) |> 
    pull(setid)
  
  cohort_post_exclusion <- temp |> 
    fsubset(setid %in% complete_sets) |> 
    ftransform(indexdate_copy=indexdate) #Make copy because variable will be deleted when used in SurvSplit
}

# create cohort with 1st March 2020: adding to the enddate var: 
create_cohort_post_exclusion_march20_enddate <- function(cohort_wide) {
  
  temp <- cohort_wide 
  
  # NB obsdate is CMD event date 
  temp <- temp |>
    ftransform(enddate=pmin( as.Date("2020-03-01"), enddate, obsdate, infection_censor_obsdate, organic_cmd_obsdate, organic_dep_obsdate, dementia_obsdate, smi_incident_obsdate, smi_historic_obsdate, ocd_obsdate , ptsd_obsdate, cmd_historic_obsdate, reg_date_of_death_ons,  na.rm = TRUE)) |> 
    ftransform(event=if_else(obsdate==enddate, 1, 0, missing=0)) |> 
    fsubset(enddate>indexdate)
  
  #Remove matched sets that have become incomplete (i.e. no longer containing at least 1 exposed and 1 unexposed person)
  complete_sets <- temp |> 
    fgroup_by(setid) |> 
    get_vars("exposed") |> 
    fmean() |> 
    fsubset(exposed!=0 & exposed !=1) |> 
    pull(setid)
  
  cohort_post_exclusion <- temp |> 
    fsubset(setid %in% complete_sets) |> 
    ftransform(indexdate_copy=indexdate) #Make copy because variable will be deleted when used in SurvSplit
}

# create cohort Incident ANXIETY as outcome: adding to the enddate var: 
create_cohort_post_exclusion_anxiety <- function(cohort_wide) {
  
  temp <- cohort_wide 
  
  # NB set obsdate to be anxiety_incident_obsdate
  temp <- temp |>
    mutate(obsdate = anxiety_incident_obsdate) |> 
    ftransform(enddate=pmin( as.Date("2020-03-01"), enddate, obsdate, infection_censor_obsdate, organic_cmd_obsdate, smi_incident_obsdate, smi_historic_obsdate, ocd_obsdate , ptsd_obsdate, cmd_historic_obsdate, reg_date_of_death_ons,  na.rm = TRUE)) |> 
    ftransform(event=if_else(obsdate==enddate, 1, 0, missing=0)) |> 
    fsubset(enddate>indexdate)
  
  #Remove matched sets that have become incomplete (i.e. no longer containing at least 1 exposed and 1 unexposed person)
  complete_sets <- temp |> 
    fgroup_by(setid) |> 
    get_vars("exposed") |> 
    fmean() |> 
    fsubset(exposed!=0 & exposed !=1) |> 
    pull(setid)
  
  cohort_post_exclusion <- temp |> 
    fsubset(setid %in% complete_sets) |> 
    ftransform(indexdate_copy=indexdate) #Make copy because variable will be deleted when used in SurvSplit
}

# create cohort Incident DEPRESSION as outcome: adding to the enddate var: 
create_cohort_post_exclusion_depression <- function(cohort_wide) {
  
  temp <- cohort_wide 
  
  # NB set obsdate to be depression_incident_obsdate
  temp <- temp |>
    mutate(obsdate = depression_incident_obsdate) |> 
    ftransform(enddate=pmin( as.Date("2020-03-01"), enddate, obsdate, infection_censor_obsdate, organic_dep_obsdate, smi_incident_obsdate, smi_historic_obsdate,  cmd_historic_obsdate, reg_date_of_death_ons,  na.rm = TRUE)) |> 
    ftransform(event=if_else(obsdate==enddate, 1, 0, missing=0)) |> 
    fsubset(enddate>indexdate)
  
  #Remove matched sets that have become incomplete (i.e. no longer containing at least 1 exposed and 1 unexposed person)
  complete_sets <- temp |> 
    fgroup_by(setid) |> 
    get_vars("exposed") |> 
    fmean() |> 
    fsubset(exposed!=0 & exposed !=1) |> 
    pull(setid)
  
  cohort_post_exclusion <- temp |> 
    fsubset(setid %in% complete_sets) |> 
    ftransform(indexdate_copy=indexdate) #Make copy because variable will be deleted when used in SurvSplit
}



# INVESTIGATING ASSOCIATION PERSISTANCE OR WANING: 
# SPLITTING IN TO THE FU PERIODS SPECIFIED IN THE PROTOCOL
# 0-6 months: 
create_cohort_post_exclusion_6mon <- function(cohort_wide) {
  
  temp <- cohort_wide 
  
  # NB obsdate is CMD event date 
  temp <- temp |>
    ftransform(enddate=pmin(enddate_6mon, obsdate, infection_censor_obsdate, organic_cmd_obsdate, organic_dep_obsdate, dementia_obsdate, smi_incident_obsdate, smi_historic_obsdate, ocd_obsdate , ptsd_obsdate, cmd_historic_obsdate, reg_date_of_death_ons,  na.rm = TRUE)) |> 
    ftransform(event=if_else(obsdate==enddate, 1, 0, missing=0)) |> 
    fsubset(enddate>indexdate)
  
  #Remove matched sets that have become incomplete (i.e. no longer containing at least 1 exposed and 1 unexposed person)
  complete_sets <- temp |> 
    fgroup_by(setid) |> 
    get_vars("exposed") |> 
    fmean() |> 
    fsubset(exposed!=0 & exposed !=1) |> 
    pull(setid)
  
  cohort_post_exclusion <- temp |> 
    fsubset(setid %in% complete_sets) |> 
    ftransform(indexdate_copy=indexdate) #Make copy because variable will be deleted when used in SurvSplit
}

# 0-1 year: 
create_cohort_post_exclusion_1yr <- function(cohort_wide) {
  
  temp <- cohort_wide 
  
  # NB obsdate is CMD event date 
  temp <- temp |>
    ftransform(enddate=pmin(enddate_1yr, obsdate, infection_censor_obsdate, organic_cmd_obsdate, organic_dep_obsdate, dementia_obsdate, smi_incident_obsdate, smi_historic_obsdate, ocd_obsdate , ptsd_obsdate, cmd_historic_obsdate, reg_date_of_death_ons,  na.rm = TRUE)) |> 
    ftransform(event=if_else(obsdate==enddate, 1, 0, missing=0)) |> 
    fsubset(enddate>indexdate)
  
  #Remove matched sets that have become incomplete (i.e. no longer containing at least 1 exposed and 1 unexposed person)
  complete_sets <- temp |> 
    fgroup_by(setid) |> 
    get_vars("exposed") |> 
    fmean() |> 
    fsubset(exposed!=0 & exposed !=1) |> 
    pull(setid)
  
  cohort_post_exclusion <- temp |> 
    fsubset(setid %in% complete_sets) |> 
    ftransform(indexdate_copy=indexdate) #Make copy because variable will be deleted when used in SurvSplit
}

# 0-2 year: 
create_cohort_post_exclusion_2yr <- function(cohort_wide) {
  
  temp <- cohort_wide 
  
  # NB obsdate is CMD event date 
  temp <- temp |>
    ftransform(enddate=pmin(enddate_2yr, obsdate, infection_censor_obsdate, organic_cmd_obsdate, organic_dep_obsdate, dementia_obsdate, smi_incident_obsdate, smi_historic_obsdate, ocd_obsdate , ptsd_obsdate, cmd_historic_obsdate, reg_date_of_death_ons,  na.rm = TRUE)) |> 
    ftransform(event=if_else(obsdate==enddate, 1, 0, missing=0)) |> 
    fsubset(enddate>indexdate)
  
  #Remove matched sets that have become incomplete (i.e. no longer containing at least 1 exposed and 1 unexposed person)
  complete_sets <- temp |> 
    fgroup_by(setid) |> 
    get_vars("exposed") |> 
    fmean() |> 
    fsubset(exposed!=0 & exposed !=1) |> 
    pull(setid)
  
  cohort_post_exclusion <- temp |> 
    fsubset(setid %in% complete_sets) |> 
    ftransform(indexdate_copy=indexdate) #Make copy because variable will be deleted when used in SurvSplit
}

# 0-3 year: 
create_cohort_post_exclusion_3yr <- function(cohort_wide) {
  
  temp <- cohort_wide 
  
  # NB obsdate is CMD event date 
  temp <- temp |>
    ftransform(enddate=pmin(enddate_3yr, obsdate, infection_censor_obsdate, organic_cmd_obsdate, organic_dep_obsdate, dementia_obsdate, smi_incident_obsdate, smi_historic_obsdate, ocd_obsdate , ptsd_obsdate, cmd_historic_obsdate, reg_date_of_death_ons,  na.rm = TRUE)) |> 
    ftransform(event=if_else(obsdate==enddate, 1, 0, missing=0)) |> 
    fsubset(enddate>indexdate)
  
  #Remove matched sets that have become incomplete (i.e. no longer containing at least 1 exposed and 1 unexposed person)
  complete_sets <- temp |> 
    fgroup_by(setid) |> 
    get_vars("exposed") |> 
    fmean() |> 
    fsubset(exposed!=0 & exposed !=1) |> 
    pull(setid)
  
  cohort_post_exclusion <- temp |> 
    fsubset(setid %in% complete_sets) |> 
    ftransform(indexdate_copy=indexdate) #Make copy because variable will be deleted when used in SurvSplit
}

# 0-4 year: 
create_cohort_post_exclusion_4yr <- function(cohort_wide) {
  
  temp <- cohort_wide 
  
  # NB obsdate is CMD event date 
  temp <- temp |>
    ftransform(enddate=pmin(enddate_4yr, obsdate, infection_censor_obsdate, organic_cmd_obsdate, organic_dep_obsdate, dementia_obsdate, smi_incident_obsdate, smi_historic_obsdate, ocd_obsdate , ptsd_obsdate, cmd_historic_obsdate, reg_date_of_death_ons,  na.rm = TRUE)) |> 
    ftransform(event=if_else(obsdate==enddate, 1, 0, missing=0)) |> 
    fsubset(enddate>indexdate)
  
  #Remove matched sets that have become incomplete (i.e. no longer containing at least 1 exposed and 1 unexposed person)
  complete_sets <- temp |> 
    fgroup_by(setid) |> 
    get_vars("exposed") |> 
    fmean() |> 
    fsubset(exposed!=0 & exposed !=1) |> 
    pull(setid)
  
  cohort_post_exclusion <- temp |> 
    fsubset(setid %in% complete_sets) |> 
    ftransform(indexdate_copy=indexdate) #Make copy because variable will be deleted when used in SurvSplit
}

# 0-5 year: 
create_cohort_post_exclusion_5yr <- function(cohort_wide) {
  
  temp <- cohort_wide 
  
  # NB obsdate is CMD event date 
  temp <- temp |>
    ftransform(enddate=pmin(enddate_5yr, obsdate, infection_censor_obsdate, organic_cmd_obsdate, organic_dep_obsdate, dementia_obsdate, smi_incident_obsdate, smi_historic_obsdate, ocd_obsdate , ptsd_obsdate, cmd_historic_obsdate, reg_date_of_death_ons,  na.rm = TRUE)) |> 
    ftransform(event=if_else(obsdate==enddate, 1, 0, missing=0)) |> 
    fsubset(enddate>indexdate)
  
  #Remove matched sets that have become incomplete (i.e. no longer containing at least 1 exposed and 1 unexposed person)
  complete_sets <- temp |> 
    fgroup_by(setid) |> 
    get_vars("exposed") |> 
    fmean() |> 
    fsubset(exposed!=0 & exposed !=1) |> 
    pull(setid)
  
  cohort_post_exclusion <- temp |> 
    fsubset(setid %in% complete_sets) |> 
    ftransform(indexdate_copy=indexdate) #Make copy because variable will be deleted when used in SurvSplit
}

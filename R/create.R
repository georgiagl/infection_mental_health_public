
# Denominator data: all eligible patients: read in txt file, open as paruet, do some mx, then collect: NB. schema is specified by order so make sure it matches the file order

denom_data <- function(file) {
  open_tsv_dataset((file), 
                   schema = schema(patid=string(), pracid=string(), gender=int64(), yob=string(), mob=string(), emis_ddate = string(), regstartdate = string() , patienttypeid = string(), regenddate = string(), acceptable=string(), cprd_ddate = string (), uts=string(), lcd=string(), region=int64()), 
                   skip = 1) %>% 
  
    # fix data formatting
    mutate(regstartdate = as.Date(regstartdate, format = "%d/%m/%Y"), 
           regenddate = as.Date(regenddate, format = "%d/%m/%Y"), 
           cprd_ddate = as.Date(cprd_ddate, format = "%d/%m/%Y"), 
           lcd = as.Date(lcd, format = "%d/%m/%Y"),
           yob = as.integer(yob), 
           patid = as.character(patid), 
           pracid = as.integer(pracid)) %>%
    
    # remove variables dont need: 
    select(-c(emis_ddate, uts, patienttypeid, mob, acceptable)) %>% 
    
    # collecting here as cant gen startdate* or enddate vars within arrow
    collect () 
}


# Denom data: do some more data management: generate study start and end dates: 

denom_data_mx <- function(data, pracid_exc) {   
  data %>% 

    # some missing cprd_ddates coded as "1860-01-01" - recode to missing (data.table::fifdate preserves date format, unlike ifelse) 
    mutate(cprd_ddate = fifelse(cprd_ddate == "1860-01-01", NA_Date_, cprd_ddate)) %>% 
    
    #dob: no month or date of birth in CPRD- assume middle of year (1st July)
    ftransform(dob = as.Date(paste0(yob, "-07-01"))) %>% 
    
    # date person turns 18: 
    ftransform(eighteen_date = as.Date(paste0(yob + 18, "-07-01"))) %>% 
    
    # Startdate: latest of  (regstartdate + 365 / eighteen_date /study start(01-01-2007))
    ftransform(startdate=pmax(regstartdate + 365, eighteen_date, as.Date("2007-01-01"))) %>%   

    # Enddate: earliest of: end of registration, date of death, date of the last data collection from the practice, or the end of the study (2024-06-15).
    ftransform(enddate=pmin(regenddate, cprd_ddate, lcd, as.Date("2024-06-15"), na.rm = TRUE)) %>% 
    
    # Remove 46 duplicated practids - as per CPRD Aurum Data Spec, v3.5, Date: 30 Aug 2024
    filter(!grepl(paste(pracid_exc, collapse='|'), pracid))

} 


## infection cases: : using the file restricted on study criteria ("*_Define_results.txt") rather than the full "Observation_001" files
inf_data <-  function(file) {
  
                open_tsv_dataset((file), 
                               schema = schema(patid=string(), indexdate=string()),
                               skip = 1) %>% 
            
                select(patid, indexdate)  %>% 
                
                distinct(patid, indexdate) %>% 
                collect() %>% 
                
                # format indexdate and generate an infection var to keep track of branches - takes name from file name
                mutate(indexdate = as.Date(indexdate, format = "%d/%m/%Y"), 
                       infection = str_extract(file, pattern = "(?<=Georgia_)(.*?)(?=_def)"))
  
} 

# Merge pracid from denominator in to infection data and drop from the infection cases the 46 duplicated practid: 

inf_pracid_exc <- function(inf_data_collect, denom_data_collect, pracid_exc) {
  
    # join in pracid
    left_join(inf_data_collect, denom_data_collect, by= "patid") %>% 
    select(c(patid, indexdate, pracid, infection)) %>% 
    
    # Remove 46 duplicated practids - as per CPRD Aurum Data Spec, v3.5, Date: 30 Aug 2024
    filter(!grepl(paste(pracid_exc, collapse='|'), pracid))
  
}

# Conduct 10% random sample of individuals, grouped (stratified) by practid (for UTI, SSTI, LRTI, and Gastro), keeping full sample for sepsis and mening/enceph
sample <- function(data) {
  
  # random sample for all 4 infections other than sepsis and mening/enceph
  if (!grepl("sepsis|mening_encep", data[1,]$infection) ) {
    
    slice_sample(data, prop = 0.1, by = "pracid")
    
  } else {
    data 
  }
}
 
# join each of the infection datasets to the denominator datasets
#  drop vars dont need any more: 

denom_inf <- function(denom, inf, pracid_exc) {
              left_join(denom, inf, by=c("patid")) %>%
               
              # tidy up the 2 pracid variables- drop one 
              rename(pracid = pracid.x) %>% 
    
              # drop denom vars dont need to make faster- 
              select(-c(regstartdate, regenddate, cprd_ddate, lcd, region, dob, eighteen_date, pracid.y)) %>% 
              arrange(patid) 
}

# linkage eligibility file: 
linkage_data <- function(file) {
  open_tsv_dataset((file), 
                   schema = schema(patid=string(), pracid=string(), linkyear=string(), 
                                   lsoa_e=int64(), sgss_e=int64(), chess_e=int64(),
                                   hes_ae_e=int64(), hes_did_e=int64(), cr_e=int64(),
                                   sact_e=int64(), rtds_e=int64(), hes_apc_e=int64(),
                                   hes_op_e=int64(), ons_death_e=int64()), 
                   skip = 1) %>% 
    
    # collecting here as cant gen startdate* or enddate vars within arrow
    collect () 
}


## COHORT matching:  CODE FROM JM adapted from Krishnan STATA code
# variables required in memory before running (correctly named):
# patid:          CPRD patient id [nb if using Aurum data, patid MUST be stored as a "double" precision variable]
# indexdate:      date of "exposure" for exposed patients (missing for potential controls)
# gender:         gender, numerically coded (e.g. 1=male, 2=female)
# startdate:      date of start of CPRD follow-up
# enddate:        date of end of follow-up as a potential control, generally = end of CPRD follow-up, but see "important note" below
# exposed:        indicator: 1 for exposed patients, 0 for potential controls
# yob:            year of birth
# IMPORTANT NOTE: in most cases it is desirable to allow exposed patients
# to be available as controls prior to their date of first exposure. Such
# patients should be included in the dataset twice (i.e. two separate rows):
# once as exposed (exposed = 1, startdate = start of CPRD follow-up,
# indexdate = date of first exposure, enddate = end of CPRD follow-up),
# and once as a potential control (exposed = 0, startdate = start of CPRD
# follow-up, indexdate = missing, enddate = date of first exposure-1)
# create_cohort_matchable <- function(cohort_matchable) {
#
#   # Get exposed patients to be available as controls prior to their date of first exposure
  create_cohort_matchable <- function(cohort_eligible) {

    # store all the data - exposed and unexposed:
    cohort_matchable <- cohort_eligible  %>%
          ftransform(exposed=if_else(is.na(indexdate), 0, 1))

    # Get exposed patients to be available as controls prior to their date of first exposure
    pre_exposure_ppl <- cohort_matchable %>%

        fsubset(exposed==1)  %>%
        mutate(enddate=indexdate-1,
               indexdate=NA_Date_,
               exposed=0) %>%
        fsubset(enddate > startdate)

      bind_rows(cohort_matchable, pre_exposure_ppl)
  }


  #' Create matched cohort using sequential trials matching
  #' @description each daily trial includes all n eligible people who
  #' become exposed on that day (exposed=1) and
  #' a sample of n eligible controls (exposed=0) who:
  #' - had not been exposed on or before that day (still at risk of becoming exposed);
  #' - still at risk of an outcome (not left the study);
  #' - had not already been selected as a control in a previous trial
  #' @param grouped_cohort_to_match A data frame with one or two rows per participant with:
  #'  1. $patid: The patient ID
  #'  2. $startdate: the date people become eligible
  #'  3. $indexdate: the date eligible people got exposed
  #'  4. $enddate: the date people leave the study
  #'  @param dayspriorreg days prior registration required for controls
  #'  GGL NB: for this study is 1 year but this is managed in the denom_data_mx function already in the generation of the studystart var
  #'
  create_cohort_matched <- function(cohort_matchable_grouped, dayspriorreg = 0) {

    #Map across every practice group
    library(future)
    plan(multisession, workers = 8)
    furrr::future_map_dfr(cohort_matchable_grouped, .progress = TRUE, \(x) {

      #Sort by indexdate, and then randomly
      cohort_matchable <- x  %>%
        ftransform(sortunique=runif(nrow(x))) %>%
        roworder(exposed, indexdate, sortunique) %>%
        fselect(-sortunique)

      exposed <- cohort_matchable %>% fsubset(exposed==1) %>% ftransform(setid=patid)
      unexposed <- cohort_matchable %>% fsubset(exposed==0)

      matched <- exposed[0,] #Make empty dataframe with same columns to be filled

      #Loop through all people that ever get exposed (each one gets matched to people who are unexposed at the same time)
      for (i in 1:nrow(exposed)) {

        exposed_pat <- exposed[i,]
        matchday <- exposed_pat$indexdate

        #Drop people that can't be matched anymore (either because they have already been matched or they have passed the study end date)
        unexposed <- unexposed %>%
          fsubset(enddate > matchday & !(patid %in% matched$patid))

        #Perform matching
        new <- unexposed %>%
          fsubset(gender==exposed_pat$gender) %>%
          fsubset((startdate+dayspriorreg) <= matchday) %>%
          ftransform(age_difference=abs(yob-exposed_pat$yob)) %>%
          fsubset(age_difference<=2) %>%
          roworder(age_difference) %>% # the closest matches are given priority
          slice(1:5) %>%
          ftransform(setid=exposed_pat$patid,
                     indexdate=as_date(matchday)) %>%  #Set the indexdate for everyone to the day the exposed individual got exposed
          fselect(-age_difference)
        if (nrow(new)>0) matched <- bind_rows(matched, exposed_pat, new)
        # cli::cli_progress_update()
      }
      matched %>%
        roworder(setid, -exposed) %>%
        select(patid, exposed, indexdate, enddate, setid)})
  }


  # Get some summary information on # of matches and if anyone unmatched:

    cohort_stats <- function(data) {
      
      setDT(data)
      
      matches <- data[, .(.N), by = .(setid)]
      table <-  table(matches$N)

    }
    
    # Save patient lists as patid text files - as needed for Extract
    
    write_delim_compressed_in_chunks_and_return_paths <- function(x, max) {
      
      size <- seq_len(nrow(x))
      chunked <- split(x, ceiling(size/max))
      
      paths <- paste0(extract_data, deparse(substitute(x)), seq_len(length(chunked)), ".txt")
      
      for (i in seq_len(length(chunked))) {
        write_delim(chunked[[i]], file=paths[[i]])
      }
      return(paths)
    }
    
    # Save patient lists as patid text files - as needed for linked data - 
    
    write_delim_compressed_in_chunks_and_return_paths_linked <- function(x, max) {
      
      size <- seq_len(nrow(x))
      chunked <- split(x, ceiling(size/max))
      

      for (i in seq_len(length(chunked))) {
        write_delim(chunked[[i]], file=paths[[i]])
      }
      return(paths)
    }
    

     # Calculate participant counts at various stages

     table_infection_n <- function(inf_data_collect) {
 
       inf_data_collect %>%
         filter(!is.na(indexdate)) %>%
         count() %>%
         mutate(step="had infection in study period")

     }
    
     # Table of infection sample size once random sample code run: 
     table_inf_sample_n <- function(random_sample) {
       random_sample %>% 
         count() %>% 
         mutate(step = "inf_n_rand_10sample")
     }
     
     
      table_cohort_matched_pop <- function(cohort_matched) {
        cohort_matched %>% 
        count() %>%
        mutate(step="matched cohort (1:5)")
      } 
      
     table_cohort_matched_exp_unexpo <- function(cohort_matched) {
       cohort_matched  %>% 
         group_by(exposed) %>% 
         count() %>% 
         mutate(step ="matched expo/unexpo")
     }
     
     
     # save matched cohort (with setid, patid, expos/unexpo status, indexdate and enddate) as parquet (to be read back in and merged with extract): 
     
     write_parquet_and_return_path <- function(x, path, infection) {
       path <- paste0(path, infection, ".parquet")
       write_parquet(x, path)
       return(path)
     }
     
     
    # read back in the matched cohort parquet files as saved in previous function: 
     
     inf_parquet_open  <- function(file) {
       open_dataset(file) %>% 
         
         # fix data formatting
         mutate(indexdate = as.Date(indexdate, format = "%d/%m/%Y"), 
                enddate = as.Date(enddate, format = "%d/%m/%Y"), 
                patid = as.character(patid), 
                setid = as.character(setid)) %>%
         
         collect () 
     }

     
     
     # merge sex and YOB back on to the matched cohort data - basis for Table 1 
     
     join_sex_age <- function(inf, denom) {
       left_join(x = inf, y = denom, by = "patid") %>% 
         # gen age at index date: 
         mutate(indexdate_year = as.integer(format(as.Date(indexdate), "%Y")), 
                indexdate_age = indexdate_year - yob, 
                # catgorise age as per catgories used in Adesanya et al 
                indexdate_age_cat = cut(indexdate_age, breaks = c(-Inf, 29, 39, 59, Inf), labels = c("18-29", "30-39", "40-59", "60+")), 
          # data management:       
                gender = factor(gender,levels = c(1,2,3),
                                labels = c("male", "female", "indeterminate"))) %>% 
         select(c(indexdate_age, indexdate_age_cat, gender))
     }


     ## Matching firstevent to patid
     # @return A of dataframes each with one row per patient, with patient ID and a TRUE/FALSE variable if they had the event in the variable name prior to index date
     create_pre_index_vars <- function(cohort_matched, eventdata, codelists) {
       pmap_dfc(list(eventdata, codelists$name), .progress=TRUE, \(eventdata, codelist_name) {
         

         # Get the first occurring event for each patient (that has an event) - 
         first_event <- eventdata |> 
           fselect(patid, obsdate) |> 
           fgroup_by(patid) |> 
           fmin()
         
         # Join the first events to the cohort of all patients
         temp <- cohort_matched[c("patid", "exposed", "indexdate")] |> 
           left_join(first_event, by=c("patid"))
         
         stopifnot(identical(temp[c("patid", "exposed")],cohort_matched[c("patid", "exposed")]))
         
         # Make a variable that is TRUE when the first event occurs before the index date
         temp[codelist_name] <- ifelse(temp$obsdate < temp$indexdate, TRUE, FALSE)
         temp[codelist_name][is.na(temp[codelist_name])] <- FALSE
         temp[codelist_name]
         # code suggeste by JM if want to keep more data- would need to also select when create "temp" for first time above
         #temp %>% select(contains(codelist_name), patid, exposed, indexdate, setid, pracid, gender, yob, indexdate_age, indexdate_age_cat)
       })
       
     }     
     
     
     ## OUTCOMES - Matching firstevent to patid - as create_pre_index_vars() but also keeps obsdate 
     # @return A of dataframes each with one row per patient, with patient ID and a TRUE/FALSE variable if they had the event in the variable name prior to index date
     create_obs_vars <- function(cohort_matched, eventdata, codelists) {
       pmap_dfc(list(eventdata, codelists$name), .progress=TRUE, \(eventdata, codelist_name) {
         
         # Get the first occurring event for each patient (that has an event) - 
         first_event <- eventdata |> 
           fselect(patid, obsdate) |> 
           fgroup_by(patid) |> 
           fmin() 
         
         # Join the first events to the cohort of all patients
         temp <- cohort_matched[c("patid", "exposed", "indexdate")] |> 
           left_join(first_event, by=c("patid"))
         
         stopifnot(identical(temp[c("patid", "exposed")],cohort_matched[c("patid", "exposed")]))
         
         # Outcome T/F flag: i.e. if any record of outcome - TRUE, else FALSE (regarldess of pre/post-indexdate)
         temp[codelist_name] <- ifelse(is.na(temp$obsdate), FALSE, TRUE)
        
         # # Make a variable that is TRUE when the first event occurs before the index date 
         # # NB TRUE means outcome before indexdate
         temp[paste0(codelist_name, "_preindex")] <- ifelse(temp$obsdate < temp$indexdate, TRUE, FALSE)
         temp[paste0(codelist_name, "_preindex")][is.na(temp[codelist_name])] <- FALSE
         temp[paste0(codelist_name, "_preindex")]

         temp[paste0(codelist_name, "_obsdate")] <- temp$obsdate

         
         # code suggested by JM if want to keep more data- would need to also select when create "temp" for first time above
         temp  |>  select(contains(codelist_name))
       })
     }
     
     
     # function to select events pre-post indexdate: 
     get_nearest_events <- function(df) {
       negs <- filter(df, !is.na(sh_indexdate_difftime) & sh_indexdate_difftime < 0) %>%
         slice_max(sh_indexdate_difftime, n = 1, with_ties = FALSE)
       
       poss <- filter(df, !is.na(sh_indexdate_difftime) & sh_indexdate_difftime >= 0) %>%
         slice_min(sh_indexdate_difftime, n = 1, with_ties = FALSE)
       
       nas <- filter(df, is.na(sh_indexdate_difftime))
       
       bind_rows(negs, poss, nas)
     }
     
     # version to CMD UPTO 5 years pre-index date (but not more than 5 years before) : 
     create_obs_vars_5yrspreindex <- function(cohort_matched, eventdata, codelists) {
       pmap_dfc(list(eventdata, codelists$name), .progress=TRUE, \(eventdata, codelist_name) {
         
         eventdata_indexdate <-  cohort_matched  |> select(patid, indexdate) |>  
           # join indexdate on to eventdata: 
           left_join(eventdata, by = "patid") |> 
           # gen time between indexdate and obsdate
           mutate(cmd_infection = as.numeric(difftime(obsdate, indexdate, units = "days"))) |> 
           # keep only events upto 5 years pre indexdate OR post-indexdate (i.e. outcome of interest): 
           filter(cmd_infection >= -1825) 
         
         
         # Get the first occurring event for each patient (that has an event within 5 yrs preindex) - just need one record, keeping first but could be any - 
         first_event <- eventdata_indexdate |> 
           select(patid, obsdate, cmd_infection) |> 
           group_by(patid) |> 
           slice_min(cmd_infection, with_ties = FALSE) 
         
         
         # Join the CMD within 5 years pre-index  to the cohort of all patients
         temp <- cohort_matched[c("patid", "exposed", "indexdate")] |> 
           left_join(first_event, by=c("patid"))
         
         # Outcome T/F flag: i.e. if any record of outcome - TRUE, else FALSE (regarldess of pre/post-indexdate)
         temp[codelist_name] <- ifelse(is.na(temp$obsdate), FALSE, TRUE)
         
         # # Make a variable that is TRUE when the first event occurs before the index date 
         # # NB TRUE means outcome before indexdate
         temp[paste0(codelist_name, "_5yrs_preindex")] <- ifelse(temp$obsdate < temp$indexdate, TRUE, FALSE)
         temp[paste0(codelist_name, "_5yrs_preindex")][is.na(temp[codelist_name])] <- FALSE
         temp[paste0(codelist_name, "_5yrs_preindex")]
         
         temp[paste0(codelist_name, "_5yrs_preindex_obsdate")] <- temp$obsdate
         
         
         # code suggested by JM if want to keep more data- would need to also select when create "temp" for first time above
         temp  |>  select(contains(codelist_name))
         
         
       })
     }
     
     
     
     # alternative dx censoring events: same as previous function but TRUE/FALSE is whether or not there is a record of the outcome, no assessment of it being before/after index date... 
     create_censor_vars <- function(cohort_matched, eventdata, codelists) {
       pmap_dfc(list(eventdata, codelists$name), .progress=TRUE, \(eventdata, codelist_name) {
         
         # Get the first occurring event for each patient (that has an event) - already filtered eventdata to be within study period (in eventdata target)
         first_event <- eventdata |> 
           fselect(patid, obsdate) |> 
           fgroup_by(patid) |> 
           fmin()
         
         # Join the first events to the cohort of all patients
         temp <- cohort_matched[c("patid", "exposed", "indexdate")] |> 
           left_join(first_event, by=c("patid"))
         
         stopifnot(identical(temp[c("patid", "exposed")],cohort_matched[c("patid", "exposed")]))
         
         temp[codelist_name] <- ifelse(!is.na(temp$obsdate), TRUE, FALSE)
         
         # and rename obsdate to reflect cenosring codelist name
         temp[paste0(codelist_name, "_obsdate")] <- temp$obsdate
         
         temp %>% select(contains(codelist_name)) 
     })
     }
     

     # save the data - Parquet - one file per infection: 
     
     list_write_parquet_and_return_path <- function(x, path, infection) {
       path <- paste0(path, infection, ".parquet")
       y <- as_tibble(x)
       write_parquet(y, path)
       return(path)
     }
     
     # save the data - csv - one file per infection: 
     
     write_and_return_path <- function(x, path, infection) {
       path <- paste0(path, infection,  ".csv")
       write_csv(x, path)
       return(path)
     }

     # save the data - RDS: 
     
     save_and_return_path_rds <- function(x, path) {
       path <- paste0(path, deparse(substitute(x)), ".RDS")
       saveRDS(x, file=path)
       return(path)
     }
     
     

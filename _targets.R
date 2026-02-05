# Setup -------------------------------------------------------------------
# This _targets.R file defines the {targets} pipeline.
# Run tar_make() to run the pipeline, tar_make(target) to run up to a defined target, and tar_read(target) to view the results.
library(targets)


# Define external paths
source("paths//0_globals_R.R")

# Set target-specific options such as packages.
tar_option_set(
  memory = "transient", garbage_collection = TRUE,
  packages = c(
    "arrow", #For reading parquet files
    "dplyr", 
    "tidyverse", 
    "data.table",
    "lubridate", # To manage dates
    "collapse", #For fast data management
    "dttr2", 
    "furrr", 
    "janitor", 
    "survival", # For survival analysis including Cox regression
    "broom",  # to clean regression output 
    "ggsurvfit")
) 

# Source all functions from the "R" folder
#sapply(list.files("R", full.names = TRUE), source, .GlobalEnv)
tar_source(files = "R")

# List of target objects.
list( 
  
  #-----------------------------------------------------------
  # Set up branching over all the infection datafiles and assign infection name as a "cohorts" variable
  #-----------------------------------------------------------
  
  tar_target(
    paths,
    list.files(path = raw_infections_files, pattern ="*.txt", full.names=TRUE)
  ),

  #-----------------------------------------------------------
  # Specifications: pracids to exclude as per CPRD Aurum Data Spec, v3.5, Date: 30 Aug 2024
  #-----------------------------------------------------------
  
  tar_target(
    name = pracid_exc,
    c("20024", "20469", "20803", "21112", "21444", "20036", "20487", "20804", "21118", "21451", "20091", "20552", "20822", "21172", "21529", "20171", "20554", "20868", "21173", "21553", "20178", "20640", "20908", "21277", "21558", "20202", "20717", "20912", "21281", "21585", "20254", "20734", "20996", "21331", "20389", "20737", "21001", "21334", "20430", "20740", "21015", "21390", "20452", "20790", "21078", "21430")
  ),
  
  # ┠ Study end date ----
  
  tar_target(study_start, as.Date("2007-01-01")),
  tar_target(study_end, as.Date("2024-06-15")),
  
  # plausible earliest record date  (for use when filter records of "ever" OUTCOME record, e.g. for CCI components)
  tar_target(study_start_plausible, as.Date("1930-01-01")),
  
  
  # Lookups
  
  tar_target(name = numunit, 
             command = read_tsv(path_numunit, col_types = cols(.default = "c"))
             ),
  
  # Codebrowsers - 
   tar_target(aurum_product, haven::read_dta(paste0(path_browsers, "\\CPRDAurumProduct.dta"))),
   tar_target(aurum_medical, haven::read_dta(paste0(path_browsers, "\\CPRDAurumMedical.dta"))),
  
  
  # Specify models for regression: 
  tar_target(
    model, 
    c(# Main models
      "matching_vars" = "",
      "A" = paste(c("", "imd"), collapse = " + ") ,
      "B" = paste(c("", "imd",  "smoking_status_simp", "alcohol_abuse", "obese", "cci_cat" ), collapse = " + "), 
      "C" = paste(c("", "imd",  "smoking_status_simp", "alcohol_abuse", "obese", "cci_cat", "ethnicity_5_na"), collapse = " + "))
  
  ),
  
  # target with just the fully adjusted model:
  tar_target(
    model_fully_adj, 
    c("C" = paste(c("", "imd",  "smoking_status_simp", "alcohol_abuse", "obese", "cci_cat", "ethnicity_5_na"), collapse = " + "))
  ),

  #-----------------------------------------------------------
  # Data Processing: Denom & infection data
  #-----------------------------------------------------------
  
  # txt files from define to parquet files for management and collect() ---------------
  
    # collect denominator data: 
  tar_target(
      name = denom_data_collect, 
      command = denom_data(file = raw_data_denominator),
      description = "Denom_read_in_data",
      format = "parquet"
  ),  

    # Denom data mx: adding study start and end dates & drop duplicated pracid: 
  tar_target(
    name = denom_mx, 
    command = denom_data_mx(data = denom_data_collect, 
                            pracid_exc = pracid_exc), 
    description = "Denom_mx_study_end_start",
    format = "parquet"
  ),
  
  tar_target(
    denom_mx_n,
    nrow(denom_mx)
  ), 

  # infection data collect():
    tar_target(
    name = inf_data_collect,
    command = inf_data(paths),
    description = "inf_read_in_data",
    pattern = map(paths),
    format = "parquet"
  ),

  # Merge in pracid (from denom data) and drop pracid to be excluded (n=46). 
  tar_target(
    name = inf_data_pracid_exc, 
    command = inf_pracid_exc(inf_data_collect = inf_data_collect, 
                             denom_data_collect = denom_data_collect,
                             pracid_exc = pracid_exc), 
    description = "merge pracid from denom and drop duplicates", 
    pattern = map(inf_data_collect), 
    format = "parquet"
  ),
  
  # Get random 10% sample for UTIs, SSTI, LRTI & Gastro (keeping whole sample for sepsis & mening/enceph)
  tar_target(
    name = random_sample, 
    command = sample(data = inf_data_pracid_exc), 
    description = "10% random sample (4 infs)", 
    pattern = map(inf_data_pracid_exc), 
    format = "parquet" 
  ), 
  
  # Join infection data to denom to get exposed & unexposed:
  tar_target(
    name = denom_inf_data,
    command = denom_inf(denom = denom_mx, inf = random_sample),
    description = "Denom+inf",
    pattern = map(random_sample),
    format = "parquet"
  ),
  


#-----------------------------------------------------------
# Create cohort of cases and comparators:
#-----------------------------------------------------------

# Get exposed patients to be available as controls prior to their date of first exposure
tar_target(
  name = cohort_matchable,
  command = create_cohort_matchable(cohort_eligible = denom_inf_data),
  description = "Cohort_mx_pre_expo",
  pattern = map(denom_inf_data),
  format = "parquet"
 ),

# Divide in to pracid groups (as matching done by pracid):
tar_target(
  cohort_matchable_grouped,
  command = split(cohort_matchable, cohort_matchable$pracid),
  description = "Split_practid",
  pattern = map(cohort_matchable),
  iteration="list"
),

# Run cohort matching algorithm: get upto 5 comparators for every 1 case
 tar_target(
  name = cohort_matched,
  command = create_cohort_matched(cohort_matchable_grouped = cohort_matchable_grouped),
  description = "Create_matched_cohort",
  pattern = map(cohort_matchable_grouped),
  format = "parquet"
),

# Generate summary statistics: # matches per case - need to make sure all >1 match
tar_target(
  name = summary_matches,
  command = cohort_stats(data = cohort_matched),
  description = "Table_matches",
  pattern = map(cohort_matched),
  iteration = "list"
), 

#-----------------------------------------------------------
# Get unique patids for extract:
#-----------------------------------------------------------

# save patids for the extracts (max nrows = 1M), seperate branches and name by cohorts-
tar_target(name = patids_for_extract_unique,
           command = cohort_matched  %>%
                     select(patid) %>%
                     filter(!duplicated(patid))

),

# 1m patids per file - as per extract instruction
tar_target(name = patids_for_extract_chunked_files,
           command = write_delim_compressed_in_chunks_and_return_paths(x = patids_for_extract_unique, max=1000000)
), 


#-----------------------------------------------------------
# Get summary stats on the numbers at each stage of cohort:
#-----------------------------------------------------------

# number in the final infection sample (once restricted to 10% for the 4 infs)
tar_target(name = table_n_inf_random_sample,
           command = table_inf_sample_n(random_sample), 
           pattern = map(random_sample)
), 

tar_target(name = table_n_cohort_matched, 
           command = table_cohort_matched_pop(cohort_matched), 
           pattern = map(cohort_matched)
), 


tar_target(name = table_n_cohort_matched_exp_unexp, 
           command = table_cohort_matched_exp_unexpo(cohort_matched), 
           pattern = map(cohort_matched)
), 


#-----------------------------------------------------------
# Saved the matched cohorts (with setid etc) as parquet files:
#-----------------------------------------------------------

tar_target(name = inf_names, 
           c("GE", "LRTI", "mening_enceph", "sepsis", "SSTI", "UTI")
), 

tar_target(name = write_parquet,
           command= write_parquet_and_return_path(x = cohort_matched, infection = inf_names, path = matched_infections),
           pattern = map(cohort_matched, inf_names)
),


#-----------------------------------------------------------
# Read back in the matched cohorts extract data (with setid etc) from parquet files:
#-----------------------------------------------------------
tar_target(name = paths_inf_parquet,
           list.files(path = matched_infections, pattern ="*.parquet", full.names=TRUE)
),

tar_target(name = inf_parquet,
           command = inf_parquet_open(paths_inf_parquet), 
           pattern = map(paths_inf_parquet), 
           format = "parquet"
), 

# SSTI enddate for extract was set to be 31-Dec rather than 31 March in the extract so correct and drop those with enddate before index date- drop accompanying matches as well
tar_target(inf_parquet_correct,
           command = inf_parquet |>    group_by(setid) %>%
                                       filter(all(indexdate <= as.Date('2023-03-31')) & all(enddate > indexdate))  |>
                      ungroup(), 
           pattern = map(inf_parquet),
           format = "parquet"
),

# FINAL INFECTION COHORT NUMBERS: 
tar_target(inf_parquet_correct_n,
           command = inf_parquet_correct |> tabyl(exposed,  show_na = TRUE) |> 
                                                 adorn_totals("row") |> 
                                                 adorn_percentages("col") |> 
                                                 adorn_pct_formatting(digits = 2)  |> 
                                                 adorn_ns(), 
           pattern = map(inf_parquet_correct),
           format = "parquet"
),


#-----------------------------------------------------------
# TYPE 2 data linkage request for CPRD data: to get additional ONS death data and HES records to classify as severe/non-severe
#-----------------------------------------------------------
# first need to read in the linkage_eligibility file- function - linkage_data
tar_target(name = read_linkage_data,
           command = linkage_data(linkage_eligibility),
           format = "parquet"
           
),

# merge in linkage eligiblity, keeping infections separate for now so can see numbers eligible per inf:  
tar_target(name = join_eligibility_flags,
           command = inf_parquet_correct  %>%
             left_join(read_linkage_data, by="patid") |>  
             filter(hes_apc_e == 1 | ons_death_e == 1 | lsoa_e == 1) ,
           pattern = map(inf_parquet_correct),
           format = "parquet"        
), 

# check numbers eligible by infection and exposure: 
tar_target(join_eligibility_flags_n,
           command = join_eligibility_flags |> tabyl(exposed,  show_na = TRUE) |> 
             adorn_totals("row") |> 
             adorn_percentages("col") |> 
             adorn_pct_formatting(digits = 2)  |> 
             adorn_ns(), 
           pattern = map(join_eligibility_flags),
           format = "parquet"
),

# look at numbers eligible for HES or ONS (ignoring LSOA): 
tar_target(join_eligibility_flags_n2,
           command = join_eligibility_flags |> filter(hes_apc_e == 1 | ons_death_e == 1) |> 
                                            tabyl(exposed,  show_na = TRUE) |> 
             adorn_totals("row") |> 
             adorn_percentages("col") |> 
             adorn_pct_formatting(digits = 2)  |> 
             adorn_ns(), 
           pattern = map(join_eligibility_flags),
           format = "parquet"
),


# remove duplicate patids (ppl duplicated between inf cohorts and included as unexpo before expo) 
# NB not mapping over infections, so get one list of patid with no duplicates
tar_target(data_for_linkage_cprd,
           command = join_eligibility_flags |> filter(!duplicated(patid)), 
           format = "parquet"
),

#  saving patid as per instructions for large requests, so reduce file size:
tar_target(name = patids_for_linkage_cprd,
           command = data_for_linkage_cprd   |>
             select(patid)
),

# check total n unique patients eligible for linkage (HES or ONS or LSOA flag)
tar_target(name = patids_for_linkage_cprd_n,
           command = patids_for_linkage_cprd   |>
                     nrow()
),
# save patids as text file - as per linked data request instruction - max 0.4M pats per file to comply with 2MB per file limit (once zipped)
tar_target(name = patids_for_linkage_chunked_files,
           command = write_delim_compressed_in_chunks_and_return_paths_linked(x = patids_for_linkage_cprd, max=400000)
),

# PRACTICE IDS FOR TYPE 2 LINKAGE-  not restricted to individuals eligible for linkage: 
# save just practice IDs 
tar_target(name = pracids_for_linkage_v2,
           command = inf_parquet_correct  |> 
             left_join(denom_mx |> select(pracid, patid), by="patid") |> 
             distinct(pracid) 
),

# save pracids as text file 
tar_target(name = pracids_for_linkage_chunked_files_v2,
           command = write_delim_compressed_in_chunks_and_return_paths_linked(x = pracids_for_linkage_v2, max=500000)
),


#-----------------------------------------------------------
# EXTRACTED DATA: Read  in the extracted Observation data 
#-----------------------------------------------------------

tar_target(name = files_obs,
           command = list.files(path = observation_files, pattern =".*?Extract_Observation.*?\\.parquet", recursive = TRUE, full.names = TRUE)
), 

tar_target(name = files_drug,
           command = list.files(path = observation_files, pattern = ".*?Extract_DrugIssue.*?\\.parquet", recursive = TRUE, full.names = TRUE)
), 

tar_target(name = files_pat,
           command = list.files(path = observation_files, pattern = ".*?Extract_Patient.*?\\.parquet", recursive = TRUE, full.names = TRUE)
), 

tar_target(name = files_prac,
           command = list.files(path = observation_files, pattern = ".*?Extract_Practice.*?\\.parquet", recursive = TRUE, full.names = TRUE)
), 


# check number of observations matches number in log- 
tar_target(name = files_obs_n, 
           command = open_dataset(files_obs) |> 
             nrow()
), 

# check number of drugs matches number in log-
tar_target(name = files_drug_n, 
           command = open_dataset(files_drug) |> 
             nrow()
), 

# check number of patients matches number in log- =
tar_target(name = files_pat_n, 
           command = open_dataset(files_pat) |> 
             nrow()
), 

# check number of practices matches number in log- 
tar_target(name = files_prac_n, 
           command = open_dataset(files_prac) |> 
             nrow()
), 

#-----------------------------------------------------------
# START COMBINING DATA:  
#-----------------------------------------------------------


tar_target( # Codelists used to extract eventdata
  codelists,
  tribble(~name, ~label, ~codevar, ~extract_from,
          
          #co-variates
          "alcohol_abuse", "Alcohol abuse", "medcodeid",	"observation",
          "alcoholism_drugs",	"Alcohol abuse drugs", "prodcodeid",	"drugissue",
          "cigarette_smoking", "Cigarette smoking", "medcodeid",	"observation", 
          "bmi", "BMI", "medcodeid",	"observation",
          "ethnicity", "Ethnicity", "medcodeid",	"observation") %>%  
    
          mutate(path=paste0(codelists_path, "\\Aurum\\", codevar, "\\", name, " aurum codelist.csv"), 
                 # store constents of each codelist csv file as list under codelist$full: 
                 full=map(path, ~read_csv(.x, col_types = cols(.default = "c"))), 
                 
                 # extract codevar (medcodeid/prodcodeid) and store under codelists$codes : 
                 codes=map2(path, codevar, ~read_csv(.x, col_types = cols(.default = "c"))[[.y]])) 
                 
  ), 

tar_target(bmi_codelist, codelists[which(codelists$name == "bmi"),]),
tar_target(alcohol_abuse_codelist, codelists[which(codelists$name == "alcohol_abuse"),]),
tar_target(alcoholism_drugs_codelist, codelists[which(codelists$name == "alcoholism_drugs"),]),
tar_target(cigarette_smoking_codelist, codelists[which(codelists$name == "cigarette_smoking"),]),
tar_target(ethnicity_codelist, codelists[which(codelists$name == "ethnicity"),]),


# Outcomes: CMD: separated by anxiety and depression and incident and historic
tar_target( # Codelists used to extract eventdata
  codelists_outcomes_cmd,
  tribble(~name, ~label, ~codevar, ~extract_from,
          
          #incident outcomes: 
          "anxiety_incident", "Anxiety Incident", "medcodeid", "observation", 
          "anxiety_historic", "Anxiety Historic", "medcodeid", "observation", 
          "depression_incident", "Depression Incident", "medcodeid", "observation", 
          "depression_historic", "Depression Historic", "medcodeid", "observation") %>%   
    
    mutate(path=paste0(codelists_path, "\\Aurum\\", codevar, "\\", name, " aurum codelist.csv"), 
           # store contents of each codelist csv file as list under codelist$full: 
           full=map(path, ~read_csv(.x, col_types = cols(.default = "c"))), 
           
           # extract codevar (medcodeid/prodcodeid) and store under codelists$codes : 
           codes=map2(path, codevar, ~read_csv(.x, col_types = cols(.default = "c"))[[.y]])) 
  
), 

tar_target(outcomes_cmd, codelists_outcomes_cmd$name),

# Outcomes: SMI: seperated by subtype and bipolar/schizophrenia/other: 

tar_target( # Codelists used to extract eventdata
  codelists_outcomes_smi,
  tribble(~name, ~label, ~codevar, ~extract_from,
          
          #incident outcomes: 
          "smi_bipolar_incident", "Bipolar Incident", "medcodeid", "observation", 
          "smi_bipolar_historic", "Bipolar Historic", "medcodeid", "observation", 
          "smi_other_psychosis_incident", "Other Psychosis Incident", "medcodeid", "observation", 
          "smi_other_psychosis_historic", "Other Psychosis Historic", "medcodeid", "observation", 
          "smi_schizo_incident", "Schizophrenia Historic", "medcodeid", "observation", 
          "smi_schizo_historic", "Schizophrenia Historic", "medcodeid", "observation") %>%   
     
    mutate(path=paste0(codelists_path, "\\Aurum\\", codevar, "\\", name, " aurum codelist.csv"), 
           # store constents of each codelist csv file as list under codelist$full: 
           full=map(path, ~read_csv(.x, col_types = cols(.default = "c"))), 
           
           # extract codevar (medcodeid/prodcodeid) and store under codelists$codes : 
           codes=map2(path, codevar, ~read_csv(.x, col_types = cols(.default = "c"))[[.y]])) 
  
), 

tar_target(outcomes_smi, codelists_outcomes_smi$name),

# Outcomes: non-fatal Self-Harm: separated by incident or historic 
tar_target( # Codelists used to extract eventdata
  codelists_outcomes_sh,
  tribble(~name, ~label, ~codevar, ~extract_from,
          
          #incident outcomes: 
          "self_harm_incident", "Self-harm Incident", "medcodeid", "observation", 
          "self_harm_historic", "Self-harm Historic", "medcodeid", "observation") %>%   
    
    mutate(path=paste0(codelists_path, "\\Aurum\\", codevar, "\\", name, " aurum codelist.csv"), 
           # store constents of each codelist csv file as list under codelist$full: 
           full=map(path, ~read_csv(.x, col_types = cols(.default = "c"))), 
           
           # extract codevar (medcodeid/prodcodeid) and store under codelists$codes : 
           codes=map2(path, codevar, ~read_csv(.x, col_types = cols(.default = "c"))[[.y]])) 
  
), 


# Suicide ICD-10 codes: 
tar_target( # Codelists used to extract eventdata
  codelists_outcomes_suicide_icd10,
  tribble(~name, ~label, ~codevar, ~extract_from,
          
          "suicide", "Suicide", "icd_10", "observation") %>%   
    
    mutate(path=paste0(codelists_path, "\\ICD_10\\", name, "_icd10.csv"), 
           # store constents of each codelist csv file as list under codelist$full: 
           full=map(path, ~read_csv(.x, col_types = cols(.default = "c"))), 
           
           # extract codevar (medcodeid/prodcodeid) and store under codelists$codes : 
           codes=map2(path, codevar, ~read_csv(.x, col_types = cols(.default = "c"))[[.y]])) 
  
), 
tar_target(suicide_hes_codes, unlist(codelists_outcomes_suicide_icd10$codes[codelists_outcomes_suicide_icd10$name=="suicide"])),

# Censors: 
tar_target( # Codelists used to extract eventdata
  codelists_censors,
  tribble(~name, ~label, ~codevar, ~extract_from,
          
          # censoring outcomes
          "organic_cmd", "Organic CMD", "medcodeid", "observation", 
          "organic_dep", "Organic Depression", "medcodeid", "observation",
          "dementia", "Dementia", "medcodeid", "observation", 
          "ocd", "OCD", "medcodeid", "observation", 
          "ptsd", "PTSD", "medcodeid", "observation", 
          "organic_psychosis", "Organic psychosis", "medcodeid", "observation") %>%  
    
    
    mutate(path=paste0(codelists_path, "\\Aurum\\", codevar, "\\", name, " aurum codelist.csv"), 
           # store constents of each codelist csv file as list under codelist$full: 
           full=map(path, ~read_csv(.x, col_types = cols(.default = "c"))), 
           
           # extract codevar (medcodeid/prodcodeid) and store under codelists$codes : 
           codes=map2(path, codevar, ~read_csv(.x, col_types = cols(.default = "c"))[[.y]])) 
  
), 

# Injecting drug use
tar_target( # Codelists used to extract eventdata for injectable drug use:
  codelist_inj_drug_use,
  tribble(~name, ~label, ~codevar, ~extract_from,
          
          #co-variates
          "injectable_drugs", "Injectable drug use", "medcodeid",	"observation") %>%  
    
    mutate(path=paste0(codelists_path, "\\Aurum\\", codevar, "\\", name, " aurum codelist.csv"), 
           # store constents of each codelist csv file as list under codelist$full: 
           full=map(path, ~read_csv(.x, col_types = cols(.default = "c"))), 
           
           # extract codevar (medcodeid/prodcodeid) and store under codelists$codes : 
           codes=map2(path, codevar, ~read_csv(.x, col_types = cols(.default = "c"))[[.y]])) 
  
), 

# Antimicrobial presciprtions
tar_target( # Codelists used to extract prodcodeid records for antimicrobial drugs:
  codelist_antimicrobials,
  tribble(~name, ~label, ~codevar, ~extract_from,
          
          #co-variates
          "antimicrobials", "Antimicrobial drug prescriptions", "prodcodeid",	"observation") %>%  
    
    mutate(path=paste0(codelists_path, "\\Aurum\\", codevar, "\\", name, " aurum codelist.csv"), 
           # store constents of each codelist csv file as list under codelist$full: 
           full=map(path, ~read_csv(.x, col_types = cols(.default = "c"))), 
           
           # extract codevar (medcodeid/prodcodeid) and store under codelists$codes : 
           codes=map2(path, codevar, ~read_csv(.x, col_types = cols(.default = "c"))[[.y]])) 
  
), 

# Charlson Comorbidities: 
tar_target(
  codelists_cci, 
  tribble(~name, ~label, ~codevar, ~extract_from,
          
          # CCI - codelists as per Quan 2011 Am J Epi paper:
          "myocardial_infarction", "Myocardial infarction","medcodeid",	"observation",
          "congestive_heart_failure", "Congestive Heart Failure", "medcodeid", "observation", 
          "peripheral_vascular_disease", "Peripheral vascular disease", "medcodeid",	"observation",
          "cerebrovascular_disease", "Cerebrovascular disease",	"medcodeid",	"observation",
          # dementia: generated separately in censoring code  
          "lung_disease", "Lung diseases", "medcodeid", "observation", # CPD 
          "rheumatoid_diseases", "Rheumatoid diseases", "medcodeid", "observation", 
          "peptic_ulcer", "Peptic ulcer", "medcodeid",	"observation",
          "mild_liver_disease", "Mild liver disease","medcodeid",	"observation",
          "diabetes_without_complications", "Diabetes: not complicated", "medcodeid", "observation", 
          "diabetes_with_complications", "Diabetes: complicated", "medcodeid", "observation", 
          "hemiplegia", "Hemiplegia  or paraplegia", "medcodeid", "observation", 
          # NB renal_disease codelist is the same as "chronic_kidney_disease" in FI so use that instead, not removing from here as has already run
          "renal_disease", "Renal Disease", "medcodeid", "observation", 
          "cancer_incl_lymphoma", "Cancer inc. lymphoma", "medcodeid",	"observation", 
          "moderate_severe_liver_disease", "Moderate or severe liver disease","medcodeid",	"observation",
          "metastatic_cancer", "Mestatic cancer", "medcodeid", "observation", 
          "hiv", "HIV/AIDs", "medcodeid", "observation") |> 
    mutate(path=paste0(codelists_path, "\\Aurum\\", codevar, "\\", name, " aurum codelist.csv"), 
           # store constents of each codelist csv file as list under codelist$full: 
           full=map(path, ~read_csv(.x, col_types = cols(.default = "c"))), 
           
           # extract codevar (medcodeid/prodcodeid) and store under codelists$codes : 
           codes=map2(path, codevar, ~read_csv(.x, col_types = cols(.default = "c"))[[.y]]))            
),



# these are the additional codlelists needed for FI (come components already captured in CCI)
tar_target( # Codelists needed to calculate Clegg Frailty Index, n=35, polypharmacy is calculate separately (not codelist)
  codelists_frailty, 
  tribble(~name, ~label, ~codevar, ~extract_from,
          #frailty - 
          "activity_limitation", "Activity limitation","medcodeid",	"observation",
          "anaemia", "Anaemia", "medcodeid",	"observation",
          "arthritis", "Arthritis", "medcodeid",	"observation",
          "atrial_fibrillation", "Atrial fibrillation", "medcodeid",	"observation",
         # already captured in CCI: 
         # "cerebrovascular_disease", "Cerebrovascular disease",	"medcodeid",	"observation",
          "chronic_kidney_disease", "Chronic Kidney Disease", "medcodeid", "observation",   
          "diabetes", "Diabetes mellitus", "medcodeid", "observation",
          "dizziness","Dizziness", "medcodeid",	"observation",
          "dyspnoea","Dyspnoea", "medcodeid",	"observation",
          "falls", "Falls", "medcodeid",	"observation",
          "foot_problems", "Foot problems", "medcodeid",	"observation",
          "fragility_fracture", "Fragility fracture", "medcodeid",	"observation",
          "hearing_loss", "Hearing loss",	"medcodeid",	"observation",
          "heart_failure", "Heart failure", "medcodeid",	"observation",
          "heart_valve_disease", "Heart valve disease", "medcodeid",	"observation",
          "housebound", "Housebound", "medcodeid",	"observation",
          "hypertension", "Hypertension", "medcodeid", "observation",
          "hypotension_syncope", "Hypotension / syncope", "medcodeid",	"observation",
          "ischaemic_heart_disease","Ischaemic heart disease","medcodeid",	"observation",
          "memory_cognitive_problems", "Memory and cognitive problems", "medcodeid", "observation", 
          "mobility_transfer_problems", "Mobility and transfer problems", "medcodeid",	"observation",
          "osteoporosis", "Osteoporosis", "medcodeid",	"observation",
          "parkinsonism_tremors", "Parkinsonism and tremor", "medcodeid",	"observation",
         # already captured in CCI: 
         # "peptic_ulcer", "Peptic ulcer", "medcodeid",	"observation",
         # "peripheral_vascular_disease", "Peripheral vascular disease", "medcodeid",	"observation",
          "requirement_for_care", "Requirement for care", "medcodeid",	"observation",
          "respiratory_disease", "Respiratory disease", "medcodeid",	"observation",
          "skin_ulcer", "Skin ulcer", "medcodeid",	"observation",
          "sleep_disorders", "Sleep Disorders", "medcodeid",	"observation", 
          "social_vulnerability", "Social vulnerability", "medcodeid",	"observation",
          "thyroid_disease", "Thyroid disease", "medcodeid",	"observation",
          "urinary_incontinence", "Urinary incontinence", "medcodeid",	"observation",
          "urinary_system_disease", "Urinary system disease", "medcodeid",	"observation",
          "visual_impairment", "Visual impairment", "medcodeid",	"observation",
          "weight_loss_anorexia", "Weight loss and anorexia", "medcodeid",	"observation") |> 
    mutate(path=paste0(codelists_path, "\\Aurum\\", codevar, "\\", name, " aurum codelist.csv"), 
           # store constents of each codelist csv file as list under codelist$full: 
           full=map(path, ~read_csv(.x, col_types = cols(.default = "c"))), 
           
           # extract codevar (medcodeid/prodcodeid) and store under codelists$codes : 
           codes=map2(path, codevar, ~read_csv(.x, col_types = cols(.default = "c"))[[.y]]))            
           
  ),


# # ┠ Extract eventdata -------------------------------------------------------

# ┠ Create analysis cohort --------------------------------------------------

## OUTCOMES -----------------------------------

# CMD outcomes: disaggregated by dep/anx and incident/historic

# extract outcome event data using codelists_outcomes 
# outcome_eventdata needs to be a list - 
tar_target(
  name = outcome_eventdata_cmd,
  command = open_dataset(files_obs)  |> 
    select(patid, obsdate, medcodeid) |> 
    filter(medcodeid %in% unlist(codelists_outcomes_cmd$codes)) |> 
    mutate(obsdate=as.Date(obsdate, "%d/%m/%Y")) |> 
    filter(obsdate >= study_start_plausible & obsdate<=study_end) |> 
    arrange(obsdate) |> 
    collect(),
  pattern = map(codelists_outcomes_cmd),
  iteration = "list",
  format = "parquet"
),

tar_target(
  name = outcomes_pre_index_var_cmd, 
  command = create_obs_vars(cohort_matched = inf_parquet_correct, 
                            eventdata = outcome_eventdata_cmd, 
                            codelists = codelists_outcomes_cmd),
  pattern = inf_parquet_correct,
  iteration = "list",
  format = "parquet",
  deployment = "main"
), 

# then need to generate CMD variable- 
tar_target(
  name = outcomes_cmd_first_event, 
  command = outcomes_pre_index_var_cmd |> 
                          mutate(cmd_incident_obsdate = pmin(anxiety_incident_obsdate, depression_incident_obsdate, na.rm = TRUE),
                                 cmd_historic_obsdate = pmin(anxiety_historic_obsdate, depression_historic_obsdate, na.rm = TRUE),
                                 cmd_preindex = ifelse((anxiety_incident_preindex == TRUE |anxiety_historic_preindex == TRUE | depression_incident_preindex == TRUE | depression_historic_preindex == TRUE ), TRUE, FALSE)),
  pattern = outcomes_pre_index_var_cmd,
  iteration = "list",
  format = "parquet",
  deployment = "main"
), 


# SMI outcomes: disaggregated by type and incident/historic

tar_target(
  name = outcome_eventdata_smi,
  command = open_dataset(files_obs)  |> 
    select(patid, obsdate, medcodeid) |> 
    filter(medcodeid %in% unlist(codelists_outcomes_smi$codes)) |> 
    mutate(obsdate=as.Date(obsdate, "%d/%m/%Y")) |> 
    filter(obsdate >= study_start_plausible & obsdate<=study_end) |> 
    arrange(obsdate) |> 
    collect(),
  pattern = map(codelists_outcomes_smi),
  iteration = "list",
  format = "parquet"
),


# SMI  observation variables: 
tar_target(
  name = outcomes_pre_index_var_smi,
  command = create_obs_vars(cohort_matched = inf_parquet_correct,
                            eventdata = outcome_eventdata_smi,
                            codelists = codelists_outcomes_smi),
  pattern = inf_parquet_correct,
  iteration = "list",
  format = "parquet",
  deployment = "main"
),

# then need to generate SMI composite variable- first of schizo/bipolar/other, with obsdate, and whether it is pre-indexdate or not
# also need an incident/historic flag for SMI composite- if have both for one day- historic trumps incident
tar_target(
  name = outcomes_smi_first_event,
  command = outcomes_pre_index_var_smi |>
    mutate(smi_incident_obsdate = pmin(smi_bipolar_incident_obsdate , smi_other_psychosis_incident_obsdate, smi_schizo_incident_obsdate, na.rm = TRUE),
           smi_historic_obsdate = pmin(smi_bipolar_historic_obsdate, smi_other_psychosis_historic_obsdate, smi_schizo_historic_obsdate, na.rm = TRUE),
           smi_preindex = ifelse((smi_bipolar_incident_preindex == TRUE |smi_bipolar_historic_preindex == TRUE | smi_other_psychosis_incident_preindex == TRUE | smi_other_psychosis_historic_preindex == TRUE | smi_schizo_incident_preindex == TRUE | smi_schizo_historic_preindex == TRUE), TRUE, FALSE)),
  pattern = outcomes_pre_index_var_smi,
  iteration = "list",
  format = "parquet",
  deployment = "main"
),


# INCIDENT SH events keeping the two eitherside of indexdate
tar_target(
  name = outcomes_pre_index_var_sh_incident,
  command = inf_parquet_correct  |> select(patid, indexdate, exposed) |>  
                # join indexdate on to eventdata: 
                left_join(outcome_eventdata_sh_incident |> select(patid, obsdate), by = "patid", relationship = "many-to-many") |> 
                
                # gen time between indexdate and obsdate
                mutate(sh_indexdate_difftime = as.numeric(obsdate - indexdate)) |> 
                
                # for each patid and each exposed status keep nearest records pre- and post-indexdate: 
                group_by(patid, exposed) |> 
                get_nearest_events()  |>       # apply function which keeps records nearest pre- and post- indexdate and keeps NA when no obsdate
                
                ungroup() , 
  pattern =  inf_parquet_correct,
  iteration = "list",
  format = "parquet",
  deployment = "main"
),

# INCIDENT SH: generate the flag variable for pre/post-obsdate and reshape to wide so one row per patid and exposed status again: 
tar_target(
  name = outcomes_pre_index_var_sh_incident_var,
  command = outcomes_pre_index_var_sh_incident  |> mutate(self_harm_incident_flag = case_when(
                                                          obsdate < indexdate ~ "self_harm_incident_preindex_obsdate", 
                                                          obsdate >= indexdate ~ "self_harm_incident_postindex_obsdate")) |> 
                                                  # need to drop the sh_indexdate_difftime var so can pivot_wider
                                                  select(-sh_indexdate_difftime) |> 
                                                       pivot_wider(names_from = self_harm_incident_flag, 
                                                                   values_from = obsdate) |> 
                                                       mutate(self_harm_incident_preindex = ifelse(!is.na(self_harm_incident_preindex_obsdate), TRUE, NA), 
                                                              self_harm_incident_postindex = ifelse(!is.na(self_harm_incident_postindex_obsdate), TRUE, NA)) |>  
                                                       select(-"NA"), 
  pattern =  map(outcomes_pre_index_var_sh_incident),
  iteration = "list",
  format = "parquet",
  deployment = "main"
),

# check no more than 1 row per patid & exposed and 2 rows per patid: 
tar_target(
  name = outcomes_pre_index_var_sh_incident_var_n,
  command = outcomes_pre_index_var_sh_incident_var  |> group_by(patid, exposed) |>  
                                                       mutate(n_patid_expo = n()) |> 
                                                       ungroup() |> 
                                                       group_by(patid) |>  
                                                       mutate(n_patid = n()) |> 
                                                       ungroup() |> 
                                                       tabyl(n_patid_expo,n_patid ), 
  pattern =  map(outcomes_pre_index_var_sh_incident_var),
  iteration = "list",
  format = "parquet",
  deployment = "main"
),


# HISTORIC SH events keeping the two eitherside of indexdate
tar_target(
  name = outcomes_pre_index_var_sh_historic,
  command = inf_parquet_correct  |> select(patid, indexdate, exposed) |>  
    # join indexdate on to eventdata: 
    left_join(outcome_eventdata_sh_historic |> select(patid, obsdate), by = "patid", relationship = "many-to-many") |> 
    
    # gen time between indexdate and obsdate
    mutate(sh_indexdate_difftime = as.numeric(obsdate - indexdate)) |> 
    
    # for each patid and each exposed status keep nearest records pre- and post-indexdate: 
    group_by(patid, exposed) |> 
    get_nearest_events()  |>       # apply function which keeps records nearest pre- and post- indexdate and keeps NA when no obsdate
    
    ungroup() , 
  pattern =  inf_parquet_correct,
  iteration = "list",
  format = "parquet",
  deployment = "main"
),



# HISTORIC SH: generate the flag variable for pre/post-obsdate and reshape to wide so one row per patid and exposed status again: 
tar_target(
  name = outcomes_pre_index_var_sh_historic_var,
  command = outcomes_pre_index_var_sh_historic  |> mutate(self_harm_historic_flag = case_when(
                                                          obsdate < indexdate ~ "self_harm_historic_preindex_obsdate", 
                                                          obsdate >= indexdate ~ "self_harm_historic_postindex_obsdate")) |> 
                                                          # need to drop the sh_indexdate_difftime var so can pivot_wider
                                                          select(-sh_indexdate_difftime) |> 
                                                          pivot_wider(names_from = self_harm_historic_flag, 
                                                                      values_from = obsdate) |> 
                                                          mutate(self_harm_historic_preindex = ifelse(!is.na(self_harm_historic_preindex_obsdate), TRUE, NA), 
                                                                 self_harm_historic_postindex = ifelse(!is.na(self_harm_historic_postindex_obsdate), TRUE, NA)) |>  
                                                          select(-"NA"), 
  pattern =  map(outcomes_pre_index_var_sh_historic),
  iteration = "list",
  format = "parquet",
  deployment = "main"
),

# Join together historic and incident SH events: 
tar_target(
  name = outcomes_pre_index_var_sh_vars_both,
  command = inf_parquet_correct |>  left_join(outcomes_pre_index_var_sh_incident_var, by = (c("patid", "exposed", "indexdate"))) |>  
                                    left_join(outcomes_pre_index_var_sh_historic_var, by = (c("patid", "exposed", "indexdate"))),
  pattern =  map(inf_parquet_correct, outcomes_pre_index_var_sh_incident_var, outcomes_pre_index_var_sh_historic_var),
  iteration = "list",
  format = "parquet",
  deployment = "main"
),

# add the ONS suicide death flag and cprd_death date from cprd to CPRD SH events
# make the key ddate variable to be earliest of CPRD and ONS (if both recorded), else take the one recorded
tar_target(
  name = outcomes_sh_suicide_cprd_ons2,
  command = bind_cols(outcomes_pre_index_var_sh_vars_both, ons_death_date_suicide, outcomes_death_cprd) |> 
    mutate(ddate = pmin(cprd_ddate, reg_date_of_death_ons, na.rm = TRUE)), 
  pattern = map(outcomes_pre_index_var_sh_vars_both, ons_death_date_suicide, outcomes_death_cprd),
  iteration = "list",
  format = "parquet",
  deployment = "main"
),


# make composite SH/suicide variable (CPRD SH/non-fatal suicide + ONS suicide):
tar_target(
  name = outcomes_sh_suicide_composite2,
  command = outcomes_sh_suicide_cprd_ons2 |> mutate(self_harm_suicide_incident_composite = ifelse((self_harm_incident_postindex == TRUE | suicide_icd10 == TRUE), TRUE, NA),  
                                                   suicide_icd10_obsdate = as.Date(ifelse(suicide_icd10== TRUE, reg_date_of_death_ons, NA)), 
                                                   self_harm_suicide_incident_composite_obsdate = pmin(self_harm_incident_postindex_obsdate, suicide_icd10_obsdate, na.rm = TRUE)) |> 
                                              # drop patid, exposed, indexdate, setid so dont have duplicate when bind in to cohort_wide
                                                  select(-c(patid, exposed, enddate, indexdate, setid)), 
  pattern = map(outcomes_sh_suicide_cprd_ons2),
  iteration = "list",
  format = "parquet"
),


# Censoring vars:  restrict to between study dates, as has to be after indexdate to be a censoring event (and therefore also after studystart)
tar_target(
  name = outcome_censors,
  command = open_dataset(files_obs)  %>%
    select(patid, obsdate, medcodeid) %>%
    filter(medcodeid %in% unlist(codelists_censors$codes)) %>%
    mutate(obsdate=as.Date(obsdate, "%d/%m/%Y")) %>%
    filter(obsdate>=study_start & obsdate<=study_end) %>%
    arrange(obsdate) %>%
    collect(),
  pattern = map(codelists_censors),
  iteration = "list",
  format = "parquet"
),


# censoring vars flagged regardless of pre/post-index date
tar_target(
  name = outcome_censors_var, 
  command = create_censor_vars(cohort_matched = inf_parquet_correct, 
                            eventdata = outcome_censors, 
                            codelists = codelists_censors), 
  pattern = inf_parquet_correct,
  iteration = "list",
  format = "parquet",
  deployment = "main"
), 




## VARIABLES -----------------------------------

# Sex and age (inc categorical age): 
tar_target(name = sex_age,
           command = join_sex_age(inf = inf_parquet_correct, denom = denom_mx),
           pattern = map(inf_parquet_correct)
),

   
#┠ ALCOHOL ABUSE -----------------------------------
# get medcodeid and prodcodeid record for study participants during study dates: 
tar_target(
  name = alc_eventdata, 
  command = alg_alcoholism(cl = bind_rows(alcohol_abuse_codelist, alcoholism_drugs_codelist), 
                           files = files_obs , 
                           drug_files = files_drug, 
                           study_end = study_end)
),


tar_target(
  name = alc_pre_index_var, 
  command = create_pre_index_vars(cohort_matched = inf_parquet_correct, 
                                  eventdata = list(alc_eventdata), 
                                  codelists = alcohol_abuse_codelist),
           pattern = inf_parquet_correct,
           iteration = "list",
           format = "parquet",
           deployment = "main"
),

#┠ Smoking status---------------------------
tar_target(smok_medcodes, unlist(cigarette_smoking_codelist$codes)),

# select observation records that match smoking medcodeid-
tar_target(
  smok_eventdata, 
  open_dataset(files_obs) |> 
    select(patid, obsdate, medcodeid) |> 
    filter(medcodeid %in% smok_medcodes) |> 
    mutate(obsdate=as.Date(obsdate, "%d/%m/%Y")) |> 
    filter(obsdate<=study_end) |> 
    arrange(obsdate) |> 
    collect() |> 
    right_join(cigarette_smoking_codelist$full[[1]][c("medcodeid", "smokstatus")]),
  format="parquet"
),

# restrict to 5 years pre-index date: looking back 5 years (1826.25 days) (and then algorithm takes nearest to indexdate)
# run smoking algorithm: 
tar_target(
  smok_var,
  alg_smoking(cohort_matched = inf_parquet_correct, 
              smok_eventdata = smok_eventdata),
  pattern = inf_parquet_correct,
  iteration = "list",
  format = "parquet"
),

#┠ Injectable drug use ---------------------------
tar_target(inj_drugs_medcodes, unlist(codelist_inj_drug_use$codes)),

tar_target(
  inj_drugs_eventdata,
  open_dataset(files_obs) |>
    select(patid, obsdate, medcodeid) |>
    filter(medcodeid %in% inj_drugs_medcodes) |>
    mutate(obsdate=as.Date(obsdate, "%d/%m/%Y")) |>
    filter(obsdate<=study_end) |>
    arrange(obsdate) |>
    collect() |>
    right_join(codelist_inj_drug_use$full[[1]]),
  format="parquet"
),

tar_target(
  name = inj_drugs_var, 
  command = create_pre_index_vars(cohort_matched = inf_parquet_correct, 
                                  eventdata = list(inj_drugs_eventdata), 
                                  codelists = codelist_inj_drug_use),
  pattern = inf_parquet_correct,
  iteration = "list",
  format = "parquet",
  deployment = "main"
),

#┠ BMI -----------------------------------
tar_target(bmi_medcodes, unlist(bmi_codelist$codes)),
tar_target(
  bmi_eventdata, 
  open_dataset(files_obs) |> 
    select(patid, obsdate, medcodeid, numunitid, value) |> 
    filter(medcodeid %in% bmi_medcodes) |> 
    mutate(obsdate=as.Date(obsdate, "%d/%m/%Y")) |> 
    filter(obsdate<=study_end) |> 
    #arrange(obsdate) |> 
    collect() |> 
    right_join(bmi_codelist$full[[1]]) |> 
    left_join(numunit, by="numunitid"),
  format="parquet"
),

tar_target(
  name = bmi_var,
  command = alg_bmi(cohort_matched = inf_parquet_correct, bmi_eventdata = bmi_eventdata),
  pattern = inf_parquet_correct,
  iteration = "list",
  format = "parquet"
), 


#┠ Ethnicity--------------------------------------
tar_target(
  ethnicity,
  alg_ethnicity(codelists_for_define = ethnicity_codelist, 
                files = files_obs)
),

tar_target(
  ethnicity_vars,
  inf_parquet_correct |> 
    left_join(ethnicity, by="patid") |> 
    select(ethnicity_5, ethnicity_16) |> 
    mutate(ethnicity_5=as_factor(ethnicity_5), 
           ethnicity_5_na=ethnicity_5,
           ethnicity_5=fct_na_value_to_level(ethnicity_5,  level="5"),
           ethnicity_5=factor(ethnicity_5, labels = c("white", "south_asian", "black", "mixed", "other", "unknown")), 
           ethnicity_16=as_factor(ethnicity_16),
           ethnicity_16_na=ethnicity_16,
           ethnicity_16=fct_na_value_to_level(ethnicity_16,  level="17"),
           ethnicity_4_na=fct_collapse(ethnicity_5_na, white="0",south_asian="1", black="2", mixed_other=c("3", "4"), unknown = "5")),
           
  pattern = inf_parquet_correct,
  iteration = "list",
  format = "parquet"
),


#Consultation in the year before cohort entry
tar_target(
  cons_in_year_pre_index,
  alg_cons_in_year_pre_index(inf_parquet_correct, files_obs, aurum_medical),
  pattern = inf_parquet_correct,
  iteration = "list",
  format = "parquet"
),

#┠ Polypharmacy--------------------------------

tar_target(
  polypharmacy,
  alg_polypharmacy(inf_parquet_correct, files_drug, aurum_product),
  pattern = inf_parquet_correct,
  iteration = "list",
  format = "parquet"
),

# ┠ Antibiotic prescribing drugs data--------------------------------
# restricting to study start and enddates as will match on indexdate
tar_target(
  antimicrobial_eventdata,
  open_dataset(files_drug) |>
    select(patid, issuedate, prodcodeid) |>
    filter(prodcodeid %in% unlist(codelist_antimicrobials$codes)) |>
    mutate(issuedate=as.Date(issuedate, "%d/%m/%Y")) |>
    filter(issuedate>=(study_start) & issuedate<= (study_end)) |>
    arrange(issuedate) |>
    collect() |>
    right_join(codelist_antimicrobials$full[[1]]),
  format="parquet"
),


tar_target(
  name = antimicrobial_var,
  command = antimicrobial_prescriptions(cohort_matched = inf_parquet_correct, 
                                        eventdata= antimicrobial_eventdata), 
  pattern = inf_parquet_correct,
  iteration = "list",
  format = "parquet",
  deployment = "main"
),

# ┠ Frailty  & CCI data--------------------------------

## all CCI data - record of "ever" ie not restricted to study dates: 
tar_target( # Get eventdata for every CCI codelist
  cci_eventdata,
  open_dataset(files_obs) |> 
    select(patid, obsdate, medcodeid) |> 
    filter(medcodeid %in% unlist(codelists_cci$codes)) |> 
    mutate(obsdate=as.Date(obsdate, "%d/%m/%Y")) |> 
    filter(obsdate<=study_end) |> 
    arrange(obsdate) |> 
    collect(),
  pattern = map(codelists_cci),
  iteration = "list",
  format = "parquet"
),

# create CCI vars - dropping cancer_incl_lymphoma as had scientific notation problems, corrected in cancer_eventdata targets: 
tar_target(
  name = cci_vars, 
  command = create_pre_index_vars(cohort_matched = inf_parquet_correct, 
                                  eventdata = cci_eventdata, 
                                  codelists = codelists_cci),
  pattern = inf_parquet_correct,
  iteration = "list",
  format = "parquet",
  deployment = "main"
  ),

## all FI data - record of "ever" ie not restricted to study dates: 
# these are the FI variables not already captured in CCI
tar_target( # Get eventdata for every frailty codelist
  frailty_eventdata,
  open_dataset(files_obs) |> 
    select(patid, obsdate, medcodeid) |> 
    filter(medcodeid %in% unlist(codelists_frailty$codes)) |> 
    mutate(obsdate=as.Date(obsdate, "%d/%m/%Y")) |> 
    filter(obsdate<=study_end) |> 
    arrange(obsdate) |> 
    collect(),
  pattern = map(codelists_frailty),
  iteration = "list",
  format = "parquet"
),

tar_target(
  name = fi_var, 
  command = create_pre_index_vars(cohort_matched = inf_parquet_correct, 
                                  eventdata = frailty_eventdata, 
                                  codelists = codelists_frailty),
  pattern = inf_parquet_correct,
  iteration = "list",
  format = "parquet",
  deployment = "main"
  ),


## Severity flag: ICD-10 code in hes_diagnosis_hosp within 28 days eitherside of indexdate (for exposed)
# join HES infection data to the CPRD infection data: 
tar_target(
  name = hes_severity,
  command = inf_parquet_correct |> left_join(hes_diagnosis_hosp_inf, by = "patid", relationship = "many-to-many") |>
                                    mutate(admidate = as.Date(admidate, format = "%Y/%m/%d"),
                                           hes_index_diff = as.numeric(difftime(admidate, indexdate, units = "days")))  |>
                                    group_by(patid, exposed) |>
                                    slice_min(hes_index_diff) |>
                                    ungroup() |>
                                    # generate severe flag 0/1 for exposed if HES event is 28 days eitherside of CPRD event, NA for unexposed:
                                      mutate(severe_flag = case_when(
                                        # HES before indexdate: 
                                        exposed == 1 & hes_index_diff >= -28 & hes_index_diff <= 0  ~ 1,
                                        # HES after indexdate: 
                                        exposed == 1 & hes_index_diff >= 1   & hes_index_diff <= 28 ~ 2,
                                        exposed == 1 & hes_index_diff > 28 ~ 0,
                                        exposed == 1 & hes_index_diff < -28 ~ 0,
                                        exposed == 1 & is.na(hes_index_diff) ~ 0, 
                                        exposed == 0 ~ NA_real_), 
                                        severe_flag = factor(
                                          severe_flag,
                                          levels = c(1, 2, 0),
                                          labels = c("hes pre-index", "hes post-index", "not severe")
                                        )) |>    
                                    # ppl with multiple HES codes for same infection type on same day: keep only one record:
                                    select(patid, exposed, ICD, severe_flag, hes_index_diff, admidate) |>

                                    group_by(patid, exposed) |>
                                    filter(row_number()==1),
  pattern = map(inf_parquet_correct,
                  cross(hes_diagnosis_hosp_inf))
),


tar_target(
  name = hes_severity_n, 
  command = hes_severity |> filter(exposed == 1) |> 
                                tabyl(exposed, severe_flag) |> 
                                adorn_percentages("row") %>%
                                adorn_pct_formatting(digits = 2), 
  pattern = map(hes_severity),  
  iteration = "list",
  format = "parquet"
  ), 

# join on to inf_parquet_correct: 
tar_target(
  name = hes_severity_var, 
  command = inf_parquet_correct |>  left_join(hes_severity, by = c("patid", "exposed")), 
   pattern = map(inf_parquet_correct,
                 cross(hes_severity)),  
   iteration = "list",
   format = "parquet"
  ), 

# doing sepsis for every branch as need sepsis flag for the other infections (except mening/enceph)
tar_target(
  name = hes_severity_sepsis,
  command = inf_parquet_correct |> left_join(hes_diagnosis_hosp_sepsis, by = "patid", multiple  = "all") |>
                                    mutate(admidate_sepsis = as.Date(admidate, format = "%Y/%m/%d"),
                                           # removed abs() so that difftime is negative if sepsis before indexdate, pos if after
                                           hes_index_diff_sepsis = as.numeric(difftime(admidate_sepsis, indexdate, units = "days")))  |>
                                    group_by(patid, exposed) |>
                                    slice_min(hes_index_diff_sepsis) |>
                                    ungroup() |>
                                    # generate severe flag 0/1 for exposed if HES event is 28 days eitherside of CPRD event, NA for unexposed:
                                    mutate(severe_flag_sepsis = case_when(
                                      # sepsis before indexdate: 
                                      exposed == 1 & hes_index_diff_sepsis >= -28 & hes_index_diff_sepsis <= 0  ~ 1,
                                      # sepsis after indexdate: 
                                      exposed == 1 & hes_index_diff_sepsis >= 1   & hes_index_diff_sepsis <= 28 ~ 2,
                                      exposed == 1 & hes_index_diff_sepsis > 28 ~ 0,
                                      exposed == 1 & hes_index_diff_sepsis < -28 ~ 0,
                                      exposed == 1 & is.na(hes_index_diff_sepsis) ~ 0, 
                                      exposed == 0 ~ NA_real_), 
                                          severe_flag_sepsis = factor(
                                            severe_flag_sepsis,
                                            levels = c(1, 2, 0),
                                            labels = c("hes pre-index", "hes post-index", "not severe")
                                      )) |>    
                                    # ppl with multiple HES codes for same infection type on same day: keep only one record:
                                    select(patid, exposed, ICD, severe_flag_sepsis, hes_index_diff_sepsis, admidate_sepsis) |>
                                    
                                    group_by(patid, exposed) |>
                                    filter(row_number()==1),
  pattern = map(inf_parquet_correct),
  iteration = "list",
  format = "parquet"
),

#join on toinf_parquet_correct: 
tar_target(
  name = hes_severity_sepsis_var, 
  command = inf_parquet_correct |>  left_join(hes_severity_sepsis, by = c("patid", "exposed")), 
  pattern = map(inf_parquet_correct,
                cross(hes_severity_sepsis)),  
  iteration = "list",
  format = "parquet"
), 

# check numbers: 
tar_target(
  name = hes_severity_sepsis_n, 
  command = hes_severity_sepsis |> filter(exposed == 1) |> 
    tabyl(exposed, severe_flag_sepsis) |> 
    adorn_percentages("row") %>%
    adorn_pct_formatting(digits = 2), 
  pattern = map(hes_severity_sepsis),  
  iteration = "list",
  format = "parquet"
), 

# join HES and infection and severity flag vars and gen composite severe_inf_composite var:
tar_target(
  name = hes_severity_composite, 
  command = hes_severity_var |>  left_join(hes_severity_sepsis_var, by = c("patid", "exposed")) |> 
                                  mutate(severe_flag_composite = case_when(
                                    exposed == 1 & severe_flag == "hes pre-index" ~ 1,
                                    exposed == 1 & severe_flag_sepsis == "hes pre-index" ~ 1,
                                    exposed == 1 & severe_flag ==  "hes post-index" ~ 2, 
                                    exposed == 1 & severe_flag_sepsis ==  "hes post-index" ~ 2, 
                                    exposed == 1 & severe_flag_sepsis == "not severe" & severe_flag == "not severe" ~ 0,
                                    exposed == 0 ~ NA_real_), 
                                    severe_flag_composite = factor(
                                      severe_flag_composite,
                                      levels = c(1, 2, 0),
                                      labels = c("hes pre-index", "hes post-index", "not severe")
                                    )) |> 
                                  select(c(ICD.x, ICD.y, starts_with("severe_flag"), starts_with("hes_index_diff"))) |> 
                                  rename(ICD_inf = ICD.x) |> 
                                  rename(ICD_sepsis = ICD.y),    
    pattern = map(hes_severity_var,
                cross(hes_severity_sepsis_var)),  
    iteration = "list",
    format = "parquet"
),

# check numbers: 
tar_target(
  name = hes_severity_composite_n, 
  command = hes_severity_composite |> filter(exposed == 1) |> 
    tabyl(exposed, severe_flag_composite) |> 
    adorn_percentages("row") %>%
    adorn_pct_formatting(digits = 2), 
  pattern = map(hes_severity_composite),  
  iteration = "list",
  format = "parquet"
), 

# need linkage eligibilty var to subset for "severe/non-severe" infection secondary analyses: 

tar_target( 
  name = severe_cohort_flag, 
  command = inf_parquet_correct |>  left_join(read_linkage_data, by="patid") |>  
                                    mutate(hes_linkage_eligible = if_else(hes_apc_e == 1, TRUE, FALSE), 
                                           # those not even appearing in linkage eligibility file: code to not eligible
                                           hes_linkage_eligible = if_else(is.na(hes_linkage_eligible), FALSE, hes_linkage_eligible)) |> 
                                    select(hes_linkage_eligible),
  pattern = map(inf_parquet_correct),  
  iteration = "list",
  format = "parquet"
  ), 

# number in severe_cohort: 
tar_target(
  name = severe_cohort_flag_n, 
  command = severe_cohort_flag |>  tabyl(hes_linkage_eligible), 
  pattern = map(severe_cohort_flag),  
  iteration = "list",
  format = "parquet"
), 


# ┠ Create analysis cohort --------------------------------------------------

# bind variables to cohort data: and run frailty & CCI algorithms 

tar_target(
  name = cohort_wide,
  command = bind_cols(inf_parquet_correct, sex_age, alc_pre_index_var, smok_var,  inj_drugs_var, bmi_var, ethnicity_vars, cons_in_year_pre_index, polypharmacy, antimicrobial_var, cci_vars, fi_var, outcomes_cmd_first_event, outcomes_smi_first_event, outcomes_sh_suicide_composite2, outcome_censors_var, imd, , hes_severity_composite, severe_cohort_flag) |> 
    alg_frailty() |>  
    alg_cci(), 
  pattern = map(inf_parquet_correct, sex_age, alc_pre_index_var, smok_var,  inj_drugs_var, bmi_var, ethnicity_vars, cons_in_year_pre_index, polypharmacy, antimicrobial_var, cci_vars, fi_var, outcomes_cmd_first_event, outcomes_smi_first_event, outcomes_sh_suicide_composite2, outcome_censors_var, imd, hes_severity_composite, severe_cohort_flag),
  iteration = "list",
  format = "parquet"
),


tar_target(
  name = cohort_wide_n, 
  command = cohort_wide %>%  tabyl(exposed) |> 
                                      adorn_totals("row"), 
  pattern = map(cohort_wide),
  iteration = "list",
  format = "parquet"
  
),

# number in severe_cohort - by exposure status: 
tar_target(
  name = severe_cohort_flag_exposed_n, 
  command = cohort_wide |>  tabyl(exposed, hes_linkage_eligible) |> 
    adorn_totals("row"), 
  pattern = map(cohort_wide),  
  iteration = "list",
  format = "parquet"
), 



# Tidy up variable labels and factor labels: focussing on those in table 1 - 
# drop those with "indeterminate" gender- there are just 2 ppl in the GE cohort: 

tar_target(
  name = cohort_wide_labels, 
  command = cohort_wide |>  filter(gender != "indeterminate") |> 
                              mutate(obese_na = factor(obese_na, labels = c("not obese" = "Not obese", 
                                                                          "obese" = "Obese")), 
                                   ethnicity_5 = factor(ethnicity_5, labels = c("white" = "White",
                                                                                "south_asian" = "South Asian", 
                                                                                "black" = "Black", 
                                                                                "mixed" = "Mixed", 
                                                                                "other" = "Other", 
                                                                                "unkown" ="Unknown")),
                                   ethnicity_5_na = ifelse(ethnicity_5_na == "5", NA, as.character(ethnicity_5_na)),
                                  ethnicity_5_na = factor(ethnicity_5_na), 
                                   imd = ordered(imd, labels =  c("1" = "1, least deprived", 
                                                                  "2" = "2", 
                                                                  "3" = "3", 
                                                                  "4" = "4", 
                                                                  "5" = "5, most deprived")), 
                                   smoking_status_simp = factor(smoking_status_simp, labels = c("non_smoker" = "Non-smoker", 
                                                                                    "smoker" = "Current- or ex-smoker"))) |> 
                                  # drop the unused levels within gender var: 
                                  droplevels() |> 
                                  mutate(gender = factor(gender, labels = c("male" = "Male", 
                                                                      "female" = "Female")))       |> 
                            labelled::set_variable_labels(ethnicity_5 = "Ethnicity", 
                                                          obese_na = "Obesity", 
                                                          smoking_status_simp = "Smoking", 
                                                          alcohol_abuse = "Harmful alcohol use", 
                                                          imd = "Index of multiple deprivation quintile", 
                                                          obese_na = "Body mass index", 
                                                          indexdate_age_cat = "Age", 
                                                          indexdate_age = "Age", 
                                                          n_cons_in_year_pre_index	= "Consultations in year before cohort entry", 
                                                          cci_cat = "Charlson Co-morbidity Index Category", 
                                                          severe_flag = "HES record of infection within 28 days (same infection)", 
                                                          severe_flag = "HES record of sepsis within 28 days", 
                                                          severe_flag_composite = "HES record of severe infection  within 28 days (same infection or sepsis"), 
  pattern = map(cohort_wide),
  iteration = "list",
  format = "parquet"
  
),


# Update enddate: 
# now have ONS death data drop those with death before infection - one is recorded incorrectly (they are typically v close in date) but cant know which so dropping these ppl entirely (not many ppl)
tar_target( 
  name = cohort_wide_enddate, 
  command = cohort_wide_labels |> mutate(enddate = pmin(ddate, enddate, na.rm = TRUE)), 
  pattern = map(cohort_wide_labels),
  iteration = "list",
  format = "parquet"
  
),

# get setids for exposed people with enddate<indexdate- need to drop whole set: 
tar_target(
  name = setids_ddate_preindexdate,
  command = cohort_wide_enddate |> filter(exposed == 1  & enddate < indexdate)  |>
    select(setid)  |>
    unique(), 
  pattern = map(cohort_wide_enddate),
  iteration = "list",
  format = "parquet"
), 

# if they are exposed: drop entire set, if they are unexposed, drop just the individual
tar_target(
  name = cohort_wide_final,
  command =  cohort_wide_enddate  |> filter(!setid %in% setids_ddate_preindexdate$setid) |>             # drop all rows with those setids
    filter(!(enddate < indexdate & exposed == 0)), # drop exposed==0 rows where enddate < indexdate
  pattern = map(cohort_wide_enddate, 
                cross(setids_ddate_preindexdate)),
  iteration = "list",
  format = "parquet"
),

tar_target(
  name = cohort_wide_final_n,
  command = cohort_wide_final %>%  tabyl(exposed)|> 
    adorn_totals("row"), 
  pattern = map(cohort_wide_final),
  iteration = "list",
  format = "parquet"
  
),

# Checking numbers of unexposed who have a record of infection (possibly from being in the 90% of infections not excluded from the denominator file - only an issue for the 4 infs took a 10% sample of):
# NB: the cohort_wide_final and inf_data_pracid_exc are in the same order of infections so can branch over them.
tar_target(
  name = unexposed_checks,
  command = cohort_wide_final |>  # join the entire infection cohort (before restricted to 10%):
                                 left_join(inf_data_pracid_exc, by = "patid") |>
                                        # flag unexposed ppl with an infection record in this original 100% of cases cohort - ie unexposed who appear in the orig cohort
                                      mutate(unexposed_inf = if_else(exposed == 0 & !is.na(infection), 1, NA),
                                            # rename the infection indexdate as it is not matching date, but infection date for the 90% not included as cases
                                            infection_censor_obsdate = indexdate.y,
                                            # keep orig indexdate variable
                                             indexdate = indexdate.x,
                                            # flag when infection date is before matching date:
                                             infection_censor_pre_matchdate = if_else(infection_censor_obsdate < indexdate, 1, 0),
                                            # recode infection_censor_obsdate to NA when exposed == 1 (as this is matching date and is fine, info captured indexdate
                                             infection_censor_obsdate = if_else(exposed == 1, NA, infection_censor_obsdate),
                                            # recode infection_censor_obsdate to NA for sepsis and mening cohorts as not a problem for them and dont want to censor on this date as it is the indexdate:
                                            infection_censor_obsdate = if_else(infection == "sepsis" | infection == "mening_encep ", NA, infection_censor_obsdate)) |>
                                      select(-c(indexdate.y, indexdate.x, infection)),
  pattern = map(cohort_wide_final, inf_data_pracid_exc),
  iteration = "list",
  format = "parquet"

),

# # tabulate numbers of controls who were in the original patient define:
tar_target(
  name = unexposed_checks_n,
  command = unexposed_checks |> tabyl(exposed, infection, show_na = TRUE),
  pattern = map(unexposed_checks),
  iteration = "list",
  format = "parquet"

),


# tabulate numbers of controls who were in the original patient define & who have record of infection pre-mathcing date: THESE ARE TO BE DROPPPED I THINK....
tar_target(
  name = unexposed_prematch_n,
  command = unexposed_checks |> filter(exposed == 0 & !is.na(infection)) |>
                             tabyl(infection_censor_pre_matchdate),
  pattern = map(unexposed_checks),
  iteration = "list",
  format = "parquet"

),

# drop unexposed who had an infection record pre-matching (indexdate) and recode infection_censor_obsdate to NA for exposed people- only need for unexposed.  
tar_target(
  name = cohort_wide_final_2, 
  command = unexposed_checks |> filter(is.na(infection_censor_pre_matchdate) | infection_censor_pre_matchdate==0) |> 
                                select(-c(infection_censor_pre_matchdate, unexposed_inf)), 
  pattern = map(unexposed_checks),
  iteration = "list",
  format = "parquet"
  
),

# Check total sample size: 
tar_target(
  name = cohort_wide_final_2_n,
  command = cohort_wide_final_2 %>%  tabyl(exposed)|> 
    adorn_totals("row"), 
  pattern = map(cohort_wide_final_2),
  iteration = "list",
  format = "parquet"
  
),


# gen total number per setid var to check how many matches ppl have now: 

tar_target(
  name = updated_set_num, 
  command = unexposed_checks |> filter(is.na(infection_censor_pre_matchdate) | infection_censor_pre_matchdate==0) |>  
    group_by(setid) |> 
    mutate(
      ppl_per_set = n(),              # Number of rows per setid
      is_first_row = row_number() == 1  # Flag TRUE for the first row within each setid
    ) |>
    ungroup(), 
  pattern = map(unexposed_checks),
  iteration = "list",
  format = "parquet"
  
),
# tabulate people in the sets: 
tar_target(
  name = updated_set_num_n, 
  command = updated_set_num |> filter(is_first_row == TRUE) |> 
    tabyl(ppl_per_set), 
  pattern = map(updated_set_num),
  iteration = "list",
  format = "parquet"
  
),

# check numbers of matches in the original cohort: 

# drop unexposed who had an infection record pre-mathcing (indexdate) and regen number of ppl per setid: 
tar_target(
  name = cohort_wide_final_set_num, 
  command = cohort_wide_final |>  group_by(setid) |> 
                                      mutate(
                                        ppl_per_set = n(),              # Number of rows per setid
                                        is_first_row = row_number() == 1  # Flag TRUE for the first row within each setid
                                      ) |>
                                      ungroup(), 
  pattern = map(cohort_wide_final),
  iteration = "list",
  format = "parquet"
  
),

# tabulate people in the sets: 
tar_target(
  name = cohort_wide_final_set_num_n, 
  command = cohort_wide_final_set_num |> filter(is_first_row == TRUE) |> 
    tabyl(ppl_per_set), 
  pattern = map(cohort_wide_final_set_num),
  iteration = "list",
  format = "parquet"
  
),


## Severity analysis: look at number people with HES record after CPRD reocrd:
tar_target(
  name = severity_timing_explore, 
  command = cohort_wide_final_2 |> filter(hes_linkage_eligible == TRUE & exposed == 1) |>  
                                   tabyl(severe_flag_composite), 
  pattern = map(cohort_wide_final_2),
  iteration = "list",
  format = "parquet"
  
),

# median gap between indexdate and HES record if HES record after indexdate: 
tar_target(
  name = severity_timing_median, 
  command = cohort_wide_final_2 |> filter(hes_linkage_eligible == TRUE & exposed == 1 & severe_flag_composite == "hes post-index") |>  
    # use HES date nearest to indexdate: 
    mutate(hes_index_diff_min = pmin(hes_index_diff, hes_index_diff_sepsis,
                                     na.rm = TRUE)) |> 
    summarize(median(hes_index_diff_min), 
              IQR (hes_index_diff_min)), 
  pattern = map(cohort_wide_final_2),
  iteration = "list",
  format = "parquet"
  
),

# see if there are any outcome events occuring before the HES records that are post index date: ie inbetween indexdate and HES record date
tar_target(
  name = severity_timing_outcomes, 
  command = cohort_wide_final_2 |> filter(hes_linkage_eligible == TRUE & exposed == 1 & severe_flag_composite == "hes post-index") |>  
    mutate( # first make the hes_index_diff & hes_index_diff_sepsis vars absolute (so can keep min without getting massive -ve value)
            hes_index_diff = abs(hes_index_diff), 
            hes_index_diff_sepsis = abs(hes_index_diff_sepsis),
            # use HES date nearest to indexdate: 
            hes_index_diff_min = pmin(hes_index_diff, hes_index_diff_sepsis,
                                     na.rm = TRUE), 
           # turn the HES record back in to a date: 
           hes_date = indexdate + hes_index_diff_min, 
           # gen a TRUE var for each outcome if outcome is before  HES record (but also make sure outcome after indexdate- ie to get at the exclusion criteria)
           suicide_sh_pre_hes = ifelse((self_harm_suicide_incident_composite_obsdate < hes_date)& (self_harm_suicide_incident_composite_obsdate > indexdate), TRUE, NA), 
           smi_pre_hes = ifelse((smi_incident_obsdate < hes_date) & (smi_incident_obsdate > indexdate), TRUE, NA),
           cmd_pre_hes = ifelse((cmd_incident_obsdate < hes_date) & (cmd_incident_obsdate > indexdate) , TRUE, NA)), 

  pattern = map(cohort_wide_final_2),
  iteration = "list",
  format = "parquet"
  
),

# Tabulate number of events pre HES record: CMD- ONLY  A FEW 

tar_target(
  name = severity_timing_outcomes_cmd, 
  command = severity_timing_outcomes |> tabyl(cmd_pre_hes), 
  
  pattern = map(severity_timing_outcomes),
  iteration = "list",
  format = "parquet"
  
),

# keep just these rows and relevant vars to check: 

tar_target(
  name = severity_timing_outcomes_cmd_dets, 
  command = severity_timing_outcomes |> filter(cmd_pre_hes == TRUE) |> 
    select(indexdate, hes_index_diff_min, hes_index_diff, hes_index_diff_sepsis, hes_date, cmd_pre_hes, cmd_incident_obsdate), 
  
  pattern = map(severity_timing_outcomes),
  iteration = "list",
  format = "parquet"
  
),

# Tabulate number of events pre HES record: SMI - NONE

tar_target(
  name = severity_timing_outcomes_smi, 
  command = severity_timing_outcomes |> tabyl(smi_pre_hes), 
  
  pattern = map(severity_timing_outcomes),
  iteration = "list",
  format = "parquet"
  
),

# Tabulate number of events pre HES record: SH - NONE

tar_target(
  name = severity_timing_outcomes_sh, 
  command = severity_timing_outcomes |> tabyl(suicide_sh_pre_hes), 
  
  pattern = map(severity_timing_outcomes),
  iteration = "list",
  format = "parquet"
  
),

# Update enddate for SUicide/SH - latest ONS linkage date is 3/04/2024 so end FU there (rather than study enddate of 15/06/2024 used for CMD and SMI) - as missing outcome data (suicide from ONS) after this date: 
tar_target( 
  name = cohort_wide_suicide_sh, 
  command = cohort_wide_final_2 |> left_join(denom_mx %>% select(patid, regenddate, lcd), by = "patid") |> 
    # update enddate var to bring study enddate forward to Apr 2024 (it is June 2024 in main cohort): 
    mutate(enddate = pmin(enddate, as.Date("2024-04-03"), na.rm = TRUE), 
    # gen reason FU ends var: 
           enddate_reason = case_when(
             enddate == ddate ~ "death",
             enddate == as.Date("2024-04-03") ~ "study end",
             enddate == regenddate ~ "end CPRD registration", 
             enddate == lcd ~ "end of CPRD collection from practice"), 
           # those with enddate that doesnt match any of the above dates- it is date of infection and becoming exposed
           enddate_reason = replace_na(enddate_reason, "moved to exposed"),
           enddate_reason = factor(enddate_reason)), 
  pattern = map(cohort_wide_final_2),
  iteration = "list",
  format = "parquet"
  
),

tar_target(
  name = cohort_wide_suicide_sh_n, 
  command = cohort_wide_suicide_sh %>%  tabyl(exposed)|> 
    adorn_totals("row"), 
  pattern = map(cohort_wide_suicide_sh),
  iteration = "list",
  format = "parquet"
  
),

# check numbers of matches per setid: 

tar_target(
  name = cohort_wide_suicide_sh_setn, 
  command = cohort_wide_suicide_sh %>%  group_by(setid) |> 
    mutate(exposed_n = sum(exposed == 1),
           unexposed_n = sum(exposed == 0)) |>  
    ungroup() |> 
    tabyl(exposed_n, unexposed_n), 
  # # keep only sets with atleast 1 exposed and 1 unexposed person: 
  # filter(unexposed_n >= 1 & exposed_n >=1) |> 
  # select(-c(unexposed_n, exposed_n)),  
  pattern = map(cohort_wide_suicide_sh),
  iteration = "list",
  format = "parquet"
  
),

# ## CMD: exclude those with CMD or SMI before indexdate 
## Adding flag of whether or not CMD event was 5 years prior to indexdate (ie. not looking back to "ever" as done in code above for main cohort)
tar_target(
  name = outcomes_pre_index_var_cmd_5yr, 
  command = create_obs_vars_5yrspreindex(cohort_matched = cohort_wide_final_2, 
                                         eventdata = outcome_eventdata_cmd, 
                                         codelists = codelists_outcomes_cmd),
  pattern = map(cohort_wide_final_2),
  iteration = "list",
  format = "parquet",
  deployment = "main"
), 
# then need to generate CMD 5 yrs pre indexdate variable- Y/N if any of anxiety of depression incident/historic and TRUE
# also need an incident/historic flag for CMD- if have both for one day- historic trumps incident
tar_target(
  name = outcomes_cmd_first_event_5yr, 
        command = outcomes_pre_index_var_cmd_5yr |> 
          mutate(cmd_preindex_5yrs = ifelse((anxiety_incident_5yrs_preindex == TRUE |anxiety_historic_5yrs_preindex == TRUE | depression_incident_5yrs_preindex == TRUE | depression_historic_5yrs_preindex == TRUE ), TRUE, FALSE), 
                 cmd_incident_obsdate_5yr = pmin(anxiety_incident_5yrs_preindex_obsdate, depression_incident_5yrs_preindex_obsdate, na.rm = TRUE),
                 cmd_historic_obsdate_5yr = pmin(anxiety_historic_5yrs_preindex_obsdate, depression_historic_5yrs_preindex_obsdate, na.rm = TRUE)),
  pattern = map(outcomes_pre_index_var_cmd_5yr),
  iteration = "list",
  format = "parquet",
  deployment = "main"
), 

# then bind on to cohort_wide_final_2
tar_target(
  name = cohort_wide_cmd5yr, 
  command = cohort_wide_final_2 |> bind_cols(outcomes_cmd_first_event_5yr), 
  pattern = map(cohort_wide_final_2, outcomes_cmd_first_event_5yr),
  iteration = "list",
  format = "parquet", 
  deployment = "main"
  
), 

# remove those with CMD within 5 years pre-infection or SMI ever before infection
# CHANIGNG OBS DATE TO BE THE DATE GENERATED ON UPTO 5 YEARS PRE-index. if want to change back to obsdate = cmd_incident_obsdate_5yr change here and also filter on ever cmd event not just 5 yrs pre-indeex
tar_target(
  name = cohort_wide_cmd_incident, 
  command = cohort_wide_cmd5yr %>%  filter((cmd_preindex_5yrs == FALSE | is.na(cmd_preindex_5yrs)) & (smi_preindex == FALSE | is.na(smi_preindex))) %>% 
    #select(-c(cmd_preindex, smi_preindex)) %>% 
    mutate(obsdate = cmd_incident_obsdate_5yr, 
           inf = inf_names), 
  pattern = map(cohort_wide_cmd5yr, 
                cross(inf_names)),
  iteration = "list",
  format = "parquet"
  
),



# need to drop sets that don't have at least one exposed and one unexposed- they come out of model anyway but dont want in table 1/cohort counts: 
tar_target(
  name = cohort_wide_cmd, 
  command = cohort_wide_cmd_incident %>%  group_by(setid) |> 
                                          mutate(exposed_n = sum(exposed == 1),
                                                 unexposed_n = sum(exposed == 0)) |>  
                                          ungroup() |>  
                                          # # keep only sets with atleast 1 exposed and 1 unexposed person: 
                                          filter(unexposed_n >= 1 & exposed_n >=1) |>
                                          select(-c(unexposed_n, exposed_n)) |> 
                                          # and combine the missing ethnicity variables in ethnicity_5_na: 
                                          mutate(ethnicity_5_na = ifelse(ethnicity_5_na == "5", NA, as.character(ethnicity_5_na)),
                                                 ethnicity_5_na = factor(ethnicity_5_na)),
  pattern = map(cohort_wide_cmd_incident),
  iteration = "list",
  format = "parquet"
  
),

tar_target(
  name = cohort_wide_cmd_n, 
  command = cohort_wide_cmd %>%  tabyl(exposed) |> 
    adorn_totals("row"), 
  pattern = map(cohort_wide_cmd),
  iteration = "list",
  format = "parquet"
  
),

# number of ppl per setid: 
tar_target(
  name = cohort_wide_cmd_set_num, 
  command = cohort_wide_cmd |>  group_by(setid) |> 
    mutate(
      ppl_per_set = n(),              # Number of rows per setid
      is_first_row = row_number() == 1  # Flag TRUE for the first row within each setid
    ) |>  ungroup() |> 
          filter(is_first_row == TRUE) |> 
    tabyl(ppl_per_set), 
  pattern = map(cohort_wide_cmd),
  iteration = "list",
  format = "parquet"
  
),


## SMI: exclude those with SMI before indexdate 
tar_target(
  name = cohort_wide_smi, 
  command = cohort_wide_final_2 %>%  filter(smi_preindex == FALSE | is.na(smi_preindex)) %>% 
    select(-smi_preindex) |> 
    mutate(obsdate = smi_incident_obsdate, 
           inf = inf_names), 
  pattern = map(cohort_wide_final_2, 
                cross(inf_names)),
    iteration = "list",
  format = "parquet"
  
), 

# remove matched sets with out exposed and unexposed: 
tar_target(
  name = cohort_wide_smi_comp, 
  command = cohort_wide_smi %>%  group_by(setid) |> 
    mutate(exposed_n = sum(exposed == 1),
           unexposed_n = sum(exposed == 0)) |>  
           ungroup() |>  
  # # keep only sets with atleast 1 exposed and 1 unexposed person: 
  filter(unexposed_n >= 1 & exposed_n >=1) |>
    select(-c(unexposed_n, exposed_n)),
  pattern = map(cohort_wide_smi),
  iteration = "list",
  format = "parquet"
  
),


tar_target(
  name = cohort_wide_smi_n, 
  command = cohort_wide_smi_comp %>%  tabyl(exposed) |> 
                                 adorn_totals("row"), 
  pattern = map(cohort_wide_smi_comp),
  iteration = "list",
  format = "parquet"
  
),



# Entire cohort size: 
tar_target(
  name = complete_cohort_numbers, 
  command = cohort_wide_final_2 %>%  tabyl(exposed) |> 
    adorn_totals("row"), 
  pattern = map(cohort_wide_final_2),
  iteration = "list",
  format = "parquet"
  
),

tar_target(
  name = complete_cohort_numbers_severe, 
  command = cohort_wide_final_2 %>%  filter(hes_linkage_eligible == TRUE) |> 
                                   tabyl(severe_flag_composite) |> 
                                    adorn_totals("row"), 
  pattern = map(cohort_wide_final_2),
  iteration = "list",
  format = "parquet"
  
),

# tabulate those with censoring event pre-indexdate

# dementia pre-indexdate:
tar_target(
  name = cmd_censoring_nums_dementia, 
  command = cohort_wide_cmd %>%  mutate(dementia_pre_index = ifelse(dementia_obsdate < indexdate, TRUE, NA)), 
  pattern = map(cohort_wide_cmd),
  iteration = "list",
  format = "parquet"
  
), 

# OCD pre-indexdate:
tar_target(
  name = cmd_censoring_nums_ocd, 
  command = cohort_wide_cmd %>%  mutate(ocd_pre_index = ifelse(ocd_obsdate < indexdate, TRUE, NA)) |> 
              select(ocd_pre_index) |> 
              tabyl(ocd_pre_index) |> 
              adorn_totals("row"), 
  pattern = map(cohort_wide_cmd),
  iteration = "list",
  format = "parquet"
  
),

# PTSD pre-indexdate:
tar_target(
  name = cmd_censoring_nums_ptsd, 
  command = cohort_wide_cmd %>%  mutate(ptsd_pre_index = ifelse(ptsd_obsdate < indexdate, TRUE, NA)) |> 
    select(ptsd_pre_index) |> 
    tabyl(ptsd_pre_index) |> 
    adorn_totals("row"), 
  pattern = map(cohort_wide_cmd),
  iteration = "list",
  format = "parquet"
  
),

#  NO dementia pre-indexdate but dementia pre-CMD obsdate:
tar_target(
  name = cmd_censoring_nums_dementia2, 
  command = cohort_wide_cmd %>%  mutate(dementia_pre_index = ifelse(dementia_obsdate < indexdate, TRUE, NA), 
                                        dementia_pre_CMD = ifelse(dementia_obsdate < obsdate, TRUE, NA))  |>
              select(dementia_pre_index, dementia_pre_CMD, dementia_obsdate, indexdate, obsdate),
  pattern = map(cohort_wide_cmd),
  iteration = "list",
  format = "parquet"
  
), 

# ┠ HRs: Main  Cohort--------------------------------------------------
tar_target(cohort_labels, 
           c("GE", 
             "LRTI", 
             "Mening", 
             "Sepsis", 
             "SSTI", 
             "UTI")
), 

# CMD analyses: 
tar_target(cohort_post_exclusion, 
           create_cohort_post_exclusion(cohort_wide = cohort_wide_cmd), 
           pattern = map(cohort_wide_cmd),
           iteration = "list",
           format = "parquet"
), 


# FU time by exposure and WITH all censoring vars: numbers to present in Table 1- with censoring
tar_target(censored_pop_fu_time_summ, 
           command = cohort_post_exclusion |> mutate(fup = (as.numeric(enddate) - as.numeric(indexdate))/365.25) |> 
             group_by(exposed) |> 
             summarize(median(enddate), 
                       pyears = sum(fup), 
                       fup_iqr = paste0(
                         round(quantile(fup, 0.25, na.rm = TRUE), digits = 1), "–",
                         round(quantile(fup, 0.75, na.rm = TRUE), digits = 1)), 
                       min(fup), 
                       mean(fup), 
                       median(fup), 
                       IQR(fup), 
                       max(fup)) |> 
             mutate(infection = cohort_labels),
           pattern = cross(map(cohort_post_exclusion,
                               cross(cohort_labels))),
           iteration = "list",
           format = "parquet"
), 

# rates: 
tar_target(
  name = results_rates,
  command = analysis_rates(outcome = "cmd", 
                 cohort_post_exclusion = cohort_post_exclusion, 
                 cohort_labels = cohort_labels,
                 model = model, 
                 exposure ="exposed"),
  pattern = cross(map(cohort_post_exclusion,
                      cross(cohort_labels)), 
                  model)	  
), 

# HRS: 
tar_target(
  name = results_hrs,
  command = analysis_hrs(outcome = "cmd", 
               cohort_post_exclusion = cohort_post_exclusion,
               cohort_labels = cohort_labels,
               exposure = "exposed",
               model = model),
  pattern = cross(map(cohort_post_exclusion,
                      cross(cohort_labels)),
                  
                  
                  model)
),


# combine HR and rate results: 
tar_target(
  results,
  results_hrs |> 
    left_join(results_rates, by=c("outcome", "model", "cohort")) |> 
    rename(term=term.x)
),


## SENSITIVITY ANALYSES: 

# restricting to those with >=1 GP consultation in 1 year pre-indexdate: 
tar_target(
  name = data_cons_preindex_check,
  command = cohort_post_exclusion |> filter(cons_in_year_pre_index == FALSE) |> 
                                     summarise(exposed1_n = sum(exposed ==1), 
                                               exposed0_n = sum(exposed ==0), 
                                               cmd= sum(!is.na(obsdate))),
  pattern = map(cohort_post_exclusion),
  iteration = "list",
  format = "parquet"
),

tar_target(
  name = data_cons_preindex,
  command = cohort_post_exclusion |> filter(cons_in_year_pre_index == TRUE),
  pattern = map(cohort_post_exclusion),
  iteration = "list",
  format = "parquet"
),

tar_target(
  name = results_hrs_cons_preindex, 
  command = analysis_hrs(outcome = "cmd", 
                         cohort_post_exclusion = data_cons_preindex,
                         cohort_labels = cohort_labels,
                         exposure = "exposed",
                         model = model_fully_adj),
  pattern = cross(map(data_cons_preindex,
                      cross(cohort_labels)),
                  model_fully_adj)
),

# rates to get N expsoed and CMD among exposed for supp infor table: 
tar_target(
  name = results_rates_cons_preindex,
  command = analysis_rates(outcome = "cmd", 
                           cohort_post_exclusion = data_cons_preindex, 
                           cohort_labels = cohort_labels,
                           model = model_fully_adj, 
                           exposure ="exposed"),
  pattern = cross(map(data_cons_preindex,
                      cross(cohort_labels)), 
                  model_fully_adj)	  
), 



# Not censoring at alternative mental illness diagnoses: just running on final fully adjusted model
tar_target(cohort_post_exclusion_no_censor, 
           create_cohort_post_exclusion_no_cmd_censor(cohort_wide = cohort_wide_cmd), 
           pattern = map(cohort_wide_cmd),
           iteration = "list",
           format = "parquet"
), 

tar_target(
  name = results_hrs_no_censor,
  command = analysis_hrs(outcome = "cmd", 
                         cohort_post_exclusion = cohort_post_exclusion_no_censor,
                         cohort_labels = cohort_labels,
                         exposure = "exposed",
                         model = model_fully_adj),
  pattern = cross(map(cohort_post_exclusion_no_censor,
                      cross(cohort_labels)),
                  model_fully_adj)
),

# rates to get N expsoed and CMD among exposed for supp infor table: 
tar_target(
  name = results_rates_no_censor,
  command = analysis_rates(outcome = "cmd", 
                           cohort_post_exclusion = cohort_post_exclusion_no_censor, 
                           cohort_labels = cohort_labels,
                           model = model_fully_adj, 
                           exposure ="exposed"),
  pattern = cross(map(cohort_post_exclusion_no_censor,
                      cross(cohort_labels)), 
                  model_fully_adj)	  
), 

# Main analysis but starting FU one year post-indexdate and dropping those with CMD in this 1st years post-indexdate: 

# update indexdate and drop those with CMD pre-indexdate or enddate now pre-indexdate

tar_target(
  name = data_indexdate_plus_yr, 
  command = cohort_wide_cmd |>  mutate(indexdate = indexdate + 365) |> 
                                filter(enddate > indexdate) |>  
                                filter(is.na(obsdate) | obsdate >= indexdate),
  pattern = map(cohort_wide_cmd),
  iteration = "list",
  format = "parquet"
),


tar_target(cohort_post_exclusion_indexdate_plus_yr, 
           create_cohort_post_exclusion(cohort_wide = data_indexdate_plus_yr), 
           pattern = map(data_indexdate_plus_yr),
           iteration = "list",
           format = "parquet"
), 

tar_target(
  name = results_hrs_indexdate_plus_yr,
  command = analysis_hrs(outcome = "cmd", 
                         cohort_post_exclusion = cohort_post_exclusion_indexdate_plus_yr,
                         cohort_labels = cohort_labels,
                         exposure = "exposed",
                         model = model_fully_adj),
  pattern = cross(map(cohort_post_exclusion_indexdate_plus_yr,
                      cross(cohort_labels)),
                  model_fully_adj)
),

# rates to get N exposed and CMD among exposed for supp infor table: 
tar_target(
  name = results_rates_indexdate_plus_yr,
  command = analysis_rates(outcome = "cmd", 
                           cohort_post_exclusion = cohort_post_exclusion_indexdate_plus_yr, 
                           cohort_labels = cohort_labels,
                           model = model_fully_adj, 
                           exposure ="exposed"),
  pattern = cross(map(cohort_post_exclusion_indexdate_plus_yr,
                      cross(cohort_labels)), 
                  model_fully_adj)	  
), 

# ENDING FOLLOW UP on 1st March 2020  (pre-covid): 
tar_target(cohort_post_exclusion_march20_enddate, 
           create_cohort_post_exclusion_march20_enddate(cohort_wide = cohort_wide_cmd), 
           pattern = map(cohort_wide_cmd),
           iteration = "list",
           format = "parquet"
), 


tar_target(
  name = results_hrs_indexdate_march20_enddate,
  command = analysis_hrs(outcome = "cmd", 
                         cohort_post_exclusion = cohort_post_exclusion_march20_enddate,
                         cohort_labels = cohort_labels,
                         exposure = "exposed",
                         model = model_fully_adj),
  pattern = cross(map(cohort_post_exclusion_march20_enddate,
                      cross(cohort_labels)),
                  model_fully_adj)
),

tar_target(
  name = results_rates_march20_enddate,
  command = analysis_rates(outcome = "cmd", 
                           cohort_post_exclusion = cohort_post_exclusion_march20_enddate, 
                           cohort_labels = cohort_labels,
                           model = model_fully_adj, 
                           exposure ="exposed"),
  pattern = cross(map(cohort_post_exclusion_march20_enddate,
                      cross(cohort_labels)), 
                  model_fully_adj)	  
), 

# ADDITIONALLY ADJUSTING FOR INJECTING DRUG USE FOR SSTI COHORT: 
# Add injecting drug use to the model: 
# NB running on all infections but only need the SSTI estimate
tar_target(
  model_fully_adj_inj_drug, 
  c("C" = paste(c("", "imd",  "smoking_status_simp", "alcohol_abuse", "obese", "cci_cat", "ethnicity_5_na", "injectable_drugs"), collapse = " + "))
),

tar_target(
  name = results_hrs_ssti_inj,
  command = analysis_hrs(outcome = "cmd", 
                         cohort_post_exclusion = cohort_post_exclusion,
                         cohort_labels = cohort_labels,
                         exposure = "exposed",
                         model = model_fully_adj_inj_drug),
  pattern = cross(map(cohort_post_exclusion,
                      cross(cohort_labels)),
                  model_fully_adj_inj_drug)
),

tar_target(
  name = results_rates_ssti_inj,
  command = analysis_rates(outcome = "cmd", 
                           cohort_post_exclusion = cohort_post_exclusion, 
                           cohort_labels = cohort_labels,
                           model = model_fully_adj_inj_drug, 
                           exposure ="exposed"),
  pattern = cross(map(cohort_post_exclusion,
                      cross(cohort_labels)), 
                  model_fully_adj_inj_drug)	  
), 

# Missing category approch to ethnicity- ie using ethnicity_5 (which has NA coded as Unknown and included in model)
tar_target(
  model_fully_adj_missing_cat, 
  c("C" = paste(c("","", "imd",  "smoking_status_simp", "alcohol_abuse", "obese", "cci_cat", "ethnicity_5"), collapse = " + "))
),

tar_target(
  name = results_hrs_missing_ethnicity,
  command = analysis_hrs(outcome = "cmd", 
                         cohort_post_exclusion = cohort_post_exclusion,
                         cohort_labels = cohort_labels,
                         exposure = "exposed",
                         model = model_fully_adj_missing_cat),
  pattern = cross(map(cohort_post_exclusion,
                      cross(cohort_labels)),
                  model_fully_adj_missing_cat)
),

tar_target(
  name = results_rates_missing_ethnicity,
  command = analysis_rates(outcome = "cmd", 
                           cohort_post_exclusion = cohort_post_exclusion, 
                           cohort_labels = cohort_labels,
                           model = model_fully_adj_missing_cat, 
                           exposure ="exposed"),
  pattern = cross(map(cohort_post_exclusion,
                      cross(cohort_labels)), 
                  model_fully_adj_missing_cat)	  
), 


# Adding just CCI to model with IMD
tar_target(
  model_imd_cci, 
  c("C" = paste(c("", "imd", "cci_cat"), collapse = " + "))
),

# HRS: 
tar_target(
  name = results_hrs_cci,
  command = analysis_hrs(outcome = "cmd", 
                         cohort_post_exclusion = cohort_post_exclusion,
                         cohort_labels = cohort_labels,
                         exposure = "exposed",
                         model = model_imd_cci),
  pattern = cross(map(cohort_post_exclusion,
                      cross(cohort_labels)),
                  
                  
                  model_imd_cci)
),

tar_target(
  name = results_rates_cci,
  command = analysis_rates(outcome = "cmd", 
                           cohort_post_exclusion = cohort_post_exclusion, 
                           cohort_labels = cohort_labels,
                           model = model_imd_cci, 
                           exposure ="exposed"),
  pattern = cross(map(cohort_post_exclusion,
                      cross(cohort_labels)), 
                  model_imd_cci)	  
), 


# HRs startified by ANTIMICROBIAL RPESCRIPTION: 
# conduct separate models for exposed + antimicrobial presc. and exposed w/out antimicro presc and compare HRs: 
# keep all unexposed and only exposed: 
tar_target(
  name = cohort_antimicro_true, 
  command = cohort_post_exclusion |> group_by(setid) |> 
                            filter(any(antimicro_presc == TRUE, na.rm = TRUE)) %>%
                            ungroup(), 
  pattern = map(cohort_post_exclusion),
  iteration = "list",
  format = "parquet"
),  





tar_target(
  name = results_hrs_antimicro_true,
  command = analysis_hrs(outcome = "cmd", 
                         cohort_post_exclusion = (cohort_antimicro_true),
                         cohort_labels = cohort_labels,
                         exposure = "exposed",
                         model = model_fully_adj),
  pattern = cross(map(cohort_antimicro_true,
                      cross(cohort_labels)),
                  model_fully_adj)
),

# exposed NOT with record of antimicrobial 
tar_target(
  name = cohort_antimicro_false, 
  command = cohort_post_exclusion |> group_by(setid) |> 
    filter(any(antimicro_presc == FALSE, na.rm = TRUE)) %>%
    ungroup(), 
  pattern = map(cohort_post_exclusion),
  iteration = "list",
  format = "parquet"
),  


tar_target(
  name = results_hrs_antimicro_false,
  command = analysis_hrs(outcome = "cmd", 
                         cohort_post_exclusion = (cohort_antimicro_false),
                         cohort_labels = cohort_labels,
                         exposure = "exposed",
                         model = model_fully_adj),
  pattern = cross(map(cohort_antimicro_false,
                      cross(cohort_labels)),
                  model_fully_adj)
),

# Tabulate severe infections and antimicroial prescriptions- 
tar_target(
  name = severity_antimicro_tab, 
  command = cohort_wide_cmd |>  filter(!is.na(severe_flag_composite) & !is.na(antimicro_presc)) |>   # remove NAs so not included in the percentages
                                tabyl(severe_flag_composite, antimicro_presc)  |> 
                                    adorn_percentages("row")  |> 
                                    adorn_totals("row")  |> 
                                    adorn_pct_formatting(), 
  pattern = map(cohort_wide_cmd),
  iteration = "list",
  format = "parquet"
),  


# SECONDARY ANALYSES: 
# Association persists/wanes over time: 
# generate first the enddate for each of these sep analyses (0-6mon, 0-1yr, 0-2,-3,-4,-5)
tar_target(
  name = cohort_wide_cmd_fu_dates, 
  command = cohort_wide_cmd |> mutate(enddate_6mon = indexdate + 182.5, 
                                      enddate_1yr = indexdate + 365, 
                                      enddate_2yr = indexdate + (365*2),
                                      enddate_3yr = indexdate + (365*3),
                                      enddate_4yr = indexdate + (365*4),
                                      enddate_5yr = indexdate + (365*5)),
    pattern = map(cohort_wide_cmd),
    iteration = "list",
    format = "parquet"
  ),

# create the post exclusion the relant function and each enddate var

# 0-6 months
tar_target(cohort_post_exclusion_6mon, 
           create_cohort_post_exclusion_6mon(cohort_wide = cohort_wide_cmd_fu_dates), 
           pattern = map(cohort_wide_cmd_fu_dates),
           iteration = "list",
           format = "parquet"
), 


tar_target(
  name = results_hrs_6mon, 
  command = analysis_hrs(outcome = "cmd", 
                         cohort_post_exclusion = cohort_post_exclusion_6mon,
                         cohort_labels = cohort_labels,
                         exposure = "exposed",
                         model = model_fully_adj),
  pattern = cross(map(cohort_post_exclusion_6mon,
                      cross(cohort_labels)),
                  model_fully_adj)
),

# 0-1 year 
tar_target(cohort_post_exclusion_1yr, 
           create_cohort_post_exclusion_1yr(cohort_wide = cohort_wide_cmd_fu_dates), 
           pattern = map(cohort_wide_cmd_fu_dates),
           iteration = "list",
           format = "parquet"
), 



tar_target(
  name = results_hrs_1yr, 
  command = analysis_hrs(outcome = "cmd", 
                         cohort_post_exclusion = cohort_post_exclusion_1yr,
                         cohort_labels = cohort_labels,
                         exposure = "exposed",
                         model = model_fully_adj),
  pattern = cross(map(cohort_post_exclusion_1yr,
                      cross(cohort_labels)),
                  model_fully_adj)
),

# 0-2 year 
tar_target(cohort_post_exclusion_2yr, 
           create_cohort_post_exclusion_2yr(cohort_wide = cohort_wide_cmd_fu_dates), 
           pattern = map(cohort_wide_cmd_fu_dates),
           iteration = "list",
           format = "parquet"
), 


tar_target(
  name = results_hrs_2yr, 
  command = analysis_hrs(outcome = "cmd", 
                         cohort_post_exclusion = cohort_post_exclusion_2yr,
                         cohort_labels = cohort_labels,
                         exposure = "exposed",
                         model = model_fully_adj),
  pattern = cross(map(cohort_post_exclusion_2yr,
                      cross(cohort_labels)),
                  model_fully_adj)
),

# 0-3 year 
tar_target(cohort_post_exclusion_3yr, 
           create_cohort_post_exclusion_3yr(cohort_wide = cohort_wide_cmd_fu_dates), 
           pattern = map(cohort_wide_cmd_fu_dates),
           iteration = "list",
           format = "parquet"
), 

tar_target(
  name = results_hrs_3yr, 
  command = analysis_hrs(outcome = "cmd", 
                         cohort_post_exclusion = cohort_post_exclusion_3yr,
                         cohort_labels = cohort_labels,
                         exposure = "exposed",
                         model = model_fully_adj),
  pattern = cross(map(cohort_post_exclusion_3yr,
                      cross(cohort_labels)),
                  model_fully_adj)
),

# 0-4 year 
tar_target(cohort_post_exclusion_4yr, 
           create_cohort_post_exclusion_4yr(cohort_wide = cohort_wide_cmd_fu_dates), 
           pattern = map(cohort_wide_cmd_fu_dates),
           iteration = "list",
           format = "parquet"
), 


tar_target(
  name = results_hrs_4yr, 
  command = analysis_hrs(outcome = "cmd", 
                         cohort_post_exclusion = cohort_post_exclusion_4yr,
                         cohort_labels = cohort_labels,
                         exposure = "exposed",
                         model = model_fully_adj),
  pattern = cross(map(cohort_post_exclusion_4yr,
                      cross(cohort_labels)),
                  model_fully_adj)
),

# 0-5 year 
tar_target(cohort_post_exclusion_5yr, 
           create_cohort_post_exclusion_5yr(cohort_wide = cohort_wide_cmd_fu_dates), 
           pattern = map(cohort_wide_cmd_fu_dates),
           iteration = "list",
           format = "parquet"
), 


tar_target(
  name = results_hrs_5yr, 
  command = analysis_hrs(outcome = "cmd", 
                         cohort_post_exclusion = cohort_post_exclusion_5yr,
                         cohort_labels = cohort_labels,
                         exposure = "exposed",
                         model = model_fully_adj),
  pattern = cross(map(cohort_post_exclusion_5yr,
                      cross(cohort_labels)),
                  model_fully_adj)
),

# join all FU period results together and format as need for table: 
tar_target(
  name = results_hrs_all_fu, 
  command = bind_rows(
    results_hrs_6mon |> mutate(period = "6mon"), 
    results_hrs_1yr |> mutate(period = "1yr"), 
    results_hrs_2yr |> mutate(period = "2yr"),
    results_hrs_3yr |> mutate(period = "3yr"),
    results_hrs_4yr |> mutate(period = "4yr"),
    results_hrs_5yr |> mutate(period = "5yr")) |> 
    filter(term == "exposed") |> 
    mutate(result = paste0(round(estimate, digits = 2), " (", round(conf.low, digits = 2), "-", round(conf.high, digits = 2), ")")) |> 
    select(result, cohort, period) |> 
    pivot_wider(names_from = cohort,
                values_from = result) |> 
    relocate(period, SSTI, UTI, GE, LRTI, Sepsis, Mening)
  
    
),

# SEPERATNG INFECTED IN TO SEVERE VS NON-SEVERE: 
# Remove those not eligible for HES linkage (cant know if they are severe or not severe):  
# simplify the severity variable-
tar_target(
  name = data_severity, 
  command = cohort_wide_cmd |>  filter(hes_linkage_eligible == TRUE) |> 
                                mutate( severity =  case_when(
                                  severe_flag_composite == "not severe" ~ "Non-severe",
                                  severe_flag_composite == "hes pre-index" ~ "Severe",   
                                  severe_flag_composite == "hes post-index" ~ "Severe",   
                                  is.na(severe_flag_composite) ~ "Uninfected"), 
                                  # set uninfected as first group so that it is the reference group in the COx regression: 
                                  severity = factor(severity, levels = c("Uninfected", "Non-severe", "Severe"))),
  pattern = map(cohort_wide_cmd),
  iteration = "list",
  format = "parquet"
),


tar_target(cohort_post_exclusion_severity, 
           create_cohort_post_exclusion(cohort_wide = data_severity), 
           pattern = map(data_severity),
           iteration = "list",
           format = "parquet"
), 


tar_target(
  name = results_hrs_severity,
  command = analysis_hrs(outcome = "cmd", 
                         cohort_post_exclusion = cohort_post_exclusion_severity,
                         cohort_labels = cohort_labels,
                         exposure = "severity",
                         model = model_fully_adj),
  pattern = cross(map(cohort_post_exclusion_severity,
                      cross(cohort_labels)),
                  model_fully_adj)
), 

tar_target(
  name = results_rates_severity,
  command = analysis_rates(outcome = "cmd", 
                           cohort_post_exclusion = cohort_post_exclusion_severity, 
                           cohort_labels = cohort_labels,
                           model = model_fully_adj, 
                           exposure ="severity"),
  pattern = cross(map(cohort_post_exclusion_severity,
                      cross(cohort_labels)), 
                  model_fully_adj)	  
), 

# # combine severity HR and rate results: 
tar_target(
  results_severity,
  results_hrs_severity |> 
    left_join(results_rates_severity, by=c("outcome", "model", "cohort")) |> 
    rename(term=term.x)
),



# tabulate number severe/non-severe/uninfected among linked cohort: 
tar_target(
  name = data_severity_n, 
  command = data_severity |> tabyl(severity), 
  pattern = map(data_severity),
  iteration = "list",
  format = "parquet"
),

# tabulate number eligible for HES linkage by exposure status: 
tar_target(
  name = cohort_wide_hes_n, 
  command = cohort_wide_final |> tabyl(exposed, hes_linkage_eligible), 
  pattern = map(cohort_wide_final),
  iteration = "list",
  format = "parquet"
),


# ANXIETY AS OUTCOME: 
tar_target(cohort_post_exclusion_anxiety, 
           create_cohort_post_exclusion_anxiety(cohort_wide = cohort_wide_cmd), 
                                                pattern = map(cohort_wide_cmd),
                                                iteration = "list",
                                                format = "parquet"
), 
           

tar_target(
  name = results_hrs_anxiety,
  command = analysis_hrs(outcome = "anxiety", 
                         cohort_post_exclusion = cohort_post_exclusion_anxiety,
                         cohort_labels = cohort_labels,
                         exposure = "exposed",
                         model = model_fully_adj),
  pattern = cross(map(cohort_post_exclusion_anxiety,
                      cross(cohort_labels)),
                  model_fully_adj)
),

# DEPRESSION AS OUTCOME: 
tar_target(cohort_post_exclusion_depression, 
           create_cohort_post_exclusion_depression(cohort_wide = cohort_wide_cmd), 
           pattern = map(cohort_wide_cmd),
           iteration = "list",
           format = "parquet"
), 


tar_target(
  name = results_hrs_depression,
  command = analysis_hrs(outcome = "depression", 
                         cohort_post_exclusion = cohort_post_exclusion_depression,
                         cohort_labels = cohort_labels,
                         exposure = "exposed",
                         model = model_fully_adj),
  pattern = cross(map(cohort_post_exclusion_depression,
                      cross(cohort_labels)),
                  model_fully_adj)
),



# PTSD AS OUTCOME: 
tar_target(cohort_post_exclusion_ptsd, 
           create_cohort_post_exclusion_ptsd(cohort_wide = cohort_wide_cmd), 
           pattern = map(cohort_wide_cmd),
           iteration = "list",
           format = "parquet"
), 

tar_target(
  name = results_hrs_ptsd,
  command = analysis_hrs(outcome = "depression", 
                         cohort_post_exclusion = cohort_post_exclusion_ptsd,
                         cohort_labels = cohort_labels,
                         exposure = "exposed",
                         model = model_fully_adj),
  pattern = cross(map(cohort_post_exclusion_ptsd,
                      cross(cohort_labels)),
                  model_fully_adj)
),

## EFEECT MODIFICATION: 
# Age as effect modifier: 
# NB using "exposed:indexdate_age_cat" as gives the HRs for each stratum (* gives the interaction terms)

tar_target(
  name = results_hrs_effect_mod_age, 
  command = analysis_hrs(outcome = "cmd", 
                         cohort_post_exclusion = cohort_post_exclusion,
                         cohort_labels = cohort_labels,
                         exposure = "exposed",
                         model =  ":indexdate_age_cat + ethnicity_5_na + smoking_status_simp + alcohol_abuse + obese + cci_cat + imd"),
  pattern = cross(map(cohort_post_exclusion,
                      cross(cohort_labels)))
),

# Sex as effect modifier: 

tar_target(
  name = results_hrs_effect_mod_sex, 
  command = analysis_hrs(outcome = "cmd", 
                                        cohort_post_exclusion = cohort_post_exclusion,
                                        cohort_labels = cohort_labels,
                                        exposure = "exposed",
                                        model =  ":gender + ethnicity_5_na + smoking_status_simp + alcohol_abuse + obese + cci_cat + imd"),
  pattern = cross(map(cohort_post_exclusion,
                      cross(cohort_labels)))
),


# Frailty as an effect modifier: 
# first filter to >65 yr 
tar_target(
  name = cohort_wide_cmd_65, 
  command = cohort_wide_cmd |> filter(indexdate_age >= 65 & indexdate_age <= 95), 
  pattern = map(cohort_wide_cmd),
  iteration = "list",
  format = "parquet"
), 


tar_target(cohort_post_exclusion_65, 
           create_cohort_post_exclusion(cohort_wide = cohort_wide_cmd_65), 
           pattern = map(cohort_wide_cmd_65),
           iteration = "list",
           format = "parquet"
), 


tar_target(
  name = results_hrs_effect_frailty, 
  command = analysis_hrs(outcome = "cmd", 
                                    cohort_post_exclusion = cohort_post_exclusion_65,
                                    cohort_labels = cohort_labels,
                                    exposure = "exposed",
                                    model =  ":frailty_score + ethnicity_5_na + smoking_status_simp + alcohol_abuse + obese + cci_cat + imd"),
  pattern = cross(map(cohort_post_exclusion_65,
                      cross(cohort_labels)))
),

tar_target(
  name = results_hrs_effect_frailty2, 
  command = analysis_hrs(outcome = "cmd", 
                         cohort_post_exclusion = cohort_post_exclusion_65,
                         cohort_labels = cohort_labels,
                         exposure = "exposed",
                         model =  ":frailty_score + frailty_score + ethnicity_5_na + smoking_status_simp + alcohol_abuse + obese + cci_cat + imd"),
  pattern = cross(map(cohort_post_exclusion_65,
                      cross(cohort_labels)))
),


tar_target(
  name = results_hrs_effect_frailty4, 
  command = analysis_hrs(outcome = "cmd", 
                         cohort_post_exclusion = cohort_post_exclusion_65,
                         cohort_labels = cohort_labels,
                         exposure = "exposed",
                         model =  "+ frailty_score + frailty_score:exposed + ethnicity_5_na + smoking_status_simp + alcohol_abuse + obese + cci_cat + imd"),
  pattern = cross(map(cohort_post_exclusion_65,
                      cross(cohort_labels)))
),


# append  together the effect mod results, tidy
tar_target(
  name = effect_mod_results, 
  command = rbind(results_hrs_effect_mod_age, results_hrs_effect_mod_sex, results_hrs_effect_frailty) |> 
            filter(str_starts(term, "exposed")) |> 
            mutate(result = paste0(round(estimate, digits = 2), " (", round(conf.low, digits = 2), "-", round(conf.high, digits = 2), ")")) |> 
            select(term, cohort, result) |> 
            spread(cohort, result) |> 
            relocate(term, SSTI, UTI, GE, LRTI, Sepsis, Mening), 
  format = "parquet"
), 


######################################
# SAVING DATA: 
# save model results for plots/tables: 
tar_target(
  results_hrs_save, 
  command = save_and_return_path_rds(path = paste0(outputs_path, "/cmd_cohort/"), 
                                     x=results_hrs)
),

# rates and HR combined: 
tar_target(
  results_save, 
  command = save_and_return_path_rds(path = paste0(outputs_path, "/cmd_cohort/"), 
                                     x=results)
),


# SMI ANALYSES: 
tar_target(cohort_post_exclusion_smi, 
           create_cohort_post_exclusion_smi(cohort_wide = cohort_wide_smi_comp), 
           pattern = map(cohort_wide_smi_comp),
           iteration = "list",
           format = "parquet"
), 


# rates: 
tar_target(
  results_rates_smi,
  analysis_rates(outcome = "smi", 
                 cohort_post_exclusion = cohort_post_exclusion_smi, 
                 cohort_labels = cohort_labels,
                 model = model_fully_adj, 
                 exposure ="exposed"),
  pattern = cross(map(cohort_post_exclusion_smi,
                      cross(cohort_labels)))	  
), 

 
# SMI hazard ratios: 
tar_target(
  results_hrs_smi,
  analysis_hrs(outcome = "smi", 
               cohort_post_exclusion = cohort_post_exclusion_smi, 
               cohort_labels = cohort_labels, 
               exposure = "exposed", 
               model = model_fully_adj),
  pattern = cross(map(cohort_post_exclusion_smi,
                      cross(cohort_labels))), 
                  	  
),


# save the data for tables/graphs: 

# save as parquet files- for use in R - one file per infection 
# CMD cohort
tar_target(
  name = write_parquet_cleaned_cmd,
  command= list_write_parquet_and_return_path(x = cohort_wide_cmd, 
                                              path = paste0(outputs_path, "/cmd_cohort/"), 
                                              infection = inf_names),
  pattern = map(cohort_wide_cmd, inf_names)
), 

# CMD severity results: save to be plotted: save as one dataset: 
tar_target(
  name = results_hrs_save_severity,
  command= save_and_return_path_rds(x = results_hrs_severity, 
                                              path = paste0(outputs_path, "/cmd_cohort/severity_sens_anal/")) 
), 


# CMD by infection severity results: 
# rates and HR combined: 
tar_target(
  results_save_severity, 
  command = save_and_return_path_rds(path = paste0(outputs_path, "/cmd_cohort/"), 
                                     x=results_severity)
),


# SMI cohort
tar_target(
  name = write_parquet_cleaned_smi,
  command= list_write_parquet_and_return_path(x = cohort_wide_smi, 
                                              path = paste0(outputs_path, "/smi_cohort/20250902/"), 
                                              infection = inf_names),
  pattern = map(cohort_wide_smi, inf_names)
), 

# all data as csv files: ie before remove those with outcome before indexdate:
tar_target(
  name = write_csv_cleaned,
  command= write_and_return_path(x = cohort_wide_final, 
                                              path =  paste0(outputs_path, "/all_data/"), 
                                              infection = inf_names),
  pattern = map(cohort_wide_final, inf_names)
),

# and all data as parquet: 
tar_target(
  name = write_parquet_cleaned,
  command= list_write_parquet_and_return_path(x = cohort_wide_final, 
                                              path = paste0(outputs_path, "/all_data/"),  
                                              infection = inf_names),
  pattern = map(cohort_wide_final, inf_names)
), 


# all data as parquet files - 
tar_target(
  name = write_parquet_cleaned_jaime,
  command= list_write_parquet_and_return_path(x = cohort_wide_suicide_sh, 
                                              path = paste0(path_jaime, "/20251127/"),  
                                              infection = inf_names),
  pattern = map(cohort_wide_suicide_sh, inf_names)
), 



########################################################
## TYPE 2 linked data as received from CPRD:
########################################################

# READ IN THE DATAFILES: 

tar_target(
  name = linked_data_t2_files,
  command = list.files(path = linked_data_t2, pattern ="*.parquet", full.names=TRUE) 
), 

# death data: 
tar_target(
  name = hes_death, 
  command = open_dataset(linked_data_t2_files[[1]]) |> 
  
  # fix data formatting
   mutate(reg_date_of_death  = as.Date(reg_date_of_death , format = "%Y-%m-%d"), 
                 s_cod_code_1 = as.character(s_cod_code_1), 
                 s_underlying_cod_icd10 = as.character(s_underlying_cod_icd10))  |> 
  
  collect ()  

), 

# join death data to infection cohorts: 
# first remove duplicate data - multiple rows per person, sometimes with same death date, other times with different dates, keep first record when dates differ: 
# this is an issue CPRD are aware of with the Nov 2024 ONS data, they offer no solution. 
tar_target(name = hes_death_clean, 
           command = hes_death  |> 
                              group_by(patid) |> 
                              slice_min(order_by = reg_date_of_death, with_ties = FALSE) |> 
                              ungroup(), 
           format = "parquet"
), 
           
tar_target( 
  name = hes_death_dups, 
  command =  hes_death  |> 
    group_by(patid) |>
    mutate(flag_disagree = n_distinct(reg_date_of_death) > 1) |>
    ungroup() |>
    filter(flag_disagree == TRUE), 
  format = "parquet"
), 

  
  
tar_target(name = ons_death_date, 
           command = inf_parquet_correct |>  
                     left_join(hes_death_clean, by = "patid"), 
           pattern = map(inf_parquet_correct),
           format = "parquet"
), 

# check inf cohort numbers: they match
tar_target(name = ons_death_date_n,
           command = nrow(ons_death_date), 
           pattern = map(ons_death_date)
), 


# identify those with a suicide ICD-10 record - (only matching on 1st 3 digits of code, as per the codelist with the ICD_10 dictionary codes, all only XXX)
tar_target(name = ons_death_date_suicide, 
           command = ons_death_date |> mutate(suicide_icd10 = ifelse(substr(s_underlying_cod_icd10, 1, 3) %in% suicide_hes_codes | substr(s_cod_code_1, 1, 3) %in% suicide_hes_codes , TRUE, NA)) |> 
             select(reg_date_of_death, suicide_icd10) |> 
             rename(reg_date_of_death_ons = reg_date_of_death), 
           pattern = map(ons_death_date),
           format = "parquet"
),

# tabulate number of ONS recorded suicides in each infection cohort: 
tar_target(name = ons_death_date_suicide_n, 
           command = ons_death_date_suicide |> tabyl(suicide_icd10, show_na = TRUE) |> 
                                                    adorn_totals("row") |> 
                                                    adorn_percentages("row") |> 
                                                    adorn_pct_formatting(digits = 2)  |> 
                                                    adorn_ns(), 
           pattern = map(ons_death_date_suicide)
),


#┠ IMD------------------------------------

# patient IMD data:
# NB: e2019_imd_5 gives quintile (1=LEAST deprived)
tar_target(
  name = imd_patient, 
  command = open_dataset(linked_data_t2_files[[8]]) |> 
    select(-pracid) |> 
    
    collect ()  
  
), 

# practice IMD data:
# NB: e2019_imd_5 gives quintile (1=LEAST deprived)
# NB: only variable need is e2019_imd_5 - others are for NI/Scot etc
# NB: got this full practice level IMD data from a 2nd request, so reading in directly
tar_target(
  name = imd_practice, 
  command = open_dataset(paste0(linked_data_t2, "\\24_003746_practiceIMD_quintiles.parquet")) |> 

        collect ()  
  
), 

# IMD data: practice and patient level
# first have to add pracid back on to infection data

tar_target(
  name = imd,
  command = inf_parquet_correct |>  left_join(denom_mx %>% select(pracid, patid), by = "patid") |> 
                                    left_join(imd_patient, by="patid") |> 
                                    rename(imd=e2019_imd_5) |> 
                                    mutate(pracid = as.character(pracid)) |> 
                                    # add practise level data: 
                                    left_join(imd_practice, by="pracid") |> 

                                    mutate(imd=ifelse(is.na(imd), e2019_imd_5, imd)) |> 
                                    select(imd), 
  pattern = map(inf_parquet_correct),
  iteration = "list",
  format = "parquet"
), 


# HES INFECTION DATA: (to generate severity flags):
# All 6 infections: simplest to branch over all 6 infections even though dont need mening/encep and sepsis needs repeating for 4 infections: UTI/LRTI/GE/SSTI - 
# order is the same as CPRD infection branches so can cross map over: 
# icd-10 codelists: 
tar_target( # Codelists used to extract eventdata
  codelists_infections_icd10,
  tribble(~name, ~label, ~codevar, ~extract_from,
          
          "gastro", "GE", "code_icd_10", "observation",
          "lrti", "LRTI", "code_icd_10", "observation",
          "mening_enceph", "Mening/Enceph", "code_icd_10", "observation",
          "sepsis", "Sepsis", "code_icd_10", "observation", 
          "ssti", "SSTI", "code_icd_10", "observation",
          "uti", "UTI", "code_icd_10", "observation") %>%   
    
    mutate(path=paste0(codelists_path, "\\ICD_10\\", name, "_icd10.csv"), 
           # store constents of each codelist csv file as list under codelist$full: 
           full=map(path, ~read_csv(.x, col_types = cols(.default = "c"))), 
           
           # extract codevar (medcodeid/prodcodeid) and store under codelists$codes : 
           codes=map2(path, codevar, ~read_csv(.x, col_types = cols(.default = "c"))[[.y]])) 
  
), 

tar_target( # Codelists used to extract eventdata
  codelists_sepsis_icd10, 
  tribble(~name, ~label, ~codevar, ~extract_from,
          
          "sepsis", "Sepsis", "code_icd_10", "observation") %>%   
    
    mutate(path=paste0(codelists_path, "\\ICD_10\\", name, "_icd10.csv"), 
           # store constents of each codelist csv file as list under codelist$full: 
           full=map(path, ~read_csv(.x, col_types = cols(.default = "c"))), 
           
           # extract codevar (medcodeid/prodcodeid) and store under codelists$codes : 
           codes=map2(path, codevar, ~read_csv(.x, col_types = cols(.default = "c"))[[.y]])) 
  
), 


# 4 infections: GE/LRTI/SSTI/UTI
tar_target(
  name = hes_diagnosis_hosp_inf, 
  command = open_dataset(linked_data_t2_files[[3]]) |> 
    filter(ICD %in% unlist(codelists_infections_icd10$codes)) %>%
    collect(),
  pattern = map(codelists_infections_icd10),
  iteration = "list",
  format = "parquet"
), 

# Sepsis: need as separate infection to run over all 4 other infection branches: 
# 4 infections: GE/LRTI/SSTI/UTI
tar_target(
  name = hes_diagnosis_hosp_sepsis, 
  command = open_dataset(linked_data_t2_files[[3]]) |> 
    filter(ICD %in% unlist(codelists_sepsis_icd10$codes)) %>%
    collect(),
  iteration = "list",
  format = "parquet"
), 

tar_target(
  name = hes_episodes, 
  command = open_dataset(linked_data_t2_files[[4]]) |> 
    collect(),
  iteration = "list",
  format = "parquet"
), 

tar_target(
  name = hes_hospital, 
  command = open_dataset(linked_data_t2_files[[5]]) |> 
    collect(),
  iteration = "list",
  format = "parquet"
), 

tar_target(
  name = hes_patient, 
  command = open_dataset(linked_data_t2_files[[6]]) |> 
    collect(),
  iteration = "list",
  format = "parquet"
), 

tar_target(
  name = hes_primary_diag_hosp, 
  command = open_dataset(linked_data_t2_files[[7]]) |> 
    filter(ICD_PRIMARY %in% unlist(codelists_infections_icd10$codes))  |> 
    collect(),
  pattern = map(codelists_infections_icd10),
  iteration = "list",
  format = "parquet"
)

)

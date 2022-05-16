# Import library

library(tidyverse)
library(lubridate)

# Load data
path <- 
  fs::path("", "Volumes", "Gillis_Research", "Christelle Colin-Leitzinger", "CH_transfer_HCT")

clinical <- 
  readxl::read_xlsx(
    paste0(path, 
           "/data/CICPT2144_cleaned vcf_01.21.22_reviewed w EHR and raw vcf_03.30.22__MUTATIONS.xlsx"),
    sheet = "ClinicalData")

calls_CH <- readxl::read_xlsx(
  paste0(path, 
         "/data/MCCdonors_CH calls cleaned_04.29.21.xlsx"),
  sheet = "MCCdonors_CH calls cleaned_04.2")

# Date cleaning
## Clinical
clinical <- clinical %>% 
  select(-c(primedx.y, bmt_date.y), primedx = primedx.x, bmt_date = bmt_date.x) %>%
  mutate(across("rec_cmv_result", ~na_if(., "NA"))) %>% 
  mutate(across(c("bmt_date", "dod", "relapse_date", "date_last_contact"), 
                ~ as.Date(.))
         ) %>% 
  mutate(across(c("rec_cmv_result", "don_cmv_res"), 
                ~ str_to_sentence(.))
  ) %>% 
  unite("sex_pairs", c(sex, don_sex), sep = "/", remove = FALSE) %>% 
  unite("race_pairs", c(race, don_race), sep = "/", remove = FALSE) %>% 
  unite("cmv_pairs", c(rec_cmv_result, don_cmv_res), sep = "/", remove = FALSE) %>% 
  
  
  
  mutate(os_event = case_when(
  is.na(dod)        ~ 0,
  !is.na(dod)       ~ 1
  )) %>% 
  mutate(os_date = coalesce(dod, date_last_contact),
         os_time = 
           interval(
           start = bmt_date, 
           end = os_date)/
           duration(n=1, units = "days")) %>% 
  mutate(relapse_event = case_when(
    is.na(relapse_date)        ~ 0,
    !is.na(relapse_date)       ~ 1
  )) %>% 
  mutate(relapse_date = coalesce(relapse_date, date_last_contact),
         relapse_time = 
           interval(
             start = bmt_date, 
             end = relapse_date)/
           duration(n=1, units = "days"))

write_rds(clinical, "clinical.rds")


## CH calls
calls_CH <- calls_CH %>% 
  mutate(CH_status = case_when(
    !is.na(CHROM)                      ~ "CH",
    is.na(CHROM)                       ~ "No CH"
  ))

write_rds(calls_CH, "calls_CH.rds")


# End Cleaning




library(touch)
library(tidyverse)
library(icd)

diag = read.csv("DiagnosesComprehensiveAll.csv")
encounter = read.csv("EncounterAll.csv")
como = read.csv("ComorbiditiesElixhauserComprehensive.csv")
race = read.csv("PatientRace.csv")
social = read.csv("ClaritySocialHistory.csv")

diag_convert = diag

# Convert to ICD-10
diag_convert[diag_convert$Lexicon == "ICD9","TermCodeMapped"] = icd_map(diag_convert[diag_convert$Lexicon == "ICD9","TermCodeMapped"])
diag_convert$TermCodeMapped = gsub(",.*", "", diag_convert$TermCodeMapped)
diag_convert$TermCodeMapped[diag_convert$Lexicon == "ICD9"] = short_to_decimal(diag_convert$TermCodeMapped[diag_convert$Lexicon == "ICD9"])

# Match the encounter date with diagnoses
# Add a new column "PE" and set it to 1
diag_date = diag_convert %>%
    select(-TermNameMapped, -Lexicon) %>%
    left_join(encounter, by = c("DeID_PatientID", "DeID_EncounterID")) %>%
    filter(!is.na(DeID_AdmitDate)) %>%
    mutate(PE = 1) 

diag_date$DeID_AdmitDate = as.Date(diag_date$DeID_AdmitDate, format = "%m/%d/%Y %H:%M")

# Como at PE's diagnosis
como_earliest = diag_date %>%
    right_join(como, by = c("DeID_PatientID", "DeID_EncounterID")) %>%
    filter(grepl("O14", TermCodeMapped)) %>%
    group_by(DeID_PatientID) %>%
    arrange(DeID_AdmitDate) %>%
    slice(1) %>%
    select(-DeID_EncounterID, -TermCodeMapped, -DeID_AdmitDate, -AgeInYears, -PE)

diag_date = diag_date %>%
    group_by(DeID_PatientID) %>%
    mutate(GestationalDiabetes = if_else(any(grepl("O24", TermCodeMapped)) & DeID_AdmitDate > min(DeID_AdmitDate[grepl("O24", TermCodeMapped)]), 1, 0))

# Filter diagnoses occurring 1~10 years after first PE
diag_date_list = list()
pop_case_list = list()

for (i in 1:10) {
    diag_date_i = diag_date %>%
        group_by(DeID_PatientID) %>%
        filter(any(grepl("O14", TermCodeMapped))) %>%
        filter(max(DeID_AdmitDate) >= (min(DeID_AdmitDate[grepl("O14", TermCodeMapped)]) + years(i))) %>%
        filter(DeID_AdmitDate < (min(DeID_AdmitDate[grepl("O14", TermCodeMapped)]) + years(i)))
    pop_case_i = length(unique(diag_date_i$DeID_PatientID))
    diag_date_list[[i]] = diag_date_i
    pop_case_list[[i]] = pop_case_i
}

# Convert ICD-10 to ElixComo
diag_elix_list = list()
for (i in 1:10) {
    diag_elix_i = as.data.frame(icd10_comorbid_elix(diag_date_list[[i]]))
    diag_elix_i = diag_elix_i %>%
        mutate(DeID_EncounterID = rownames(diag_elix_i)) %>%
        gather(key = "key", value = "value", -DeID_EncounterID) %>%
        filter(value) %>%
        group_by(DeID_EncounterID) %>%
        summarise(Comorbidity = list(key)) %>%
        ungroup()
    new_data = diag_date_list[[i]] %>% 
        left_join(diag_elix_i, by = "DeID_EncounterID") %>%
        unnest(Comorbidity) %>%
        select(-TermCodeMapped) %>%
        group_by(DeID_PatientID) %>%
        arrange(DeID_AdmitDate) %>%
        distinct(Comorbidity, .keep_all = TRUE) %>%
        relocate(Comorbidity, .before = "AgeInYears")
    diag_elix_list[[i]] = new_data
}

# Social demographic history at PE's diagnosis
social_hist = social %>%
    select(-TobaccoPacksPerDay, -TobaccoUsedYears) %>%
    mutate(SmokingStatusMapped = ifelse(SmokingStatusMapped %in% c("Current", "Former"), 1, 0)) %>%
    mutate(AlcoholUseStatusMapped = ifelse(AlcoholUseStatusMapped == "Yes", 1, 0)) %>%
    mutate(IllegalDrugUserStatusMapped = ifelse(IllegalDrugUserStatusMapped == "Yes", 1, 0)) %>%
    mutate(SexuallyActiveStatusMapped = ifelse(SexuallyActiveStatusMapped == "Yes", 1, 0)) %>%
    mutate(CigarettesYN = ifelse(CigarettesYN == "Y", 1, 0)) %>%
    mutate(PipesYN = ifelse(PipesYN == "Y", 1, 0)) %>%
    mutate(CigarsYN = ifelse(CigarsYN == "Y", 1, 0)) %>%
    mutate(SnuffYN = ifelse(SnuffYN == "Y", 1, 0)) %>%
    mutate(ChewYN = ifelse(ChewYN == "Y", 1, 0)) %>%
    mutate(IVDrugUserYN = ifelse(IVDrugUserYN == "Y", 1, 0)) %>%
    group_by(DeID_PatientID) %>%
    mutate(across(SmokingStatusMapped:IVDrugUserYN, ~ifelse(any(.x == 1), 1, 0 ))) %>%
    distinct(DeID_PatientID, .keep_all = TRUE) %>%
    select(-DeID_EncounterID)

# Race
race_hist = race %>%
    mutate(RaceName = if_else(!(RaceName %in% c("Caucasian", "African American")), "Other", RaceName)) %>%
    mutate(value = 1) %>%
    spread(key = RaceName, value = value, fill = 0) %>%
    select(-Line, -RaceCode) %>%
    distinct(DeID_PatientID, .keep_all = TRUE)

# Combine comorbidities, social hist, race
case_diag_all_elix_list = list()
for (i in 1:10) {
    diag_all_elix_i = diag_elix_list[[i]] %>% 
        left_join(como_earliest, by = "DeID_PatientID") %>%
        left_join(social_hist, by = "DeID_PatientID") %>%
        left_join(race_hist, by = "DeID_PatientID")
    
    case_diag_all_elix_list[[i]] = diag_all_elix_i
}

save(case_diag_all_elix_list, file = "case_diag_all_elix_list.RData")


# Basic statistics
case_stats_elix_list = list()
for (i in 1:10) {
    case_stats_elix_i = diag_elix_list[[i]] %>%
        group_by(Comorbidity) %>%
        summarise(count = n()) %>%
        mutate(Rate = count / pop_case_list[[i]])
    case_stats_elix_list[[i]] = case_stats_elix_i
}

save(case_stats_elix_list, file = "case_stats_elix_list.RData")

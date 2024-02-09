library(tidyverse)
library(touch)
library(icd)

diag = read.csv("DiagnosesComprehensiveAll.csv")
encounter = read.csv("EncounterAll.csv")
como = read.csv("ComorbiditiesElixhauserComprehensive.csv")
race = read.csv("PatientRace.csv")
social = read.csv("ClaritySocialHistory.csv")
death = read.csv("MichiganDeathIndex.csv")

diag_convert = diag

# Convert to ICD-10
diag_convert[diag_convert$Lexicon == "ICD9","TermCodeMapped"] = icd_map(diag_convert[diag_convert$Lexicon == "ICD9","TermCodeMapped"])
diag_convert$TermCodeMapped = gsub(",.*", "", diag_convert$TermCodeMapped)
diag_convert$TermCodeMapped[diag_convert$Lexicon == "ICD9"] = short_to_decimal(diag_convert$TermCodeMapped[diag_convert$Lexicon == "ICD9"])

# Match the encounter date with diagnoses
# Add a new column "PE" and set it to 0
diag_date = diag_convert %>%
    select(-TermNameMapped, -Lexicon) %>%
    left_join(encounter, by = c("DeID_PatientID", "DeID_EncounterID")) %>%
    filter(!is.na(DeID_AdmitDate)) %>%
    mutate(PE = 0) 

diag_date$DeID_AdmitDate = as.Date(diag_date$DeID_AdmitDate, format = "%m/%d/%Y %H:%M")

diag_date = diag_date %>%
    group_by(DeID_PatientID) %>%
    filter(any(grepl("Z34", TermCodeMapped))) %>%
    mutate(start_time = min(DeID_AdmitDate[grepl("Z34", TermCodeMapped)])) %>%
    mutate(GestationalDiabetes = if_else(any(grepl("O24", TermCodeMapped)) & DeID_AdmitDate > min(DeID_AdmitDate[grepl("O24", TermCodeMapped)]), 1, 0))

# Como at pregnancy diagnosis
como_earliest = diag_date %>%
    right_join(como, by = c("DeID_PatientID", "DeID_EncounterID")) %>%
    group_by(DeID_PatientID) %>%
    filter(DeID_AdmitDate <= start_time) %>%
    mutate(across(AidsHIV:WeightLoss, ~max(.))) %>%
    slice(1) %>%
    select(-DeID_EncounterID, -TermCodeMapped, -DeID_AdmitDate, -AgeInYears, -PE, -start_time, -GestationalDiabetes)

# Convert ICD-10 to ElixComo
diag_elix = as.data.frame(icd10_comorbid_elix(diag_date))
diag_elix = diag_elix %>%
    mutate(DeID_EncounterID = rownames(diag_elix)) %>%
    gather(key = "key", value = "value", -DeID_EncounterID) %>%
    filter(value) %>%
    group_by(DeID_EncounterID) %>%
    summarise(Comorbidity = list(key)) %>%
    ungroup()
new_data = diag_date %>% 
    left_join(diag_elix, by = "DeID_EncounterID") %>%
    unnest(Comorbidity) %>%
    select(-TermCodeMapped) %>%
    group_by(DeID_PatientID) %>%
    arrange(DeID_AdmitDate) %>%
    distinct(Comorbidity, .keep_all = TRUE) %>%
    relocate(start_time, .before = "DeID_AdmitDate")
diag_elix = new_data

death$DeID_MDIDeceasedDate = as.Date(death$DeID_MDIDeceasedDate, format = "%m/%d/%Y %H:%M")
death$DeID_RDWDeceasedDate = as.Date(death$DeID_RDWDeceasedDate, format = "%m/%d/%Y %H:%M")
death = death %>%
    group_by(DeID_PatientID) %>%
    mutate(DeceaseDate = min(DeID_RDWDeceasedDate ,DeID_MDIDeceasedDate, na.rm = TRUE)) %>%
    select(-DeID_MDIDeceasedDate, -DeID_RDWDeceasedDate)
death_como = as.data.frame(icd10_comorbid_elix(death))
death_como = death_como %>%
    mutate(DeID_PatientID = rownames(death_como)) %>%
    gather(key = "key", value = "value", -DeID_PatientID) %>%
    filter(value) %>%
    group_by(DeID_PatientID) %>%
    summarise(Comorbidity = list(key))
death = death %>%
    left_join(death_como, by = "DeID_PatientID") %>%
    select(-UnderlyingCOD_ICD10)

# Target diseases
target = c("HTN", "DMcx", "CHF", "Renal", "Obesity", "Hypothyroid")
diag_date_list = list()

for (disease in target) {
    diag_date_list[[disease]] = diag_elix %>%
        group_by(DeID_PatientID) %>%
        filter(DeID_AdmitDate > start_time) %>%
        mutate(status = ifelse(any(Comorbidity == disease), 1, 0)) %>%
        mutate(time = ifelse(any(Comorbidity == disease), min(DeID_AdmitDate[Comorbidity == disease]) - start_time, max(DeID_AdmitDate - start_time))) %>%
        mutate(AgeInYears = ifelse(any(Comorbidity == disease), min(AgeInYears[Comorbidity == disease]), max(AgeInYears))) %>%
        distinct(DeID_PatientID, .keep_all = TRUE) %>%
        select(-Comorbidity)
    
    temp = death %>%
        inner_join(diag_date_list[[disease]], by = "DeID_PatientID")%>%
        select(start_time, DeID_PatientID, Comorbidity, DeceaseDate) %>%
        mutate(time = as.double(DeceaseDate - start_time)) %>%
        mutate(status = ifelse(Comorbidity == disease, 1, 2)) %>%
        select(-Comorbidity, -DeceaseDate, -start_time)
    
    diag_date_list[[disease]] = diag_date_list[[disease]] %>%
        left_join(temp, by = "DeID_PatientID") 
    diag_date_list[[disease]] = diag_date_list[[disease]] %>%
        mutate(status = coalesce(status.y, status.x)) %>%
        mutate(time = coalesce(time.y, time.x)) %>%
        select(-status.x, -status.y, -time.x, -time.y) %>%
        relocate(status, .before = "AgeInYears") %>%
        relocate(time, .before = "AgeInYears")
}

# Social demographic history at first pregnancy's diagnosis
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


# Add comorbidities, social hist, race

control_survival_list = list()

# Combine data
for (disease in target) {
    control_survival_list[[disease]] = diag_date_list[[disease]] %>% 
        left_join(como_earliest, by = "DeID_PatientID") %>%
        left_join(social_hist, by = "DeID_PatientID") %>%
        left_join(race_hist, by = "DeID_PatientID")
}

save(control_survival_list, file = "control_survival_list.RData")


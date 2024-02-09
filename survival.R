# Import libraries
library(survival)
library(tidyverse)
library(survminer)
library(ggsurvfit)
library(tidycmprsk)

### Regression analysis

load("case_survival_list.RData")
load("control_survival_list.RData")

target = c("HTN", "DMcx", "CHF", "Renal", "Obesity", "Hypothyroid")

survival_list = list()
for (disease in target) {
    # Merge case & control
    survival_list[[disease]] = rbind(case_survival_list[[disease]], control_survival_list[[disease]])
    survival_list[[disease]] = survival_list[[disease]] %>%
        ungroup() %>%
        mutate(across(6:ncol(survival_list[[disease]]), as.double))
    
    # Find the effective variables
    sums = case_survival_list[[disease]] %>%
        ungroup() %>%
        distinct(DeID_PatientID, .keep_all = TRUE) %>%
        select(DeID_PatientID, 9:ncol(case_survival_list[[disease]])) %>%
        mutate(across(everything(), ~ifelse(is.na(.)|.x == 0, 0, 1))) %>%
        select(-DeID_PatientID) %>%
        colSums()
    na_prop = colSums(is.na(survival_list[[disease]][,9:ncol(survival_list[[disease]])])) / nrow(survival_list[[disease]])
    
    # Remove variables that have more than 20% NA
    ineffective_var = names(na_prop[na_prop >= 0.2])
    
    # Final variables
    survival_list[[disease]] = survival_list[[disease]] %>%
        select(-all_of(ineffective_var))

}



# Basic statistics
pop_case = length(unique(case_survival_list[[1]]$DeID_PatientID))
pop_control = length(unique(control_survival_list[[1]]$DeID_PatientID))

case_stats = case_survival_list[[1]] %>%
    distinct(DeID_PatientID, .keep_all = TRUE) %>%
    ungroup() %>%
    select(8:last_col()) %>%
    replace(is.na(.), 0) %>%
    colSums()
case_stats = data.frame(Comorbidity = names(case_stats), case = case_stats/pop_case)

control_stats = control_survival_list[[1]] %>%
    distinct(DeID_PatientID, .keep_all = TRUE) %>%
    ungroup() %>%
    select(8:last_col()) %>%
    replace(is.na(.), 0) %>%
    colSums()
control_stats = data.frame(Comorbidity = names(control_stats), control = control_stats/pop_control)

sur_stats = case_stats %>%
    left_join(control_stats, by = "Comorbidity") %>%
    group_by(Comorbidity) %>%
    mutate(pvalue = prop.test(c(case*pop_case, control*pop_control), c(pop_case, pop_control))$p.value)



# Survival analysis
survival_model_list = list()
result_survival = data.frame()
for (disease in target) {
    survival_list[[disease]] = survival_list[[disease]] %>%
        mutate(status = as.factor(status))
    formula = as.formula("Surv(time, ifelse(status == 1, 1, 0)) ~  AgeInYears + PE + GestationalDiabetes + DiabetesComplicated + HypertensionUncomplicated + Hypothyroidism + Obesity + RenalFailure + CongestiveHeartFailure + SmokingStatusMapped + AlcoholUseStatusMapped + `African American` + Caucasian")  
    model = coxph(formula, data = survival_list[[disease]])
    survival_model_list[[disease]] = model
    result_survival = result_survival %>% 
        bind_rows(data.frame(Comorbidity = disease,
                             coef = summary(model)$coefficients["PE", "coef"],
                             stderr = summary(model)$coefficients["PE", "se(coef)"],
                             pvalue = summary(model)$coefficients["PE", "Pr(>|z|)"]
        )
        )
}

save(survival_model_list, file = "survival_model_list.RData")
save(result_survival, file = "result_survival.RData")


# Competing risk +
survival_model_cmprsk_list = list()
result_survival_cmprsk = data.frame()
for (disease in target) {
    survival_list[[disease]] = survival_list[[disease]] %>%
        mutate(status = as.factor(status))
    formula = as.formula("Surv(time, status) ~  AgeInYears + PE + GestationalDiabetes + DiabetesComplicated + HypertensionUncomplicated + Hypothyroidism + Obesity + RenalFailure + CongestiveHeartFailure + SmokingStatusMapped + AlcoholUseStatusMapped + `African American` + Caucasian")  
    model = crr(formula, data = survival_list[[disease]])
    survival_model_cmprsk_list[[disease]] = model
    temp = model$tidy
    result_survival_cmprsk = result_survival_cmprsk %>%
        bind_rows(data.frame(Comorbidity = disease,
                             coef = temp$estimate[temp$term == "PE"],
                             stderr = temp$std.error[temp$term == "PE"],
                             pvalue = temp$p.value[temp$term == "PE"]
        )
        )
}

save(result_survival_cmprsk, file = "result_survival_cmprsk.RData")
save(survival_model_cmprsk_list, file = "survival_model_cmprsk_list.RData")

# Stratify races
result_survival_race = data.frame()
for (disease in target) {
    temp = survival_list[[disease]] %>% 
        filter(Caucasian == 1) %>%
        mutate(status = as.factor(status))
    formula = as.formula("Surv(time, ifelse(status == 1, 1, 0)) ~  AgeInYears + PE + GestationalDiabetes + DiabetesComplicated + HypertensionUncomplicated + Hypothyroidism + Obesity + RenalFailure + CongestiveHeartFailure + SmokingStatusMapped + AlcoholUseStatusMapped + `African American` + Caucasian")  
    model = coxph(formula, data = temp)
    result_survival_race = result_survival_race %>% 
        bind_rows(data.frame(Comorbidity = disease,
                             coef = summary(model)$coefficients["PE", "coef"],
                             stderr = summary(model)$coefficients["PE", "se(coef)"],
                             pvalue = summary(model)$coefficients["PE", "Pr(>|z|)"],
                             race = "Caucasian"
        )
        )
    
    temp = survival_list[[disease]] %>% 
        filter(`African American` == 1) %>%
        mutate(status = as.factor(status))
    formula = as.formula("Surv(time, ifelse(status == 1, 1, 0)) ~  AgeInYears + PE + GestationalDiabetes + DiabetesComplicated + HypertensionUncomplicated + Hypothyroidism + Obesity + RenalFailure + CongestiveHeartFailure + SmokingStatusMapped + AlcoholUseStatusMapped + `African American` + Caucasian")  
    model = coxph(formula, data = temp)
    result_survival_race = result_survival_race %>% 
        bind_rows(data.frame(Comorbidity = disease,
                             coef = summary(model)$coefficients["PE", "coef"],
                             stderr = summary(model)$coefficients["PE", "se(coef)"],
                             pvalue = summary(model)$coefficients["PE", "Pr(>|z|)"],
                             race = "African American"
        )
        )
}


result_survival = result_survival %>%
    mutate(xmin=exp(coef-1.96*stderr), xmax=exp(coef+1.96*stderr), x = exp(coef))

result_survival_cmprsk = result_survival_cmprsk %>%
    mutate(xmin=exp(coef-1.96*stderr), xmax=exp(coef+1.96*stderr), x = exp(coef))

result_survival_race = result_survival_race %>%
    mutate(xmin=exp(coef-1.96*stderr), xmax=exp(coef+1.96*stderr), x = exp(coef))

# Full names
replacement_rules <- c(
    "CHF" = "Congestive Heart Failure",
    "DMcx" = "Diabetes Complicated",
    "HTN" = "Hypertension Uncomplicated",
    "Hypothyroid" = "Hypothyroidism",
    "Renal" = "Renal Disease",
    "Obesity" = "Obesity"
)
result_survival$Comorbidity = replacement_rules[result_survival$Comorbidity]
result_survival_race$Comorbidity = replacement_rules[result_survival_race$Comorbidity]
result_survival_cmprsk$Comorbidity = replacement_rules[result_survival_cmprsk$Comorbidity]
fullname =c("Hypertension Uncomplicated", "Diabetes Complicated", "Congestive Heart Failure", "Renal Disease", "Obesity", "Hypothyroidism")



### Plotting

# Survival curves
i = 1
for (disease in target) {
    sdata = survival_list[[disease]] %>%
        mutate(time = time %/% 7)
    sfit = survfit(Surv(time, ifelse(status == 1, 1, 0), type = "right")~PE, data = sdata)
    temp = ggsurvplot(sfit, pval = TRUE, pval.coord = c(0, max(min(sfit$lower), 0)), 
                      conf.int=TRUE, ylim = c(max(min(sfit$lower)-0.05, 0), 1),
                      legend.title = "Group", legend.labs = c("Control", "Case"), palette=c("#2E9FDF", "#E7B800"),
                      title = paste0("Survival Curve of ", fullname[i]), xlab = "Time (weeks)", ylab = "Survival Probability")
    png(height = 4000, width = 4000, file = paste0("./figures/Survival_", disease, ".png"), res = 600)
    print(temp)
    dev.off()
    i = i + 1
}


# Survival curves by race
i = 1
for (disease in target) {
    sdata = survival_list[[disease]] %>%
        filter(Caucasian == 1 | `African American` == 1) %>%
        mutate(time = time %/% 7) %>%
        mutate(Race = ifelse(Caucasian == 1, "Caucasian", "African American")) %>%
        mutate(Group = ifelse(PE == 1, "Case", "Control"))
    sdata = as.data.frame(sdata)
    sfit = survfit(Surv(time, ifelse(status == 1, 1, 0), type = "right")~Group+Race, data = sdata)
 
    temp = ggsurvplot(sfit, sdata, pval = TRUE,
                      conf.int=TRUE,
                      legend.title = "Group", 
                      title = paste0("Survival Curve of ", fullname[i]), xlab = "Time (weeks)", ylab = "Survival Probability",
                      # palette = c("#D0717C","#6DA6CF"), 
                      font.tickslab = c(12),
                      facet.by = "Group",
    )
    
    png(height = 3600, width = 5000, file = paste0("./figures/Survival_race_", disease, "_2.png"), res = 600)
    print(temp)
    dev.off()
    i = i + 1
}

# Harzard Ratios
temp = ggplot(data=result_survival, aes(y= Comorbidity, x=exp(coef), xmin=exp(coef-1.96*stderr), xmax=exp(coef+1.96*stderr), color = Comorbidity)) +
    theme_minimal()+
    geom_errorbarh(height=.1)+
    geom_point() + 
    geom_vline(xintercept = 1, linetype = "dashed") +
    labs(title = "Hazard Ratios of PE", x = "Hazard Ratio", y = "Comorbidity") +
    theme(legend.position = "none")+
    scale_x_continuous(breaks = c(1,2,4,6))
ggsave("./figures/HazardRatio.png", height = 4.8, width = 7.2, dpi = 600)

# Harzard Ratios + competing risks
temp = ggplot(data=result_survival_cmprsk, aes(y= Comorbidity, x=exp(coef), xmin=exp(coef-1.96*stderr), xmax=exp(coef+1.96*stderr), color = Comorbidity)) +
    theme_minimal()+
    geom_errorbarh(height=.1)+
    geom_point() + 
    geom_vline(xintercept = 1, linetype = "dashed") +
    labs(title = "Hazard Ratios, UM Medicine", x = "Hazard Ratio", y = "Comorbidity") +
    theme(legend.position = "none")+
    scale_x_continuous(breaks = c(1,2,4,6))
ggsave("./figures/HazardRatio_CMPRSK.png", height = 4, width = 6)

# Harzard Ratios by race
temp = ggplot(data=result_survival_race, aes(y= Comorbidity, x=exp(coef), xmin=exp(coef-1.96*stderr), xmax=exp(coef+1.96*stderr), color = race, group = race)) +
    theme_minimal()+
    geom_errorbarh(height=0.2, position = position_dodge(width = 0.5))+
    geom_pointrange(position = position_dodge(width = 0.5), size = 0.15) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    labs(title = "Hazard Ratios of PE", x = "Hazard Ratio", y = "Comorbidity") +
    scale_x_continuous(breaks = c(1,2,4,6)) 
ggsave("./figures/HazardRatio_Race.png", height = 4.8, width = 7.2, dpi = 600)



# Import libraries
library(tidyverse)
library(reshape2)
library(gridExtra)
library(corrplot)
library(qacReg)

### Odds Ratios

trend_odds = data.frame(Comorbidity = result_elix_list[[1]]$Comorbidity)
trend_omin = data.frame(Comorbidity = result_elix_list[[1]]$Comorbidity)
trend_omax = data.frame(Comorbidity = result_elix_list[[1]]$Comorbidity)

for (i in 1:10) {
    temp = result_elix_list[[i]] %>% select(Comorbidity, odds)
    trend_odds = trend_odds %>%
        full_join(temp, by = "Comorbidity")
    colnames(trend_odds)[i+1] = paste0(i)
    
    temp = result_elix_list[[i]] %>% select(Comorbidity, omin)
    trend_omin = trend_omin %>%
        full_join(temp, by = "Comorbidity")
    colnames(trend_omin)[i+1] = paste0(i)
    
    temp = result_elix_list[[i]] %>% select(Comorbidity, omax)
    trend_omax = trend_omax %>%
        full_join(temp, by = "Comorbidity")
    colnames(trend_omax)[i+1] = paste0(i)
}

trend_odds_long <- melt(trend_odds, id.vars = "Comorbidity", variable.name = "Year", value.name = "Odds")
trend_odds_long$Year = as.integer(trend_odds_long$Year)

trend_omin_long <- melt(trend_omin, id.vars = "Comorbidity", variable.name = "Year", value.name = "omin")
trend_omin_long$Year = as.integer(trend_omin_long$Year)

trend_omax_long <- melt(trend_omax, id.vars = "Comorbidity", variable.name = "Year", value.name = "omax")
trend_omax_long$Year = as.integer(trend_omax_long$Year)

trend_odds = trend_odds_long %>% 
    left_join(trend_omin_long, by = c("Comorbidity", "Year")) %>%
    left_join(trend_omax_long, by = c("Comorbidity", "Year"))

temp = ggplot(data = trend_odds, aes(x = Year, y = Odds, color = Comorbidity)) +
    ylim(0, max(trend_odds$omax)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "black") +
    geom_point() +
    geom_line(linewidth = 0.5) +
    geom_errorbar(aes(ymin = omin, ymax = omax), width = 0.1) +
    geom_ribbon(aes(ymin = omin, ymax = omax, fill = Comorbidity), alpha = 0.2, color = NA) +
    facet_wrap(~ Comorbidity, scales = "free_y", nrow = 5) +
    xlab("Year") +
    ylab("Odds Ratio") +
    ggtitle("Trends of Elixhauser Comorbidities Odds Ratio During 1~10 Years") +
    theme(legend.position = "none") +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10))
ggsave("./figures/TrendOdds.png", temp, height = 12, width = 18)


### Ages

ages = all_final_elix_list[[1]] %>%
    mutate(Group = if_else(PE == 1, "Case", "Control")) %>%
    group_by(DeID_PatientID) %>%
    reframe(AgeInYears = min(AgeInYears), Group = Group[1])
t.test(ages$AgeInYears[ages$Group == "Case"], ages$AgeInYears[ages$Group == "Control"])

plot_ages = ggplot(ages, aes(x = Group, y = AgeInYears, fill = Group)) +
    geom_boxplot() +
    labs(title = "Age Distribution, 5 Years", x = "Group", y = "AgeInYears") +
    theme_minimal()
ggsave("./figures/Ages.png", plot_ages)


### VIF

for (i in 1:10) {
    plot_i = vif_plot(model_list[[i]]) +
        ggtitle(paste0("Variance Inflation Plot, ", i, " Years"))
    ggsave(paste0("./figures/VIFPlot_", i, ".png"), plot_i, height = 6, width = 6)
}


### Populations

pop_case = unlist(pop_case_list)
pop_control = unlist(pop_control_list)
pop_all = data.frame(Year = 1:10, case = pop_case, control = pop_control)

temp = ggplot(data = pop_all, aes(x = Year)) +
    geom_point(aes(y = case, color = "case")) +
    geom_line(aes(y = case, color = "case")) +
    geom_point(aes(y = control, color = "control")) +
    geom_line(aes(y = control, color = "control")) +
    labs(title = "Populations in Case and Control Group after 1~10 Years", x = "Year", y = "Population", color = "Group") +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10))

ggsave("./figures/Population.png", temp, height = 6, width = 8)

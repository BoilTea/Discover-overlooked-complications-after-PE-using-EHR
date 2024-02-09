# Discover-overlooked-complications-after-PE-using-EHR
R codes for the manuscript "Discover overlooked complications after preeclampsia using electronic health records"

## Author
Haoming Zhu (University of Michigan)

## Description for each R script
case/control_elix.R: preprocess cases and controls  data for Logistic Regression

logit_regression_elix.R: perform Logistic Regression on each Elixhauser Comorbidity (complications)

logit_plot.R: plot the Logistic Regression results

case/control_survival.R: preprocess cases and controls  data for survival analysis

survival.R: perform survival analysis on significant Elixhauser Comorbidities (complications)

### Install necessary packages (ONLY if not previously installed)

#install.packages("dplyr")
#install.packages("tidyr")
#install.packages("ggplot2")

### Access necessary packages

library(tidyr)
library(dplyr)
library(ggplot2)

### Reference the plaquant.R script as the source for functions

source("plaquant.R")


### Establish names of columns containing count data and all possible names of control and treated samples.

vis_data <- read.csv("EC virus summary.csv")

for(i in 1:dim(vis_data)[1]){
  if(vis_data$logR[i] < 0){
    vis_data$logRstack[i] <- 0
  }
  else(vis_data$logRstack[i] <- vis_data$logR[i])
}

#
##
### Statistical analyses
##
#


# Find differences in mean and cumulatve standard error for each test compared to baseline

for(i in 1: dim(vis_data)[1]){
  corresponding_baseline <- filter(vis_data, Pathogen == vis_data$Pathogen[i], Reduction == vis_data$Reduction[i], Test == "Baseline")
  vis_data$Meandiff[i] <- vis_data[i, "logR"] - corresponding_baseline$logR
  vis_data$CumulativeSE[i] <- sqrt(vis_data[i, "logR_SE"]^2 + corresponding_baseline$logR_SE^2)
}

# Calculate t-values for independent t-test (Ho: mu1 - mu2 = 0)

Ho <- 0
vis_data$t_value <- abs((vis_data$Meandiff - Ho) / vis_data$CumulativeSE)

#all tests in triplicate, assume df = 3 + 3 - 2

df = 4

#calculate p-values

vis_data$p_value <- dt(vis_data$t_value, df)

# compare to required alpha value for a 2-tailed test with Bonferroni adjustment

confidence <- 0.95

req_alpha <- (1-confidence)/2

# Bonferroni correction -- adjust as appropriate for number of comparisons

bonferroni_n <- 3

corr_alpha <- req_alpha / bonferroni_n


# if p-value based on calculated t-value is less than the corrected alpha, 
#       signify significant difference with an asterisk (*) 

vis_data$t_decision <- ""

for(i in 1:dim(vis_data)[1]){
  if(is.na(vis_data$p_value[i])){vis_data$t_decision[i] <- ""}
  else{
    if(vis_data$p_value[i] <= corr_alpha){
    vis_data$t_decision[i] <- "*"
      }
    }
}

write.csv(vis_data, "virus output.csv")

#
##
### Data visualization
##
#


# Select the tests, pathogens and treatments of interest

tests <- c("Baseline", "NOM", "Turbidity", "Chloride")
pathogens <- c("ADV", "ECV", "FCV")
treatments <- c("beef broth", "filtered")

# Create a stacked bar chart

stacked <- logR_bar(vis_data, tests, pathogens, treatments, layout = "stacked", significance = TRUE, y_offset = 0.6)


# Change to Marquette colors and change the axis labels and tick marks

stacked +
  scale_fill_manual(labels = c("Coagulation", "Inactivation"), values=c("#003366", "#FFCC00")) + 
  scale_color_manual(values=c("#003366", "#FFCC00"), guide = "none") +
  labs(x = "", y = "Log Reduction") + 
  ylim(0,6) +
  scale_x_discrete(limits = tests, labels=c("Baseline", "33 mg/L NOM", "50 NTU", "115 mg/L Cl")) +
  theme(axis.title = element_text(face="bold", size=12), axis.text = element_text(face="bold", size=11), axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.5))


### pH

# Select the tests, pathogens and treatments of interest

tests <- c("pH 6", "Baseline", "pH 8")
pathogens <- c("ADV", "ECV", "FCV")
treatments <- c("beef broth", "filtered")

# Create a staggered bar chart using logR_bar()

stacked <- logR_bar(vis_data, tests, pathogens, treatments, layout = "stacked", significance = FALSE)


# Change to Marquette colors and change the axis labels and tick marks

stacked +
  scale_fill_manual(labels = c("Coagulation", "Inactivation"), values=c("#003366", "#FFCC00")) + 
  scale_color_manual(values=c("#003366", "#FFCC00")) +
  labs(x = "pH", y = 'Log Reduction') + 
  ylim(0,6) +
  scale_x_discrete(limits = tests, labels=c("6", "7", "8")) +
  theme(axis.title = element_text(face="bold", size=12), axis.text = element_text(face="bold", size=11))
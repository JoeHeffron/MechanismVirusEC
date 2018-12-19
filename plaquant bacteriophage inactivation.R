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

count_cols <- c("Count1", "Count2", "Count3", "Count4", "Count5", "Count6", "Count7", "Count8", "Count9", "Count10")
control_names <- c("control")
treated_names <- c("1", "2", "3")


#
##
### Data organization
##
#

tidy <- data_tidy("Inactivation Master Treat 20180722.csv", count_cols, control_names, treated_names)

stat_data <- stat_transform(tidy)

total_recovery <- c("beef broth")

vis_data <- vis_transform(tidy, total_recovery)

#
##
### Statistical analyses
##
#


## Effect of pH on inactivation

pH_tests <- c("pH 6", "Baseline", "pH 8")

pH_data <- filter(stat_data, Sample == "treated")

pH_data <- filter(pH_data, Test %in% pH_tests)

pH_values <- c(6, 7, 8)

pH_data$pH <- pH_values[match(pH_data$Test, pH_tests)]

pH_data <- mutate(pH_data, pH2 = pH^2)

pH_data <- mutate(pH_data, OH = 10^-(14-pH))

pH_data <- mutate(pH_data, invOH = OH^-1)

pH_data <- mutate(pH_data, invOH2 = OH^-2)

pH_data <- filter(pH_data, Treatment == "beef broth")

pHMS2 <- filter(pH_data, Pathogen == "MS2")

pHfr <- filter(pH_data, Pathogen == "fr")

pHphi <- filter(pH_data, Pathogen == "phiX174")

pHP22 <- filter(pH_data, Pathogen == "P22")


Mlm <- lm(logR ~ pH + pH2, data = pHMS2)

flm <- lm(logR ~ pH + pH2, data = pHfr)

philm <- lm(logR ~ pH + pH2, data = pHphi)

Plm <- lm(logR ~ pH + pH2, data = pHP22)


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

# Bonferroni correction, assuming all Pathogens have the same number and name of tests
bonferroni_n <- length(unique(vis_data$Test)) - 1


corr_alpha <- req_alpha / bonferroni_n


# if p-value based on calculated t-value is less than the corrected alpha, 
#       signify significant difference with an asterisk (*) 

vis_data$t_decision <- ""

for(i in 1:dim(vis_data)[1]){
  if(vis_data$p_value[i] <= corr_alpha){
    vis_data$t_decision[i] <- "*"
  }
}

write.csv(vis_data, "phage output.csv")


data_tidy()

#
##
### Data visualization
##
#

# Select the tests, pathogens and treatments of interest

tests <- c("Baseline", "NOM", "Turbidity", "Chloride")
pathogens <- c("MS2", "phiX174", "fr", "P22")
treatments <- c("beef broth", "filtered")

# Create a stacked bar chart using logR_bar()


stacked <- logR_bar(vis_data, tests, pathogens, treatments, layout = "stacked", significance = TRUE )


# Change to Marquette colors and change the axis labels and tick marks

stacked +
  scale_fill_manual(labels = c("Coagulation", "Inactivation"), values=c("#003366", "#FFCC00"), guide = "none") +
  scale_color_manual(values=c("#003366", "#FFCC00"), guide = "none") +
  labs(x = "", y = bquote("Log Reduction")) + 
  scale_x_discrete(limits =  tests, labels=c("Baseline", "33 mg/L NOM", "50 NTU", "115 mg/L Cl")) +
  ylim(0,6) +
  theme(axis.title = element_text(face="bold", size=12), axis.text = element_text(face="bold", size=11), axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.5))



### pH

# Select the tests, pathogens and treatments of interest

# List tests in order of appearance on graph
tests <- c("pH 6", "Baseline", "pH 8")
pathogens <- c("MS2", "P22", "phiX174", "fr")
treatments <- c("beef broth", "filtered")

# Create a stacked bar chart using logR_bar()

stacked <- logR_bar(vis_data, tests, pathogens, treatments, layout = "stacked", significance = FALSE)

# Change to Marquette colors and change the axis labels and tick marks

stacked +
  scale_fill_manual(labels = c("Coagulation", "Inactivation"), values=c("#003366", "#FFCC00"), guide = "none") + 
  scale_color_manual(values=c("#003366", "#FFCC00"), guide = "none") +
  labs(x = "pH", y = "Log Reduction") + 
  scale_x_discrete(limits = c("pH 6", "Baseline", "pH 8"), labels=c("6", "7", "8")) +
  theme(axis.title = element_text(face="bold", size=12), axis.text = element_text(face="bold", size=11))

### chem coag
tests <- c("Baseline", "CCC-FeCl3-2mg", "CCC-FeCl2-2mg", "pre-formed floc", "Titanium")
pathogens <- c("MS2", "P22", "phiX174", "fr")
treatments <- c("beef broth", "filtered")

levels(vis_data$Test) <- tests

# Create a stacked bar chart using logR_bar()

stacked <- logR_bar(vis_data, tests, pathogens, treatments, layout = "stacked", significance = TRUE, y_offset = 1)


# Change to Marquette colors and change the axis labels and tick marks

stacked +
  scale_fill_manual(labels = c("Coagulation", "Inactivation"), values=c("#003366", "#FFCC00")) + 
  scale_color_manual(values=c("#003366", "#FFCC00"), guide = "none") +
  labs(x = "", y = "Log Reduction") + 
  ylim(0, 6)+
  scale_x_discrete(limits = tests, labels=c("EC", "FeCl3", "FeCl2", "Pre-formed floc", "Titanium")) +
  theme(axis.title = element_text(face="bold", size=12), axis.text = element_text(face="bold", size=11), axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.5))

### "plaquant" 
# Plaque/colony count data organization, analysis and visualization
# Author: Joe Heffron 
# Copyright 2017, Joe Heffron


#
##
### Data organization
##
#

data_tidy <- function(filename, count_cols, control_names = c("control"), treated_names = c("treated")){
  
  require(tidyr)
  require(dplyr)

  # open file
  data <- read.csv(filename, stringsAsFactors = F)
  
  # Remove rows with no count data ("NA" for Dilution)
  
  for(i in 1:dim(data)[1]){
    if(is.na(data[i, "Dilution"])){
    data <- data[-i, ]
    }
  }
  
  # remove all but lowest dilution for each sample
  
  by_Sample <- group_by(data, Test, Sample, Treatment, Pathogen)
  
  clean_data <- filter(by_Sample, min_rank(Dilution) <= 1)
  
  
  ### convert observations with zeros to ones
  
  clean_data[clean_data == 0] <- 1
  
  ### convert sample names to standard names ("control" and "treated")

  clean_data$Sample[which(clean_data$Sample %in% control_names)] <- "control"
  
  clean_data$Sample[which(clean_data$Sample %in% treated_names)] <- "treated"


  ### Melt data to convert counts from wide to long format 
  
  data.tidy <- gather(clean_data, key, Count, one_of(count_cols), na.rm = TRUE, convert = TRUE)
  data.tidy <- select(data.tidy, -key)
    
  return(data.tidy)
  }

#### Generic function for selecting tests based on Test, Pathogen and Treatment

data_selector <- function(data, tests, pathogens, treatments){
  data_temp <- filter(data, Test %in% tests & Pathogen %in% pathogens & Treatment %in% treatments)
  return(data_temp)
}

  
#
##
### Data transformation - statistics
##
#

stat_transform <- function(data.tidy){

  ### Derive a concentration in PFU/mL from the count and dilution
  
  data.trans <- mutate(data.tidy, Concentration  = Count * 10^Dilution)
  
  
  ### take log of all concentrations
  data.trans <- mutate(data.trans, logC = log10(Concentration))
  
  ### remove count, dilution and key columns
  data.trans <- select(data.trans, -(Dilution:Concentration))
  
  
  ### Summarize control data for each set of replicates
  
  controls <- filter(data.trans, Sample == "control")
  control.sum <- aggregate(logC ~ Treatment + Sample + Pathogen + Test, controls, function(x) c(M = mean(x), Weight = 1/var(x)))
  control.sum[,5:6] <- control.sum[,5]
  colnames(control.sum)[5:6] <- c("logC", "Weight")
  
  ### Find corresponding control mean log C for each treated logC
  treateds <- filter(data.trans, Sample == "treated")
  
  for(i in 1:dim(data.trans)[1]){
    if(data.trans[i, "Sample"] == "treated"){
    corresponding_control <- filter(control.sum, Pathogen == data.trans$Pathogen[i], Test == data.trans$Test[i], Treatment == data.trans$Treatment[i])
    data.trans[i,"logR"] <- corresponding_control$logC - data.trans[i, "logC"]
    }
  }
  
  return(data.trans)
}

#
##
### Data transformation - visualization
##
#


vis_transform <- function(data.tidy, recovery_treatment){
  
  ### Define basic functions
  
  std.err <- function(x) {
    sd(x)/sqrt(length(x))
  }
  
  ### Derive a concentration in PFU/mL from the count and dilution
  data.trans <- mutate(data.tidy, Concentration  = Count * 10^Dilution)
  
  ### take log of all concentrations
  data.trans <- mutate(data.trans, logC = log10(Concentration))
  
  
  ### Summarize data for each set of replicates
  data.sum <- aggregate(logC ~ Treatment + Sample + Pathogen + Test, data.trans, function(x) c(M = mean(x), SE = std.err(x), n = length(x)))
  data.sum[,5:7] <- data.sum[,5]
  agNames <- c("logSE", "n")
  colnames(data.sum)[6:7] <- agNames
  
  
  ### determine log Reductions and SE of log Reductions
  
  data.sum <- data.sum %>%
    gather(temp, stats, logC, logSE, n) %>%
    unite(temp2, Sample, temp) %>%
    spread(temp2, stats)
  
  # Compute log Reduction
  data.sum <- mutate(data.sum, logR = control_logC - treated_logC)
  
  # Compute SE of log Reduction
  data.sum <- mutate(data.sum, logR_SE = sqrt(control_logSE^2 + treated_logSE^2))
  
  ### determine reversible and irreversible log Reductions
  
  data.reduction <- select(data.sum, -contains("control"), -contains("treated"))
  
  
  for(i in 1:dim(data.reduction)[1]){
    if(data.reduction$Treatment[i] %in% recovery_treatment){
      data.reduction$Reduction[i] <- "Irreversible"
    }
    else{data.reduction$Reduction[i] <- "Total"}
  }
  
  
  data.reduction <- data.reduction %>%
    gather(temp, value, logR, logR_SE) %>%
    unite(temp2, Reduction, temp) %>%
    spread(temp2, value)
  
  # Separate Irreversible and Total reduction data
  
  Total <- filter(data.reduction, !is.na(Total_logR))
  
  data.reduction <- filter(data.reduction, !is.na(Irreversible_logR))
  
  # Fill total log reduction in to irreversible data (fewer total reduction stats due to different elution methods used for irreversible reduction)
  
  for(i in 1:dim(data.reduction)[1]){
    tot_rem <- filter(Total, Pathogen == data.reduction$Pathogen[i], Test == data.reduction$Test[i])
    data.reduction[i,"Total_logR"] <- tot_rem$Total_logR
    data.reduction[i, "Total_logR_SE"] <- tot_rem$Total_logR_SE
  }
  
  # Calculate Reversible reduction from Total and Irreversible
  
  data.reduction <- mutate(data.reduction, Reversible_logR = Total_logR - Irreversible_logR)
  
  ## Calculate Reversible reduction standard error based on additive standard error
  data.reduction <- mutate(data.reduction, Reversible_logR_SE = sqrt(Total_logR_SE^2 + Irreversible_logR_SE^2))
  

  # Convert log reduction and SE values back to long format, with Reversible, Irreversible or Total as a new "key" column
  
  data.reduction <- data.reduction %>%
    gather(key, value, -Treatment, -Pathogen, -Test) %>%
    separate(key, into = c("Reduction", "stat"), sep = "_", extra = "merge") %>%
    spread(stat, value)
  
  # For stacked bar charts, convert negative log removal to zero and adjust SE
  
  for(i in 1:dim(data.reduction)[1]){
    if(data.reduction$logR[i] < 0){
      data.reduction$logRstack[i] <- 0
    }
    else{data.reduction$logRstack[i] <- data.reduction$logR[i]}
    }
    
  
  return(data.reduction)
}

#
##
### Data Visualization
##
#


logR_bar<- function(data, tests, pathogens, treatments, layout = "stagger", significance = TRUE, y_offset = 0.5){
  require(dplyr)
  require(ggplot2)
  # limit results to those tests and treatments the user is interested in
  results_temp <- data_selector(data, tests, pathogens, treatments)
  

  if(layout == "stagger"){
    
    # set values for error bars
    
    results_temp$LRmax <- with(results_temp, logR + logR_SE)
    results_temp$LRmin <- with(results_temp, logR - logR_SE)
    
    # limit data to total Reduction and irreversible Reduction
    results_temp <- filter(results_temp, Reduction == "Total" | Reduction == "Irreversible")
    # reorder the factors so that Irreversible Reduction is printed after Total Reduction
    results_temp$Reduction <- ordered(results_temp$Reduction, c("Irreversible", "Total"))
    ## compile staggered bar chart
    staggered.bar <- ggplot(results_temp, aes(x = Test, y = logR, fill = Reduction)) +
      geom_bar(stat = "identity", position = "dodge") +
      facet_grid(.~Pathogen) +
      geom_errorbar(aes(ymax=(LRmax),
                        ymin=(LRmin)),
                    position=position_dodge(0.9),
                    data=results_temp)
    return(staggered.bar)
  }
  else{
    if(layout == "stacked"){
      # set values for error bars
      
      results_temp$LRmax <- with(results_temp, logRstack + logR_SE)
      results_temp$LRmin <- with(results_temp, logRstack - logR_SE)
      
      # limit data to reversible Reduction and irreversible Reduction
      results_temp <- filter(results_temp, Reduction == "Reversible" | Reduction == "Irreversible")
      # alter the LRmin and LRmax values so that they account for stacking height
      results_temp[results_temp$Reduction == "Reversible", ] <- transform(results_temp[results_temp$Reduction == "Reversible", ],
                                                                          LRmax = LRmax + results_temp[results_temp$Reduction == "Irreversible", "logRstack"],
                                                                          LRmin = LRmin + results_temp[results_temp$Reduction == "Irreversible", "logRstack"])

      # reorder the factors so that Reversible Reduction is shown below Irreversible Reduction
      results_temp$Reduction <- ordered(results_temp$Reduction, c("Reversible", "Irreversible"))
      
      ## compile stacked bar chart
      stacked.bar <- ggplot(results_temp, aes(x = Test, y = logRstack, fill = Reduction)) +
        geom_bar(stat = "identity") +
        facet_grid(.~Pathogen) +
        geom_errorbar(aes(ymax=(LRmax),
                          ymin=(LRmin)),
                      data=results_temp)
      
      if(significance == TRUE){
        stacked.bar <- stacked.bar + 
          geom_text(data = results_temp, size = 8, aes(x=Test,y= logRstack + y_offset, label=t_decision, color = Reduction), position = position_stack())}    
      
      return(stacked.bar)
      }
    else{print('Layout must be either \"staggered" or \"stacked\".', quote = F)}
  }
}
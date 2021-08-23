# partitioning of variance

# import all data and combine to 1 dataframe#####
library(dplyr)
library(stringr)
library(broom)

Merge_Demographic_data <- function() {
  Combine_data <- function(mypath, mypattern, ...) {
    setwd(mypath)
    import.multiple.files <- function(mypath, mypattern, ...)
    {
      setwd(mypath)
      tmp.list.1 <- list.files(mypath, pattern = mypattern)
      tmp.list.2 <- list(length = length(tmp.list.1))
      for (i in 1:length(tmp.list.1)) {
        tmp.list.2[[i]] <- read.csv(tmp.list.1[i], ...)
      }
      names(tmp.list.2) <- tmp.list.1
      tmp.list.2
    }
    files <-
      import.multiple.files(mypath = mypath,
                            mypattern = mypattern,
                            sep = ",")
    
    dataframe <- bind_rows(files)
    return(dataframe)
  }
  
  
  Indi <-
    Combine_data(mypath = "D:/Models_Aug2020/data_analysis/Individual_level",
                 mypattern = "indi.csv",
                 sep = ",")
  Isolated_patch <-
    Combine_data(mypath = "D:/Models_Aug2020/data_analysis/Isolated_patch",
                 mypattern = ".csv",
                 sep = ",")
  pop_size <-
    Combine_data(mypath = "D:/Models_Aug2020/data_analysis/pop_size",
                 mypattern = ".csv",
                 sep = ",")
  
  
  Batch <-
    read.csv("D:/Models_Aug2020/data_analysis/Models/B.csv", header = T)
  Land <-
    read.csv("D:/Models_Aug2020/data_analysis/Models/Land.csv",
             header = T)
  Landscape_analysis <-
    read.csv("D:/Models_Aug2020/data_analysis/Models/Landscape_analysis.csv",
             header = T)
  Simulation <-
    read.csv("D:/Models_Aug2020/data_analysis/Models/sims.csv",
             header = T)
  
  
  
  mergeCols <- c("Simulation", "Batch", "Land", "Rep", "Year")
  Master_data <-
    merge(Isolated_patch, Indi, by = mergeCols, all = TRUE)
  Master_data <-
    merge(Master_data, pop_size, by = mergeCols, all = TRUE)
  Master_data <- merge(Master_data, Batch, by = "Batch", all = TRUE)
  Master_data <-
    merge(Master_data, Simulation, by = "Simulation", all = TRUE)
  Master_data <- merge(Master_data, Land, by = "Land", all = TRUE)
  Master_data <-
    merge(Master_data,
          Landscape_analysis,
          by = "grid",
          all = TRUE)
  return(Master_data)
}
Demographic <- Merge_Demographic_data()
Demographic$Tree_disease_factor <- factor(Demographic$Tree_disease)
Demographic$Management_factor <- factor(Demographic$Management)

Merge_genetic_data <- function(IBD) {
  Combine_data <- function(mypath, mypattern, ...) {
    setwd(mypath)
    library(dplyr)
    import.multiple.files <- function(mypath, mypattern, ...)
    {
      setwd(mypath)
      tmp.list.1 <- list.files(mypath, pattern = mypattern)
      tmp.list.2 <- list(length = length(tmp.list.1))
      for (i in 1:length(tmp.list.1)) {
        tmp.list.2[[i]] <- read.csv(tmp.list.1[i], ...)
      }
      names(tmp.list.2) <- tmp.list.1
      tmp.list.2
    }
    files <-
      import.multiple.files(mypath = mypath,
                            mypattern = mypattern,
                            sep = ",")
    
    dataframe <- bind_rows(files)
    return(dataframe)
  }
  
  
  Genetic_data <-
    Combine_data(mypath = "D:/Models_Aug2020/data_analysis/Genetics",
                 mypattern = "Genetics.csv",
                 sep = ",") #combines genetic outputs from RS
  IBD <-
    Combine_data(mypath = IBD,
                 mypattern = "IBD.csv",
                 sep = ",") #combines IBD estimates
  
  Batch <-
    read.csv("D:/Models_Aug2020/data_analysis/Models/Batch.csv", header = T) #file with details of batches - each batch is different landscape
  Land <-
    read.csv("D:/Models_Aug2020/data_analysis/Models/Land.csv",
             header = T) #landscape treatment information - described by land number 
  Landscape_analysis <-
    read.csv("D:/Models_Aug2020/data_analysis/Models/Landscape_analysis.csv",
             header = T) #analysis of landscape traits 
  Simulation <-
    read.csv("D:/Models_Aug2020/data_analysis/Models/sims.csv",
             header = T) #species traits - change in
  
  
  mergeCols <- c("Simulation", "Batch", "Land", "Rep", "Year")
  Master_data <-
    merge(Genetic_data, CC, by = mergeCols, all = TRUE)
  Master_data <- merge(Master_data, Batch, by = "Batch", all = TRUE)
  Master_data <-
    merge(Master_data, Simulation, by = "Simulation", all = TRUE)
  Master_data <- merge(Master_data, Land, by = "Land", all = TRUE)
  Master_data <-
    merge(Master_data,
          Landscape_analysis,
          by = "grid",
          all = TRUE)
  return(Master_data)
}
Genetic <- Merge_genetic_data(IBD = "D:/Models_Outputs/data_analysis/IBD_fst")
Genetic_Josts <- Merge_genetic_data(IBD = "D:/Models_Outputs/data_analysis/IBD_Jost") # specifiy folder where IBD estimates are stored
  
### proportion of successful disperses ####

Demographic_dispersers <-
  subset(Demographic, Demographic$DistMoved > 0)
Demographic_dispersers <-
  Demographic_dispersers %>% group_by(Simulation, Batch, Land, Rep, Year) %>%
  mutate(proportion = Nind / sum(Nind)) %>% #calculates this by group then ungroups dataframe
  ungroup
Demographic_dispersers_summary <-
  Demographic_dispersers %>% group_by(Simulation, Batch, Land, Rep, Year, Success) %>%
  summarise(proportion_success = sum(proportion))
Demographic_dispersers_summary <-
  as.data.frame(Demographic_dispersers_summary)
Demographic_dispersers <-
  merge(
    Demographic_dispersers,
    Demographic_dispersers_summary,
    by = c("Simulation", "Batch", "Land", "Rep", "Year", "Success"),
    all = TRUE
  )
Demographic_dispersers <-
  subset(Demographic_dispersers,
         Demographic_dispersers$Success == "Yes") #data frame shows proportion successful dispersers by group



##  partitioning of variance ####

# Fit a single anova for the response variable, and the fixed effects should be (in order):
#   Disease level (D)
# Mgt response level (M)
# D * M
# K
# DP
# HM
# K * DP
# K * HM
# DP * HM
# D * K
# D * DP
# D * HM
# M * K
# M * DP
# M * HM
# 10km ID (Square)
# 5km ID (grid)
# 10km ID * D
# 10km ID * M
# 5km ID * D
# 5km ID * M
# 5km ID * Landscape replicate 
# 5km ID * Removal replicate
# 
# The terms from '10km ID' onwards are in reality the random effects in the model, but for the purposes of partitioning variance they should be fitted as fixed effects.
# 


### disperses and isolated patches ####


Demographic_dispersers$K <- factor(Demographic_dispersers$K)
Demographic_dispersers$DP <- factor(Demographic_dispersers$DP)
Demographic_dispersers$HM <- factor(Demographic_dispersers$HM)
Demographic_dispersers$grid <- factor(Demographic_dispersers$grid)
Demographic_dispersers$Square <- factor(Demographic_dispersers$Square)
Demographic_dispersers$Ash_rep <- factor(Demographic_dispersers$Ash_rep)
Demographic_dispersers$Removal_rep <- factor(Demographic_dispersers$Removal_rep)

is_lm <-  anova(lm(proportion_isolated_NFI ~ Tree_disease + Management + Tree_disease * Management + K + DP + HM + K*DP + K*HM + HM*DP + Tree_disease * K + Tree_disease * DP + Tree_disease * HM + Management * K + Management * DP + Management * HM + grid + Square + grid *Tree_disease + Square *Tree_disease + grid *Management + Square*Management + Square * Ash_rep + Square * Removal_rep, data = Demographic_dispersers))

dis_lm <-  anova(lm(proportion_success ~ Tree_disease + Management + Tree_disease * Management + K + DP + HM + K*DP + K*HM + HM*DP + Tree_disease * K + Tree_disease * DP + Tree_disease * HM + Management * K + Management * DP + Management * HM + grid + Square + grid *Tree_disease + Square *Tree_disease + grid *Management + Square*Management + Square * Ash_rep + Square * Removal_rep, data = Demographic_dispersers))




# genetic data ####

### Josts D
gen_summary <- subset(Genetic_Josts, Genetic_Josts$model == "lm")
gen_summary <- subset(gen_summary, gen_summary$term == "Distance_km")
gen_summary$K <- factor(gen_summary$K)
gen_summary$DP <- factor(gen_summary$DP)
gen_summary$HM <- factor(gen_summary$HM)
gen_summary$grid <- factor(gen_summary$grid)
gen_summary$Square <- factor(gen_summary$Square)
gen_summary$Ash_rep <- factor(gen_summary$Ash_rep)
gen_summary$Removal_rep <- factor(gen_summary$Removal_rep)
gen_summary <- gen_summary %>% filter(!is.na(estimate))
gen_summary <- subset(gen_summary, gen_summary$df.residual > 4)

gen_lm <-
  lm(estimate ~ Tree_disease + Management + Tree_disease * Management + K + DP + HM + K*DP + K*HM + HM*DP + Tree_disease * K + Tree_disease * DP + Tree_disease * HM + Management * K + Management * DP + Management * HM + grid + Square + grid *Tree_disease + Square *Tree_disease + grid *Management + Square*Management + Square * Ash_rep + Square * Removal_rep, data = gen_summary)
anova(gen_lm) #lm across all landscape - reponse is the IBD estimate 

gen_lm_5k  <- gen_summary %>%
  group_by(grid) %>%            # group by column z
  do(tidy(anova(lm(estimate ~ Tree_disease + Management + Tree_disease * Management + K + DP + HM + K*DP + K*HM + HM*DP + Tree_disease * K + Tree_disease * DP + Tree_disease * HM + Management * K + Management * DP + Management * HM + bearing *Tree_disease +  bearing *Management + bearing *Ash_rep + bearing * Removal_rep, data =.)))) #linear model for each 5 km grid 

gen_lm_5k <-
  gen_lm_5k %>% group_by(grid) %>%
  mutate(percentage_variance = sumsq / sum(sumsq)) %>% 
  ungroup
gen_lm_5k$percentage_variance <- gen_lm_5k$percentage_variance *100 # calculate % variance 


### FST
gen_summary <- subset(Genetic, Genetic$model == "lm")
gen_summary <- subset(gen_summary, gen_summary$term == "Distance")
gen_summary$K <- factor(gen_summary$K)
gen_summary$DP <- factor(gen_summary$DP)
gen_summary$HM <- factor(gen_summary$HM)
gen_summary$grid <- factor(gen_summary$grid)
gen_summary$Square <- factor(gen_summary$Square)
gen_summary$Ash_rep <- factor(gen_summary$Ash_rep)
gen_summary$Removal_rep <- factor(gen_summary$Removal_rep)
gen_summary <- gen_summary %>% filter(!is.na(estimate))
gen_summary <- subset(gen_summary, gen_summary$df.residual > 4)

gen_lm <-
  lm(estimate ~ Tree_disease + Management + Tree_disease * Management + K + DP + HM + K*DP + K*HM + HM*DP + Tree_disease * K + Tree_disease * DP + Tree_disease * HM + Management * K + Management * DP + Management * HM + grid + Square + grid *Tree_disease + Square *Tree_disease + grid *Management + Square*Management + Square * Ash_rep + Square * Removal_rep, data = gen_summary)
anova(gen_lm) #lm across all landscape - reponse is the IBD estimate 

gen_lm_5k  <- gen_summary %>%
  group_by(grid) %>%            # group by column z
  do(tidy(anova(lm(estimate ~ Tree_disease + Management + Tree_disease * Management + K + DP + HM + K*DP + K*HM + HM*DP + Tree_disease * K + Tree_disease * DP + Tree_disease * HM + Management * K + Management * DP + Management * HM + bearing *Tree_disease +  bearing *Management + bearing *Ash_rep + bearing * Removal_rep, data =.)))) #linear model for each 5 km grid 

gen_lm_5k <-
  gen_lm_5k %>% group_by(grid) %>%
  mutate(percentage_variance = sumsq / sum(sumsq)) %>% 
  ungroup
gen_lm_5k$percentage_variance <- gen_lm_5k$percentage_variance *100 # calculate % variance 

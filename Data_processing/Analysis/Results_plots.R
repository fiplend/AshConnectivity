# code for creating summary figures for results sections
library(ggplot2)
library(dplyr)
library(plyr)
library(ggpubr)
library(data.table)

# Dispersal data - impacts to functional connectvity ####
# Combine demographic summary data into 1 dataset 
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
          by = "Square",
          all = TRUE)
  return(Master_data)
}
Demographic <- Merge_Demographic_data()

Demographic$Tree_disease_factor <- factor(Demographic$Tree_disease)
Demographic$Management_factor <- factor(Demographic$Management)

# population size is calculated before demographic deaths therefore we need to summarize number of dispersers and use this to calculate percentages

Demographic_dispersers <-
  subset(Demographic, Demographic$DistMoved > 0)
Demographic_dispersers_summary <-
  Demographic_dispersers %>% 
  group_by(Simulation, Batch, Land, Rep, Year) %>%
  dplyr::summarise(Sum_Nind = sum(Nind))
Demographic_dispersers_sum <-
  Demographic_dispersers %>% 
  group_by(Simulation, Batch, Land, Rep, Year, Success) %>%
  dplyr::summarise(Sum_success = sum(Nind))
Demographic_dispersers_sums <-
  merge(
    Demographic_dispersers_sum,
    Demographic_dispersers_summary,
    by = c("Simulation", "Batch", "Land", "Rep", "Year"),
    all = TRUE
  )
Demographic_dispersers <-
  merge(
    Demographic_dispersers,
    Demographic_dispersers_sums,
    by = c("Simulation", "Batch", "Land", "Rep", "Year", "Success"),
    all = TRUE
  )
Demographic_dispersers$proportion_success <- Demographic_dispersers$Sum_success / Demographic_dispersers$Sum_Nind

Demographic_dispersers <-
  as.data.frame(Demographic_dispersers)
Demographic_dispersers <-
  subset(Demographic_dispersers,
         Demographic_dispersers$Success == "Yes")
Demographic_dispersers$proportion_isolated_NFI[is.na(Demographic_dispersers$proportion_isolated_NFI)] <-
  0
Demographic_dispersers$landrep <-
  paste(Demographic_dispersers$Ash_rep,
        Demographic_dispersers$Removal_rep,
        sep = "_")
Demographic_dispersers$K <- factor(Demographic_dispersers$K)
Demographic_dispersers$HM <- factor(Demographic_dispersers$HM)
Demographic_dispersers$DP <- factor(Demographic_dispersers$DP)
Demographic_dispersers$Rep <- factor(Demographic_dispersers$Rep)
Demographic_dispersers$Year <- factor(Demographic_dispersers$Year)


# data needs to be summarised in order to make plots 
#### Isolated plots ####
Demographic_dispersers_NA <-
  Demographic_dispersers %>% filter(!is.na(proportion_isolated_NFI))
Demographic_dispersers_NA$DP <- factor(Demographic_dispersers_NA$DP)
Demographic_dispersers_NA$K <- factor(Demographic_dispersers_NA$K)
Demographic_dispersers_NA$HM <- factor(Demographic_dispersers_NA$HM)
Demographic_dispersers_NA$Management_factor <-
  factor(Demographic_dispersers_NA$Management)

DP_iso <-
  ddply(
    Demographic_dispersers_NA ,
    c("Tree_disease", "DP", "Management_factor"),
    summarise,
    N_iso    = length(proportion_isolated_NFI),
    mean_proportion_isolated_NFI = mean(proportion_isolated_NFI),
    sd_iso   = sd(proportion_isolated_NFI),
    se_iso   = sd_iso / sqrt(N_iso)
  )
HM_iso <-
  ddply(
    Demographic_dispersers_NA ,
    c("Tree_disease", "HM", "Management_factor"),
    summarise,
    N_iso    = length(proportion_isolated_NFI),
    mean_proportion_isolated_NFI = mean(proportion_isolated_NFI),
    sd_iso   = sd(proportion_isolated_NFI),
    se_iso   = sd_iso / sqrt(N_iso)
  )
K_iso <-
  ddply(
    Demographic_dispersers_NA ,
    c("Tree_disease", "K", "Management_factor"),
    summarise,
    N_iso    = length(proportion_isolated_NFI),
    mean_proportion_isolated_NFI = mean(proportion_isolated_NFI),
    sd_iso   = sd(proportion_isolated_NFI),
    se_iso   = sd_iso / sqrt(N_iso)
  )

#### dispersal proportion ####

Demographic_dispersers_NA <-
  Demographic_dispersers %>% filter(!is.na(proportion_success))
Demographic_dispersers_NA$DP <- factor(Demographic_dispersers_NA$DP)
Demographic_dispersers_NA$K <- factor(Demographic_dispersers_NA$K)
Demographic_dispersers_NA$HM <- factor(Demographic_dispersers_NA$HM)
Demographic_dispersers_NA$Management_factor <-
  factor(Demographic_dispersers_NA$Management)

DP_dispersal <-
  ddply(
    Demographic_dispersers_NA ,
    c("Tree_disease", "DP", "Management_factor"),
    summarise,
    N_disp    = length(proportion_success),
    mean_proportion_success = mean(proportion_success),
    sd_disp   = sd(proportion_success),
    se_disp   = sd_disp / sqrt(N_disp)
  )

HM_dispersal <-
  ddply(
    Demographic_dispersers_NA ,
    c("Tree_disease", "HM", "Management_factor"),
    summarise,
    N_disp    = length(proportion_success),
    mean_proportion_success = mean(proportion_success),
    sd_disp   = sd(proportion_success),
    se_disp   = sd_disp / sqrt(N_disp)
  )
K_dispersal <-
  ddply(
    Demographic_dispersers_NA ,
    c("Tree_disease", "K", "Management_factor"),
    summarise,
    N_disp    = length(proportion_success),
    mean_proportion_success = mean(proportion_success),
    sd_disp   = sd(proportion_success),
    se_disp   = sd_disp / sqrt(N_disp)
  )


#### summary data - species traits data ####
HM_dispersal$Species_trait <- "Habitat mortality"
HM_iso$Species_trait <- "Habitat mortality"
DP_dispersal$Species_trait <- "Directional persistence"
DP_iso$Species_trait <- "Directional persistence"
K_dispersal$Species_trait <- "Carrying capacity"
K_iso$Species_trait <- "Carrying capacity"

names(HM_dispersal)[2] <- "Trait"
names(HM_iso)[2] <- "Trait"
names(DP_dispersal)[2] <- "Trait"
names(DP_iso)[2] <- "Trait"
names(K_dispersal)[2] <- "Trait"
names(K_iso)[2] <- "Trait"

HM_dispersal$Trait_value <- NA
HM_dispersal$Trait_value <- ifelse(
  HM_dispersal$Trait == 0.02,
  'low',
  ifelse(
    HM_dispersal$Trait == 0.035,
    'medium',
    ifelse(HM_dispersal$Trait == 0.05, 'high', 'NA')
  )
)
HM_iso$Trait_value <- NA
HM_iso$Trait_value <- ifelse(
  HM_iso$Trait == 0.02,
  'low',
  ifelse(
    HM_iso$Trait == 0.035,
    'medium',
    ifelse(HM_iso$Trait == 0.05, 'high', 'NA')
  )
)    

DP_dispersal$Trait_value <- NA
DP_dispersal$Trait_value <- ifelse(
  DP_dispersal$Trait == 5,
  'low',
  ifelse(
    DP_dispersal$Trait == 7,
    'medium',
    ifelse(DP_dispersal$Trait == 9, 'high', 'NA')
  )
)
DP_iso$Trait_value <- NA
DP_iso$Trait_value <- ifelse(
  DP_iso$Trait == 5,
  'low',
  ifelse(
    DP_iso$Trait == 7,
    'medium',
    ifelse(DP_iso$Trait == 9, 'high', 'NA')
  )
)     

K_dispersal$Trait_value <- NA
K_dispersal$Trait_value <- ifelse(
  K_dispersal$Trait == 500,
  'low',
  ifelse(
    K_dispersal$Trait == 750,
    'medium',
    ifelse(K_dispersal$Trait == 1000, 'high', 'NA')
  )
)
K_iso$Trait_value <- NA
K_iso$Trait_value <- ifelse(
  K_iso$Trait == 500,
  'low',
  ifelse(
    K_iso$Trait == 750,
    'medium',
    ifelse(K_iso$Trait == 1000, 'high', 'NA')
  )
)     

dispersal_summary <- rbind(HM_dispersal, DP_dispersal, K_dispersal)
iso_summary <- rbind(HM_iso, DP_iso, K_iso)
Disp_data <-
  merge(
    iso_summary,
    dispersal_summary,
    by = c("Trait", "Tree_disease", "Trait_value", "Management_factor", "Species_trait")
  )

Disp_data$Trait_value <- factor(Disp_data$Trait_value, levels = c("low", "medium", "high"))

# averaged plots ####

Summary <-
  ddply(
    Disp_data ,
    c("Tree_disease", "Management_factor"),
    summarise,
    N_disp    = length(N_disp),
    mean_proportion_success = mean(mean_proportion_success),
    sd_disp   = sd(sd_disp),
    se_disp   = sd_disp / sqrt(N_disp),
    N_iso    = length(N_iso),
    mean_proportion_isolated_NFI = mean(mean_proportion_isolated_NFI),
    sd_iso   = sd(sd_iso),
    se_iso   = sd_iso / sqrt(N_iso)
  )

Dispersal_succes_sum <-
  ggplot(Summary,
         aes(x = Tree_disease, y = mean_proportion_success)) +
  geom_point(aes(shape = Management_factor), size = 1.5) +
  geom_line(aes(linetype = Management_factor), size = 1) +
  ylab("Proportion succesful dispersers") +
  xlab("Level of tree disease (%)") +
  guides(shape=guide_legend("Management effort (%)"),
         linetype=guide_legend("Management effort (%)"))+
  plot_theme

Isolated_woods_sum  <-
  ggplot(Summary,
         aes(x = Tree_disease, y = mean_proportion_isolated_NFI)) +
  geom_point(aes(shape = Management_factor), size = 1.5) +
  geom_line(aes(linetype = Management_factor), size = 1) +
  ylab("Proportion isolated woodlands") +
  xlab("Level of tree disease (%)") +
  guides(shape=guide_legend("Management effort (%)"),
         linetype=guide_legend("Management effort (%)"))+
  plot_theme

ggarrange(
  Dispersal_succes_sum,
  Isolated_woods_sum,
  nrow =  1,
  common.legend = TRUE,
  legend = "right",
  heights = c(1, 1),
  labels = "AUTO"
)

# species traits plots ####

Dispersal_succes <-
  ggplot(Disp_data,
         aes(x = Tree_disease, y = mean_proportion_success, colour = Trait_value)) +
  geom_point(aes(shape = Management_factor), size = 1.5) +
  geom_line(aes(linetype = Management_factor), size = 1) +
  facet_wrap( ~ Species_trait) +
  ylab("Proportion succesful dispersers") +
  xlab("Level of tree disease (%)") +
  scale_color_manual(values = c("#BA7B1D", "#1E9C2E", "#306CA1"), name = "Trait value")+
  guides(col=guide_legend("Trait value"),
         shape=guide_legend("Management effort (%)"),
         linetype=guide_legend("Management effort (%)"))+
  plot_theme

Isolated_woods <-
  ggplot(Disp_data,
         aes(x = Tree_disease, y = mean_proportion_isolated_NFI, colour = Trait_value)) +
  geom_point(aes(shape = Management_factor), size = 1.5) +
  geom_line(aes(linetype = Management_factor), size = 1) +
  facet_wrap( ~ Species_trait) +
  ylab("Proportion isolated woodlands") +
  xlab("Level of tree disease (%)") +
  scale_color_manual(values = c("#BA7B1D", "#1E9C2E", "#306CA1"), name = "Trait value")+
  guides(col=guide_legend("Trait value"),
         shape=guide_legend("Management effort (%)"),
         linetype=guide_legend("Management effort (%)"))+
  plot_theme


ggarrange(
  Dispersal_succes + rremove("x.text") + rremove("xlab") + rremove("x.ticks"),
  Isolated_woods,
  nrow =  2,
  common.legend = TRUE,
  legend = "right",
  heights = c(1, 1, 1.25),
  labels = "AUTO"
)

### woodland disperses graph ####
Demographic_disp <-
  subset(Demographic, Demographic$DistMoved > 0)
#create column to store woodland info
Demographic_disp$Woodland_type <- ifelse(
  Demographic_disp$Natal_patch_type == "Woodland" & Demographic_disp$Focal_patch_type == "Woodland",
  'Woodland',
  ifelse(
    Demographic_disp$Focal_patch_type == "Matrix", "Matrix", 'Steppingstone')
)
Demographic_disp <-
  subset(Demographic_disp,
         Demographic_disp$Success == "Yes")
Demographic_disp_summary <-
  Demographic_disp %>% group_by(Simulation, Batch, Land, Rep, Year, Woodland_type, Success) %>%
  dplyr::summarise(mean_dist = mean(DistMoved), 
                   Nind = sum(Nind))
Demographic_disp_summary <-
  as.data.frame(Demographic_disp_summary)


# combine summary data with batch data 
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
Dispersal_distance <- merge(Demographic_disp_summary, Batch, by = "Batch", all = TRUE)
Dispersal_distance <-
  merge(Dispersal_distance, Simulation, by = "Simulation", all = TRUE)
Dispersal_distance <- merge(Dispersal_distance, Land, by = "Land", all = TRUE)
Dispersal_distance <-
  merge(Dispersal_distance,
        Landscape_analysis,
        by = "Square",
        all = TRUE)
Dispersal_distance$Management_factor <- factor(Dispersal_distance$Management)

#summarise data for plots 
Dispersal_distance_summary <-
  ddply(
    Dispersal_distance,
    c("Tree_disease", "Management_factor", "Woodland_type", "K", "HM"),
    summarise,
    mean_proportion_success = mean(proportion_success),
    mean_distance  = mean(mean_dist),
    mean_Nind = mean(Nind), 
  )

Dispersal_distance_summary_woodland <- subset(Dispersal_distance_summary, Dispersal_distance_summary$Woodland_type == "Woodland")

Plot_1 <-
  ggplot(
    Dispersal_distance_summary_woodland,
    aes(x = Tree_disease, y = mean_proportion_success)
  ) +
  geom_point(aes(shape = Management_factor), size = 1.5) +
  geom_line(aes(linetype = Management_factor), size = 1) +
  ylab("Mean number of individuals") +
  xlab("Level of tree disease (%)") +
  scale_color_manual(values = c("#BA7B1D", "#1E9C2E"), name = "Patch type") +
  guides(
    col = guide_legend("Patch type"),
    shape = guide_legend("Management effort (%)"),
    linetype = guide_legend("Management effort (%)")
  ) +
  plot_theme
Plot_2 <-
  ggplot(
    Dispersal_distance_summary,
    aes(x = Tree_disease, y = mean_distance)
  ) +
  geom_point(aes(shape = Management_factor), size = 1.5) +
  geom_line(aes(linetype = Management_factor), size = 1) +
  ylab("Mean distance traveled (m)") +
  xlab("Level of tree disease (%)") +
  scale_color_manual(values = c("#BA7B1D", "#1E9C2E"), name = "Patch type") +
  guides(
    col = guide_legend("Patch type"),
    shape = guide_legend("Management effort (%)"),
    linetype = guide_legend("Management effort (%)")
  ) +
  plot_theme

ggarrange(
  Plot_1 + rremove("x.text") + rremove("xlab") + rremove("x.ticks"),
  Plot_2,
  nrow = 2,
  common.legend = T,
  legend = "right",
  labels = "AUTO",
  heights = c(.48, .52))



# supplementary fig - tree cover ####

Demographic_dispersers_NA <-
  Demographic_dispersers %>% filter(!is.na(proportion_isolated_NFI))
Demographic_dispersers_NA$DP <- factor(Demographic_dispersers_NA$DP)
Demographic_dispersers_NA$K <- factor(Demographic_dispersers_NA$K)
Demographic_dispersers_NA$HM <- factor(Demographic_dispersers_NA$HM)
Demographic_dispersers_NA$Management_factor <-
  factor(Demographic_dispersers_NA$Management_factor)

Dispersal <-
  ddply(Demographic_dispersers_NA,
        c("Tree_disease", "Ash_rep", "Management_factor", "grid", "bearing"),
        summarise,
        N_disp    = length(proportion_success),
        mean_proportion_success = mean(proportion_success),
        sd_disp   = sd(proportion_success),
        se_disp   = sd_disp / sqrt(N_disp),
        mean_isolated = mean(proportion_isolated_NFI),
        mean_tree = mean(percent_tree_cover),
        mean_woodland = mean(percent_woodland_tree)
  )
plot_1 <- ggplot(Dispersal,
                 aes(x = Tree_disease, y = mean_isolated, colour = mean_tree)) +
  geom_point(size = 3) +
  scale_colour_viridis()+
  plot_theme+
  labs(colour = "Tree cover (%)")+
  ylab("Proportion isolated woodlands") +
  xlab("Level of tree disease (%)")

Demographic_dispersers_NA_2 <-
  Demographic_dispersers %>% filter(!is.na(proportion_success))
Demographic_dispersers_NA_2$DP <- factor(Demographic_dispersers_NA$DP)
Demographic_dispersers_NA_2$K <- factor(Demographic_dispersers_NA$K)
Demographic_dispersers_NA_2$HM <- factor(Demographic_dispersers_NA$HM)
Demographic_dispersers_NA_2$Management_factor <-
  factor(Demographic_dispersers_NA$Management_factor)

Dispersal_2 <-
  ddply(Demographic_dispersers_NA_2,
        c("Tree_disease", "Ash_rep", "Management_factor", "grid", "bearing"),
        summarise,
        N_disp    = length(proportion_success),
        mean_proportion_success = mean(proportion_success),
        sd_disp   = sd(proportion_success),
        se_disp   = sd_disp / sqrt(N_disp),
        mean_isolated = mean(proportion_isolated_NFI),
        mean_tree = mean(percent_tree_cover),
        mean_woodland = mean(percent_woodland_tree)
  )
plot_2 <- ggplot(Dispersal_2,
                 aes(x = Tree_disease, y = mean_proportion_success, colour = mean_tree)) +
  geom_point(size = 3) +
  scale_colour_viridis()+
  plot_theme+
  labs(colour = "Tree cover (%)")+
  ylab("Proportion succesful dispersers") +
  xlab("Level of tree disease (%)")+ 
  ylim(0.1,0.5)

ggarrange(plot_2 + rremove("x.text") + rremove("xlab") + rremove("x.ticks"),
          plot_1,
          nrow =  2,
          common.legend = TRUE,
          legend = "right",
          heights = c(1, 1, 1.25),
          labels = "AUTO")



# genetic data #####

Genetic <-
  read.csv("D:/Models_Aug2020/data_analysis/Genetic_data.csv", header = T)
Genetic_Jost <-
  read.csv("D:/Models_Aug2020/data_analysis/Genetic_Josts.csv", header = T)


Merge_genetic_data_Jost <- function() {
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
  
  
  Genetic_data <-
    Combine_data(mypath = "D:/Models_Aug2020/data_analysis/Genetics",
                 mypattern = "Genetics.csv",
                 sep = ",")
  CC <-
    Combine_data(mypath = "D:/Models_Aug2020/data_analysis/Jost",
                 mypattern = "CC.csv",
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
  
  Master_data <- merge(Genetic_data, CC, by = mergeCols, all = TRUE)
  Master_data <- merge(Master_data, Batch, by = "Batch", all = TRUE)
  Master_data <-
    merge(Master_data, Simulation, by = "Simulation", all = TRUE)
  Master_data <- merge(Master_data, Land, by = "Land", all = TRUE)
  Master_data <-
    merge(Master_data,
          Landscape_analysis,
          by = "Square",
          all = TRUE)
  return(Master_data)
}
Genetic_Jost <- Merge_genetic_data_Jost()
Merge_genetic_data <- function() {
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
  
  
  Genetic_data <-
    Combine_data(mypath = "D:/Models_Aug2020/data_analysis/Genetics",
                 mypattern = "Genetics.csv",
                 sep = ",")
  iso_plot <-
    Combine_data(mypath = "D:/Models_Aug2020/data_analysis/iso_plot",
                 mypattern = "iso_plot.csv",
                 sep = ",")
  CC <-
    Combine_data(mypath = "D:/Models_Aug2020/data_analysis/CC",
                 mypattern = "CC.csv",
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
    merge(Genetic_data, iso_plot, by = mergeCols, all = TRUE)
  Master_data <- merge(Master_data, CC, by = mergeCols, all = TRUE)
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
Genetic <- Merge_genetic_data()

Genetic$IBD <- "FST"
Genetic_Jost$IBD <- "Jost"
Genetic<- subset(Genetic, Genetic$model == "lm")
Genetic <- subset(Genetic, Genetic$term == "Distance")
Genetic_Jost<- subset(Genetic_Jost, Genetic_Jost$model == "lm")
Genetic_Jost <- subset(Genetic_Jost, Genetic_Jost$term == "Distance_km")
Genetic <- Genetic %>% filter(!is.na(estimate))
Genetic <- subset(Genetic, Genetic$df.residual > 4)
Genetic$estimate <- Genetic$estimate * 1000
Genetic_Jost <- Genetic_Jost %>% filter(!is.na(estimate))
Genetic_Jost <- subset(Genetic_Jost, Genetic_Jost$df.residual > 4)
Genetic$COR_FST <- NULL
Genetic <- rbind(Genetic, Genetic_Jost)

Genetic$K <- factor(Genetic$K)
Genetic$HM <- factor(Genetic$HM)
Genetic$DP <- factor(Genetic$DP)
Genetic$Management_factor <- factor(Genetic$Management)
Genetic_Jost$Management_factor <- factor(Genetic_Jost$Management)


DP_IBD <-
  ddply(
    Genetic,
    c("Tree_disease", "DP", "Management_factor", "IBD"),
    summarise,
    N_cor    = length(estimate),
    mean_cor = mean(estimate),
    sd_cor   = sd(estimate),
    se_cor   = sd_cor / sqrt(N_cor)
  )
HM_IBD <-
  ddply(
    Genetic,
    c("Tree_disease", "HM", "Management_factor", "IBD"),
    summarise,
    N_cor    = length(estimate),
    mean_cor = mean(estimate),
    sd_cor   = sd(estimate),
    se_cor   = sd_cor / sqrt(N_cor)
  )
K_IBD <-
  ddply(
    Genetic,
    c("Tree_disease", "K", "Management_factor", "IBD"),
    summarise,
    N_cor    = length(estimate),
    mean_cor = mean(estimate),
    sd_cor   = sd(estimate),
    se_cor   = sd_cor / sqrt(N_cor)
  )

HM_IBD$Species_trait <- "Habitat mortality"
DP_IBD$Species_trait <- "Directional persistence"
K_IBD$Species_trait <- "Carrying capacity"

names(HM_IBD)[2] <- "Trait"
names(DP_IBD)[2] <- "Trait"
names(K_IBD)[2] <- "Trait"
HM_IBD$Trait_value <- NA
HM_IBD$Trait_value <- ifelse(
  HM_IBD$Trait == 0.02,
  'low',
  ifelse(
    HM_IBD$Trait == 0.035,
    'medium',
    ifelse(HM_IBD$Trait == 0.05, 'high', 'NA')
  )
)    
DP_IBD$Trait_value <- NA
DP_IBD$Trait_value <- ifelse(
  DP_IBD$Trait == 5,
  'low',
  ifelse(
    DP_IBD$Trait == 7,
    'medium',
    ifelse(DP_IBD$Trait == 9, 'high', 'NA')
  )
)    
K_IBD$Trait_value <- NA
K_IBD$Trait_value <- ifelse(
  K_IBD$Trait == 500,
  'low',
  ifelse(
    K_IBD$Trait == 750,
    'medium',
    ifelse(K_IBD$Trait == 1000, 'high', 'NA')
  )
)    
IBD_summary <- rbind(HM_IBD, DP_IBD, K_IBD)
IBD_summary$Trait_value <- factor(IBD_summary$Trait_value, levels = c("low", "medium", "high"))

IBD_Jost <- subset(IBD_summary, IBD_summary$IBD == "Jost")
IBD_jost <-
  ggplot(IBD_Jost, aes(x = Tree_disease, y = mean_cor, colour = Trait_value)) +
  geom_point(aes(shape = Management_factor), size = 1.5) +
  geom_line(aes(linetype = Management_factor), size = 1) +
  facet_wrap(~ Species_trait, scales = "fixed") +
  ylab("Isolation by distance (Jost D)") +
  xlab("Level of tree disease (%)") +
  scale_color_manual(values = c("#BA7B1D", "#1E9C2E", "#306CA1"), name = "Trait value")+
  guides(col=guide_legend("Trait value"),
         shape=guide_legend("Management effort (%)"),
         linetype=guide_legend("Management effort (%)"))+
  plot_theme
IBD_FST <- subset(IBD_summary, IBD_summary$IBD == "FST")
IBD_F <-
  ggplot(IBD_FST, aes(x = Tree_disease, y = mean_cor, colour = Trait_value)) +
  geom_point(aes(shape = Management_factor), size = 1.5) +
  geom_line(aes(linetype = Management_factor), size = 1) +
  facet_wrap(~ Species_trait, scales = "fixed") +
  ylab("Isolation by distance (FST)") +
  xlab("Level of tree disease (%)") +
  scale_color_manual(values = c("#BA7B1D", "#1E9C2E", "#306CA1"), name = "Trait value")+
  guides(col=guide_legend("Trait value"),
         shape=guide_legend("Management effort (%)"),
         linetype=guide_legend("Management effort (%)"))+
  plot_theme

ggarrange(IBD_jost + rremove("x.text") + rremove("xlab") + rremove("x.ticks"),
          IBD_F,
          nrow =  2,
          common.legend = TRUE,
          legend = "right",
          heights = c(1, 1, 1.25),
          labels = "AUTO")

landscape_IBD <-
  ddply(
    Genetic,
    c("Tree_disease", "Square", "grid", "bearing", "IBD"),
    summarise,
    N_cor    = length(estimate),
    mean_cor = mean(estimate),
    sd_cor   = sd(estimate),
    se_cor   = sd_cor / sqrt(N_cor)
  )



landscape_fst <- subset(landscape_IBD, landscape_IBD$IBD == "FST")
landscape_jost <- subset(landscape_IBD, landscape_IBD$IBD == "Jost")
landscape_fst_TL90 <- subset(landscape_fst, landscape_fst$grid == "TL90")
landscape_jost_TL90 <- subset(landscape_jost, landscape_jost$grid == "TL90")
landscape_fst <- subset(landscape_fst, landscape_fst$grid != "TL90")
landscape_jost <- subset(landscape_jost, landscape_jost$grid != "TL90")



#landscape_IBD <- subset(landscape_IBD, landscape_IBD$grid != "TL90")

landscape_fst$x <- -1
landscape_fst$y <- 0.08
landscape_fst$y2 <- -0.01
landscape_fst$y[landscape_fst$grid =="TL90"] <- 0.5
landscape_fst$grid <- factor(landscape_fst$grid , levels = c("TL90", "TL54", "TL74", "TL96", "TM15", "TM17"))

blank_data <- read.csv("D:/Models_Aug2020/limits.csv", header = TRUE)

Landscape_Fst <- ggplot(landscape_fst,
                        aes(x = Tree_disease, y = mean_cor, colour = bearing)) +
  geom_blank(data = landscape_fst, aes(x = x, y = y)) +
  geom_blank(data = landscape_fst, aes(x = x, y = y2))+
  geom_line(aes(), size = 1) +
  geom_point(aes(), size = 1.5) +
  facet_wrap(~ grid, nrow = 6, scales = "free_y") +
  ylab("Isolation by distance (FST)") +
  xlab("Level of tree disease (%)")+
  plot_theme+
  scale_color_manual(values = c("#BA7B1D", "#1E9C2E", "#306CA1", "#B22222"), name = "Quadrant")


landscape_jost$x <- -1
landscape_jost$y <- 0.06
landscape_jost$y2 <- -0.01
landscape_jost$y[landscape_jost$grid =="TL90"] <- 0.5
landscape_jost$grid <- factor(landscape_jost$grid , levels = c("TL90", "TL54", "TL74", "TL96", "TM15", "TM17"))
Landscape_Jost <-
  ggplot(landscape_jost,
         aes(x = Tree_disease, y = mean_cor, colour = bearing)) +
  geom_blank(data = landscape_jost, aes(x = x, y = y)) +
  geom_blank(data = landscape_jost, aes(x = x, y = y2))+
  geom_line(aes(), size = 1) +
  geom_point(aes(), size = 1.5) +
  facet_wrap(~ grid, nrow = 6, scales = "free_y") +
  ylab("Isolation by distance (Jost D)") +
  xlab("Level of tree disease (%)")+
  plot_theme+
  scale_color_manual(values = c("#BA7B1D", "#1E9C2E", "#306CA1", "#B22222"), name = "Quadrant")



ggarrange(Landscape_Fst,
          Landscape_Jost,
          common.legend = T,
          legend = "right", 
          ncol =2)



# landscape plots - genetic ####


landscape_IBD <-
  ddply(
    Genetic,
    c("Tree_disease", "Square", "grid", "bearing", "IBD"),
    summarise,
    N_cor    = length(estimate),
    mean_cor = mean(estimate),
    sd_cor   = sd(estimate),
    se_cor   = sd_cor / sqrt(N_cor)
  )



landscape_fst <- subset(landscape_IBD, landscape_IBD$IBD == "FST")
landscape_jost <- subset(landscape_IBD, landscape_IBD$IBD == "Jost")
landscape_fst_TL90 <- subset(landscape_fst, landscape_fst$grid == "TL90")
landscape_jost_TL90 <- subset(landscape_jost, landscape_jost$grid == "TL90")
landscape_fst <- subset(landscape_fst, landscape_fst$grid != "TL90")
landscape_jost <- subset(landscape_jost, landscape_jost$grid != "TL90")



#landscape_IBD <- subset(landscape_IBD, landscape_IBD$grid != "TL90")

landscape_fst$x <- -1
landscape_fst$y <- 0.08
landscape_fst$y2 <- -0.01
landscape_fst$y[landscape_fst$grid =="TL90"] <- 0.5
landscape_fst$grid <- factor(landscape_fst$grid , levels = c("TL90", "TL54", "TL74", "TL96", "TM15", "TM17"))

blank_data <- read.csv("D:/Models_Aug2020/limits.csv", header = TRUE)

Landscape_Fst <- ggplot(landscape_fst,
                        aes(x = Tree_disease, y = mean_cor, colour = bearing)) +
  geom_blank(data = landscape_fst, aes(x = x, y = y)) +
  geom_blank(data = landscape_fst, aes(x = x, y = y2))+
  geom_line(aes(), size = 1) +
  geom_point(aes(), size = 1.5) +
  facet_wrap(~ grid, nrow = 6, scales = "free_y") +
  ylab("Isolation by distance (FST)") +
  xlab("Level of tree disease (%)")+
  plot_theme+
  scale_color_manual(values = c("#BA7B1D", "#1E9C2E", "#306CA1", "#B22222"), name = "Quadrant")


landscape_jost$x <- -1
landscape_jost$y <- 0.06
landscape_jost$y2 <- -0.01
landscape_jost$y[landscape_jost$grid =="TL90"] <- 0.5
landscape_jost$grid <- factor(landscape_jost$grid , levels = c("TL90", "TL54", "TL74", "TL96", "TM15", "TM17"))
Landscape_Jost <-
  ggplot(landscape_jost,
         aes(x = Tree_disease, y = mean_cor, colour = bearing)) +
  geom_blank(data = landscape_jost, aes(x = x, y = y)) +
  geom_blank(data = landscape_jost, aes(x = x, y = y2))+
  geom_line(aes(), size = 1) +
  geom_point(aes(), size = 1.5) +
  facet_wrap(~ grid, nrow = 6, scales = "free_y") +
  ylab("Isolation by distance (Jost D)") +
  xlab("Level of tree disease (%)")+
  plot_theme+
  scale_color_manual(values = c("#BA7B1D", "#1E9C2E", "#306CA1", "#B22222"), name = "Quadrant")



ggarrange(Landscape_Fst,
          Landscape_Jost,
          common.legend = T,
          legend = "right", 
          ncol =2)



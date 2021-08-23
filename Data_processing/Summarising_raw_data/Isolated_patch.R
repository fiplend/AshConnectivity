# Number of isolated patches ####
# conditions - if number ind incoming = 0 for both years
# and then for NFI so size of patch <200 cells

library(stringr)
library(raster)
Isolated_patch <- function(inputpath, patch, mypath, mypattern, ...) {
  setwd(mypath)
  library(dplyr)
  import.multiple.files<-function(mypath,mypattern,...)
  {
    setwd(mypath)
    tmp.list.1<-list.files(mypath, pattern=mypattern)
    tmp.list.2<-list(length=length(tmp.list.1))
    for (i in 1:length(tmp.list.1)){tmp.list.2[[i]]<-read.csv(tmp.list.1[i],...)}
    names(tmp.list.2)<-tmp.list.1
    tmp.list.2
  }
  files <-
    import.multiple.files(mypath = mypath,
                          mypattern = mypattern,
                          sep = "\t")
  
  filenames <-
    list.files(mypath, pattern = mypattern)
  filenames <- as.data.frame(filenames)
  filenames <-
    as.data.frame(str_split_fixed(filenames$filenames, "_", 4))
  Batch <- unname(unlist(filenames["V1"]))
  Sim <- unname(unlist(filenames["V2"]))
  Land <- unname(unlist(filenames["V3"]))
  Batch_update <- Map(cbind, files, Simulation = Sim)
  Batch_update <- Map(cbind, Batch_update, Batch = Batch)
  Batch_update <- Map(cbind, Batch_update, Land = Land)
  combined <- bind_rows(Batch_update, .id = "column_label")
  
  
  
  combined <- subset(combined, combined$Ninds == 0) # only patches receiving 0 immigrants
  combined <- subset(combined, combined$EndPatch != -999) #remove failed attempts
  number_isolated_patches <- combined %>%
    group_by(Simulation, Batch, Land, Rep, Year) %>%
    summarise(
      isolated_patch_count = n_distinct(EndPatch))
  
  
  # filter for only large patches ####
  library(dplyr)
  
  
  setwd(inputpath)
  patch <- raster(patch) #import patch map
  patch_frequency <- as.data.frame(freq(patch)) #calculate freq
  patch_frequency_NFI <- subset(patch_frequency, patch_frequency$count>200) #subset to Ha
  patch_frequency_NFI <- subset(patch_frequency_NFI, patch_frequency_NFI$value != 0) #remove matric cells
  
  NFIpatches <- combined %>%
    filter(EndPatch %in% patch_frequency_NFI$value) #filter isolated patches by NFI only
  
  number_isolated_NFIpatches<- NFIpatches %>%
    group_by(Simulation, Batch, Land, Rep, Year) %>%
    summarise(isolated_NFIpatch_count = n_distinct(EndPatch))
  
  mergeCols <- c("Simulation", "Batch", "Land", "Rep", "Year")
  isolated_patch <- merge(number_isolated_NFIpatches, number_isolated_patches, by = mergeCols, all = TRUE)
  isolated_patch$totalNFI <- length(patch_frequency_NFI$value)
  isolated_patch$proportion_isolated_NFI  <- isolated_patch$isolated_NFIpatch_count / isolated_patch$totalNFI
  isolated_patch$total_patch <- length(patch_frequency$value)
  isolated_patch$proportion_isolated  <- isolated_patch$isolated_patch_count / isolated_patch$total_patch
  
  isolated_patch <- subset(isolated_patch, select = - c(total_patch, totalNFI))
  
  return(isolated_patch)
}

# Batch1 ####
Batch1001_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch1001/Inputs", patch = "L1_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch1001/Outputs", mypattern = "Connect.txt", sep="\t")
Batch1002_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch1002/Inputs", patch = "L1_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch1002/Outputs", mypattern = "Connect.txt", sep="\t")
Batch1003_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch1003/Inputs", patch = "L1_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch1003/Outputs", mypattern = "Connect.txt", sep="\t")
Batch1004_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch1004/Inputs", patch = "L2_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch1004/Outputs", mypattern = "Connect.txt", sep="\t")
Batch1005_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch1005/Inputs", patch = "L2_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch1005/Outputs", mypattern = "Connect.txt", sep="\t")
Batch1006_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch1006/Inputs", patch = "L2_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch1006/Outputs", mypattern = "Connect.txt", sep="\t")
Batch1007_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch1007/Inputs", patch = "L3_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch1007/Outputs", mypattern = "Connect.txt", sep="\t")
Batch1008_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch1008/Inputs", patch = "L3_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch1008/Outputs", mypattern = "Connect.txt", sep="\t")
Batch1009_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch1009/Inputs", patch = "L3_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch1009/Outputs", mypattern = "Connect.txt", sep="\t")
Batch1 <- rbind(Batch1001_isolated_patch, Batch1002_isolated_patch, Batch1003_isolated_patch,Batch1004_isolated_patch, Batch1005_isolated_patch, Batch1006_isolated_patch, Batch1007_isolated_patch, Batch1008_isolated_patch, Batch1009_isolated_patch)
write.csv(Batch1, file = "E:/Models_Outputs/data_analysis/Isolated_patch/isolated_land1.csv", row.names = F)








# Batch2 ####
Batch2001_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch2001/Inputs", patch = "L1_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch2001/Outputs", mypattern = "Connect.txt", sep="\t")
Batch2002_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch2002/Inputs", patch = "L1_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch2002/Outputs", mypattern = "Connect.txt", sep="\t")
Batch2003_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch2003/Inputs", patch = "L1_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch2003/Outputs", mypattern = "Connect.txt", sep="\t")
Batch2004_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch2004/Inputs", patch = "L2_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch2004/Outputs", mypattern = "Connect.txt", sep="\t")
Batch2005_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch2005/Inputs", patch = "L2_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch2005/Outputs", mypattern = "Connect.txt", sep="\t")
Batch2006_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch2006/Inputs", patch = "L2_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch2006/Outputs", mypattern = "Connect.txt", sep="\t")
Batch2007_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch2007/Inputs", patch = "L3_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch2007/Outputs", mypattern = "Connect.txt", sep="\t")
Batch2008_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch2008/Inputs", patch = "L3_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch2008/Outputs", mypattern = "Connect.txt", sep="\t")
Batch2009_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch2009/Inputs", patch = "L3_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch2009/Outputs", mypattern = "Connect.txt", sep="\t")
Batch2 <- rbind(Batch2001_isolated_patch, Batch2002_isolated_patch, Batch2003_isolated_patch,Batch2004_isolated_patch, Batch2005_isolated_patch, Batch2006_isolated_patch, Batch2007_isolated_patch, Batch2008_isolated_patch, Batch2009_isolated_patch)
write.csv(Batch2, file = "E:/Models_Outputs/data_analysis/Isolated_patch/isolated_land2.csv", row.names = F)








# Batch3 ####
Batch3001_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch3001/Inputs", patch = "L1_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch3001/Outputs", mypattern = "Connect.txt", sep="\t")
Batch3002_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch3002/Inputs", patch = "L1_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch3002/Outputs", mypattern = "Connect.txt", sep="\t")
Batch3003_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch3003/Inputs", patch = "L1_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch3003/Outputs", mypattern = "Connect.txt", sep="\t")
Batch3004_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch3004/Inputs", patch = "L2_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch3004/Outputs", mypattern = "Connect.txt", sep="\t")
Batch3005_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch3005/Inputs", patch = "L2_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch3005/Outputs", mypattern = "Connect.txt", sep="\t")
Batch3006_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch3006/Inputs", patch = "L2_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch3006/Outputs", mypattern = "Connect.txt", sep="\t")
Batch3007_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch3007/Inputs", patch = "L3_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch3007/Outputs", mypattern = "Connect.txt", sep="\t")
Batch3008_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch3008/Inputs", patch = "L3_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch3008/Outputs", mypattern = "Connect.txt", sep="\t")
Batch3009_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch3009/Inputs", patch = "L3_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch3009/Outputs", mypattern = "Connect.txt", sep="\t")
Batch3 <- rbind(Batch3001_isolated_patch, Batch3002_isolated_patch, Batch3003_isolated_patch,Batch3004_isolated_patch, Batch3005_isolated_patch, Batch3006_isolated_patch, Batch3007_isolated_patch, Batch3008_isolated_patch, Batch3009_isolated_patch)
write.csv(Batch3, file = "E:/Models_Outputs/data_analysis/Isolated_patch/isolated_land3.csv", row.names = F)








# Batch4 ####
Batch4001_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch4001/Inputs", patch = "L1_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch4001/Outputs", mypattern = "Connect.txt", sep="\t")
Batch4002_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch4002/Inputs", patch = "L1_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch4002/Outputs", mypattern = "Connect.txt", sep="\t")
Batch4003_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch4003/Inputs", patch = "L1_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch4003/Outputs", mypattern = "Connect.txt", sep="\t")
Batch4004_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch4004/Inputs", patch = "L2_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch4004/Outputs", mypattern = "Connect.txt", sep="\t")
Batch4005_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch4005/Inputs", patch = "L2_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch4005/Outputs", mypattern = "Connect.txt", sep="\t")
Batch4006_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch4006/Inputs", patch = "L2_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch4006/Outputs", mypattern = "Connect.txt", sep="\t")
Batch4007_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch4007/Inputs", patch = "L3_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch4007/Outputs", mypattern = "Connect.txt", sep="\t")
Batch4008_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch4008/Inputs", patch = "L3_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch4008/Outputs", mypattern = "Connect.txt", sep="\t")
Batch4009_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch4009/Inputs", patch = "L3_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch4009/Outputs", mypattern = "Connect.txt", sep="\t")
Batch4 <- rbind(Batch4001_isolated_patch, Batch4002_isolated_patch, Batch4003_isolated_patch,Batch4004_isolated_patch, Batch4005_isolated_patch, Batch4006_isolated_patch, Batch4007_isolated_patch, Batch4008_isolated_patch, Batch4009_isolated_patch)
write.csv(Batch4, file = "E:/Models_Outputs/data_analysis/Isolated_patch/isolated_land4.csv", row.names = F)








# Batch5 ####
Batch5001_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch5001/Inputs", patch = "L1_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch5001/Outputs", mypattern = "Connect.txt", sep="\t")
Batch5002_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch5002/Inputs", patch = "L1_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch5002/Outputs", mypattern = "Connect.txt", sep="\t")
Batch5003_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch5003/Inputs", patch = "L1_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch5003/Outputs", mypattern = "Connect.txt", sep="\t")
Batch5004_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch5004/Inputs", patch = "L2_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch5004/Outputs", mypattern = "Connect.txt", sep="\t")
Batch5005_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch5005/Inputs", patch = "L2_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch5005/Outputs", mypattern = "Connect.txt", sep="\t")
Batch5006_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch5006/Inputs", patch = "L2_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch5006/Outputs", mypattern = "Connect.txt", sep="\t")
Batch5007_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch5007/Inputs", patch = "L3_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch5007/Outputs", mypattern = "Connect.txt", sep="\t")
Batch5008_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch5008/Inputs", patch = "L3_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch5008/Outputs", mypattern = "Connect.txt", sep="\t")
Batch5009_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch5009/Inputs", patch = "L3_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch5009/Outputs", mypattern = "Connect.txt", sep="\t")
Batch5 <- rbind(Batch5001_isolated_patch, Batch5002_isolated_patch, Batch5003_isolated_patch,Batch5004_isolated_patch, Batch5005_isolated_patch, Batch5006_isolated_patch, Batch5007_isolated_patch, Batch5008_isolated_patch, Batch5009_isolated_patch)
write.csv(Batch5, file = "E:/Models_Outputs/data_analysis/Isolated_patch/isolated_land5.csv", row.names = F)


# Batch6 ####
Batch6001_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch6001/Inputs", patch = "L1_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch6001/Outputs", mypattern = "Connect.txt", sep="\t")
Batch6002_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch6002/Inputs", patch = "L1_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch6002/Outputs", mypattern = "Connect.txt", sep="\t")
Batch6003_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch6003/Inputs", patch = "L1_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch6003/Outputs", mypattern = "Connect.txt", sep="\t")
Batch6004_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch6004/Inputs", patch = "L2_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch6004/Outputs", mypattern = "Connect.txt", sep="\t")
Batch6005_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch6005/Inputs", patch = "L2_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch6005/Outputs", mypattern = "Connect.txt", sep="\t")
Batch6006_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch6006/Inputs", patch = "L2_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch6006/Outputs", mypattern = "Connect.txt", sep="\t")
Batch6007_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch6007/Inputs", patch = "L3_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch6007/Outputs", mypattern = "Connect.txt", sep="\t")
Batch6008_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch6008/Inputs", patch = "L3_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch6008/Outputs", mypattern = "Connect.txt", sep="\t")
Batch6009_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch6009/Inputs", patch = "L3_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch6009/Outputs", mypattern = "Connect.txt", sep="\t")
Batch6 <- rbind(Batch6001_isolated_patch, Batch6002_isolated_patch, Batch6003_isolated_patch,Batch6004_isolated_patch, Batch6005_isolated_patch, Batch6006_isolated_patch, Batch6007_isolated_patch, Batch6008_isolated_patch, Batch6009_isolated_patch)
write.csv(Batch6, file = "E:/Models_Outputs/data_analysis/Isolated_patch/isolated_land6.csv", row.names = F)



# Batch7 ####
Batch7001_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch7001/Inputs", patch = "L1_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch7001/Outputs", mypattern = "Connect.txt", sep="\t")
Batch7002_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch7002/Inputs", patch = "L1_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch7002/Outputs", mypattern = "Connect.txt", sep="\t")
Batch7003_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch7003/Inputs", patch = "L1_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch7003/Outputs", mypattern = "Connect.txt", sep="\t")
Batch7004_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch7004/Inputs", patch = "L2_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch7004/Outputs", mypattern = "Connect.txt", sep="\t")
Batch7005_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch7005/Inputs", patch = "L2_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch7005/Outputs", mypattern = "Connect.txt", sep="\t")
Batch7006_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch7006/Inputs", patch = "L2_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch7006/Outputs", mypattern = "Connect.txt", sep="\t")
Batch7007_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch7007/Inputs", patch = "L3_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch7007/Outputs", mypattern = "Connect.txt", sep="\t")
Batch7008_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch7008/Inputs", patch = "L3_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch7008/Outputs", mypattern = "Connect.txt", sep="\t")
Batch7009_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch7009/Inputs", patch = "L3_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch7009/Outputs", mypattern = "Connect.txt", sep="\t")
Batch7 <- rbind(Batch7001_isolated_patch, Batch7002_isolated_patch, Batch7003_isolated_patch,Batch7004_isolated_patch, Batch7005_isolated_patch, Batch7006_isolated_patch, Batch7007_isolated_patch, Batch7008_isolated_patch, Batch7009_isolated_patch)
write.csv(Batch7, file = "E:/Models_Outputs/data_analysis/Isolated_patch/isolated_land7.csv", row.names = F)



# Batch8 ####
Batch8001_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch8001/Inputs", patch = "L1_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch8001/Outputs", mypattern = "Connect.txt", sep="\t")
Batch8002_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch8002/Inputs", patch = "L1_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch8002/Outputs", mypattern = "Connect.txt", sep="\t")
Batch8003_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch8003/Inputs", patch = "L1_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch8003/Outputs", mypattern = "Connect.txt", sep="\t")
Batch8004_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch8004/Inputs", patch = "L2_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch8004/Outputs", mypattern = "Connect.txt", sep="\t")
Batch8005_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch8005/Inputs", patch = "L2_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch8005/Outputs", mypattern = "Connect.txt", sep="\t")
Batch8006_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch8006/Inputs", patch = "L2_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch8006/Outputs", mypattern = "Connect.txt", sep="\t")
Batch8007_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch8007/Inputs", patch = "L3_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch8007/Outputs", mypattern = "Connect.txt", sep="\t")
Batch8008_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch8008/Inputs", patch = "L3_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch8008/Outputs", mypattern = "Connect.txt", sep="\t")
Batch8009_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch8009/Inputs", patch = "L3_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch8009/Outputs", mypattern = "Connect.txt", sep="\t")
Batch8 <- rbind(Batch8001_isolated_patch, Batch8002_isolated_patch, Batch8003_isolated_patch,Batch8004_isolated_patch, Batch8005_isolated_patch, Batch8006_isolated_patch, Batch8007_isolated_patch, Batch8008_isolated_patch, Batch8009_isolated_patch)
write.csv(Batch8, file = "E:/Models_Outputs/data_analysis/Isolated_patch/isolated_land8.csv", row.names = F)




# Batch9 ####
Batch9001_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch9001/Inputs", patch = "L1_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch9001/Outputs", mypattern = "Connect.txt", sep="\t")
Batch9002_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch9002/Inputs", patch = "L1_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch9002/Outputs", mypattern = "Connect.txt", sep="\t")
Batch9003_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch9003/Inputs", patch = "L1_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch9003/Outputs", mypattern = "Connect.txt", sep="\t")
Batch9004_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch9004/Inputs", patch = "L2_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch9004/Outputs", mypattern = "Connect.txt", sep="\t")
Batch9005_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch9005/Inputs", patch = "L2_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch9005/Outputs", mypattern = "Connect.txt", sep="\t")
Batch9006_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch9006/Inputs", patch = "L2_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch9006/Outputs", mypattern = "Connect.txt", sep="\t")
Batch9007_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch9007/Inputs", patch = "L3_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch9007/Outputs", mypattern = "Connect.txt", sep="\t")
Batch9008_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch9008/Inputs", patch = "L3_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch9008/Outputs", mypattern = "Connect.txt", sep="\t")
Batch9009_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch9009/Inputs", patch = "L3_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch9009/Outputs", mypattern = "Connect.txt", sep="\t")
Batch9 <- rbind(Batch9001_isolated_patch, Batch9002_isolated_patch, Batch9003_isolated_patch,Batch9004_isolated_patch, Batch9005_isolated_patch, Batch9006_isolated_patch, Batch9007_isolated_patch, Batch9008_isolated_patch, Batch9009_isolated_patch)
write.csv(Batch9, file = "E:/Models_Outputs/data_analysis/Isolated_patch/isolated_land9.csv", row.names = F)



# Batch10 ####
Batch10001_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch10001/Inputs", patch = "L1_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch10001/Outputs", mypattern = "Connect.txt", sep="\t")
Batch10002_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch10002/Inputs", patch = "L1_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch10002/Outputs", mypattern = "Connect.txt", sep="\t")
Batch10003_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch10003/Inputs", patch = "L1_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch10003/Outputs", mypattern = "Connect.txt", sep="\t")
Batch10004_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch10004/Inputs", patch = "L2_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch10004/Outputs", mypattern = "Connect.txt", sep="\t")
Batch10005_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch10005/Inputs", patch = "L2_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch10005/Outputs", mypattern = "Connect.txt", sep="\t")
Batch10006_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch10006/Inputs", patch = "L2_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch10006/Outputs", mypattern = "Connect.txt", sep="\t")
Batch10007_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch10007/Inputs", patch = "L3_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch10007/Outputs", mypattern = "Connect.txt", sep="\t")
Batch10008_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch10008/Inputs", patch = "L3_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch10008/Outputs", mypattern = "Connect.txt", sep="\t")
Batch10009_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch10009/Inputs", patch = "L3_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch10009/Outputs", mypattern = "Connect.txt", sep="\t")
Batch10 <- rbind(Batch10001_isolated_patch, Batch10002_isolated_patch, Batch10003_isolated_patch,Batch10004_isolated_patch, Batch10005_isolated_patch, Batch10006_isolated_patch, Batch10007_isolated_patch, Batch10008_isolated_patch, Batch10009_isolated_patch)
write.csv(Batch10, file = "E:/Models_Outputs/data_analysis/Isolated_patch/isolated_land10.csv", row.names = F)

# Batch11 ####
Batch11001_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch11001/Inputs", patch = "L1_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch11001/Outputs", mypattern = "Connect.txt", sep="\t")
Batch11002_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch11002/Inputs", patch = "L1_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch11002/Outputs", mypattern = "Connect.txt", sep="\t")
Batch11003_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch11003/Inputs", patch = "L1_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch11003/Outputs", mypattern = "Connect.txt", sep="\t")
Batch11004_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch11004/Inputs", patch = "L2_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch11004/Outputs", mypattern = "Connect.txt", sep="\t")
Batch11005_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch11005/Inputs", patch = "L2_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch11005/Outputs", mypattern = "Connect.txt", sep="\t")
Batch11006_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch11006/Inputs", patch = "L2_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch11006/Outputs", mypattern = "Connect.txt", sep="\t")
Batch11007_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch11007/Inputs", patch = "L3_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch11007/Outputs", mypattern = "Connect.txt", sep="\t")
Batch11008_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch11008/Inputs", patch = "L3_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch11008/Outputs", mypattern = "Connect.txt", sep="\t")
Batch11009_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch11009/Inputs", patch = "L3_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch11009/Outputs", mypattern = "Connect.txt", sep="\t")
Batch11 <- rbind(Batch11001_isolated_patch, Batch11002_isolated_patch, Batch11003_isolated_patch,Batch11004_isolated_patch, Batch11005_isolated_patch, Batch11006_isolated_patch, Batch11007_isolated_patch, Batch11008_isolated_patch, Batch11009_isolated_patch)
write.csv(Batch11, file = "E:/Models_Outputs/data_analysis/Isolated_patch/isolated_land11.csv", row.names = F)


# Batch12 ####
Batch12001_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch12001/Inputs", patch = "L1_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch12001/Outputs", mypattern = "Connect.txt", sep="\t")
Batch12002_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch12002/Inputs", patch = "L1_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch12002/Outputs", mypattern = "Connect.txt", sep="\t")
Batch12003_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch12003/Inputs", patch = "L1_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch12003/Outputs", mypattern = "Connect.txt", sep="\t")
Batch12004_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch12004/Inputs", patch = "L2_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch12004/Outputs", mypattern = "Connect.txt", sep="\t")
Batch12005_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch12005/Inputs", patch = "L2_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch12005/Outputs", mypattern = "Connect.txt", sep="\t")
Batch12006_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch12006/Inputs", patch = "L2_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch12006/Outputs", mypattern = "Connect.txt", sep="\t")
Batch12007_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch12007/Inputs", patch = "L3_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch12007/Outputs", mypattern = "Connect.txt", sep="\t")
Batch12008_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch12008/Inputs", patch = "L3_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch12008/Outputs", mypattern = "Connect.txt", sep="\t")
Batch12009_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch12009/Inputs", patch = "L3_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch12009/Outputs", mypattern = "Connect.txt", sep="\t")
Batch12 <- rbind(Batch12001_isolated_patch, Batch12002_isolated_patch, Batch12003_isolated_patch,Batch12004_isolated_patch, Batch12005_isolated_patch, Batch12006_isolated_patch, Batch12007_isolated_patch, Batch12008_isolated_patch, Batch12009_isolated_patch)
write.csv(Batch12, file = "E:/Models_Outputs/data_analysis/Isolated_patch/isolated_land12.csv", row.names = F)

































# Batch13 ####
Batch13001_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch13001/Inputs", patch = "L1_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch13001/Outputs", mypattern = "Connect.txt", sep="\t")
Batch13002_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch13002/Inputs", patch = "L1_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch13002/Outputs", mypattern = "Connect.txt", sep="\t")
Batch13003_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch13003/Inputs", patch = "L1_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch13003/Outputs", mypattern = "Connect.txt", sep="\t")
Batch13004_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch13004/Inputs", patch = "L2_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch13004/Outputs", mypattern = "Connect.txt", sep="\t")
Batch13005_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch13005/Inputs", patch = "L2_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch13005/Outputs", mypattern = "Connect.txt", sep="\t")
Batch13006_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch13006/Inputs", patch = "L2_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch13006/Outputs", mypattern = "Connect.txt", sep="\t")
Batch13007_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch13007/Inputs", patch = "L3_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch13007/Outputs", mypattern = "Connect.txt", sep="\t")
Batch13008_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch13008/Inputs", patch = "L3_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch13008/Outputs", mypattern = "Connect.txt", sep="\t")
Batch13009_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch13009/Inputs", patch = "L3_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch13009/Outputs", mypattern = "Connect.txt", sep="\t")
Batch13 <- rbind(Batch13001_isolated_patch, Batch13002_isolated_patch, Batch13003_isolated_patch,Batch13004_isolated_patch, Batch13005_isolated_patch, Batch13006_isolated_patch, Batch13007_isolated_patch, Batch13008_isolated_patch, Batch13009_isolated_patch)
write.csv(Batch13, file = "E:/Models_Outputs/data_analysis/Isolated_patch/isolated_land13.csv", row.names = F)
# Batch14 ####
Batch14001_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch14001/Inputs", patch = "L1_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch14001/Outputs", mypattern = "Connect.txt", sep="\t")
Batch14002_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch14002/Inputs", patch = "L1_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch14002/Outputs", mypattern = "Connect.txt", sep="\t")
Batch14003_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch14003/Inputs", patch = "L1_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch14003/Outputs", mypattern = "Connect.txt", sep="\t")
Batch14004_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch14004/Inputs", patch = "L2_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch14004/Outputs", mypattern = "Connect.txt", sep="\t")
Batch14005_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch14005/Inputs", patch = "L2_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch14005/Outputs", mypattern = "Connect.txt", sep="\t")
Batch14006_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch14006/Inputs", patch = "L2_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch14006/Outputs", mypattern = "Connect.txt", sep="\t")
Batch14007_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch14007/Inputs", patch = "L3_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch14007/Outputs", mypattern = "Connect.txt", sep="\t")
Batch14008_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch14008/Inputs", patch = "L3_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch14008/Outputs", mypattern = "Connect.txt", sep="\t")
Batch14009_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch14009/Inputs", patch = "L3_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch14009/Outputs", mypattern = "Connect.txt", sep="\t")
Batch14 <- rbind(Batch14001_isolated_patch, Batch14002_isolated_patch, Batch14003_isolated_patch,Batch14004_isolated_patch, Batch14005_isolated_patch, Batch14006_isolated_patch, Batch14007_isolated_patch, Batch14008_isolated_patch, Batch14009_isolated_patch)
write.csv(Batch14, file = "E:/Models_Outputs/data_analysis/Isolated_patch/isolated_land14.csv", row.names = F)
# Batch15 ####
Batch15001_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch15001/Inputs", patch = "L1_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch15001/Outputs", mypattern = "Connect.txt", sep="\t")
Batch15002_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch15002/Inputs", patch = "L1_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch15002/Outputs", mypattern = "Connect.txt", sep="\t")
Batch15003_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch15003/Inputs", patch = "L1_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch15003/Outputs", mypattern = "Connect.txt", sep="\t")
Batch15004_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch15004/Inputs", patch = "L2_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch15004/Outputs", mypattern = "Connect.txt", sep="\t")
Batch15005_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch15005/Inputs", patch = "L2_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch15005/Outputs", mypattern = "Connect.txt", sep="\t")
Batch15006_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch15006/Inputs", patch = "L2_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch15006/Outputs", mypattern = "Connect.txt", sep="\t")
Batch15007_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch15007/Inputs", patch = "L3_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch15007/Outputs", mypattern = "Connect.txt", sep="\t")
Batch15008_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch15008/Inputs", patch = "L3_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch15008/Outputs", mypattern = "Connect.txt", sep="\t")
Batch15009_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch15009/Inputs", patch = "L3_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch15009/Outputs", mypattern = "Connect.txt", sep="\t")
Batch15 <- rbind(Batch15001_isolated_patch, Batch15002_isolated_patch, Batch15003_isolated_patch,Batch15004_isolated_patch, Batch15005_isolated_patch, Batch15006_isolated_patch, Batch15007_isolated_patch, Batch15008_isolated_patch, Batch15009_isolated_patch)
write.csv(Batch15, file = "E:/Models_Outputs/data_analysis/Isolated_patch/isolated_land15.csv", row.names = F)
# Batch16 ####
Batch16001_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch16001/Inputs", patch = "L1_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch16001/Outputs", mypattern = "Connect.txt", sep="\t")
Batch16002_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch16002/Inputs", patch = "L1_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch16002/Outputs", mypattern = "Connect.txt", sep="\t")
Batch16003_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch16003/Inputs", patch = "L1_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch16003/Outputs", mypattern = "Connect.txt", sep="\t")
Batch16004_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch16004/Inputs", patch = "L2_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch16004/Outputs", mypattern = "Connect.txt", sep="\t")
Batch16005_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch16005/Inputs", patch = "L2_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch16005/Outputs", mypattern = "Connect.txt", sep="\t")
Batch16006_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch16006/Inputs", patch = "L2_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch16006/Outputs", mypattern = "Connect.txt", sep="\t")
Batch16007_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch16007/Inputs", patch = "L3_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch16007/Outputs", mypattern = "Connect.txt", sep="\t")
Batch16008_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch16008/Inputs", patch = "L3_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch16008/Outputs", mypattern = "Connect.txt", sep="\t")
Batch16009_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch16009/Inputs", patch = "L3_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch16009/Outputs", mypattern = "Connect.txt", sep="\t")
Batch16 <- rbind(Batch16001_isolated_patch, Batch16002_isolated_patch, Batch16003_isolated_patch,Batch16004_isolated_patch, Batch16005_isolated_patch, Batch16006_isolated_patch, Batch16007_isolated_patch, Batch16008_isolated_patch, Batch16009_isolated_patch)
write.csv(Batch16, file = "E:/Models_Outputs/data_analysis/Isolated_patch/isolated_land16.csv", row.names = F)
# Batch17 ####
Batch17001_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch17001/Inputs", patch = "L1_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch17001/Outputs", mypattern = "Connect.txt", sep="\t")
Batch17002_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch17002/Inputs", patch = "L1_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch17002/Outputs", mypattern = "Connect.txt", sep="\t")
Batch17003_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch17003/Inputs", patch = "L1_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch17003/Outputs", mypattern = "Connect.txt", sep="\t")
Batch17004_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch17004/Inputs", patch = "L2_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch17004/Outputs", mypattern = "Connect.txt", sep="\t")
Batch17005_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch17005/Inputs", patch = "L2_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch17005/Outputs", mypattern = "Connect.txt", sep="\t")
Batch17006_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch17006/Inputs", patch = "L2_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch17006/Outputs", mypattern = "Connect.txt", sep="\t")
Batch17007_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch17007/Inputs", patch = "L3_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch17007/Outputs", mypattern = "Connect.txt", sep="\t")
Batch17008_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch17008/Inputs", patch = "L3_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch17008/Outputs", mypattern = "Connect.txt", sep="\t")
Batch17009_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch17009/Inputs", patch = "L3_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch17009/Outputs", mypattern = "Connect.txt", sep="\t")
Batch17 <- rbind(Batch17001_isolated_patch, Batch17002_isolated_patch, Batch17003_isolated_patch,Batch17004_isolated_patch, Batch17005_isolated_patch, Batch17006_isolated_patch, Batch17007_isolated_patch, Batch17008_isolated_patch, Batch17009_isolated_patch)
write.csv(Batch17, file = "E:/Models_Outputs/data_analysis/Isolated_patch/isolated_land17.csv", row.names = F)


# Batch18 ####
Batch18001_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch18001/Inputs", patch = "L1_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch18001/Outputs", mypattern = "Connect.txt", sep="\t")
Batch18002_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch18002/Inputs", patch = "L1_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch18002/Outputs", mypattern = "Connect.txt", sep="\t")
Batch18003_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch18003/Inputs", patch = "L1_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch18003/Outputs", mypattern = "Connect.txt", sep="\t")
Batch18004_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch18004/Inputs", patch = "L2_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch18004/Outputs", mypattern = "Connect.txt", sep="\t")
Batch18005_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch18005/Inputs", patch = "L2_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch18005/Outputs", mypattern = "Connect.txt", sep="\t")
Batch18006_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch18006/Inputs", patch = "L2_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch18006/Outputs", mypattern = "Connect.txt", sep="\t")
Batch18007_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch18007/Inputs", patch = "L3_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch18007/Outputs", mypattern = "Connect.txt", sep="\t")
Batch18008_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch18008/Inputs", patch = "L3_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch18008/Outputs", mypattern = "Connect.txt", sep="\t")
Batch18009_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch18009/Inputs", patch = "L3_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch18009/Outputs", mypattern = "Connect.txt", sep="\t")
Batch18 <- rbind(Batch18001_isolated_patch, Batch18002_isolated_patch, Batch18003_isolated_patch,Batch18004_isolated_patch, Batch18005_isolated_patch, Batch18006_isolated_patch, Batch18007_isolated_patch, Batch18008_isolated_patch, Batch18009_isolated_patch)
write.csv(Batch18, file = "E:/Models_Outputs/data_analysis/Isolated_patch/isolated_land18.csv", row.names = F)

# Batch19 ####
Batch19001_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch19001/Inputs", patch = "L1_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch19001/Outputs", mypattern = "Connect.txt", sep="\t")
Batch19002_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch19002/Inputs", patch = "L1_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch19002/Outputs", mypattern = "Connect.txt", sep="\t")
Batch19003_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch19003/Inputs", patch = "L1_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch19003/Outputs", mypattern = "Connect.txt", sep="\t")
Batch19004_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch19004/Inputs", patch = "L2_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch19004/Outputs", mypattern = "Connect.txt", sep="\t")
Batch19005_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch19005/Inputs", patch = "L2_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch19005/Outputs", mypattern = "Connect.txt", sep="\t")
Batch19006_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch19006/Inputs", patch = "L2_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch19006/Outputs", mypattern = "Connect.txt", sep="\t")
Batch19007_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch19007/Inputs", patch = "L3_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch19007/Outputs", mypattern = "Connect.txt", sep="\t")
Batch19008_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch19008/Inputs", patch = "L3_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch19008/Outputs", mypattern = "Connect.txt", sep="\t")
Batch19009_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch19009/Inputs", patch = "L3_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch19009/Outputs", mypattern = "Connect.txt", sep="\t")
Batch19 <- rbind(Batch19001_isolated_patch, Batch19002_isolated_patch, Batch19003_isolated_patch,Batch19004_isolated_patch, Batch19005_isolated_patch, Batch19006_isolated_patch, Batch19007_isolated_patch, Batch19008_isolated_patch, Batch19009_isolated_patch)
write.csv(Batch19, file = "E:/Models_Outputs/data_analysis/Isolated_patch/isolated_land19.csv", row.names = F)



# Batch20 ####
Batch20001_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch20001/Inputs", patch = "L1_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch20001/Outputs", mypattern = "Connect.txt", sep="\t")
Batch20002_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch20002/Inputs", patch = "L1_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch20002/Outputs", mypattern = "Connect.txt", sep="\t")
Batch20003_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch20003/Inputs", patch = "L1_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch20003/Outputs", mypattern = "Connect.txt", sep="\t")
Batch20004_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch20004/Inputs", patch = "L2_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch20004/Outputs", mypattern = "Connect.txt", sep="\t")
Batch20005_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch20005/Inputs", patch = "L2_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch20005/Outputs", mypattern = "Connect.txt", sep="\t")
Batch20006_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch20006/Inputs", patch = "L2_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch20006/Outputs", mypattern = "Connect.txt", sep="\t")
Batch20007_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch20007/Inputs", patch = "L3_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch20007/Outputs", mypattern = "Connect.txt", sep="\t")
Batch20008_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch20008/Inputs", patch = "L3_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch20008/Outputs", mypattern = "Connect.txt", sep="\t")
Batch20009_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch20009/Inputs", patch = "L3_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch20009/Outputs", mypattern = "Connect.txt", sep="\t")
Batch20 <- rbind(Batch20001_isolated_patch, Batch20002_isolated_patch, Batch20003_isolated_patch,Batch20004_isolated_patch, Batch20005_isolated_patch, Batch20006_isolated_patch, Batch20007_isolated_patch, Batch20008_isolated_patch, Batch20009_isolated_patch)
write.csv(Batch20, file = "E:/Models_Outputs/data_analysis/Isolated_patch/isolated_land20.csv", row.names = F)
# Batch21 ####
Batch21001_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch21001/Inputs", patch = "L1_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch21001/Outputs", mypattern = "Connect.txt", sep="\t")
Batch21002_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch21002/Inputs", patch = "L1_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch21002/Outputs", mypattern = "Connect.txt", sep="\t")
Batch21003_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch21003/Inputs", patch = "L1_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch21003/Outputs", mypattern = "Connect.txt", sep="\t")
Batch21004_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch21004/Inputs", patch = "L2_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch21004/Outputs", mypattern = "Connect.txt", sep="\t")
Batch21005_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch21005/Inputs", patch = "L2_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch21005/Outputs", mypattern = "Connect.txt", sep="\t")
Batch21006_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch21006/Inputs", patch = "L2_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch21006/Outputs", mypattern = "Connect.txt", sep="\t")
Batch21007_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch21007/Inputs", patch = "L3_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch21007/Outputs", mypattern = "Connect.txt", sep="\t")
Batch21008_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch21008/Inputs", patch = "L3_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch21008/Outputs", mypattern = "Connect.txt", sep="\t")
Batch21009_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch21009/Inputs", patch = "L3_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch21009/Outputs", mypattern = "Connect.txt", sep="\t")
Batch21 <- rbind(Batch21001_isolated_patch, Batch21002_isolated_patch, Batch21003_isolated_patch,Batch21004_isolated_patch, Batch21005_isolated_patch, Batch21006_isolated_patch, Batch21007_isolated_patch, Batch21008_isolated_patch, Batch21009_isolated_patch)
write.csv(Batch21, file = "E:/Models_Outputs/data_analysis/Isolated_patch/isolated_land21.csv", row.names = F)
# Batch22 ####
Batch22001_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch22001/Inputs", patch = "L1_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch22001/Outputs", mypattern = "Connect.txt", sep="\t")
Batch22002_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch22002/Inputs", patch = "L1_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch22002/Outputs", mypattern = "Connect.txt", sep="\t")
Batch22003_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch22003/Inputs", patch = "L1_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch22003/Outputs", mypattern = "Connect.txt", sep="\t")
Batch22004_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch22004/Inputs", patch = "L2_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch22004/Outputs", mypattern = "Connect.txt", sep="\t")
Batch22005_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch22005/Inputs", patch = "L2_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch22005/Outputs", mypattern = "Connect.txt", sep="\t")
Batch22006_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch22006/Inputs", patch = "L2_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch22006/Outputs", mypattern = "Connect.txt", sep="\t")
Batch22007_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch22007/Inputs", patch = "L3_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch22007/Outputs", mypattern = "Connect.txt", sep="\t")
Batch22008_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch22008/Inputs", patch = "L3_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch22008/Outputs", mypattern = "Connect.txt", sep="\t")
Batch22009_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch22009/Inputs", patch = "L3_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch22009/Outputs", mypattern = "Connect.txt", sep="\t")
Batch22 <- rbind(Batch22001_isolated_patch, Batch22002_isolated_patch, Batch22003_isolated_patch,Batch22004_isolated_patch, Batch22005_isolated_patch, Batch22006_isolated_patch, Batch22007_isolated_patch, Batch22008_isolated_patch, Batch22009_isolated_patch)
write.csv(Batch22, file = "E:/Models_Outputs/data_analysis/Isolated_patch/isolated_land22.csv", row.names = F)
# Batch23 ####
Batch23001_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch23001/Inputs", patch = "L1_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch23001/Outputs", mypattern = "Connect.txt", sep="\t")
Batch23002_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch23002/Inputs", patch = "L1_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch23002/Outputs", mypattern = "Connect.txt", sep="\t")
Batch23003_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch23003/Inputs", patch = "L1_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch23003/Outputs", mypattern = "Connect.txt", sep="\t")
Batch23004_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch23004/Inputs", patch = "L2_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch23004/Outputs", mypattern = "Connect.txt", sep="\t")
Batch23005_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch23005/Inputs", patch = "L2_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch23005/Outputs", mypattern = "Connect.txt", sep="\t")
Batch23006_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch23006/Inputs", patch = "L2_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch23006/Outputs", mypattern = "Connect.txt", sep="\t")
Batch23007_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch23007/Inputs", patch = "L3_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch23007/Outputs", mypattern = "Connect.txt", sep="\t")
Batch23008_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch23008/Inputs", patch = "L3_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch23008/Outputs", mypattern = "Connect.txt", sep="\t")
Batch23009_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch23009/Inputs", patch = "L3_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch23009/Outputs", mypattern = "Connect.txt", sep="\t")
Batch23 <- rbind(Batch23001_isolated_patch, Batch23002_isolated_patch, Batch23003_isolated_patch,Batch23004_isolated_patch, Batch23005_isolated_patch, Batch23006_isolated_patch, Batch23007_isolated_patch, Batch23008_isolated_patch, Batch23009_isolated_patch)
write.csv(Batch23, file = "E:/Models_Outputs/data_analysis/Isolated_patch/isolated_land23.csv", row.names = F)
# Batch24 ####
Batch24001_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch24001/Inputs", patch = "L1_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch24001/Outputs", mypattern = "Connect.txt", sep="\t")
Batch24002_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch24002/Inputs", patch = "L1_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch24002/Outputs", mypattern = "Connect.txt", sep="\t")
Batch24003_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch24003/Inputs", patch = "L1_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch24003/Outputs", mypattern = "Connect.txt", sep="\t")
Batch24004_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch24004/Inputs", patch = "L2_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch24004/Outputs", mypattern = "Connect.txt", sep="\t")
Batch24005_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch24005/Inputs", patch = "L2_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch24005/Outputs", mypattern = "Connect.txt", sep="\t")
Batch24006_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch24006/Inputs", patch = "L2_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch24006/Outputs", mypattern = "Connect.txt", sep="\t")
Batch24007_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch24007/Inputs", patch = "L3_R1_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch24007/Outputs", mypattern = "Connect.txt", sep="\t")
Batch24008_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch24008/Inputs", patch = "L3_R2_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch24008/Outputs", mypattern = "Connect.txt", sep="\t")
Batch24009_isolated_patch <-  Isolated_patch(inputpath = "E:/Models_Outputs/Batch24009/Inputs", patch = "L3_R3_Patch_Scenario0.asc",mypath = "E:/Models_Outputs/Batch24009/Outputs", mypattern = "Connect.txt", sep="\t")
Batch24 <- rbind(Batch24001_isolated_patch, Batch24002_isolated_patch, Batch24003_isolated_patch,Batch24004_isolated_patch, Batch24005_isolated_patch, Batch24006_isolated_patch, Batch24007_isolated_patch, Batch24008_isolated_patch, Batch24009_isolated_patch)
write.csv(Batch24, file = "E:/Models_Outputs/data_analysis/Isolated_patch/isolated_land24.csv", row.names = F)
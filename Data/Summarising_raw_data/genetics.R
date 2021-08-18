

# this function imports population files and genetic files, assigns pop size to all patches and filters genetic data by min pop size

Filter_genetics <- function(min_pop_size, mypath, ...) {
  setwd(mypath)
  library(dplyr)
  import.multiple.files<-function(mypath,pattern,...)
  {
    setwd(mypath)
    tmp.list.1<-list.files(mypath, pattern)
    tmp.list.2<-list(length=length(tmp.list.1))
    for (i in 1:length(tmp.list.1)){tmp.list.2[[i]]<-read.csv(tmp.list.1[i],...)}
    names(tmp.list.2)<-tmp.list.1
    tmp.list.2
  }
  pop_files <-
    import.multiple.files(mypath = mypath,
                          pattern = "Pop.txt",
                          sep = "\t")
  pop_filenames <-
    list.files(mypath, pattern = "Pop.txt")
  pop_filenames <- as.data.frame(pop_filenames)
  pop_filenames <-
    as.data.frame(str_split_fixed(pop_filenames$pop_filenames, "_", 4))
  pop_Batch <- unname(unlist(pop_filenames["V1"]))
  pop_Sim <- unname(unlist(pop_filenames["V2"]))
  pop_Land <- unname(unlist(pop_filenames["V3"]))
  Pop <- Map(cbind, pop_files, Simulation = pop_Sim)
  Pop <- Map(cbind, Pop, Batch = pop_Batch)
  Pop <- Map(cbind, Pop, Land = pop_Land)
  Population <- bind_rows(Pop, .id = "column_label")
  
  gen_files <-
    import.multiple.files(mypath = mypath,
                          pattern = "LandGen.txt",
                          sep = "\t")
  gen_filenames <-
    list.files(mypath, pattern = "LandGen.txt")
  gen_filenames <- as.data.frame(gen_filenames)
  gen_filenames <-
    as.data.frame(str_split_fixed(gen_filenames$gen_filenames, "_", 4))
  gen_Batch <- unname(unlist(gen_filenames["V1"]))
  gen_Sim <- unname(unlist(gen_filenames["V2"]))
  gen_Land <- unname(unlist(gen_filenames["V3"]))
  gen <- Map(cbind, gen_files, Simulation = gen_Sim)
  gen <- Map(cbind, gen, Batch = gen_Batch)
  gen <- Map(cbind, gen, Land = gen_Land)
  Genetics <- bind_rows(gen, .id = "column_label")
  
  Genetics$PatchID0_size <- Population[match(paste(Genetics$Rep, Genetics$Year, Genetics$Simulation, Genetics$Batch, Genetics$Land, Genetics$PatchID0),paste(Population$Rep, Population$Year, Population$Simulation, Population$Batch, Population$Land, Population$PatchID)),"NInd"]
  Genetics$PatchID1_size <- Population[match(paste(Genetics$Rep, Genetics$Year, Genetics$Simulation, Genetics$Batch, Genetics$Land, Genetics$PatchID1),paste(Population$Rep, Population$Year, Population$Simulation, Population$Batch, Population$Land, Population$PatchID)),"NInd"]
  
  
  Genetics <- subset(Genetics, Genetics$PatchID0_size > min_pop_size)
  Genetics <- subset(Genetics, Genetics$PatchID1_size > min_pop_size)
  
  Genetics$FST[Genetics$FST[] == -999] <- "NA"
  Genetics$FST[Genetics$FST < 0] <- 0
  Genetics$FST <- as.numeric(Genetics$FST)
  
  Genetics_summary <- Genetics %>%
    group_by(Simulation, Batch, Land, Rep, Year) %>%
    summarise(
      Distance = mean(Distance),
      MeanNAlleles = mean(MeanNAlleles),
      FIS = mean(FIS),
      FST = mean(FST),
      MeanD = mean(MeanD),
      HarmMeanD = mean(HarmMeanD)
    )
  
  
  return(Genetics_summary)
}



Batch1001_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch1001/Outputs")
Batch1002_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch1002/Outputs")
Batch1003_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch1003/Outputs")
Batch1004_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch1004/Outputs")
Batch1005_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch1005/Outputs")
Batch1006_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch1006/Outputs")
Batch1007_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch1007/Outputs")
Batch1008_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch1008/Outputs")
Batch1009_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch1009/Outputs")
Batch1_Genetics <-
  rbind(
    Batch1001_gen_Genetics,
    Batch1002_gen_Genetics,
    Batch1003_gen_Genetics,
    Batch1004_gen_Genetics,
    Batch1005_gen_Genetics,
    Batch1006_gen_Genetics,
    Batch1007_gen_Genetics,
    Batch1008_gen_Genetics,
    Batch1009_gen_Genetics
  )
write.csv(Batch1_Genetics, file = "E:/Models_outputs/data_analysis/Genetics/Batch1_Genetics.csv", row.names = F)





Batch2001_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch2001/Outputs")
Batch2002_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch2002/Outputs")
Batch2003_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch2003/Outputs")
Batch2004_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch2004/Outputs")
Batch2005_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch2005/Outputs")
Batch2006_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch2006/Outputs")
Batch2007_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch2007/Outputs")
Batch2008_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch2008/Outputs")
Batch2009_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch2009/Outputs")
Batch2_Genetics <-
  rbind(
    Batch2001_gen_Genetics,
    Batch2002_gen_Genetics,
    Batch2003_gen_Genetics,
    Batch2004_gen_Genetics,
    Batch2005_gen_Genetics,
    Batch2006_gen_Genetics,
    Batch2007_gen_Genetics,
    Batch2008_gen_Genetics,
    Batch2009_gen_Genetics
  )
write.csv(Batch2_Genetics, file = "E:/Models_outputs/data_analysis/Genetics/Batch2_Genetics.csv", row.names = F)



Batch3001_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch3001/Outputs")
Batch3002_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch3002/Outputs")
Batch3003_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch3003/Outputs")
Batch3004_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch3004/Outputs")
Batch3005_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch3005/Outputs")
Batch3006_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch3006/Outputs")
Batch3007_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch3007/Outputs")
Batch3008_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch3008/Outputs")
Batch3009_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch3009/Outputs")
Batch3_Genetics <-
  rbind(
    Batch3001_gen_Genetics,
    Batch3002_gen_Genetics,
    Batch3003_gen_Genetics,
    Batch3004_gen_Genetics,
    Batch3005_gen_Genetics,
    Batch3006_gen_Genetics,
    Batch3007_gen_Genetics,
    Batch3008_gen_Genetics,
    Batch3009_gen_Genetics
  )
write.csv(Batch3_Genetics, file = "E:/Models_outputs/data_analysis/Genetics/Batch3_Genetics.csv", row.names = F)


Batch4001_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch4001/Outputs")
Batch4002_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch4002/Outputs")
Batch4003_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch4003/Outputs")
Batch4004_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch4004/Outputs")
Batch4005_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch4005/Outputs")
Batch4006_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch4006/Outputs")
Batch4007_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch4007/Outputs")
Batch4008_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch4008/Outputs")
Batch4009_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch4009/Outputs")
Batch4_Genetics <-
  rbind(
    Batch4001_gen_Genetics,
    Batch4002_gen_Genetics,
    Batch4003_gen_Genetics,
    Batch4004_gen_Genetics,
    Batch4005_gen_Genetics,
    Batch4006_gen_Genetics,
    Batch4007_gen_Genetics,
    Batch4008_gen_Genetics,
    Batch4009_gen_Genetics
  )
write.csv(Batch4_Genetics, file = "E:/Models_outputs/data_analysis/Genetics/Batch4_Genetics.csv", row.names = F)


Batch5001_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch5001/Outputs")
Batch5002_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch5002/Outputs")
Batch5003_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch5003/Outputs")
Batch5004_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch5004/Outputs")
Batch5005_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch5005/Outputs")
Batch5006_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch5006/Outputs")
Batch5007_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch5007/Outputs")
Batch5008_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch5008/Outputs")
Batch5009_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch5009/Outputs")
Batch5_Genetics <-
  rbind(
    Batch5001_gen_Genetics,
    Batch5002_gen_Genetics,
    Batch5003_gen_Genetics,
    Batch5004_gen_Genetics,
    Batch5005_gen_Genetics,
    Batch5006_gen_Genetics,
    Batch5007_gen_Genetics,
    Batch5008_gen_Genetics,
    Batch5009_gen_Genetics
  )
write.csv(Batch5_Genetics, file = "E:/Models_outputs/data_analysis/Genetics/Batch5_Genetics.csv", row.names = F)


Batch6001_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch6001/Outputs")
Batch6002_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch6002/Outputs")
Batch6003_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch6003/Outputs")
Batch6004_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch6004/Outputs")
Batch6005_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch6005/Outputs")
Batch6006_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch6006/Outputs")
Batch6007_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch6007/Outputs")
Batch6008_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch6008/Outputs")
Batch6009_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch6009/Outputs")
Batch6_Genetics <-
  rbind(
    Batch6001_gen_Genetics,
    Batch6002_gen_Genetics,
    Batch6003_gen_Genetics,
    Batch6004_gen_Genetics,
    Batch6005_gen_Genetics,
    Batch6006_gen_Genetics,
    Batch6007_gen_Genetics,
    Batch6008_gen_Genetics,
    Batch6009_gen_Genetics
  )
write.csv(Batch6_Genetics, file = "E:/Models_outputs/data_analysis/Genetics/Batch6_Genetics.csv", row.names = F)


Batch7001_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch7001/Outputs")
Batch7002_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch7002/Outputs")
Batch7003_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch7003/Outputs")
Batch7004_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch7004/Outputs")
Batch7005_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch7005/Outputs")
Batch7006_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch7006/Outputs")
Batch7007_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch7007/Outputs")
Batch7008_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch7008/Outputs")
Batch7009_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch7009/Outputs")
Batch7_Genetics <-
  rbind(
    Batch7001_gen_Genetics,
    Batch7002_gen_Genetics,
    Batch7003_gen_Genetics,
    Batch7004_gen_Genetics,
    Batch7005_gen_Genetics,
    Batch7006_gen_Genetics,
    Batch7007_gen_Genetics,
    Batch7008_gen_Genetics,
    Batch7009_gen_Genetics
  )
write.csv(Batch7_Genetics, file = "E:/Models_outputs/data_analysis/Genetics/Batch7_Genetics.csv", row.names = F)


Batch8001_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch8001/Outputs")
Batch8002_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch8002/Outputs")
Batch8003_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch8003/Outputs")
Batch8004_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch8004/Outputs")
Batch8005_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch8005/Outputs")
Batch8006_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch8006/Outputs")
Batch8007_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch8007/Outputs")
Batch8008_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch8008/Outputs")
Batch8009_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch8009/Outputs")
Batch8_Genetics <-
  rbind(
    Batch8001_gen_Genetics,
    Batch8002_gen_Genetics,
    Batch8003_gen_Genetics,
    Batch8004_gen_Genetics,
    Batch8005_gen_Genetics,
    Batch8006_gen_Genetics,
    Batch8007_gen_Genetics,
    Batch8008_gen_Genetics,
    Batch8009_gen_Genetics
  )
write.csv(Batch8_Genetics, file = "E:/Models_outputs/data_analysis/Genetics/Batch8_Genetics.csv", row.names = F)



Batch9001_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch9001/Outputs")
Batch9002_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch9002/Outputs")
Batch9003_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch9003/Outputs")
Batch9004_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch9004/Outputs")
Batch9005_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch9005/Outputs")
Batch9006_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch9006/Outputs")
Batch9007_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch9007/Outputs")
Batch9008_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch9008/Outputs")
Batch9009_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch9009/Outputs")
Batch9_Genetics <-
  rbind(
    Batch9001_gen_Genetics,
    Batch9002_gen_Genetics,
    Batch9003_gen_Genetics,
    Batch9004_gen_Genetics,
    Batch9005_gen_Genetics,
    Batch9006_gen_Genetics,
    Batch9007_gen_Genetics,
    Batch9008_gen_Genetics,
    Batch9009_gen_Genetics
  )
write.csv(Batch9_Genetics, file = "E:/Models_outputs/data_analysis/Genetics/Batch9_Genetics.csv", row.names = F)



Batch10001_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch10001/Outputs")
Batch10002_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch10002/Outputs")
Batch10003_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch10003/Outputs")
Batch10004_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch10004/Outputs")
Batch10005_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch10005/Outputs")
Batch10006_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch10006/Outputs")
Batch10007_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch10007/Outputs")
Batch10008_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch10008/Outputs")
Batch10009_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch10009/Outputs")
Batch10_Genetics <-
  rbind(
    Batch10001_gen_Genetics,
    Batch10002_gen_Genetics,
    Batch10003_gen_Genetics,
    Batch10004_gen_Genetics,
    Batch10005_gen_Genetics,
    Batch10006_gen_Genetics,
    Batch10007_gen_Genetics,
    Batch10008_gen_Genetics,
    Batch10009_gen_Genetics
  )
write.csv(Batch10_Genetics, file = "E:/Models_outputs/data_analysis/Genetics/Batch10_Genetics.csv", row.names = F)



Batch11001_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch11001/Outputs")
Batch11002_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch11002/Outputs")
Batch11003_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch11003/Outputs")
Batch11004_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch11004/Outputs")
Batch11005_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch11005/Outputs")
Batch11006_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch11006/Outputs")
Batch11007_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch11007/Outputs")
Batch11008_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch11008/Outputs")
Batch11009_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch11009/Outputs")
Batch11_Genetics <-
  rbind(
    Batch11001_gen_Genetics,
    Batch11002_gen_Genetics,
    Batch11003_gen_Genetics,
    Batch11004_gen_Genetics,
    Batch11005_gen_Genetics,
    Batch11006_gen_Genetics,
    Batch11007_gen_Genetics,
    Batch11008_gen_Genetics,
    Batch11009_gen_Genetics
  )
write.csv(Batch11_Genetics, file = "E:/Models_outputs/data_analysis/Genetics/Batch11_Genetics.csv", row.names = F)



Batch12001_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch12001/Outputs")
Batch12002_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch12002/Outputs")
Batch12003_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch12003/Outputs")
Batch12004_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch12004/Outputs")
Batch12005_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch12005/Outputs")
Batch12006_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch12006/Outputs")
Batch12007_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch12007/Outputs")
Batch12008_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch12008/Outputs")
Batch12009_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch12009/Outputs")
Batch12_Genetics <-
  rbind(
    Batch12001_gen_Genetics,
    Batch12002_gen_Genetics,
    Batch12003_gen_Genetics,
    Batch12004_gen_Genetics,
    Batch12005_gen_Genetics,
    Batch12006_gen_Genetics,
    Batch12007_gen_Genetics,
    Batch12008_gen_Genetics,
    Batch12009_gen_Genetics
  )
write.csv(Batch12_Genetics, file = "E:/Models_outputs/data_analysis/Genetics/Batch12_Genetics.csv", row.names = F)



Batch13001_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch13001/Outputs")
Batch13002_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch13002/Outputs")
Batch13003_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch13003/Outputs")
Batch13004_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch13004/Outputs")
Batch13005_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch13005/Outputs")
Batch13006_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch13006/Outputs")
Batch13007_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch13007/Outputs")
Batch13008_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch13008/Outputs")
Batch13009_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch13009/Outputs")
Batch13_Genetics <-
  rbind(
    Batch13001_gen_Genetics,
    Batch13002_gen_Genetics,
    Batch13003_gen_Genetics,
    Batch13004_gen_Genetics,
    Batch13005_gen_Genetics,
    Batch13006_gen_Genetics,
    Batch13007_gen_Genetics,
    Batch13008_gen_Genetics,
    Batch13009_gen_Genetics
  )
write.csv(Batch13_Genetics, file = "E:/Models_outputs/data_analysis/Genetics/Batch13_Genetics.csv", row.names = F)



Batch14001_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch14001/Outputs")
Batch14002_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch14002/Outputs")
Batch14003_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch14003/Outputs")
Batch14004_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch14004/Outputs")
Batch14005_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch14005/Outputs")
Batch14006_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch14006/Outputs")
Batch14007_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch14007/Outputs")
Batch14008_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch14008/Outputs")
Batch14009_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch14009/Outputs")
Batch14_Genetics <-
  rbind(
    Batch14001_gen_Genetics,
    Batch14002_gen_Genetics,
    Batch14003_gen_Genetics,
    Batch14004_gen_Genetics,
    Batch14005_gen_Genetics,
    Batch14006_gen_Genetics,
    Batch14007_gen_Genetics,
    Batch14008_gen_Genetics,
    Batch14009_gen_Genetics
  )
write.csv(Batch14_Genetics, file = "E:/Models_outputs/data_analysis/Genetics/Batch14_Genetics.csv", row.names = F)


Batch15001_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch15001/Outputs")
Batch15002_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch15002/Outputs")
Batch15003_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch15003/Outputs")
Batch15004_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch15004/Outputs")
Batch15005_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch15005/Outputs")
Batch15006_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch15006/Outputs")
Batch15007_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch15007/Outputs")
Batch15008_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch15008/Outputs")
Batch15009_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch15009/Outputs")
Batch15_Genetics <-
  rbind(
    Batch15001_gen_Genetics,
    Batch15002_gen_Genetics,
    Batch15003_gen_Genetics,
    Batch15004_gen_Genetics,
    Batch15005_gen_Genetics,
    Batch15006_gen_Genetics,
    Batch15007_gen_Genetics,
    Batch15008_gen_Genetics,
    Batch15009_gen_Genetics
  )
write.csv(Batch15_Genetics, file = "E:/Models_outputs/data_analysis/Genetics/Batch15_Genetics.csv", row.names = F)



Batch16001_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch16001/Outputs")
Batch16002_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch16002/Outputs")
Batch16003_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch16003/Outputs")
Batch16004_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch16004/Outputs")
Batch16005_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch16005/Outputs")
Batch16006_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch16006/Outputs")
Batch16007_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch16007/Outputs")
Batch16008_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch16008/Outputs")
Batch16009_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch16009/Outputs")
Batch16_Genetics <-
  rbind(
    Batch16001_gen_Genetics,
    Batch16002_gen_Genetics,
    Batch16003_gen_Genetics,
    Batch16004_gen_Genetics,
    Batch16005_gen_Genetics,
    Batch16006_gen_Genetics,
    Batch16007_gen_Genetics,
    Batch16008_gen_Genetics,
    Batch16009_gen_Genetics
  )
write.csv(Batch16_Genetics, file = "E:/Models_outputs/data_analysis/Genetics/Batch16_Genetics.csv", row.names = F)



Batch17001_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch17001/Outputs")
Batch17002_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch17002/Outputs")
Batch17003_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch17003/Outputs")
Batch17004_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch17004/Outputs")
Batch17005_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch17005/Outputs")
Batch17006_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch17006/Outputs")
Batch17007_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch17007/Outputs")
Batch17008_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch17008/Outputs")
Batch17009_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch17009/Outputs")
Batch17_Genetics <-
  rbind(
    Batch17001_gen_Genetics,
    Batch17002_gen_Genetics,
    Batch17003_gen_Genetics,
    Batch17004_gen_Genetics,
    Batch17005_gen_Genetics,
    Batch17006_gen_Genetics,
    Batch17007_gen_Genetics,
    Batch17008_gen_Genetics,
    Batch17009_gen_Genetics
  )
write.csv(Batch17_Genetics, file = "E:/Models_outputs/data_analysis/Genetics/Batch17_Genetics.csv", row.names = F)



Batch18001_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch18001/Outputs")
Batch18002_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch18002/Outputs")
Batch18003_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch18003/Outputs")
Batch18004_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch18004/Outputs")
Batch18005_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch18005/Outputs")
Batch18006_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch18006/Outputs")
Batch18007_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch18007/Outputs")
Batch18008_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch18008/Outputs")
Batch18009_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch18009/Outputs")
Batch18_Genetics <-
  rbind(
    Batch18001_gen_Genetics,
    Batch18002_gen_Genetics,
    Batch18003_gen_Genetics,
    Batch18004_gen_Genetics,
    Batch18005_gen_Genetics,
    Batch18006_gen_Genetics,
    Batch18007_gen_Genetics,
    Batch18008_gen_Genetics,
    Batch18009_gen_Genetics
  )
write.csv(Batch18_Genetics, file = "E:/Models_outputs/data_analysis/Genetics/Batch18_Genetics.csv", row.names = F)


Batch19001_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch19001/Outputs")
Batch19002_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch19002/Outputs")
Batch19003_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch19003/Outputs")
Batch19004_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch19004/Outputs")
Batch19005_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch19005/Outputs")
Batch19006_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch19006/Outputs")
Batch19007_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch19007/Outputs")
Batch19008_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch19008/Outputs")
Batch19009_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch19009/Outputs")
Batch19_Genetics <-
  rbind(
    Batch19001_gen_Genetics,
    Batch19002_gen_Genetics,
    Batch19003_gen_Genetics,
    Batch19004_gen_Genetics,
    Batch19005_gen_Genetics,
    Batch19006_gen_Genetics,
    Batch19007_gen_Genetics,
    Batch19008_gen_Genetics,
    Batch19009_gen_Genetics
  )
write.csv(Batch19_Genetics, file = "E:/Models_outputs/data_analysis/Genetics/Batch19_Genetics.csv", row.names = F)


Batch20001_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch20001/Outputs")
Batch20002_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch20002/Outputs")
Batch20003_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch20003/Outputs")
Batch20004_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch20004/Outputs")
Batch20005_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch20005/Outputs")
Batch20006_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch20006/Outputs")
Batch20007_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch20007/Outputs")
Batch20008_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch20008/Outputs")
Batch20009_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch20009/Outputs")
Batch20_Genetics <-
  rbind(
    Batch20001_gen_Genetics,
    Batch20002_gen_Genetics,
    Batch20003_gen_Genetics,
    Batch20004_gen_Genetics,
    Batch20005_gen_Genetics,
    Batch20006_gen_Genetics,
    Batch20007_gen_Genetics,
    Batch20008_gen_Genetics,
    Batch20009_gen_Genetics
  )
write.csv(Batch20_Genetics, file = "E:/Models_outputs/data_analysis/Genetics/Batch20_Genetics.csv", row.names = F)


Batch21001_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch21001/Outputs")
Batch21002_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch21002/Outputs")
Batch21003_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch21003/Outputs")
Batch21004_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch21004/Outputs")
Batch21005_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch21005/Outputs")
Batch21006_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch21006/Outputs")
Batch21007_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch21007/Outputs")
Batch21008_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch21008/Outputs")
Batch21009_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch21009/Outputs")
Batch21_Genetics <-
  rbind(
    Batch21001_gen_Genetics,
    Batch21002_gen_Genetics,
    Batch21003_gen_Genetics,
    Batch21004_gen_Genetics,
    Batch21005_gen_Genetics,
    Batch21006_gen_Genetics,
    Batch21007_gen_Genetics,
    Batch21008_gen_Genetics,
    Batch21009_gen_Genetics
  )
write.csv(Batch21_Genetics, file = "E:/Models_outputs/data_analysis/Genetics/Batch21_Genetics.csv", row.names = F)



Batch22001_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch22001/Outputs")
Batch22002_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch22002/Outputs")
Batch22003_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch22003/Outputs")
Batch22004_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch22004/Outputs")
Batch22005_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch22005/Outputs")
Batch22006_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch22006/Outputs")
Batch22007_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch22007/Outputs")
Batch22008_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch22008/Outputs")
Batch22009_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch22009/Outputs")
Batch22_Genetics <-
  rbind(
    Batch22001_gen_Genetics,
    Batch22002_gen_Genetics,
    Batch22003_gen_Genetics,
    Batch22004_gen_Genetics,
    Batch22005_gen_Genetics,
    Batch22006_gen_Genetics,
    Batch22007_gen_Genetics,
    Batch22008_gen_Genetics,
    Batch22009_gen_Genetics
  )
write.csv(Batch22_Genetics, file = "E:/Models_outputs/data_analysis/Genetics/Batch22_Genetics.csv", row.names = F)


Batch23001_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch23001/Outputs")
Batch23002_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch23002/Outputs")
Batch23003_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch23003/Outputs")
Batch23004_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch23004/Outputs")
Batch23005_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch23005/Outputs")
Batch23006_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch23006/Outputs")
Batch23007_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch23007/Outputs")
Batch23008_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch23008/Outputs")
Batch23009_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch23009/Outputs")
Batch23_Genetics <-
  rbind(
    Batch23001_gen_Genetics,
    Batch23002_gen_Genetics,
    Batch23003_gen_Genetics,
    Batch23004_gen_Genetics,
    Batch23005_gen_Genetics,
    Batch23006_gen_Genetics,
    Batch23007_gen_Genetics,
    Batch23008_gen_Genetics,
    Batch23009_gen_Genetics
  )
write.csv(Batch23_Genetics, file = "E:/Models_outputs/data_analysis/Genetics/Batch23_Genetics.csv", row.names = F)



Batch24001_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch24001/Outputs")
Batch24002_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch24002/Outputs")
Batch24003_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch24003/Outputs")
Batch24004_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch24004/Outputs")
Batch24005_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch24005/Outputs")
Batch24006_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch24006/Outputs")
Batch24007_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch24007/Outputs")
Batch24008_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch24008/Outputs")
Batch24009_gen_Genetics <-
  Filter_genetics(min_pop_size = 5, mypath = "E:/Models_outputs/Batch24009/Outputs")
Batch24_Genetics <-
  rbind(
    Batch24001_gen_Genetics,
    Batch24002_gen_Genetics,
    Batch24003_gen_Genetics,
    Batch24004_gen_Genetics,
    Batch24005_gen_Genetics,
    Batch24006_gen_Genetics,
    Batch24007_gen_Genetics,
    Batch24008_gen_Genetics,
    Batch24009_gen_Genetics
  )
write.csv(Batch24_Genetics, file = "E:/Models_outputs/data_analysis/Genetics/Batch24_Genetics.csv", row.names = F)


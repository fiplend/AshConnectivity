# Pop size


Pop_size <- function(mypath, mypattern, ...) {
  setwd(mypath)
  library(dplyr)
  library(stringr)
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
  
  
  population_size <- combined %>%
    group_by(Simulation, Batch, Land, Rep, Year) %>%
    summarise(
      pop_size = sum(NInd))
  
  return(population_size)
}

### Land 1 analysis ####
Batch1001<-Pop_size(mypath = "E:/Models_Aug2020/Batch1001/Outputs", mypattern = "Pop.txt", sep="\t")
Batch1002<-Pop_size(mypath = "E:/Models_Aug2020/Batch1002/Outputs", mypattern = "Pop.txt", sep="\t")
Batch1003<-Pop_size(mypath = "E:/Models_Aug2020/Batch1003/Outputs", mypattern = "Pop.txt", sep="\t")
Batch1004<-Pop_size(mypath = "E:/Models_Aug2020/Batch1004/Outputs", mypattern = "Pop.txt", sep="\t")
Batch1005<-Pop_size(mypath = "E:/Models_Aug2020/Batch1005/Outputs", mypattern = "Pop.txt", sep="\t")
Batch1006<-Pop_size(mypath = "E:/Models_Aug2020/Batch1006/Outputs", mypattern = "Pop.txt", sep="\t")
Batch1007<-Pop_size(mypath = "E:/Models_Aug2020/Batch1007/Outputs", mypattern = "Pop.txt", sep="\t")
Batch1008<-Pop_size(mypath = "E:/Models_Aug2020/Batch1008/Outputs", mypattern = "Pop.txt", sep="\t")
Batch1009<-Pop_size(mypath = "E:/Models_Aug2020/Batch1009/Outputs", mypattern = "Pop.txt", sep="\t")
Land1_pop <- rbind(Batch1001, Batch1002, Batch1003, Batch1004, Batch1005, Batch1006, Batch1007, Batch1008, Batch1009)

write.csv(Land1_pop, file = "E:/Models_Aug2020/data_analysis/pop_size/population_size_land1.csv", row.names = F)


### Land 2 analysis ####
Batch2001<-Pop_size(mypath = "E:/Models_Aug2020/Batch2001/Outputs", mypattern = "Pop.txt", sep="\t")
Batch2002<-Pop_size(mypath = "E:/Models_Aug2020/Batch2002/Outputs", mypattern = "Pop.txt", sep="\t")
Batch2003<-Pop_size(mypath = "E:/Models_Aug2020/Batch2003/Outputs", mypattern = "Pop.txt", sep="\t")
Batch2004<-Pop_size(mypath = "E:/Models_Aug2020/Batch2004/Outputs", mypattern = "Pop.txt", sep="\t")
Batch2005<-Pop_size(mypath = "E:/Models_Aug2020/Batch2005/Outputs", mypattern = "Pop.txt", sep="\t")
Batch2006<-Pop_size(mypath = "E:/Models_Aug2020/Batch2006/Outputs", mypattern = "Pop.txt", sep="\t")
Batch2007<-Pop_size(mypath = "E:/Models_Aug2020/Batch2007/Outputs", mypattern = "Pop.txt", sep="\t")
Batch2008<-Pop_size(mypath = "E:/Models_Aug2020/Batch2008/Outputs", mypattern = "Pop.txt", sep="\t")
Batch2009<-Pop_size(mypath = "E:/Models_Aug2020/Batch2009/Outputs", mypattern = "Pop.txt", sep="\t")
Land2_pop <- rbind(Batch2001, Batch2002, Batch2003, Batch2004, Batch2005, Batch2006, Batch2007, Batch2008, Batch2009)

write.csv(Land2_pop, file = "E:/Models_Aug2020/data_analysis/pop_size/population_size_land2.csv", row.names = F)


### Land 3 analysis ####
Batch3001<-Pop_size(mypath = "E:/Models_Aug2020/Batch3001/Outputs", mypattern = "Pop.txt", sep="\t")
Batch3002<-Pop_size(mypath = "E:/Models_Aug2020/Batch3002/Outputs", mypattern = "Pop.txt", sep="\t")
Batch3003<-Pop_size(mypath = "E:/Models_Aug2020/Batch3003/Outputs", mypattern = "Pop.txt", sep="\t")
Batch3004<-Pop_size(mypath = "E:/Models_Aug2020/Batch3004/Outputs", mypattern = "Pop.txt", sep="\t")
Batch3005<-Pop_size(mypath = "E:/Models_Aug2020/Batch3005/Outputs", mypattern = "Pop.txt", sep="\t")
Batch3006<-Pop_size(mypath = "E:/Models_Aug2020/Batch3006/Outputs", mypattern = "Pop.txt", sep="\t")
Batch3007<-Pop_size(mypath = "E:/Models_Aug2020/Batch3007/Outputs", mypattern = "Pop.txt", sep="\t")
Batch3008<-Pop_size(mypath = "E:/Models_Aug2020/Batch3008/Outputs", mypattern = "Pop.txt", sep="\t")
Batch3009<-Pop_size(mypath = "E:/Models_Aug2020/Batch3009/Outputs", mypattern = "Pop.txt", sep="\t")
Land3_pop <- rbind(Batch3001, Batch3002, Batch3003, Batch3004, Batch3005, Batch3006, Batch3007, Batch3008, Batch3009)

write.csv(Land3_pop, file = "E:/Models_Aug2020/data_analysis/pop_size/population_size_land3.csv", row.names = F)


### Land 4 analysis ####
Batch4001<-Pop_size(mypath = "E:/Models_Aug2020/Batch4001/Outputs", mypattern = "Pop.txt", sep="\t")
Batch4002<-Pop_size(mypath = "E:/Models_Aug2020/Batch4002/Outputs", mypattern = "Pop.txt", sep="\t")
Batch4003<-Pop_size(mypath = "E:/Models_Aug2020/Batch4003/Outputs", mypattern = "Pop.txt", sep="\t")
Batch4004<-Pop_size(mypath = "E:/Models_Aug2020/Batch4004/Outputs", mypattern = "Pop.txt", sep="\t")
Batch4005<-Pop_size(mypath = "E:/Models_Aug2020/Batch4005/Outputs", mypattern = "Pop.txt", sep="\t")
Batch4006<-Pop_size(mypath = "E:/Models_Aug2020/Batch4006/Outputs", mypattern = "Pop.txt", sep="\t")
Batch4007<-Pop_size(mypath = "E:/Models_Aug2020/Batch4007/Outputs", mypattern = "Pop.txt", sep="\t")
Batch4008<-Pop_size(mypath = "E:/Models_Aug2020/Batch4008/Outputs", mypattern = "Pop.txt", sep="\t")
Batch4009<-Pop_size(mypath = "E:/Models_Aug2020/Batch4009/Outputs", mypattern = "Pop.txt", sep="\t")
Land4_pop <- rbind(Batch4001, Batch4002, Batch4003, Batch4004, Batch4005, Batch4006, Batch4007, Batch4008, Batch4009)

write.csv(Land4_pop, file = "E:/Models_Aug2020/data_analysis/pop_size/population_size_land4.csv", row.names = F)


### Land 5 analysis ####
Batch5001<-Pop_size(mypath = "E:/Models_Aug2020/Batch5001/Outputs", mypattern = "Pop.txt", sep="\t")
Batch5002<-Pop_size(mypath = "E:/Models_Aug2020/Batch5002/Outputs", mypattern = "Pop.txt", sep="\t")
Batch5003<-Pop_size(mypath = "E:/Models_Aug2020/Batch5003/Outputs", mypattern = "Pop.txt", sep="\t")
Batch5004<-Pop_size(mypath = "E:/Models_Aug2020/Batch5004/Outputs", mypattern = "Pop.txt", sep="\t")
Batch5005<-Pop_size(mypath = "E:/Models_Aug2020/Batch5005/Outputs", mypattern = "Pop.txt", sep="\t")
Batch5006<-Pop_size(mypath = "E:/Models_Aug2020/Batch5006/Outputs", mypattern = "Pop.txt", sep="\t")
Batch5007<-Pop_size(mypath = "E:/Models_Aug2020/Batch5007/Outputs", mypattern = "Pop.txt", sep="\t")
Batch5008<-Pop_size(mypath = "E:/Models_Aug2020/Batch5008/Outputs", mypattern = "Pop.txt", sep="\t")
Batch5009<-Pop_size(mypath = "E:/Models_Aug2020/Batch5009/Outputs", mypattern = "Pop.txt", sep="\t")
Land5_pop <- rbind(Batch5001, Batch5002, Batch5003, Batch5004, Batch5005, Batch5006, Batch5007, Batch5008, Batch5009)

write.csv(Land5_pop, file = "E:/Models_Aug2020/data_analysis/pop_size/population_size_land5.csv", row.names = F)


### Land 6 analysis ####
Batch6001<-Pop_size(mypath = "E:/Models_Aug2020/Batch6001/Outputs", mypattern = "Pop.txt", sep="\t")
Batch6002<-Pop_size(mypath = "E:/Models_Aug2020/Batch6002/Outputs", mypattern = "Pop.txt", sep="\t")
Batch6003<-Pop_size(mypath = "E:/Models_Aug2020/Batch6003/Outputs", mypattern = "Pop.txt", sep="\t")
Batch6004<-Pop_size(mypath = "E:/Models_Aug2020/Batch6004/Outputs", mypattern = "Pop.txt", sep="\t")
Batch6005<-Pop_size(mypath = "E:/Models_Aug2020/Batch6005/Outputs", mypattern = "Pop.txt", sep="\t")
Batch6006<-Pop_size(mypath = "E:/Models_Aug2020/Batch6006/Outputs", mypattern = "Pop.txt", sep="\t")
Batch6007<-Pop_size(mypath = "E:/Models_Aug2020/Batch6007/Outputs", mypattern = "Pop.txt", sep="\t")
Batch6008<-Pop_size(mypath = "E:/Models_Aug2020/Batch6008/Outputs", mypattern = "Pop.txt", sep="\t")
Batch6009<-Pop_size(mypath = "E:/Models_Aug2020/Batch6009/Outputs", mypattern = "Pop.txt", sep="\t")
Land6_pop <- rbind(Batch6001, Batch6002, Batch6003, Batch6004, Batch6005, Batch6006, Batch6007, Batch6008, Batch6009)

write.csv(Land6_pop, file = "E:/Models_Aug2020/data_analysis/pop_size/population_size_land6.csv", row.names = F)


### Land 7 analysis ####
Batch7001<-Pop_size(mypath = "E:/Models_Aug2020/Batch7001/Outputs", mypattern = "Pop.txt", sep="\t")
Batch7002<-Pop_size(mypath = "E:/Models_Aug2020/Batch7002/Outputs", mypattern = "Pop.txt", sep="\t")
Batch7003<-Pop_size(mypath = "E:/Models_Aug2020/Batch7003/Outputs", mypattern = "Pop.txt", sep="\t")
Batch7004<-Pop_size(mypath = "E:/Models_Aug2020/Batch7004/Outputs", mypattern = "Pop.txt", sep="\t")
Batch7005<-Pop_size(mypath = "E:/Models_Aug2020/Batch7005/Outputs", mypattern = "Pop.txt", sep="\t")
Batch7006<-Pop_size(mypath = "E:/Models_Aug2020/Batch7006/Outputs", mypattern = "Pop.txt", sep="\t")
Batch7007<-Pop_size(mypath = "E:/Models_Aug2020/Batch7007/Outputs", mypattern = "Pop.txt", sep="\t")
Batch7008<-Pop_size(mypath = "E:/Models_Aug2020/Batch7008/Outputs", mypattern = "Pop.txt", sep="\t")
Batch7009<-Pop_size(mypath = "E:/Models_Aug2020/Batch7009/Outputs", mypattern = "Pop.txt", sep="\t")
Land7_pop <- rbind(Batch7001, Batch7002, Batch7003, Batch7004, Batch7005, Batch7006, Batch7007, Batch7008, Batch7009)

write.csv(Land7_pop, file = "E:/Models_Aug2020/data_analysis/pop_size/population_size_land7.csv", row.names = F)

### Land 8 analysis ####
Batch8001<-Pop_size(mypath = "E:/Models_Aug2020/Batch8001/Outputs", mypattern = "Pop.txt", sep="\t")
Batch8002<-Pop_size(mypath = "E:/Models_Aug2020/Batch8002/Outputs", mypattern = "Pop.txt", sep="\t")
Batch8003<-Pop_size(mypath = "E:/Models_Aug2020/Batch8003/Outputs", mypattern = "Pop.txt", sep="\t")
Batch8004<-Pop_size(mypath = "E:/Models_Aug2020/Batch8004/Outputs", mypattern = "Pop.txt", sep="\t")
Batch8005<-Pop_size(mypath = "E:/Models_Aug2020/Batch8005/Outputs", mypattern = "Pop.txt", sep="\t")
Batch8006<-Pop_size(mypath = "E:/Models_Aug2020/Batch8006/Outputs", mypattern = "Pop.txt", sep="\t")
Batch8007<-Pop_size(mypath = "E:/Models_Aug2020/Batch8007/Outputs", mypattern = "Pop.txt", sep="\t")
Batch8008<-Pop_size(mypath = "E:/Models_Aug2020/Batch8008/Outputs", mypattern = "Pop.txt", sep="\t")
Batch8009<-Pop_size(mypath = "E:/Models_Aug2020/Batch8009/Outputs", mypattern = "Pop.txt", sep="\t")
Land8_pop <- rbind(Batch8001, Batch8002, Batch8003, Batch8004, Batch8005, Batch8006, Batch8007, Batch8008, Batch8009)

write.csv(Land8_pop, file = "E:/Models_Aug2020/data_analysis/pop_size/population_size_land8.csv", row.names = F)

### Land 9 analysis ####
Batch9001<-Pop_size(mypath = "E:/Models_Aug2020/Batch9001/Outputs", mypattern = "Pop.txt", sep="\t")
Batch9002<-Pop_size(mypath = "E:/Models_Aug2020/Batch9002/Outputs", mypattern = "Pop.txt", sep="\t")
Batch9003<-Pop_size(mypath = "E:/Models_Aug2020/Batch9003/Outputs", mypattern = "Pop.txt", sep="\t")
Batch9004<-Pop_size(mypath = "E:/Models_Aug2020/Batch9004/Outputs", mypattern = "Pop.txt", sep="\t")
Batch9005<-Pop_size(mypath = "E:/Models_Aug2020/Batch9005/Outputs", mypattern = "Pop.txt", sep="\t")
Batch9006<-Pop_size(mypath = "E:/Models_Aug2020/Batch9006/Outputs", mypattern = "Pop.txt", sep="\t")
Batch9007<-Pop_size(mypath = "E:/Models_Aug2020/Batch9007/Outputs", mypattern = "Pop.txt", sep="\t")
Batch9008<-Pop_size(mypath = "E:/Models_Aug2020/Batch9008/Outputs", mypattern = "Pop.txt", sep="\t")
Batch9009<-Pop_size(mypath = "E:/Models_Aug2020/Batch9009/Outputs", mypattern = "Pop.txt", sep="\t")
Land9_pop <- rbind(Batch9001, Batch9002, Batch9003, Batch9004, Batch9005, Batch9006, Batch9007, Batch9008, Batch9009)

write.csv(Land9_pop, file = "E:/Models_Aug2020/data_analysis/pop_size/population_size_land9.csv", row.names = F)



### Land 10 analysis ####
Batch10001<-Pop_size(mypath = "E:/Models_Aug2020/Batch10001/Outputs", mypattern = "Pop.txt", sep="\t")
Batch10002<-Pop_size(mypath = "E:/Models_Aug2020/Batch10002/Outputs", mypattern = "Pop.txt", sep="\t")
Batch10003<-Pop_size(mypath = "E:/Models_Aug2020/Batch10003/Outputs", mypattern = "Pop.txt", sep="\t")
Batch10004<-Pop_size(mypath = "E:/Models_Aug2020/Batch10004/Outputs", mypattern = "Pop.txt", sep="\t")
Batch10005<-Pop_size(mypath = "E:/Models_Aug2020/Batch10005/Outputs", mypattern = "Pop.txt", sep="\t")
Batch10006<-Pop_size(mypath = "E:/Models_Aug2020/Batch10006/Outputs", mypattern = "Pop.txt", sep="\t")
Batch10007<-Pop_size(mypath = "E:/Models_Aug2020/Batch10007/Outputs", mypattern = "Pop.txt", sep="\t")
Batch10008<-Pop_size(mypath = "E:/Models_Aug2020/Batch10008/Outputs", mypattern = "Pop.txt", sep="\t")
Batch10009<-Pop_size(mypath = "E:/Models_Aug2020/Batch10009/Outputs", mypattern = "Pop.txt", sep="\t")
Land10_pop <- rbind(Batch10001, Batch10002, Batch10003, Batch10004, Batch10005, Batch10006, Batch10007, Batch10008, Batch10009)

write.csv(Land10_pop, file = "E:/Models_Aug2020/data_analysis/pop_size/population_size_land10.csv", row.names = F)

### Land 11 analysis ####
Batch11001<-Pop_size(mypath = "E:/Models_Aug2020/Batch11001/Outputs", mypattern = "Pop.txt", sep="\t")
Batch11002<-Pop_size(mypath = "E:/Models_Aug2020/Batch11002/Outputs", mypattern = "Pop.txt", sep="\t")
Batch11003<-Pop_size(mypath = "E:/Models_Aug2020/Batch11003/Outputs", mypattern = "Pop.txt", sep="\t")
Batch11004<-Pop_size(mypath = "E:/Models_Aug2020/Batch11004/Outputs", mypattern = "Pop.txt", sep="\t")
Batch11005<-Pop_size(mypath = "E:/Models_Aug2020/Batch11005/Outputs", mypattern = "Pop.txt", sep="\t")
Batch11006<-Pop_size(mypath = "E:/Models_Aug2020/Batch11006/Outputs", mypattern = "Pop.txt", sep="\t")
Batch11007<-Pop_size(mypath = "E:/Models_Aug2020/Batch11007/Outputs", mypattern = "Pop.txt", sep="\t")
Batch11008<-Pop_size(mypath = "E:/Models_Aug2020/Batch11008/Outputs", mypattern = "Pop.txt", sep="\t")
Batch11009<-Pop_size(mypath = "E:/Models_Aug2020/Batch11009/Outputs", mypattern = "Pop.txt", sep="\t")
Land11_pop <- rbind(Batch11001, Batch11002, Batch11003, Batch11004, Batch11005, Batch11006, Batch11007, Batch11008, Batch11009)

write.csv(Land11_pop, file = "E:/Models_Aug2020/data_analysis/pop_size/population_size_land11.csv", row.names = F)


### Land 12 analysis ####
Batch12001<-Pop_size(mypath = "E:/Models_Aug2020/Batch12001/Outputs", mypattern = "Pop.txt", sep="\t")
Batch12002<-Pop_size(mypath = "E:/Models_Aug2020/Batch12002/Outputs", mypattern = "Pop.txt", sep="\t")
Batch12003<-Pop_size(mypath = "E:/Models_Aug2020/Batch12003/Outputs", mypattern = "Pop.txt", sep="\t")
Batch12004<-Pop_size(mypath = "E:/Models_Aug2020/Batch12004/Outputs", mypattern = "Pop.txt", sep="\t")
Batch12005<-Pop_size(mypath = "E:/Models_Aug2020/Batch12005/Outputs", mypattern = "Pop.txt", sep="\t")
Batch12006<-Pop_size(mypath = "E:/Models_Aug2020/Batch12006/Outputs", mypattern = "Pop.txt", sep="\t")
Batch12007<-Pop_size(mypath = "E:/Models_Aug2020/Batch12007/Outputs", mypattern = "Pop.txt", sep="\t")
Batch12008<-Pop_size(mypath = "E:/Models_Aug2020/Batch12008/Outputs", mypattern = "Pop.txt", sep="\t")
Batch12009<-Pop_size(mypath = "E:/Models_Aug2020/Batch12009/Outputs", mypattern = "Pop.txt", sep="\t")
Land12_pop <- rbind(Batch12001, Batch12002, Batch12003, Batch12004, Batch12005, Batch12006, Batch12007, Batch12008, Batch12009)

write.csv(Land12_pop, file = "E:/Models_Aug2020/data_analysis/pop_size/population_size_land12.csv", row.names = F)

### Land 13 analysis ####
Batch13001<-Pop_size(mypath = "E:/Models_Aug2020/Batch13001/Outputs", mypattern = "Pop.txt", sep="\t")
Batch13002<-Pop_size(mypath = "E:/Models_Aug2020/Batch13002/Outputs", mypattern = "Pop.txt", sep="\t")
Batch13003<-Pop_size(mypath = "E:/Models_Aug2020/Batch13003/Outputs", mypattern = "Pop.txt", sep="\t")
Batch13004<-Pop_size(mypath = "E:/Models_Aug2020/Batch13004/Outputs", mypattern = "Pop.txt", sep="\t")
Batch13005<-Pop_size(mypath = "E:/Models_Aug2020/Batch13005/Outputs", mypattern = "Pop.txt", sep="\t")
Batch13006<-Pop_size(mypath = "E:/Models_Aug2020/Batch13006/Outputs", mypattern = "Pop.txt", sep="\t")
Batch13007<-Pop_size(mypath = "E:/Models_Aug2020/Batch13007/Outputs", mypattern = "Pop.txt", sep="\t")
Batch13008<-Pop_size(mypath = "E:/Models_Aug2020/Batch13008/Outputs", mypattern = "Pop.txt", sep="\t")
Batch13009<-Pop_size(mypath = "E:/Models_Aug2020/Batch13009/Outputs", mypattern = "Pop.txt", sep="\t")
Land13_pop <- rbind(Batch13001, Batch13002, Batch13003, Batch13004, Batch13005, Batch13006, Batch13007, Batch13008, Batch13009)

write.csv(Land13_pop, file = "E:/Models_Aug2020/data_analysis/pop_size/population_size_land13.csv", row.names = F)

### Land 14 analysis ####
Batch14001<-Pop_size(mypath = "E:/Models_Aug2020/Batch14001/Outputs", mypattern = "Pop.txt", sep="\t")
Batch14002<-Pop_size(mypath = "E:/Models_Aug2020/Batch14002/Outputs", mypattern = "Pop.txt", sep="\t")
Batch14003<-Pop_size(mypath = "E:/Models_Aug2020/Batch14003/Outputs", mypattern = "Pop.txt", sep="\t")
Batch14004<-Pop_size(mypath = "E:/Models_Aug2020/Batch14004/Outputs", mypattern = "Pop.txt", sep="\t")
Batch14005<-Pop_size(mypath = "E:/Models_Aug2020/Batch14005/Outputs", mypattern = "Pop.txt", sep="\t")
Batch14006<-Pop_size(mypath = "E:/Models_Aug2020/Batch14006/Outputs", mypattern = "Pop.txt", sep="\t")
Batch14007<-Pop_size(mypath = "E:/Models_Aug2020/Batch14007/Outputs", mypattern = "Pop.txt", sep="\t")
Batch14008<-Pop_size(mypath = "E:/Models_Aug2020/Batch14008/Outputs", mypattern = "Pop.txt", sep="\t")
Batch14009<-Pop_size(mypath = "E:/Models_Aug2020/Batch14009/Outputs", mypattern = "Pop.txt", sep="\t")
Land14_pop <- rbind(Batch14001, Batch14002, Batch14003, Batch14004, Batch14005, Batch14006, Batch14007, Batch14008, Batch14009)

write.csv(Land14_pop, file = "E:/Models_Aug2020/data_analysis/pop_size/population_size_land14.csv", row.names = F)

### Land 15 analysis ####
Batch15001<-Pop_size(mypath = "E:/Models_Aug2020/Batch15001/Outputs", mypattern = "Pop.txt", sep="\t")
Batch15002<-Pop_size(mypath = "E:/Models_Aug2020/Batch15002/Outputs", mypattern = "Pop.txt", sep="\t")
Batch15003<-Pop_size(mypath = "E:/Models_Aug2020/Batch15003/Outputs", mypattern = "Pop.txt", sep="\t")
Batch15004<-Pop_size(mypath = "E:/Models_Aug2020/Batch15004/Outputs", mypattern = "Pop.txt", sep="\t")
Batch15005<-Pop_size(mypath = "E:/Models_Aug2020/Batch15005/Outputs", mypattern = "Pop.txt", sep="\t")
Batch15006<-Pop_size(mypath = "E:/Models_Aug2020/Batch15006/Outputs", mypattern = "Pop.txt", sep="\t")
Batch15007<-Pop_size(mypath = "E:/Models_Aug2020/Batch15007/Outputs", mypattern = "Pop.txt", sep="\t")
Batch15008<-Pop_size(mypath = "E:/Models_Aug2020/Batch15008/Outputs", mypattern = "Pop.txt", sep="\t")
Batch15009<-Pop_size(mypath = "E:/Models_Aug2020/Batch15009/Outputs", mypattern = "Pop.txt", sep="\t")
Land15_pop <- rbind(Batch15001, Batch15002, Batch15003, Batch15004, Batch15005, Batch15006, Batch15007, Batch15008, Batch15009)

write.csv(Land15_pop, file = "E:/Models_Aug2020/data_analysis/pop_size/population_size_land15.csv", row.names = F)


### Land 16 analysis ####
Batch16001<-Pop_size(mypath = "E:/Models_Aug2020/Batch16001/Outputs", mypattern = "Pop.txt", sep="\t")
Batch16002<-Pop_size(mypath = "E:/Models_Aug2020/Batch16002/Outputs", mypattern = "Pop.txt", sep="\t")
Batch16003<-Pop_size(mypath = "E:/Models_Aug2020/Batch16003/Outputs", mypattern = "Pop.txt", sep="\t")
Batch16004<-Pop_size(mypath = "E:/Models_Aug2020/Batch16004/Outputs", mypattern = "Pop.txt", sep="\t")
Batch16005<-Pop_size(mypath = "E:/Models_Aug2020/Batch16005/Outputs", mypattern = "Pop.txt", sep="\t")
Batch16006<-Pop_size(mypath = "E:/Models_Aug2020/Batch16006/Outputs", mypattern = "Pop.txt", sep="\t")
Batch16007<-Pop_size(mypath = "E:/Models_Aug2020/Batch16007/Outputs", mypattern = "Pop.txt", sep="\t")
Batch16008<-Pop_size(mypath = "E:/Models_Aug2020/Batch16008/Outputs", mypattern = "Pop.txt", sep="\t")
Batch16009<-Pop_size(mypath = "E:/Models_Aug2020/Batch16009/Outputs", mypattern = "Pop.txt", sep="\t")
Land16_pop <- rbind(Batch16001, Batch16002, Batch16003, Batch16004, Batch16005, Batch16006, Batch16007, Batch16008, Batch16009)

write.csv(Land16_pop, file = "E:/Models_Aug2020/data_analysis/pop_size/population_size_land16.csv", row.names = F)



### Land 17 analysis ####
Batch17001<-Pop_size(mypath = "E:/Models_Aug2020/Batch17001/Outputs", mypattern = "Pop.txt", sep="\t")
Batch17002<-Pop_size(mypath = "E:/Models_Aug2020/Batch17002/Outputs", mypattern = "Pop.txt", sep="\t")
Batch17003<-Pop_size(mypath = "E:/Models_Aug2020/Batch17003/Outputs", mypattern = "Pop.txt", sep="\t")
Batch17004<-Pop_size(mypath = "E:/Models_Aug2020/Batch17004/Outputs", mypattern = "Pop.txt", sep="\t")
Batch17005<-Pop_size(mypath = "E:/Models_Aug2020/Batch17005/Outputs", mypattern = "Pop.txt", sep="\t")
Batch17006<-Pop_size(mypath = "E:/Models_Aug2020/Batch17006/Outputs", mypattern = "Pop.txt", sep="\t")
Batch17007<-Pop_size(mypath = "E:/Models_Aug2020/Batch17007/Outputs", mypattern = "Pop.txt", sep="\t")
Batch17008<-Pop_size(mypath = "E:/Models_Aug2020/Batch17008/Outputs", mypattern = "Pop.txt", sep="\t")
Batch17009<-Pop_size(mypath = "E:/Models_Aug2020/Batch17009/Outputs", mypattern = "Pop.txt", sep="\t")
Land17_pop <- rbind(Batch17001, Batch17002, Batch17003, Batch17004, Batch17005, Batch17006, Batch17007, Batch17008, Batch17009)

write.csv(Land17_pop, file = "E:/Models_Aug2020/data_analysis/pop_size/population_size_land17.csv", row.names = F)

### Land 18 analysis ####
Batch18001<-Pop_size(mypath = "E:/Models_Aug2020/Batch18001/Outputs", mypattern = "Pop.txt", sep="\t")
Batch18002<-Pop_size(mypath = "E:/Models_Aug2020/Batch18002/Outputs", mypattern = "Pop.txt", sep="\t")
Batch18003<-Pop_size(mypath = "E:/Models_Aug2020/Batch18003/Outputs", mypattern = "Pop.txt", sep="\t")
Batch18004<-Pop_size(mypath = "E:/Models_Aug2020/Batch18004/Outputs", mypattern = "Pop.txt", sep="\t")
Batch18005<-Pop_size(mypath = "E:/Models_Aug2020/Batch18005/Outputs", mypattern = "Pop.txt", sep="\t")
Batch18006<-Pop_size(mypath = "E:/Models_Aug2020/Batch18006/Outputs", mypattern = "Pop.txt", sep="\t")
Batch18007<-Pop_size(mypath = "E:/Models_Aug2020/Batch18007/Outputs", mypattern = "Pop.txt", sep="\t")
Batch18008<-Pop_size(mypath = "E:/Models_Aug2020/Batch18008/Outputs", mypattern = "Pop.txt", sep="\t")
Batch18009<-Pop_size(mypath = "E:/Models_Aug2020/Batch18009/Outputs", mypattern = "Pop.txt", sep="\t")
Land18_pop <- rbind(Batch18001, Batch18002, Batch18003, Batch18004, Batch18005, Batch18006, Batch18007, Batch18008, Batch18009)

write.csv(Land18_pop, file = "E:/Models_Aug2020/data_analysis/pop_size/population_size_land18.csv", row.names = F)


### Land 19 analysis ####
Batch19001<-Pop_size(mypath = "E:/Models_Aug2020/Batch19001/Outputs", mypattern = "Pop.txt", sep="\t")
Batch19002<-Pop_size(mypath = "E:/Models_Aug2020/Batch19002/Outputs", mypattern = "Pop.txt", sep="\t")
Batch19003<-Pop_size(mypath = "E:/Models_Aug2020/Batch19003/Outputs", mypattern = "Pop.txt", sep="\t")
Batch19004<-Pop_size(mypath = "E:/Models_Aug2020/Batch19004/Outputs", mypattern = "Pop.txt", sep="\t")
Batch19005<-Pop_size(mypath = "E:/Models_Aug2020/Batch19005/Outputs", mypattern = "Pop.txt", sep="\t")
Batch19006<-Pop_size(mypath = "E:/Models_Aug2020/Batch19006/Outputs", mypattern = "Pop.txt", sep="\t")
Batch19007<-Pop_size(mypath = "E:/Models_Aug2020/Batch19007/Outputs", mypattern = "Pop.txt", sep="\t")
Batch19008<-Pop_size(mypath = "E:/Models_Aug2020/Batch19008/Outputs", mypattern = "Pop.txt", sep="\t")
Batch19009<-Pop_size(mypath = "E:/Models_Aug2020/Batch19009/Outputs", mypattern = "Pop.txt", sep="\t")
Land19_pop <- rbind(Batch19001, Batch19002, Batch19003, Batch19004, Batch19005, Batch19006, Batch19007, Batch19008, Batch19009)

write.csv(Land19_pop, file = "E:/Models_Aug2020/data_analysis/pop_size/population_size_land19.csv", row.names = F)



### Land 20 analysis ####
Batch20001<-Pop_size(mypath = "E:/Models_Aug2020/Batch20001/Outputs", mypattern = "Pop.txt", sep="\t")
Batch20002<-Pop_size(mypath = "E:/Models_Aug2020/Batch20002/Outputs", mypattern = "Pop.txt", sep="\t")
Batch20003<-Pop_size(mypath = "E:/Models_Aug2020/Batch20003/Outputs", mypattern = "Pop.txt", sep="\t")
Batch20004<-Pop_size(mypath = "E:/Models_Aug2020/Batch20004/Outputs", mypattern = "Pop.txt", sep="\t")
Batch20005<-Pop_size(mypath = "E:/Models_Aug2020/Batch20005/Outputs", mypattern = "Pop.txt", sep="\t")
Batch20006<-Pop_size(mypath = "E:/Models_Aug2020/Batch20006/Outputs", mypattern = "Pop.txt", sep="\t")
Batch20007<-Pop_size(mypath = "E:/Models_Aug2020/Batch20007/Outputs", mypattern = "Pop.txt", sep="\t")
Batch20008<-Pop_size(mypath = "E:/Models_Aug2020/Batch20008/Outputs", mypattern = "Pop.txt", sep="\t")
Batch20009<-Pop_size(mypath = "E:/Models_Aug2020/Batch20009/Outputs", mypattern = "Pop.txt", sep="\t")
Land20_pop <- rbind(Batch20001, Batch20002, Batch20003, Batch20004, Batch20005, Batch20006, Batch20007, Batch20008, Batch20009)

write.csv(Land20_pop, file = "E:/Models_Aug2020/data_analysis/pop_size/population_size_land20.csv", row.names = F)


### Land 21 analysis ####
Batch21001<-Pop_size(mypath = "E:/Models_Aug2020/Batch21001/Outputs", mypattern = "Pop.txt", sep="\t")
Batch21002<-Pop_size(mypath = "E:/Models_Aug2020/Batch21002/Outputs", mypattern = "Pop.txt", sep="\t")
Batch21003<-Pop_size(mypath = "E:/Models_Aug2020/Batch21003/Outputs", mypattern = "Pop.txt", sep="\t")
Batch21004<-Pop_size(mypath = "E:/Models_Aug2020/Batch21004/Outputs", mypattern = "Pop.txt", sep="\t")
Batch21005<-Pop_size(mypath = "E:/Models_Aug2020/Batch21005/Outputs", mypattern = "Pop.txt", sep="\t")
Batch21006<-Pop_size(mypath = "E:/Models_Aug2020/Batch21006/Outputs", mypattern = "Pop.txt", sep="\t")
Batch21007<-Pop_size(mypath = "E:/Models_Aug2020/Batch21007/Outputs", mypattern = "Pop.txt", sep="\t")
Batch21008<-Pop_size(mypath = "E:/Models_Aug2020/Batch21008/Outputs", mypattern = "Pop.txt", sep="\t")
Batch21009<-Pop_size(mypath = "E:/Models_Aug2020/Batch21009/Outputs", mypattern = "Pop.txt", sep="\t")
Land21_pop <- rbind(Batch21001, Batch21002, Batch21003, Batch21004, Batch21005, Batch21006, Batch21007, Batch21008, Batch21009)

write.csv(Land21_pop, file = "E:/Models_Aug2020/data_analysis/pop_size/population_size_land21.csv", row.names = F)



### Land 22 analysis ####
Batch22001<-Pop_size(mypath = "E:/Models_Aug2020/Batch22001/Outputs", mypattern = "Pop.txt", sep="\t")
Batch22002<-Pop_size(mypath = "E:/Models_Aug2020/Batch22002/Outputs", mypattern = "Pop.txt", sep="\t")
Batch22003<-Pop_size(mypath = "E:/Models_Aug2020/Batch22003/Outputs", mypattern = "Pop.txt", sep="\t")
Batch22004<-Pop_size(mypath = "E:/Models_Aug2020/Batch22004/Outputs", mypattern = "Pop.txt", sep="\t")
Batch22005<-Pop_size(mypath = "E:/Models_Aug2020/Batch22005/Outputs", mypattern = "Pop.txt", sep="\t")
Batch22006<-Pop_size(mypath = "E:/Models_Aug2020/Batch22006/Outputs", mypattern = "Pop.txt", sep="\t")
Batch22007<-Pop_size(mypath = "E:/Models_Aug2020/Batch22007/Outputs", mypattern = "Pop.txt", sep="\t")
Batch22008<-Pop_size(mypath = "E:/Models_Aug2020/Batch22008/Outputs", mypattern = "Pop.txt", sep="\t")
Batch22009<-Pop_size(mypath = "E:/Models_Aug2020/Batch22009/Outputs", mypattern = "Pop.txt", sep="\t")
Land22_pop <- rbind(Batch22001, Batch22002, Batch22003, Batch22004, Batch22005, Batch22006, Batch22007, Batch22008, Batch22009)

write.csv(Land22_pop, file = "E:/Models_Aug2020/data_analysis/pop_size/population_size_land22.csv", row.names = F)


### Land 23 analysis ####
Batch23001<-Pop_size(mypath = "E:/Models_Aug2020/Batch23001/Outputs", mypattern = "Pop.txt", sep="\t")
Batch23002<-Pop_size(mypath = "E:/Models_Aug2020/Batch23002/Outputs", mypattern = "Pop.txt", sep="\t")
Batch23003<-Pop_size(mypath = "E:/Models_Aug2020/Batch23003/Outputs", mypattern = "Pop.txt", sep="\t")
Batch23004<-Pop_size(mypath = "E:/Models_Aug2020/Batch23004/Outputs", mypattern = "Pop.txt", sep="\t")
Batch23005<-Pop_size(mypath = "E:/Models_Aug2020/Batch23005/Outputs", mypattern = "Pop.txt", sep="\t")
Batch23006<-Pop_size(mypath = "E:/Models_Aug2020/Batch23006/Outputs", mypattern = "Pop.txt", sep="\t")
Batch23007<-Pop_size(mypath = "E:/Models_Aug2020/Batch23007/Outputs", mypattern = "Pop.txt", sep="\t")
Batch23008<-Pop_size(mypath = "E:/Models_Aug2020/Batch23008/Outputs", mypattern = "Pop.txt", sep="\t")
Batch23009<-Pop_size(mypath = "E:/Models_Aug2020/Batch23009/Outputs", mypattern = "Pop.txt", sep="\t")
Land23_pop <- rbind(Batch23001, Batch23002, Batch23003, Batch23004, Batch23005, Batch23006, Batch23007, Batch23008, Batch23009)

write.csv(Land23_pop, file = "E:/Models_Aug2020/data_analysis/pop_size/population_size_land23.csv", row.names = F)



### Land 24 analysis ####
Batch24001<-Pop_size(mypath = "E:/Models_Aug2020/Batch24001/Outputs", mypattern = "Pop.txt", sep="\t")
Batch24002<-Pop_size(mypath = "E:/Models_Aug2020/Batch24002/Outputs", mypattern = "Pop.txt", sep="\t")
Batch24003<-Pop_size(mypath = "E:/Models_Aug2020/Batch24003/Outputs", mypattern = "Pop.txt", sep="\t")
Batch24004<-Pop_size(mypath = "E:/Models_Aug2020/Batch24004/Outputs", mypattern = "Pop.txt", sep="\t")
Batch24005<-Pop_size(mypath = "E:/Models_Aug2020/Batch24005/Outputs", mypattern = "Pop.txt", sep="\t")
Batch24006<-Pop_size(mypath = "E:/Models_Aug2020/Batch24006/Outputs", mypattern = "Pop.txt", sep="\t")
Batch24007<-Pop_size(mypath = "E:/Models_Aug2020/Batch24007/Outputs", mypattern = "Pop.txt", sep="\t")
Batch24008<-Pop_size(mypath = "E:/Models_Aug2020/Batch24008/Outputs", mypattern = "Pop.txt", sep="\t")
Batch24009<-Pop_size(mypath = "E:/Models_Aug2020/Batch24009/Outputs", mypattern = "Pop.txt", sep="\t")
Land24_pop <- rbind(Batch24001, Batch24002, Batch24003, Batch24004, Batch24005, Batch24006, Batch24007, Batch24008, Batch24009)

write.csv(Land24_pop, file = "E:/Models_Aug2020/data_analysis/pop_size/population_size_land24.csv", row.names = F)


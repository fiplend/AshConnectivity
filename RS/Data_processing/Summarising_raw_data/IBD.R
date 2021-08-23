# Landscape genetics 

## function for filtering data, running model and extracting data

#IBD
Isolation_by_distance  <- function(min_pop_size, mypath, ...) {
  setwd(mypath)
  library(dplyr)
  library(stringr)
  library(broom)
  import.multiple.files <- function(mypath, pattern, ...)
  {
    setwd(mypath)
    tmp.list.1 <- list.files(mypath, pattern)
    tmp.list.2 <- list(length = length(tmp.list.1))
    for (i in 1:length(tmp.list.1)) {
      tmp.list.2[[i]] <- read.csv(tmp.list.1[i], ...)
    }
    names(tmp.list.2) <- tmp.list.1
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
  
  Genetics$PatchID0_size <-
    Population[match(
      paste(
        Genetics$Rep,
        Genetics$Year,
        Genetics$Simulation,
        Genetics$Batch,
        Genetics$Land,
        Genetics$PatchID0
      ),
      paste(
        Population$Rep,
        Population$Year,
        Population$Simulation,
        Population$Batch,
        Population$Land,
        Population$PatchID
      )
    ), "NInd"]
  Genetics$PatchID1_size <-
    Population[match(
      paste(
        Genetics$Rep,
        Genetics$Year,
        Genetics$Simulation,
        Genetics$Batch,
        Genetics$Land,
        Genetics$PatchID1
      ),
      paste(
        Population$Rep,
        Population$Year,
        Population$Simulation,
        Population$Batch,
        Population$Land,
        Population$PatchID
      )
    ), "NInd"]
  
  
  Genetics <-
    subset(Genetics, Genetics$PatchID0_size > min_pop_size)
  Genetics <-
    subset(Genetics, Genetics$PatchID1_size > min_pop_size)
  
  Genetics$FST <- as.numeric(Genetics$FST)
  Genetics$FST[Genetics$FST[] == -999] <- "NA"
  Genetics$FST[Genetics$FST < 0] <- 0
  Genetics$FST <- as.numeric(Genetics$FST)
  Genetics$Distance_km <- Genetics$Distance / 1000
  
  analysis_estimate <- Genetics %>%
    group_by(Batch, Land, Simulation, Rep, Year) %>%            # group by column z
    do(tidy(lm(FST ~ Distance_km, data = .)))             # group by column z
  # for each group create model using corresponding x and y values
  analysis_estimate$model <- "lm"
  
  
  analysis_glance <- Genetics %>%
    group_by(Batch, Land, Simulation, Rep, Year) %>%            # group by column z
    do(glance(lm(FST ~ Distance_km, data = .)))
  analysis_glance$model <- "lm"
  
  analysis_estimate_quad <- Genetics %>%
    group_by(Batch, Land, Simulation, Rep, Year) %>%            # group by column z
    do(tidy(lm(FST ~ Distance_km + I(Distance_km ^ 2), data = .)))             # group by column z
  # for each group create model using corresponding x and y values
  analysis_estimate_quad$model <- "quadratic"
  
  
  analysis_glance_quad  <- Genetics %>%
    group_by(Batch, Land, Simulation, Rep, Year) %>%            # group by column z
    do(glance(lm(FST ~ Distance_km + I(Distance_km ^ 2), data = .)))
  analysis_glance_quad$model <- "quadratic"
 
  
  Genetic_analysis_lm <-
    merge(
      analysis_estimate,
      analysis_glance,
      by = c("Simulation", "Batch", "Land", "Rep", "Year", "model")
    )
  Genetic_analysis_quad <-
    merge(
      analysis_glance_quad,
      analysis_estimate_quad,
      by = c("Simulation", "Batch", "Land", "Rep", "Year", "model")
    )
  Genetic_analysis <-
    rbind(Genetic_analysis_lm,
          Genetic_analysis_quad)
  
  return(Genetic_analysis)
}
Isolation_by_distance_jost <- function(min_pop_size, mypath, ...) {
  setwd(mypath)
  library(dplyr)
  library(stringr)
  library(broom)
  import.multiple.files <- function(mypath, pattern, ...)
  {
    setwd(mypath)
    tmp.list.1 <- list.files(mypath, pattern)
    tmp.list.2 <- list(length = length(tmp.list.1))
    for (i in 1:length(tmp.list.1)) {
      tmp.list.2[[i]] <- read.csv(tmp.list.1[i], ...)
    }
    names(tmp.list.2) <- tmp.list.1
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
  
  Genetics$PatchID0_size <-
    Population[match(
      paste(
        Genetics$Rep,
        Genetics$Year,
        Genetics$Simulation,
        Genetics$Batch,
        Genetics$Land,
        Genetics$PatchID0
      ),
      paste(
        Population$Rep,
        Population$Year,
        Population$Simulation,
        Population$Batch,
        Population$Land,
        Population$PatchID
      )
    ), "NInd"]
  Genetics$PatchID1_size <-
    Population[match(
      paste(
        Genetics$Rep,
        Genetics$Year,
        Genetics$Simulation,
        Genetics$Batch,
        Genetics$Land,
        Genetics$PatchID1
      ),
      paste(
        Population$Rep,
        Population$Year,
        Population$Simulation,
        Population$Batch,
        Population$Land,
        Population$PatchID
      )
    ), "NInd"]
  
  
  Genetics <-
    subset(Genetics, Genetics$PatchID0_size > min_pop_size)
  Genetics <-
    subset(Genetics, Genetics$PatchID1_size > min_pop_size)
  
  Genetics$MeanD	 <- as.numeric(Genetics$MeanD)
  Genetics$MeanD[Genetics$MeanD == -999] <- "NA"
  Genetics$MeanD[Genetics$MeanD < 0] <- 0
  Genetics$MeanD	 <- as.numeric(Genetics$MeanD)
  Genetics$Distance_km <- Genetics$Distance / 1000
  
  analysis_estimate <- Genetics %>%
    group_by(Batch, Land, Simulation, Rep, Year) %>%            # group by column z
    do(tidy(lm(MeanD	 ~ Distance_km, data = .)))             # group by column z
  # for each group create model using corresponding x and y values
  analysis_estimate$model <- "lm"
  
  
  analysis_glance <- Genetics %>%
    group_by(Batch, Land, Simulation, Rep, Year) %>%            # group by column z
    do(glance(lm(MeanD	 ~ Distance_km, data = .)))
  analysis_glance$model <- "lm"
  
  analysis_estimate_quad <- Genetics %>%
    group_by(Batch, Land, Simulation, Rep, Year) %>%            # group by column z
    do(tidy(lm(MeanD ~ Distance_km + I(Distance_km ^ 2), data = .)))             # group by column z
  # for each group create model using corresponding x and y values
  analysis_estimate_quad$model <- "quadratic"
  
  
  analysis_glance_quad  <- Genetics %>%
    group_by(Batch, Land, Simulation, Rep, Year) %>%            # group by column z
    do(glance(lm(MeanD ~ Distance_km + I(Distance_km ^ 2), data = .)))
  analysis_glance_quad$model <- "quadratic"
  
  
  Genetic_analysis_lm <-
    merge(
      analysis_estimate,
      analysis_glance,
      by = c("Simulation", "Batch", "Land", "Rep", "Year", "model")
    )
  Genetic_analysis_quad <-
    merge(
      analysis_glance_quad,
      analysis_estimate_quad,
      by = c("Simulation", "Batch", "Land", "Rep", "Year", "model")
    )
  Genetic_analysis <-
    rbind(Genetic_analysis_lm,
          Genetic_analysis_quad)
  
  return(Genetic_analysis)
}


Batch1001_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch1001/Outputs")
Batch1002_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch1002/Outputs")
Batch1003_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch1003/Outputs")
Batch1004_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch1004/Outputs")
Batch1005_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch1005/Outputs")
Batch1006_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch1006/Outputs")
Batch1007_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch1007/Outputs")
Batch1008_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch1008/Outputs")
Batch1009_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch1009/Outputs")
Batch1_IBD <-
  rbind(
    Batch1001_gen_IBD,
    Batch1002_gen_IBD,
    Batch1003_gen_IBD,
    Batch1004_gen_IBD,
    Batch1005_gen_IBD,
    Batch1006_gen_IBD,
    Batch1007_gen_IBD,
    Batch1008_gen_IBD,
    Batch1009_gen_IBD
  )
write.csv(Batch1_IBD, file = "E:/Models_Outputs/data_analysis/IBD/Batch1_IBD.csv", row.names = F)





Batch2001_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch2001/Outputs")
Batch2002_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch2002/Outputs")
Batch2003_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch2003/Outputs")
Batch2004_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch2004/Outputs")
Batch2005_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch2005/Outputs")
Batch2006_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch2006/Outputs")
Batch2007_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch2007/Outputs")
Batch2008_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch2008/Outputs")
Batch2009_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch2009/Outputs")
Batch2_IBD <-
  rbind(
    Batch2001_gen_IBD,
    Batch2002_gen_IBD,
    Batch2003_gen_IBD,
    Batch2004_gen_IBD,
    Batch2005_gen_IBD,
    Batch2006_gen_IBD,
    Batch2007_gen_IBD,
    Batch2008_gen_IBD,
    Batch2009_gen_IBD
  )
write.csv(Batch2_IBD, file = "E:/Models_Outputs/data_analysis/IBD/Batch2_IBD.csv", row.names = F)



Batch3001_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch3001/Outputs")
Batch3002_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch3002/Outputs")
Batch3003_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch3003/Outputs")
Batch3004_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch3004/Outputs")
Batch3005_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch3005/Outputs")
Batch3006_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch3006/Outputs")
Batch3007_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch3007/Outputs")
Batch3008_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch3008/Outputs")
Batch3009_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch3009/Outputs")
Batch3_IBD <-
  rbind(
    Batch3001_gen_IBD,
    Batch3002_gen_IBD,
    Batch3003_gen_IBD,
    Batch3004_gen_IBD,
    Batch3005_gen_IBD,
    Batch3006_gen_IBD,
    Batch3007_gen_IBD,
    Batch3008_gen_IBD,
    Batch3009_gen_IBD
  )
write.csv(Batch3_IBD, file = "E:/Models_Outputs/data_analysis/IBD/Batch3_IBD.csv", row.names = F)


Batch4001_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch4001/Outputs")
Batch4002_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch4002/Outputs")
Batch4003_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch4003/Outputs")
Batch4004_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch4004/Outputs")
Batch4005_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch4005/Outputs")
Batch4006_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch4006/Outputs")
Batch4007_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch4007/Outputs")
Batch4008_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch4008/Outputs")
Batch4009_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch4009/Outputs")
Batch4_IBD <-
  rbind(
    Batch4001_gen_IBD,
    Batch4002_gen_IBD,
    Batch4003_gen_IBD,
    Batch4004_gen_IBD,
    Batch4005_gen_IBD,
    Batch4006_gen_IBD,
    Batch4007_gen_IBD,
    Batch4008_gen_IBD,
    Batch4009_gen_IBD
  )
write.csv(Batch4_IBD, file = "E:/Models_Outputs/data_analysis/IBD/Batch4_IBD.csv", row.names = F)


Batch5001_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch5001/Outputs")
Batch5002_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch5002/Outputs")
Batch5003_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch5003/Outputs")
Batch5004_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch5004/Outputs")
Batch5005_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch5005/Outputs")
Batch5006_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch5006/Outputs")
Batch5007_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch5007/Outputs")
Batch5008_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch5008/Outputs")
Batch5009_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch5009/Outputs")
Batch5_IBD <-
  rbind(
    Batch5001_gen_IBD,
    Batch5002_gen_IBD,
    Batch5003_gen_IBD,
    Batch5004_gen_IBD,
    Batch5005_gen_IBD,
    Batch5006_gen_IBD,
    Batch5007_gen_IBD,
    Batch5008_gen_IBD,
    Batch5009_gen_IBD
  )
write.csv(Batch5_IBD, file = "E:/Models_Outputs/data_analysis/IBD/Batch5_IBD.csv", row.names = F)


Batch6001_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch6001/Outputs")
Batch6002_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch6002/Outputs")
Batch6003_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch6003/Outputs")
Batch6004_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch6004/Outputs")
Batch6005_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch6005/Outputs")
Batch6006_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch6006/Outputs")
Batch6007_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch6007/Outputs")
Batch6008_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch6008/Outputs")
Batch6009_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch6009/Outputs")
Batch6_IBD <-
  rbind(
    Batch6001_gen_IBD,
    Batch6002_gen_IBD,
    Batch6003_gen_IBD,
    Batch6004_gen_IBD,
    Batch6005_gen_IBD,
    Batch6006_gen_IBD,
    Batch6007_gen_IBD,
    Batch6008_gen_IBD,
    Batch6009_gen_IBD
  )
write.csv(Batch6_IBD, file = "E:/Models_Outputs/data_analysis/IBD/Batch6_IBD.csv", row.names = F)


Batch7001_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch7001/Outputs")
Batch7002_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch7002/Outputs")
Batch7003_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch7003/Outputs")
Batch7004_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch7004/Outputs")
Batch7005_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch7005/Outputs")
Batch7006_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch7006/Outputs")
Batch7007_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch7007/Outputs")
Batch7008_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch7008/Outputs")
Batch7009_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch7009/Outputs")
Batch7_IBD <-
  rbind(
    Batch7001_gen_IBD,
    Batch7002_gen_IBD,
    Batch7003_gen_IBD,
    Batch7004_gen_IBD,
    Batch7005_gen_IBD,
    Batch7006_gen_IBD,
    Batch7007_gen_IBD,
    Batch7008_gen_IBD,
    Batch7009_gen_IBD
  )
write.csv(Batch7_IBD, file = "E:/Models_Outputs/data_analysis/IBD/Batch7_IBD.csv", row.names = F)


Batch8001_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch8001/Outputs")
Batch8002_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch8002/Outputs")
Batch8003_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch8003/Outputs")
Batch8004_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch8004/Outputs")
Batch8005_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch8005/Outputs")
Batch8006_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch8006/Outputs")
Batch8007_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch8007/Outputs")
Batch8008_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch8008/Outputs")
Batch8009_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch8009/Outputs")
Batch8_IBD <-
  rbind(
    Batch8001_gen_IBD,
    Batch8002_gen_IBD,
    Batch8003_gen_IBD,
    Batch8004_gen_IBD,
    Batch8005_gen_IBD,
    Batch8006_gen_IBD,
    Batch8007_gen_IBD,
    Batch8008_gen_IBD,
    Batch8009_gen_IBD
  )
write.csv(Batch8_IBD, file = "E:/Models_Outputs/data_analysis/IBD/Batch8_IBD.csv", row.names = F)



Batch9001_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch9001/Outputs")
Batch9002_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch9002/Outputs")
Batch9003_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch9003/Outputs")
Batch9004_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch9004/Outputs")
Batch9005_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch9005/Outputs")
Batch9006_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch9006/Outputs")
Batch9007_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch9007/Outputs")
Batch9008_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch9008/Outputs")
Batch9009_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch9009/Outputs")
Batch9_IBD <-
  rbind(
    Batch9001_gen_IBD,
    Batch9002_gen_IBD,
    Batch9003_gen_IBD,
    Batch9004_gen_IBD,
    Batch9005_gen_IBD,
    Batch9006_gen_IBD,
    Batch9007_gen_IBD,
    Batch9008_gen_IBD,
    Batch9009_gen_IBD
  )
write.csv(Batch9_IBD, file = "E:/Models_Outputs/data_analysis/IBD/Batch9_IBD.csv", row.names = F)



Batch10001_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch10001/Outputs")
Batch10002_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch10002/Outputs")
Batch10003_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch10003/Outputs")
Batch10004_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch10004/Outputs")
Batch10005_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch10005/Outputs")
Batch10006_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch10006/Outputs")
Batch10007_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch10007/Outputs")
Batch10008_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch10008/Outputs")
Batch10009_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch10009/Outputs")
Batch10_IBD <-
  rbind(
    Batch10001_gen_IBD,
    Batch10002_gen_IBD,
    Batch10003_gen_IBD,
    Batch10004_gen_IBD,
    Batch10005_gen_IBD,
    Batch10006_gen_IBD,
    Batch10007_gen_IBD,
    Batch10008_gen_IBD,
    Batch10009_gen_IBD
  )
write.csv(Batch10_IBD, file = "E:/Models_Outputs/data_analysis/IBD/Batch10_IBD.csv", row.names = F)



Batch11001_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch11001/Outputs")
Batch11002_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch11002/Outputs")
Batch11003_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch11003/Outputs")
Batch11004_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch11004/Outputs")
Batch11005_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch11005/Outputs")
Batch11006_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch11006/Outputs")
Batch11007_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch11007/Outputs")
Batch11008_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch11008/Outputs")
Batch11009_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch11009/Outputs")
Batch11_IBD <-
  rbind(
    Batch11001_gen_IBD,
    Batch11002_gen_IBD,
    Batch11003_gen_IBD,
    Batch11004_gen_IBD,
    Batch11005_gen_IBD,
    Batch11006_gen_IBD,
    Batch11007_gen_IBD,
    Batch11008_gen_IBD,
    Batch11009_gen_IBD
  )
write.csv(Batch11_IBD, file = "E:/Models_Outputs/data_analysis/IBD/Batch11_IBD.csv", row.names = F)



Batch12001_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch12001/Outputs")
Batch12002_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch12002/Outputs")
Batch12003_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch12003/Outputs")
Batch12004_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch12004/Outputs")
Batch12005_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch12005/Outputs")
Batch12006_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch12006/Outputs")
Batch12007_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch12007/Outputs")
Batch12008_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch12008/Outputs")
Batch12009_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch12009/Outputs")
Batch12_IBD <-
  rbind(
    Batch12001_gen_IBD,
    Batch12002_gen_IBD,
    Batch12003_gen_IBD,
    Batch12004_gen_IBD,
    Batch12005_gen_IBD,
    Batch12006_gen_IBD,
    Batch12007_gen_IBD,
    Batch12008_gen_IBD,
    Batch12009_gen_IBD
  )
write.csv(Batch12_IBD, file = "E:/Models_Outputs/data_analysis/IBD/Batch12_IBD.csv", row.names = F)



Batch13001_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch13001/Outputs")
Batch13002_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch13002/Outputs")
Batch13003_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch13003/Outputs")
Batch13004_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch13004/Outputs")
Batch13005_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch13005/Outputs")
Batch13006_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch13006/Outputs")
Batch13007_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch13007/Outputs")
Batch13008_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch13008/Outputs")
Batch13009_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch13009/Outputs")
Batch13_IBD <-
  rbind(
    Batch13001_gen_IBD,
    Batch13002_gen_IBD,
    Batch13003_gen_IBD,
    Batch13004_gen_IBD,
    Batch13005_gen_IBD,
    Batch13006_gen_IBD,
    Batch13007_gen_IBD,
    Batch13008_gen_IBD,
    Batch13009_gen_IBD
  )
write.csv(Batch13_IBD, file = "E:/Models_Outputs/data_analysis/IBD/Batch13_IBD.csv", row.names = F)



Batch14001_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch14001/Outputs")
Batch14002_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch14002/Outputs")
Batch14003_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch14003/Outputs")
Batch14004_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch14004/Outputs")
Batch14005_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch14005/Outputs")
Batch14006_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch14006/Outputs")
Batch14007_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch14007/Outputs")
Batch14008_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch14008/Outputs")
Batch14009_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch14009/Outputs")
Batch14_IBD <-
  rbind(
    Batch14001_gen_IBD,
    Batch14002_gen_IBD,
    Batch14003_gen_IBD,
    Batch14004_gen_IBD,
    Batch14005_gen_IBD,
    Batch14006_gen_IBD,
    Batch14007_gen_IBD,
    Batch14008_gen_IBD,
    Batch14009_gen_IBD
  )
write.csv(Batch14_IBD, file = "E:/Models_Outputs/data_analysis/IBD/Batch14_IBD.csv", row.names = F)


Batch15001_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch15001/Outputs")
Batch15002_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch15002/Outputs")
Batch15003_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch15003/Outputs")
Batch15004_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch15004/Outputs")
Batch15005_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch15005/Outputs")
Batch15006_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch15006/Outputs")
Batch15007_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch15007/Outputs")
Batch15008_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch15008/Outputs")
Batch15009_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch15009/Outputs")
Batch15_IBD <-
  rbind(
    Batch15001_gen_IBD,
    Batch15002_gen_IBD,
    Batch15003_gen_IBD,
    Batch15004_gen_IBD,
    Batch15005_gen_IBD,
    Batch15006_gen_IBD,
    Batch15007_gen_IBD,
    Batch15008_gen_IBD,
    Batch15009_gen_IBD
  )
write.csv(Batch15_IBD, file = "E:/Models_Outputs/data_analysis/IBD/Batch15_IBD.csv", row.names = F)



Batch16001_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch16001/Outputs")
Batch16002_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch16002/Outputs")
Batch16003_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch16003/Outputs")
Batch16004_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch16004/Outputs")
Batch16005_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch16005/Outputs")
Batch16006_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch16006/Outputs")
Batch16007_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch16007/Outputs")
Batch16008_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch16008/Outputs")
Batch16009_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch16009/Outputs")
Batch16_IBD <-
  rbind(
    Batch16001_gen_IBD,
    Batch16002_gen_IBD,
    Batch16003_gen_IBD,
    Batch16004_gen_IBD,
    Batch16005_gen_IBD,
    Batch16006_gen_IBD,
    Batch16007_gen_IBD,
    Batch16008_gen_IBD,
    Batch16009_gen_IBD
  )
write.csv(Batch16_IBD, file = "E:/Models_Outputs/data_analysis/IBD/Batch16_IBD.csv", row.names = F)



Batch17001_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch17001/Outputs")
Batch17002_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch17002/Outputs")
Batch17003_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch17003/Outputs")
Batch17004_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch17004/Outputs")
Batch17005_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch17005/Outputs")
Batch17006_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch17006/Outputs")
Batch17007_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch17007/Outputs")
Batch17008_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch17008/Outputs")
Batch17009_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch17009/Outputs")
Batch17_IBD <-
  rbind(
    Batch17001_gen_IBD,
    Batch17002_gen_IBD,
    Batch17003_gen_IBD,
    Batch17004_gen_IBD,
    Batch17005_gen_IBD,
    Batch17006_gen_IBD,
    Batch17007_gen_IBD,
    Batch17008_gen_IBD,
    Batch17009_gen_IBD
  )
write.csv(Batch17_IBD, file = "E:/Models_Outputs/data_analysis/IBD/Batch17_IBD.csv", row.names = F)



Batch18001_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch18001/Outputs")
Batch18002_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch18002/Outputs")
Batch18003_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch18003/Outputs")
Batch18004_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch18004/Outputs")
Batch18005_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch18005/Outputs")
Batch18006_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch18006/Outputs")
Batch18007_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch18007/Outputs")
Batch18008_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch18008/Outputs")
Batch18009_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch18009/Outputs")
Batch18_IBD <-
  rbind(
    Batch18001_gen_IBD,
    Batch18002_gen_IBD,
    Batch18003_gen_IBD,
    Batch18004_gen_IBD,
    Batch18005_gen_IBD,
    Batch18006_gen_IBD,
    Batch18007_gen_IBD,
    Batch18008_gen_IBD,
    Batch18009_gen_IBD
  )
write.csv(Batch18_IBD, file = "E:/Models_Outputs/data_analysis/IBD/Batch18_IBD.csv", row.names = F)


Batch19001_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch19001/Outputs")
Batch19002_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch19002/Outputs")
Batch19003_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch19003/Outputs")
Batch19004_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch19004/Outputs")
Batch19005_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch19005/Outputs")
Batch19006_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch19006/Outputs")
Batch19007_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch19007/Outputs")
Batch19008_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch19008/Outputs")
Batch19009_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Batch19009/Outputs")
Batch19_IBD <-
  rbind(
    Batch19001_gen_IBD,
    Batch19002_gen_IBD,
    Batch19003_gen_IBD,
    Batch19004_gen_IBD,
    Batch19005_gen_IBD,
    Batch19006_gen_IBD,
    Batch19007_gen_IBD,
    Batch19008_gen_IBD,
    Batch19009_gen_IBD
  )
write.csv(Batch19_IBD, file = "E:/Models_Outputs/data_analysis/IBD/Batch19_IBD.csv", row.names = F)


Batch20001_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Land20/Batch20001/Outputs")
Batch20002_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Land20/Batch20002/Outputs")
Batch20003_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Land20/Batch20003/Outputs")
Batch20004_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Land20/Batch20004/Outputs")
Batch20005_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Land20/Batch20005/Outputs")
Batch20006_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Land20/Batch20006/Outputs")
Batch20007_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Land20/Batch20007/Outputs")
Batch20008_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Land20/Batch20008/Outputs")
Batch20009_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Land20/Batch20009/Outputs")
Batch20_IBD <-
  rbind(
    Batch20001_gen_IBD,
    Batch20002_gen_IBD,
    Batch20003_gen_IBD,
    Batch20004_gen_IBD,
    Batch20005_gen_IBD,
    Batch20006_gen_IBD,
    Batch20007_gen_IBD,
    Batch20008_gen_IBD,
    Batch20009_gen_IBD
  )
write.csv(Batch20_IBD, file = "E:/Models_Outputs/data_analysis/IBD/Batch20_IBD.csv", row.names = F)


Batch21001_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Land21/Batch21001/Outputs")
Batch21002_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Land21/Batch21002/Outputs")
Batch21003_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Land21/Batch21003/Outputs")
Batch21004_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Land21/Batch21004/Outputs")
Batch21005_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Land21/Batch21005/Outputs")
Batch21006_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Land21/Batch21006/Outputs")
Batch21007_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Land21/Batch21007/Outputs")
Batch21008_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Land21/Batch21008/Outputs")
Batch21009_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Land21/Batch21009/Outputs")
Batch21_IBD <-
  rbind(
    Batch21001_gen_IBD,
    Batch21002_gen_IBD,
    Batch21003_gen_IBD,
    Batch21004_gen_IBD,
    Batch21005_gen_IBD,
    Batch21006_gen_IBD,
    Batch21007_gen_IBD,
    Batch21008_gen_IBD,
    Batch21009_gen_IBD
  )
write.csv(Batch21_IBD, file = "E:/Models_Outputs/data_analysis/IBD/Batch21_IBD.csv", row.names = F)



Batch22001_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Land22/Batch22001/Outputs")
Batch22002_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Land22/Batch22002/Outputs")
Batch22003_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Land22/Batch22003/Outputs")
Batch22004_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Land22/Batch22004/Outputs")
Batch22005_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Land22/Batch22005/Outputs")
Batch22006_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Land22/Batch22006/Outputs")
Batch22007_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Land22/Batch22007/Outputs")
Batch22008_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Land22/Batch22008/Outputs")
Batch22009_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Land22/Batch22009/Outputs")
Batch22_IBD <-
  rbind(
    Batch22001_gen_IBD,
    Batch22002_gen_IBD,
    Batch22003_gen_IBD,
    Batch22004_gen_IBD,
    Batch22005_gen_IBD,
    Batch22006_gen_IBD,
    Batch22007_gen_IBD,
    Batch22008_gen_IBD,
    Batch22009_gen_IBD
  )
write.csv(Batch22_IBD, file = "E:/Models_Outputs/data_analysis/IBD/Batch22_IBD.csv", row.names = F)


Batch23001_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Land23/Batch23001/Outputs")
Batch23002_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Land23/Batch23002/Outputs")
Batch23003_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Land23/Batch23003/Outputs")
Batch23004_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Land23/Batch23004/Outputs")
Batch23005_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Land23/Batch23005/Outputs")
Batch23006_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Land23/Batch23006/Outputs")
Batch23007_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Land23/Batch23007/Outputs")
Batch23008_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Land23/Batch23008/Outputs")
Batch23009_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Land23/Batch23009/Outputs")
Batch23_IBD <-
  rbind(
    Batch23001_gen_IBD,
    Batch23002_gen_IBD,
    Batch23003_gen_IBD,
    Batch23004_gen_IBD,
    Batch23005_gen_IBD,
    Batch23006_gen_IBD,
    Batch23007_gen_IBD,
    Batch23008_gen_IBD,
    Batch23009_gen_IBD
  )
write.csv(Batch23_IBD, file = "E:/Models_Outputs/data_analysis/IBD/Batch23_IBD.csv", row.names = F)



Batch24001_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Land24/Batch24001/Outputs")
Batch24002_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Land24/Batch24002/Outputs")
Batch24003_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Land24/Batch24003/Outputs")
Batch24004_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Land24/Batch24004/Outputs")
Batch24005_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Land24/Batch24005/Outputs")
Batch24006_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Land24/Batch24006/Outputs")
Batch24007_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Land24/Batch24007/Outputs")
Batch24008_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Land24/Batch24008/Outputs")
Batch24009_gen_IBD <-
  Isolation_by_distance (min_pop_size = 5, mypath = "E:/Models_Outputs/Land24/Batch24009/Outputs")
Batch24_IBD <-
  rbind(
    Batch24001_gen_IBD,
    Batch24002_gen_IBD,
    Batch24003_gen_IBD,
    Batch24004_gen_IBD,
    Batch24005_gen_IBD,
    Batch24006_gen_IBD,
    Batch24007_gen_IBD,
    Batch24008_gen_IBD,
    Batch24009_gen_IBD
  )
write.csv(Batch24_IBD, file = "E:/Models_Outputs/data_analysis/IBD/Batch24_IBD.csv", row.names = F)


IBD <-
  rbind(
    Batch1_IBD,
    Batch2_IBD,
    Batch3_IBD,
    Batch4_IBD,
    Batch5_IBD,
    Batch6_IBD,
    Batch7_IBD,
    Batch8_IBD,
    Batch9_IBD,
    Batch10_IBD,
    Batch11_IBD,
    Batch12_IBD,
    Batch13_IBD,
    Batch14_IBD,
    Batch15_IBD,
    Batch16_IBD,
    Batch17_IBD,
    Batch18_IBD,
    Batch19_IBD,
    Batch20_IBD,
    Batch21_IBD,
    Batch22_IBD,
    Batch23_IBD,
    Batch24_IBD
  )


# correltion co ####

CC <- function(min_pop_size, mypath, ...) {
  setwd(mypath)
  library(dplyr)
  library(stringr)
  library(broom)
  import.multiple.files <- function(mypath, pattern, ...)
  {
    setwd(mypath)
    tmp.list.1 <- list.files(mypath, pattern)
    tmp.list.2 <- list(length = length(tmp.list.1))
    for (i in 1:length(tmp.list.1)) {
      tmp.list.2[[i]] <- read.csv(tmp.list.1[i], ...)
    }
    names(tmp.list.2) <- tmp.list.1
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
  
  Genetics$PatchID0_size <-
    Population[match(
      paste(
        Genetics$Rep,
        Genetics$Year,
        Genetics$Simulation,
        Genetics$Batch,
        Genetics$Land,
        Genetics$PatchID0
      ),
      paste(
        Population$Rep,
        Population$Year,
        Population$Simulation,
        Population$Batch,
        Population$Land,
        Population$PatchID
      )
    ), "NInd"]
  Genetics$PatchID1_size <-
    Population[match(
      paste(
        Genetics$Rep,
        Genetics$Year,
        Genetics$Simulation,
        Genetics$Batch,
        Genetics$Land,
        Genetics$PatchID1
      ),
      paste(
        Population$Rep,
        Population$Year,
        Population$Simulation,
        Population$Batch,
        Population$Land,
        Population$PatchID
      )
    ), "NInd"]
  
  
  Genetics <-
    subset(Genetics, Genetics$PatchID0_size > min_pop_size)
  Genetics <-
    subset(Genetics, Genetics$PatchID1_size > min_pop_size)
  
  Genetics$FST <- as.numeric(Genetics$FST)
  Genetics$FST[Genetics$FST[] == -999] <- "NA"
  Genetics$FST[Genetics$FST < 0] <- 0
  Genetics$FST <- as.numeric(Genetics$FST)
  
  
  analysis_glance <- Genetics %>%
    group_by(Batch, Land, Simulation, Rep, Year) %>%
    summarize(COR_FST=cor(Distance,FST))
 
    return(analysis_glance)
}




Batch1001_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch1001/Outputs")
Batch1002_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch1002/Outputs")
Batch1003_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch1003/Outputs")
Batch1004_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch1004/Outputs")
Batch1005_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch1005/Outputs")
Batch1006_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch1006/Outputs")
Batch1007_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch1007/Outputs")
Batch1008_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch1008/Outputs")
Batch1009_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch1009/Outputs")
Batch1_CC <-
  rbind(
    Batch1001_gen_CC,
    Batch1002_gen_CC,
    Batch1003_gen_CC,
    Batch1004_gen_CC,
    Batch1005_gen_CC,
    Batch1006_gen_CC,
    Batch1007_gen_CC,
    Batch1008_gen_CC,
    Batch1009_gen_CC
  )
write.csv(Batch1_CC, file = "E:/Models_Outputs/data_analysis/CC/Batch1_CC.csv", row.names = F)





Batch2001_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch2001/Outputs")
Batch2002_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch2002/Outputs")
Batch2003_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch2003/Outputs")
Batch2004_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch2004/Outputs")
Batch2005_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch2005/Outputs")
Batch2006_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch2006/Outputs")
Batch2007_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch2007/Outputs")
Batch2008_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch2008/Outputs")
Batch2009_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch2009/Outputs")
Batch2_CC <-
  rbind(
    Batch2001_gen_CC,
    Batch2002_gen_CC,
    Batch2003_gen_CC,
    Batch2004_gen_CC,
    Batch2005_gen_CC,
    Batch2006_gen_CC,
    Batch2007_gen_CC,
    Batch2008_gen_CC,
    Batch2009_gen_CC
  )
write.csv(Batch2_CC, file = "E:/Models_Outputs/data_analysis/CC/Batch2_CC.csv", row.names = F)



Batch3001_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch3001/Outputs")
Batch3002_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch3002/Outputs")
Batch3003_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch3003/Outputs")
Batch3004_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch3004/Outputs")
Batch3005_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch3005/Outputs")
Batch3006_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch3006/Outputs")
Batch3007_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch3007/Outputs")
Batch3008_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch3008/Outputs")
Batch3009_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch3009/Outputs")
Batch3_CC <-
  rbind(
    Batch3001_gen_CC,
    Batch3002_gen_CC,
    Batch3003_gen_CC,
    Batch3004_gen_CC,
    Batch3005_gen_CC,
    Batch3006_gen_CC,
    Batch3007_gen_CC,
    Batch3008_gen_CC,
    Batch3009_gen_CC
  )
write.csv(Batch3_CC, file = "E:/Models_Outputs/data_analysis/CC/Batch3_CC.csv", row.names = F)


Batch4001_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch4001/Outputs")
Batch4002_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch4002/Outputs")
Batch4003_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch4003/Outputs")
Batch4004_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch4004/Outputs")
Batch4005_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch4005/Outputs")
Batch4006_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch4006/Outputs")
Batch4007_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch4007/Outputs")
Batch4008_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch4008/Outputs")
Batch4009_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch4009/Outputs")
Batch4_CC <-
  rbind(
    Batch4001_gen_CC,
    Batch4002_gen_CC,
    Batch4003_gen_CC,
    Batch4004_gen_CC,
    Batch4005_gen_CC,
    Batch4006_gen_CC,
    Batch4007_gen_CC,
    Batch4008_gen_CC,
    Batch4009_gen_CC
  )
write.csv(Batch4_CC, file = "E:/Models_Outputs/data_analysis/CC/Batch4_CC.csv", row.names = F)


Batch5001_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch5001/Outputs")
Batch5002_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch5002/Outputs")
Batch5003_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch5003/Outputs")
Batch5004_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch5004/Outputs")
Batch5005_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch5005/Outputs")
Batch5006_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch5006/Outputs")
Batch5007_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch5007/Outputs")
Batch5008_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch5008/Outputs")
Batch5009_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch5009/Outputs")
Batch5_CC <-
  rbind(
    Batch5001_gen_CC,
    Batch5002_gen_CC,
    Batch5003_gen_CC,
    Batch5004_gen_CC,
    Batch5005_gen_CC,
    Batch5006_gen_CC,
    Batch5007_gen_CC,
    Batch5008_gen_CC,
    Batch5009_gen_CC
  )
write.csv(Batch5_CC, file = "E:/Models_Outputs/data_analysis/CC/Batch5_CC.csv", row.names = F)


Batch6001_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch6001/Outputs")
Batch6002_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch6002/Outputs")
Batch6003_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch6003/Outputs")
Batch6004_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch6004/Outputs")
Batch6005_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch6005/Outputs")
Batch6006_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch6006/Outputs")
Batch6007_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch6007/Outputs")
Batch6008_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch6008/Outputs")
Batch6009_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch6009/Outputs")
Batch6_CC <-
  rbind(
    Batch6001_gen_CC,
    Batch6002_gen_CC,
    Batch6003_gen_CC,
    Batch6004_gen_CC,
    Batch6005_gen_CC,
    Batch6006_gen_CC,
    Batch6007_gen_CC,
    Batch6008_gen_CC,
    Batch6009_gen_CC
  )
write.csv(Batch6_CC, file = "E:/Models_Outputs/data_analysis/CC/Batch6_CC.csv", row.names = F)


Batch7001_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch7001/Outputs")
Batch7002_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch7002/Outputs")
Batch7003_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch7003/Outputs")
Batch7004_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch7004/Outputs")
Batch7005_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch7005/Outputs")
Batch7006_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch7006/Outputs")
Batch7007_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch7007/Outputs")
Batch7008_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch7008/Outputs")
Batch7009_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch7009/Outputs")
Batch7_CC <-
  rbind(
    Batch7001_gen_CC,
    Batch7002_gen_CC,
    Batch7003_gen_CC,
    Batch7004_gen_CC,
    Batch7005_gen_CC,
    Batch7006_gen_CC,
    Batch7007_gen_CC,
    Batch7008_gen_CC,
    Batch7009_gen_CC
  )
write.csv(Batch7_CC, file = "E:/Models_Outputs/data_analysis/CC/Batch7_CC.csv", row.names = F)


Batch8001_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch8001/Outputs")
Batch8002_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch8002/Outputs")
Batch8003_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch8003/Outputs")
Batch8004_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch8004/Outputs")
Batch8005_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch8005/Outputs")
Batch8006_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch8006/Outputs")
Batch8007_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch8007/Outputs")
Batch8008_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch8008/Outputs")
Batch8009_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch8009/Outputs")
Batch8_CC <-
  rbind(
    Batch8001_gen_CC,
    Batch8002_gen_CC,
    Batch8003_gen_CC,
    Batch8004_gen_CC,
    Batch8005_gen_CC,
    Batch8006_gen_CC,
    Batch8007_gen_CC,
    Batch8008_gen_CC,
    Batch8009_gen_CC
  )
write.csv(Batch8_CC, file = "E:/Models_Outputs/data_analysis/CC/Batch8_CC.csv", row.names = F)



Batch9001_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch9001/Outputs")
Batch9002_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch9002/Outputs")
Batch9003_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch9003/Outputs")
Batch9004_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch9004/Outputs")
Batch9005_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch9005/Outputs")
Batch9006_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch9006/Outputs")
Batch9007_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch9007/Outputs")
Batch9008_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch9008/Outputs")
Batch9009_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch9009/Outputs")
Batch9_CC <-
  rbind(
    Batch9001_gen_CC,
    Batch9002_gen_CC,
    Batch9003_gen_CC,
    Batch9004_gen_CC,
    Batch9005_gen_CC,
    Batch9006_gen_CC,
    Batch9007_gen_CC,
    Batch9008_gen_CC,
    Batch9009_gen_CC
  )
write.csv(Batch9_CC, file = "E:/Models_Outputs/data_analysis/CC/Batch9_CC.csv", row.names = F)



Batch10001_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch10001/Outputs")
Batch10002_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch10002/Outputs")
Batch10003_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch10003/Outputs")
Batch10004_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch10004/Outputs")
Batch10005_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch10005/Outputs")
Batch10006_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch10006/Outputs")
Batch10007_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch10007/Outputs")
Batch10008_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch10008/Outputs")
Batch10009_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch10009/Outputs")
Batch10_CC <-
  rbind(
    Batch10001_gen_CC,
    Batch10002_gen_CC,
    Batch10003_gen_CC,
    Batch10004_gen_CC,
    Batch10005_gen_CC,
    Batch10006_gen_CC,
    Batch10007_gen_CC,
    Batch10008_gen_CC,
    Batch10009_gen_CC
  )
write.csv(Batch10_CC, file = "E:/Models_Outputs/data_analysis/CC/Batch10_CC.csv", row.names = F)



Batch11001_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch11001/Outputs")
Batch11002_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch11002/Outputs")
Batch11003_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch11003/Outputs")
Batch11004_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch11004/Outputs")
Batch11005_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch11005/Outputs")
Batch11006_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch11006/Outputs")
Batch11007_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch11007/Outputs")
Batch11008_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch11008/Outputs")
Batch11009_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch11009/Outputs")
Batch11_CC <-
  rbind(
    Batch11001_gen_CC,
    Batch11002_gen_CC,
    Batch11003_gen_CC,
    Batch11004_gen_CC,
    Batch11005_gen_CC,
    Batch11006_gen_CC,
    Batch11007_gen_CC,
    Batch11008_gen_CC,
    Batch11009_gen_CC
  )
write.csv(Batch11_CC, file = "E:/Models_Outputs/data_analysis/CC/Batch11_CC.csv", row.names = F)



Batch12001_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch12001/Outputs")
Batch12002_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch12002/Outputs")
Batch12003_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch12003/Outputs")
Batch12004_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch12004/Outputs")
Batch12005_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch12005/Outputs")
Batch12006_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch12006/Outputs")
Batch12007_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch12007/Outputs")
Batch12008_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch12008/Outputs")
Batch12009_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch12009/Outputs")
Batch12_CC <-
  rbind(
    Batch12001_gen_CC,
    Batch12002_gen_CC,
    Batch12003_gen_CC,
    Batch12004_gen_CC,
    Batch12005_gen_CC,
    Batch12006_gen_CC,
    Batch12007_gen_CC,
    Batch12008_gen_CC,
    Batch12009_gen_CC
  )
write.csv(Batch12_CC, file = "E:/Models_Outputs/data_analysis/CC/Batch12_CC.csv", row.names = F)



Batch13001_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch13001/Outputs")
Batch13002_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch13002/Outputs")
Batch13003_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch13003/Outputs")
Batch13004_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch13004/Outputs")
Batch13005_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch13005/Outputs")
Batch13006_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch13006/Outputs")
Batch13007_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch13007/Outputs")
Batch13008_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch13008/Outputs")
Batch13009_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch13009/Outputs")
Batch13_CC <-
  rbind(
    Batch13001_gen_CC,
    Batch13002_gen_CC,
    Batch13003_gen_CC,
    Batch13004_gen_CC,
    Batch13005_gen_CC,
    Batch13006_gen_CC,
    Batch13007_gen_CC,
    Batch13008_gen_CC,
    Batch13009_gen_CC
  )
write.csv(Batch13_CC, file = "E:/Models_Outputs/data_analysis/CC/Batch13_CC.csv", row.names = F)



Batch14001_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch14001/Outputs")
Batch14002_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch14002/Outputs")
Batch14003_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch14003/Outputs")
Batch14004_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch14004/Outputs")
Batch14005_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch14005/Outputs")
Batch14006_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch14006/Outputs")
Batch14007_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch14007/Outputs")
Batch14008_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch14008/Outputs")
Batch14009_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch14009/Outputs")
Batch14_CC <-
  rbind(
    Batch14001_gen_CC,
    Batch14002_gen_CC,
    Batch14003_gen_CC,
    Batch14004_gen_CC,
    Batch14005_gen_CC,
    Batch14006_gen_CC,
    Batch14007_gen_CC,
    Batch14008_gen_CC,
    Batch14009_gen_CC
  )
write.csv(Batch14_CC, file = "E:/Models_Outputs/data_analysis/CC/Batch14_CC.csv", row.names = F)


Batch15001_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch15001/Outputs")
Batch15002_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch15002/Outputs")
Batch15003_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch15003/Outputs")
Batch15004_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch15004/Outputs")
Batch15005_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch15005/Outputs")
Batch15006_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch15006/Outputs")
Batch15007_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch15007/Outputs")
Batch15008_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch15008/Outputs")
Batch15009_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch15009/Outputs")
Batch15_CC <-
  rbind(
    Batch15001_gen_CC,
    Batch15002_gen_CC,
    Batch15003_gen_CC,
    Batch15004_gen_CC,
    Batch15005_gen_CC,
    Batch15006_gen_CC,
    Batch15007_gen_CC,
    Batch15008_gen_CC,
    Batch15009_gen_CC
  )
write.csv(Batch15_CC, file = "E:/Models_Outputs/data_analysis/CC/Batch15_CC.csv", row.names = F)



Batch16001_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch16001/Outputs")
Batch16002_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch16002/Outputs")
Batch16003_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch16003/Outputs")
Batch16004_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch16004/Outputs")
Batch16005_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch16005/Outputs")
Batch16006_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch16006/Outputs")
Batch16007_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch16007/Outputs")
Batch16008_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch16008/Outputs")
Batch16009_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch16009/Outputs")
Batch16_CC <-
  rbind(
    Batch16001_gen_CC,
    Batch16002_gen_CC,
    Batch16003_gen_CC,
    Batch16004_gen_CC,
    Batch16005_gen_CC,
    Batch16006_gen_CC,
    Batch16007_gen_CC,
    Batch16008_gen_CC,
    Batch16009_gen_CC
  )
write.csv(Batch16_CC, file = "E:/Models_Outputs/data_analysis/CC/Batch16_CC.csv", row.names = F)



Batch17001_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch17001/Outputs")
Batch17002_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch17002/Outputs")
Batch17003_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch17003/Outputs")
Batch17004_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch17004/Outputs")
Batch17005_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch17005/Outputs")
Batch17006_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch17006/Outputs")
Batch17007_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch17007/Outputs")
Batch17008_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch17008/Outputs")
Batch17009_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch17009/Outputs")
Batch17_CC <-
  rbind(
    Batch17001_gen_CC,
    Batch17002_gen_CC,
    Batch17003_gen_CC,
    Batch17004_gen_CC,
    Batch17005_gen_CC,
    Batch17006_gen_CC,
    Batch17007_gen_CC,
    Batch17008_gen_CC,
    Batch17009_gen_CC
  )
write.csv(Batch17_CC, file = "E:/Models_Outputs/data_analysis/CC/Batch17_CC.csv", row.names = F)



Batch18001_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch18001/Outputs")
Batch18002_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch18002/Outputs")
Batch18003_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch18003/Outputs")
Batch18004_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch18004/Outputs")
Batch18005_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch18005/Outputs")
Batch18006_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch18006/Outputs")
Batch18007_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch18007/Outputs")
Batch18008_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch18008/Outputs")
Batch18009_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch18009/Outputs")
Batch18_CC <-
  rbind(
    Batch18001_gen_CC,
    Batch18002_gen_CC,
    Batch18003_gen_CC,
    Batch18004_gen_CC,
    Batch18005_gen_CC,
    Batch18006_gen_CC,
    Batch18007_gen_CC,
    Batch18008_gen_CC,
    Batch18009_gen_CC
  )
write.csv(Batch18_CC, file = "E:/Models_Outputs/data_analysis/CC/Batch18_CC.csv", row.names = F)


Batch19001_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch19001/Outputs")
Batch19002_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch19002/Outputs")
Batch19003_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch19003/Outputs")
Batch19004_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch19004/Outputs")
Batch19005_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch19005/Outputs")
Batch19006_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch19006/Outputs")
Batch19007_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch19007/Outputs")
Batch19008_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch19008/Outputs")
Batch19009_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch19009/Outputs")
Batch19_CC <-
  rbind(
    Batch19001_gen_CC,
    Batch19002_gen_CC,
    Batch19003_gen_CC,
    Batch19004_gen_CC,
    Batch19005_gen_CC,
    Batch19006_gen_CC,
    Batch19007_gen_CC,
    Batch19008_gen_CC,
    Batch19009_gen_CC
  )
write.csv(Batch19_CC, file = "E:/Models_Outputs/data_analysis/CC/Batch19_CC.csv", row.names = F)


Batch20001_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land20/Batch20001/Outputs")
Batch20002_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land20/Batch20002/Outputs")
Batch20003_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land20/Batch20003/Outputs")
Batch20004_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land20/Batch20004/Outputs")
Batch20005_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land20/Batch20005/Outputs")
Batch20006_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land20/Batch20006/Outputs")
Batch20007_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land20/Batch20007/Outputs")
Batch20008_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land20/Batch20008/Outputs")
Batch20009_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land20/Batch20009/Outputs")
Batch20_CC <-
  rbind(
    Batch20001_gen_CC,
    Batch20002_gen_CC,
    Batch20003_gen_CC,
    Batch20004_gen_CC,
    Batch20005_gen_CC,
    Batch20006_gen_CC,
    Batch20007_gen_CC,
    Batch20008_gen_CC,
    Batch20009_gen_CC
  )
write.csv(Batch20_CC, file = "E:/Models_Outputs/data_analysis/CC/Batch20_CC.csv", row.names = F)


Batch21001_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land21/Batch21001/Outputs")
Batch21002_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land21/Batch21002/Outputs")
Batch21003_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land21/Batch21003/Outputs")
Batch21004_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land21/Batch21004/Outputs")
Batch21005_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land21/Batch21005/Outputs")
Batch21006_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land21/Batch21006/Outputs")
Batch21007_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land21/Batch21007/Outputs")
Batch21008_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land21/Batch21008/Outputs")
Batch21009_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land21/Batch21009/Outputs")
Batch21_CC <-
  rbind(
    Batch21001_gen_CC,
    Batch21002_gen_CC,
    Batch21003_gen_CC,
    Batch21004_gen_CC,
    Batch21005_gen_CC,
    Batch21006_gen_CC,
    Batch21007_gen_CC,
    Batch21008_gen_CC,
    Batch21009_gen_CC
  )
write.csv(Batch21_CC, file = "E:/Models_Outputs/data_analysis/CC/Batch21_CC.csv", row.names = F)



Batch22001_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land22/Batch22001/Outputs")
Batch22002_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land22/Batch22002/Outputs")
Batch22003_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land22/Batch22003/Outputs")
Batch22004_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land22/Batch22004/Outputs")
Batch22005_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land22/Batch22005/Outputs")
Batch22006_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land22/Batch22006/Outputs")
Batch22007_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land22/Batch22007/Outputs")
Batch22008_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land22/Batch22008/Outputs")
Batch22009_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land22/Batch22009/Outputs")
Batch22_CC <-
  rbind(
    Batch22001_gen_CC,
    Batch22002_gen_CC,
    Batch22003_gen_CC,
    Batch22004_gen_CC,
    Batch22005_gen_CC,
    Batch22006_gen_CC,
    Batch22007_gen_CC,
    Batch22008_gen_CC,
    Batch22009_gen_CC
  )
write.csv(Batch22_CC, file = "E:/Models_Outputs/data_analysis/CC/Batch22_CC.csv", row.names = F)


Batch23001_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land23/Batch23001/Outputs")
Batch23002_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land23/Batch23002/Outputs")
Batch23003_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land23/Batch23003/Outputs")
Batch23004_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land23/Batch23004/Outputs")
Batch23005_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land23/Batch23005/Outputs")
Batch23006_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land23/Batch23006/Outputs")
Batch23007_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land23/Batch23007/Outputs")
Batch23008_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land23/Batch23008/Outputs")
Batch23009_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land23/Batch23009/Outputs")
Batch23_CC <-
  rbind(
    Batch23001_gen_CC,
    Batch23002_gen_CC,
    Batch23003_gen_CC,
    Batch23004_gen_CC,
    Batch23005_gen_CC,
    Batch23006_gen_CC,
    Batch23007_gen_CC,
    Batch23008_gen_CC,
    Batch23009_gen_CC
  )
write.csv(Batch23_CC, file = "E:/Models_Outputs/data_analysis/CC/Batch23_CC.csv", row.names = F)



Batch24001_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land24/Batch24001/Outputs")
Batch24002_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land24/Batch24002/Outputs")
Batch24003_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land24/Batch24003/Outputs")
Batch24004_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land24/Batch24004/Outputs")
Batch24005_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land24/Batch24005/Outputs")
Batch24006_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land24/Batch24006/Outputs")
Batch24007_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land24/Batch24007/Outputs")
Batch24008_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land24/Batch24008/Outputs")
Batch24009_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land24/Batch24009/Outputs")
Batch24_CC <-
  rbind(
    Batch24001_gen_CC,
    Batch24002_gen_CC,
    Batch24003_gen_CC,
    Batch24004_gen_CC,
    Batch24005_gen_CC,
    Batch24006_gen_CC,
    Batch24007_gen_CC,
    Batch24008_gen_CC,
    Batch24009_gen_CC
  )
write.csv(Batch24_CC, file = "E:/Models_Outputs/data_analysis/CC/Batch24_CC.csv", row.names = F)



#### Josts D ####



Batch1001_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch1001/Outputs")
Batch1002_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch1002/Outputs")
Batch1003_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch1003/Outputs")
Batch1004_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch1004/Outputs")
Batch1005_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch1005/Outputs")
Batch1006_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch1006/Outputs")
Batch1007_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch1007/Outputs")
Batch1008_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch1008/Outputs")
Batch1009_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch1009/Outputs")
Batch1_CC <-
  rbind(
    Batch1001_gen_CC,
    Batch1002_gen_CC,
    Batch1003_gen_CC,
    Batch1004_gen_CC,
    Batch1005_gen_CC,
    Batch1006_gen_CC,
    Batch1007_gen_CC,
    Batch1008_gen_CC,
    Batch1009_gen_CC
  )
write.csv(Batch1_CC, file = "E:/Models_Outputs/data_analysis/Josts/Batch1_CC.csv", row.names = F)





Batch2001_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch2001/Outputs")
Batch2002_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch2002/Outputs")
Batch2003_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch2003/Outputs")
Batch2004_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch2004/Outputs")
Batch2005_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch2005/Outputs")
Batch2006_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch2006/Outputs")
Batch2007_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch2007/Outputs")
Batch2008_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch2008/Outputs")
Batch2009_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch2009/Outputs")
Batch2_CC <-
  rbind(
    Batch2001_gen_CC,
    Batch2002_gen_CC,
    Batch2003_gen_CC,
    Batch2004_gen_CC,
    Batch2005_gen_CC,
    Batch2006_gen_CC,
    Batch2007_gen_CC,
    Batch2008_gen_CC,
    Batch2009_gen_CC
  )
write.csv(Batch2_CC, file = "E:/Models_Outputs/data_analysis/Josts/Batch2_CC.csv", row.names = F)



Batch3001_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch3001/Outputs")
Batch3002_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch3002/Outputs")
Batch3003_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch3003/Outputs")
Batch3004_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch3004/Outputs")
Batch3005_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch3005/Outputs")
Batch3006_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch3006/Outputs")
Batch3007_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch3007/Outputs")
Batch3008_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch3008/Outputs")
Batch3009_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch3009/Outputs")
Batch3_CC <-
  rbind(
    Batch3001_gen_CC,
    Batch3002_gen_CC,
    Batch3003_gen_CC,
    Batch3004_gen_CC,
    Batch3005_gen_CC,
    Batch3006_gen_CC,
    Batch3007_gen_CC,
    Batch3008_gen_CC,
    Batch3009_gen_CC
  )
write.csv(Batch3_CC, file = "E:/Models_Outputs/data_analysis/Josts/Batch3_CC.csv", row.names = F)


Batch4001_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch4001/Outputs")
Batch4002_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch4002/Outputs")
Batch4003_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch4003/Outputs")
Batch4004_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch4004/Outputs")
Batch4005_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch4005/Outputs")
Batch4006_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch4006/Outputs")
Batch4007_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch4007/Outputs")
Batch4008_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch4008/Outputs")
Batch4009_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch4009/Outputs")
Batch4_CC <-
  rbind(
    Batch4001_gen_CC,
    Batch4002_gen_CC,
    Batch4003_gen_CC,
    Batch4004_gen_CC,
    Batch4005_gen_CC,
    Batch4006_gen_CC,
    Batch4007_gen_CC,
    Batch4008_gen_CC,
    Batch4009_gen_CC
  )
write.csv(Batch4_CC, file = "E:/Models_Outputs/data_analysis/Josts/Batch4_CC.csv", row.names = F)


Batch5001_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch5001/Outputs")
Batch5002_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch5002/Outputs")
Batch5003_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch5003/Outputs")
Batch5004_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch5004/Outputs")
Batch5005_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch5005/Outputs")
Batch5006_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch5006/Outputs")
Batch5007_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch5007/Outputs")
Batch5008_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch5008/Outputs")
Batch5009_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch5009/Outputs")
Batch5_CC <-
  rbind(
    Batch5001_gen_CC,
    Batch5002_gen_CC,
    Batch5003_gen_CC,
    Batch5004_gen_CC,
    Batch5005_gen_CC,
    Batch5006_gen_CC,
    Batch5007_gen_CC,
    Batch5008_gen_CC,
    Batch5009_gen_CC
  )
write.csv(Batch5_CC, file = "E:/Models_Outputs/data_analysis/Josts/Batch5_CC.csv", row.names = F)


Batch6001_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch6001/Outputs")
Batch6002_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch6002/Outputs")
Batch6003_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch6003/Outputs")
Batch6004_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch6004/Outputs")
Batch6005_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch6005/Outputs")
Batch6006_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch6006/Outputs")
Batch6007_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch6007/Outputs")
Batch6008_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch6008/Outputs")
Batch6009_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch6009/Outputs")
Batch6_CC <-
  rbind(
    Batch6001_gen_CC,
    Batch6002_gen_CC,
    Batch6003_gen_CC,
    Batch6004_gen_CC,
    Batch6005_gen_CC,
    Batch6006_gen_CC,
    Batch6007_gen_CC,
    Batch6008_gen_CC,
    Batch6009_gen_CC
  )
write.csv(Batch6_CC, file = "E:/Models_Outputs/data_analysis/Josts/Batch6_CC.csv", row.names = F)


Batch7001_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch7001/Outputs")
Batch7002_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch7002/Outputs")
Batch7003_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch7003/Outputs")
Batch7004_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch7004/Outputs")
Batch7005_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch7005/Outputs")
Batch7006_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch7006/Outputs")
Batch7007_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch7007/Outputs")
Batch7008_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch7008/Outputs")
Batch7009_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch7009/Outputs")
Batch7_CC <-
  rbind(
    Batch7001_gen_CC,
    Batch7002_gen_CC,
    Batch7003_gen_CC,
    Batch7004_gen_CC,
    Batch7005_gen_CC,
    Batch7006_gen_CC,
    Batch7007_gen_CC,
    Batch7008_gen_CC,
    Batch7009_gen_CC
  )
write.csv(Batch7_CC, file = "E:/Models_Outputs/data_analysis/Josts/Batch7_CC.csv", row.names = F)


Batch8001_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch8001/Outputs")
Batch8002_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch8002/Outputs")
Batch8003_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch8003/Outputs")
Batch8004_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch8004/Outputs")
Batch8005_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch8005/Outputs")
Batch8006_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch8006/Outputs")
Batch8007_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch8007/Outputs")
Batch8008_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch8008/Outputs")
Batch8009_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch8009/Outputs")
Batch8_CC <-
  rbind(
    Batch8001_gen_CC,
    Batch8002_gen_CC,
    Batch8003_gen_CC,
    Batch8004_gen_CC,
    Batch8005_gen_CC,
    Batch8006_gen_CC,
    Batch8007_gen_CC,
    Batch8008_gen_CC,
    Batch8009_gen_CC
  )
write.csv(Batch8_CC, file = "E:/Models_Outputs/data_analysis/Josts/Batch8_CC.csv", row.names = F)



Batch9001_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch9001/Outputs")
Batch9002_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch9002/Outputs")
Batch9003_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch9003/Outputs")
Batch9004_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch9004/Outputs")
Batch9005_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch9005/Outputs")
Batch9006_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch9006/Outputs")
Batch9007_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch9007/Outputs")
Batch9008_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch9008/Outputs")
Batch9009_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch9009/Outputs")
Batch9_CC <-
  rbind(
    Batch9001_gen_CC,
    Batch9002_gen_CC,
    Batch9003_gen_CC,
    Batch9004_gen_CC,
    Batch9005_gen_CC,
    Batch9006_gen_CC,
    Batch9007_gen_CC,
    Batch9008_gen_CC,
    Batch9009_gen_CC
  )
write.csv(Batch9_CC, file = "E:/Models_Outputs/data_analysis/Josts/Batch9_CC.csv", row.names = F)



Batch10001_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch10001/Outputs")
Batch10002_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch10002/Outputs")
Batch10003_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch10003/Outputs")
Batch10004_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch10004/Outputs")
Batch10005_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch10005/Outputs")
Batch10006_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch10006/Outputs")
Batch10007_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch10007/Outputs")
Batch10008_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch10008/Outputs")
Batch10009_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch10009/Outputs")
Batch10_CC <-
  rbind(
    Batch10001_gen_CC,
    Batch10002_gen_CC,
    Batch10003_gen_CC,
    Batch10004_gen_CC,
    Batch10005_gen_CC,
    Batch10006_gen_CC,
    Batch10007_gen_CC,
    Batch10008_gen_CC,
    Batch10009_gen_CC
  )
write.csv(Batch10_CC, file = "E:/Models_Outputs/data_analysis/Josts/Batch10_CC.csv", row.names = F)



Batch11001_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch11001/Outputs")
Batch11002_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch11002/Outputs")
Batch11003_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch11003/Outputs")
Batch11004_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch11004/Outputs")
Batch11005_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch11005/Outputs")
Batch11006_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch11006/Outputs")
Batch11007_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch11007/Outputs")
Batch11008_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch11008/Outputs")
Batch11009_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch11009/Outputs")
Batch11_CC <-
  rbind(
    Batch11001_gen_CC,
    Batch11002_gen_CC,
    Batch11003_gen_CC,
    Batch11004_gen_CC,
    Batch11005_gen_CC,
    Batch11006_gen_CC,
    Batch11007_gen_CC,
    Batch11008_gen_CC,
    Batch11009_gen_CC
  )
write.csv(Batch11_CC, file = "E:/Models_Outputs/data_analysis/Josts/Batch11_CC.csv", row.names = F)



Batch12001_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch12001/Outputs")
Batch12002_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch12002/Outputs")
Batch12003_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch12003/Outputs")
Batch12004_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch12004/Outputs")
Batch12005_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch12005/Outputs")
Batch12006_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch12006/Outputs")
Batch12007_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch12007/Outputs")
Batch12008_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch12008/Outputs")
Batch12009_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch12009/Outputs")
Batch12_CC <-
  rbind(
    Batch12001_gen_CC,
    Batch12002_gen_CC,
    Batch12003_gen_CC,
    Batch12004_gen_CC,
    Batch12005_gen_CC,
    Batch12006_gen_CC,
    Batch12007_gen_CC,
    Batch12008_gen_CC,
    Batch12009_gen_CC
  )
write.csv(Batch12_CC, file = "E:/Models_Outputs/data_analysis/Josts/Batch12_CC.csv", row.names = F)



Batch13001_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch13001/Outputs")
Batch13002_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch13002/Outputs")
Batch13003_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch13003/Outputs")
Batch13004_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch13004/Outputs")
Batch13005_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch13005/Outputs")
Batch13006_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch13006/Outputs")
Batch13007_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch13007/Outputs")
Batch13008_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch13008/Outputs")
Batch13009_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch13009/Outputs")
Batch13_CC <-
  rbind(
    Batch13001_gen_CC,
    Batch13002_gen_CC,
    Batch13003_gen_CC,
    Batch13004_gen_CC,
    Batch13005_gen_CC,
    Batch13006_gen_CC,
    Batch13007_gen_CC,
    Batch13008_gen_CC,
    Batch13009_gen_CC
  )
write.csv(Batch13_CC, file = "E:/Models_Outputs/data_analysis/Josts/Batch13_CC.csv", row.names = F)



Batch14001_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch14001/Outputs")
Batch14002_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch14002/Outputs")
Batch14003_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch14003/Outputs")
Batch14004_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch14004/Outputs")
Batch14005_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch14005/Outputs")
Batch14006_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch14006/Outputs")
Batch14007_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch14007/Outputs")
Batch14008_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch14008/Outputs")
Batch14009_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch14009/Outputs")
Batch14_CC <-
  rbind(
    Batch14001_gen_CC,
    Batch14002_gen_CC,
    Batch14003_gen_CC,
    Batch14004_gen_CC,
    Batch14005_gen_CC,
    Batch14006_gen_CC,
    Batch14007_gen_CC,
    Batch14008_gen_CC,
    Batch14009_gen_CC
  )
write.csv(Batch14_CC, file = "E:/Models_Outputs/data_analysis/Josts/Batch14_CC.csv", row.names = F)


Batch15001_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch15001/Outputs")
Batch15002_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch15002/Outputs")
Batch15003_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch15003/Outputs")
Batch15004_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch15004/Outputs")
Batch15005_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch15005/Outputs")
Batch15006_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch15006/Outputs")
Batch15007_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch15007/Outputs")
Batch15008_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch15008/Outputs")
Batch15009_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch15009/Outputs")
Batch15_CC <-
  rbind(
    Batch15001_gen_CC,
    Batch15002_gen_CC,
    Batch15003_gen_CC,
    Batch15004_gen_CC,
    Batch15005_gen_CC,
    Batch15006_gen_CC,
    Batch15007_gen_CC,
    Batch15008_gen_CC,
    Batch15009_gen_CC
  )
write.csv(Batch15_CC, file = "E:/Models_Outputs/data_analysis/Josts/Batch15_CC.csv", row.names = F)



Batch16001_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch16001/Outputs")
Batch16002_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch16002/Outputs")
Batch16003_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch16003/Outputs")
Batch16004_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch16004/Outputs")
Batch16005_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch16005/Outputs")
Batch16006_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch16006/Outputs")
Batch16007_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch16007/Outputs")
Batch16008_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch16008/Outputs")
Batch16009_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch16009/Outputs")
Batch16_CC <-
  rbind(
    Batch16001_gen_CC,
    Batch16002_gen_CC,
    Batch16003_gen_CC,
    Batch16004_gen_CC,
    Batch16005_gen_CC,
    Batch16006_gen_CC,
    Batch16007_gen_CC,
    Batch16008_gen_CC,
    Batch16009_gen_CC
  )
write.csv(Batch16_CC, file = "E:/Models_Outputs/data_analysis/Josts/Batch16_CC.csv", row.names = F)



Batch17001_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch17001/Outputs")
Batch17002_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch17002/Outputs")
Batch17003_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch17003/Outputs")
Batch17004_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch17004/Outputs")
Batch17005_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch17005/Outputs")
Batch17006_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch17006/Outputs")
Batch17007_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch17007/Outputs")
Batch17008_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch17008/Outputs")
Batch17009_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch17009/Outputs")
Batch17_CC <-
  rbind(
    Batch17001_gen_CC,
    Batch17002_gen_CC,
    Batch17003_gen_CC,
    Batch17004_gen_CC,
    Batch17005_gen_CC,
    Batch17006_gen_CC,
    Batch17007_gen_CC,
    Batch17008_gen_CC,
    Batch17009_gen_CC
  )
write.csv(Batch17_CC, file = "E:/Models_Outputs/data_analysis/Josts/Batch17_CC.csv", row.names = F)



Batch18001_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch18001/Outputs")
Batch18002_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch18002/Outputs")
Batch18003_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch18003/Outputs")
Batch18004_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch18004/Outputs")
Batch18005_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch18005/Outputs")
Batch18006_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch18006/Outputs")
Batch18007_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch18007/Outputs")
Batch18008_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch18008/Outputs")
Batch18009_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch18009/Outputs")
Batch18_CC <-
  rbind(
    Batch18001_gen_CC,
    Batch18002_gen_CC,
    Batch18003_gen_CC,
    Batch18004_gen_CC,
    Batch18005_gen_CC,
    Batch18006_gen_CC,
    Batch18007_gen_CC,
    Batch18008_gen_CC,
    Batch18009_gen_CC
  )
write.csv(Batch18_CC, file = "E:/Models_Outputs/data_analysis/Josts/Batch18_CC.csv", row.names = F)


Batch19001_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch19001/Outputs")
Batch19002_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch19002/Outputs")
Batch19003_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch19003/Outputs")
Batch19004_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch19004/Outputs")
Batch19005_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch19005/Outputs")
Batch19006_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch19006/Outputs")
Batch19007_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch19007/Outputs")
Batch19008_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch19008/Outputs")
Batch19009_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Batch19009/Outputs")
Batch19_CC <-
  rbind(
    Batch19001_gen_CC,
    Batch19002_gen_CC,
    Batch19003_gen_CC,
    Batch19004_gen_CC,
    Batch19005_gen_CC,
    Batch19006_gen_CC,
    Batch19007_gen_CC,
    Batch19008_gen_CC,
    Batch19009_gen_CC
  )
write.csv(Batch19_CC, file = "E:/Models_Outputs/data_analysis/Josts/Batch19_CC.csv", row.names = F)


Batch20001_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land20/Batch20001/Outputs")
Batch20002_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land20/Batch20002/Outputs")
Batch20003_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land20/Batch20003/Outputs")
Batch20004_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land20/Batch20004/Outputs")
Batch20005_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land20/Batch20005/Outputs")
Batch20006_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land20/Batch20006/Outputs")
Batch20007_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land20/Batch20007/Outputs")
Batch20008_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land20/Batch20008/Outputs")
Batch20009_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land20/Batch20009/Outputs")
Batch20_CC <-
  rbind(
    Batch20001_gen_CC,
    Batch20002_gen_CC,
    Batch20003_gen_CC,
    Batch20004_gen_CC,
    Batch20005_gen_CC,
    Batch20006_gen_CC,
    Batch20007_gen_CC,
    Batch20008_gen_CC,
    Batch20009_gen_CC
  )
write.csv(Batch20_CC, file = "E:/Models_Outputs/data_analysis/Josts/Batch20_CC.csv", row.names = F)


Batch21001_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land21/Batch21001/Outputs")
Batch21002_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land21/Batch21002/Outputs")
Batch21003_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land21/Batch21003/Outputs")
Batch21004_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land21/Batch21004/Outputs")
Batch21005_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land21/Batch21005/Outputs")
Batch21006_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land21/Batch21006/Outputs")
Batch21007_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land21/Batch21007/Outputs")
Batch21008_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land21/Batch21008/Outputs")
Batch21009_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land21/Batch21009/Outputs")
Batch21_CC <-
  rbind(
    Batch21001_gen_CC,
    Batch21002_gen_CC,
    Batch21003_gen_CC,
    Batch21004_gen_CC,
    Batch21005_gen_CC,
    Batch21006_gen_CC,
    Batch21007_gen_CC,
    Batch21008_gen_CC,
    Batch21009_gen_CC
  )
write.csv(Batch21_CC, file = "E:/Models_Outputs/data_analysis/Josts/Batch21_CC.csv", row.names = F)



Batch22001_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land22/Batch22001/Outputs")
Batch22002_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land22/Batch22002/Outputs")
Batch22003_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land22/Batch22003/Outputs")
Batch22004_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land22/Batch22004/Outputs")
Batch22005_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land22/Batch22005/Outputs")
Batch22006_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land22/Batch22006/Outputs")
Batch22007_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land22/Batch22007/Outputs")
Batch22008_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land22/Batch22008/Outputs")
Batch22009_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land22/Batch22009/Outputs")
Batch22_CC <-
  rbind(
    Batch22001_gen_CC,
    Batch22002_gen_CC,
    Batch22003_gen_CC,
    Batch22004_gen_CC,
    Batch22005_gen_CC,
    Batch22006_gen_CC,
    Batch22007_gen_CC,
    Batch22008_gen_CC,
    Batch22009_gen_CC
  )
write.csv(Batch22_CC, file = "E:/Models_Outputs/data_analysis/Josts/Batch22_CC.csv", row.names = F)


Batch23001_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land23/Batch23001/Outputs")
Batch23002_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land23/Batch23002/Outputs")
Batch23003_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land23/Batch23003/Outputs")
Batch23004_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land23/Batch23004/Outputs")
Batch23005_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land23/Batch23005/Outputs")
Batch23006_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land23/Batch23006/Outputs")
Batch23007_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land23/Batch23007/Outputs")
Batch23008_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land23/Batch23008/Outputs")
Batch23009_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land23/Batch23009/Outputs")
Batch23_CC <-
  rbind(
    Batch23001_gen_CC,
    Batch23002_gen_CC,
    Batch23003_gen_CC,
    Batch23004_gen_CC,
    Batch23005_gen_CC,
    Batch23006_gen_CC,
    Batch23007_gen_CC,
    Batch23008_gen_CC,
    Batch23009_gen_CC
  )
write.csv(Batch23_CC, file = "E:/Models_Outputs/data_analysis/Josts/Batch23_CC.csv", row.names = F)



Batch24001_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land24/Batch24001/Outputs")
Batch24002_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land24/Batch24002/Outputs")
Batch24003_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land24/Batch24003/Outputs")
Batch24004_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land24/Batch24004/Outputs")
Batch24005_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land24/Batch24005/Outputs")
Batch24006_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land24/Batch24006/Outputs")
Batch24007_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land24/Batch24007/Outputs")
Batch24008_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land24/Batch24008/Outputs")
Batch24009_gen_CC <-
  Isolation_by_distance_jost(min_pop_size = 5, mypath = "E:/Models_Outputs/Land24/Batch24009/Outputs")
Batch24_CC <-
  rbind(
    Batch24001_gen_CC,
    Batch24002_gen_CC,
    Batch24003_gen_CC,
    Batch24004_gen_CC,
    Batch24005_gen_CC,
    Batch24006_gen_CC,
    Batch24007_gen_CC,
    Batch24008_gen_CC,
    Batch24009_gen_CC
  )
write.csv(Batch24_CC, file = "E:/Models_Outputs/data_analysis/Josts/Batch24_CC.csv", row.names = F)




# saving ####

write.csv(Batch1_CC, file = "E:/Models_Outputs/data_analysis/Josts/Batch1_CC.csv", row.names = F)
write.csv(Batch2_CC, file = "E:/Models_Outputs/data_analysis/Josts/Batch2_CC.csv", row.names = F)
write.csv(Batch3_CC, file = "E:/Models_Outputs/data_analysis/Josts/Batch3_CC.csv", row.names = F)
write.csv(Batch4_CC, file = "E:/Models_Outputs/data_analysis/Josts/Batch4_CC.csv", row.names = F)
write.csv(Batch5_CC, file = "E:/Models_Outputs/data_analysis/Josts/Batch5_CC.csv", row.names = F)
write.csv(Batch6_CC, file = "E:/Models_Outputs/data_analysis/Josts/Batch6_CC.csv", row.names = F)
write.csv(Batch7_CC, file = "E:/Models_Outputs/data_analysis/Josts/Batch7_CC.csv", row.names = F)
write.csv(Batch8_CC, file = "E:/Models_Outputs/data_analysis/Josts/Batch8_CC.csv", row.names = F)
write.csv(Batch9_CC, file = "E:/Models_Outputs/data_analysis/Josts/Batch9_CC.csv", row.names = F)
write.csv(Batch10_CC, file = "E:/Models_Outputs/data_analysis/Josts/Batch10_CC.csv", row.names = F)
write.csv(Batch11_CC, file = "E:/Models_Outputs/data_analysis/Josts/Batch11_CC.csv", row.names = F)
write.csv(Batch12_CC, file = "E:/Models_Outputs/data_analysis/Josts/Batch12_CC.csv", row.names = F)
write.csv(Batch13_CC, file = "E:/Models_Outputs/data_analysis/Josts/Batch13_CC.csv", row.names = F)
write.csv(Batch14_CC, file = "E:/Models_Outputs/data_analysis/Josts/Batch14_CC.csv", row.names = F)
write.csv(Batch15_CC, file = "E:/Models_Outputs/data_analysis/Josts/Batch15_CC.csv", row.names = F)




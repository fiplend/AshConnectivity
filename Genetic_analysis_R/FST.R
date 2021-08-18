## reformating allele data from R

setwd("E:/Subset_data")
library(adegenet)
library(hierfstat)

#Test 1 ####
Test1_samples <- read.table("Batch13_Sim1_Land1_GenSamples.txt", header=TRUE)
T1_subset <- subset(Test1_samples, Test1_samples$Year == 90)
T1_subset <- subset(T1_subset, T1_subset$Rep == 1)
location_1 <- unname(unlist(T1_subset["PatchID"]))

T1_formatted <- T1_subset
T1_formatted <- T1_formatted[,6:35]+150
T1_formatted$loc0 <- paste(T1_formatted$Chr0Loc0Allele0, "/", T1_formatted$Chr0Loc0Allele1)
T1_formatted$loc1 <- paste(T1_formatted$Chr0Loc1Allele0, "/", T1_formatted$Chr0Loc1Allele1)
T1_formatted$loc2 <- paste(T1_formatted$Chr0Loc2Allele0, "/", T1_formatted$Chr0Loc2Allele1)
T1_formatted$loc3 <- paste(T1_formatted$Chr0Loc3Allele0, "/", T1_formatted$Chr0Loc3Allele1)
T1_formatted$loc4 <- paste(T1_formatted$Chr0Loc4Allele0, "/", T1_formatted$Chr0Loc4Allele1)
T1_formatted$loc5 <- paste(T1_formatted$Chr0Loc5Allele0, "/", T1_formatted$Chr0Loc5Allele1)
T1_formatted$loc6 <- paste(T1_formatted$Chr0Loc6Allele0, "/", T1_formatted$Chr0Loc6Allele1)
T1_formatted$loc7 <- paste(T1_formatted$Chr0Loc7Allele0, "/", T1_formatted$Chr0Loc7Allele1)
T1_formatted$loc8 <- paste(T1_formatted$Chr0Loc8Allele0, "/", T1_formatted$Chr0Loc8Allele1)
T1_formatted$loc9 <- paste(T1_formatted$Chr0Loc9Allele0, "/", T1_formatted$Chr0Loc9Allele1)
T1_formatted$loc10 <- paste(T1_formatted$Chr0Loc10Allele0, "/", T1_formatted$Chr0Loc10Allele1)
T1_formatted$loc11 <- paste(T1_formatted$Chr0Loc11Allele0, "/", T1_formatted$Chr0Loc11Allele1)
T1_formatted$loc12 <- paste(T1_formatted$Chr0Loc12Allele0, "/", T1_formatted$Chr0Loc12Allele1)
T1_formatted$loc13 <- paste(T1_formatted$Chr0Loc13Allele0, "/", T1_formatted$Chr0Loc13Allele1)
T1_formatted$loc14 <- paste(T1_formatted$Chr0Loc14Allele0, "/", T1_formatted$Chr0Loc14Allele1)
T1_formatted <- T1_formatted[,31:45]
T_genind_1 <- df2genind(T1_formatted, pop = location_1, ncode = 3, sep = "/")

alleles(T_genind_1)

T1_raw_gen <- genet.dist(T_genind_1, method = "WC84")

##### extracting genetic data

library(data.table)

gen_raw1 <- as.matrix(genet.dist(T_genind_1, method = "WC84"))
gen_raw1 <- as.data.frame(gen_raw1)
setDT(gen_raw1, keep.rownames = TRUE)[]

# prepare data
# csv data has useful information in the row.names
# this has been loaded as the first column
# rename as patch_start

colnames(gen_raw1)[1] <- "patch_start"


# convert patch info to factor
gen_raw1$patch_start <- as.factor(gen_raw1$patch_start)


# data is a triangular matrix
# let's just keep the lower half & convert to lower triangular matrix
# note: we need the diagonal
# (could alternatively use upper triangular matrix)


gen_ind1 <- lower.tri(gen_raw1, diag = TRUE)

# select values of interest and replace rest with NAs

gen_tri1 <- gen_raw1
gen_tri1[gen_ind1 == FALSE] <- NA


#melt data

gen_long1 <-
  melt(
    data = gen_tri1,
    id.vars = "patch_start",
    variable.name = "patch_end",
    value.name = "gen_distance",
    na.rm = TRUE
  )


gen_long1$Tree_disease <- 0
gen_long1$Management <- 0
gen_long1$Replicate <- 1

#Test 2 ####
Test2_samples <- read.table("Batch13_Sim1_Land2_GenSamples.txt", header=TRUE)
T2_subset <- subset(Test2_samples, Test2_samples$Year == 90)
T2_subset <- subset(T2_subset, T2_subset$Rep == 1)
location_2 <- unname(unlist(T2_subset["PatchID"]))

T2_formatted <- T2_subset
T2_formatted <- T2_formatted[,6:35]+150
T2_formatted$loc0 <- paste(T2_formatted$Chr0Loc0Allele0, "/", T2_formatted$Chr0Loc0Allele1)
T2_formatted$loc1 <- paste(T2_formatted$Chr0Loc1Allele0, "/", T2_formatted$Chr0Loc1Allele1)
T2_formatted$loc2 <- paste(T2_formatted$Chr0Loc2Allele0, "/", T2_formatted$Chr0Loc2Allele1)
T2_formatted$loc3 <- paste(T2_formatted$Chr0Loc3Allele0, "/", T2_formatted$Chr0Loc3Allele1)
T2_formatted$loc4 <- paste(T2_formatted$Chr0Loc4Allele0, "/", T2_formatted$Chr0Loc4Allele1)
T2_formatted$loc5 <- paste(T2_formatted$Chr0Loc5Allele0, "/", T2_formatted$Chr0Loc5Allele1)
T2_formatted$loc6 <- paste(T2_formatted$Chr0Loc6Allele0, "/", T2_formatted$Chr0Loc6Allele1)
T2_formatted$loc7 <- paste(T2_formatted$Chr0Loc7Allele0, "/", T2_formatted$Chr0Loc7Allele1)
T2_formatted$loc8 <- paste(T2_formatted$Chr0Loc8Allele0, "/", T2_formatted$Chr0Loc8Allele1)
T2_formatted$loc9 <- paste(T2_formatted$Chr0Loc9Allele0, "/", T2_formatted$Chr0Loc9Allele1)
T2_formatted$loc10 <- paste(T2_formatted$Chr0Loc10Allele0, "/", T2_formatted$Chr0Loc10Allele1)
T2_formatted$loc11 <- paste(T2_formatted$Chr0Loc11Allele0, "/", T2_formatted$Chr0Loc11Allele1)
T2_formatted$loc12 <- paste(T2_formatted$Chr0Loc12Allele0, "/", T2_formatted$Chr0Loc12Allele1)
T2_formatted$loc13 <- paste(T2_formatted$Chr0Loc13Allele0, "/", T2_formatted$Chr0Loc13Allele1)
T2_formatted$loc14 <- paste(T2_formatted$Chr0Loc14Allele0, "/", T2_formatted$Chr0Loc14Allele1)
T2_formatted <- T2_formatted[,31:45]
T_genind_2 <- df2genind(T2_formatted, pop = location_2, ncode = 3, sep = "/")

T2_raw_gen <- genet.dist(T_genind_2, method = "WC84")

gen_raw2 <- as.matrix(genet.dist(T_genind_2, method = "WC84"))
gen_raw2 <- as.data.frame(gen_raw2)
setDT(gen_raw2, keep.rownames = TRUE)[]

# prepare data
# csv data has useful information in the row.names
# this has been loaded as the first column
# rename as patch_start

colnames(gen_raw2)[1] <- "patch_start"


# convert patch info to factor
gen_raw2$patch_start <- as.factor(gen_raw2$patch_start)


# data is a triangular matrix
# let's just keep the lower half & convert to lower triangular matrix
# note: we need the diagonal
# (could alternatively use upper triangular matrix)


gen_ind2 <- lower.tri(gen_raw2, diag = TRUE)

# select values of interest and replace rest with NAs

gen_tri2 <- gen_raw2
gen_tri2[gen_ind2 == FALSE] <- NA


#melt data

gen_long2 <-
  melt(
    data = gen_tri2,
    id.vars = "patch_start",
    variable.name = "patch_end",
    value.name = "gen_distance",
    na.rm = TRUE
  )

gen_long2$Tree_disease <- 70
gen_long2$Management <- 0
gen_long2$Replicate <- 1



#Test 3 ####
Test3_samples <- read.table("Batch13_Sim1_Land3_GenSamples.txt", header=TRUE)
T3_subset <- subset(Test3_samples, Test3_samples$Year == 90)
T3_subset <- subset(T3_subset, T3_subset$Rep == 1)
location_3 <- unname(unlist(T3_subset["PatchID"]))

T3_formatted <- T3_subset
T3_formatted <- T3_formatted[,6:35]+150
T3_formatted$loc0 <- paste(T3_formatted$Chr0Loc0Allele0, "/", T3_formatted$Chr0Loc0Allele1)
T3_formatted$loc1 <- paste(T3_formatted$Chr0Loc1Allele0, "/", T3_formatted$Chr0Loc1Allele1)
T3_formatted$loc2 <- paste(T3_formatted$Chr0Loc2Allele0, "/", T3_formatted$Chr0Loc2Allele1)
T3_formatted$loc3 <- paste(T3_formatted$Chr0Loc3Allele0, "/", T3_formatted$Chr0Loc3Allele1)
T3_formatted$loc4 <- paste(T3_formatted$Chr0Loc4Allele0, "/", T3_formatted$Chr0Loc4Allele1)
T3_formatted$loc5 <- paste(T3_formatted$Chr0Loc5Allele0, "/", T3_formatted$Chr0Loc5Allele1)
T3_formatted$loc6 <- paste(T3_formatted$Chr0Loc6Allele0, "/", T3_formatted$Chr0Loc6Allele1)
T3_formatted$loc7 <- paste(T3_formatted$Chr0Loc7Allele0, "/", T3_formatted$Chr0Loc7Allele1)
T3_formatted$loc8 <- paste(T3_formatted$Chr0Loc8Allele0, "/", T3_formatted$Chr0Loc8Allele1)
T3_formatted$loc9 <- paste(T3_formatted$Chr0Loc9Allele0, "/", T3_formatted$Chr0Loc9Allele1)
T3_formatted$loc10 <- paste(T3_formatted$Chr0Loc10Allele0, "/", T3_formatted$Chr0Loc10Allele1)
T3_formatted$loc11 <- paste(T3_formatted$Chr0Loc11Allele0, "/", T3_formatted$Chr0Loc11Allele1)
T3_formatted$loc12 <- paste(T3_formatted$Chr0Loc12Allele0, "/", T3_formatted$Chr0Loc12Allele1)
T3_formatted$loc13 <- paste(T3_formatted$Chr0Loc13Allele0, "/", T3_formatted$Chr0Loc13Allele1)
T3_formatted$loc14 <- paste(T3_formatted$Chr0Loc14Allele0, "/", T3_formatted$Chr0Loc14Allele1)
T3_formatted <- T3_formatted[,31:45]
T_genind_3 <- df2genind(T3_formatted, pop = location_3, ncode = 3, sep = "/")

T3_raw_gen <- genet.dist(T_genind_3, method = "WC84")

gen_raw3 <- as.matrix(genet.dist(T_genind_3, method = "WC84"))
gen_raw3 <- as.data.frame(gen_raw3)
setDT(gen_raw3, keep.rownames = TRUE)[]

# prepare data
# csv data has useful information in the row.names
# this has been loaded as the first column
# rename as patch_start

colnames(gen_raw3)[1] <- "patch_start"


# convert patch info to factor
gen_raw3$patch_start <- as.factor(gen_raw3$patch_start)


# data is a triangular matrix
# let's just keep the lower half & convert to lower triangular matrix
# note: we need the diagonal
# (could alternatively use upper triangular matrix)


gen_ind3 <- lower.tri(gen_raw3, diag = TRUE)

# select values of interest and replace rest with NAs

gen_tri3 <- gen_raw3
gen_tri3[gen_ind3 == FALSE] <- NA


#melt data

gen_long3 <-
  melt(
    data = gen_tri3,
    id.vars = "patch_start",
    variable.name = "patch_end",
    value.name = "gen_distance",
    na.rm = TRUE
  )

gen_long3$Tree_disease <- 70
gen_long3$Management <- 40
gen_long3$Replicate <- 1




#Test 4 ####
Test4_samples <- read.table("Batch13_Sim1_Land4_GenSamples.txt", header=TRUE)
T4_subset <- subset(Test4_samples, Test4_samples$Year == 90)
T4_subset <- subset(T4_subset, T4_subset$Rep == 1)
location_4 <- unname(unlist(T4_subset["PatchID"]))

T4_formatted <- T4_subset
T4_formatted <- T4_formatted[,6:35]+150
T4_formatted$loc0 <- paste(T4_formatted$Chr0Loc0Allele0, "/", T4_formatted$Chr0Loc0Allele1)
T4_formatted$loc1 <- paste(T4_formatted$Chr0Loc1Allele0, "/", T4_formatted$Chr0Loc1Allele1)
T4_formatted$loc2 <- paste(T4_formatted$Chr0Loc2Allele0, "/", T4_formatted$Chr0Loc2Allele1)
T4_formatted$loc3 <- paste(T4_formatted$Chr0Loc3Allele0, "/", T4_formatted$Chr0Loc3Allele1)
T4_formatted$loc4 <- paste(T4_formatted$Chr0Loc4Allele0, "/", T4_formatted$Chr0Loc4Allele1)
T4_formatted$loc5 <- paste(T4_formatted$Chr0Loc5Allele0, "/", T4_formatted$Chr0Loc5Allele1)
T4_formatted$loc6 <- paste(T4_formatted$Chr0Loc6Allele0, "/", T4_formatted$Chr0Loc6Allele1)
T4_formatted$loc7 <- paste(T4_formatted$Chr0Loc7Allele0, "/", T4_formatted$Chr0Loc7Allele1)
T4_formatted$loc8 <- paste(T4_formatted$Chr0Loc8Allele0, "/", T4_formatted$Chr0Loc8Allele1)
T4_formatted$loc9 <- paste(T4_formatted$Chr0Loc9Allele0, "/", T4_formatted$Chr0Loc9Allele1)
T4_formatted$loc10 <- paste(T4_formatted$Chr0Loc10Allele0, "/", T4_formatted$Chr0Loc10Allele1)
T4_formatted$loc11 <- paste(T4_formatted$Chr0Loc11Allele0, "/", T4_formatted$Chr0Loc11Allele1)
T4_formatted$loc12 <- paste(T4_formatted$Chr0Loc12Allele0, "/", T4_formatted$Chr0Loc12Allele1)
T4_formatted$loc13 <- paste(T4_formatted$Chr0Loc13Allele0, "/", T4_formatted$Chr0Loc13Allele1)
T4_formatted$loc14 <- paste(T4_formatted$Chr0Loc14Allele0, "/", T4_formatted$Chr0Loc14Allele1)
T4_formatted <- T4_formatted[,31:45]
T_genind_4 <- df2genind(T4_formatted, pop = location_4, ncode = 3, sep = "/")

T4_raw_gen <- genet.dist(T_genind_4, method = "WC84")

gen_raw4 <- as.matrix(genet.dist(T_genind_4, method = "WC84"))
gen_raw4 <- as.data.frame(gen_raw4)
setDT(gen_raw4, keep.rownames = TRUE)[]

# prepare data
# csv data has useful information in the row.names
# this has been loaded as the first column
# rename as patch_start

colnames(gen_raw4)[1] <- "patch_start"


# convert patch info to factor
gen_raw4$patch_start <- as.factor(gen_raw4$patch_start)


# data is a triangular matrix
# let's just keep the lower half & convert to lower triangular matrix
# note: we need the diagonal
# (could alternatively use upper triangular matrix)


gen_ind4 <- lower.tri(gen_raw4, diag = TRUE)

# select values of interest and replace rest with NAs

gen_tri4 <- gen_raw4
gen_tri4[gen_ind4 == FALSE] <- NA


#melt data

gen_long4 <-
  melt(
    data = gen_tri4,
    id.vars = "patch_start",
    variable.name = "patch_end",
    value.name = "gen_distance",
    na.rm = TRUE
  )

gen_long4$Tree_disease <- 70
gen_long4$Management <- 80
gen_long4$Replicate <- 1





#Test 5 ####
Test5_samples <- read.table("Batch13_Sim1_Land5_GenSamples.txt", header=TRUE)
T5_subset <- subset(Test5_samples, Test5_samples$Year == 90)
T5_subset <- subset(T5_subset, T5_subset$Rep == 1)
location_5 <- unname(unlist(T5_subset["PatchID"]))

T5_formatted <- T5_subset
T5_formatted <- T5_formatted[,6:35]+150
T5_formatted$loc0 <- paste(T5_formatted$Chr0Loc0Allele0, "/", T5_formatted$Chr0Loc0Allele1)
T5_formatted$loc1 <- paste(T5_formatted$Chr0Loc1Allele0, "/", T5_formatted$Chr0Loc1Allele1)
T5_formatted$loc2 <- paste(T5_formatted$Chr0Loc2Allele0, "/", T5_formatted$Chr0Loc2Allele1)
T5_formatted$loc3 <- paste(T5_formatted$Chr0Loc3Allele0, "/", T5_formatted$Chr0Loc3Allele1)
T5_formatted$loc4 <- paste(T5_formatted$Chr0Loc4Allele0, "/", T5_formatted$Chr0Loc4Allele1)
T5_formatted$loc5 <- paste(T5_formatted$Chr0Loc5Allele0, "/", T5_formatted$Chr0Loc5Allele1)
T5_formatted$loc6 <- paste(T5_formatted$Chr0Loc6Allele0, "/", T5_formatted$Chr0Loc6Allele1)
T5_formatted$loc7 <- paste(T5_formatted$Chr0Loc7Allele0, "/", T5_formatted$Chr0Loc7Allele1)
T5_formatted$loc8 <- paste(T5_formatted$Chr0Loc8Allele0, "/", T5_formatted$Chr0Loc8Allele1)
T5_formatted$loc9 <- paste(T5_formatted$Chr0Loc9Allele0, "/", T5_formatted$Chr0Loc9Allele1)
T5_formatted$loc10 <- paste(T5_formatted$Chr0Loc10Allele0, "/", T5_formatted$Chr0Loc10Allele1)
T5_formatted$loc11 <- paste(T5_formatted$Chr0Loc11Allele0, "/", T5_formatted$Chr0Loc11Allele1)
T5_formatted$loc12 <- paste(T5_formatted$Chr0Loc12Allele0, "/", T5_formatted$Chr0Loc12Allele1)
T5_formatted$loc13 <- paste(T5_formatted$Chr0Loc13Allele0, "/", T5_formatted$Chr0Loc13Allele1)
T5_formatted$loc14 <- paste(T5_formatted$Chr0Loc14Allele0, "/", T5_formatted$Chr0Loc14Allele1)
T5_formatted <- T5_formatted[,31:45]
T_genind_5 <- df2genind(T5_formatted, pop = location_5, ncode = 3, sep = "/")

T5_raw_gen <- genet.dist(T_genind_5, method = "WC84")

gen_raw5 <- as.matrix(genet.dist(T_genind_5, method = "WC84"))
gen_raw5 <- as.data.frame(gen_raw5)
setDT(gen_raw5, keep.rownames = TRUE)[]

# prepare data
# csv data has useful information in the row.names
# this has been loaded as the first column
# rename as patch_start

colnames(gen_raw5)[1] <- "patch_start"


# convert patch info to factor
gen_raw5$patch_start <- as.factor(gen_raw5$patch_start)


# data is a triangular matrix
# let's just keep the lower half & convert to lower triangular matrix
# note: we need the diagonal
# (could alternatively use upper triangular matrix)


gen_ind5 <- lower.tri(gen_raw5, diag = TRUE)

# select values of interest and replace rest with NAs

gen_tri5 <- gen_raw5
gen_tri5[gen_ind5 == FALSE] <- NA


#melt data

gen_long5 <-
  melt(
    data = gen_tri5,
    id.vars = "patch_start",
    variable.name = "patch_end",
    value.name = "gen_distance",
    na.rm = TRUE
  )

gen_long5$Tree_disease <- 80
gen_long5$Management <- 0
gen_long5$Replicate <- 1


#Test 6 ####
Test6_samples <- read.table("Batch13_Sim1_Land6_GenSamples.txt", header=TRUE)
T6_subset <- subset(Test6_samples, Test6_samples$Year == 90)
T6_subset <- subset(T6_subset, T6_subset$Rep == 1)
location_6 <- unname(unlist(T6_subset["PatchID"]))

T6_formatted <- T6_subset
T6_formatted <- T6_formatted[,6:35]+150
T6_formatted$loc0 <- paste(T6_formatted$Chr0Loc0Allele0, "/", T6_formatted$Chr0Loc0Allele1)
T6_formatted$loc1 <- paste(T6_formatted$Chr0Loc1Allele0, "/", T6_formatted$Chr0Loc1Allele1)
T6_formatted$loc2 <- paste(T6_formatted$Chr0Loc2Allele0, "/", T6_formatted$Chr0Loc2Allele1)
T6_formatted$loc3 <- paste(T6_formatted$Chr0Loc3Allele0, "/", T6_formatted$Chr0Loc3Allele1)
T6_formatted$loc4 <- paste(T6_formatted$Chr0Loc4Allele0, "/", T6_formatted$Chr0Loc4Allele1)
T6_formatted$loc5 <- paste(T6_formatted$Chr0Loc5Allele0, "/", T6_formatted$Chr0Loc5Allele1)
T6_formatted$loc6 <- paste(T6_formatted$Chr0Loc6Allele0, "/", T6_formatted$Chr0Loc6Allele1)
T6_formatted$loc7 <- paste(T6_formatted$Chr0Loc7Allele0, "/", T6_formatted$Chr0Loc7Allele1)
T6_formatted$loc8 <- paste(T6_formatted$Chr0Loc8Allele0, "/", T6_formatted$Chr0Loc8Allele1)
T6_formatted$loc9 <- paste(T6_formatted$Chr0Loc9Allele0, "/", T6_formatted$Chr0Loc9Allele1)
T6_formatted$loc10 <- paste(T6_formatted$Chr0Loc10Allele0, "/", T6_formatted$Chr0Loc10Allele1)
T6_formatted$loc11 <- paste(T6_formatted$Chr0Loc11Allele0, "/", T6_formatted$Chr0Loc11Allele1)
T6_formatted$loc12 <- paste(T6_formatted$Chr0Loc12Allele0, "/", T6_formatted$Chr0Loc12Allele1)
T6_formatted$loc13 <- paste(T6_formatted$Chr0Loc13Allele0, "/", T6_formatted$Chr0Loc13Allele1)
T6_formatted$loc14 <- paste(T6_formatted$Chr0Loc14Allele0, "/", T6_formatted$Chr0Loc14Allele1)
T6_formatted <- T6_formatted[,31:45]
T_genind_6 <- df2genind(T6_formatted, pop = location_6, ncode = 3, sep = "/")

T6_raw_gen <- genet.dist(T_genind_6, method = "WC84")

gen_raw6 <- as.matrix(genet.dist(T_genind_6, method = "WC84"))
gen_raw6 <- as.data.frame(gen_raw6)
setDT(gen_raw6, keep.rownames = TRUE)[]

# prepare data
# csv data has useful information in the row.names
# this has been loaded as the first column
# rename as patch_start

colnames(gen_raw6)[1] <- "patch_start"


# convert patch info to factor
gen_raw6$patch_start <- as.factor(gen_raw6$patch_start)


# data is a triangular matrix
# let's just keep the lower half & convert to lower triangular matrix
# note: we need the diagonal
# (could alternatively use upper triangular matrix)


gen_ind6 <- lower.tri(gen_raw6, diag = TRUE)

# select values of interest and replace rest with NAs

gen_tri6 <- gen_raw6
gen_tri6[gen_ind6 == FALSE] <- NA


#melt data

gen_long6 <-
  melt(
    data = gen_tri6,
    id.vars = "patch_start",
    variable.name = "patch_end",
    value.name = "gen_distance",
    na.rm = TRUE
  )

gen_long6$Tree_disease <- 80
gen_long6$Management <- 40
gen_long6$Replicate <- 1



#Test 7 ####
Test7_samples <- read.table("Batch13_Sim1_Land7_GenSamples.txt", header=TRUE)
T7_subset <- subset(Test7_samples, Test7_samples$Year == 90)
T7_subset <- subset(T7_subset, T7_subset$Rep == 1)
location_7 <- unname(unlist(T7_subset["PatchID"]))

T7_formatted <- T7_subset
T7_formatted <- T7_formatted[,6:35]+150
T7_formatted$loc0 <- paste(T7_formatted$Chr0Loc0Allele0, "/", T7_formatted$Chr0Loc0Allele1)
T7_formatted$loc1 <- paste(T7_formatted$Chr0Loc1Allele0, "/", T7_formatted$Chr0Loc1Allele1)
T7_formatted$loc2 <- paste(T7_formatted$Chr0Loc2Allele0, "/", T7_formatted$Chr0Loc2Allele1)
T7_formatted$loc3 <- paste(T7_formatted$Chr0Loc3Allele0, "/", T7_formatted$Chr0Loc3Allele1)
T7_formatted$loc4 <- paste(T7_formatted$Chr0Loc4Allele0, "/", T7_formatted$Chr0Loc4Allele1)
T7_formatted$loc5 <- paste(T7_formatted$Chr0Loc5Allele0, "/", T7_formatted$Chr0Loc5Allele1)
T7_formatted$loc6 <- paste(T7_formatted$Chr0Loc6Allele0, "/", T7_formatted$Chr0Loc6Allele1)
T7_formatted$loc7 <- paste(T7_formatted$Chr0Loc7Allele0, "/", T7_formatted$Chr0Loc7Allele1)
T7_formatted$loc8 <- paste(T7_formatted$Chr0Loc8Allele0, "/", T7_formatted$Chr0Loc8Allele1)
T7_formatted$loc9 <- paste(T7_formatted$Chr0Loc9Allele0, "/", T7_formatted$Chr0Loc9Allele1)
T7_formatted$loc10 <- paste(T7_formatted$Chr0Loc10Allele0, "/", T7_formatted$Chr0Loc10Allele1)
T7_formatted$loc11 <- paste(T7_formatted$Chr0Loc11Allele0, "/", T7_formatted$Chr0Loc11Allele1)
T7_formatted$loc12 <- paste(T7_formatted$Chr0Loc12Allele0, "/", T7_formatted$Chr0Loc12Allele1)
T7_formatted$loc13 <- paste(T7_formatted$Chr0Loc13Allele0, "/", T7_formatted$Chr0Loc13Allele1)
T7_formatted$loc14 <- paste(T7_formatted$Chr0Loc14Allele0, "/", T7_formatted$Chr0Loc14Allele1)
T7_formatted <- T7_formatted[,31:45]
T_genind_7 <- df2genind(T7_formatted, pop = location_7, ncode = 3, sep = "/")

T7_raw_gen <- genet.dist(T_genind_7, method = "WC84")

gen_raw7 <- as.matrix(genet.dist(T_genind_7, method = "WC84"))
gen_raw7 <- as.data.frame(gen_raw7)
setDT(gen_raw7, keep.rownames = TRUE)[]

# prepare data
# csv data has useful information in the row.names
# this has been loaded as the first column
# rename as patch_start

colnames(gen_raw7)[1] <- "patch_start"


# convert patch info to factor
gen_raw7$patch_start <- as.factor(gen_raw7$patch_start)


# data is a triangular matrix
# let's just keep the lower half & convert to lower triangular matrix
# note: we need the diagonal
# (could alternatively use upper triangular matrix)


gen_ind7 <- lower.tri(gen_raw7, diag = TRUE)

# select values of interest and replace rest with NAs

gen_tri7 <- gen_raw7
gen_tri7[gen_ind7 == FALSE] <- NA


#melt data

gen_long7 <-
  melt(
    data = gen_tri7,
    id.vars = "patch_start",
    variable.name = "patch_end",
    value.name = "gen_distance",
    na.rm = TRUE
  )

gen_long7$Tree_disease <- 80
gen_long7$Management <- 80
gen_long7$Replicate <- 1


##### bind data ####

gen_data <- rbind(gen_long1, gen_long2, gen_long3, gen_long4, gen_long5, gen_long6, gen_long7)








### RangeShifter ####

Test1_RS <- read.table("Batch13_Sim1_Land1_LandGen.txt", header=TRUE)
Test2_RS <- read.table("Batch13_Sim1_Land2_LandGen.txt", header=TRUE)
Test3_RS<- read.table("Batch13_Sim1_Land3_LandGen.txt", header=TRUE)
Test4_RS <- read.table("Batch13_Sim1_Land4_LandGen.txt", header=TRUE)
Test5_RS <- read.table("Batch13_Sim1_Land5_LandGen.txt", header=TRUE)
Test6_RS <- read.table("Batch13_Sim1_Land6_LandGen.txt", header=TRUE)
Test7_RS<- read.table("Batch13_Sim1_Land7_LandGen.txt", header=TRUE)

T1 <- subset(Test1_RS, Test1_RS$Year == 90)
T1 <- subset(T1, T1$Rep == 1)
T1$Tree_disease <- 0
T1$Management <- 0
T1$Replicate <- 1
T2 <- subset(Test2_RS, Test2_RS$Year == 90)
T2 <- subset(T2, T2$Rep == 1)
T2$Tree_disease <- 70
T2$Management <- 0
T2$Replicate <- 1
T3 <- subset(Test3_RS, Test3_RS$Year == 90)
T3 <- subset(T3, T3$Rep == 1)
T3$Tree_disease <- 70
T3$Management <- 40
T3$Replicate <- 1
T4 <- subset(Test4_RS, Test4_RS$Year == 90)
T4 <- subset(T4, T4$Rep == 1)
T4$Tree_disease <- 70
T4$Management <- 80
T4$Replicate <- 1
T5 <- subset(Test5_RS, Test5_RS$Year == 90)
T5 <- subset(T5, T5$Rep == 1)
T5$Tree_disease <- 80
T5$Management <- 0
T5$Replicate <- 1
T6 <- subset(Test6_RS, Test6_RS$Year == 90)
T6 <- subset(T6, T6$Rep == 1)
T6$Tree_disease <- 80
T6$Management <- 40
T6$Replicate <- 1
T7 <- subset(Test7_RS, Test7_RS$Year == 90)
T7 <- subset(T7, T7$Rep == 1)
T7$Tree_disease <- 80
T7$Management <- 80
T7$Replicate <- 1

Data_RS <- rbind(T1, T2, T3, T4, T5, T6, T7)
Data_RS <- subset(Data_RS, Data_RS$PatchID0 > 0)
Data_RS <- subset(Data_RS, Data_RS$PatchID1 > 0)

Data_RS_subset <- Data_RS[, 3:14]
Data_RS_subset <- Data_RS[, -c(1, 2, 5, 6, 7, 8, 10, 11)]
library(plyr) 
gen_data_names <- rename(gen_data, c("patch_start" = "PatchID1", "patch_end" = "PatchID0", "gen_distance" = "FST"))

Data_RS_subset$method <- "RS"
gen_data_names$method <- "R"

finaldata <- rbind(Data_RS_subset, gen_data_names)

library(ggplot2)

finaldata$PatchID0 <- factor(finaldata$PatchID0)
finaldata$PatchID1 <- factor(finaldata$PatchID1)
finaldata$method <- factor(finaldata$method)
finaldata$Tree_disease <- factor(finaldata$Tree_disease)
finaldata$Management <- factor(finaldata$Management)


ggplot(finaldata, aes(x=Tree_disease, y=FST, fill=method)) +  geom_boxplot()+
  facet_wrap(~Management)

write.csv(finaldata, file = "comparison_data.csv")

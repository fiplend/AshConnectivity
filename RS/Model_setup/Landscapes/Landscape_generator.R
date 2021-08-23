# Code for Ash maps and tree removal
# Author: FP
# Date created: 30_03_2020  

##function that removes ash trees
## tree disease treatment: each ash tree has a xx% of ash dieback
## Management scenario: Ash surrounding infected ash trees on roadsides are removed until XX% of infected roadside trees have been infected
## Woodland ash once infected is converted to woodland cell
## Matrix and roadside ash are both converted to matrix cells once removed or felled

library(dplyr)
library(raster)
library(igraph)
library(rgdal)

AssigningAsh <- function(hab_map, prob_ash) {
  hab_map[hab_map[] == 3] <- sample(
    c(3, 6),
    size = length(hab_map[hab_map[] ==
                            3]),
    #create vector of the same lengtht to replace current values with
    replace  = T,
    prob = c(1 - prob_ash, prob_ash)
  )  #90% of being a 3s (broadleaf) in raster with 3 and 10% of being a 6 (Ash)
  hab_map[hab_map[] == 4] <- sample(
    c(4, 7),
    size = length(hab_map[hab_map[] ==
                            4]),
    replace  = T,
    prob = c(1 - prob_ash, prob_ash)
  )
  hab_map[hab_map[] == 5] <- sample(
    c(5, 8),
    size = length(hab_map[hab_map[] ==
                            5]),
    replace  = T,
    prob = c(1 - prob_ash, prob_ash)
  )
  BasilineMap <- hab_map
  return(BasilineMap)
}
killingAsh <-
  # function that applies treatment to baseline maps - produces 6 maps in a raster brick
  function(hab_map,
           hab_baseline,
           prob_infection1,
           prob_infection2,
           prob_infection3,
           prob_infection4,
           prob_infection5,
           prob_infection6,
           prob_infection7,
           prob_infection8,
           prob_fell,
           prob_fell2,
           buffer_distance) {
    # set seed for reproducibility
    set.seed(0)
    
    #control treatment - no disease 
    hab_scenario0 <- hab_baseline
    
    #### corresponding patch map
    hab_baseline_patch <- hab_baseline
    hab_baseline_patch[hab_baseline_patch[] == 1] <-
      0 ### converting all matrix cells to 0
    hab_baseline_patch[hab_baseline_patch[] == 2] <- 0
    Patch_raster <- clump(hab_baseline_patch)
    ## Make an NA-value raster based on the LC raster attributes
    formask <- setValues(raster(Patch_raster), NA)
    ## Assign 1 to formask to all cells corresponding to the forest class
    formask[hab_baseline == 2] <- 1
    clumpFreq <- freq(Patch_raster)
    ## Coerce freq table to data.frame
    clumpFreq <- as.data.frame(clumpFreq)
    ## which rows of the data.frame are only represented by xx cell?
    str(which(clumpFreq$count < 5)) 
    ## which values do these correspond to?
    str(clumpFreq$value[which(clumpFreq$count < 5)])
    ## Put these into a vector of clump ID's to be removed
    excludeID <- clumpFreq$value[which(clumpFreq$count < 5)]
    ## Make a new forest mask to be sieved
    formaskSieve <- Patch_raster
    ## Assign NA to all clumps whose IDs are found in excludeID
    formaskSieve[Patch_raster %in% excludeID] <- NA
    
    ashcountstats <-
      data.frame(zonal(
        match(hab_scenario0, c(6, 7, 8)),
        formaskSieve,
        #each forest patch is a zone
        fun = 'count',
        digits = 0,
        na.rm = TRUE
      ))
    cellcount <-
      data.frame(freq(formaskSieve, useNA = "no", merge = TRUE))
    
    #create replicate raster to work with
    Patch_scenario0 <- formaskSieve
    excludeIDash <-
      ashcountstats$zone[which(ashcountstats$count < 1)]  ## count of ash trees and update ID list
    
    ## Assign NA to all clumps whose IDs are found in excludeID
    Patch_scenario0[formaskSieve %in% excludeIDash] <- NA
    
    
    
    ## Tree disease 10 ####
    
    hab_scenario1 <- hab_baseline
    # replace elements of raster with dead ash trees
    hab_scenario1[hab_scenario1[] == 6] <-
      sample(
        c(6, 9),
        size = length(hab_scenario1[hab_scenario1[] ==
                                      6]),
        replace  = T,
        prob = c(1 - prob_infection1, prob_infection1)
      )
    hab_scenario1[hab_scenario1[] == 7] <-
      sample(
        c(7, 10),
        size = length(hab_scenario1[hab_scenario1[] ==
                                      7]),
        replace  = T,
        prob = c(1 - prob_infection1, prob_infection1)
      )
    hab_scenario1[hab_scenario1[] == 8] <-
      sample(
        c(8, 9),
        size = length(hab_scenario1[hab_scenario1[] ==
                                      8]),
        replace  = T,
        prob = c(1 - prob_infection1, prob_infection1)
      )
    
    
    ashcountstats1 <-
      data.frame(zonal(
        match(hab_scenario1, c(6, 7, 8)),
        formaskSieve,
        #each forest patch is a zone
        fun = 'count',
        digits = 0,
        na.rm = TRUE
      ))
    cellcount <-
      data.frame(freq(formaskSieve, useNA = "no", merge = TRUE))
    
    #create replicate raster to work with
    Patch_scenario1 <- formaskSieve
    excludeIDash1 <-
      ashcountstats1$zone[which(ashcountstats1$count < 1)]  ## count of ash trees and update ID list
    
    ## Assign NA to all clumps whose IDs are found in excludeID
    Patch_scenario1[formaskSieve %in% excludeIDash1] <- NA
    
    
    ################ with management 1 
    
    # create raster which shows just the roadside tree patches and the infected trees
    hab_scenario2 <- hab_scenario1
    while (length(hab_scenario2[hab_scenario2[] == 11]) < prob_fell * length(hab_scenario1[hab_scenario1[] == 10])) {
      hab_road <- hab_map
      hab_road[hab_road[] != 4] <- 0
      hab_road[hab_road[] == 4] <- 1
      Habitat_infect_road <- hab_road    #create copy of raod trees
      Habitat_infect_road[hab_scenario2[] == 10] <-
        2  #mark all infected road trees
      Habitat_infect_road[hab_road[] != 1] <-
        0   # remove and non-roadside infected trees
      Habitat_road_infected_buf <- Habitat_infect_road
      Habitat_road_infected_buf[Habitat_infect_road[] != 2] <-
        NA   #all cells except NA become NA
      Habitat_road_infected_buf[Habitat_infect_road[] == 2] <-
        1  # infected cells now 1
      Habitat_road_infected_buffer <-
        Habitat_road_infected_buf #copy
      Habitat_road_infected_buffer[Habitat_road_infected_buf[] == 1] <-
        NA  #empty version of map
      Habitat_road_infected_buffer[sample(which(values(Habitat_road_infected_buf) ==
                                                  1), 1)] <-
        1 #converts random tree infected cell to one
      Habitat_road_infected_buffer <-
        buffer(Habitat_road_infected_buffer, width = buffer_distance) # add a 10 m buffer
      #now remove all ash trees that fall within specified buffer and create new raster holding this treatment
      Habitat_road_infected_markedtrees <-
        Habitat_road_infected_buffer
      Habitat_road_infected_markedtrees[hab_scenario2[] == 7] <-
        2  #mark all ash trees within the buffer
      Habitat_road_infected_markedtrees[hab_scenario2[] == 10] <-
        3
      Habitat_road_infected_markedtrees[is.na(Habitat_road_infected_buffer[])] <-
        0
      hab_scenario2[Habitat_road_infected_markedtrees[] == 3] <-
        11 # 11 = felled and infected
      hab_scenario2[Habitat_road_infected_markedtrees[] == 2] <-
        12 # 12 = felled and healthy
    }
    
    #### corresponding patch map
    
    ashcountstats2  <-
      #count of ash trees in each forest patch that are 6, 7 and 8 (ash trees)
      data.frame(zonal(
        match(hab_scenario2, c(6, 7, 8)),
        formaskSieve,
        #each forest patch is a zone
        fun = 'count',
        digits = 0,
        na.rm = TRUE
      ))
    #create replicate raster to work with
    Patch_scenario2 <- formaskSieve
    excludeIDash2 <-
      ashcountstats2$zone[which(ashcountstats2$count < 1)]  ## count of ash trees and update ID list
    
    ## Assign NA to all clumps whose IDs are found in excludeID
    Patch_scenario2[formaskSieve %in% excludeIDash2] <- NA
    
    ################ with management 2
    
    # create raster which shows just the roadside tree patches and the infected trees
    hab_scenario3 <- hab_scenario1
    while (length(hab_scenario3[hab_scenario3[] == 11]) < prob_fell2 * length(hab_scenario1[hab_scenario1[] == 10])) {
      hab_road3 <- hab_map
      hab_road3[hab_road3[] != 4] <- 0
      hab_road3[hab_road3[] == 4] <- 1
      Habitat_infect_road3 <-
        hab_road3    #create copy of raod trees
      Habitat_infect_road3[hab_scenario3[] == 10] <-
        2  #mark all infected road trees
      Habitat_infect_road3[hab_road3[] != 1] <-
        0   # remove and non-roadside infected trees
      Habitat_road_infected_buf3 <- Habitat_infect_road3
      Habitat_road_infected_buf3[Habitat_infect_road3[] != 2] <-
        NA   #all cells except NA become NA
      Habitat_road_infected_buf3[Habitat_infect_road3[] == 2] <-
        1  # infected cells now 1
      Habitat_road_infected_buffer3 <-
        Habitat_road_infected_buf3 #copy
      Habitat_road_infected_buffer3[Habitat_road_infected_buf3[] == 1] <-
        NA  #empty version of map
      Habitat_road_infected_buffer3[sample(which(values(Habitat_road_infected_buf3) ==
                                                   1), 1)] <-
        1 #converts random tree infected cell to one
      Habitat_road_infected_buffer3 <-
        buffer(Habitat_road_infected_buffer3, width = buffer_distance) # add a 10 m buffer
      #now remove all ash trees that fall within specified buffer and create new raster holding this treatment
      Habitat_road_infected_markedtrees3 <-
        Habitat_road_infected_buffer3
      Habitat_road_infected_markedtrees3[hab_scenario3[] == 7] <-
        2  #mark all ash trees within the buffer
      Habitat_road_infected_markedtrees3[hab_scenario3[] == 10] <-
        3
      Habitat_road_infected_markedtrees3[is.na(Habitat_road_infected_buffer3[])] <-
        0
      hab_scenario3[Habitat_road_infected_markedtrees3[] == 3] <-
        11 # 11 = felled and infected
      hab_scenario3[Habitat_road_infected_markedtrees3[] == 2] <-
        12 # 12 = felled and healthy
    }
    
    #### corresponding patch map
    
    ashcountstats3  <-
      #count of ash trees in each forest patch that are 6, 7 and 8 (ash trees)
      data.frame(zonal(
        match(hab_scenario3, c(6, 7, 8)),
        formaskSieve,
        #each forest patch is a zone
        fun = 'count',
        digits = 0,
        na.rm = TRUE
      ))
    #create replicate raster to work with
    Patch_scenario3 <- formaskSieve
    excludeIDash3 <-
      ashcountstats3$zone[which(ashcountstats3$count < 1)]  ## count of ash trees and update ID list
    
    ## Assign NA to all clumps whose IDs are found in excludeID
    Patch_scenario3[formaskSieve %in% excludeIDash3] <- NA
    
    
     
    ##### Tree disease 20 ####
    
    
    hab_scenario4 <- hab_baseline
    # replace elements of raster with dead ash trees
    hab_scenario4[hab_scenario4[] == 6] <-
      sample(
        c(6, 9),
        size = length(hab_scenario4[hab_scenario4[] ==
                                      6]),
        replace  = T,
        prob = c(1 - prob_infection2, prob_infection2)
      )
    hab_scenario4[hab_scenario4[] == 7] <-
      sample(
        c(7, 10),
        size = length(hab_scenario4[hab_scenario4[] ==
                                      7]),
        replace  = T,
        prob = c(1 - prob_infection2, prob_infection2)
      )
    hab_scenario4[hab_scenario4[] == 8] <-
      sample(
        c(8, 9),
        size = length(hab_scenario4[hab_scenario4[] ==
                                      8]),
        replace  = T,
        prob = c(1 - prob_infection2, prob_infection2)
      )
    
    
    ashcountstats4 <-
      data.frame(zonal(
        match(hab_scenario4, c(6, 7, 8)),
        formaskSieve,
        #each forest patch is a zone
        fun = 'count',
        digits = 0,
        na.rm = TRUE
      ))
    cellcount <-
      data.frame(freq(formaskSieve, useNA = "no", merge = TRUE))
    
    #create replicate raster to work with
    Patch_scenario4 <- formaskSieve
    excludeIDash4 <-
      ashcountstats4$zone[which(ashcountstats4$count < 1)]  ## count of ash trees and update ID list
    
    ## Assign NA to all clumps whose IDs are found in excludeID
    Patch_scenario4[formaskSieve %in% excludeIDash4] <- NA
    
    ################ with management 1
    
    # create raster which shows just the roadside tree patches and the infected trees
    hab_scenario5 <- hab_scenario4
    while (length(hab_scenario5[hab_scenario5[] == 11]) < prob_fell * length(hab_scenario4[hab_scenario4[] == 10])) {
      hab_road5 <- hab_map
      hab_road5[hab_road5[] != 4] <- 0
      hab_road5[hab_road5[] == 4] <- 1
      Habitat_infect_road5 <-
        hab_road5    #create copy of raod trees
      Habitat_infect_road5[hab_scenario5[] == 10] <-
        2  #mark all infected road trees
      Habitat_infect_road5[hab_road5[] != 1] <-
        0   # remove and non-roadside infected trees
      Habitat_road_infected_buf5 <- Habitat_infect_road5
      Habitat_road_infected_buf5[Habitat_infect_road5[] != 2] <-
        NA   #all cells except NA become NA
      Habitat_road_infected_buf5[Habitat_infect_road5[] == 2] <-
        1  # infected cells now 1
      Habitat_road_infected_buffer5 <-
        Habitat_road_infected_buf5 #copy
      Habitat_road_infected_buffer5[Habitat_road_infected_buf5[] == 1] <-
        NA  #empty version of map
      Habitat_road_infected_buffer5[sample(which(values(Habitat_road_infected_buf5) ==
                                                   1), 1)] <-
        1 #converts random tree infected cell to one
      Habitat_road_infected_buffer5 <-
        buffer(Habitat_road_infected_buffer5, width = buffer_distance) # add a 10 m buffer
      #now remove all ash trees that fall within specified buffer and create new raster holding this treatment
      Habitat_road_infected_markedtrees5 <-
        Habitat_road_infected_buffer5
      Habitat_road_infected_markedtrees5[hab_scenario5[] == 7] <-
        2  #mark all ash trees within the buffer
      Habitat_road_infected_markedtrees5[hab_scenario5[] == 10] <-
        3
      Habitat_road_infected_markedtrees5[is.na(Habitat_road_infected_buffer5[])] <-
        0
      hab_scenario5[Habitat_road_infected_markedtrees5[] == 3] <-
        11 # 11 = felled and infected
      hab_scenario5[Habitat_road_infected_markedtrees5[] == 2] <-
        12 # 12 = felled and healthy
    }
    
    #### corresponding patch map
    
    ashcountstats5  <-
      #count of ash trees in each forest patch that are 6, 7 and 8 (ash trees)
      data.frame(zonal(
        match(hab_scenario5, c(6, 7, 8)),
        formaskSieve,
        #each forest patch is a zone
        fun = 'count',
        digits = 0,
        na.rm = TRUE
      ))
    #create replicate raster to work with
    Patch_scenario5 <- formaskSieve
    excludeIDash5 <-
      ashcountstats5$zone[which(ashcountstats5$count < 1)]  ## count of ash trees and update ID list
    
    ## Assign NA to all clumps whose IDs are found in excludeID
    Patch_scenario5[formaskSieve %in% excludeIDash5] <- NA
    
    ################ with management 2
    
    # create raster which shows just the roadside tree patches and the infected trees
    hab_scenario6 <- hab_scenario4
    while (length(hab_scenario6[hab_scenario6[] == 11]) < prob_fell2 * length(hab_scenario4[hab_scenario4[] == 10])) {
      hab_road6 <- hab_map
      hab_road6[hab_road6[] != 4] <- 0
      hab_road6[hab_road6[] == 4] <- 1
      Habitat_infect_road6 <-
        hab_road6    #create copy of raod trees
      Habitat_infect_road6[hab_scenario6[] == 10] <-
        2  #mark all infected road trees
      Habitat_infect_road6[hab_road6[] != 1] <-
        0   # remove and non-roadside infected trees
      Habitat_road_infected_buf6 <- Habitat_infect_road6
      Habitat_road_infected_buf6[Habitat_infect_road6[] != 2] <-
        NA   #all cells except NA become NA
      Habitat_road_infected_buf6[Habitat_infect_road6[] == 2] <-
        1  # infected cells now 1
      Habitat_road_infected_buffer6 <-
        Habitat_road_infected_buf6 #copy
      Habitat_road_infected_buffer6[Habitat_road_infected_buf6[] == 1] <-
        NA  #empty version of map
      Habitat_road_infected_buffer6[sample(which(values(Habitat_road_infected_buf6) ==
                                                   1), 1)] <-
        1 #converts random tree infected cell to one
      Habitat_road_infected_buffer6 <-
        buffer(Habitat_road_infected_buffer6, width = buffer_distance) # add a 10 m buffer
      #now remove all ash trees that fall within specified buffer and create new raster holding this treatment
      Habitat_road_infected_markedtrees6 <-
        Habitat_road_infected_buffer6
      Habitat_road_infected_markedtrees6[hab_scenario6[] == 7] <-
        2  #mark all ash trees within the buffer
      Habitat_road_infected_markedtrees6[hab_scenario6[] == 10] <-
        3
      Habitat_road_infected_markedtrees6[is.na(Habitat_road_infected_buffer6[])] <-
        0
      hab_scenario6[Habitat_road_infected_markedtrees6[] == 3] <-
        11 # 11 = felled and infected
      hab_scenario6[Habitat_road_infected_markedtrees6[] == 2] <-
        12 # 12 = felled and healthy
    }
    
    #### corresponding patch map
    
    ashcountstats6  <-
      #count of ash trees in each forest patch that are 6, 7 and 8 (ash trees)
      data.frame(zonal(
        match(hab_scenario6, c(6, 7, 8)),
        formaskSieve,
        #each forest patch is a zone
        fun = 'count',
        digits = 0,
        na.rm = TRUE
      ))
    #create replicate raster to work with
    Patch_scenario6 <- formaskSieve
    excludeIDash6 <-
      ashcountstats6$zone[which(ashcountstats6$count < 1)]  ## count of ash trees and update ID list
    
    ## Assign NA to all clumps whose IDs are found in excludeID
    Patch_scenario6[formaskSieve %in% excludeIDash6] <- NA
    
    ##### Tree disease 30 ####
    
    
    hab_scenario7 <- hab_baseline
    # replace elements of raster with dead ash trees
    hab_scenario7[hab_scenario7[] == 6] <-
      sample(
        c(6, 9),
        size = length(hab_scenario7[hab_scenario7[] ==
                                      6]),
        replace  = T,
        prob = c(1 - prob_infection3, prob_infection3)
      )
    hab_scenario7[hab_scenario7[] == 7] <-
      sample(
        c(7, 10),
        size = length(hab_scenario7[hab_scenario7[] ==
                                      7]),
        replace  = T,
        prob = c(1 - prob_infection3, prob_infection3)
      )
    hab_scenario7[hab_scenario7[] == 8] <-
      sample(
        c(8, 9),
        size = length(hab_scenario7[hab_scenario7[] ==
                                      8]),
        replace  = T,
        prob = c(1 - prob_infection3, prob_infection3)
      )
    
    
    ashcountstats7 <-
      data.frame(zonal(
        match(hab_scenario7, c(6, 7, 8)),
        formaskSieve,
        #each forest patch is a zone
        fun = 'count',
        digits = 0,
        na.rm = TRUE
      ))
    cellcount <-
      data.frame(freq(formaskSieve, useNA = "no", merge = TRUE))
    
    #create replicate raster to work with
    Patch_scenario7 <- formaskSieve
    excludeIDash7 <-
      ashcountstats7$zone[which(ashcountstats7$count < 1)]  ## count of ash trees and update ID list
    
    ## Assign NA to all clumps whose IDs are found in excludeID
    Patch_scenario7[formaskSieve %in% excludeIDash7] <- NA
    
    ################ with management 1
    
    # create raster which shows just the roadside tree patches and the infected trees
    hab_scenario8 <- hab_scenario7
    while (length(hab_scenario8[hab_scenario8[] == 11]) < prob_fell * length(hab_scenario7[hab_scenario7[] == 10])) {
      hab_road8 <- hab_map
      hab_road8[hab_road8[] != 4] <- 0
      hab_road8[hab_road8[] == 4] <- 1
      Habitat_infect_road8 <-
        hab_road8    #create copy of raod trees
      Habitat_infect_road8[hab_scenario8[] == 10] <-
        2  #mark all infected road trees
      Habitat_infect_road8[hab_road8[] != 1] <-
        0   # remove and non-roadside infected trees
      Habitat_road_infected_buf8 <- Habitat_infect_road8
      Habitat_road_infected_buf8[Habitat_infect_road8[] != 2] <-
        NA   #all cells except NA become NA
      Habitat_road_infected_buf8[Habitat_infect_road8[] == 2] <-
        1  # infected cells now 1
      Habitat_road_infected_buffer8 <-
        Habitat_road_infected_buf8 #copy
      Habitat_road_infected_buffer8[Habitat_road_infected_buf8[] == 1] <-
        NA  #empty version of map
      Habitat_road_infected_buffer8[sample(which(values(Habitat_road_infected_buf8) ==
                                                   1), 1)] <-
        1 #converts random tree infected cell to one
      Habitat_road_infected_buffer8 <-
        buffer(Habitat_road_infected_buffer8, width = buffer_distance) # add a 10 m buffer
      #now remove all ash trees that fall within specified buffer and create new raster holding this treatment
      Habitat_road_infected_markedtrees8 <-
        Habitat_road_infected_buffer8
      Habitat_road_infected_markedtrees8[hab_scenario8[] == 7] <-
        2  #mark all ash trees within the buffer
      Habitat_road_infected_markedtrees8[hab_scenario8[] == 10] <-
        3
      Habitat_road_infected_markedtrees8[is.na(Habitat_road_infected_buffer8[])] <-
        0
      hab_scenario8[Habitat_road_infected_markedtrees8[] == 3] <-
        11 # 11 = felled and infected
      hab_scenario8[Habitat_road_infected_markedtrees8[] == 2] <-
        12 # 12 = felled and healthy
    }
    
    #### corresponding patch map
    
    ashcountstats8  <-
      #count of ash trees in each forest patch that are 6, 7 and 8 (ash trees)
      data.frame(zonal(
        match(hab_scenario8, c(6, 7, 8)),
        formaskSieve,
        #each forest patch is a zone
        fun = 'count',
        digits = 0,
        na.rm = TRUE
      ))
    #create replicate raster to work with
    Patch_scenario8 <- formaskSieve
    excludeIDash8 <-
      ashcountstats8$zone[which(ashcountstats8$count < 1)]  ## count of ash trees and update ID list
    
    ## Assign NA to all clumps whose IDs are found in excludeID
    Patch_scenario8[formaskSieve %in% excludeIDash8] <- NA
    
    ################ with management 2
    
    # create raster which shows just the roadside tree patches and the infected trees
    hab_scenario9 <- hab_scenario7
    while (length(hab_scenario9[hab_scenario9[] == 11]) < prob_fell2 * length(hab_scenario7[hab_scenario7[] == 10])) {
      hab_road9 <- hab_map
      hab_road9[hab_road9[] != 4] <- 0
      hab_road9[hab_road9[] == 4] <- 1
      Habitat_infect_road9 <-
        hab_road9   #create copy of raod trees
      Habitat_infect_road9[hab_scenario9[] == 10] <-
        2  #mark all infected road trees
      Habitat_infect_road9[hab_road9[] != 1] <-
        0   # remove and non-roadside infected trees
      Habitat_road_infected_buf9 <- Habitat_infect_road9
      Habitat_road_infected_buf9[Habitat_infect_road9[] != 2] <-
        NA   #all cells except NA become NA
      Habitat_road_infected_buf9[Habitat_infect_road9[] == 2] <-
        1  # infected cells now 1
      Habitat_road_infected_buffer9 <-
        Habitat_road_infected_buf9 #copy
      Habitat_road_infected_buffer9[Habitat_road_infected_buf9[] == 1] <-
        NA  #empty version of map
      Habitat_road_infected_buffer9[sample(which(values(Habitat_road_infected_buf9) ==
                                                   1), 1)] <-
        1 #converts random tree infected cell to one
      Habitat_road_infected_buffer9 <-
        buffer(Habitat_road_infected_buffer9, width = buffer_distance) # add a 10 m buffer
      #now remove all ash trees that fall within specified buffer and create new raster holding this treatment
      Habitat_road_infected_markedtrees9 <-
        Habitat_road_infected_buffer9
      Habitat_road_infected_markedtrees9[hab_scenario9[] == 7] <-
        2  #mark all ash trees within the buffer
      Habitat_road_infected_markedtrees9[hab_scenario9[] == 10] <-
        3
      Habitat_road_infected_markedtrees9[is.na(Habitat_road_infected_buffer9[])] <-
        0
      hab_scenario9[Habitat_road_infected_markedtrees9[] == 3] <-
        11 # 11 = felled and infected
      hab_scenario9[Habitat_road_infected_markedtrees9[] == 2] <-
        12 # 12 = felled and healthy
    }
    
    #### corresponding patch map
    
    ashcountstats9  <-
      #count of ash trees in each forest patch that are 6, 7 and 8 (ash trees)
      data.frame(zonal(
        match(hab_scenario9, c(6, 7, 8)),
        formaskSieve,
        #each forest patch is a zone
        fun = 'count',
        digits = 0,
        na.rm = TRUE
      ))
    #create replicate raster to work with
    Patch_scenario9 <- formaskSieve
    excludeIDash9 <-
      ashcountstats9$zone[which(ashcountstats9$count < 1)]  ## count of ash trees and update ID list
    
    ## Assign NA to all clumps whose IDs are found in excludeID
    Patch_scenario9[formaskSieve %in% excludeIDash9] <- NA
    
    
    ##### Tree disease 40 ####
    
    
    hab_scenario10 <- hab_baseline
    # replace elements of raster with dead ash trees
    hab_scenario10[hab_scenario10[] == 6] <-
      sample(
        c(6, 9),
        size = length(hab_scenario10[hab_scenario10[] ==
                                      6]),
        replace  = T,
        prob = c(1 - prob_infection4, prob_infection4)
      )
    hab_scenario10[hab_scenario10[] == 7] <-
      sample(
        c(7, 10),
        size = length(hab_scenario10[hab_scenario10[] ==
                                      7]),
        replace  = T,
        prob = c(1 - prob_infection4, prob_infection4)
      )
    hab_scenario10[hab_scenario10[] == 8] <-
      sample(
        c(8, 9),
        size = length(hab_scenario10[hab_scenario10[] ==
                                      8]),
        replace  = T,
        prob = c(1 - prob_infection4, prob_infection4)
      )
    
    
    ashcountstats10 <-
      data.frame(zonal(
        match(hab_scenario10, c(6, 7, 8)),
        formaskSieve,
        #each forest patch is a zone
        fun = 'count',
        digits = 0,
        na.rm = TRUE
      ))
    cellcount <-
      data.frame(freq(formaskSieve, useNA = "no", merge = TRUE))
    
    #create replicate raster to work with
    Patch_scenario10 <- formaskSieve
    excludeIDash10 <-
      ashcountstats10$zone[which(ashcountstats10$count < 1)]  ## count of ash trees and update ID list
    
    ## Assign NA to all clumps whose IDs are found in excludeID
    Patch_scenario10[formaskSieve %in% excludeIDash10] <- NA
    
    ################ with management 1
    
    # create raster which shows just the roadside tree patches and the infected trees
    hab_scenario11 <- hab_scenario10
    while (length(hab_scenario11[hab_scenario11[] == 11]) < prob_fell * length(hab_scenario11[hab_scenario11[] == 10])) {
      hab_road11 <- hab_map
      hab_road11[hab_road11[] != 4] <- 0
      hab_road11[hab_road11[] == 4] <- 1
      Habitat_infect_road11 <-
        hab_road11    #create copy of raod trees
      Habitat_infect_road11[hab_scenario11[] == 10] <-
        2  #mark all infected road trees
      Habitat_infect_road11[hab_road11[] != 1] <-
        0   # remove and non-roadside infected trees
      Habitat_road_infected_buf11 <- Habitat_infect_road11
      Habitat_road_infected_buf11[Habitat_infect_road11[] != 2] <-
        NA   #all cells except NA become NA
      Habitat_road_infected_buf11[Habitat_infect_road11[] == 2] <-
        1  # infected cells now 1
      Habitat_road_infected_buffer11 <-
        Habitat_road_infected_buf11 #copy
      Habitat_road_infected_buffer11[Habitat_road_infected_buf11[] == 1] <-
        NA  #empty version of map
      Habitat_road_infected_buffer11[sample(which(values(Habitat_road_infected_buf11) ==
                                                   1), 1)] <-
        1 #converts random tree infected cell to one
      Habitat_road_infected_buffer11 <-
        buffer(Habitat_road_infected_buffer11, width = buffer_distance) # add a 10 m buffer
      #now remove all ash trees that fall within specified buffer and create new raster holding this treatment
      Habitat_road_infected_markedtrees11 <-
        Habitat_road_infected_buffer11
      Habitat_road_infected_markedtrees11[hab_scenario11[] == 7] <-
        2  #mark all ash trees within the buffer
      Habitat_road_infected_markedtrees11[hab_scenario11[] == 10] <-
        3
      Habitat_road_infected_markedtrees11[is.na(Habitat_road_infected_buffer11[])] <-
        0
      hab_scenario11[Habitat_road_infected_markedtrees11[] == 3] <-
        11 # 11 = felled and infected
      hab_scenario11[Habitat_road_infected_markedtrees11[] == 2] <-
        12 # 12 = felled and healthy
    }
    
    #### corresponding patch map
    
    ashcountstats11  <-
      #count of ash trees in each forest patch that are 6, 7 and 8 (ash trees)
      data.frame(zonal(
        match(hab_scenario11, c(6, 7, 8)),
        formaskSieve,
        #each forest patch is a zone
        fun = 'count',
        digits = 0,
        na.rm = TRUE
      ))
    #create replicate raster to work with
    Patch_scenario11 <- formaskSieve
    excludeIDash11 <-
      ashcountstats11$zone[which(ashcountstats11$count < 1)]  ## count of ash trees and update ID list
    
    ## Assign NA to all clumps whose IDs are found in excludeID
    Patch_scenario11[formaskSieve %in% excludeIDash11] <- NA
    
    ################ with management 2
    
    # create raster which shows just the roadside tree patches and the infected trees
    hab_scenario12 <- hab_scenario10
    while (length(hab_scenario12[hab_scenario12[] == 11]) < prob_fell2 * length(hab_scenario10[hab_scenario10[] == 10])) {
      hab_road12 <- hab_map
      hab_road12[hab_road12[] != 4] <- 0
      hab_road12[hab_road12[] == 4] <- 1
      Habitat_infect_road12 <-
        hab_road12  #create copy of raod trees
      Habitat_infect_road12[hab_scenario12[] == 10] <-
        2  #mark all infected road trees
      Habitat_infect_road12[hab_road12[] != 1] <-
        0   # remove and non-roadside infected trees
      Habitat_road_infected_buf12 <- Habitat_infect_road12
      Habitat_road_infected_buf12[Habitat_infect_road12[] != 2] <-
        NA   #all cells except NA become NA
      Habitat_road_infected_buf12[Habitat_infect_road12[] == 2] <-
        1  # infected cells now 1
      Habitat_road_infected_buffer12 <-
        Habitat_road_infected_buf12 #copy
      Habitat_road_infected_buffer12[Habitat_road_infected_buf12[] == 1] <-
        NA  #empty version of map
      Habitat_road_infected_buffer12[sample(which(values(Habitat_road_infected_buf12) ==
                                                   1), 1)] <-
        1 #converts random tree infected cell to one
      Habitat_road_infected_buffer12 <-
        buffer(Habitat_road_infected_buffer12, width = buffer_distance) # add a 10 m buffer
      #now remove all ash trees that fall within specified buffer and create new raster holding this treatment
      Habitat_road_infected_markedtrees12 <-
        Habitat_road_infected_buffer12
      Habitat_road_infected_markedtrees12[hab_scenario12[] == 7] <-
        2  #mark all ash trees within the buffer
      Habitat_road_infected_markedtrees12[hab_scenario12[] == 10] <-
        3
      Habitat_road_infected_markedtrees12[is.na(Habitat_road_infected_buffer12[])] <-
        0
      hab_scenario12[Habitat_road_infected_markedtrees12[] == 3] <-
        11 # 11 = felled and infected
      hab_scenario12[Habitat_road_infected_markedtrees12[] == 2] <-
        12 # 12 = felled and healthy
    }
    
    #### corresponding patch map
    
    ashcountstats12  <-
      #count of ash trees in each forest patch that are 6, 7 and 8 (ash trees)
      data.frame(zonal(
        match(hab_scenario12, c(6, 7, 8)),
        formaskSieve,
        #each forest patch is a zone
        fun = 'count',
        digits = 0,
        na.rm = TRUE
      ))
    #create replicate raster to work with
    Patch_scenario12 <- formaskSieve
    excludeIDash12 <-
      ashcountstats12$zone[which(ashcountstats12$count < 1)]  ## count of ash trees and update ID list
    
    ## Assign NA to all clumps whose IDs are found in excludeID
    Patch_scenario12[formaskSieve %in% excludeIDash12] <- NA
    
    
    
    
    ##### Tree disease 50 ####
    
    
    hab_scenario13 <- hab_baseline
    # replace elements of raster with dead ash trees
    hab_scenario13[hab_scenario13[] == 6] <-
      sample(
        c(6, 9),
        size = length(hab_scenario13[hab_scenario13[] ==
                                       6]),
        replace  = T,
        prob = c(1 - prob_infection5, prob_infection5)
      )
    hab_scenario13[hab_scenario13[] == 7] <-
      sample(
        c(7, 10),
        size = length(hab_scenario13[hab_scenario13[] ==
                                       7]),
        replace  = T,
        prob = c(1 - prob_infection5, prob_infection5)
      )
    hab_scenario13[hab_scenario13[] == 8] <-
      sample(
        c(8, 9),
        size = length(hab_scenario13[hab_scenario13[] ==
                                       8]),
        replace  = T,
        prob = c(1 - prob_infection5, prob_infection5)
      )
    
    
    ashcountstats13 <-
      data.frame(zonal(
        match(hab_scenario13, c(6, 7, 8)),
        formaskSieve,
        #each forest patch is a zone
        fun = 'count',
        digits = 0,
        na.rm = TRUE
      ))
    cellcount <-
      data.frame(freq(formaskSieve, useNA = "no", merge = TRUE))
    
    #create replicate raster to work with
    Patch_scenario13 <- formaskSieve
    excludeIDash13 <-
      ashcountstats13$zone[which(ashcountstats13$count < 1)]  ## count of ash trees and update ID list
    
    ## Assign NA to all clumps whose IDs are found in excludeID
    Patch_scenario13[formaskSieve %in% excludeIDash13] <- NA
    
    ################ with management 1
    
    # create raster which shows just the roadside tree patches and the infected trees
    hab_scenario14 <- hab_scenario13
    while (length(hab_scenario14[hab_scenario14[] == 11]) < prob_fell * length(hab_scenario14[hab_scenario14[] == 10])) {
      hab_road14 <- hab_map
      hab_road14[hab_road14[] != 4] <- 0
      hab_road14[hab_road14[] == 4] <- 1
      Habitat_infect_road14 <-
        hab_road14    #create copy of raod trees
      Habitat_infect_road14[hab_scenario14[] == 10] <-
        2  #mark all infected road trees
      Habitat_infect_road14[hab_road14[] != 1] <-
        0   # remove and non-roadside infected trees
      Habitat_road_infected_buf14 <- Habitat_infect_road14
      Habitat_road_infected_buf14[Habitat_infect_road14[] != 2] <-
        NA   #all cells except NA become NA
      Habitat_road_infected_buf14[Habitat_infect_road14[] == 2] <-
        1  # infected cells now 1
      Habitat_road_infected_buffer14 <-
        Habitat_road_infected_buf14 #copy
      Habitat_road_infected_buffer14[Habitat_road_infected_buf14[] == 1] <-
        NA  #empty version of map
      Habitat_road_infected_buffer14[sample(which(values(Habitat_road_infected_buf14) ==
                                                    1), 1)] <-
        1 #converts random tree infected cell to one
      Habitat_road_infected_buffer14 <-
        buffer(Habitat_road_infected_buffer14, width = buffer_distance) # add a 10 m buffer
      #now remove all ash trees that fall within specified buffer and create new raster holding this treatment
      Habitat_road_infected_markedtrees14 <-
        Habitat_road_infected_buffer14
      Habitat_road_infected_markedtrees14[hab_scenario14[] == 7] <-
        2  #mark all ash trees within the buffer
      Habitat_road_infected_markedtrees14[hab_scenario14[] == 10] <-
        3
      Habitat_road_infected_markedtrees14[is.na(Habitat_road_infected_buffer14[])] <-
        0
      hab_scenario14[Habitat_road_infected_markedtrees14[] == 3] <-
        11 # 11 = felled and infected
      hab_scenario11[Habitat_road_infected_markedtrees11[] == 2] <-
        12 # 12 = felled and healthy
    }
    
    #### corresponding patch map
    
    ashcountstats14  <-
      #count of ash trees in each forest patch that are 6, 7 and 8 (ash trees)
      data.frame(zonal(
        match(hab_scenario14, c(6, 7, 8)),
        formaskSieve,
        #each forest patch is a zone
        fun = 'count',
        digits = 0,
        na.rm = TRUE
      ))
    #create replicate raster to work with
    Patch_scenario14 <- formaskSieve
    excludeIDash14 <-
      ashcountstats14$zone[which(ashcountstats14$count < 1)]  ## count of ash trees and update ID list
    
    ## Assign NA to all clumps whose IDs are found in excludeID
    Patch_scenario14[formaskSieve %in% excludeIDash14] <- NA
    
    ################ with management 2
    
    # create raster which shows just the roadside tree patches and the infected trees
    hab_scenario15 <- hab_scenario13
    while (length(hab_scenario15[hab_scenario15[] == 11]) < prob_fell2 * length(hab_scenario13[hab_scenario13[] == 10])) {
      hab_road15 <- hab_map
      hab_road15[hab_road15[] != 4] <- 0
      hab_road15[hab_road15[] == 4] <- 1
      Habitat_infect_road15 <-
        hab_road15  #create copy of raod trees
      Habitat_infect_road15[hab_scenario15[] == 10] <-
        2  #mark all infected road trees
      Habitat_infect_road15[hab_road15[] != 1] <-
        0   # remove and non-roadside infected trees
      Habitat_road_infected_buf15 <- Habitat_infect_road15
      Habitat_road_infected_buf15[Habitat_infect_road15[] != 2] <-
        NA   #all cells except NA become NA
      Habitat_road_infected_buf15[Habitat_infect_road15[] == 2] <-
        1  # infected cells now 1
      Habitat_road_infected_buffer15 <-
        Habitat_road_infected_buf15 #copy
      Habitat_road_infected_buffer15[Habitat_road_infected_buf15[] == 1] <-
        NA  #empty version of map
      Habitat_road_infected_buffer15[sample(which(values(Habitat_road_infected_buf15) ==
                                                    1), 1)] <-
        1 #converts random tree infected cell to one
      Habitat_road_infected_buffer15 <-
        buffer(Habitat_road_infected_buffer15, width = buffer_distance) # add a 10 m buffer
      #now remove all ash trees that fall within specified buffer and create new raster holding this treatment
      Habitat_road_infected_markedtrees15 <-
        Habitat_road_infected_buffer15
      Habitat_road_infected_markedtrees15[hab_scenario15[] == 7] <-
        2  #mark all ash trees within the buffer
      Habitat_road_infected_markedtrees15[hab_scenario15[] == 10] <-
        3
      Habitat_road_infected_markedtrees15[is.na(Habitat_road_infected_buffer15[])] <-
        0
      hab_scenario15[Habitat_road_infected_markedtrees15[] == 3] <-
        11 # 11 = felled and infected
      hab_scenario15[Habitat_road_infected_markedtrees15[] == 2] <-
        12 # 12 = felled and healthy
    }
    
    #### corresponding patch map
    
    ashcountstats15  <-
      #count of ash trees in each forest patch that are 6, 7 and 8 (ash trees)
      data.frame(zonal(
        match(hab_scenario15, c(6, 7, 8)),
        formaskSieve,
        #each forest patch is a zone
        fun = 'count',
        digits = 0,
        na.rm = TRUE
      ))
    #create replicate raster to work with
    Patch_scenario15 <- formaskSieve
    excludeIDash15 <-
      ashcountstats15$zone[which(ashcountstats15$count < 1)]  ## count of ash trees and update ID list
    
    ## Assign NA to all clumps whose IDs are found in excludeID
    Patch_scenario15[formaskSieve %in% excludeIDash15] <- NA
    
    
    ##### Tree disease 60 ####
    
    
    hab_scenario16 <- hab_baseline
    # replace elements of raster with dead ash trees
    hab_scenario16[hab_scenario16[] == 6] <-
      sample(
        c(6, 9),
        size = length(hab_scenario16[hab_scenario16[] ==
                                       6]),
        replace  = T,
        prob = c(1 - prob_infection6, prob_infection6)
      )
    hab_scenario16[hab_scenario16[] == 7] <-
      sample(
        c(7, 10),
        size = length(hab_scenario16[hab_scenario16[] ==
                                       7]),
        replace  = T,
        prob = c(1 - prob_infection6, prob_infection6)
      )
    hab_scenario16[hab_scenario16[] == 8] <-
      sample(
        c(8, 9),
        size = length(hab_scenario16[hab_scenario16[] ==
                                       8]),
        replace  = T,
        prob = c(1 - prob_infection6, prob_infection6)
      )
    
    
    ashcountstats16 <-
      data.frame(zonal(
        match(hab_scenario16, c(6, 7, 8)),
        formaskSieve,
        #each forest patch is a zone
        fun = 'count',
        digits = 0,
        na.rm = TRUE
      ))
    cellcount <-
      data.frame(freq(formaskSieve, useNA = "no", merge = TRUE))
    
    #create replicate raster to work with
    Patch_scenario16 <- formaskSieve
    excludeIDash16 <-
      ashcountstats16$zone[which(ashcountstats16$count < 1)]  ## count of ash trees and update ID list
    
    ## Assign NA to all clumps whose IDs are found in excludeID
    Patch_scenario16[formaskSieve %in% excludeIDash16] <- NA
    
    ################ with management 1
    
    # create raster which shows just the roadside tree patches and the infected trees
    hab_scenario17 <- hab_scenario16
    while (length(hab_scenario17[hab_scenario17[] == 11]) < prob_fell * length(hab_scenario16[hab_scenario16[] == 10])) {
      hab_road17 <- hab_map
      hab_road17[hab_road17[] != 4] <- 0
      hab_road17[hab_road17[] == 4] <- 1
      Habitat_infect_road17 <-
        hab_road17    #create copy of raod trees
      Habitat_infect_road17[hab_scenario17[] == 10] <-
        2  #mark all infected road trees
      Habitat_infect_road17[hab_road17[] != 1] <-
        0   # remove and non-roadside infected trees
      Habitat_road_infected_buf17 <- Habitat_infect_road17
      Habitat_road_infected_buf17[Habitat_infect_road17[] != 2] <-
        NA   #all cells except NA become NA
      Habitat_road_infected_buf17[Habitat_infect_road17[] == 2] <-
        1  # infected cells now 1
      Habitat_road_infected_buffer17 <-
        Habitat_road_infected_buf17 #copy
      Habitat_road_infected_buffer17[Habitat_road_infected_buf17[] == 1] <-
        NA  #empty version of map
      Habitat_road_infected_buffer17[sample(which(values(Habitat_road_infected_buf17) ==
                                                    1), 1)] <-
        1 #converts random tree infected cell to one
      Habitat_road_infected_buffer17 <-
        buffer(Habitat_road_infected_buffer17, width = buffer_distance) # add a 10 m buffer
      #now remove all ash trees that fall within specified buffer and create new raster holding this treatment
      Habitat_road_infected_markedtrees17 <-
        Habitat_road_infected_buffer17
      Habitat_road_infected_markedtrees17[hab_scenario17[] == 7] <-
        2  #mark all ash trees within the buffer
      Habitat_road_infected_markedtrees17[hab_scenario17[] == 10] <-
        3
      Habitat_road_infected_markedtrees17[is.na(Habitat_road_infected_buffer17[])] <-
        0
      hab_scenario17[Habitat_road_infected_markedtrees17[] == 3] <-
        11 # 11 = felled and infected
      hab_scenario17[Habitat_road_infected_markedtrees17[] == 2] <-
        12 # 12 = felled and healthy
    }
    
    #### corresponding patch map
    
    ashcountstats17  <-
      #count of ash trees in each forest patch that are 6, 7 and 8 (ash trees)
      data.frame(zonal(
        match(hab_scenario17, c(6, 7, 8)),
        formaskSieve,
        #each forest patch is a zone
        fun = 'count',
        digits = 0,
        na.rm = TRUE
      ))
    #create replicate raster to work with
    Patch_scenario17 <- formaskSieve
    excludeIDash17 <-
      ashcountstats17$zone[which(ashcountstats17$count < 1)]  ## count of ash trees and update ID list
    
    ## Assign NA to all clumps whose IDs are found in excludeID
    Patch_scenario17[formaskSieve %in% excludeIDash17] <- NA
    
    ################ with management 2
    
    # create raster which shows just the roadside tree patches and the infected trees
    hab_scenario18 <- hab_scenario16
    while (length(hab_scenario18[hab_scenario18[] == 11]) < prob_fell2 * length(hab_scenario16[hab_scenario16[] == 10])) {
      hab_road18 <- hab_map
      hab_road18[hab_road18[] != 4] <- 0
      hab_road18[hab_road18[] == 4] <- 1
      Habitat_infect_road18 <-
        hab_road18  #create copy of raod trees
      Habitat_infect_road18[hab_scenario18[] == 10] <-
        2  #mark all infected road trees
      Habitat_infect_road18[hab_road18[] != 1] <-
        0   # remove and non-roadside infected trees
      Habitat_road_infected_buf18 <- Habitat_infect_road18
      Habitat_road_infected_buf18[Habitat_infect_road18[] != 2] <-
        NA   #all cells except NA become NA
      Habitat_road_infected_buf18[Habitat_infect_road18[] == 2] <-
        1  # infected cells now 1
      Habitat_road_infected_buffer18 <-
        Habitat_road_infected_buf18 #copy
      Habitat_road_infected_buffer18[Habitat_road_infected_buf18[] == 1] <-
        NA  #empty version of map
      Habitat_road_infected_buffer18[sample(which(values(Habitat_road_infected_buf18) ==
                                                    1), 1)] <-
        1 #converts random tree infected cell to one
      Habitat_road_infected_buffer18 <-
        buffer(Habitat_road_infected_buffer18, width = buffer_distance) # add a 10 m buffer
      #now remove all ash trees that fall within specified buffer and create new raster holding this treatment
      Habitat_road_infected_markedtrees18 <-
        Habitat_road_infected_buffer18
      Habitat_road_infected_markedtrees18[hab_scenario18[] == 7] <-
        2  #mark all ash trees within the buffer
      Habitat_road_infected_markedtrees18[hab_scenario18[] == 10] <-
        3
      Habitat_road_infected_markedtrees18[is.na(Habitat_road_infected_buffer18[])] <-
        0
      hab_scenario18[Habitat_road_infected_markedtrees18[] == 3] <-
        11 # 11 = felled and infected
      hab_scenario18[Habitat_road_infected_markedtrees18[] == 2] <-
        12 # 12 = felled and healthy
    }
    
    #### corresponding patch map
    
    ashcountstats18  <-
      #count of ash trees in each forest patch that are 6, 7 and 8 (ash trees)
      data.frame(zonal(
        match(hab_scenario18, c(6, 7, 8)),
        formaskSieve,
        #each forest patch is a zone
        fun = 'count',
        digits = 0,
        na.rm = TRUE
      ))
    #create replicate raster to work with
    Patch_scenario18 <- formaskSieve
    excludeIDash18 <-
      ashcountstats18$zone[which(ashcountstats18$count < 1)]  ## count of ash trees and update ID list
    
    ## Assign NA to all clumps whose IDs are found in excludeID
    Patch_scenario18[formaskSieve %in% excludeIDash18] <- NA
    
    
    
    ##### Tree disease 70 ####
    
    
    hab_scenario19 <- hab_baseline
    # replace elements of raster with dead ash trees
    hab_scenario19[hab_scenario19[] == 6] <-
      sample(
        c(6, 9),
        size = length(hab_scenario19[hab_scenario19[] ==
                                       6]),
        replace  = T,
        prob = c(1 - prob_infection7, prob_infection7)
      )
    hab_scenario19[hab_scenario19[] == 7] <-
      sample(
        c(7, 10),
        size = length(hab_scenario19[hab_scenario19[] ==
                                       7]),
        replace  = T,
        prob = c(1 - prob_infection7, prob_infection7)
      )
    hab_scenario19[hab_scenario19[] == 8] <-
      sample(
        c(8, 9),
        size = length(hab_scenario19[hab_scenario19[] ==
                                       8]),
        replace  = T,
        prob = c(1 - prob_infection7, prob_infection7)
      )
    
    
    ashcountstats19 <-
      data.frame(zonal(
        match(hab_scenario19, c(6, 7, 8)),
        formaskSieve,
        #each forest patch is a zone
        fun = 'count',
        digits = 0,
        na.rm = TRUE
      ))
    cellcount <-
      data.frame(freq(formaskSieve, useNA = "no", merge = TRUE))
    
    #create replicate raster to work with
    Patch_scenario19 <- formaskSieve
    excludeIDash19 <-
      ashcountstats19$zone[which(ashcountstats19$count < 1)]  ## count of ash trees and update ID list
    
    ## Assign NA to all clumps whose IDs are found in excludeID
    Patch_scenario19[formaskSieve %in% excludeIDash19] <- NA
    
    ################ with management 1
    
    # create raster which shows just the roadside tree patches and the infected trees
    hab_scenario20 <- hab_scenario19
    while (length(hab_scenario20[hab_scenario20[] == 11]) < prob_fell * length(hab_scenario20[hab_scenario20[] == 10])) {
      hab_road20 <- hab_map
      hab_road20[hab_road20[] != 4] <- 0
      hab_road20[hab_road20[] == 4] <- 1
      Habitat_infect_road20 <-
        hab_road20    #create copy of raod trees
      Habitat_infect_road20[hab_scenario20[] == 10] <-
        2  #mark all infected road trees
      Habitat_infect_road20[hab_road20[] != 1] <-
        0   # remove and non-roadside infected trees
      Habitat_road_infected_buf20 <- Habitat_infect_road20
      Habitat_road_infected_buf20[Habitat_infect_road20[] != 2] <-
        NA   #all cells except NA become NA
      Habitat_road_infected_buf20[Habitat_infect_road20[] == 2] <-
        1  # infected cells now 1
      Habitat_road_infected_buffer20 <-
        Habitat_road_infected_buf20 #copy
      Habitat_road_infected_buffer20[Habitat_road_infected_buf20[] == 1] <-
        NA  #empty version of map
      Habitat_road_infected_buffer20[sample(which(values(Habitat_road_infected_buf20) ==
                                                    1), 1)] <-
        1 #converts random tree infected cell to one
      Habitat_road_infected_buffer20 <-
        buffer(Habitat_road_infected_buffer20, width = buffer_distance) # add a 10 m buffer
      #now remove all ash trees that fall within specified buffer and create new raster holding this treatment
      Habitat_road_infected_markedtrees20 <-
        Habitat_road_infected_buffer20
      Habitat_road_infected_markedtrees20[hab_scenario20[] == 7] <-
        2  #mark all ash trees within the buffer
      Habitat_road_infected_markedtrees20[hab_scenario20[] == 10] <-
        3
      Habitat_road_infected_markedtrees20[is.na(Habitat_road_infected_buffer20[])] <-
        0
      hab_scenario20[Habitat_road_infected_markedtrees20[] == 3] <-
        11 # 11 = felled and infected
      hab_scenario20[Habitat_road_infected_markedtrees20[] == 2] <-
        12 # 12 = felled and healthy
    }
    
    #### corresponding patch map
    
    ashcountstats20  <-
      #count of ash trees in each forest patch that are 6, 7 and 8 (ash trees)
      data.frame(zonal(
        match(hab_scenario20, c(6, 7, 8)),
        formaskSieve,
        #each forest patch is a zone
        fun = 'count',
        digits = 0,
        na.rm = TRUE
      ))
    #create replicate raster to work with
    Patch_scenario20 <- formaskSieve
    excludeIDash20 <-
      ashcountstats20$zone[which(ashcountstats20$count < 1)]  ## count of ash trees and update ID list
    
    ## Assign NA to all clumps whose IDs are found in excludeID
    Patch_scenario20[formaskSieve %in% excludeIDash20] <- NA
    
    ################ with management 2
    
    # create raster which shows just the roadside tree patches and the infected trees
    hab_scenario21 <- hab_scenario19
    while (length(hab_scenario21[hab_scenario21[] == 11]) < prob_fell2 * length(hab_scenario19[hab_scenario19[] == 10])) {
      hab_road21 <- hab_map
      hab_road21[hab_road21[] != 4] <- 0
      hab_road21[hab_road21[] == 4] <- 1
      Habitat_infect_road21 <-
        hab_road21  #create copy of raod trees
      Habitat_infect_road21[hab_scenario21[] == 10] <-
        2  #mark all infected road trees
      Habitat_infect_road21[hab_road21[] != 1] <-
        0   # remove and non-roadside infected trees
      Habitat_road_infected_buf21 <- Habitat_infect_road21
      Habitat_road_infected_buf21[Habitat_infect_road21[] != 2] <-
        NA   #all cells except NA become NA
      Habitat_road_infected_buf21[Habitat_infect_road21[] == 2] <-
        1  # infected cells now 1
      Habitat_road_infected_buffer21 <-
        Habitat_road_infected_buf21 #copy
      Habitat_road_infected_buffer21[Habitat_road_infected_buf21[] == 1] <-
        NA  #empty version of map
      Habitat_road_infected_buffer21[sample(which(values(Habitat_road_infected_buf21) ==
                                                    1), 1)] <-
        1 #converts random tree infected cell to one
      Habitat_road_infected_buffer21 <-
        buffer(Habitat_road_infected_buffer21, width = buffer_distance) # add a 10 m buffer
      #now remove all ash trees that fall within specified buffer and create new raster holding this treatment
      Habitat_road_infected_markedtrees21 <-
        Habitat_road_infected_buffer21
      Habitat_road_infected_markedtrees21[hab_scenario21[] == 7] <-
        2  #mark all ash trees within the buffer
      Habitat_road_infected_markedtrees21[hab_scenario21[] == 10] <-
        3
      Habitat_road_infected_markedtrees21[is.na(Habitat_road_infected_buffer21[])] <-
        0
      hab_scenario21[Habitat_road_infected_markedtrees21[] == 3] <-
        11 # 11 = felled and infected
      hab_scenario21[Habitat_road_infected_markedtrees21[] == 2] <-
        12 # 12 = felled and healthy
    }
    
    #### corresponding patch map
    
    ashcountstats21  <-
      #count of ash trees in each forest patch that are 6, 7 and 8 (ash trees)
      data.frame(zonal(
        match(hab_scenario21, c(6, 7, 8)),
        formaskSieve,
        #each forest patch is a zone
        fun = 'count',
        digits = 0,
        na.rm = TRUE
      ))
    #create replicate raster to work with
    Patch_scenario21 <- formaskSieve
    excludeIDash21 <-
      ashcountstats21$zone[which(ashcountstats21$count < 1)]  ## count of ash trees and update ID list
    
    ## Assign NA to all clumps whose IDs are found in excludeID
    Patch_scenario21[formaskSieve %in% excludeIDash21] <- NA
    
    
    
    ##### Tree disease 80 ####
    
    
    hab_scenario22 <- hab_baseline
    # replace elements of raster with dead ash trees
    hab_scenario22[hab_scenario22[] == 6] <-
      sample(
        c(6, 9),
        size = length(hab_scenario22[hab_scenario22[] ==
                                       6]),
        replace  = T,
        prob = c(1 - prob_infection8, prob_infection8)
      )
    hab_scenario22[hab_scenario22[] == 7] <-
      sample(
        c(7, 10),
        size = length(hab_scenario22[hab_scenario22[] ==
                                       7]),
        replace  = T,
        prob = c(1 - prob_infection8, prob_infection8)
      )
    hab_scenario22[hab_scenario22[] == 8] <-
      sample(
        c(8, 9),
        size = length(hab_scenario22[hab_scenario22[] ==
                                       8]),
        replace  = T,
        prob = c(1 - prob_infection8, prob_infection8)
      )
    
    
    ashcountstats22 <-
      data.frame(zonal(
        match(hab_scenario22, c(6, 7, 8)),
        formaskSieve,
        #each forest patch is a zone
        fun = 'count',
        digits = 0,
        na.rm = TRUE
      ))
    cellcount <-
      data.frame(freq(formaskSieve, useNA = "no", merge = TRUE))
    
    #create replicate raster to work with
    Patch_scenario22 <- formaskSieve
    excludeIDash22 <-
      ashcountstats22$zone[which(ashcountstats22$count < 1)]  ## count of ash trees and update ID list
    
    ## Assign NA to all clumps whose IDs are found in excludeID
    Patch_scenario22[formaskSieve %in% excludeIDash22] <- NA
    
    ################ with management 1
    
    # create raster which shows just the roadside tree patches and the infected trees
    hab_scenario23 <- hab_scenario22
    while (length(hab_scenario23[hab_scenario23[] == 11]) < prob_fell * length(hab_scenario23[hab_scenario23[] == 10])) {
      hab_road23 <- hab_map
      hab_road23[hab_road23[] != 4] <- 0
      hab_road23[hab_road23[] == 4] <- 1
      Habitat_infect_road23 <-
        hab_road23    #create copy of raod trees
      Habitat_infect_road23[hab_scenario23[] == 10] <-
        2  #mark all infected road trees
      Habitat_infect_road23[hab_road23[] != 1] <-
        0   # remove and non-roadside infected trees
      Habitat_road_infected_buf23 <- Habitat_infect_road23
      Habitat_road_infected_buf23[Habitat_infect_road23[] != 2] <-
        NA   #all cells except NA become NA
      Habitat_road_infected_buf23[Habitat_infect_road23[] == 2] <-
        1  # infected cells now 1
      Habitat_road_infected_buffer23 <-
        Habitat_road_infected_buf23 #copy
      Habitat_road_infected_buffer23[Habitat_road_infected_buf23[] == 1] <-
        NA  #empty version of map
      Habitat_road_infected_buffer23[sample(which(values(Habitat_road_infected_buf23) ==
                                                    1), 1)] <-
        1 #converts random tree infected cell to one
      Habitat_road_infected_buffer23 <-
        buffer(Habitat_road_infected_buffer23, width = buffer_distance) # add a 10 m buffer
      #now remove all ash trees that fall within specified buffer and create new raster holding this treatment
      Habitat_road_infected_markedtrees23 <-
        Habitat_road_infected_buffer23
      Habitat_road_infected_markedtrees23[hab_scenario23[] == 7] <-
        2  #mark all ash trees within the buffer
      Habitat_road_infected_markedtrees23[hab_scenario23[] == 10] <-
        3
      Habitat_road_infected_markedtrees23[is.na(Habitat_road_infected_buffer23[])] <-
        0
      hab_scenario23[Habitat_road_infected_markedtrees23[] == 3] <-
        11 # 11 = felled and infected
      hab_scenario23[Habitat_road_infected_markedtrees23[] == 2] <-
        12 # 12 = felled and healthy
    }
    
    #### corresponding patch map
    
    ashcountstats23  <-
      #count of ash trees in each forest patch that are 6, 7 and 8 (ash trees)
      data.frame(zonal(
        match(hab_scenario23, c(6, 7, 8)),
        formaskSieve,
        #each forest patch is a zone
        fun = 'count',
        digits = 0,
        na.rm = TRUE
      ))
    #create replicate raster to work with
    Patch_scenario23 <- formaskSieve
    excludeIDash23 <-
      ashcountstats23$zone[which(ashcountstats23$count < 1)]  ## count of ash trees and update ID list
    
    ## Assign NA to all clumps whose IDs are found in excludeID
    Patch_scenario23[formaskSieve %in% excludeIDash23] <- NA
    
    ################ with management 2
    
    # create raster which shows just the roadside tree patches and the infected trees
    hab_scenario24 <- hab_scenario22
    while (length(hab_scenario24[hab_scenario24[] == 11]) < prob_fell2 * length(hab_scenario22[hab_scenario22[] == 10])) {
      hab_road24 <- hab_map
      hab_road24[hab_road24[] != 4] <- 0
      hab_road24[hab_road24[] == 4] <- 1
      Habitat_infect_road24 <-
        hab_road24  #create copy of raod trees
      Habitat_infect_road24[hab_scenario24[] == 10] <-
        2  #mark all infected road trees
      Habitat_infect_road24[hab_road24[] != 1] <-
        0   # remove and non-roadside infected trees
      Habitat_road_infected_buf24 <- Habitat_infect_road24
      Habitat_road_infected_buf24[Habitat_infect_road24[] != 2] <-
        NA   #all cells except NA become NA
      Habitat_road_infected_buf24[Habitat_infect_road24[] == 2] <-
        1  # infected cells now 1
      Habitat_road_infected_buffer24 <-
        Habitat_road_infected_buf24 #copy
      Habitat_road_infected_buffer24[Habitat_road_infected_buf24[] == 1] <-
        NA  #empty version of map
      Habitat_road_infected_buffer24[sample(which(values(Habitat_road_infected_buf24) ==
                                                    1), 1)] <-
        1 #converts random tree infected cell to one
      Habitat_road_infected_buffer24 <-
        buffer(Habitat_road_infected_buffer24, width = buffer_distance) # add a 10 m buffer
      #now remove all ash trees that fall within specified buffer and create new raster holding this treatment
      Habitat_road_infected_markedtrees24 <-
        Habitat_road_infected_buffer24
      Habitat_road_infected_markedtrees24[hab_scenario24[] == 7] <-
        2  #mark all ash trees within the buffer
      Habitat_road_infected_markedtrees24[hab_scenario24[] == 10] <-
        3
      Habitat_road_infected_markedtrees24[is.na(Habitat_road_infected_buffer24[])] <-
        0
      hab_scenario24[Habitat_road_infected_markedtrees24[] == 3] <-
        11 # 11 = felled and infected
      hab_scenario24[Habitat_road_infected_markedtrees24[] == 2] <-
        12 # 12 = felled and healthy
    }
    
    #### corresponding patch map
    
    ashcountstats24  <-
      #count of ash trees in each forest patch that are 6, 7 and 8 (ash trees)
      data.frame(zonal(
        match(hab_scenario24, c(6, 7, 8)),
        formaskSieve,
        #each forest patch is a zone
        fun = 'count',
        digits = 0,
        na.rm = TRUE
      ))
    #create replicate raster to work with
    Patch_scenario24 <- formaskSieve
    excludeIDash24 <-
      ashcountstats24$zone[which(ashcountstats24$count < 1)]  ## count of ash trees and update ID list
    
    ## Assign NA to all clumps whose IDs are found in excludeID
    Patch_scenario24[formaskSieve %in% excludeIDash24] <- NA
    
    
    
    ############ removing excess categories
    hab_scenario1[hab_scenario1[]==10] <- 9
    hab_scenario2[hab_scenario2[]==10] <- 9
    hab_scenario3[hab_scenario3[]==10] <- 9
    hab_scenario4[hab_scenario4[]==10] <- 9
    hab_scenario5[hab_scenario5[]==10] <- 9
    hab_scenario6[hab_scenario6[]==10] <- 9
    hab_scenario7[hab_scenario7[]==10] <- 9
    hab_scenario8[hab_scenario8[]==10] <- 9
    hab_scenario9[hab_scenario9[]==10] <- 9
    hab_scenario10[hab_scenario10[]==10] <- 9
    hab_scenario11[hab_scenario11[]==10] <- 9
    hab_scenario12[hab_scenario12[]==10] <- 9
    hab_scenario13[hab_scenario13[]==10] <- 9
    hab_scenario14[hab_scenario14[]==10] <- 9
    hab_scenario15[hab_scenario15[]==10] <- 9
    hab_scenario16[hab_scenario16[]==10] <- 9
    hab_scenario17[hab_scenario17[]==10] <- 9
    hab_scenario18[hab_scenario18[]==10] <- 9
    hab_scenario19[hab_scenario19[]==10] <- 9
    hab_scenario20[hab_scenario20[]==10] <- 9
    hab_scenario21[hab_scenario21[]==10] <- 9
    hab_scenario22[hab_scenario22[]==10] <- 9
    hab_scenario23[hab_scenario23[]==10] <- 9
    hab_scenario24[hab_scenario24[]==10] <- 9
    
    hab_scenario1[hab_scenario1[]==11] <- 1
    hab_scenario2[hab_scenario2[]==11] <- 1
    hab_scenario3[hab_scenario3[]==11] <- 1
    hab_scenario4[hab_scenario4[]==11] <- 1
    hab_scenario5[hab_scenario5[]==11] <- 1
    hab_scenario6[hab_scenario6[]==11] <- 1
    hab_scenario7[hab_scenario7[]==11] <- 1
    hab_scenario8[hab_scenario8[]==11] <- 1
    hab_scenario9[hab_scenario9[]==11] <- 1
    hab_scenario10[hab_scenario10[]==11] <- 1
    hab_scenario11[hab_scenario11[]==11] <- 1
    hab_scenario12[hab_scenario12[]==11] <- 1
    hab_scenario13[hab_scenario13[]==11] <- 1
    hab_scenario14[hab_scenario14[]==11] <- 1
    hab_scenario15[hab_scenario15[]==11] <- 1
    hab_scenario16[hab_scenario16[]==11] <- 1
    hab_scenario17[hab_scenario17[]==11] <- 1
    hab_scenario18[hab_scenario18[]==11] <- 1
    hab_scenario11[hab_scenario19[]==11] <- 1
    hab_scenario20[hab_scenario20[]==11] <- 1
    hab_scenario21[hab_scenario21[]==11] <- 1
    hab_scenario22[hab_scenario22[]==11] <- 1
    hab_scenario23[hab_scenario23[]==11] <- 1
    hab_scenario24[hab_scenario24[]==11] <- 1
    hab_scenario1[hab_scenario1[]==12] <- 1
    hab_scenario2[hab_scenario2[]==12] <- 1
    hab_scenario3[hab_scenario3[]==12] <- 1
    hab_scenario4[hab_scenario4[]==12] <- 1
    hab_scenario5[hab_scenario5[]==12] <- 1
    hab_scenario6[hab_scenario6[]==12] <- 1
    hab_scenario7[hab_scenario7[]==12] <- 1
    hab_scenario8[hab_scenario8[]==12] <- 1
    hab_scenario9[hab_scenario9[]==12] <- 1
    hab_scenario10[hab_scenario10[]==12] <- 1
    hab_scenario11[hab_scenario11[]==12] <- 1
    hab_scenario12[hab_scenario12[]==12] <- 1
    hab_scenario13[hab_scenario13[]==12] <- 1
    hab_scenario14[hab_scenario14[]==12] <- 1
    hab_scenario15[hab_scenario15[]==12] <- 1
    hab_scenario16[hab_scenario16[]==12] <- 1
    hab_scenario17[hab_scenario17[]==12] <- 1
    hab_scenario18[hab_scenario18[]==12] <- 1
    hab_scenario11[hab_scenario19[]==12] <- 1
    hab_scenario20[hab_scenario20[]==12] <- 1
    hab_scenario21[hab_scenario21[]==12] <- 1
    hab_scenario22[hab_scenario22[]==12] <- 1
    hab_scenario23[hab_scenario23[]==12] <- 1
    hab_scenario24[hab_scenario24[]==12] <- 1
    
    Patch_scenario0[is.na(Patch_scenario0[])] <- 0 
    Patch_scenario1[is.na(Patch_scenario1[])] <- 0 
    Patch_scenario2[is.na(Patch_scenario2[])] <- 0 
    Patch_scenario3[is.na(Patch_scenario3[])] <- 0 
    Patch_scenario4[is.na(Patch_scenario4[])] <- 0 
    Patch_scenario5[is.na(Patch_scenario5[])] <- 0 
    Patch_scenario6[is.na(Patch_scenario6[])] <- 0 
    Patch_scenario7[is.na(Patch_scenario7[])] <- 0 
    Patch_scenario8[is.na(Patch_scenario8[])] <- 0 
    Patch_scenario9[is.na(Patch_scenario9[])] <- 0 
    Patch_scenario10[is.na(Patch_scenario10[])] <- 0 
    Patch_scenario11[is.na(Patch_scenario11[])] <- 0 
    Patch_scenario12[is.na(Patch_scenario12[])] <- 0 
    Patch_scenario13[is.na(Patch_scenario13[])] <- 0 
    Patch_scenario14[is.na(Patch_scenario14[])] <- 0 
    Patch_scenario15[is.na(Patch_scenario15[])] <- 0 
    Patch_scenario16[is.na(Patch_scenario16[])] <- 0 
    Patch_scenario17[is.na(Patch_scenario17[])] <- 0 
    Patch_scenario18[is.na(Patch_scenario18[])] <- 0 
    Patch_scenario19[is.na(Patch_scenario19[])] <- 0 
    Patch_scenario20[is.na(Patch_scenario20[])] <- 0 
    Patch_scenario21[is.na(Patch_scenario21[])] <- 0 
    Patch_scenario22[is.na(Patch_scenario22[])] <- 0 
    Patch_scenario23[is.na(Patch_scenario23[])] <- 0 
    Patch_scenario24[is.na(Patch_scenario24[])] <- 0 
    
    Tree <-
      brick(
        hab_scenario0,
        Patch_scenario0,
        hab_scenario1,
        Patch_scenario1,
        hab_scenario2,
        Patch_scenario2,
        hab_scenario3,
        Patch_scenario3,
        hab_scenario4,
        Patch_scenario4,
        hab_scenario5,
        Patch_scenario5,
        hab_scenario6,
        Patch_scenario6,
        hab_scenario7,
        Patch_scenario7,
        hab_scenario8,
        Patch_scenario8,
        hab_scenario9,
        Patch_scenario9,
        hab_scenario10,
        Patch_scenario10,
        hab_scenario11,
        Patch_scenario11,
        hab_scenario12,
        Patch_scenario12,
        hab_scenario13,
        Patch_scenario13,
        hab_scenario14,
        Patch_scenario14,
        hab_scenario15,
        Patch_scenario15,
        hab_scenario16,
        Patch_scenario16,
        hab_scenario17,
        Patch_scenario17,
        hab_scenario18,
        Patch_scenario18,
        hab_scenario19,
        Patch_scenario19,
        hab_scenario20,
        Patch_scenario20,
        hab_scenario21,
        Patch_scenario21,
        hab_scenario22,
        Patch_scenario22,
        hab_scenario23,
        Patch_scenario23,
        hab_scenario24,
        Patch_scenario24
      )
    
    return(Tree)
    
  }

setwd("C:/Baseline_maps") 
habitat_raster <- raster("Hab_TM17SW.tif")
hab_baseline1 <- AssigningAsh(hab_map = habitat_raster, prob_ash = .13)
hab_baseline2 <- AssigningAsh(hab_map = habitat_raster, prob_ash = .13)
hab_baseline3 <- AssigningAsh(hab_map = habitat_raster, prob_ash = .13)

### L1_R1 ####
Landrep1_R1 <- killingAsh(
  hab_map = habitat_raster,
  hab_baseline = hab_baseline1,
  prob_infection1 = .1,
  prob_infection2 = .2,
  prob_infection3 = .3,
  prob_infection4 = .4,
  prob_infection5 = .5,
  prob_infection6 = .6,
  prob_infection7 = .7,
  prob_infection8 = .8,
  prob_fell = .4,
  prob_fell2 = .8,
  buffer_distance = 50
)
L1_R1_Habitat_Scenario0 <- subset(Landrep1_R1, 1)
L1_R1_Patch_Scenario0 <- subset(Landrep1_R1, 2)
L1_R1_Habitat_Scenario1 <- subset(Landrep1_R1, 3)
L1_R1_Patch_Scenario1 <- subset(Landrep1_R1, 4)
L1_R1_Habitat_Scenario2 <- subset(Landrep1_R1, 5)
L1_R1_Patch_Scenario2 <- subset(Landrep1_R1, 6)
L1_R1_Habitat_Scenario3 <- subset(Landrep1_R1, 7)
L1_R1_Patch_Scenario3 <- subset(Landrep1_R1, 8)
L1_R1_Habitat_Scenario4 <- subset(Landrep1_R1, 9)
L1_R1_Patch_Scenario4 <- subset(Landrep1_R1, 10)
L1_R1_Habitat_Scenario5 <- subset(Landrep1_R1, 11)
L1_R1_Patch_Scenario5 <- subset(Landrep1_R1, 12)
L1_R1_Habitat_Scenario6 <- subset(Landrep1_R1, 13)
L1_R1_Patch_Scenario6 <- subset(Landrep1_R1, 14)
L1_R1_Habitat_Scenario7 <- subset(Landrep1_R1, 15)
L1_R1_Patch_Scenario7 <- subset(Landrep1_R1, 16)
L1_R1_Habitat_Scenario8 <- subset(Landrep1_R1, 17)
L1_R1_Patch_Scenario8 <- subset(Landrep1_R1, 18)
L1_R1_Habitat_Scenario9 <- subset(Landrep1_R1, 19)
L1_R1_Patch_Scenario9 <- subset(Landrep1_R1, 20)
L1_R1_Habitat_Scenario10 <- subset(Landrep1_R1, 21)
L1_R1_Patch_Scenario10 <- subset(Landrep1_R1, 22)
L1_R1_Habitat_Scenario11 <- subset(Landrep1_R1, 23)
L1_R1_Patch_Scenario11 <- subset(Landrep1_R1, 24)
L1_R1_Habitat_Scenario12 <- subset(Landrep1_R1, 25)
L1_R1_Patch_Scenario12 <- subset(Landrep1_R1, 26)
L1_R1_Habitat_Scenario13 <- subset(Landrep1_R1, 27)
L1_R1_Patch_Scenario13 <- subset(Landrep1_R1, 28)
L1_R1_Habitat_Scenario14 <- subset(Landrep1_R1, 29)
L1_R1_Patch_Scenario14 <- subset(Landrep1_R1, 30)
L1_R1_Habitat_Scenario15 <- subset(Landrep1_R1, 31)
L1_R1_Patch_Scenario15 <- subset(Landrep1_R1, 32)
L1_R1_Habitat_Scenario16 <- subset(Landrep1_R1, 33)
L1_R1_Patch_Scenario16 <- subset(Landrep1_R1, 34)
L1_R1_Habitat_Scenario17 <- subset(Landrep1_R1, 35)
L1_R1_Patch_Scenario17 <- subset(Landrep1_R1, 36)
L1_R1_Habitat_Scenario18 <- subset(Landrep1_R1, 37)
L1_R1_Patch_Scenario18 <- subset(Landrep1_R1, 38)
L1_R1_Habitat_Scenario19 <- subset(Landrep1_R1, 39)
L1_R1_Patch_Scenario19 <- subset(Landrep1_R1, 40)
L1_R1_Habitat_Scenario20 <- subset(Landrep1_R1, 41)
L1_R1_Patch_Scenario20 <- subset(Landrep1_R1, 42)
L1_R1_Habitat_Scenario21 <- subset(Landrep1_R1, 43)
L1_R1_Patch_Scenario21 <- subset(Landrep1_R1, 44)
L1_R1_Habitat_Scenario22 <- subset(Landrep1_R1, 45)
L1_R1_Patch_Scenario22 <- subset(Landrep1_R1, 46)
L1_R1_Habitat_Scenario23 <- subset(Landrep1_R1, 47)
L1_R1_Patch_Scenario23 <- subset(Landrep1_R1, 48)
L1_R1_Habitat_Scenario24 <- subset(Landrep1_R1, 49)
L1_R1_Patch_Scenario24 <- subset(Landrep1_R1, 50)


#### exporting maps as ascii files ###

writeRaster(L1_R1_Habitat_Scenario0, "C:/FionaPlenderleith/Batch1001/Inputs/L1_R1_Habitat_Scenario0.asc", NAflag= -9999, format="ascii", overwrite=TRUE)
writeRaster(L1_R1_Patch_Scenario0, "C:/FionaPlenderleith/Batch1001/Inputs/L1_R1_Patch_Scenario0.asc", NAflag= -9999, format="ascii", overwrite=TRUE)
writeRaster(L1_R1_Habitat_Scenario1, "C:/FionaPlenderleith/Batch1001/Inputs/L1_R1_Habitat_Scenario1.asc", NAflag= -9999,  format="ascii", overwrite=TRUE)
writeRaster(L1_R1_Patch_Scenario1, "C:/FionaPlenderleith/Batch1001/Inputs/L1_R1_Patch_Scenario1.asc", NAflag= -9999,  format="ascii", overwrite=TRUE)
writeRaster(L1_R1_Habitat_Scenario2, "C:/FionaPlenderleith/Batch1001/Inputs/L1_R1_Habitat_Scenario2.asc", NAflag= -9999,  format="ascii", overwrite=TRUE)
writeRaster(L1_R1_Patch_Scenario2, "C:/FionaPlenderleith/Batch1001/Inputs/L1_R1_Patch_Scenario2.asc", NAflag= -9999, format="ascii", overwrite=TRUE)
writeRaster(L1_R1_Habitat_Scenario3, "C:/FionaPlenderleith/Batch1001/Inputs/L1_R1_Habitat_Scenario3.asc", NAflag= -9999, format="ascii", overwrite=TRUE)
writeRaster(L1_R1_Patch_Scenario3, "C:/FionaPlenderleith/Batch1001/Inputs/L1_R1_Patch_Scenario3.asc", NAflag= -9999, format="ascii", overwrite=TRUE)
writeRaster(L1_R1_Habitat_Scenario4, "C:/FionaPlenderleith/Batch1001/Inputs/L1_R1_Habitat_Scenario4.asc", NAflag= -9999, format="ascii",   overwrite=TRUE)
writeRaster(L1_R1_Patch_Scenario4, "C:/FionaPlenderleith/Batch1001/Inputs/L1_R1_Patch_Scenario4.asc", NAflag= -9999, format="ascii",   overwrite=TRUE)
writeRaster(L1_R1_Habitat_Scenario5, "C:/FionaPlenderleith/Batch1001/Inputs/L1_R1_Habitat_Scenario5.asc", NAflag= -9999, format="ascii",   overwrite=TRUE)
writeRaster(L1_R1_Patch_Scenario5, "C:/FionaPlenderleith/Batch1001/Inputs/L1_R1_Patch_Scenario5.asc", NAflag= -9999, format="ascii",   overwrite=TRUE)
writeRaster(L1_R1_Habitat_Scenario6, "C:/FionaPlenderleith/Batch1001/Inputs/L1_R1_Habitat_Scenario6.asc", NAflag= -9999, format="ascii",   overwrite=TRUE)
writeRaster(L1_R1_Patch_Scenario6, "C:/FionaPlenderleith/Batch1001/Inputs/L1_R1_Patch_Scenario6.asc", NAflag= -9999, format="ascii",   overwrite=TRUE)
writeRaster(L1_R1_Habitat_Scenario7, "C:/FionaPlenderleith/Batch1001/Inputs/L1_R1_Habitat_Scenario7.asc", NAflag= -9999, format="ascii",   overwrite=TRUE)
writeRaster(L1_R1_Patch_Scenario7, "C:/FionaPlenderleith/Batch1001/Inputs/L1_R1_Patch_Scenario7.asc", NAflag= -9999, format="ascii",   overwrite=TRUE)
writeRaster(L1_R1_Habitat_Scenario8, "C:/FionaPlenderleith/Batch1001/Inputs/L1_R1_Habitat_Scenario8.asc", NAflag= -9999, format="ascii",   overwrite=TRUE)
writeRaster(L1_R1_Patch_Scenario8, "C:/FionaPlenderleith/Batch1001/Inputs/L1_R1_Patch_Scenario8.asc", NAflag= -9999, format="ascii",   overwrite=TRUE)
writeRaster(L1_R1_Habitat_Scenario9, "C:/FionaPlenderleith/Batch1001/Inputs/L1_R1_Habitat_Scenario9.asc", NAflag= -9999, format="ascii",   overwrite=TRUE)
writeRaster(L1_R1_Patch_Scenario9, "C:/FionaPlenderleith/Batch1001/Inputs/L1_R1_Patch_Scenario9.asc", NAflag= -9999, format="ascii",   overwrite=TRUE)
writeRaster(L1_R1_Habitat_Scenario10, "C:/FionaPlenderleith/Batch1001/Inputs/L1_R1_Habitat_Scenario10.asc", NAflag= -9999, format="ascii",   overwrite=TRUE)
writeRaster(L1_R1_Patch_Scenario10, "C:/FionaPlenderleith/Batch1001/Inputs/L1_R1_Patch_Scenario10.asc", NAflag= -9999, format="ascii", overwrite=TRUE)
writeRaster(L1_R1_Habitat_Scenario11, "C:/FionaPlenderleith/Batch1001/Inputs/L1_R1_Habitat_Scenario11.asc", NAflag= -9999, format="ascii", overwrite=TRUE)
writeRaster(L1_R1_Patch_Scenario11, "C:/FionaPlenderleith/Batch1001/Inputs/L1_R1_Patch_Scenario11.asc", NAflag= -9999, format="ascii", overwrite=TRUE)
writeRaster(L1_R1_Habitat_Scenario12, "C:/FionaPlenderleith/Batch1001/Inputs/L1_R1_Habitat_Scenario12.asc", NAflag= -9999, format="ascii", overwrite=TRUE)
writeRaster(L1_R1_Patch_Scenario12, "C:/FionaPlenderleith/Batch1001/Inputs/L1_R1_Patch_Scenario12.asc", NAflag= -9999, format="ascii",  overwrite=TRUE)
writeRaster(L1_R1_Habitat_Scenario13, "C:/FionaPlenderleith/Batch1001/Inputs/L1_R1_Habitat_Scenario13.asc", NAflag= -9999, format="ascii",  overwrite=TRUE)
writeRaster(L1_R1_Patch_Scenario13, "C:/FionaPlenderleith/Batch1001/Inputs/L1_R1_Patch_Scenario13.asc", NAflag= -9999, format="ascii", overwrite=TRUE)
writeRaster(L1_R1_Habitat_Scenario14, "C:/FionaPlenderleith/Batch1001/Inputs/L1_R1_Habitat_Scenario14.asc", NAflag= -9999, format="ascii",   overwrite=TRUE)
writeRaster(L1_R1_Patch_Scenario14, "C:/FionaPlenderleith/Batch1001/Inputs/L1_R1_Patch_Scenario14.asc", NAflag= -9999, format="ascii",   overwrite=TRUE)
writeRaster(L1_R1_Habitat_Scenario15, "C:/FionaPlenderleith/Batch1001/Inputs/L1_R1_Habitat_Scenario15.asc", NAflag= -9999, format="ascii",   overwrite=TRUE)
writeRaster(L1_R1_Patch_Scenario15, "C:/FionaPlenderleith/Batch1001/Inputs/L1_R1_Patch_Scenario15.asc", NAflag= -9999, format="ascii", overwrite=TRUE)
writeRaster(L1_R1_Habitat_Scenario16, "C:/FionaPlenderleith/Batch1001/Inputs/L1_R1_Habitat_Scenario16.asc", NAflag= -9999, format="ascii", overwrite=TRUE)
writeRaster(L1_R1_Patch_Scenario16, "C:/FionaPlenderleith/Batch1001/Inputs/L1_R1_Patch_Scenario16.asc", NAflag= -9999, format="ascii", overwrite=TRUE)
writeRaster(L1_R1_Habitat_Scenario17, "C:/FionaPlenderleith/Batch1001/Inputs/L1_R1_Habitat_Scenario17.asc", NAflag= -9999, format="ascii", overwrite=TRUE)
writeRaster(L1_R1_Patch_Scenario17, "C:/FionaPlenderleith/Batch1001/Inputs/L1_R1_Patch_Scenario17.asc", NAflag= -9999, format="ascii",   overwrite=TRUE)
writeRaster(L1_R1_Habitat_Scenario18, "C:/FionaPlenderleith/Batch1001/Inputs/L1_R1_Habitat_Scenario18.asc", NAflag= -9999, format="ascii", overwrite=TRUE)
writeRaster(L1_R1_Patch_Scenario18, "C:/FionaPlenderleith/Batch1001/Inputs/L1_R1_Patch_Scenario18.asc", NAflag= -9999, format="ascii", overwrite=TRUE)
writeRaster(L1_R1_Habitat_Scenario19, "C:/FionaPlenderleith/Batch1001/Inputs/L1_R1_Habitat_Scenario19.asc", NAflag= -9999, format="ascii", overwrite=TRUE)
writeRaster(L1_R1_Patch_Scenario19, "C:/FionaPlenderleith/Batch1001/Inputs/L1_R1_Patch_Scenario19.asc", NAflag= -9999, format="ascii", overwrite=TRUE)
writeRaster(L1_R1_Habitat_Scenario20, "C:/FionaPlenderleith/Batch1001/Inputs/L1_R1_Habitat_Scenario20.asc", NAflag= -9999, format="ascii", overwrite=TRUE)
writeRaster(L1_R1_Patch_Scenario20, "C:/FionaPlenderleith/Batch1001/Inputs/L1_R1_Patch_Scenario20.asc", NAflag= -9999, format="ascii", overwrite=TRUE)
writeRaster(L1_R1_Habitat_Scenario21, "C:/FionaPlenderleith/Batch1001/Inputs/L1_R1_Habitat_Scenario21.asc", NAflag= -9999, format="ascii", overwrite=TRUE)
writeRaster(L1_R1_Patch_Scenario21, "C:/FionaPlenderleith/Batch1001/Inputs/L1_R1_Patch_Scenario21.asc", NAflag= -9999, format="ascii", overwrite=TRUE)
writeRaster(L1_R1_Habitat_Scenario22, "C:/FionaPlenderleith/Batch1001/Inputs/L1_R1_Habitat_Scenario22.asc", NAflag= -9999, format="ascii", overwrite=TRUE)
writeRaster(L1_R1_Patch_Scenario22, "C:/FionaPlenderleith/Batch1001/Inputs/L1_R1_Patch_Scenario22.asc", NAflag= -9999, format="ascii", overwrite=TRUE)
writeRaster(L1_R1_Habitat_Scenario23, "C:/FionaPlenderleith/Batch1001/Inputs/L1_R1_Habitat_Scenario23.asc", NAflag= -9999, format="ascii", overwrite=TRUE)
writeRaster(L1_R1_Patch_Scenario23, "C:/FionaPlenderleith/Batch1001/Inputs/L1_R1_Patch_Scenario23.asc", NAflag= -9999, format="ascii", overwrite=TRUE)
writeRaster(L1_R1_Habitat_Scenario24, "C:/FionaPlenderleith/Batch1001/Inputs/L1_R1_Habitat_Scenario24.asc", NAflag= -9999, format="ascii", overwrite=TRUE)
writeRaster(L1_R1_Patch_Scenario24, "C:/FionaPlenderleith/Batch1001/Inputs/L1_R1_Patch_Scenario24.asc", NAflag= -9999, format="ascii", overwrite=TRUE)

patch_frequency_Batch1001 <- as.data.frame(freq(L1_R1_Patch_Scenario0))
patch_frequency_Batch1001 <- subset(patch_frequency_Batch1001, patch_frequency_Batch1001$count>400)
patch_frequency_Batch1001 <- subset(patch_frequency_Batch1001, patch_frequency_Batch1001$value != 0)
patch_frequency_Batch1001 <- sample_n(patch_frequency_Batch1001, 15)
patch_frequency_Batch1001 <- patch_frequency_Batch1001[order(patch_frequency_Batch1001$value),]
patches_Batch1001 <- unname(unlist(patch_frequency_Batch1001["value"]))
write.table(patches_Batch1001, sep = " ", row.names = F, col.names = F, "C:/FionaPlenderleith/Batch1001/Inputs/patches.txt")



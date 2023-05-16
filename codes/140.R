#miR-140-5p
#microRNA-140-5p
miR1405p_Targets <- getTargets(microRNA = "miR-140-5p")
Day1Targets140 <- confidentTargets(DEres1 = hpc1_199a1, DEres2 = hpc1_199b1, 
                                   targets1 = miR1405p_Targets,
                                   targets2 = miR1405p_Targets)
Day0Targets140 <- confidentTargets(DEres1 = hpc0_199a0, DEres2 = hpc0_199b0, 
                                   targets1 = miR1405p_Targets,
                                   targets2 = miR1405p_Targets)
list140 <- listOfGenes(D0 = Day0Targets140, D1 = Day1Targets140)
MSC_miR1405p_targets_199a <- getLogfc(D0 = hpc0_199a0, D1 = hpc1_199a1, 
                                      Targets = list140)
MSC_miR1405p_targets_199b <- getLogfc(D0 = hpc0_199b0, D1 = hpc1_199b1, 
                                      Targets = list140)
colnames(MSC_miR1405p_targets_199a) <- c("24H", "48H")
colnames(MSC_miR1405p_targets_199b) <- c("24H", "48H")
a_140 <- createHeatMap(targets = MSC_miR1405p_targets_199a, microRNA = "miR1405p")
b_140 <- createHeatMap(targets = MSC_miR1405p_targets_199b, microRNA = "miR1405p")
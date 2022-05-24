if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("ballgown")

library(ballgown)
bg_hed = ballgown(dataDir="/scratch/henrygeo/Results/stringtie/", samplePattern='\\w{2}\\d+', meas='all')
saveRDS(bg_hed, file='/scratch/henrygeo/Results/ballgown/bg_hed.rds')
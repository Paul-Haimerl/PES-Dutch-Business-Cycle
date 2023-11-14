setwd("C:/Users/phaim/OneDrive/PES Research Div/Dutch-Business-Cycle")
source("R/0_Init.R")
plotPath <- "C:/Users/phaim/OneDrive/PES Research Div/Figures"
outputPath <- paste0(getwd(), "/Output/")
InitProject(outputPath = outputPath)

nRandom <- 1e3
UC_result <- UC_Wrapper(y = y, det.type = "trend", corr.ind = FALSE, d = 1, nRandom = nRandom, outputPath = outputPath, fileName = "UC")
saveRDS(object = UC_result, file = "UC_result.rds")
UC_Corr_result <- UC_Wrapper(y = y, det.type = "trend", corr.ind = TRUE, d = 1, nRandom = nRandom, outputPath = outputPath, fileName = "UC_Corr")
saveRDS(object = UC_Corr_result, file = "UC_Corr_result.rds")
UC_Frac_result <- UC_Wrapper(y = y, det.type = "frac", corr.ind = TRUE, d = FALSE, nRandom = nRandom, outputPath = outputPath, fileName = "UC_Frac")
saveRDS(object = UC_Frac_result, file = "UC_Frac_result.rds")

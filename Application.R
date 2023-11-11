nRandom <- 100
UC_result <- UC_Wrapper(y = y, det.type = "trend", corr.ind = FALSE, d = 1, nRandom = nRandom, outputPath = outputPath, fileName = "UC")

UC_Corr_result <- UC_Wrapper(y = y, det.type = "trend", corr.ind = TRUE, d = 1, nRandom = nRandom, outputPath = outputPath, fileName = "UC_Corr")

UC_Frac_result <- UC_Wrapper(y = y, det.type = "frac", corr.ind = TRUE, d = FALSE, nRandom = nRandom, outputPath = outputPath, fileName = "UC_Frac")
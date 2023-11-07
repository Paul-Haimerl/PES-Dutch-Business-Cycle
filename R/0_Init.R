#' @description Function that loads required packages or installs them if necessary
#' @param string of a package name

LoadPackages <- function(pckg) {
  if (!require(pckg, character.only = TRUE)) {
    install.packages(pckg, dep = TRUE, quiet = T)
  }
  require(pckg, character.only = TRUE, quietly = T)
}


#' @description Function that sets up the project
#' @param outputPath to the folder where the regression results are to be stored in
#' @return tibble with the data set

InitProject <- function(outputPath) {
  # Load packages
  pckgs <- c(
    "openxlsx", "readxl", "tidyverse", "fUCpack", "bHP", "pbapply", "lubridate"
  )
  suppressWarnings(sapply(pckgs, LoadPackages))

  # Create the output folder
  if (!dir.exists(outputPath)) {
    dir.create(outputPath)
  }

  # Load other scripts
  scripts <- paste0("R/", list.files(path = paste0(getwd(), "/R")))
  sapply(scripts, source)

  # Load data
  suppressMessages({
    yTib <<- read_csv("GDP_NL.csv") %>%
      mutate(logGDP = log(GDP),
             Date = mdy(Date))
  })
  y <<- yTib$logGDP
}

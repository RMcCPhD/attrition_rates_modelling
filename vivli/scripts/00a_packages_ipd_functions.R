
# Libraries 
library(tidyverse)
library(readxl)
library(styler)
library(haven)
library(broom)
library(survival)
library(muhaz)
library(flexsurv)
library(gt)
library(RDP)
library(gridExtra)

# Data extraction and version control ------------------------------------------

options(scipen = 999)

ListZip <- function(ipdname = "V:/SYR-322-INS-011/SYR-322-INS-001_IPD.zip",
                    separate = TRUE) {
  a <- unzip(zipfile = ipdname, list = TRUE)
  names(a) <- str_to_lower(names(a))
  if (separate) {
    a <- a %>%
      separate(
        col = "name", into = c("folder", "file"),
        sep = fixed("/"), remove = FALSE
      )
  } else {
    a <- a %>%
      mutate(file = name)
  }

  a %>%
    mutate(file_simple = str_remove(file, fixed(".sas7bdat")))
}

ExtractLabel <- function(x, return_object = FALSE, ncharprint = 40) {
  lbls <- map(x, ~ attr(.x, "label"))
  names(lbls) <- names(x)
  lbls <- unlist(lbls)

  # Following code adds in name if no label
  y <- names(x)
  names(y) <- y
  y[names(lbls)] <- lbls
  lbls <- y

  if (return_object) {
    lbls
  } else {
    lbls[] <- str_sub(lbls, 1, pmin(ncharprint, str_length(lbls)))
    print(as.data.frame(lbls))
  }
}

VersionControl <- function(desc) {
  # Load tidyverse
  library(tidyverse)

  # Create versions if it does not exist, otherwise load in
  if (!file.exists("versions.Rds")) {
    versions <- list()
  } else {
    versions <- readRDS("versions.Rds")
  }

  # Get data
  datetime <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  scripts <- map(list.files("Scripts", full.names = TRUE), read_lines)
  script_names <- list.files("Scripts", full.names = TRUE)
  names(scripts) <- script_names

  # Add new entry
  new_entry <- list(
    desc = desc,
    scripts = scripts
  )
  versions[[datetime]] <- new_entry

  # Save
  saveRDS(versions, "versions.rds")
}


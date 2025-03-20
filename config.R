
rm(list = ls())
gc()

options(scipen = 999)
set.seed(123)

# Unload MASS after use to avoid dplyr conflicts
detach("package:MASS", unload = TRUE)
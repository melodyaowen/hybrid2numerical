# Run this before any other code to make sure all packages are installed and loaded

# Package names
packages <- c("tidyverse", "lme4", "MASS", "Matrix", "matrixcalc",
              "ICC", "nlme", "bindata", "gee", "crt2power", "mosaic",
              "tmvtnorm", "ggplot2", "latex2exp")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))

## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(gdalraster)

set_config_option("GDAL_NUM_THREADS", "ALL_CPUS")
# unset:
set_config_option("GDAL_NUM_THREADS", "")


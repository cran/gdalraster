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

## -----------------------------------------------------------------------------
# bilinear interpolation (2x2 neighborhood of pixels)
set_config_option("GDAL_RASTERIO_RESAMPLING", "BILINEAR")

## -----------------------------------------------------------------------------
set_config_option("CPL_TMPDIR", "<dirname>") # tmpdir to use

## -----------------------------------------------------------------------------
# set to a specific size in MB
set_config_option("GDAL_CACHEMAX", "800")

# or percent of physical RAM
set_config_option("GDAL_CACHEMAX", "10%")

## -----------------------------------------------------------------------------
# default is 100
set_config_option("GDAL_MAX_DATASET_POOL_SIZE", "450")

## -----------------------------------------------------------------------------
# use COPY for inserting to PostGIS
set_config_option("PG_USE_COPY", "YES")

## -----------------------------------------------------------------------------
# SQLite: GPKG (.gpkg) and Spatialite (.sqlite)
# enable extra buffering/caching by the GDAL/OGR I/O layer
set_config_option("SQLITE_USE_OGR_VFS", "YES")

## -----------------------------------------------------------------------------
# configure SQLite to store the rollback journal in RAM
set_config_option("OGR_SQLITE_JOURNAL", "MEMORY")

## -----------------------------------------------------------------------------
# SFSQL/WKT1_SIMPLE/WKT1/WKT1_GDAL/WKT1_ESRI/WKT2_2015/WKT2_2018/WKT2/DEFAULT
set_config_option("OSR_WKT_FORMAT", "WKT2")

## -----------------------------------------------------------------------------
# note this also affects several other parts of GDAL
set_config_option("GDAL_NUM_THREADS", "4") # number of threads or ALL_CPUS

## -----------------------------------------------------------------------------
# specify the number of worker threads or ALL_CPUS
# note this also affects several other parts of GDAL
set_config_option("GDAL_NUM_THREADS", "ALL_CPUS")

## -----------------------------------------------------------------------------
# applies to external overviews (.ovr), and internal overviews if GDAL >= 3.6
# LZW is a good default but several other compression algorithms are available
set_config_option("COMPRESS_OVERVIEW", "LZW")

## -----------------------------------------------------------------------------
# horizontal differencing
set_config_option("PREDICTOR_OVERVIEW", "2")

## -----------------------------------------------------------------------------
# public bucket no AWS account required
set_config_option("AWS_NO_SIGN_REQUEST", "YES")

## -----------------------------------------------------------------------------
set_config_option("AWS_ACCESS_KEY_ID", "<value>") # key ID
set_config_option("AWS_SECRET_ACCESS_KEY", "<value>") # secret access key
# used for validation if using temporary credentials:
set_config_option("AWS_SESSION_TOKEN", "<value>") # session token
# if requester pays:
set_config_option("AWS_REQUEST_PAYER", "<value>") # requester

## -----------------------------------------------------------------------------
# specify region
set_config_option("AWS_REGION", "us-west-2")

## -----------------------------------------------------------------------------
# SOZip optimization defaults to AUTO
set_config_option("CPL_SOZIP_ENABLED", "YES")

## -----------------------------------------------------------------------------
# SOZip minimum file size
set_config_option("CPL_SOZIP_MIN_FILE_SIZE", "100K")


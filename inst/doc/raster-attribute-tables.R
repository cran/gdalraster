## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(gdalraster)

# LANDFIRE Existing Vegetation Type (EVT)
evt_file <- system.file("extdata/storml_evt.tif", package="gdalraster")

# make a copy to modify
f <- paste0(tempdir(), "/", "storml_evt_tmp.tif")
file.copy(evt_file,  f)
ds <- new(GDALRaster, f, read_only=FALSE)
ds$getDefaultRAT(band=1)

# get the full attribute table for LANDFIRE EVT from its CSV file
evt_csv <- system.file("extdata/LF20_EVT_220.csv", package="gdalraster")
evt_df <- read.csv(evt_csv)
nrow(evt_df)
head(evt_df)

# keep just the R, G, B fields (0-255) and drop RED, GREEN, BLUE
evt_df <- evt_df[,1:7]

# build a RAT for the EVT raster, attaching additional columns from evt_df
tbl <- buildRAT(ds,
                table_type = "thematic",
                na_value = -9999,
                join_df = evt_df)

nrow(tbl)
head(tbl)

# attributes on the returned data frame and its columns define RAT metadata
attr(tbl, "GDALRATTableType")
attributes(tbl$VALUE)     # GFU_MinMax for column of discrete pixel values
attributes(tbl$COUNT)     # pixel counts
attributes(tbl$EVT_NAME)  # the class names
attributes(tbl$EVT_LF)    # ancillary attribute
attributes(tbl$EVT_PHYS)  # ancillary attribute
attributes(tbl$R)         # red 0-255
attributes(tbl$G)         # green 0-255
attributes(tbl$B)         # blue 0-255

# set as default RAT on the EVT raster
ds$setDefaultRAT(band=1, tbl)
ds$flushCache()

# it can now be read from the raster dataset
rm(tbl)
tbl <- ds$getDefaultRAT(band=1)
nrow(tbl)

## ----fig.width=6, fig.height=4, dev="png"-------------------------------------
bb <- ds$bbox()
plot_raster(data = ds,
            col_tbl = tbl[,c(1,6:8)],
            maxColorValue = 255,
            interpolate = FALSE,
            main = "Storm Lake LANDFIRE EVT")

## -----------------------------------------------------------------------------
displayRAT(tbl, title = "Raster Attribute Table for Storm Lake EVT")

## -----------------------------------------------------------------------------
ds$close()

## ----echo=FALSE, fig.cap="LANDFIRE EVT in the Raster Attribute Table QGIS Plugin", out.width = '90%'----
knitr::include_graphics("qgis_rat_classify.png")


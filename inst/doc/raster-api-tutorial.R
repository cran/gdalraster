## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(gdalraster)

tcc_file <- system.file("extdata/storml_tcc.tif", package="gdalraster")
ds <- new(GDALRaster, tcc_file, read_only = TRUE)

## -----------------------------------------------------------------------------
ds
str(ds)

## -----------------------------------------------------------------------------
gt <- ds$getGeoTransform()
gt[1]  # x-coordinate of upper-left corner of the upper-left pixel
gt[2]  # pixel width (w-e resolution)
gt[3]  # 0 for north-up
gt[4]  # y-coordinate of upper-left corner of the upper-left pixel
gt[5]  # 0 for north-up
gt[6]  # pixel height (n-s resolution, negative value)

## -----------------------------------------------------------------------------
ds$bbox()  # xmin, ymin, xmax, ymax
ds$res()   # pixel width, pixel height as positive values

## -----------------------------------------------------------------------------
# GDAL format driver
ds$getDriverShortName()
ds$getDriverLongName()

# raster size in pixels, number of bands
ds$getRasterXSize()
ds$getRasterYSize()
ds$getRasterCount()
ds$dim()

# coordinate reference system as WKT string
ds$getProjectionRef()

# origin and pixel size from the geotransform
print(paste("Origin:", gt[1], gt[4]))
print(paste("Pixel size:", gt[2], gt[6]))

## -----------------------------------------------------------------------------
# block	size
ds$getBlockSize(band = 1)

# data type
ds$getDataTypeName(band = 1)

# nodata value
ds$getNoDataValue(band = 1)

# min, max, mean, sd of pixel values in the band
ds$getStatistics(band = 1, approx_ok = FALSE, force = TRUE)

# does this band have overviews? (aka "pyramids")
ds$getOverviewCount(band = 1)

# does this band have a color table?
col_tbl <- ds$getColorTable(band = 1)
if (!is.null(col_tbl))
  head(col_tbl)

## -----------------------------------------------------------------------------
# read the first row of pixel values
ncols <- ds$getRasterXSize()
rowdata <- ds$read(band = 1,
                   xoff = 0,
                   yoff = 0,
                   xsize = ncols,
                   ysize = 1,
                   out_xsize = ncols,
                   out_ysize = 1)

length(rowdata)
typeof(rowdata)
head(rowdata)

## ----fig.width=6, fig.height=4, dev="png"-------------------------------------
plot_raster(ds, legend = TRUE, main = "Storm Lake Tree Canopy Cover (%)")

## -----------------------------------------------------------------------------
# close the dataset for proper cleanup
ds$close()

## -----------------------------------------------------------------------------
lcp_file <- system.file("extdata/storm_lake.lcp", package="gdalraster")
tif_file <- file.path(tempdir(), "storml_lndscp.tif")
ds <- createCopy(format = "GTiff",
                 dst_filename = tif_file,
                 src_filename = lcp_file,
                 options = "COMPRESS=LZW",
                 return_obj = TRUE)

# band = 0 for dataset-level metadata:
ds$getMetadata(band = 0, domain = "IMAGE_STRUCTURE")

# set nodata value for all bands
for (band in 1:ds$getRasterCount())
  ds$setNoDataValue(band, -9999)

# band 2 of an LCP file is slope degrees
ds$getStatistics(band = 2, approx_ok = FALSE, force = TRUE)
ds$close()

vsi_stat(lcp_file, "size")
vsi_stat(tif_file, "size")

## -----------------------------------------------------------------------------
getCreationOptions("GTiff", "COMPRESS")

getCreationOptions("GTiff", "SPARSE_OK")

## -----------------------------------------------------------------------------
new_file <- file.path(tempdir(), "newdata.tif")
ds <- create(format = "GTiff",
             dst_filename = new_file,
             xsize = 143,
             ysize = 107,
             nbands = 1, 
             dataType = "Int16",
             return_obj = TRUE)

## -----------------------------------------------------------------------------
# EPSG:26912 - NAD83 / UTM zone 12N
ds$setProjection(epsg_to_wkt(26912))

gt <- c(323476.1, 30, 0, 5105082.0, 0, -30)
ds$setGeoTransform(gt)

ds$setNoDataValue(band = 1, -9999)
ds$fillRaster(band = 1, -9999, 0)

# ...

# close the dataset when done
ds$close()


## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(gdalraster)

elev_file <- system.file("extdata/storml_elev.tif", package="gdalraster")
ds <- new(GDALRaster, filename = elev_file, read_only = TRUE)

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

# coordinate reference system
ds$getProjectionRef()

# origin and pixel size from the geotransform
print(paste("Origin:", gt[1], gt[4]))
print(paste("Pixel size:", gt[2], gt[6]))

## -----------------------------------------------------------------------------
# block	size
ds$getBlockSize(band=1)

# data type
ds$getDataTypeName(band=1)

# nodata value
ds$getNoDataValue(band=1)

# min, max, mean, sd of pixel values in the band
ds$getStatistics(band=1, approx_ok = FALSE, force = TRUE)

# does this band have overviews? (aka "pyramids")
ds$getOverviewCount(band=1)

# gdalraster currently does not support access to color tables

## -----------------------------------------------------------------------------
ncols <- ds$getRasterXSize()
rowdata <- ds$read(band=1, 
                   xoff=0, yoff=0,
                   xsize=ncols, ysize=1,
                   out_xsize=ncols, out_ysize=1)

length(rowdata)
typeof(rowdata)
head(rowdata)

## -----------------------------------------------------------------------------
# close the dataset for proper cleanup
ds$close()

## -----------------------------------------------------------------------------
lcp_file <- system.file("extdata/storm_lake.lcp", package="gdalraster")
tif_file <- paste0(tempdir(), "/", "storml_lndscp.tif")
options <- c("COMPRESS=LZW")
createCopy(format="GTiff", dst_filename=tif_file, src_filename=lcp_file, 
           options=options)

file.size(lcp_file)
file.size(tif_file)

ds <- new(GDALRaster, tif_file, read_only=FALSE)

# band=0 for dataset-level metadata:
ds$getMetadata(band=0, domain="IMAGE_STRUCTURE")

# set nodata value for all bands
for (band in 1:ds$getRasterCount())
    ds$setNoDataValue(band, -9999)

# band 2 of an LCP file is slope degrees
ds$getStatistics(band=2, approx_ok=FALSE, force=TRUE)
ds$close()

## -----------------------------------------------------------------------------
new_file <- paste0(tempdir(), "/", "newdata.tif")
create(format="GTiff", dst_filename=new_file, xsize=143, ysize=107, nbands=1, 
       dataType="Int16")

## -----------------------------------------------------------------------------
ds <- new(GDALRaster, new_file, read_only=FALSE)

# EPSG:26912 - NAD83 / UTM zone 12N
ds$setProjection(epsg_to_wkt(26912))

gt <- c(323476.1, 30, 0, 5105082.0, 0, -30)
ds$setGeoTransform(gt)

ds$setNoDataValue(band=1, -9999)
ds$fillRaster(band=1, -9999, 0)

# ...

# close the dataset when done
ds$close()


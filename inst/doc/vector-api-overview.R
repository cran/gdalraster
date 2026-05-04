## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

library(gdalraster)
gdal_3_6_0_min <- (gdal_version_num() >= gdal_compute_version(3, 6, 0))
gdal_3_7_0_min <- (gdal_version_num() >= gdal_compute_version(3, 7, 0))

## -----------------------------------------------------------------------------
library(gdalraster)

# get path to the Yellowstone National Park sample dataset
f <- system.file("extdata/ynp_features.zip", package = "gdalraster")

# add the VSI prefix 
(zf <- file.path("/vsizip", f))

# list files in the zip archive
vsi_read_dir(zf)

# VSI path to the GPKG file
(zf_gpkg <- file.path(zf, "ynp_features.gpkg"))

## -----------------------------------------------------------------------------
if (gdal_version_num() >= gdal_compute_version(3, 7, 0)) {
  cat("SOZip metadata for ynp_features.gpkg:\n")
  vsi_get_file_metadata(zf_gpkg, domain = "ZIP") |> print()
} else {
  cat("SOZip support requires GDAL >= 3.7\n")
}

## -----------------------------------------------------------------------------
inspectDataset(zf_gpkg)

## -----------------------------------------------------------------------------
# test for existence of a vector data source with at least read access
ogr_ds_exists(zf_gpkg)
# note that update of an existing zipped .gpkg file is not supported
ogr_ds_exists(zf_gpkg, with_update = TRUE)

# list the vector layers
ogr_ds_layer_names(zf_gpkg)

## ----eval=gdal_3_7_0_min------------------------------------------------------
# list the layers in a data source
ogrinfo(zf_gpkg)

# detailed information about a specific layer
ogrinfo(zf_gpkg, "ynp_bnd", cl_arg = c("-geom=SUMMARY", "-wkt_format", "WKT2"))

## -----------------------------------------------------------------------------
# copy ynp_features.gpkg from the zip file to an in-memory file
mem_gpkg <- "/vsimem/tmp/ynp_features.gpkg"
ogr2ogr(zf_gpkg, mem_gpkg)

vsi_read_dir("/vsimem/tmp")

# confirm it's writable
ogr_ds_exists(mem_gpkg, with_update = TRUE)

rd_layer <- "roads"

# rename a layer (requires GDAL >= 3.5)
if (gdal_version_num() < gdal_compute_version(3, 5, 0)) {
  message("ogr_layer_rename() requires GDAL >= 3.5")
} else if (ogr_layer_test_cap(mem_gpkg, rd_layer)$Rename) {
  ogr_layer_rename(mem_gpkg, rd_layer, "roads2")
  rd_layer <- "roads2"
} else {
  message("layer does not have 'Rename' capability")
}

ogr_ds_layer_names(mem_gpkg)

# delete a layer
if (ogr_ds_test_cap(mem_gpkg)$DeleteLayer) {
  ogr_layer_delete(mem_gpkg, rd_layer)
} else {
  message("dataset does not have 'DeleteLayer' capability")
}

ogr_ds_layer_names(mem_gpkg)

# manage fields on a layer
ogr_layer_field_names(mem_gpkg, "points_of_interest")

# delete a field
if (ogr_layer_test_cap(mem_gpkg, "points_of_interest")$DeleteField) {
  ogr_field_delete(mem_gpkg, "points_of_interest", "createdate")
} else {
  message("layer does not have 'DeleteField' capability")
}

# rename fields
if (ogr_layer_test_cap(mem_gpkg, "points_of_interest")$AlterFieldDefn) {
  ogr_field_rename(mem_gpkg, "points_of_interest", "poiname", "poi_name")
  ogr_field_rename(mem_gpkg, "points_of_interest", "poitype", "poi_type")
} else {
  message("layer does not have 'AlterFieldDefn' capability")
}

## -----------------------------------------------------------------------------
# create a new field
if (ogr_layer_test_cap(mem_gpkg, "points_of_interest")$CreateField) {
  ogr_field_create(mem_gpkg, "points_of_interest",
                   fld_name = "is_geothermal",
                   fld_type = "OFTInteger",
                   fld_subtype = "OFSTBoolean")
} else {
  message("layer does not have 'CreateField' capability")
}

# edit data with SQL
sql <- "UPDATE points_of_interest SET is_geothermal = 
          CASE
            WHEN poi_type IN ('Basin', 'Geyser') THEN TRUE
            ELSE FALSE
          END"
ogr_execute_sql(mem_gpkg, sql)

ogr_layer_field_names(mem_gpkg, "points_of_interest")

## ----fig.alt = "A plot of the Yellowstone National Park (YNP) boundary in geographic coordinate system showing point locations of geothermal features. The YNP boundary polygon has background R color 'wheat'. The locations of geothermal features are solid circles with R color 'steelblue1. The x-axis label is 'longitude' and the y-axis label is 'latitude'. The plot title is 'YNP Geothermal Features."----
# read and display the geothermal features
sql <- "SELECT poi_name, geom
          FROM points_of_interest
            WHERE is_geothermal = TRUE"
(lyr <- ogr_execute_sql(mem_gpkg, sql))

lyr$getFeatureCount()

# retrieve all features in the result set (cf. DBI::dbFetch())
feat_set <- lyr$fetch(-1)
head(feat_set)

# plot the park boundary
# the layer contains a single polygon feature which is piped directly to plot()
GDALVector$new(mem_gpkg, "ynp_bnd")$getNextFeature() |>
  plot(col = "wheat", xlab = "longitude", ylab = "latitude",
       main = "YNP Geothermal Features")

plot(feat_set, pch = 19, col = "steelblue1", add = TRUE)

## -----------------------------------------------------------------------------
lyr$close()

# delete a data source
vsi_unlink(mem_gpkg)  # delete a single file
# or, deleteDataset(mem_gpkg)

## -----------------------------------------------------------------------------
f <- system.file("extdata/ynp_features.zip", package = "gdalraster")
ynp_dsn <- file.path("/vsizip", f, "ynp_features.gpkg")

# the park boundary layer containing a single feature
(bnd <- new(GDALVector, ynp_dsn, "ynp_bnd"))

bnd$getFeatureCount()

# spatial reference definition as WKT string
bnd$getSpatialRef()

# xmin, ymin, xmax, ymax
bnd$bbox()

# structure of the layer definition (a.k.a. feature class definition)
bnd$getLayerDefn() |> str()

## -----------------------------------------------------------------------------
bnd_feat <- bnd$getNextFeature()
str(bnd_feat)

# no more features
bnd$getNextFeature()

bnd$close()

## ----fig.alt = "A plot of the Yellowstone National Park (YNP) boundary in geographic coordinate system showing public roads as LineString features. The YNP boundary polygon has background R color 'wheat', and the road features are shown as double-width lines with R color 'slategray'. The x-axis label is 'longitude' and the y-axis label is 'latitude'. The plot title is 'YNP Public Roads'."----
# SQL layer for public roads
sql <- "SELECT rdname, opentopubl, geom FROM roads WHERE opentopubl = 'Yes'"
(roads <- new(GDALVector, ynp_dsn, sql))

roads$getFeatureCount()

roads$getFieldNames()

roads_featset <- roads$fetch(-1)
nrow(roads_featset)

head(roads_featset)

plot(bnd_feat, col = "wheat", xlab = "longitude", ylab = "latitude",
     main = "YNP Public Roads")

plot(roads_featset, col = "slategray", lwd = 2, add = TRUE)

roads$close()

## -----------------------------------------------------------------------------
poi <- new(GDALVector, ynp_dsn, "points_of_interest")

poi$getFeatureCount()

# read progressively in batches
batch_size <- 500
batch <- 0
while (TRUE) {
    poi_featset <- poi$fetch(batch_size)
    if (nrow(poi_featset) == 0) break
    cat("batch", batch <- batch + 1, ":", nrow(poi_featset), "features\n")
    # process batch
    # ...
}

poi$close()

## ----eval=gdal_3_6_0_min, warning=FALSE---------------------------------------
# Expose an ArrowArrayStream (requires GDAL >= 3.6)

# re-open the roads layer with the required argument for type of access
roads$open(read_only = TRUE)
roads$resetReading()

# does the layer have a specialized implementation
roads$testCapability()$FastGetArrowStream

# optionally set ArrowStream options as character vector of "NAME=VALUE", e.g.,
roads$arrowStreamOptions <- "INCLUDE_FID=NO"
# the default batch size of 65,536 could also be configured with
# MAX_FEATURES_IN_BATCH=num

(stream <- roads$getArrowStream())

# get the array schema if needed
schema <- stream$get_schema()

# optionally read by batch, NULL is returned when no more features are available
# batch <- stream$get_next()
# if (!is.null(batch))
#   d_batch <- as.data.frame(batch)

# or, pull all the batches into a data frame
d <- as.data.frame(stream)
nrow(d)
head(d)

# release the stream when finished
stream$release()

# 'geom' is a list column of WKB raw vectors which can be passed to the
# Geometry API functions
geom_utm <- g_transform(d$geom,
                        srs_from = roads$getSpatialRef(),
                        srs_to = "EPSG:26912")

# add a column with road lengths in meters
d$rdlength <- g_length(geom_utm)
head(d)

roads$close()

## ----fig.alt = "A plot of the Yellowstone National Park (YNP) boundary in a projected coordinate system showing the perimeter of the 2016 Maple Fire, along with three points of interest located within the fire polygon. The YNP boundary polygon has background R color 'wheat'. The Maple Fire perimeter is shown as a filled polygon in R color 'orangered'. The three points of interest are shown as filled square symbols in R color 'royalblue'. The x-axis label is 'x' and the y-axis label is 'y'. The plot is untitled."----
# MTBS fire perimeters in Yellowstone National Park 1984-2022
f <- system.file("extdata/ynp_fires_1984_2022.gpkg", package = "gdalraster")

# copy to a temporary writable file
mtbs_dsn <- "/vsimem/tmp/ynp_fires_1984_2022.gpkg"
ogr2ogr(f, mtbs_dsn)

(fires <- new(GDALVector, mtbs_dsn, "mtbs_perims"))

srs_mtsp <- fires$getSpatialRef()  # Montana state plane metric definition

# reproject the boundary in ynp_features.gpkg to match the MTBS layer,
# returning a GDALVector object on the output layer by default
(bnd <- ogr_reproject(src_dsn = ynp_dsn, src_layer = "ynp_bnd",
                      out_dsn = mtbs_dsn, out_srs = srs_mtsp))

(bnd_feat <- bnd$getNextFeature())

bnd$close()

# reproject points_of_interest
poi <- ogr_reproject(ynp_dsn, "points_of_interest", mtbs_dsn, srs_mtsp)

# set an attribute filter to select the Maple Fire
fires$setAttributeFilter("incid_name = 'MAPLE'")
fires$getFeatureCount()

maple_fire <- fires$getNextFeature()

# use the fire polygon as a spatial filter for points_of_interest
# setSpatialFilter() expects WKT input, so convert the WKB geometry
g_wk2wk(maple_fire$geom) |> poi$setSpatialFilter()
poi$getFeatureCount()

poi$setSelectedFields(c("poiname", "poitype", "geom"))
(maple_pois <- poi$fetch(-1))

plot(bnd_feat, col = "wheat")
plot(maple_fire, col = "orangered", border = NA, add = TRUE)
plot(maple_pois, pch = 15, col = "royalblue", add = TRUE)

fires$close()
poi$close()

## ----fig.alt = "A plot of the Yellowstone National Park (YNP) boundary in a projected coordinate system showing the centroid of the park boundary polygon as a single point. The YNP boundary polygon has background R color 'wheat'. The centroid point is shown as an R 'circle plus' symbol, a circle with a plus sign inside resembling crosshairs. The x-axis label is 'x' and the y-axis label is 'y'. The plot is untitled."----
# create a feature object for the YNP centroid as a point of interest
(bnd_centroid_xy <- g_centroid(bnd_feat$geom))

feat <- list()
feat$poiname <- "YNP centroid"
feat$poitype <- "Information"
feat$createdate <- Sys.Date()
feat$editdate <- Sys.Date()
feat$geom <- g_create("POINT", bnd_centroid_xy)

# re-open the "points_of_interest" layer with update access
poi$open(read_only = FALSE)
poi$testCapability()$SequentialWrite

# create and write the new feature on the layer
poi$createFeature(feat)

# be sure pending writes are flushed
poi$syncToDisk()

# read back
fid <- poi$getLastWriteFID()
(ynp_centroid <- poi$getFeature(fid))

plot(bnd_feat, col = "wheat")
plot(ynp_centroid, pch = 10, cex = 1.5, add = TRUE)

## -----------------------------------------------------------------------------
# rewrite a feature in the "point_of_interest" layer updating the feature name
# verify the layer has random write capability
poi$testCapability()$RandomWrite

# FID 3809 is missing the trailhead name
(feat <- poi$getFeature(3809))

# update the name field and the date of the edit
feat$poiname <- "Ice Lake Trailhead"
feat$editdate <- Sys.Date()

# rewrite the feature
poi$setFeature(feat)
poi$syncToDisk()

(feat <- poi$getFeature(3809))

## -----------------------------------------------------------------------------
# delete the "YNP centroid" feature that was created above
# verify the layer has delete feature capability
poi$testCapability()$DeleteFeature

# the feature ID was obtained above as: fid <- poi$getLastWriteFID()
poi$getFeature(fid)

poi$deleteFeature(fid)
poi$syncToDisk()

poi$getFeature(fid)

poi$close()

## -----------------------------------------------------------------------------
# create a layer definition for random_points
# the spatial ref was obtained above as: srs_mtsp <- fires$getSpatialRef()
defn <- ogr_def_layer("POINT", srs = srs_mtsp)
defn$pt_desc <- ogr_def_field("OFTString")
defn$create_time <- ogr_def_field("OFTDateTime")

ogr_layer_create(mtbs_dsn, "random_points", layer_defn = defn)

lyr <- new(GDALVector, mtbs_dsn, "random_points", read_only = FALSE)

bb <- g_wk2wk(bnd_feat$geom) |> bbox_from_wkt()

## -----------------------------------------------------------------------------
batch_size <- as.integer(1e5)

# create a batch of features
rndX <- sample((bb[1] + 1):(bb[3] - 1), batch_size, replace = TRUE)
rndY <- sample((bb[2] + 1):(bb[4] - 1), batch_size, replace = TRUE)
pts <- cbind(rndX, rndY)
pts_geom <- g_create("POINT", pts)
d <- data.frame(pt_desc = rep("random points batch 1", batch_size),
                create_time = rep(Sys.time(), batch_size))
d$geom <- pts_geom

# write the batch (no transaction)
system.time(res <- lyr$batchCreateFeature(d))

(all(res))

lyr$syncToDisk()

## -----------------------------------------------------------------------------
rndX <- sample((bb[1] + 1):(bb[3] - 1), batch_size, replace = TRUE)
rndY <- sample((bb[2] + 1):(bb[4] - 1), batch_size, replace = TRUE)
pts <- cbind(rndX, rndY)
pts_geom <- g_create("POINT", pts)
d <- data.frame(pt_desc = rep("random points batch 2", batch_size),
                create_time = rep(Sys.time(), batch_size))
d$geom <- pts_geom

# write the batch using a transaction
system.time({
  lyr$startTransaction()
  res2 <- lyr$batchCreateFeature(d)
  if (all(res2))
    lyr$commitTransaction()
  else
    lyr$rollbackTransaction()
})

(all(res2))

# check the output data
d_out <- lyr$fetch(-1)
(nrow(d_out) == batch_size * 2)

head(d_out)

tail(d_out)

lyr$close()

## -----------------------------------------------------------------------------
# write the Maple Fire AOI bounding box as GeoJSON in EPSG 3857
json_file <- file.path(tempdir(), "maple_fire_aoi.geojson")

lyr <- ogr_ds_create("GeoJSON", json_file, layer = "maple_fire_aoi",
                     geom_type = "POLYGON", srs = "EPSG:3857",
                     fld_name = "id", fld_type = "OFTString",
                     overwrite = TRUE, return_obj = TRUE)

# The Maple Fire feature object and spatial reference were obtained above in
# the section on "Reproject vector layers".
# Here we extend the minimum bounding box by 500 m in each direction.
feat <- list()
feat$id <- "dataDownloadBox"
feat$geom <- g_transform(maple_fire$geom, srs_from = srs_mtsp,
                         srs_to = "EPSG:3857", as_wkb = FALSE) |>
               bbox_from_wkt(extend_x = 500, extend_y = 500) |>
               bbox_to_wkt()

lyr$createFeature(feat)

lyr$close()

readLines(json_file) |> writeLines()

## ----fig.alt = "A plot of the 1988 North Fork fire perimeter showing areas within the North Fork burn scar that have subsequently re-burned as of 2022. The North Fork fire is shown as an unfilled polygon with black outline. Areas within the North Fork polygon that have re-burned are shown as filled polygons in R color 'orangered'. The x-axis label is 'x' and the y-axis label is 'y'. The plot title is '1988 North Fork fire perimeter showing re-burned areas in red'."----
# layer filtered to fires after 1988
lyr1 <- new(GDALVector, mtbs_dsn, "mtbs_perims")
lyr1$setAttributeFilter("ig_year > 1988")
lyr1$getFeatureCount()

# second layer for the 1988 North Fork fire perimeter
sql <- "SELECT incid_name, ig_year, geom FROM mtbs_perims
        WHERE incid_name = 'NORTH FORK'"
lyr2 <- new(GDALVector, mtbs_dsn, sql)
lyr2$getFeatureCount()

north_fork_feat <- lyr2$getNextFeature()

# set mode options for the intersection
opt <- c("INPUT_PREFIX=layer1_",
         "METHOD_PREFIX=layer2_",
         "PROMOTE_TO_MULTI=YES")

# intersect to obtain areas re-burned since 2000
lyr_out <- ogr_proc(mode = "Intersection",
                    input_lyr = lyr1,
                    method_lyr = lyr2,
                    out_dsn = mtbs_dsn,
                    out_lyr_name = "north_fork_reburned",
                    out_geom_type = "MULTIPOLYGON",
                    mode_opt = opt)

# the output layer has attributes of both the input and method layers
(reburn_feat_set <- lyr_out$fetch(-1))

plot(north_fork_feat)
plot(reburn_feat_set, col = "orangered", border = NA, add = TRUE,
     main = "1988 North Fork fire perimeter showing re-burned areas in red")

# clean up
lyr1$close()
lyr2$close()
lyr_out$close()

## ----include = FALSE----------------------------------------------------------
# clean up
unlink(json_file)
vsi_unlink(mtbs_dsn)
vsi_rmdir("/vsimem/tmp/")


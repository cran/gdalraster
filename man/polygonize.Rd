% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gdalraster_proc.R
\name{polygonize}
\alias{polygonize}
\title{Create a polygon feature layer from raster data}
\usage{
polygonize(
  raster_file,
  out_dsn,
  out_layer,
  fld_name = "DN",
  out_fmt = NULL,
  connectedness = 4,
  src_band = 1,
  mask_file = NULL,
  nomask = FALSE,
  overwrite = FALSE,
  dsco = NULL,
  lco = NULL,
  quiet = FALSE
)
}
\arguments{
\item{raster_file}{Filename of the source raster.}

\item{out_dsn}{The destination vector filename to which the polygons will be
written (or database connection string).}

\item{out_layer}{Name of the layer for writing the polygon features. For
single-layer file formats such as \code{"ESRI Shapefile"}, the layer name is the
same as the filename without the path or extension (e.g., \code{out_dsn = "path_to_file/polygon_output.shp"}, the layer name is \code{"polygon_output"}).}

\item{fld_name}{Name of an integer attribute field in \code{out_layer} to which
the pixel values will be written. Will be created if necessary when using an
existing layer.}

\item{out_fmt}{GDAL short name of the output vector format. If unspecified,
the function will attempt to guess the format from the filename/connection
string.}

\item{connectedness}{Integer scalar. Must be either \code{4} or \code{8}. For the
default 4-connectedness, pixels with the same value are considered connected
only if they touch along one of the four sides, while 8-connectedness
also includes pixels that touch at one of the corners.}

\item{src_band}{The band on \code{raster_file} to build the polygons from
(default is \code{1}).}

\item{mask_file}{Use the first band of the specified raster as a
validity mask (zero is invalid, non-zero is valid). If not specified, the
default validity mask for the input band (such as nodata, or alpha masks)
will be used (unless \code{nomask} is set to \code{TRUE}).}

\item{nomask}{Logical scalar. If \code{TRUE}, do not use the default validity
mask for the input band (such as nodata, or alpha masks).
Default is \code{FALSE}.}

\item{overwrite}{Logical scalar. If \code{TRUE}, overwrite \code{out_layer} if it
already exists. Default is \code{FALSE}.}

\item{dsco}{Optional character vector of format-specific creation options
for \code{out_dsn} (\code{"NAME=VALUE"} pairs).}

\item{lco}{Optional character vector of format-specific creation options
for \code{out_layer} (\code{"NAME=VALUE"} pairs).}

\item{quiet}{Logical scalar. If \code{TRUE}, a progress bar will not be
displayed. Defaults to \code{FALSE}.}
}
\description{
\code{polygonize()} creates vector polygons for all connected regions of pixels
in a source raster sharing a common pixel value. Each polygon is created
with an attribute indicating the pixel value of that polygon. A raster mask
may also be provided to determine which pixels are eligible for processing.
The function will create the output vector layer if it does not already
exist, otherwise it will try to append to an existing one.
This function is a wrapper of \code{GDALPolygonize} in the GDAL Algorithms API.
It provides essentially the same functionality as the \code{gdal_polygonize.py}
command-line program (\url{https://gdal.org/en/stable/programs/gdal_polygonize.html}).
}
\details{
Polygon features will be created on the output layer, with polygon
geometries representing the polygons. The polygon geometries will be in the
georeferenced coordinate system of the raster (based on the geotransform of
the source dataset). It is acceptable for the output layer to already have
features. If the output layer does not already exist, it will be created
with coordinate system matching the source raster.

The algorithm attempts to minimize memory use so that very large rasters can
be processed. However, if the raster has many polygons or very large/complex
polygons, the memory use for holding polygon enumerations and active polygon
geometries may grow to be quite large.

The algorithm will generally produce very dense polygon geometries, with
edges that follow exactly on pixel boundaries for all non-interior pixels.
For non-thematic raster data (such as satellite images) the result will
essentially be one small polygon per pixel, and memory and output layer
sizes will be substantial. The algorithm is primarily intended for
relatively simple thematic rasters, masks, and classification results.
}
\note{
The source pixel band values are read into a signed 64-bit integer buffer
(\code{Int64}) by \code{GDALPolygonize}, so floating point or complex bands will be
implicitly truncated before processing.

When 8-connectedness is used, many of the resulting polygons will likely be
invalid due to ring self-intersection (in the strict OGC definition of
polygon validity). They may be suitable as-is for certain purposes such as
calculating geometry attributes (area, perimeter). Package \strong{sf} has
\code{st_make_valid()}, PostGIS has \code{ST_MakeValid()}, and QGIS has vector
processing utility "Fix geometries" (single polygons can become MultiPolygon
in the case of self-intersections).

If writing to a SQLite database format as either \code{GPKG} (GeoPackage
vector) or \code{SQLite} (Spatialite vector), setting the
\code{SQLITE_USE_OGR_VFS} and \code{OGR_SQLITE_JOURNAL} configuration options may
increase performance substantially. If writing to \code{PostgreSQL}
(PostGIS vector), setting \code{PG_USE_COPY=YES} is faster:

\if{html}{\out{<div class="sourceCode">}}\preformatted{# SQLite: GPKG (.gpkg) and Spatialite (.sqlite)
# enable extra buffering/caching by the GDAL/OGR I/O layer
set_config_option("SQLITE_USE_OGR_VFS", "YES")
# set the journal mode for the SQLite database to MEMORY
set_config_option("OGR_SQLITE_JOURNAL", "MEMORY")

# PostgreSQL / PostGIS
# use COPY for inserting data rather than INSERT
set_config_option("PG_USE_COPY", "YES")
}\if{html}{\out{</div>}}
}
\examples{
evt_file <- system.file("extdata/storml_evt.tif", package="gdalraster")
dsn <- file.path(tempdir(), "storm_lake.gpkg")
layer <- "lf_evt"
fld <- "evt_value"
set_config_option("SQLITE_USE_OGR_VFS", "YES")
set_config_option("OGR_SQLITE_JOURNAL", "MEMORY")
polygonize(evt_file, dsn, layer, fld)
set_config_option("SQLITE_USE_OGR_VFS", "")
set_config_option("OGR_SQLITE_JOURNAL", "")
\dontshow{deleteDataset(dsn)}
}
\seealso{
\code{\link[=rasterize]{rasterize()}}

\code{vignette("gdal-config-quick-ref")}
}

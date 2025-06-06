% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gdal_create.R
\name{createCopy}
\alias{createCopy}
\title{Create a copy of a raster}
\usage{
createCopy(
  format,
  dst_filename,
  src_filename,
  strict = FALSE,
  options = NULL,
  quiet = FALSE,
  return_obj = FALSE
)
}
\arguments{
\item{format}{Character string giving the format short name for the
output raster (e.g., \code{"GTiff"}).}

\item{dst_filename}{Character string giving the filename to create.}

\item{src_filename}{Either a character string giving the filename of the
source raster, or object of class \code{GDALRaster} for the source.}

\item{strict}{Logical. \code{TRUE} if the copy must be strictly equivalent,
or more normally \code{FALSE} (the default) indicating that the copy may adapt
as needed for the output format.}

\item{options}{Optional list of format-specific creation options in a
vector of \code{"NAME=VALUE"} pairs
(e.g., \code{options = c("COMPRESS=LZW")} to set \code{LZW}
compression during creation of a GTiff file).
The APPEND_SUBDATASET=YES option can be
specified to avoid prior destruction of existing dataset.}

\item{quiet}{Logical scalar. If \code{TRUE}, a progress bar will be not be
displayed. Defaults to \code{FALSE}.}

\item{return_obj}{Logical scalar. If \code{TRUE}, an object of class
\code{\link{GDALRaster}} opened on the newly created dataset will be
returned. Defaults to \code{FALSE}.}
}
\value{
By default, returns a logical value indicating success (invisible
\code{TRUE}, output written to \code{dst_filename}). An error is raised if the
operation fails. An object of class \code{\link{GDALRaster}} open on the
output dataset will be returned if \code{return_obj = TRUE}.
}
\description{
\code{createCopy()} copies a raster dataset, optionally changing the format.
The extent, cell size, number of bands, data type, projection, and
geotransform are all copied from the source raster.
}
\note{
\code{dst_filename} may be an empty string (\code{""}) with \code{format = "MEM"} and
\code{return_obj = TRUE} to create an In-memory Raster
(\url{https://gdal.org/en/stable/drivers/raster/mem.html}).
}
\examples{
lcp_file <- system.file("extdata/storm_lake.lcp", package="gdalraster")
tif_file <- file.path(tempdir(), "storml_lndscp.tif")
ds <- createCopy(format = "GTiff",
                 dst_filename = tif_file,
                 src_filename = lcp_file,
                 options = "COMPRESS=LZW",
                 return_obj = TRUE)

ds$getMetadata(band = 0, domain = "IMAGE_STRUCTURE")

for (band in 1:ds$getRasterCount())
    ds$setNoDataValue(band, -9999)
ds$getStatistics(band = 1, approx_ok = FALSE, force = TRUE)

ds$close()
\dontshow{deleteDataset(tif_file)}
}
\seealso{
\code{\link[=GDALRaster]{GDALRaster-class}}, \code{\link[=create]{create()}}, \code{\link[=getCreationOptions]{getCreationOptions()}},
\code{\link[=rasterFromRaster]{rasterFromRaster()}}, \code{\link[=translate]{translate()}}
}

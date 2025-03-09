#' Generar capas geográficas basadas en una grilla y datos raster
#'
#' Esta función toma una grilla de tipo `spread.grid`, un modelo digital de elevación (DEM),
#' y una capa de uso del suelo para calcular matrices representativas de las celdas de la grilla.
#'
#' @param sgrid Un objeto de clase `spread.grid`, que define la grilla espacial.
#' @param dem Un objeto SpatRaster que representa un modelo digital de elevación (DEM).
#' @param land Un objeto SpatRaster que representa el uso del suelo o alguna otra variable categórica.
#' @param cellra Un valor numérico opcional que define el radio del buffer en metros; si no se especifica,
#' se usará la mitad del tamaño de celda de la grilla (\code{sgrid@dgrid / 2}).
#'
#' @return Una lista con dos elementos:
#' \itemize{
#'   \item \code{dem}: Una matriz (nx x ny) con los valores promedio del DEM dentro de cada celda.
#'   \item \code{land}: Una matriz (nx x ny) con los valores más comunes del uso del suelo dentro de cada celda.
#' }
#'
#' @examples
#' # Crear un objeto spread.grid
#' library(terra)
#' dominio <- spread_grid(xorig = 298811, yorig = 4106479, nx = 100,
#'                        ny = 100, dgrid = 25, epsg = 26911)
#'
#' # Datos geográficos
#' path_dem <- system.file("extdata", "dem.tif", package="rspread")
#' terreno <- terra::rast(path_dem)/3.2808 # Conversión de pies a metros
#' path_land <- system.file("extdata", "land.tif", package="rspread")
#' cobertura <- terra::rast(path_land)
#' geodatos <- spreadgeo(dominio, terreno, cobertura, fast=TRUE, type="NLCD")
#'
#' @export
spreadgeo <- function(sgrid, dem, land, cellra=NULL, fast=FALSE, type=NULL){

  if (type == "ESA") {
    link_df <- data.frame("ESA" = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 100),
                          "value" = c(4,5,3,3,6,1,7,7,3,4,3))
    land <- terra::classify(land, link_df)
  } else if (type == "NLCD") {
    link_df <- data.frame("NLCD" = c(11,12,21,22,23,24,31,32,41,42,43,51,52,71,
                                     72,73,74,81,82,90,91,92,93,94,95,96,97,98,99),
                          "value" = c(7,7,6,6,6,6,1,1,4,2,2,5,5,3,3,3,3,3,3,4,4,
                                      5,4,5,3,3,3,3,3))
    land <- terra::classify(land, link_df)
  } else if (is.null(type)) {
    u_val <- unique(terra::values(land, mat=F))
    if (!all(u_val %in% 1:7)) {
      stop("La capa de uso del suelo no tiene valores categóricos válidos\n")
    }
  }

  if(fast){
    r <- to_rast(sgrid)
    values(r) <- 0
    e <- terra::resample(dem, r, method = "mode")
    em <- as.matrix(e, wide=T)
    l <- terra::resample(land, r, method = "mode")
    lc <- as.matrix(l, wide=T)
  } else {
    p <- terra::as.points(to_rast(sgrid), values=F)
    if(is.null(cellra)){
      cellra <- sgrid@dgrid/2
    }
    b <- terra::buffer(p, cellra, capstyle="square")
    values(b) <- 1:sgrid@total_points
    e <- terra::zonal(dem, b, mean, na.rm=T)
    em <- matrix(e[, 1], nrow=sgrid@ny, ncol=sgrid@nx, byrow=T)

    getmode <- function(v) {
      uniqv <- unique(v)
      uniqv[which.max(tabulate(match(v, uniqv)))]
    }

    l <- terra::extract(land, b, getmode, ID=F)
    lc <- matrix(l[, 1], nrow=sgrid@ny, ncol=sgrid@nx, byrow=T)
  }

  return(list("dem"=em, "land"=lc))
}

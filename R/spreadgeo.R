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
#' geodatos <- spreadgeo(dominio, terreno, cobertura)
#'
#' @export
spreadgeo <- function(sgrid, dem, land, cellra=NULL){

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
  lm <- matrix(l[, 1], nrow=sgrid@ny, ncol=sgrid@nx, byrow=T)

  return(list("dem"=em, "land"=lm))
}

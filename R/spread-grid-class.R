#' spread.grid: Clase para representar una grilla espacial
#'
#' @slot xorig Coordenada X de origen.
#' @slot yorig Coordenada Y de origen.
#' @slot nx Número de celdas en el eje X.
#' @slot ny Número de celdas en el eje Y.
#' @slot dgrid Tamaño de las celdas.
#' @slot epsg Código EPSG del sistema de referencia.
#' @slot total_points Número total de puntos.
#'
#' @exportClass spread.grid

# Definir la clase S4
setClass("spread.grid",
         slots = list(xorig = "numeric", yorig = "numeric",
                      nx = "numeric", ny = "numeric",
                      dgrid = "numeric", epsg = "numeric",
                      total_points = "numeric"))

# Definir el método initialize para la clase spread.grid
setMethod("initialize", "spread.grid", function(.Object, xorig, yorig, nx, ny, dgrid, epsg) {
  .Object@xorig <- xorig
  .Object@yorig <- yorig
  .Object@nx <- nx
  .Object@ny <- ny
  .Object@dgrid <- dgrid
  .Object@total_points <- nx * ny
  .Object@epsg <- epsg
  .Object
})

# Definir el método show para la clase spread.grid
setMethod("show", "spread.grid",
          function(object) {
            cat("X origin:", object@xorig, "\n")
            cat("Y Origin:", object@yorig, "\n")
            cat("Number of X cells:", object@nx, "\n")
            cat("Number of Y cells:", object@ny, "\n")
            cat("Grid spacing:", object@dgrid, "\n")
            cat("Total points:", object@total_points, "\n")
            cat("EPSG code:", object@epsg, "\n")
          })

# Definir el método to_rast para la clase spread.grid
setGeneric("to_rast", function(object) standardGeneric("to_rast"))

setMethod("to_rast", "spread.grid",
          function(object) {
            # Crear un objeto rast con las especificaciones de spread.grid
            rast_obj <- rast(nrows = object@ny, ncols = object@nx,
                             xmin = object@xorig, xmax = object@xorig + object@nx * object@dgrid,
                             ymin = object@yorig, ymax = object@yorig + object@ny * object@dgrid,
                             crs = paste0("epsg:", object@epsg))
            return(rast_obj)
          })

# Definir el método to_coords para generar las coordenas en x e y
setGeneric("to_coords", function(object) standardGeneric("to_coords"))

setMethod("to_coords", "spread.grid",
          function(object) {
            r <- terra::rast(to_rast(object), val=0)
            df <- terra::as.data.frame(r, xy=T)
            mx <- matrix(df$x, nrow = object@ny, ncol = object@nx, byrow = T)
            my <- matrix(df$y, nrow = object@ny, ncol = object@nx, byrow = T)
            return(list("xcoords" = mx, "ycoords" = my))
          })


# Definir el método to_coords para generar las coordenas en x e y
setGeneric("to_points", function(object) standardGeneric("to_points"))

setMethod("to_points", "spread.grid",
          function(object) {
            r <- terra::rast(to_rast(object), val=0)
            pts <- terra::as.data.frame(r, xy=T)[,1:2] |>
              st_as_sf(coords = c("x", "y"), crs = paste0("epsg:", object@epsg))
            return(pts)
          })

# Crear la función constructora
spread_grid <- function(xorig, yorig, nx, ny, dgrid, epsg) {
  new("spread.grid", xorig = xorig, yorig = yorig, nx = nx, ny = ny, dgrid = dgrid, epsg = epsg)
}

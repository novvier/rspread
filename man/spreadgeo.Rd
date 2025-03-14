% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spreadgeo.R
\name{spreadgeo}
\alias{spreadgeo}
\title{Generar capas geográficas basadas en una grilla y datos raster}
\usage{
spreadgeo(sgrid, dem, land, cellra = NULL, fast = FALSE, type = NULL)
}
\arguments{
\item{sgrid}{Un objeto de clase `spread.grid`, que define la grilla espacial.}

\item{dem}{Un objeto SpatRaster que representa un modelo digital de elevación (DEM).}

\item{land}{Un objeto SpatRaster que representa el uso del suelo o alguna otra variable categórica.}

\item{cellra}{Un valor numérico opcional que define el radio del buffer en metros; si no se especifica,
se usará la mitad del tamaño de celda de la grilla (\code{sgrid@dgrid / 2}).}
}
\value{
Una lista con dos elementos:
\itemize{
  \item \code{dem}: Una matriz (nx x ny) con los valores promedio del DEM dentro de cada celda.
  \item \code{land}: Una matriz (nx x ny) con los valores más comunes del uso del suelo dentro de cada celda.
}
}
\description{
Esta función toma una grilla de tipo `spread.grid`, un modelo digital de elevación (DEM),
y una capa de uso del suelo para calcular matrices representativas de las celdas de la grilla.
}
\examples{
# Crear un objeto spread.grid
library(terra)
dominio <- spread_grid(xorig = 298811, yorig = 4106479, nx = 100,
                       ny = 100, dgrid = 25, epsg = 26911)

# Datos geográficos
path_dem <- system.file("extdata", "dem.tif", package="rspread")
terreno <- terra::rast(path_dem)/3.2808 # Conversión de pies a metros
path_land <- system.file("extdata", "land.tif", package="rspread")
cobertura <- terra::rast(path_land)
geodatos <- spreadgeo(dominio, terreno, cobertura, fast=TRUE, type="NLCD")

}

#' Propagación del sonido para múltiples fuentes geográficas
#'
#' Esta función calcula la propagación del sonido desde múltiples fuentes geográficas (puntos o polígonos)
#' en un dominio definido, considerando factores como condiciones atmosféricas, características geográficas
#' y la posible presencia de barreras.
#'
#' @param sgrid Un objeto de clase `spread.grid` que define el dominio espacial (grilla).
#' @param src_all Fuentes sonoras. Puede ser un objeto `sf` que contiene puntos o polígonos con una columna `ids`
#' que identifica cada fuente, o una lista de estos objetos.
#' @param data_src Un `data.frame` con las características de cada fuente. Debe contener al menos las siguientes columnas:
#' \itemize{
#'   \item \code{ids}: Identificador único de la fuente (debe coincidir con los identificadores en `src_all`).
#'   \item \code{sfreq}: Frecuencias de la fuente.
#'   \item \code{db}: Nivel de presión sonora en decibelios (dB).
#'   \item \code{sdist}: Distancia de medición del nivel de presión sonora desde la fuente (en metros).
#' }
#' @param atmos Condiciones atmosféricas. Un listado que incluye:
#' \itemize{
#'   \item \code{temp}: Temperatura en grados Celsius (°C).
#'   \item \code{hum}: Humedad relativa en porcentaje (%).
#'   \item \code{ws}: Velocidad del viento en m/s.
#'   \item \code{wd}: Dirección del viento en grados.
#'   \item \code{seas_cond}: Condición estacional (valores del 1 al 10, donde 1=verano despejado y calmado, 10=nublado y ventoso).
#' }
#' @param geo Un objeto que contiene datos de elevación (DEM en pies) y cobertura terrestre.
#' @param new_barrier Lógico. Si es \code{TRUE}, se utiliza un nuevo algoritmo de barreras; de lo contrario, se usa el método basado en cuencas hidrográficas.
#' @param all_basins Opcional. Un objeto que representa las cuencas hidrográficas, requerido si \code{new_barrier=FALSE}. Si es \code{NULL}, las cuencas se calculan a partir del DEM (no recomendado).
#' @param source_offset Distancia en metros entre la fuente y el punto de emisión inicial del sonido. Por defecto, \code{0.34}.
#' @param receiver_offset Distancia en metros entre el receptor y el nivel del suelo. Por defecto, \code{1}.
#'
#' @return Una lista con el dominio definido por \code{sgrid} y los resultados de la propagación del sonido para cada fuente.
#' Cada elemento de la lista contiene los cálculos realizados para una fuente específica.
#'
#' @details
#' - La función verifica que las fuentes sean objetos de clase `sf` y que contengan una columna `ids` con identificadores únicos.
#' - Si se proporciona una lista de fuentes, estas se procesan y combinan automáticamente.
#' - El cálculo de barreras puede realizarse mediante un nuevo algoritmo o basado en cuencas hidrográficas.
#'
#' @examples
#' # Crear un objeto spread.grid
#' dominio <- spread_grid(xorig = 298811, yorig = 4106479, nx = 100,
#'                        ny = 100, dgrid = 25, epsg = 26911)
#'
#' # Crear fuentes como ejemplo
#' library(sf)
#' library(dplyr)
#' library(terra)
#' fuente1 <- data.frame("ids" = c("sr1", "sr2"),
#'                       "geometry" = c("POINT(300158.3 4107751.4)",
#'                                      "POINT(300370.4 4107881.9)")) |>
#'                                        st_as_sf(wkt = c("geometry"), crs = 26911)
#'
#' fuente2 <- data.frame("ids" = c("sr3"), "geometry" = c(
#'   "POLYGON((300229.1 4107824.1, 300310.3 4107813.1, 300321.6 4107932.7,
#'             300225.8 4107949.2, 300229.1 4107824.1))")) |>
#'               st_as_sf(wkt = c("geometry"), crs = 26911)
#'
#' fuentes <- list(fuente1, fuente2)
#'
#' # Condiciones atmosféricas
#' atmosfera <- list("temp"=28, "rh"=50, "wd"=105, "ws"=12.8, "seas_cond"=5)
#'
#' # Datos geográficos
#' path_dem <- system.file("extdata", "dem.tif", package="rspread")
#' terreno <- terra::rast(path_dem)/3.2808 # Conversión de pies a metros
#' path_land <- system.file("extdata", "land.tif", package="rspread")
#' cobertura <- terra::rast(path_land)
#' geodatos <- spreadgeo(dominio, terreno, cobertura)
#'
#' # Cargar frecuencias
#' path_freq <- system.file("extdata", "all_freqs.csv", package="rspread")
#' frecuencias <- read.csv(path_freq)
#'
#' # Ejecutar la función
#' resultado <- soundprog_all(
#'   sgrid = dominio,
#'   src_all = fuentes,
#'   data_src = frecuencias,
#'   atmos = atmosfera,
#'   geo = geodatos)
#'
#' @export
soundprog_all <- function(sgrid,
                          src_all,
                          data_src,
                          atmos,
                          geo,
                          new_barrier = TRUE,
                          all_basins = NULL,
                          source_offset = 0.34,
                          receiver_offset = 1){

  verificate_sfobj <- function(x){
    # Verificar si el onjeto es un sf
    if(class(x)[1] != "sf"){
      stop("The source must be an sf object")
    }
    # Verificar identificador
    ncols <- colnames(x)
    if(!"ids" %in% ncols){
      stop("The source data must have an 'ids' column")
    }
    x_spit <- split(x, ~ids)
    # Verificaciones
    lapply(x_spit, function(x){
      id <- unique(x$ids)
      # Verificar si hay más de una fuente con un solo indicador
      if(nrow(x) > 1){
        stop("The source data must have a single source per id:", id)
      }
      # Verificar si hay la fuente es punto o polígono
      geom_type <- unique(st_geometry_type(x))
      if(geom_type != "POINT" & geom_type != "POLYGON"){
        stop("The source must be points or polygons")
      }
    })
    return(x_spit)
  }

  if(class(src_all) == "list"){
    src_all_list <- lapply(src_all, function(x){
      x_split <- verificate_sfobj(x)
      return(x_split)
    })
    # Unir las listas
    src_all_list <- do.call(c, src_all_list)
  } else {
    src_all_list <- verificate_sfobj(src_all)
  }

  # Verificar si hay indicadores duplicados
  src_id <- lapply(src_all_list, function(x) x$ids[1])
  src_id_u <- src_id[duplicated(src_id) | duplicated(src_id, fromLast = T)]

  if(length(src_id_u) > 0){
    src_id_u <- unique(src_id_u)
    stop("The source data must have a single source per id: ", src_id_u)
  } else {
    names(src_all_list) <- as.character(src_id)
  }

  # Recortar la tabla de fuentes
  data_all_list <- split(data_src, ~ids)

  # Empezando los cálculos por fuente

  results_all <- rep(list(NULL), length(src_id))
  names(results_all) <- src_id

  for(i in src_id){
    cat("Processing source: ", i, "\n")
    src_one <- src_all_list[[i]]
    data_one <- data_all_list[[i]]

    results_all[[i]] <- soundprog_one(
      sgrid = sgrid,
      src_one = src_one,
      data_one = data_one,
      atmos = atmos,
      geo = geo,
      new_barrier = new_barrier,
      all_basins = all_basins,
      source_offset = source_offset,
      receiver_offset = receiver_offset
    )
  }
  return(c("domain"=sgrid, results_all))
}

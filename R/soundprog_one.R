#' Calculate Noise Propagation Patterns
#'
#' This function computes noise propagation patterns considering atmospheric,
#' geographic, and barrier effects, with options for using a new or old
#' barrier approach.
#'
#' @param sgrid A spatial grid configuration representing the domain of analysis.
#' @param src_one A source of sound, either as a point or polygon (sf object).
#' @param data_one A data.frame containing sound frequency (\code{freq}),
#'   sound pressure levels (\code{db}), and measurement distance (\code{dist}).
#' @param atmos Atmospheric conditions including:
#'   \itemize{
#'     \item \code{temp}: Temperature in degrees Celsius (°C).
#'     \item \code{hum}: Relative humidity as a percentage (%).
#'     \item \code{ws}: Wind speed in meters per second (m/s).
#'     \item \code{wd}: Wind direction in degrees.
#'     \item \code{seas_cond}: Seasonal conditions, where:
#'       \enumerate{
#'         \item Clear, calm summer day.
#'         \item Clear, calm winter day.
#'         \item Clear, calm summer night.
#'         \item Clear, calm winter night.
#'         \item Clear, windy summer day.
#'         \item Clear, windy winter day.
#'         \item Clear, windy summer night.
#'         \item Clear, windy winter night.
#'         \item Cloudy, calm.
#'         \item Cloudy, windy.
#'       }
#'   }
#' @param geo Geographic data, including:
#'   \itemize{
#'     \item \code{dem}: Digital Elevation Model (DEM) raster.
#'     \item \code{land}: Land cover data.
#'   }
#' @param old_barrier Logical. If \code{TRUE}, the old barrier algorithm is used.
#'   Defaults to \code{TRUE}.
#' @param all_basins A \code{SpatRaster} representing hydrological basins.
#'   Required only if \code{old_barrier = TRUE}.
#' @param source_offset Numeric. The vertical offset of the source above the
#'   ground in meters. Defaults to \code{0.34}.
#' @param receiver_offset Numeric. The vertical offset of the receiver above
#'   the ground in meters. Defaults to \code{1}.
#'
#' @return A list containing:
#'   \itemize{
#'     \item \code{noise_propagation}: A list of noise propagation patterns for
#'       each frequency.
#'     \item \code{loss}: A list of losses, including:
#'       \itemize{
#'         \item \code{ssl}: Spherical spreading loss.
#'         \item \code{aal}: Atmospheric absorption loss.
#'         \item \code{veg}: Vegetation and ground cover loss (if
#'           \code{old_barrier = TRUE}).
#'         \item \code{wind}: Wind loss.
#'         \item \code{bar}: Barrier loss.
#'         \item \code{barwind}: Combined barrier and wind loss (if
#'           \code{old_barrier = FALSE}).
#'         \item \code{topo_zones}: Topographic zones (if \code{old_barrier = TRUE}).
#'       }
#'   }
#'
#' @examples
#' # Example usage of soundprog_one
#' # Crear un objeto spread.grid
#' dominio <- spread_grid(xorig = 298811, yorig = 4106479, nx = 100,
#'                        ny = 100, dgrid = 25, epsg = 26911)
#'
#' # Crear fuentes como ejemplo
#' library(sf)
#' library(dplyr)
#' library(terra)
#' fuente <- data.frame("ids" = c("sr1"),
#'                       "geometry" = c("POINT(300158.3 4107751.4)")) |>
#'   st_as_sf(wkt = c("geometry"), crs = 26911)
#'
#' # Datos geográficos
#' path_dem <- system.file("extdata", "dem.tif", package="rspread")
#' terreno <- terra::rast(path_dem)/3.2808 # Conversión de pies a metros
#' path_land <- system.file("extdata", "land.tif", package="rspread")
#' cobertura <- terra::rast(path_land)
#' geodatos <- spreadgeo(dominio, terreno, cobertura, fast=TRUE, type="NLCD")
#'
#' # Cargar frecuencias
#' path_freq <- system.file("extdata", "all_freqs.csv", package="rspread")
#' frecuencias <- read.csv(path_freq)
#' frecuencias <- split(frecuencias, ~ids)[[1]]
#'
#' results <- soundprog_one(
#'   sgrid = dominio,
#'   src_one = fuente,
#'   data_one = frecuencias,
#'   atmos = list("temp"=28, "rh"=50, "wd"=105, "ws"=12.8, "seas_cond"=5),
#'   geo = geodatos,
#'   old_barrier = FALSE
#' )
#'
#' @export
soundprog_one <- function(sgrid,
                          src_one,
                          data_one,
                          atmos,
                          geo,
                          old_barrier = FALSE,
                          all_basins = NULL,
                          source_offset = 0.34,
                          receiver_offset = 1,
                          all_results = FALSE){

  seascond <- data.frame(
    id = as.integer(1:10),
    cond = c("clear, calm summer day", "clear, calm winter day",
             "clear, calm summer night", "clear, calm winter night",
             "clear, windy summer day", "clear, windy winter day",
             "clear, windy summer night", "clear, windy winter night",
             "cloudy, calm", "cloudy, windy"))

  # Determinar dirección y distancia euclidiana
  euc <- euclidean_dist_dir(sgrid, src_one)

  # Centroide de la
  src_centroid <- st_centroid(src_one) |> st_coordinates() |> as.vector()

  # Calcular distancia en pies
  eucdist_ft <- euc$dist * 3.28084

  # Convertir sf a SaptVector
  if(class(src_one)[1] != "SpatVector"){
    src_spvec <- vect(src_one)
  } else {
    src_spvec <- src_one
  }

  # Extraer elevación de la fuente
  if(geomtype(src_spvec) == "points") {
    demr <- terra::rast(to_rast(sgrid), val=geo$dem)
    src_elev <- terra::extract(demr, src_spvec, ID=F, method="bilinear")[1,1]
  } else if (geomtype(src_spvec) == "polygons"){
    demr <- terra::rast(to_rast(sgrid), val=geo$dem)
    src_elev <- terra::zonal(demr, src_spvec, fun="mean")[1,1]
  } else {
    stop("The geometry type of the source must be points or polygons")
  }

  if(old_barrier){

    cat("\nEl uso del antiguo método de barreras para SPreAD-GIS ha quedado obsoleto.\n",
      "Se ha mantenido como predeterminado para mantener la compatibilidad con versiones\n",
      "anteriores del código. Sin embargo, se recomienda el uso de los nuevos cálculos de barreras..\n")

    if(is.null(all_basins)){
      cat("All basins don't exist. Calculating basins from DEM...\n")
      all_basins = delineate_basins(sgrid, geo$dem)
    } else {
      if(class(all_basins)[1] != "SpatRaster"){
        stop("The basin must be a SpatRaster")
      }
    }
  }

  # Verificar si la distancia medida para todas las frecuencias es la misma
  sdist <- unique(data_one$dist)
  if(length(sdist) > 1){
    stop("The measurement distance must be the same for all frequencies")
  }

  # CONVERTIR DEL SISTEMA MÉTRICO AL INGLÉS ANTIGUO POR RAZONES HISTÓRICAS
  temp_fr = atmos$temp * 1.8 + 32 # convert to F
  ws_mph = atmos$ws * 0.621371 # convert to mph
  mdist_ft = sdist * 3.28084 # convert from meters to ft
  wind_dir = convert_wind_dir(atmos$wd)
  src_elev_ft <- src_elev * 3.28084 # convert from meters to ft

  # CALCULAR LA PÉRDIDA POR DISPERSIÓN ESFÉRICA
  ssl = spherical_spreading_loss(eucdist_ft, mdist_ft)

  # Verficar frecuencias (REVISAR)
  sfreqs <- as.integer(data_one$freq)
  # sfreqs_def <- c(125L, 160L, 200L, 250L, 315L, 400L, 500L, 630L,
  #                 800L, 1000L, 1250L, 1600L, 2000L)
  # if(!all(sfreqs %in% sfreqs_def)){
  #   stop("Las frecuencias deben estar en el rango:\n",
  #        "125, 160, 200, 250, 315, 400, 500, 630, 800, 1000, 1250, 1600, 2000")
  # }
  name_freqs <- paste0("f", sfreqs)
  # CALCULAR LA PÉRDIDA POR ABSORCIÓN ATMOSFÉRICA
  aal <- rep(list(NULL), length(sfreqs))
  names(aal) <- name_freqs
  for(i in 1:length(sfreqs)){
    aal[[i]] <- atmospheric_absorption_loss(src_elev_ft, eucdist_ft, atmos$rh,
                                            temp_fr, sfreqs[i])
  }

  if(old_barrier){
    # CALCULATE FOLIAGE AND GROUND COVER LOSS
    veg <- foliage_groundcover_loss(sgrid, eucdist_ft, geo$land, src_centroid)
  }

  # CALCULAR LA PÉRDIDA A FAVOR Y A CONTRAVIENTO
  wind <- rep(list(NULL), length(sfreqs))
  names(wind) <- name_freqs
  for(i in 1:length(sfreqs)){
    wind[[i]] <- windloss(sgrid, wind_dir, ws_mph, atmos$seas_cond,
                          eucdist_ft, euc$dir, sfreqs[i])
  }

  # Convertir elevacion a pies
  dem_ft <- geo$dem * 3.28084

  if(old_barrier){
    # Delinear barrera para calcular pérdida de barrera y zonas topográficas
    barrier_ground = delineate_barrier(src_one, all_basins)

    # CALCULAR EFECTOS TOPOGRÁFICOS Y PÉRDIDA DE BARRERA
    bar <- rep(list(NULL), length(sfreqs))
    names(bar) <- name_freqs
    for(i in 1:length(sfreqs)){
      bar[[i]] <- topographic_barrier_effects(
        sgrid, barrier_ground$barrier, src_elev*3.28084, eucdist_ft, euc$dir,
        dem_ft, sfreqs[i])
    }

    # CALCULAR UBICACIONES DONDE DOMINAN LOS EFECTOS TIERRA, ATMOSFÉRICOS Y DE BARRERA
    topo_zones <- calculate_topozones(sgrid, src_one, dem_ft, barrier_ground$ground)

    # CALCULAR PATRONES RESUMIDOS DE PROPAGACIÓN DE RUIDO
    noise_propagation <- rep(list(NULL), length(sfreqs))
    names(noise_propagation) <- name_freqs
    for(i in 1:length(sfreqs)){
      noise_propagation[[i]] <- compute_noise_propagation(sgrid,
        data_one$db[i], ssl, aal[[i]], veg, wind[[i]], bar[[i]], topo_zones)
    }

  } else {

    # Cálculo de la distancia de la trayectoria de la barrera y la pérdida máxima de vegetación
    BPD_max_veg_loss <- calculate_barrier_path_distance_and_vegmax(
      sgrid, src_centroid, src_elev_ft, dem_ft, geo$land, eucdist_ft,
      source_offset, receiver_offset)

    # Calcular el efecto barrera
    bar <- rep(list(NULL), length(sfreqs))
    names(bar) <- name_freqs
    for(i in 1:length(sfreqs)){
      bar[[i]] <- barrier_effects_v2(BPD_max_veg_loss$BPD, sfreqs[i])
    }

    # Verificar si la barrera + viento supera los 25 dB, si es así, limitar a 25 dB
    barwind <- rep(list(NULL), length(sfreqs))
    names(barwind) <- name_freqs
    for(i in 1:length(sfreqs)){
      barwind[[i]] <- check_barrier_wind(bar[[i]], wind[[i]])
    }

    # CALCULAR PATRONES FINALES DE PROPAGACIÓN DE RUIDO
    noise_propagation <- rep(list(NULL), length(sfreqs))
    names(noise_propagation) <- name_freqs
    for(i in 1:length(sfreqs)){
      noise_propagation[[i]] <- compute_noise_propagation_v2(
        sgrid, data_one$db[i], ssl, aal[[i]], BPD_max_veg_loss$vegmax, barwind[[i]])
    }
  }

  if(old_barrier){
    if(all_results){
      results <- list("noise_propagation" = noise_propagation,
                      "loss" = list("ssl" = ssl, "aal" = aal, "veg" = veg, "wind" = wind,
                                    "bar" = bar, "topo_zones" = topo_zones))
    } else {
      results <- noise_propagation
    }

  } else {
    if(all_results){
      results <- list("noise_propagation" = noise_propagation,
                      "loss" = list("ssl" = ssl, "aal" = aal, "veg" = BPD_max_veg_loss,
                                    "wind" = wind, "bar" = bar, "barwind" = barwind))
    } else {
      results <- noise_propagation
    }
  }

  return(results)
}

# R/utils.R

delineate_basins <- function(sgrid, dem, out_folder = F){

  demr <- terra::rast(to_rast(sgrid), val=dem)
  if(!wbt_init()){
    stop("WhiteboxTools could not be initialized, use install_whitebox() to install the tools")
  }

  if(oit_folder){
    dem_file <- file.path(out_folder, "dem.tif")
    fil_file <- file.path(out_folder, "fil.tif")
    d8p_file <- file.path(out_folder, "d8p.tif")
    bsn_file <- file.path(out_folder, "bsn.tif")
  } else {
    cat("The output folder is not defined, the output will be stored in a temporal directory\n")
    dem_file <- tempfile(fileext = ".tif")
    fil_file <- tempfile(fileext = ".tif")
    d8p_file <- tempfile(fileext = ".tif")
    bsn_file <- tempfile(fileext = ".tif")
  }

  writeRaster(demr, dem_file)
  wbt_fill_depressions_wang_and_liu(dem_file, fil_file)
  wbt_d8_pointer(fil_file, d8p_file)
  wbt_basins(d8p_file, bsn_file)
  basins <- rast(bsn_file)
  return(basins)
}

euclidean_dist_dir <-  function(sgrid, vec){
  grid_pts <- to_points(sgrid)
  eucdist <- sf::st_distance(grid_pts, vec) |> as.vector()
  eucdist <- matrix(eucdist, nrow=sgrid@ny, ncol=sgrid@nx, byrow=T)

  # Calculate euclidian direction
  r <- terra::rasterize(vec, rast(to_rast(sgrid)))
  if(all(is.na(values(r)))){
    cat("Polygon minus to resolution\n")
    cent <- st_centroid(vec)
    r <- terra::rasterize(cent, rast(to_rast(sgrid)))
  }
  re <- terra::direction(r,degrees=TRUE)
  eucdir <- as.matrix(re, wide=T)
  eucdir <- ifelse(eucdir == 0, 360, eucdir)
  eucdir <- ifelse(is.na(eucdir), 0, eucdir)
  eucdir <- round(eucdir, 0)
  return(list("dist" = eucdist, "dir" = eucdir))
}

PathDistance <- function(points, surface = NULL, friction = NULL){

  if(is.null(surface)){
    stop("The surface raster is required")
  }

  # Convert to rasterLayer
  if(as.character(class(surface)) != "RasterLayer"){
    if(as.character(class(surface)) == "SpatRaster"){
      surface <- raster(surface)
    } else {
      stop("The surface raster is not a RasterLayer or SpatRaster")
    }
  }

  # Obtener resolución
  resol <- raster::res(surface)
  if(resol[1] - resol[2] > 0.001){
    stop("The resolution of the raster is not square")
  } else {
    resol <- mean(resol)
  }

  # Determinar si la friccion existe
  if(is.null(friction)){
    cat("Creating a friction raster with all values equal to 1\n")
    friction <- surface
    values(friction) <- 1
  } else {
    if(as.character(class(friction)) != "RasterLayer"){
      if(as.character(class(friction)) == "SpatRaster"){
        friction <- raster(friction)
      } else {
        stop("The surface raster is not a RasterLayer or SpatRaster")
      }
    }
    # Verificar si los raster tienen la misma dimensión y resolución
    compareRaster(friction, surface)
  }
  # Determinar el costo o friccion
  tc <- transition(friction, transitionFunction = function(x) 1/mean(x), directions = 8)
  transitionMatrix(tc)
  # tcf <- geoCorrection(tc*resol)
  tcf <- tc
  transitionMatrix(tcf)
  # Determinar la distancia vertical
  dv <- transition(surface, transitionFunction = function(x) 1/(abs(x[2] - x[1])), directions = 8)
  transitionMatrix(dv)
  # Distancia horizontal = 1
  # dh = 1
  u <- surface
  values(u) <- 1
  dh <- transition(u, transitionFunction = function(x) mean(x), directions = 8)
  # Aplicar factor a diagonales
  dh <- geoCorrection(dh)
  transitionMatrix(dh)
  # Extraer matriz de adyacencia
  adj <- adjacent(tc, cells=1:ncell(tc), directions = 8, pairs = TRUE)
  # Determinar distancia real
  dr <- tc
  dr[adj] <- (1/(sqrt(((1/dh[adj])**2 + (1/dv[adj])**2))))
  transitionMatrix(dr)
  tf <- tcf*dr
  transitionMatrix(tf)
  ac <- accCost(tf, points)
  return(ac)
}

FocalAllocation <- function(r){

  FindValue <- function(x){
    if(is.na(x[5])){
      if(all(is.na(x))){
        y = NA
      } else if(all(is.na(x[c(2, 4, 6, 8)]))) {
        y = min(x[c(1, 3, 7, 9)], na.rm = T)
      } else {
        y = min(x[c(2, 4, 6, 8)], na.rm = T)
      }
    } else {
      y = x[5]
    }
    return(y)
  }

  hasNA <- any(is.na(values(r)))
  # i = 1
  while(hasNA){
    # cat("Iteration", i, "\n")
    r <- terra::focal(r, matrix(1, nrow=3, ncol=3), fun = FindValue)
    hasNA <- any(is.na(values(r)))
    # i = i + 1
  }
  return(r)
}

convert_wind_dir <- function(wind_dir){
  if (wind_dir < 180 & wind_dir >= 0){
    wind_dir_out = wind_dir + 180
  } else if (wind_dir >= 180 & wind_dir <=360){
    wind_dir_out = wind_dir - 180
  } else {
    stop("wind value is not valid. Wind direction must be in range 0 - 360 degrees")
  }
  return(wind_dir_out)
}

spherical_spreading_loss <- function(eucdist_ft, measurement_distance){
  # Divide Euclidean distance by measurement distance (X/y)
  divide_result <- eucdist_ft / measurement_distance
  # Reclassify 0 cell at center to be Spherical Spreading Loss at 1 m
  reclass_value = 3.28084 / measurement_distance # convert ft to meters
  ssl_patch <- ifelse(divide_result > 0, divide_result, reclass_value)
  # Calculate spherical spreading loss
  ssl = 20 * log10(ssl_patch) # NMSIMGIS & ISO 9613-2
  return(ssl)
}

atmospheric_absorption_loss_core <- function(elev_m, rh, temp_k, sfreq){
  # Convert elevation to atmospheric pressure
  p_a = 101.325 * (1 - (2.25577 * (10 ** (-5)) * elev_m)) ** 5.25588

  # Convert relative humidity to molar concentration of water vapor
  C = (-6.8346 * ((273.16 / temp_k) ** 1.261)) + 4.6151
  psat_pr = 10 ** C
  h = (rh) * (psat_pr) * ((p_a / 101.325) ** (-1))

  # Calculate derived values for subsequent equations
  pa_pr = p_a / 101.325
  T_Tr = temp_k / 293.15
  e = 2.7182818284

  # Calculate frO (equation 3)
  frO = ((pa_pr) * ((24 + (4.04 * 10000)) * h) * (0.02 + h)) / (0.391 + h)

  # Calculate frN (equation 4)
  frN = pa_pr * (T_Tr ** (-0.5)) * (9 + (280 * h * (e ** (-4.170 * ((T_Tr ** (-0.33333)) - 1)))))

  # Calculate alpha (equation 5)
  term1 = 1.84 * (10 ** (-11)) * (pa_pr ** (-1)) * (T_Tr ** 0.5)
  term2 = (T_Tr ** (-2.5)) * (0.01275 * (e ** (-2239.1 / temp_k)) * (frO / ((frO ** 2) + (sfreq ** 2))))
  term3 = 0.1068 * (e ** (-3352 / temp_k)) * (frN / ((frN ** 2) + (sfreq ** 2)))
  alpha = 8.686 * (sfreq ** 2)*(term1 + term2 + term3)

  return(alpha)
}

atmospheric_absorption_loss <- function(elev_ft, eucdist_ft, rh, temp_fr, sfreq){

  # Conversions for calcualtions
  elev_m = elev_ft / 3.28084 # Converto to meters
  temp_c = (temp_fr - 32) * 5 / 9 # Convert to Celsius
  temp_k = temp_c + 273.15    # Convert to Kelvins
  # Calculate atmospheric absorption coefficient
  alpha = atmospheric_absorption_loss_core(elev_m, rh, temp_k, sfreq)
  alpha_ft = alpha / 3.28084
  # Calculate atmospheric absorption loss
  aal <- eucdist_ft * alpha_ft
  return(aal)
}

windloss <- function(sgrid, wd, ws, seas_cond, eucdist_ft, eucdir, sfreq){
  # windloss(sgrid, wind_dir, ws_mph, atmos$seas_cond, eucdist_ft, euc$dir, sfreqs[i])

  seascond <- data.frame(
    id = as.integer(1:10),
    phi = c(180, 180, 0, 180, 144, 144, 62, 70, 90, 90))

  if(ws < 0){
    stop("\nWind speed cannot be less than 0! Wind speed was", wind_sp_f)
  }

  if(ws == 0){
    wind <- matrix(0, nrow=sgrid@ny, ncol=sgrid@nx)
  } else {
    # Convert seasonal conditions to phi
    phi <- seascond$phi[seascond$id == seas_cond]

    # Compute sound propagation angles away from the source
    prop_angle <- ifelse(eucdir < 180, eucdir + 180, eucdir - 180)

    # Subtract prevailing wind direction from 180 and add to sound propagation angles
    wd_f = as.numeric(wd)
    Value1 <- prop_angle + (180 - wd_f)

    # Re-classify propagation values that are less than zero or greater than 360
    Value2 <- ifelse(Value1 >= 360, Value1 - 360, ifelse(Value1 < 0, Value1 + 360, Value1))

    # Re-classify wind angle values to range between 0 and 180
    wind_ang <- ifelse(Value2 > 180, 360 - Value2, Value2)

    # Identify upwind and downwind areas
    updownwind <- phi - wind_ang
    upwind <- ifelse(updownwind > 0, updownwind, 0)
    downwind <- ifelse(updownwind <= 0, 1, 0)

    # Calculate upwind loss
    upwind_loss <- ifelse(upwind <= 0, 0, ifelse(upwind >= 50, 25, 5.7642 * log(upwind) + 2.5664))

    # Calculate shadow zone correction to upwind loss
    if(ws > 0){
      d = 375 * (ws ** (-0.85))
      x_d <- eucdist_ft / d
      if(sfreq <= 125){
        freq_w = 125
      } else if(sfreq >= 2000){
        freq_w = 2000
      } else {
        freq_w = sfreq
      }
      path <- system.file("extdata", "table_13.dbf", package = "rspread")
      table_13 <- foreign::read.dbf(path)
      table_freq <- table_13[table_13$FREQ == sfreq, ]
      upwind_szf <- table_freq$SZF_100[findInterval(x_d, table_freq$FROM)]
      upwind_loss_c <- upwind_loss * upwind_szf / 100
    } else {
      upwind_loss_c <- upwind_loss
    }

    # Calculate downwind loss
    freq_dist <- eucdist_ft * sfreq
    downwind_loss <- ifelse(freq_dist <= 406237, 0, downwind * (4.2598 * log(freq_dist) - 55.014))

    # Combine upwind and downwind loss
    wind_loss <- upwind_loss_c + downwind_loss

    # Smooth transition between upwind and downwind areas and calculate wind loss
    windr <- terra::focal(rast(to_rast(sgrid), val=wind_loss), matrix(1,nrow=9,ncol=9),
                          fun = "mean", na.rm = T, pad = T)
    wind <- terra::as.matrix(windr, wide=T)
  }
  return(wind)
}

calculate_barrier_path_distance_and_vegmax <- function(
    sgrid, srcxy, srcz, dem, land, eucdist_ft, source_offset, receiver_offset){

  # Extract paremeters
  cellsize <- sgrid@dgrid

  # dll_path <- system.file("libs", "new_barrier.dll", package = "rspread")
  # dyn.load(dll_path)
  # check <- is.loaded("new_barrier")
  # if(!check) {
  #   stop("The new_barrier.dll did not load")
  # }

  if(!is.matrix(dem) && !is.matrix(lnd)) {
    stop("The dem and landcover are not matrices")
  }

  dim_dem <- dim(dem)
  dim_land <- dim(land)

  if(dim_dem[1] != dim_land[1] || dim_dem[2] != dim_land[2]) {
    stop("The dimensions of the dem and landcover do not match")
  }

  if(!is.numeric(srcxy) && !length(srcxy) == 3) {
    stop("The source coordinates are not numeric or do not have length 3")
  }

  if(!is.numeric(cellsize) && !length(cellsize) == 1) {
    stop("The size is not numeric or does not have length 1")
  }

  if(!is.numeric(source_offset) && !length(source_offset) == 1) {
    stop("The offset is not numeric or does not have length 1")
  }

  if(!is.numeric(receiver_offset) && !length(receiver_offset) == 1) {
    stop("The receiver offset is not numeric or does not have length 1")
  }

  if(!is.integer(land)){
    cat("The landcover is not an integer, converting to integer (truncation)")
    lnd_v <- as.integer(land)
    land <- matrix(lnd_v, nrow=dim_land[1], ncol=dim_land[2])
  }

  if(!is.double(dem)){
    cat("The dem is not a double, converting to double")
    dem_v <- as.double(dem)
    dem <- matrix(dem_v, nrow=dd[1], ncol=dd[2])
  }

  rows <- dim_dem[1]
  cols <- dim_dem[2]

  srcxyz <- c((srcxy[1] - sgrid@xorig)/cellsize + 0.5,
              rows - (srcxy[2] - sgrid@yorig)/cellsize + 0.5,
              srcz)

  barheight <- matrix(rep(0, rows*cols), nrow=rows, ncol=cols)
  bardist <- matrix(rep(0, rows*cols), nrow=rows, ncol=cols)
  vegmax <- matrix(rep(0, rows*cols), nrow=rows, ncol=cols)

  rsl <- .Fortran("new_barrier", srcxyz, rows, cols, cellsize, dem, land,
                  source_offset, receiver_offset, barheight, bardist, vegmax)

  # dyn.unload(dll_path)

  # Calculate Barrier Path Distance (ft)
  term1 <- sqrt(barheight**2 + bardist**2)
  term2 <- sqrt(barheight**2 + (eucdist_ft - bardist)**2)
  BPD <- term1 + term2 - eucdist_ft

  BPD2 <- ifelse(barheight == 0, barheight, BPD)

  # return(list("BPD" = BPD2, "vegmax"=vegmax,
  #             "BPD_pre" = BPD, "term1"=term1, "term2" = term2,
  #             "barrier_distance" = barrier_distance, "barrier_height" = barrier_height))
  return(list("BPD" = BPD2, "vegmax"=vegmax))
}

barrier_effects_v2 <- function(barrier_path_distance, freq_s){
  # Calculate barrier factor (N)
  L = ((0.0000000000005*freq_s**4) - (0.000000001*freq_s**3) - (0.0000004*freq_s**2) + (0.0028*freq_s) - (0.3051))
  bar_factor <- barrier_path_distance * L
  # Calculate barrier loss
  bar = (13.573) * (bar_factor**0.2299)
  return(bar)
}

check_barrier_wind <- function(bar, wind){

  #wind_raster = wind # wind is already a raster
  barwind_prelim <- bar + wind

  # Check if barrier plus wind exceeds 25 dB, if so, cap to 25 dB
  barwind <- ifelse(barwind_prelim > 25, 25, barwind_prelim)

  return(barwind)
}

compute_noise_propagation_v2 <- function(sgrid, source_level, ssl, aal, veg, barwind){
  # Subtract spherical spreading loss from source sound level
  sslloss <- source_level - ssl

  # Calculate cumulative spherical spreading and atmospheric absorption loss
  sslaal <- sslloss + aal

  # Calculate cumulative spherical spreading, atmospheric absorption, and vegetation loss
  salveg <- sslaal + veg

  # Calculate cumulative spherical spreading, atmospheric absorption, vegetation, wind, and barrier loss
  salvgwnbr <- salveg + barwind

  # Smooth noise propagation patterns
  smoothed <- terra::focal(rast(to_rast(sgrid), val = salvgwnbr), matrix(1,nrow=3,ncol=3),
                           fun = "mean", na.rm = T)

  smoothed <- as.matrix(smoothed, wide=T)

  source_level_m <- matrix(source_level, nrow=sgrid@ny, ncol=sgrid@nx)

  # Patch to prevent smoothing at cell of origin
  pr <- ifelse(source_level_m > 0, salvgwnbr, smoothed)

  return(list("sslloss" = sslloss, "sslaal" = sslaal, "salveg" = salveg,
              "salvgwnbr" = salvgwnbr, "pr" = pr))
}

delineate_barrier <- function(src_vec, all_basins){

  # Identify sound source basin
  sound_src_basin <- terra::extract(all_basins, src_vec, method = "simple", ID=F)[,1]
  basin <- app(all_basins, function(x) ifelse(x != sound_src_basin, NA, x))

  # Delineate barrier (ie, ridgeline) around sound source basin
  basin_exp <- terra::focal(basin, matrix(1, 3, 3), fun = "sum", na.rm = T)
  basin_exp <- app(basin_exp, function(x) {x*0+as.numeric(sound_src_basin)})
  basin_shr <- basin_exp*0+1
  basin_shr <- basin_shr |> terra::focal(matrix(1, 3, 3), function(x) ifelse(sum(x, na.rm=T) < 9, NA, x[5]))
  basin_shr <- app(basin_shr, function(x) {x*0+as.numeric(sound_src_basin)})

  basin_shr_rc <- app(basin_shr, function(x) ifelse(is.na(x), as.numeric(sound_src_basin) + 1, NA))
  barrier <- basin_shr_rc - basin_exp

  # Define areas where ground effects dominate
  ground <- app(basin_exp, function(x) ifelse(is.na(x), 0, ifelse(x == sound_src_basin, 3, x)))

  barrier_ground <- list("barrier" = barrier, "ground" = ground)

  return(barrier_ground)
}

topographic_barrier_effects <- function(sgrid, barrier, elev, eucdist_ft, eucdir,
                                        dem_ft, freq_s){

  barrier_pts <- terra::as.data.frame(barrier, xy = TRUE) |>
    st_as_sf(coords = c("x", "y"), crs = st_crs(sgrid@epsg))

  raster_vec <- c(rast(to_rast(sgrid), val = eucdir, name = "eucdir"),
                  dem_ft,
                  rast(to_rast(sgrid), val = eucdist_ft, name = "eucdist_ft"))
  barrier_pts_new <- terra::extract(raster_vec, barrier_pts, method = "simple", ID=F)

  elev_dist_mean <- barrier_pts_new |>
    dplyr::mutate(eucdir = as.integer(eucdir)) |>
    group_by(eucdir) |>
    dplyr::summarise(FREQUENCY = n(),
                     mean_dem_f = mean(dem_ft, na.rm=T),
                     min_eucdis = min(eucdist_ft, na.rm=T),
                     .groups = "drop") |>
    dplyr::mutate(elev = round(mean_dem_f, 0),
                  dist = round(min_eucdis, 0))

  # Assign barrier elevation value by Euclidean direction
  elev_bar_dir <- elev_dist_mean$elev[match(eucdir, elev_dist_mean$eucdir)]
  elev_bar_dir_r <- rast(to_rast(sgrid), val = elev_bar_dir, names = "elev_bar_dir")
  elev_barrier <- FocalAllocation(elev_bar_dir_r)

  # Assign barrier distance value by Euclidean direction
  dist_bar_dir <- elev_dist_mean$dist[match(eucdir, elev_dist_mean$eucdir)]
  dist_bar_dir_r <- rast(to_rast(sgrid), val = dist_bar_dir, names = "dist_bar_dir")
  dist_barrier <- FocalAllocation(dist_bar_dir_r)

  # Calculate slope between source and receiver
  slope <- (as.vector(dem_ft) - elev) / eucdist_ft
  slope <- ifelse(is.infinite(slope), NA, slope)

  # Calculate elevation of source-receiver line under barrier
  elev_sr <- slope * as.vector(dist_barrier) + elev

  # Calculate barrier height
  h_b <- as.vector(elev_barrier) - elev_sr

  # Reclassify negative barrier height values to zero
  h_b_rc <- ifelse(h_b >= -10000 & h_b <= 0, 0, h_b) # Debe ser solo <= 0

  # Calculate barrier path distance (BPD)
  term1 <- sqrt(h_b_rc**2 + as.vector(dist_barrier)**2)
  term2 <- sqrt(h_b_rc**2 + (eucdist_ft - as.vector(dist_barrier))**2)
  calc_result <- term1 + term2 - eucdist_ft

  # Calculate barrier factor (N)
  L <- 0.0000000000005*freq_s**4 - 0.000000001*freq_s**3 - 0.0000004*freq_s**2 + 0.0028*freq_s - 0.3051
  bar_factor <- calc_result * L

  # The coefficients appear to be from a power regression based on Table 14.
  # Note that I got 13.451 and 0.2427 for my coefficients based on medians.
  # Calculate barrier loss
  bar = (13.573) * (bar_factor**0.2299)

  return(bar)
}

calculate_topozones <- function(srcv_vec, dem_ft, ground){

  # Alternate approach using observer points code
  # potential for efficiency gain as multiple points can be calculated at once.
  dem_file <- tempfile(fileext = ".tif")
  pnt_file <- tempfile(fileext = ".shp")
  view_file <- tempfile(fileext = ".tif")

  writeRaster(dem_ft*3.28084, dem_file)
  st_write(srcv_vec, pnt_file)

  wbt_viewshed(dem_file, pnt_file, view_file, height = 0)
  viewshed <- rast(view_file)

  # Define areas where ground effects and atmospheric effects dominate
  ground_atmos <- max(ground, viewshed)

  # Define areas where ground, atmospheric, and barrier effects dominate
  topo_zones <- app(ground_atmos, function(x) ifelse(is.na(x), 2, ifelse(x == 0, 2, x)))

  return(as.vector(topo_zones))
}

compute_noise_propagation <- function(source_level, ssl, aal, veg, wind, bar,
                                      topo_zones){
  # Subtract spherical spreading loss from source sound level
  sslloss <- source_level - ssl

  # Calculate cumulative spherical spreading and atmospheric absorption loss
  sslaal <- sslloss + aal

  # Calculate cumulative spherical spreading, atmospheric absorption, and vegetation loss
  salveg <- sslaal + veg

  # Calculate cumulative spherical spreading, atmospheric absorption, vegetation, and wind loss
  salvegwin <- salveg + wind

  # Calculate cumulative spherical spreading, atmospheric absorption, vegetation, wind, and barrier loss
  salvgwnbr <- salvegwin + bar

  # Pick noise propagation values for areas where ground, atmospheric,
  # or barrier effects dominate
  combined <- mapply(function(x, y, z, w){
    df <- data.frame(x = x, y = y, z = z)
    val <- rep(NA, length(x))
    for(i in 1:length(x)){
      val[i] <- df[i, w[i]]
    }
    return(val)
  },x=sslaal, y=salvgwnbr, z=salvegwin, w=as.vector(topo_zones))

  # Smooth noise propagation patterns
  smoothed <- terra::focal(rast(grid_rst, val = combined), matrix(1,nrow=3,ncol=3),
                           fun = "mean", na.rm = T)

  # Patch to prevent smoothing at cell of origin
  pr <- ifelse(source_level > 0, combined, as.vector(smoothed))

  return(list("sslloss" = sslloss, "sslaal" = sslaal, "salveg" = salveg,
              "salvegwin" = salvegwin, "salvgwnbr" = salvgwnbr, "combined" = combined,
              "smoothed" = smoothed, "pr" = pr))
}

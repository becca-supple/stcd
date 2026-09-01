#' Simulate spatially-confounded systems of time series from a linear structural
#' causal model
#'
#' @param mask SpatVector of polygons of countries/boundaries
#' @param area_ID Name of object that stores the area IDs in the mask
#' @param t_steps Number of time steps
#' @param autocorrelation Autocorrelation between years, rho
#' @param vgm_data Variogram model type for data - "Exp", "Sph", or "Gau"
#' @param nugget Variance at 0 distance
#' @param sill Asymptotic max variance
#' @param range_percent Percent of mask radius at which 95% of max variance
#' reached
#' @param spat_conf Coefficient for effect of spatial confounding
#' @param k Neighbors to consider when simulating
#' @param n_sample Number of points to sample for simulation
#' @param strength Effect size strength
#' @param noise Variance of non-spatial random Gaussian noise with mean 0
#' @param return_rast If TRUE, return raw raster data
#' @param discrete Set = TRUE for discrete space, FALSE for continuous
#' @param n_sites Number of sites to simulate. Set to NULL to use number of
#' sub-polygons in mask.
#'
#' @returns Data frame of time series at locations, and rasters if return_rast = TRUE
#' @export
#'
#' @examples
sim_system_linear <- function(mask,
                              area_ID,
                              t_steps,
                              autocorrelation = 0.1,
                              vgm_data = "Gau",
                              nugget = 0,
                              sill = 1,
                              range_percent = 50,
                              spat_conf = 1,
                              k = 4,
                              n_sample = 2500,
                              strength = c("none", "low", "med", "high"),
                              noise = 0.1,
                              return_rast = FALSE,
                              discrete = TRUE,
                              n_sites = NULL
){

  ###
  # Prepare geodata
  ###

  #sample points from the vector
  pointsample <- terra::spatSample(mask, size = n_sample, method = "regular")

  #get coordinates
  xy <- as.data.frame(crds(pointsample, df = TRUE))

  #transform range percent value
  area <- sum(expanse(mask, unit = "m"))
  range <- (range_percent/100)*(sqrt(area/pi))

  ###
  # Simulate U
  ###

  #Set up blank list
  U <- list()

  U_field <- sim_field(xy = xy,
                       vgm_model = vgm_data,
                       nugget = nugget,
                       sill = sill+5,
                       range = range,
                       k = k,
                       n_sample = n_sample)

  n_locations <- nrow(U_field)
  xy <- U_field[, c(1, 2)]

  U[[1]] <- U_field[["sim1"]]

  for(t in seq(2, t_steps + 1)){

    U[[t]] <- numeric(n_locations)

    epsilon <- sim_field(xy = xy,
                         vgm_model = vgm_data,
                         nugget = nugget,
                         sill = sill,
                         range = range,
                         k = k,
                         n_sample = n_sample)[["sim1"]]

    U[[t]] <- autocorrelation * U[[t-1]] + sqrt(1 - autocorrelation^2) * epsilon

  }

  ###
  # Simulate Agriculture
  ###

  #Set up blank list
  agr <- list()

  #Initialize time series
  agr[[1]] <- spat_conf * U[[1]] +
    rnorm(n_locations, 0, noise)

  for(t in seq(2, t_steps + 1)){

    agr[[t]] <- numeric(n_locations)

    agr[[t]] <- autocorrelation * agr[[t-1]] + spat_conf * U[[t]] + rnorm(n_locations, 0, noise)

  }

  #Transform data
  agr <- lapply(agr, pnorm) #probit link to make data in [0,1]


  ###
  # Simulate Fencing
  ###

  #Set up blank list
  fence <- list()

  #Vector of parameters
  beta <- ifelse(strength == "none", 0,
                 ifelse(strength == "low", 0.5,
                        ifelse(strength == "med", 1, 2)))

  #Initialize time series
  fence[[1]] <- beta * agr[[1]] + spat_conf * U[[1]] +
    rnorm(n_locations, 0, noise)

  for(t in seq(2, t_steps + 1)){

    fence[[t]] <- numeric(n_locations)

    fence[[t]] <- autocorrelation * fence[[t-1]] +
      beta * agr[[t]] + spat_conf * U[[t]] + rnorm(n_locations, 0, noise)

  }

  #center and scale data
  fence <- lapply(fence, function(i){
    (i - mean(i))/(sd(i))
  })

  ###
  # Simulate Fertilizer
  ###

  #Set up blank list
  fert <- list()

  #Initialize time series
  fert[[1]] <- beta * agr[[1]] + spat_conf * U[[1]] +
    rnorm(n_locations, 0, noise)

  for(t in seq(2, t_steps + 1)){

    fert[[t]] <- numeric(n_locations)

    fert[[t]] <- autocorrelation * fert[[t-1]] +
      beta * agr[[t]] + spat_conf * U[[t]] + rnorm(n_locations, 0, noise)

  }

  #center and scale data
  fert <- lapply(fert, function(i){
    (i - mean(i))/(sd(i))
  })

  ###
  # Simulate Occupancy
  ###

  #Set up blank list
  occ <- list()

  #Initialize field
  occ[[1]] <- beta * agr[[1]] + spat_conf * U[[1]] +
    rnorm(n_locations, 0, noise)

  #Sim forward
  for(t in seq(2, t_steps + 1)){

    occ[[t]] <- numeric(n_locations)

    occ[[t]] <- autocorrelation * occ[[t-1]] +
      beta * agr[[t]] +
      (-1) * beta * fert[[t-1]] + spat_conf * U[[1]] +
      rnorm(n_locations, 0, noise)
  }

  #transform
  occ <- lapply(occ, pnorm)


  ###
  # Summarize by Area if Discrete
  ###

  if(discrete){

    if(return_rast){
      rast_stack <- list()
    }

    data_short <- data.frame(area_ID = unique(mask[[area_ID]]))


    ##Summarize agriculture##

    for(t in seq(2, t_steps + 1)){

      #attach data to coordinates
      spat_sim <- cbind(xy, agr[[t]])

      #make raster
      rast_sim <- rast(spat_sim, type = "xyz", crs = crs(mask))
      rast_sim <- mask(rast_sim, mask)

      if(return_rast){
        rast_stack[[t-1]] <- c(rast_sim)
      }

      #extract mean by country
      country_means <- terra::extract(rast_sim, mask, fun = mean, na.rm = TRUE)

      data_short[[paste0("T", t - 1)]] <- country_means[,2]

    }

    #Pivot data
    data_long <- pivot_longer(data_short,
                              cols = paste0("T", seq(1, t_steps)),
                              names_to = "Time",
                              values_to = "Agr")

    ##Summarize fencing##

    for(t in seq(2, t_steps + 1)){

      #attach data to coordinates
      spat_sim <- cbind(xy, fence[[t]])

      #make raster
      rast_sim <- rast(spat_sim, type = "xyz", crs = crs(mask))
      rast_sim <- mask(rast_sim, mask)

      if(return_rast){
        rast_stack[[t-1]] <- c(rast_stack[[t-1]], rast_sim)
      }

      #extract mean by country
      country_means <- terra::extract(rast_sim, mask, fun = mean, na.rm = TRUE)

      data_short[[paste0("T", t - 1)]] <- country_means[,2]

    }

    #Pivot data
    fence_long <- pivot_longer(data_short,
                               cols = paste0("T", seq(1, t_steps)),
                               names_to = "Time",
                               values_to = "Fence")

    data_long[["Fence"]] <- fence_long[["Fence"]]

    ##Summarize fertilizer##

    for(t in seq(2, t_steps + 1)){

      #attach data to coordinates
      spat_sim <- cbind(xy, fert[[t]])

      #make raster
      rast_sim <- rast(spat_sim, type = "xyz", crs = crs(mask))
      rast_sim <- mask(rast_sim, mask)

      if(return_rast){
        rast_stack[[t-1]] <- c(rast_stack[[t-1]], rast_sim)
      }

      #extract mean by country
      country_means <- terra::extract(rast_sim, mask, fun = mean, na.rm = TRUE)

      data_short[[paste0("T", t - 1)]] <- country_means[,2]

    }

    #Pivot data
    fert_long <- pivot_longer(data_short,
                              cols = paste0("T", seq(1, t_steps)),
                              names_to = "Time",
                              values_to = "Fert")

    data_long[["Fert"]] <- fert_long[["Fert"]]

    ##Summarize occupancy##

    for(t in seq(2, t_steps + 1)){

      #attach data to coordinates
      spat_sim <- cbind(xy, occ[[t]])

      #make raster
      rast_sim <- rast(spat_sim, type = "xyz", crs = crs(mask))
      rast_sim <- mask(rast_sim, mask)

      if(return_rast){
        rast_stack[[t-1]] <- c(rast_stack[[t-1]], rast_sim)
        names(rast_stack[[t-1]]) <- c("Agr", "Fence", "Fert", "Occ")
      }

      #extract mean by country
      country_means <- terra::extract(rast_sim, mask, fun = mean, na.rm = TRUE)

      data_short[[paste0("T", t - 1)]] <- country_means[,2]

    }

    #Pivot data
    occ_long <- pivot_longer(data_short,
                             cols = paste0("T", seq(1, t_steps)),
                             names_to = "Time",
                             values_to = "Occ")

    data_long[["Occ"]] <- occ_long[["Occ"]]

    ###
    # Format and return
    ###

    data_long <- data_long |>
      filter(!is.na(Agr)) |> #remove countries with NA (out of frame)
      mutate(Time = as.numeric(gsub("T", "", Time))) #numeric time steps

    if(return_rast){
      list2return <- list(data = data_long, rasters = rast_stack)
      return(list2return)
    }else{return(data_long)}

  }

  ###
  # Or sample sites and format for coordinate space
  ###

  else{

    if(return_rast){
      rast_stack <- list()

      for(t in seq(2, t_steps + 1)){

        #attach data to coordinates
        spat_sim <- cbind(xy, agr[[t]])

        #make raster
        rast_sim <- rast(spat_sim, type = "xyz", crs = crs(mask))
        rast_sim <- mask(rast_sim, mask)

        #add to stack
        rast_stack[[t-1]] <- c(rast_sim)

        #attach data to coordinates
        spat_sim <- cbind(xy, fence[[t]])

        #make raster
        rast_sim <- rast(spat_sim, type = "xyz", crs = crs(mask))
        rast_sim <- mask(rast_sim, mask)

        #add to stack
        rast_stack[[t-1]] <- c(rast_stack[[t-1]], rast_sim)

        #attach data to coordinates
        spat_sim <- cbind(xy, fert[[t]])

        #make raster
        rast_sim <- rast(spat_sim, type = "xyz", crs = crs(mask))
        rast_sim <- mask(rast_sim, mask)

        #add to stack
        rast_stack[[t-1]] <- c(rast_stack[[t-1]], rast_sim)

        #attach data to coordinates
        spat_sim <- cbind(xy, occ[[t]])

        #make raster
        rast_sim <- rast(spat_sim, type = "xyz", crs = crs(mask))
        rast_sim <- mask(rast_sim, mask)

        #add to stack
        rast_stack[[t-1]] <- c(rast_stack[[t-1]], rast_sim)

      }

    }

    if(is.null(n_sites)){
      n_sites <- length(unlist(unique(mask[[area_ID]])))
    }

    sites <- sample(seq(1, length(agr[[1]])), n_sites, replace = FALSE)

    long_data <- as.data.frame(expand.grid(sites, seq(1, t_steps)))
    colnames(long_data) <- c("site_ID", "Time")

    long_data$Agr <- numeric(nrow(long_data))
    long_data$Fence <- numeric(nrow(long_data))
    long_data$Fert <- numeric(nrow(long_data))
    long_data$Occ <- numeric(nrow(long_data))
    long_data$x <- numeric(nrow(long_data))
    long_data$y <- numeric(nrow(long_data))

    for(i in seq(1, nrow(long_data))){

      site <- long_data$site_ID[i]
      y <- long_data$Time[i] + 1

      long_data$Agr[i] <- agr[[y]][site]
      long_data$Fence[i] <- fence[[y]][site]
      long_data$Fert[i] <- fert[[y]][site]
      long_data$Occ[i] <- occ[[y]][site]
      long_data$x[i] <- xy$x[site]
      long_data$y[i] <- xy$y[site]

    }

    if(return_rast){
      list2return <- list(data = long_data, rasters = rast_stack)
      return(list2return)
    }else{return(long_data)}

  }

}

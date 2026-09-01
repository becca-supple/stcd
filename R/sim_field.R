#' Wrapper for gstat to simulate from a variogram model
#'
#' @param xy Data frame of xy coordinates of point samples
#' @param vgm_model Variogram model type - "Exp", "Sph", or "Gau"
#' @param nugget Variance at 0 distance
#' @param sill Asymptotic max variance
#' @param range Distance at which 95% of max variance reached (m)
#' @param k Number of nearest neighbors to consider when simulating
#' @param n_sample Number of points to sample for simulation
#'
#' @returns A data frame of samples from the spatial model with values in a
#' column named "sim1"
#' @export
#'
#' @examples
#'
#' points <- as.data.frame(unique(example_coord_data[,c("x", "y")]))
#' ex_field <- sim_field(xy = points,
#'              vgm_model = "Gau",
#'              nugget = 0,
#'              sill = 1,
#'              range = 10000,
#'              k = 4,
#'              n_sample = 1000)
#' head(ex_field)
#'
sim_field <- function(xy,
                      vgm_model = "Gau",
                      nugget = 0,
                      sill = 1,
                      range,
                      k = 4,
                      n_sample = 1000

){

  ###
  # Define spatial model and predict
  ###

  spat_model <- gstat::gstat(formula = z ~ 1, #single variable no dependence
                             locations = ~x+y,
                             dummy = TRUE, #unconditional simulation
                             beta = nugget,
                             model = gstat::vgm(model = vgm_model,
                                         psill = sill,
                                         range = range
                             ),
                             nmax = k
  )

  spat_sim <- predict(spat_model,
                      newdata = xy,
                      nsim = 1,
                      debug.level = 0)

  return(spat_sim)
}

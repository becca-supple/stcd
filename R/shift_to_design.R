#' Create shifted time series design matrices from data
#'
#' @param data `data.frame` or `list` of \link[terra:SpatRaster-class]{SpatRaster} objects
#'
#'    * `data.frame` objects must have columns for each time series variable in the system, (a) column(s) for spatial structure, and a column indexing the time.
#'
#'    * A `list` of `SpatRaster` objects must be in time order, with identically named data layers for each `SpatRaster` object
#' @param time_col String naming the column indexing time
#' @param max_lag Integer value representing the maximum time lag at which you expect relationships to exist. In practice, the function will shift time series by 2 times the `max_lag` for lagged relationship discovery.
#' @param ... other arguments to specify spatial structure, see instructions below.
#'
#' @returns a `dsgn` object (either `discrete_dsgn` or `coord_dsgn` depending on spatial structure)
#'
#' @section Creating discrete space design matrices:
#'
#' Time series indexed by discrete space (sites, counties, countries, etc) should be shifted to
#' discrete design matrices. Discrete space design matrices (objects of class `discrete_dsgn`)
#' can only be created from `data.frame` objects. This is done by setting the `site_col` argument to
#' a string naming the column containing site information.
#'
#' @section Creating coordinate/continuous space design matrices:
#'
#' Time series indexed by coordinate/continuous space in 2 dimensions (x-y, latitude-longitude,
#' easting-northing, etc) should be shifted to coordinate design matrices. Coordinate space design
#' matrices (objects of class `coord_dsgn`) can be created from either a `data.frame` object or a
#' `list` of \link[terra:SpatRaster-class]{SpatRaster} objects:
#'
#' * Create a `coord_dsgn` object from a `data.frame` by setting the `coords` argument to a vector
#' of strings indicating the two column names containing coordinate information
#' (e.g. `coords = c("x", "y")`).
#'
#' * Creation of a design matrix from a `list` of `SpatRaster` objects will only ever create a
#' `coord_dsgn` object. You will need to define the number of samples to take from the rasters with
#' the `n_sample` argument (by default, `n_samples = 100`). This is a random sample, but you could set
#' a seed to ensure reproducibility.
#'
#'
#' @export
#'
#' @examples
#' # The Chloris_chloris_data set has a discrete spatial structure with a single
#' # column designating the site ID (`CountryGroup`)
#' head(Chloris_chloris_data)
#'
#' # Make sure to set the `site_col` argument to create a `discrete_dsgn` object
#' CC_design <- shift_to_design(data = Chloris_chloris_data,
#'                 time_col = "Year",
#'                 max_lag = 2,
#'                 site_col = "CountryGroup")
#' head(CC_design)
#' class(CC_design)
shift_to_design <- function(data,
                            time_col = "Time",
                            max_lag = 2,
                            ...){

  UseMethod("shift_to_design", data)

}

#' @export
shift_to_design.data.frame <- function(data,
                                       site_col = NULL,
                                       coords = NULL,
                                       time_col = "Time",
                                       max_lag = 2){

  if(is.null(coords) && is.null(site_col)){
    stop("\nBoth site_col and coords missing. Specify either site_col for discrete or coords for coordinate/continuous space.")
  }

  if(is.null(coords) == FALSE){

    #unique site IDs from coordinates
    data$site_ID <- stringr::str_c(data[[coords[1]]], data[[coords[2]]], sep = "_")
    n <- sum(data$site_ID == unique(data$site_ID[1]))
    num_rows <- n - 2 * max_lag #number of rows in design per site_col level

    if(num_rows < 1){
      stop(paste0("\nTime series of length ", n,
                  " too short to shift with maximum lag ", max_lag))
    }

    n_X <- ncol(data) - 4
    X <- seq(1, n_X)

    #make sure time and site columns are at end
    data <- data |>
      dplyr::select(-c(all_of("site_ID"), all_of(time_col), all_of(coords)), c(all_of(coords), all_of("site_ID"), all_of(time_col)))

    design <- data.frame(matrix(
      nrow = (num_rows)*length(unique(data$site_ID)),
      ncol = (2*max_lag + 1)*n_X + 2))

    colnames(design)[((2*max_lag + 1)*n_X + 1):((2*max_lag + 1)*n_X + 2)] <- coords
    col_index <- 0

    #Name the columns
    for(d in X){

      variable <- colnames(data)[d]
      for(k in seq(0, 2*max_lag)){
        if(k > 0){
          colnames(design)[col_index + k + 1] <- paste(variable, "_tminus", k, sep = "")
        } else{
          colnames(design)[col_index + k + 1] <- variable
        }
      }
      col_index <- col_index + (2*max_lag + 1)
    }

    #Fill in data values
    start <- 1
    for(fac in unique(data$site_ID)){
      filtered <- dplyr::filter(data, .data$site_ID == fac) #filter time series
      x_value <- unique(dplyr::filter(data, .data$site_ID == fac)[[coords[1]]])
      y_value <- unique(dplyr::filter(data, .data$site_ID == fac)[[coords[2]]])
      col_index <- 0

      for(d in X){
        for(k in seq(0, 2*max_lag)){ #fill in data
          design[start:(start + num_rows - 1), col_index + k + 1] <- filtered[seq(
            2*max_lag + 1 - k, n - k), d]
        }
        col_index <- col_index + (2*max_lag + 1)
      }

      #label by coordinates
      design[[coords[1]]][start:(start + num_rows - 1)] <- rep(x_value, num_rows)
      design[[coords[2]]][start:(start + num_rows - 1)] <- rep(y_value, num_rows)

      #shift the data entry to account for previous entered values
      start <- start + num_rows
    }

    #Set class
    design <- coord_dsgn(design)

  }

  if(is.null(site_col) == FALSE){
    #Create a version of data with shifted versions of all time series
    n <- sum(data[[site_col]] == unique(data[[site_col]][1]))
    num_rows <- n - 2 * max_lag #number of rows in design per site_col level
    n_X <- ncol(data) - 2
    X <- seq(1, n_X)

    #make sure time and site columns are at end
    data <- data |>
      dplyr::select(-c(all_of(site_col), all_of(time_col)), c(all_of(site_col), all_of(time_col)))

    design <- data.frame(matrix(
      nrow = (num_rows)*length(unique(data[[site_col]])),
      ncol = (2*max_lag + 1)*n_X + 1))

    colnames(design)[(2*max_lag + 1)*n_X + 1] <- site_col
    col_index <- 0

    #Name the columns
    for(d in X){

      variable <- colnames(data)[d]
      for(k in seq(0, 2*max_lag)){
        if(k > 0){
          colnames(design)[col_index + k + 1] <- paste(variable, "_tminus", k, sep = "")
        } else{
          colnames(design)[col_index + k + 1] <- variable
        }
      }
      col_index <- col_index + (2*max_lag + 1)
    }

    #Fill in data values
    start <- 1
    for(fac in unique(data[[site_col]])){
      filtered <- dplyr::filter(data, .data[[site_col]] == fac) #filter time series
      col_index <- 0

      for(d in X){
        for(k in seq(0, 2*max_lag)){ #fill in data
          design[start:(start + num_rows - 1), col_index + k + 1] <- filtered[seq(
            2*max_lag + 1 - k, n - k), d]
        }
        col_index <- col_index + (2*max_lag + 1)
      }

      #label by factor value
      design[[site_col]][start:(start + num_rows - 1)] <- rep(fac, num_rows)

      #shift the data entry to account for previous entered values
      start <- start + num_rows
    }

    #site_col variable as a factor
    design[[site_col]] <- factor(design[[site_col]])

    #Set class
    design <- discrete_dsgn(design)
  }

  return(design)

}

#' @export
shift_to_design.list <- function(data,
                                 time_col = "Time",
                                 max_lag = 2,
                                 n_samples = 100){

  #data frame from first year rasters
  df <- terra::as.data.frame(data[[1]], xy = TRUE)
  df$Time <- 1

  #sample random sites & filter
  sites <- sample(seq(1, nrow(df)), size = n_samples, replace = FALSE)
  df <- df[sites,]

  #add other years
  for(t in seq(2, length(data))){

    temp <- terra::as.data.frame(data[[t]], xy = TRUE)
    temp$Time <- t
    temp <- temp[sites,]

    df <- rbind(df, temp)

  }

  design <- shift_to_design.data.frame(df,
                                       coords = c("x", "y"),
                                       time_col = "Time")

}
